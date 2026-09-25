"""LCU/block-encoding assembly for quantum-chemistry low-rank decomposition (DF/THC) Hamiltonians.

Implements the "low-rank tensor -> LCU -> BE" pipeline of the double factorization (DF)
of Berry et al. 2019 and von Burg et al. 2021, and of the tensor hypercontraction
(THC) of Lee et al. 2021, corresponding to item C1 of the algorithm coverage board.
Integral tensors are given as classical input data; this module performs no real
chemical integral computation. The explicit small-matrix path targets small
validation instances; at large scale the spectra and rotation angle tables of
``U_r``/``G_r`` should instead go through QRAM data binding (see the three-layer
qram_prepare/alias_prepare paradigm of prepare_select).
"""

from __future__ import annotations

import cmath
import contextlib
import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from oracq.algorithms.input_model.block_encoding import adjoint_be, lcu, matrix_pauli_encoding
from oracq.algorithms.input_model.contracts import finite_real, require_instance
from oracq.algorithms.input_model.operators import BlockEncoding, _name, identity, product
from oracq.algorithms.input_model.oracles import annotate, invoke, resources_for
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits, Ref, ValidationError, fuse

_SYNTHESIS_TOLERANCE = 1e-15
_MATRIX_LIMIT_QUBITS = 5


def _as_complex_matrix(
    value: Iterable[Iterable[complex]], path: str
) -> tuple[tuple[complex, ...], ...]:
    """Normalize the input into a tuple-of-rows square complex matrix with finite entries, with dimension restricted to powers of two in 2..32."""
    try:
        matrix = tuple(tuple(complex(v) for v in row) for row in value)
    except TypeError as exc:
        raise ValidationError(f"{path} entries must be numeric") from exc
    if not matrix or any(len(row) != len(matrix) for row in matrix):
        raise ValidationError(f"{path} must be a non-empty square matrix")
    d = len(matrix)
    if d < 2 or d & (d - 1) or d > (1 << _MATRIX_LIMIT_QUBITS):
        raise ValidationError(f"{path} dimension must be a power of two in 2..32")
    if not all(math.isfinite(v.real) and math.isfinite(v.imag) for row in matrix for v in row):
        raise ValidationError(f"{path} entries must be finite")
    return matrix


def _matmul(
    a: Sequence[Sequence[complex]], b: Sequence[Sequence[complex]]
) -> tuple[tuple[complex, ...], ...]:
    """Compute the matrix product of two square complex matrices, returned as a row-major nested tuple."""
    return tuple(
        tuple(sum(ar[k] * b[k][j] for k in range(len(b))) for j in range(len(b[0]))) for ar in a
    )


def _check_unitary(
    matrix: Sequence[Sequence[complex]], path: str, *, tolerance: float = 1e-9
) -> None:
    """Check the column orthonormality of a square matrix; raise ``ValidationError`` when it fails."""
    d = len(matrix)
    for i in range(d):
        for j in range(d):
            value = sum(matrix[k][i].conjugate() * matrix[k][j] for k in range(d))
            if abs(value - (1 if i == j else 0)) > tolerance:
                raise ValidationError(f"{path} must be a unitary matrix with orthonormal columns")


def diagonalize_symmetric(
    matrix: Iterable[Iterable[float]], *, tolerance: float = 1e-12, max_sweeps: int = 100
) -> tuple[tuple[float, ...], tuple[tuple[float, ...], ...]]:
    """Jacobi eigendecomposition of a real symmetric matrix: returns eigenvalues and the eigenvector matrix satisfying ``G = V diag(λ) Vᵀ``.

    This is the classical preprocessing of the DF input model: after the symmetric
    matrix ``G_r`` is diagonalized, the eigenvector matrix is folded into the
    rotation ``U_r`` by ``from_symmetric``. The columns of the eigenvector matrix
    are the eigenvectors.

    Args:
        matrix: Real symmetric square matrix with finite real entries.
        tolerance: Tolerance shared by the symmetry check and the convergence test.
        max_sweeps: Maximum number of Jacobi sweeps; an error is raised if convergence is not reached in time.

    Returns:
        tuple[tuple[float, ...], tuple[tuple[float, ...], ...]]:
        ``(eigenvalue tuple, matrix whose columns are the eigenvectors)``, satisfying ``G = V diag(λ) Vᵀ``.
    """
    raw = tuple(tuple(row) for row in matrix)
    if not raw or any(len(row) != len(raw) for row in raw):
        raise ValidationError("diagonalize_symmetric.matrix must be a non-empty square matrix")
    d = len(raw)
    for row in raw:
        for value in row:
            finite_real(value, "diagonalize_symmetric.matrix")
    a = [list(map(float, row)) for row in raw]
    for p in range(d):
        for q in range(p + 1, d):
            if abs(a[p][q] - a[q][p]) > tolerance:
                raise ValidationError("diagonalize_symmetric.matrix must be a real symmetric matrix")
    v = [[1.0 if i == j else 0.0 for j in range(d)] for i in range(d)]
    for _ in range(max_sweeps):
        off = max((abs(a[p][q]) for p in range(d) for q in range(p + 1, d)), default=0.0)
        if off <= tolerance:
            break
        for p in range(d - 1):
            for q in range(p + 1, d):
                if abs(a[p][q]) <= tolerance:
                    continue
                theta = (a[q][q] - a[p][p]) / (2 * a[p][q])
                t = (1 if theta >= 0 else -1) / (abs(theta) + math.hypot(1.0, theta))
                c = 1.0 / math.sqrt(1.0 + t * t)
                s = t * c
                for k in range(d):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p] = c * akp - s * akq
                    a[k][q] = s * akp + c * akq
                for k in range(d):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k] = c * apk - s * aqk
                    a[q][k] = s * apk + c * aqk
                for k in range(d):
                    vkp, vkq = v[k][p], v[k][q]
                    v[k][p] = c * vkp - s * vkq
                    v[k][q] = s * vkp + c * vkq
    else:
        raise ValidationError("Jacobi eigendecomposition did not converge within max_sweeps")
    return tuple(a[k][k] for k in range(d)), tuple(tuple(row) for row in v)


@dataclass(frozen=True)
class DoubleFactorization:
    """DF input model (CP input model): small-instance representation of ``H = scalar·I + Σ_r U_r diag(g_r) U_r†``.

    rotations holds explicit small unitary matrices and spectra holds the real
    eigenvalues obtained by classically diagonalizing the symmetric matrices
    ``G_r``; in physical DF the ``G_r`` are received by ``from_symmetric`` and
    folded into the rotations. The explicit-matrix path targets small validation
    instances only; at large scale the spectral norm tables and rotation angles
    go through QRAM data binding.
    """

    scalar: float
    rotations: tuple[tuple[tuple[complex, ...], ...], ...]
    spectra: tuple[tuple[float, ...], ...]

    def __post_init__(self) -> None:
        """Validate the scalar, the unitarity of each rank rotation, and the spectrum lengths, and normalize the storage format."""
        finite_real(self.scalar, "DoubleFactorization.scalar")
        rotations = tuple(
            _as_complex_matrix(r, f"DoubleFactorization.rotations[{i}]")
            for i, r in enumerate(self.rotations)
        )
        if not rotations:
            raise ValidationError("DoubleFactorization requires at least one rank term")
        if len(rotations) != len(self.spectra):
            raise ValidationError("DoubleFactorization has mismatched counts of rotations and spectra")
        d = len(rotations[0])
        if any(len(r) != d for r in rotations):
            raise ValidationError("DoubleFactorization rank rotations have inconsistent dimensions")
        for i, rotation in enumerate(rotations):
            _check_unitary(rotation, f"DoubleFactorization.rotations[{i}]")
        spectra: list[tuple[float, ...]] = []
        for i, spectrum in enumerate(self.spectra):
            values = tuple(spectrum)
            if len(values) != d:
                raise ValidationError(f"DoubleFactorization.spectra[{i}] length must equal dimension {d}")
            for value in values:
                finite_real(value, f"DoubleFactorization.spectra[{i}]")
            spectra.append(tuple(float(v) for v in values))
        object.__setattr__(self, "scalar", float(self.scalar))
        object.__setattr__(self, "rotations", rotations)
        object.__setattr__(self, "spectra", tuple(spectra))

    @classmethod
    def from_symmetric(
        cls,
        scalar: float,
        rotations: Iterable[Iterable[Iterable[complex]]],
        factors: Iterable[Iterable[Iterable[float]]],
    ) -> DoubleFactorization:
        """Construct from explicit unitaries ``U_r`` and real symmetric ``G_r``: ``G_r = V_r diag(g_r) V_rᵀ`` is folded into ``U_r V_r``.

        Args:
            scalar: Coefficient of the identity term.
            rotations: Explicit unitary matrices ``U_r`` of each rank.
            factors: Real symmetric matrices ``G_r`` paired one-to-one with the rotations, of matching dimension.

        Returns:
            DoubleFactorization: DF input model with the spectrum unchanged after folding ``U_r V_r``.
        """
        rotations = tuple(rotations)
        factors = tuple(factors)
        if len(rotations) != len(factors):
            raise ValidationError("from_symmetric has mismatched counts of rotations and factors")
        combined: list[tuple[tuple[complex, ...], ...]] = []
        spectra: list[tuple[float, ...]] = []
        for rotation, factor in zip(rotations, factors, strict=True):
            unitary = _as_complex_matrix(rotation, "from_symmetric.rotations")
            eigenvalues, vectors = diagonalize_symmetric(factor)
            if len(vectors) != len(unitary):
                raise ValidationError("from_symmetric has inconsistent rotation and symmetric matrix dimensions")
            combined.append(_matmul(unitary, vectors))
            spectra.append(eigenvalues)
        return cls(scalar, tuple(combined), tuple(spectra))

    @property
    def width(self) -> int:
        """Number of target qubits."""
        return (len(self.rotations[0]) - 1).bit_length()

    @property
    def rank(self) -> int:
        """Number of DF rank terms."""
        return len(self.rotations)


@dataclass(frozen=True)
class THCDecomposition:
    """THC input model: ``H = Σ_{μν} ζ_{μν} L_μ L_ν†`` where ``L_μ`` are explicit small-matrix leaf operators.

    coefficients is the real symmetric ``ζ`` matrix and leaves holds explicit
    small matrices of matching dimension; the leaf operators need not be unitary
    or Hermitian, and the block encoding of the product ``L_μ L_ν†`` is assembled
    from matrix_pauli_encoding and BE products.
    """

    coefficients: tuple[tuple[float, ...], ...]
    leaves: tuple[tuple[tuple[complex, ...], ...], ...]

    def __post_init__(self) -> None:
        """Validate the real symmetry of the ζ matrix, dimensional consistency, and the leaf operator format."""
        zeta = tuple(tuple(row) for row in self.coefficients)
        if not zeta or any(len(row) != len(zeta) for row in zeta):
            raise ValidationError("THCDecomposition.coefficients must be a non-empty square matrix")
        for row in zeta:
            for value in row:
                finite_real(value, "THCDecomposition.coefficients")
        size = len(zeta)
        for mu in range(size):
            for nu in range(mu + 1, size):
                if abs(zeta[mu][nu] - zeta[nu][mu]) > 1e-12:
                    raise ValidationError("THCDecomposition.coefficients must be a real symmetric matrix")
        if len(self.leaves) != size:
            raise ValidationError("THCDecomposition leaves count must match the ζ dimension")
        leaves = tuple(
            _as_complex_matrix(leaf, f"THCDecomposition.leaves[{i}]")
            for i, leaf in enumerate(self.leaves)
        )
        if len({len(leaf) for leaf in leaves}) != 1:
            raise ValidationError("THCDecomposition leaf operators have inconsistent dimensions")
        object.__setattr__(self, "coefficients", zeta)
        object.__setattr__(self, "leaves", leaves)

    @property
    def width(self) -> int:
        """Number of target qubits."""
        return (len(self.leaves[0]) - 1).bit_length()

    @property
    def leaf_count(self) -> int:
        """Number of THC leaf operators."""
        return len(self.leaves)


def _emit_transposition(builder: Builder, ref: Ref, first: int, second: int) -> None:
    """Swap ``|first>`` and ``|second>`` along a Gray-code path with multi-controlled X gates, leaving other basis states unchanged."""
    path = [first]
    for bit in range(ref.width):
        if ((first ^ second) >> bit) & 1:
            path.append(path[-1] ^ (1 << bit))
    edges = list(zip(path, path[1:], strict=False))
    for left, _right in edges + list(reversed(edges[:-1])):
        bit = (left ^ _right).bit_length() - 1
        controls = fuse(ref[:bit], ref[bit + 1 :])
        value = (left & ((1 << bit) - 1)) | ((left >> (bit + 1)) << bit)
        with builder.control(controls, value) if controls.width else contextlib.nullcontext():
            builder.x(ref[bit])


def _zyz(matrix: Sequence[Sequence[complex]]) -> tuple[float, float, float, float]:
    """Exact ``e^{iφ} Rz(α) Ry(β) Rz(γ)`` decomposition of a second-order unitary (simulator convention in execution.gate_matrix)."""
    a00, a01, a10, a11 = matrix[0][0], matrix[0][1], matrix[1][0], matrix[1][1]
    beta = 2 * math.atan2(abs(a10), abs(a00))
    if abs(a10) <= _SYNTHESIS_TOLERANCE:
        phi = (cmath.phase(a00) + cmath.phase(a11)) / 2
        return phi, cmath.phase(a11) - cmath.phase(a00), 0.0, 0.0
    if abs(a00) <= _SYNTHESIS_TOLERANCE:
        lower, upper = cmath.phase(a10), cmath.phase(-a01)
        return (lower + upper) / 2, lower - upper, beta, 0.0
    phi = (cmath.phase(a00) + cmath.phase(a11)) / 2
    alpha = cmath.phase(a11) + cmath.phase(a10) - 2 * phi
    gamma = cmath.phase(a11) - cmath.phase(a10)
    return phi, alpha, beta, gamma


def _emit_controlled_two_by_two(
    builder: Builder, ref: Ref, bit: int, anchor: int, matrix: Sequence[Sequence[complex]]
) -> None:
    """Apply ``matrix`` on the subspace where the remaining bits equal ``anchor`` (basis order bit=0,1)."""
    phi, alpha, beta, gamma = _zyz(matrix)
    controls = fuse(ref[:bit], ref[bit + 1 :])
    value = (anchor & ((1 << bit) - 1)) | ((anchor >> (bit + 1)) << bit)
    with builder.control(controls, value) if controls.width else contextlib.nullcontext():
        if abs(gamma) > _SYNTHESIS_TOLERANCE:
            builder.rz(ref[bit], gamma)
        if abs(beta) > _SYNTHESIS_TOLERANCE:
            builder.ry(ref[bit], beta)
        if abs(alpha) > _SYNTHESIS_TOLERANCE:
            builder.rz(ref[bit], alpha)
        if abs(phi) > _SYNTHESIS_TOLERANCE:
            builder.global_phase(phi)


def _emit_two_level(
    builder: Builder, ref: Ref, p: int, q: int, matrix: Sequence[Sequence[complex]]
) -> None:
    """Apply a two-level unitary acting only on ``span{|p>, |q>}``, with basis order (``|p>``, ``|q>``)."""
    bit = ((p ^ q) & -(p ^ q)).bit_length() - 1
    moved = p ^ (1 << bit)
    if moved != q:
        _emit_transposition(builder, ref, q, moved)
    if (p >> bit) & 1:
        swapped = ((matrix[1][1], matrix[1][0]), (matrix[0][1], matrix[0][0]))
        _emit_controlled_two_by_two(builder, ref, bit, p, swapped)
    else:
        _emit_controlled_two_by_two(builder, ref, bit, p, matrix)
    if moved != q:
        _emit_transposition(builder, ref, q, moved)


def _emit_unitary(builder: Builder, ref: Ref, matrix: Sequence[Sequence[complex]]) -> None:
    """Synthesize an explicit small unitary by two-level decomposition: eliminate column by column down to diagonal phases, then replay in reverse order."""
    d = len(matrix)
    u = [list(row) for row in matrix]
    steps: list[tuple[int, int, tuple[tuple[complex, complex], tuple[complex, complex]]]] = []
    for column in range(d - 1):
        for row in range(d - 1, column, -1):
            value = u[row][column]
            if abs(value) <= _SYNTHESIS_TOLERANCE:
                continue
            pivot = row - 1
            head = u[pivot][column]
            norm = math.sqrt(abs(head) ** 2 + abs(value) ** 2)
            block = (
                (head.conjugate() / norm, value.conjugate() / norm),
                (-value / norm, head / norm),
            )
            steps.append((pivot, row, block))
            for k in range(d):
                top, bottom = u[pivot][k], u[row][k]
                u[pivot][k] = block[0][0] * top + block[0][1] * bottom
                u[row][k] = block[1][0] * top + block[1][1] * bottom
    for basis in range(d):
        phase = cmath.phase(u[basis][basis])
        if abs(phase) > _SYNTHESIS_TOLERANCE:
            with builder.control(ref, basis):
                builder.global_phase(phase)
    for pivot, row, block in reversed(steps):
        inverse = (
            (block[0][0].conjugate(), block[1][0].conjugate()),
            (block[0][1].conjugate(), block[1][1].conjugate()),
        )
        _emit_two_level(builder, ref, pivot, row, inverse)


def _rotation_operation(matrix: Sequence[Sequence[complex]]) -> Operation:
    """Synthesize an explicit small unitary matrix into a quantum operation via two-level decomposition."""
    n = (len(matrix) - 1).bit_length()
    b = Builder(_name("df_rotation", matrix), {"target": Bits(n)})
    _emit_unitary(b, b["target"], matrix)
    return annotate(b.finish(), "unitary", algorithm="two_level_synthesis")


def _diagonal_encoding(spectrum: Sequence[float]) -> BlockEncoding:
    """Controlled-rotation block encoding of a diagonal matrix: single-bit signal, ``cos(θ_t/2) = g_t/α``, ``α = Σ_p |g_p|``.

    A controlled ``Ry(θ_t)`` is applied for each basis state ``|t>`` so the (0,0)
    block is exactly ``diag(g)/α``; ``α`` is the 1-norm of the spectrum, matching
    the reporting convention of the outer DF/THC assembly.
    """
    d = len(spectrum)
    n = (d - 1).bit_length()
    alpha = sum(abs(g) for g in spectrum)
    if not alpha:
        raise ValidationError("diagonal spectrum cannot be all zero")
    b = Builder(
        _name("df_diagonal", spectrum),
        {"target": Bits(n), "signal": Bits(1)},
    )
    for t, g in enumerate(spectrum):
        theta = 2 * math.acos(max(-1.0, min(1.0, g / alpha)))
        if theta:
            with b.control(b["target"], t):
                b.ry(b["signal"], theta)
    return BlockEncoding(
        annotate(b.finish(), "block_encoding", be_alpha=alpha, be_form="diagonal_controlled_rotation")
    )


def _conjugated_encoding(
    rotation: Sequence[Sequence[complex]], spectrum: Sequence[float]
) -> BlockEncoding:
    """Block encoding of ``U diag(g) U†``: apply the two-level-synthesized ``U`` conjugately on both sides of the diagonal block encoding."""
    n = (len(rotation) - 1).bit_length()
    rotation_op = _rotation_operation(rotation)
    diagonal = _diagonal_encoding(spectrum)
    b = Builder(
        _name("df_term", rotation, spectrum),
        {"target": Bits(n), "signal": Bits(diagonal.signal_qubits)},
        resources_for(("rot", rotation_op), ("diag", diagonal.operation)),
    )
    with b.adjoint():
        invoke(b, rotation_op, "rot", target=b["target"])
    invoke(b, diagonal.operation, "diag", target=b["target"], signal=b["signal"])
    invoke(b, rotation_op, "rot", target=b["target"])
    return BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=diagonal.alpha))


def double_factorized_encoding(
    df: DoubleFactorization, *, name: str | None = None
) -> BlockEncoding:
    """LCU block encoding of a DF Hamiltonian: the outer PREPARE runs over the rank index ``r`` with weights ``∝ ‖g_r‖₁``.

    SELECT applies the term block encoding of ``U_r diag(g_r)/‖g_r‖₁ U_r†`` under
    control, where the diagonal ``G_r`` evolution uses an explicit PREPARE–SELECT
    and ``U_r`` is synthesized by two-level decomposition; the ``scalar`` term
    joins the same outer LCU. alpha is ``|scalar| + Σ_r ‖g_r‖₁`` and is reported
    in the ``df_lambda`` attribute. The result can be fed directly to
    transforms.qubitization_walk.

    Args:
        df: DF input model instance.
        name: Unused by the current implementation; kept only for interface compatibility.

    Returns:
        BlockEncoding: LCU block encoding of ``H = scalar·I + Σ_r U_r diag(g_r) U_r†``.
    """
    require_instance(df, DoubleFactorization, "double_factorized_encoding.df")
    del name
    n = df.width
    terms: list[tuple[float, BlockEncoding]] = []
    if df.scalar:
        terms.append((df.scalar, identity(n)))
    for rotation, spectrum in zip(df.rotations, df.spectra, strict=True):
        if sum(abs(g) for g in spectrum):
            terms.append((1.0, _conjugated_encoding(rotation, spectrum)))
    if not terms:
        raise ValidationError("all coefficients of the DF Hamiltonian are zero")
    out = lcu(terms)
    return BlockEncoding(
        annotate(
            out.operation,
            "block_encoding",
            be_alpha=out.alpha,
            be_form="double_factorization",
            df_rank=df.rank,
            df_lambda=out.alpha,
            lcu_terms=len(terms),
        )
    )


def thc_encoding(thc: THCDecomposition) -> BlockEncoding:
    """LCU block encoding of a THC Hamiltonian: the outer PREPARE runs over ``(μ, ν)`` pairs with weights ``∝ |ζ_{μν}|·α_μ α_ν``.

    Each ``(μ, ν)`` term is a block-encoding product of ``L_μ L_ν†`` (leaf
    operators encoded via matrix_pauli_encoding, with ``α_μ`` their Pauli l1
    bound). alpha is ``Σ_{μν} |ζ_{μν}| α_μ α_ν`` and is reported in the
    ``thc_lambda`` attribute. The result can be fed directly to
    transforms.qubitization_walk.

    Args:
        thc: THC input model instance.

    Returns:
        BlockEncoding: LCU block encoding of ``Σ_{μν} ζ_{μν} L_μ L_ν†``.
    """
    require_instance(thc, THCDecomposition, "thc_encoding.thc")
    encodings = [matrix_pauli_encoding(leaf) for leaf in thc.leaves]
    adjoints = [adjoint_be(be) for be in encodings]
    terms: list[tuple[float, BlockEncoding]] = []
    for mu, row in enumerate(thc.coefficients):
        for nu, zeta in enumerate(row):
            if zeta:
                terms.append((zeta, product(encodings[mu], adjoints[nu])))
    if not terms:
        raise ValidationError("all coefficients of the THC Hamiltonian are zero")
    out = lcu(terms)
    return BlockEncoding(
        annotate(
            out.operation,
            "block_encoding",
            be_alpha=out.alpha,
            be_form="thc",
            thc_leaves=thc.leaf_count,
            thc_lambda=out.alpha,
            lcu_terms=len(terms),
        )
    )
