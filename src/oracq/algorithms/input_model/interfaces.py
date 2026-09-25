"Python structural protocols shared by the algorithm library; applications may keep defining protocols in their own modules."

from __future__ import annotations

from inspect import signature
from typing import TYPE_CHECKING, Protocol, TypeVar, runtime_checkable

from oracq.algorithms.input_model.contracts import (
    AlgorithmContract,
    ContractError,
    InputRequirement,
    fail,
    require_instance,
)

if TYPE_CHECKING:
    from oracq.algorithms.input_model.operators import BlockEncoding
    from oracq.algorithms.input_model.oracles import (
        SparseAccess,
        StateOracle,
        StatePreparation,
        XorDatabase,
    )
    from oracq.infrastructure.builder import Operation


@runtime_checkable
class StatePreparationProtocol(Protocol):
    """Requires a ``state_preparation()`` method; anything satisfying it can serve as an algorithm's initial-state input."""

    def state_preparation(self) -> StatePreparation:
        """Return the ``StatePreparation`` access view corresponding to this input.

        Returns:
            StatePreparation: View that prepares the target state on this input's zero-input registers.
        """
        ...


@runtime_checkable
class UnitaryProtocol(StatePreparationProtocol, Protocol):
    """Adds the full unitary interface ``unitary()`` on top of the state-preparation protocol."""

    def unitary(self) -> Operation:
        """Return the unitary ``Operation`` on the full register space.

        Returns:
            Operation: Unitary operation acting on all registers of this input.
        """
        ...


@runtime_checkable
class BlockEncodingProtocol(Protocol):
    """Requires a ``block_encoding()`` method; anything satisfying it can serve as a block-encoding input."""

    def block_encoding(self) -> BlockEncoding:
        """Return the ``BlockEncoding`` access view corresponding to this input.

        Returns:
            BlockEncoding: Access view encoding this input as a normalized sub-block of a unitary.
        """
        ...


@runtime_checkable
class CKSSparseProtocol(Protocol):
    """Requires a ``sparse_access()`` method; anything satisfying it can serve as a CKS sparse-matrix input."""

    def sparse_access(self) -> SparseAccess:
        """Return the ``SparseAccess`` composed of the CKS location and entry operations.

        Returns:
            SparseAccess: Sparse access view composed of a location oracle and an entry oracle.
        """
        ...


@runtime_checkable
class XorDatabaseProtocol(Protocol):
    """Requires a ``xor_database()`` method; anything satisfying it can serve as an XOR database input."""

    def xor_database(self) -> XorDatabase:
        """Return the ``XorDatabase`` access view corresponding to this input.

        Returns:
            XorDatabase: Database access view that answers queries by XOR.
        """
        ...


@runtime_checkable
class StateOracleProtocol(Protocol):
    """Requires a ``state_oracle()`` method; anything satisfying it can serve as an output-state oracle."""

    def state_oracle(self) -> StateOracle:
        """Return the ``StateOracle`` carrying a success signal.

        Returns:
            StateOracle: Oracle that marks the success subspace within the output state.
        """
        ...


_View = TypeVar("_View")


def _view(value: object, method_name: str, cls: type[_View]) -> _View:
    """Obtain the access view exactly once; call-boundary errors and provider-internal exceptions are handled separately."""
    if isinstance(value, cls):
        return value
    method = getattr(value, method_name, None)
    if method is None:
        fail("INPUT_PROTOCOL", method_name, "provider method", type(value).__name__, "missing access method")
    if not callable(method):
        fail("INPUT_CALLABLE", method_name, "callable", type(method).__name__, "access attribute is not callable")
    try:
        parameters = signature(method)
    except (TypeError, ValueError):
        # Native callables may not expose a signature; exceptions from the actual call keep the original traceback.
        parameters = None
    if parameters is not None:
        try:
            parameters.bind()
        except TypeError as exc:
            from oracq.algorithms.input_model.contracts import ContractIssue

            raise ContractError((ContractIssue(
                "INPUT_SIGNATURE", method_name, "zero-argument call", str(parameters),
                "access method must support a zero-argument call",
            ),)) from exc
    return require_instance(method(), cls, method_name + ".output")


def as_state_preparation(value: StatePreparationProtocol) -> StatePreparation:
    """Convert a protocol-satisfying input into a ``StatePreparation``.

    Args:
        value: Object satisfying ``StatePreparationProtocol``.

    Returns:
        StatePreparation: The concrete view obtained by calling ``state_preparation()``.

    Raises:
        ContractError: Reported as ``INPUT_PROTOCOL`` when the interface is missing;
            as ``INPUT_TYPE`` when the return value is not a ``StatePreparation``.
    """
    from oracq.algorithms.input_model.oracles import StatePreparation

    return _view(value, "state_preparation", StatePreparation)


def checked_state_preparation(
    value: StatePreparationProtocol,
    *,
    adjoint: bool = False,
    controlled: bool = False,
    path: str = "preparation",
) -> StatePreparation:
    """Obtain a zero-input, clean-work preparation; the calling algorithm selects the transformation capabilities it needs.

    Args:
        value: Input object satisfying ``StatePreparationProtocol``.
        adjoint: Whether the preparation must offer a Hermitian-adjoint version.
        controlled: Whether the preparation must support a controlled version.
        path: Path name of this input in the contract report.

    Returns:
        StatePreparation: Preparation view that passed the zero-input and clean-work checks.
    """
    result = as_state_preparation(value)
    requirement = InputRequirement(
        path,
        (StatePreparationProtocol,),
        zero_input=True,
        clean_work=True,
        adjoint=adjoint,
        controlled=controlled,
    )
    _, issues = requirement.inspect(result)
    if issues:
        raise ContractError(issues)
    return result


def as_block_encoding(value: BlockEncodingProtocol) -> BlockEncoding:
    """Convert a protocol-satisfying input into a ``BlockEncoding``.

    Args:
        value: Object satisfying ``BlockEncodingProtocol``.

    Returns:
        BlockEncoding: The concrete view obtained by calling ``block_encoding()``.

    Raises:
        ContractError: Reported as ``INPUT_PROTOCOL`` when the interface is missing;
            as ``INPUT_TYPE`` when the return value is not a ``BlockEncoding``.
    """
    from oracq.algorithms.input_model.operators import BlockEncoding

    return _view(value, "block_encoding", BlockEncoding)


def as_sparse_access(value: CKSSparseProtocol) -> SparseAccess:
    """Convert a protocol-satisfying input into a ``SparseAccess``.

    Args:
        value: Object satisfying ``CKSSparseProtocol``.

    Returns:
        SparseAccess: The concrete view obtained by calling ``sparse_access()``.

    Raises:
        ContractError: Reported as ``INPUT_PROTOCOL`` when the interface is missing;
            as ``INPUT_TYPE`` when the return value is not a ``SparseAccess``.
    """
    from oracq.algorithms.input_model.oracles import SparseAccess

    return _view(value, "sparse_access", SparseAccess)


def as_qlss_matrix(
    value: BlockEncodingProtocol | CKSSparseProtocol,
) -> BlockEncoding | SparseAccess:
    """Convert a matrix input into a ``BlockEncoding`` or ``SparseAccess`` according to the available interface.

    When ``BlockEncodingProtocol`` is satisfied the block-encoding view is
    returned preferentially; otherwise the sparse access view is taken per
    ``CKSSparseProtocol``.

    Args:
        value: Matrix input object.

    Returns:
        A ``BlockEncoding`` when the block-encoding protocol is satisfied, otherwise a ``SparseAccess``.

    Raises:
        ContractError: Neither protocol is satisfied, or the obtained return value has the wrong type.
    """
    return (
        as_block_encoding(value)
        if isinstance(value, BlockEncodingProtocol)
        else as_sparse_access(value)
    )


def operator_state_contract(
    name: str,
    *,
    matrix: str = "generator",
    state: str = "initial",
    composable: bool = True,
) -> AlgorithmContract:
    """Requirements of the current BE-type evolution algorithms; other algorithms may define entirely different contracts.

    Args:
        name: Contract name, used in reports and diagnostics.
        matrix: Path name of the block-encoding matrix input in the contract.
        state: Path name of the initial-state preparation input in the contract.
        composable: Whether the initial state requires adjoint and controlled capabilities for use by composable algorithms.

    Returns:
        AlgorithmContract: Contract requiring the matrix to be a block encoding
        with adjoint and controlled capabilities, the initial state to be a
        zero-input and clean-work preparation, and constraining both to equal
        register widths.
    """
    return AlgorithmContract(
        name,
        (
            InputRequirement(
                matrix,
                (BlockEncodingProtocol,),
                adapter=as_block_encoding,
                adjoint=True,
                controlled=True,
            ),
            InputRequirement(
                state,
                (StatePreparationProtocol,),
                adapter=as_state_preparation,
                adjoint=composable,
                controlled=composable,
                zero_input=True,
                clean_work=True,
            ),
        ),
        same_width=((matrix, state),),
    )
