# General QHAM: derivation from PDEs to quantum-adapted linear systems

**English** · [简体中文](../zh/reference/qham-derivation.html)

This document first fixes the mathematical rules and then implements the generator accordingly. The primary reference is [QHAM v2](https://arxiv.org/html/2411.06759v2) II.1–II.3; its [v1](https://arxiv.org/html/2411.06759v1) uses the name secondary linearization. "Secondary" here means a second linearization, not first reducing every PDE to a quadratic polynomial. For the engineering description see the manual chapter [General QHAM automatic generation: PDE → HAM → QCL → QODE](../manual/qham.md); for a runnable example see the tutorial [generating QHAM inputs from PDE expressions](../tutorials/qham.md).

## 1. Input scope and conventions

Consider first-order time-evolution systems, written uniformly after spatial discretization as

```text
u' = f + L u + sum_{r=2}^D B_r(u,...,u),    u(0)=u_in.
```

u may contain multiple components and multidimensional grids. Each B_r is a known r-linear map V^tensor(r) → V, for example -u*u_x in Burgers, -6*u*u_x in KdV, and cubic reaction terms. Spatial derivatives, known spatial coefficients, component selection, and same-point products are folded into these maps, not treated as nonlinear gates on the quantum state.

The first batch of automatic representation rules covers autonomous, finite-order polynomial evolution PDEs, with known spatial forcing and coefficients allowed. Unknown field denominators, arbitrary transcendental nonlinearity, implicit time equations, and uneliminated constraints are outside the rules. The polynomial degree, HAM truncation order, component count, and spatial dimension within the rules are not pinned to any particular case. Time-dependent L, f, B_r, and eta can follow the algebra below, but they need time-dependent QODE inputs; freezing them directly and calling it the same problem is not allowed.

Choose L_H = partial_t - L, U0'=L U0+f, U0(0)=u_in, Ui(0)=0 (i>0). Let eta=hH, taken as a constant in the current quantum assembly.

## 2. Deriving the higher-order recursion from the homotopy equation

Starting from (1-q)L_H(U-U0)-q*eta*(partial_t U-N(U))=0, write U(q)=sum_i q^i Ui and define

```text
C_l = sum_{r=2}^D sum_{a1+...+ar=l} B_r(U_a1,...,U_ar),
D_i = (partial_t-L) Ui.
```

Comparing coefficients of q gives D_1=-eta*C_0 and D_i=(1+eta)D_(i-1)-eta*C_(i-1) (i>1). Hence

```text
Ui' = L Ui - eta * sum_{l=0}^{i-1} (1+eta)^(i-1-l) C_l.       (A)
```

This identity is the basis of the generator. In particular, for m=2:

```text
U0' = L U0 + f
U1' = L U1 - eta B_2(U0,U0)
U2' = L U2 - eta [B_2(U0,U1)+B_2(U1,U0)]
              - eta(1+eta) B_2(U0,U0)
```

The last term vanishes exactly at eta=-1; testing only with that parameter would miss errors at general h. The -eta in the paper's local-coupling schematic must not simply be copied onto every higher-order edge.

Forcing appears explicitly only in the U0 equation; the other Ui are affected by the forcing indirectly through U0. This depends on the choice of initial guess above.

## 3. Tensor variables on independent coordinates

For an ordered word of non-negative integers a=(a0,...,a_(k-1)), define

```text
Y_a(x0,...,x_(k-1)) = product_j U_aj(xj).
```

These factors have independent spatial coordinates. Y_(0,1) and Y_(1,0) must not be merged directly without a coordinate permutation. Same-point multiplication must apply the spatial operators to the corresponding factors first and only then perform the diagonal contraction. In the discrete case, B_r is exactly such a rectangular linear map. In the implementation, {obj}`structured_fd_bindings <oracq.applications.qham.stencils.structured_fd_bindings>` generates the elementary matrices from shifts and contractions without materializing port matrices.

Apply the product rule once for the time derivative:

```text
Y_a' = sum_p I ⊗ ... ⊗ (U_ap') ⊗ ... ⊗ I.
```

Substituting into (A) yields three classes of linear edges:

1. Linear terms: same word to same word, applying L on the p-th coordinate.
2. Forcing terms with a_p=0: starting from the word with that factor deleted, insert the known vector f at the p-th position.
3. Nonlinear terms with a_p>0: replace that factor by an ordered r-tuple b with sum(b)=l<a_p; B_r acts on these r input coordinates and contracts to one output coordinate with coefficient -eta(1+eta)^(a_p-1-l).

Consequently, the quantum state of U_i need not serve as an input that another solver copies repeatedly. The whole system needs only access to the initial data and to fixed linear operators/multilinear maps.

## 4. The finite closure proof

Let p=max(1,D-1) and define the weight of a word as

```text
weight(a)=p*sum(a)+len(a),    W=p*m+1.
```

The initial single-factor Ui (i≤m) all satisfy weight≤W. After a nonlinear replacement

```text
weight(new)-weight(old)
 = p*(l-a_p)+(r-1) <= -p+(D-1) <= 0.
```

Moreover the order sum strictly decreases. A forcing edge deletes a zero-order factor and reduces the weight by 1; a linear edge leaves it unchanged. All edges therefore stay inside the finite set weight(a)≤W, with at most W dynamical factors needed. The all-zero word produces only linear terms and forcing step-down terms, so the closure terminates.

This closure is exact for the given HAM-truncated system; it is not another Carleman order truncation. The generator-side {obj}`QHAMPlan <oracq.applications.qham.linearization.QHAMPlan>` lazily enumerates the closure under the same weight criterion and does not materialize matrices. The differences from the original nonlinear PDE still come from the HAM truncation, the spatial discretization, and the subsequent numerical/quantum solving.

For a quadratic PDE, p=1 and the word count is 2^(m+1)-1. Adding the physical output block gives the paper's 2^(m+1) function blocks. When each single-coordinate space has dimension N, the original total dimension is

```text
N + sum_{k=1}^{m+1} binom(m+1,k) N^k
  = (N+1)^(m+1)+N-1.
```

Higher-degree cases use the canonical superset of the finite closure above; no claim is made that every block is necessary or that the paper's optimal resource scale is reached.

## 5. The physical solution and forcing homogenization

Additionally introduce the physical block Y_phys=sum_{i=0}^m Ui. Summing (A) gives

```text
Y_phys' = L Y_phys + f
        + sum_{l=0}^{m-1} [1-(1+eta)^(m-l)] C_l.              (B)
```

The relation also holds at eta=0; no division by eta is needed. The initial value is Y_phys(0)=u_in. Among the tensor words only the all-zero word has a nonzero initial value: Y_(0,...,0)(0)=u_in^tensor(k); all other initial values are zero.

With forcing, introduce a constant component c'=0, c(0)=1 and write the f injection from the empty word into the linear generator as well, obtaining a homogeneous augmented QODE. The action of the forcing on the tensor blocks is a rectangular injection from lower-rank words to higher-rank words; one must not merely add a right-hand-side vector to U0 and omit the forcing terms of the tensor equations.

## 6. Validation method and boundaries

Mathematical validation does not depend on QLSS outputs. For arbitrary random Ui, compute separately:

```text
chain_rule = derivative_of_Lift(U) applied to HAM_rhs(U)
linear_rhs = G * Lift(U).
```

The two should agree component by component. The tests must cover m>1, eta other than -1, non-zero f, cubic terms, mixed polynomials, and coupled components. With m=1 and no forcing, the existing [u,U0,U1,U0⊗U0] model should also be recovered.

The closure proof does not prove HAM convergence, nor does it guarantee that an arbitrary generator satisfies the dissipativity assumptions of LCHS/CBMD. An implementation that explicitly generates all blocks has combinatorial size; lazier row-by-row rules can later be used to design more efficient sparse oracles. Keeping modules in MIR/RIR does not automatically confer the paper's complexity advantages.
