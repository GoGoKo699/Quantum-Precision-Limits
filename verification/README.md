# Reproduce the finite checks

The written proofs establish the general tomography, stability, and resource
claims. This directory checks specific algebraic and numerical consequences
that could expose an implementation or derivation error. The checks cover
three distinct failure modes: accumulated Bures error, validity of classical
certificates at a singular Gram, and incompatibility in a concrete circuit.

## Run

The recorded environment uses Python 3 with NumPy 2.3.5, SciPy 1.17.0, and
SymPy 1.14.0. Install the [pinned dependencies](requirements.txt) in an environment
of your choice, then run from the repository root:

```bash
python3 -m pip install -r verification/requirements.txt
python3 verification/run_all.py
```

The runner executes all three checks, checks their assertions, and compares
them with the committed reference reports. Newly generated reports go to the
ignored `verification/actual/` directory. It does not replace the reference
reports. Floating-point comparisons allow small platform variation; exact
rational certificate assertions use no numerical tolerance.

| Check | What is verified | Arithmetic and limits |
|---|---|---|
| [bures_verify.py](bures_verify.py) | Residual mass removal, weighted error accounting, exact normalization identity, factor rounding, and rank preservation under phase sequences | Deterministic seeded small matrices with floating-point residual bounds; the actual quantum Frobenius learner is not simulated |
| [verify_certified_optimization.py](verify_certified_optimization.py) | Physical-fit constraints, rational LDL factors, Hermitian action constraints, primal/dual PSD witnesses, and objective equality | Exact complex rational arithmetic; two analytic models, including a singular tangent Gram; not a generic optimizer |
| [verify_holevo_circuit.py](verify_holevo_circuit.py) | Four-qubit state preparation, exact rank/weights, locally unbiased witnesses, certified value interval, Fisher floor, and numerical tangent compression | Rational primal/dual certification after floating-point proposals; one fixed circuit, with a separate numerical compression comparison |

## Bures geometry

The test replaces each stochastic quantum Frobenius estimate with a controlled
finite matrix perturbation that satisfies the residual-estimation lemma. It
includes rank-one states, mixed states, degeneracy, and eigenvalues across many
scales. Dense, diagonal, and coherent perturbations exercise different error
locations. The checks inspect the deterministic implications used by the
[Bures proof](../science/bures_tomography.md).

Additional cases apply rounded generalized phase rotations to a purification
whose environment is larger than its system. They verify the two-dimensional
invariant span and the absence of a numerical singular-value tail beyond the
input rank. These cases test the algebraic rank argument. They do not show that
arbitrary gate noise preserves rank; the proof handles compilation errors with
a separate global distribution bound.

Reference: [bures_verification.json](expected/bures_verification.json).

## Exact optimization certificates

The two analytic models have compatibility parameters $`c=1`$ and $`c=3/5`$.
The first has a tangent Gram of order 3 and rank 2. The exact upper and lower
Holevo values are respectively $`4`$ and $`20/9`$; both SLD costs are $`2`$.
Each case includes a rational physical-fit witness against a regularized
reference. At $`\gamma=1/100`$, its squared-distance objective is exactly
$`1/640000`$, below the acceptance threshold $`1/40000`$.

The checker uses recursive rational Schur complements, including exact zero
pivots, to establish positive semidefiniteness. It then verifies every original
affine equation, the reduced independent equations, dual stationarity, and the
primal–dual gap. See [the finite optimization proof](../science/certified_optimization.md)
for the general formulation and the distinction between a certifiable finite
algorithm and the implemented example generators.

Reference: [certified_optimization_checks.json](expected/certified_optimization_checks.json).

## Four-qubit circuit

The [application account](../science/applications_and_limits.md) specifies all
nine gates as $`cI-isP`$ with rational $`c,s`$. Parity gives two exactly orthogonal
probe states with weights $`3/4,1/4`$. The sensed generators are $`X_1/2,Y_1/2`$.
The program constructs an exact rational complex unitary basis completion and
works on off-diagonal Hermitian estimator coordinates in the state's eigenbasis.

Floating-point optimization proposes a scalar dual parameter, multipliers,
and estimator coordinates. The program rounds them to rationals, repairs the
unbiasedness equations exactly, and checks the primal and dual values using
exact rational arithmetic. Finding the exact numerical optimizer is unnecessary
for the validity of the resulting interval:

```math
2.250762525\le C_H(I_2)\le2.250762539.
```

The SLD cost is enclosed by $`[2.233360231,2.233360232]`$. Rational principal
minors certify $`J\succeq(4/5)I_2`$. The incompatibility penalty is at least
$`0.017402294`$. The recorded interval endpoints are rounded outward; their
displayed difference can exceed the separately rounded bound on the exact
certificate gap.

A separate floating-point calculation compresses the 16-dimensional model to
its six-dimensional tangent support. It checks the projection residual and
agreement of the computed value. The exact interval is certified on the full
model, so the floating-point compression check is not needed to validate its
endpoints. The reference witness hash records the proposed rational witnesses;
platform-dependent changes to those proposals are acceptable when the new
exact checks and reference interval comparisons pass.

Reference: [holevo_circuit.json](expected/holevo_circuit.json).

## Interpret the evidence at its stated scope

No check here implements the general quantum acquisition circuit, a generic
SDP optimizer, or fault-tolerant hardware. The four-qubit example is
classically easy. It validates a nontrivial certificate and dimension reduction;
it does not establish classical computational hardness or a hardware crossover.
The general query bounds, conditioning obstructions, and asymptotic claims
depend on the complete scientific arguments and their primary-source inputs.
