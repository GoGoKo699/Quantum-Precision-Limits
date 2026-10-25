# Question, model, and guarantees

How accurately can coherent access to a quantum circuit determine a mixed-state
precision limit, when the physical state is too large to reconstruct classically?

The project connects two tasks. The first is a general tomography problem:
recover a small reduced state in Bures distance even when its purification lives
in a much larger workspace. The second uses that capability to evaluate the
Holevo bound of a low-rank circuit model. The output of the second task is a
classical certified interval for a local precision benchmark.

## Reading route

| Question | Authoritative account |
|---|---|
| How does coherent tomography achieve the Bures guarantee? | [Bures tomography](bures_tomography.md) |
| Why can a large sensing model be represented by a small Gram matrix? | [Compression and stability](holevo_compression_and_stability.md) |
| How are physical fitting and optimization certified? | [Certified optimization](certified_optimization.md) |
| What must the circuit input provide, and what does the algorithm cost? | [Circuit access and algorithm](circuit_access.md) |
| Which examples, baselines, and noise assumptions matter? | [Applications and limits](applications_and_limits.md) |
| Why is state-weight conditioning substantive? | [Conditioning obstruction](conditioning_limits.md) |
| Which ingredients are inherited? | [Related work](../references/RELATED_WORK.md) |

## Two different dimensions

Use $`d`$ for the dimension of the state to be learned by the tomography primitive.
Use $`D`$ for the physical dimension of a sensing model. These need not be equal.
The application will learn a label state of dimension $`m+1`$, where
$`m=r(p+1)`$, while the $`D`$-dimensional sensor becomes part of its purification.

| Symbol | Meaning |
|---|---|
| $`r`$ | Rank bound in tomography; occupied state rank in the sensing application |
| $`p`$ | Number of real parameters being estimated |
| $`\delta`$ | Failure probability of the acquisition procedure |
| $`\varepsilon_B`$ | Requested Bures-distance accuracy |
| $`\varepsilon`$ | Requested absolute error in the Holevo value |
| $`\nu`$ | Promised lower bound on every occupied state eigenvalue |
| $`j`$ | Promised lower bound on the smallest SLD Fisher eigenvalue |
| $`L`$ | Bound on the combined derivative-factor norm |
| $`\alpha`$ | Known normalization of the coherent factor/derivative preparation |

The notation $`\widetilde O`$ suppresses logarithms in the dimensions, inverse
accuracy, and inverse failure probability. Polynomial dependence on
$`\alpha,1/\nu,1/j,L,p,r`$ is never hidden in that notation.

## Distances and their operational role

For normalized positive semidefinite matrices, define fidelity amplitude and
Bures distance by

```math
f(\rho,\sigma)=\|\sqrt\rho\sqrt\sigma\|_1,
\qquad
D_B(\rho,\sigma)=\sqrt{2-2f(\rho,\sigma)}.
```

The half trace distance obeys $`\tfrac12\|\rho-\sigma\|_1\le D_B(\rho,\sigma)`$.
For arbitrary positive semidefinite matrices, use the homogeneous extension

```math
d_B(K,H)^2=\mathrm{Tr}K+\mathrm{Tr}H
 -2\|\sqrt K\sqrt H\|_1.
```

This equals the minimum squared Frobenius distance between compatible factors,
after aligning their row spaces by an isometry. It therefore controls the
physical perturbation of state and derivative factors even when the full Gram
matrix has zero eigenvalues. These metric facts are standard; the
[tomography proof](bures_tomography.md) identifies the particular weighted
quadratic estimate used here.

## General tomography theorem

Let a unitary $`U`$ and its inverse prepare a purification of a
$`d`$-dimensional state $`\rho`$ of rank at most $`r`$. The purifying register can
be arbitrarily large. For $`d\ge2`$ and $`0<\varepsilon_B\le1/4`$, the algorithm
outputs a classical density matrix $`\widehat\rho`$ satisfying

```math
\Pr\{D_B(\rho,\widehat\rho)\le\varepsilon_B\}\ge1-\delta,
\qquad
Q=\widetilde O\!\left(\frac{dr}{\varepsilon_B}\right).
```

The query count includes every use of $`U`$ and $`U^\dagger`$ inside conditioning,
amplitude estimation, and the invoked Frobenius tomography. The preparation
gate cost and workspace width remain charged. No positive eigenvalue floor,
Schmidt basis, or minimal purification is supplied.

The matching worst-case query lower bound, for sufficiently small
constant-bounded accuracy, follows from the published trace-norm tomography
lower bound of van Apeldoorn, Cornelissen, Gilyén, and Nannicini. The upper proof
combines their Frobenius tomography with the adaptive spectral geometry of
Flammia and O'Donnell. See the [proof and attribution](bures_tomography.md) for
the exact relationship. With a purification of dimension exactly $`r`$,
pure-state tomography already gives the same scaling; the arbitrary-purifier
guarantee is the additional capability studied here.

## Local sensing model

At an operating point, let $`\Psi,\Phi_i\in\mathbb C^{D\times r}`$ satisfy

```math
\rho=\Psi\Psi^\dagger,\qquad
D_i=\Phi_i\Psi^\dagger+\Psi\Phi_i^\dagger,
\qquad
\|\Psi\|_F=1,\quad
\mathrm{Re}\mathrm{Tr}(\Psi^\dagger\Phi_i)=0.
```

Here $`D_i`$ denotes the derivative of the density matrix, not its dimension.
The factor $`\Psi`$ has exactly $`r`$ independent columns. Define

```math
V=[\Psi,\Phi_1,\ldots,\Phi_p],\qquad
K=V^\dagger V,\qquad m=r(p+1),\qquad
\sum_i\|\Phi_i\|_F^2\le L^2.
```

The occupied state weights are at least $`\nu>0`$. The symmetric logarithmic
derivatives $`L_i^{\rm SLD}`$ solve
$`D_i=(\rho L_i^{\rm SLD}+L_i^{\rm SLD}\rho)/2`$. Their real covariance defines
the SLD Fisher matrix

```math
J_{ij}=\mathrm{Re}\mathrm{Tr}
 (\rho L_i^{\rm SLD}L_j^{\rm SLD}),\qquad J\succeq jI_p,
\quad j>0.
```

Kernel-kernel SLD blocks may be set to zero. No lower bound is imposed on the
positive eigenvalues of the **whole** tangent Gram matrix $`K`$. The state-weight
promise concerns $`K_{00}=\Psi^\dagger\Psi`$ and is a different assumption.

The local data describe a valid smooth fixed-rank model: near the origin one can
normalize $`\Psi+\sum_i\theta_i\Phi_i`$ and multiply it by its adjoint.
Unitary sensing is an important special case,
$`\Phi_i=-iG_i\Psi`$, with known bounded Hermitian generators.

## The quantity computed

For equal parameter weights, the Holevo value is

```math
C_H(I_p)=\min_{X_1,\ldots,X_p}
 \left\{\mathrm{Tr}(\mathrm{Re}Z)
       +\|\mathrm{Im}Z\|_1\right\},
\qquad Z_{ij}=\mathrm{Tr}(\rho X_iX_j),
```

where $`X_i=X_i^\dagger`$,
$`\mathrm{Tr}(\rho X_i)=0`$, and
$`\mathrm{Tr}(D_jX_i)=\delta_{ij}`$. Matrix real and imaginary parts are
entrywise. The imaginary part is a real antisymmetric matrix; its trace norm
accounts for incompatibility between the parameters.

Under the usual regularity conditions of finite-dimensional local asymptotic
quantum estimation, this is the asymptotically attainable coefficient of the
weighted mean-square error with collective measurements. The algorithm
evaluates the local optimization. The algorithmic promises here do not by
themselves establish an attainment theorem for every rank-deficient model;
boundary models require a separate applicable result. The selected
[estimation bridge](../learning/ESTIMATION_BRIDGE.md) distinguishes these
conditions from the review's full-rank QLAN setting.

The algorithm does not output a collective measurement,
a finite-sample estimator, or a guarantee for unrestricted laboratory noise.
The [related-work account](../references/RELATED_WORK.md) gives the metrology
sources and the distinction from the SLD bound $`\mathrm{Tr}J^{-1}`$.

For a real positive definite weight $`W`$, replace the objective by
$`\mathrm{Tr}(W\mathrm{Re}Z)+
\|\sqrt W\mathrm{Im}Z\sqrt W\|_1`$.
The weighted stability theorem retains the extreme eigenvalues of $`W`$.
The headline algorithm uses $`W=I_p`$ to keep these extra conditioning factors
visible rather than absorb them into a constant.

## Scope of the algorithmic claim

The input includes coherent preparation of all factor and derivative columns
with known normalization $`\alpha`$, and an inverse of that preparation. The
[access construction](circuit_access.md) specifies how a short, hardcoded
linear combination of known generator unitaries supplies it without QRAM.
The small label state can have rank $`m+1`$, even though the physical state has
rank $`r`$.

At fixed $`p,r,\nu,j,L,\alpha`$, acquisition has
$`\widetilde O(1/\varepsilon)`$ query dependence. The complete resource formula
also includes the conditioning threshold, preparation gates, small-register
operations, and classical fitting/optimization. A matching precision lower
bound holds in the stated opaque preparation-oracle model. A known toy circuit
with an elementary closed form can of course be evaluated classically.

The exact finite checks in [verification](../verification/README.md) inspect
specific algebraic and numerical consequences. General optimality and
stability are supplied by the proofs and their cited dependencies.
