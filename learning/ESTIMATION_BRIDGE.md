# From multiparameter estimation to certified precision limits

A probe circuit is intended to sense two components of a weak field. Its designer can vary the preparation, include a specified preparation error, and compare the precision available near an operating point. The relevant question is the best joint estimation cost. Two large individual sensitivities do not settle it, because the measurements that reveal the two components can compete.

The selected educational anchor is **Rafał Demkowicz-Dobrzański, Wojciech Górecki and Mădălin Guţă, _Multi-parameter estimation beyond Quantum Fisher Information_**, Journal of Physics A 53, 363001 (2020). All section and equation references below use [arXiv:2001.11742v2, 1 September 2020](https://arxiv.org/html/2001.11742v2). The [open PDF](https://arxiv.org/pdf/2001.11742v2) provides the same version.

The review supplies the estimation framework. This bridge connects that framework to the project's factor representation, coherent acquisition, and finite certificates. The [30-milestone map](MILESTONES.md) organizes the complete route. The [scientific model](../science/model.md) remains the authoritative statement of the theorem assumptions.

## The target and the two kinds of estimation

There are two statistical tasks in this project. Keeping them separate prevents a misleading interpretation of the algorithm's precision scaling.

| Task | Input and unknown | Output and resource |
|---|---|---|
| Sensing | A statistical model is specified; the encoded parameter $`\theta`$ is unknown. | An estimate of $`\theta`$, obtained from $`N_{\mathrm{sense}}`$ independently prepared probes, possibly with collective measurements. |
| Computing a precision benchmark | The operating point and preparation/derivative circuits are specified; their many-body amplitudes need not be known classically. | An interval for the local Holevo value $`C_H(W)`$, using $`Q`$ circuit and inverse-circuit calls. |

For a locally unbiased sensing strategy at the operating point, the bound is


```math
\mathrm{Tr}(W\Sigma_{N_{\mathrm{sense}}})
\ge \frac{C_H(W)}{N_{\mathrm{sense}}}.
```


Interpreting this coefficient as attainable requires an applicable statistical attainment theorem. The computational task estimates the coefficient to additive error $`\varepsilon`$. Its inverse-$`\varepsilon`$ query law is not a claim that a physical sensor's mean-square error scales as $`N_{\mathrm{sense}}^{-2}`$.

The project evaluates a specified model at a specified operating point. Searching over probe circuits, locating an unknown operating point, constructing an attaining measurement, and characterizing an uncontrolled environment are additional tasks. A comparison between two candidate circuits is meaningful only after their parameter units, interrogation time, and weight matrix are fixed.

## A focused route through the anchor

These are targeted readings within the single selected source.

| Connection | Anchor location | Project continuation |
|---|---|---|
| Local estimation and weighted risk | §1.2, Eqs. (1)–(5); §2.1, Eqs. (7)–(10) | [M02](MILESTONES.md#m02-local-unbiasedness-and-estimation-risk) and the calibration calculation below |
| SLD information and the Holevo optimization | §§2.2–2.5, especially Eqs. (12), (27), and (30) | [Model](../science/model.md#the-quantity-computed) and [certified optimization](../science/certified_optimization.md) |
| Compatibility | §2.6, Eq. (40); §2.8, Eq. (44) | The complete qubit example below |
| Copies and attainment | §§2.10–2.11, Eqs. (55) and (61) | The distinction between value evaluation and an attainable sensing protocol |
| A small exact example | §3.1.1, Eqs. (67)–(75) | The rotated-coordinate calculation below |
| Asymptotic interpretation and its domain | §4.5; §§4.6.1–4.6.3 | The qualifications on applying QLAN to a particular model |

The general QLAN presentation in §4.6 begins with distinct, strictly positive eigenvalues. Section 4.6.3 expressly separates rank-deficient boundary models from that statement. The project's low-rank computational promises do not establish QLAN attainability. Pure-state attainment has its own treatment in §2.11; other restricted models need an applicable result.

## A notation and assumption map

Several symbols change meaning between the review and the repository.

| Anchor | Project | Meaning |
|---|---|---|
| $`\mathfrak p`$ | $`p`$ | Number of real parameters |
| $`C`$ | $`W`$ | Real cost matrix in parameter space |
| $`F_Q`$ | $`J`$ | SLD quantum Fisher information matrix |
| $`\mathcal C^{\mathrm H}`$ | $`C_H(W)`$ | Holevo value |
| $`\mathcal C^{\mathrm{SLD}}`$ | $`\mathrm{Tr}(WJ^{-1})`$ | SLD scalar benchmark |
| $`V`$ in the covariance optimization | A real covariance majorant, called $`C`$ in the certificate SDP | This is not the project's factor array $`V=[\Psi,\Phi_1,\ldots,\Phi_p]`$. |
| $`J`$ in §2.12 | A parameter-change Jacobian | This is not the project's Fisher matrix $`J`$. |
| $`L_i`$ | $`L_i^{\mathrm{SLD}}`$ | SLD operators; the unindexed project constant $`L`$ instead bounds derivative-factor norms. |

The project uses $`D`$ for physical Hilbert-space dimension and $`D_i`$ for a state derivative. The subscript distinguishes them. The label dimension is at most $`d_{\mathrm{lab}}=r(p+1)+1`$. The probe rank $`r`$, label-state rank $`R`$, and number of parameters $`p`$ need not coincide.

The project's Bures convention is explicit:


```math
f(\rho,\sigma)=\|\sqrt\rho\sqrt\sigma\|_1,
\qquad
D_B(\rho,\sigma)=\sqrt{2-2f(\rho,\sigma)}.
```


For unnormalized positive matrices it uses


```math
d_B(K,H)^2=\mathrm{Tr}K+\mathrm{Tr}H
-2\|\sqrt K\sqrt H\|_1.
```


A quadratic Bures loss and Bures distance have different precision parameters. The review discusses a quadratic Bures cost in §4.6.4; the tomography theorem here bounds the square-root distance above. When fidelity means $`f^2`$, small infidelity is of order $`D_B^2`$. An inverse-infidelity sample bound therefore has inverse-square dependence on Bures-distance error.

The headline Holevo algorithm takes $`W=I_p`$. General positive definite weights are treated with their extreme eigenvalues retained in the stability bound. A singular weight introduces an estimand/nuisance-parameter problem; it is not silently covered by the same conditioning constants.

## Local unbiasedness as a calibration condition

Write the local expansion at the operating point as


```math
\rho_\theta=\rho+\sum_{j=1}^{p}\theta_jD_j+O(\|\theta\|^2),
\qquad \mathrm{Tr}D_j=0.
```


The origin is a coordinate choice. Suppose a POVM has effects $`M_s`$ and an outcome is assigned a real estimate vector $`t(s)`$. The corresponding first-moment operators are


```math
X_i=\sum_s t_i(s)M_s.
```


Then


```math
\mathbb E_\theta[t_i]
=\mathrm{Tr}(\rho X_i)
+\sum_j\theta_j\mathrm{Tr}(D_jX_i)
+O(\|\theta\|^2).
```


The requirements


```math
\mathrm{Tr}(\rho X_i)=0,
\qquad
\mathrm{Tr}(D_jX_i)=\delta_{ij}
```


say that the instrument has zero offset and unit first-order gain at the calibration point. They do not demand unbiasedness throughout parameter space. These are the constraints used throughout the [local model](../science/model.md).

The matrix


```math
Z_{ij}=\mathrm{Tr}(\rho X_iX_j)
```


is Hermitian and positive semidefinite, but can have complex entries. A measured covariance matrix is real. The Holevo optimization finds the least weighted trace of a real matrix that can dominate $`Z`$, while retaining the calibration constraints. Equivalently,


```math
C_H(W)=\min_{\{X_i\}}
\left[
\mathrm{Tr}(W\,\mathrm{Re}Z)
+\|\sqrt W\,\mathrm{Im}Z\,\sqrt W\|_1
\right].
```


Here real and imaginary parts are entrywise. The imaginary term keeps information that would disappear on taking only the real trace. The SLD benchmark retains only the first part of the optimization.

This is also why the operators $`X_i`$ in a feasible tuple must not automatically be identified with jointly measurable observables. They summarize first moments of a candidate estimation strategy. A real covariance majorant can exceed their real covariance, and the next example shows the difference exactly.

## A complete two-parameter qubit calculation

This is the pure-qubit example of anchor §3.1.1 in rotated coordinates. The calculation adds the factor and label-state description used by the project.

Let $`\sigma_x,\sigma_y,\sigma_z`$ denote the Pauli matrices. Prepare $`|0\rangle`$ and encode small dimensionless rotation angles through


```math
|\psi_\theta\rangle
=e^{-i\theta_2\sigma_y/2}e^{-i\theta_1\sigma_x/2}|0\rangle.
```


At $`\theta=0`$,


```math
\rho=\frac{I+\sigma_z}{2},
\qquad
D_1=-\frac{\sigma_y}{2},
\qquad
D_2=\frac{\sigma_x}{2}.
```


Both tangent directions are nonzero. Using rotation coordinates avoids the polar-angle singularity at the north pole of the Bloch sphere.

### Fisher information

The SLD equations are


```math
D_i=\frac{\rho L_i^{\mathrm{SLD}}
+L_i^{\mathrm{SLD}}\rho}{2}.
```


A solution with zero kernel–kernel blocks is


```math
L_1^{\mathrm{SLD}}=-\sigma_y,\qquad
L_2^{\mathrm{SLD}}=\sigma_x.
```


Their real covariance is $`J=I_2`$. For equal weights the SLD scalar benchmark is therefore $`2`$.

### The complete Holevo minimization

Write an arbitrary Hermitian qubit operator as a real combination of $`I,\sigma_x,\sigma_y,\sigma_z`$. The centering and derivative constraints force every feasible pair to have the form


```math
X_1=-\sigma_y+a|1\rangle\langle1|,
\qquad
X_2=\sigma_x+b|1\rangle\langle1|,
\qquad a,b\in\mathbb R.
```


The remaining terms annihilate the occupied state. Thus every feasible pair has the same covariance,


```math
Z=
\begin{pmatrix}
1&i\\
-i&1
\end{pmatrix}.
```


Its real trace is $`2`$. Its imaginary part has two singular values equal to one. Consequently


```math
C_H(I_2)=2+2=4.
```


This proves both the lower and upper value: all feasible tuples have that objective, and $`a=b=0`$ is feasible. The average SLD commutator is $`2i`$, so it does not satisfy the compatibility criterion.

### An actual locally unbiased measurement

Use four effects indexed by a setting $`k\in\{1,2\}`$ and a sign $`s\in\{+1,-1\}`$:


```math
M_{1,s}=\frac{I-s\sigma_y}{4},
\qquad
M_{2,s}=\frac{I+s\sigma_x}{4}.
```


They are positive and sum to $`I`$. Operationally, a fair coin chooses between measuring $`-\sigma_y`$ and measuring $`\sigma_x`$. Retain both the setting and the sign.

Assign estimates


```math
t(1,s)=(2s,0),\qquad t(2,s)=(0,2s).
```


At the origin every outcome has probability $`1/4`$. To first order,


```math
p_{1,s}(\theta)=\frac14+\frac{s\theta_1}{4}
+O(\|\theta\|^2),
\qquad
p_{2,s}(\theta)=\frac14+\frac{s\theta_2}{4}
+O(\|\theta\|^2).
```


Hence $`\mathbb E_\theta t=\theta+O(\|\theta\|^2)`$, and the covariance at the origin is $`2I_2`$. The classical Fisher matrix of these four outcomes is $`I_2/2`$. Averaging independent outcome estimates over $`N_{\mathrm{sense}}`$ probes gives covariance $`2I_2/N_{\mathrm{sense}}`$ at the origin, and total variance $`4/N_{\mathrm{sense}}`$. This realizes the local Holevo value.

The operator covariance $`Z`$ has real part $`I_2`$, whereas the actual outcome covariance is $`2I_2`$. Indeed,


```math
2I_2-Z=
\begin{pmatrix}
1&-i\\
i&1
\end{pmatrix}\succeq0.
```


The extra real covariance is exactly what the joint estimation problem requires. The construction is local to the specified operating point; it is not a globally unbiased estimator for arbitrary rotation angles.

### A singular tangent Gram with good conditioning

The factor columns for the same model are


```math
\Psi=|0\rangle,\qquad
\Phi_1=-\frac{i}{2}|1\rangle,\qquad
\Phi_2=\frac12|1\rangle.
```


With $`V=[\Psi,\Phi_1,\Phi_2]`$,


```math
K=V^\dagger V=
\begin{pmatrix}
1&0&0\\
0&1/4&i/4\\
0&-i/4&1/4
\end{pmatrix}.
```


Its eigenvalues are $`1,1/2,0`$. Yet the occupied probe eigenvalue is one, $`J=I_2`$, and $`\|\Phi_1\|_F^2+\|\Phi_2\|_F^2=1/2`$. Thus $`\nu=j=1`$ are valid promises. The tangent Gram is singular because the two derivative columns are complex-linearly dependent; the two real parameter directions remain identifiable.

A normalized joint preparation is


```math
|\Omega\rangle=
\sqrt{\frac23}|0\rangle_L|0\rangle_S
-\frac{i}{\sqrt6}|1\rangle_L|1\rangle_S
+\frac1{\sqrt6}|2\rangle_L|1\rangle_S.
```


Tracing out the sensor gives $`K^{\mathsf T}/(3/2)`$, whose eigenvalues are $`2/3,1/3,0`$. No failure label is needed in this exact Pauli example. It can be added with zero probability to match the general construction.

All these quantities are elementary to compute classically. The example establishes the distinctions between Fisher information, Holevo cost, occupied weights, and tangent rank; it provides no computational advantage claim.

## From a large probe to a small tangent Gram

For a mixed probe, write


```math
\rho=\Psi\Psi^\dagger,\qquad
D_i=\Phi_i\Psi^\dagger+\Psi\Phi_i^\dagger,
```


where every matrix has $`D`$ rows and $`r`$ columns. Normalization and a trace-preserving tangent require


```math
\|\Psi\|_F^2=1,\qquad
\mathrm{Re}\,\mathrm{Tr}(\Psi^\dagger\Phi_i)=0.
```


For unitary encoding with known generator $`G_i`$, the choice $`\Phi_i=-iG_i\Psi`$ meets the tangent equation. A general factor need not list orthogonal state eigenvectors, but the theorem requires its $`r`$ columns to be independent.

Form $`V=[\Psi,\Phi_1,\ldots,\Phi_p]`$. Its column span $`S`$ has dimension at most $`m=r(p+1)`$, even when $`D=2^N`$. Both $`\rho`$ and every $`D_i`$ are supported in $`S`$. If $`P`$ projects onto $`S`$, replacing $`X_i`$ by $`PX_iP`$ preserves the calibration constraints and gives


```math
Z(X)=Z(PXP)+Q,\qquad
Q_{ij}=\mathrm{Tr}\bigl[\rho X_i(I-P)X_j\bigr].
```


The matrix $`Q`$ is a positive Gram matrix. The Holevo objective is monotone under such additions, so the compression preserves the minimum. This is the [exact compression theorem](../science/holevo_compression_and_stability.md#2-exact-compression).

The Gram $`K=V^\dagger V`$ determines that compressed local model up to a common isometry. Learning all $`D^2`$ state entries is unnecessary. Learning only the support of $`\rho`$ would be insufficient: derivatives can point out of that support, as they do in the pure-qubit example.

This reduction has close metrology precedents. In particular, [Zhou's purification/nuisance-parameter reduction](https://arxiv.org/html/2605.03975v2) already supplies an overlap-based route to a small local problem. The project's additional task is quantitative acquisition and certification of the value from imperfect Gram data.

## Why Bures distance is the acquired accuracy

A small absolute error in a nearly zero matrix eigenvalue can be a much larger error in its square root. The scalar calculation


```math
|\eta-0|=\eta,\qquad
|\sqrt\eta-\sqrt0|=\sqrt\eta
```


already explains the issue. Entrywise Gram accuracy does not automatically give equally accurate factors.

For positive Grams, the homogeneous Bures distance is exactly the minimum Frobenius distance between compatible aligned factors. After embedding into a common row space,


```math
d_B(K,\widehat K)
=\min_U\|V-U\widehat V\|_F.
```


The alignment changes the physical representation, not the model's intrinsic value. It controls the state factor and all derivative factors together. If the aligned joint factor error is at most $`\gamma`$, the normalized state error obeys


```math
\|\rho-\widehat\rho\|_1\le2\gamma.
```


The derivative error is also linear in $`\gamma`$, with the derivative-norm constant retained. This is the connection between the metric learned by the quantum procedure and the stability theorem for the precision value.

The general Bures learner has no occupied-eigenvalue floor. The Holevo application does require an occupied-state floor $`\nu`$ and a Fisher floor $`j`$. Those promises bound estimator norms and make the calibration constraints stable. They do not impose a lower bound on every positive eigenvalue of $`K`$.

## Coherent access and paid conditioning

The input is a circuit and its inverse, not an unlimited supply of unspecified state vectors. For the factor array, a successful component has the form


```math
\frac1\alpha\sum_{i=0}^{p}\sum_{a=1}^{r}
|i,a\rangle_L|V_{i,a}\rangle_S|\mathrm{good}\rangle_F,
\qquad
\alpha^2\ge\mathrm{Tr}K.
```


The columns are unnormalized. In a short linear-combination implementation of known generators, preparation, branch selection, coefficients, controls, ancillas, and reversal are all charged.

A failed branch is retained as one additional label. After reversible cleanup and tracing out the physical workspace, the label state is


```math
\omega=
\frac{K^{\mathsf T}}{\alpha^2}
\oplus
\left(1-\frac{\mathrm{Tr}K}{\alpha^2}\right).
```


Its dimension is $`d_{\mathrm{lab}}=m+1`$, and its rank $`R`$ can be as large as $`m+1`$. Substituting the smaller probe rank $`r`$ for $`R`$ would undercount the tomography cost. The transpose is the result of the partial trace: the $`(j,k)`$ label entry is $`\langle V_k|V_j\rangle`$.

The label state is small, but its purification includes the entire physical sensor. The [Bures theorem](../science/bures_tomography.md) learns the reduced state without reconstructing that environment. An already supplied minimal purification would permit a simpler pure-state tomography baseline; the theorem does not assume that stronger representation.

### Why inverse access changes the accuracy cost

With $`U|0\rangle=|\Omega\rangle`$, a reflection about the prepared state can be implemented as


```math
U(2|0\rangle\langle0|-I)U^\dagger.
```


The inverse permits coherent amplification of a known system subspace. It is not supplied by an ordinary laboratory source merely because that source emits copies. The zero-state phase also acts on the charged workspace; its cost does not disappear when the environment is large.

The adaptive learner resolves successively smaller residual subspaces. If a residual has probability mass $`\tau`$, preparing its normalized component costs approximately $`\tau^{-1/2}`$ original queries. To obtain a subnormalized Frobenius error of order $`\sqrt{\tau h/R}`$, the normalized residual needs error of order $`\sqrt{h/(R\tau)}`$. The inherited coherent Frobenius learner then uses approximately


```math
\frac{d_{\mathrm{lab}}\sqrt R}
{\sqrt{h/(R\tau)}}
=d_{\mathrm{lab}}R\sqrt{\frac{\tau}{h}}
```


calls to the residual preparation. Multiplying by its amplification cost cancels the residual-mass dependence:


```math
d_{\mathrm{lab}}R\sqrt{\frac{\tau}{h}}
\times\tau^{-1/2}
=\frac{d_{\mathrm{lab}}R}{\sqrt h}.
```


At each successful round, the procedure retains directions carrying a constant fraction of the remaining mass. Weighted error accounting controls the coherences between directions retained in different rounds. Taking $`h`$ of order $`\zeta^2`$ divided by a logarithmic round count gives


```math
Q_{\mathrm{Bures}}
=\widetilde O(d_{\mathrm{lab}}R/\zeta).
```


These relations explain the resource cancellation; the [proof](../science/bures_tomography.md#3-a-subnormalized-residual-estimation-lemma) handles the stopping rule, finite amplification error, exact rank preservation, and confidence allocation. Conditional preparation is never supplied at unit cost.

By comparison, optimal copy-access Bures tomography uses order $`d_{\mathrm{lab}}R/\zeta^2`$ copies. Published coherent trace tomography followed by the generic trace-to-Bures conversion also has inverse-square Bures-accuracy dependence. Direct overlap methods remain competitive when a model has additional structure. The comparison is between explicitly different access models and guarantees.

## Turning acquired data into a certified value

Acquired data do not automatically describe an exactly normalized state and trace-preserving tangent. The algorithm first fits the learned Gram to a physical set with explicit slack:


```math
\widehat K\succeq0,\qquad
\mathrm{Tr}\widehat K_{00}=1,\qquad
\widehat K_{00}\succeq\frac\nu2 I_r,
```



```math
\mathrm{Re}\,\mathrm{Tr}\widehat K_{0i}=0,\qquad
\sum_i\mathrm{Tr}\widehat K_{ii}\le2L^2.
```


The true model satisfies the stronger promises $`\nu`$ and $`L^2`$. The relaxed fit makes finite certified computation possible near a boundary. It does not relax the input assumptions without accounting for the change.

The [physical-fit certificate](../science/certified_optimization.md#2-physical-fitting-without-a-hidden-boundary-assumption) uses a positive block matrix and a squared Bures-distance bound. A rational feasible witness establishes proximity; a solver's success flag alone does not.

The stability mechanism is a calibration repair. For a feasible estimator tuple of the true model, form


```math
M_{ij}=\mathrm{Tr}(\widehat D_jX_i),\qquad
\mu_i=\mathrm{Tr}(\widehat\rho X_i).
```


When $`M`$ is sufficiently close to $`I`$, the tuple


```math
\widehat X_i
=\sum_k(M^{-1})_{ik}(X_k-\mu_kI)
```


is exactly centered and locally unbiased for the fitted model. The state-weight and Fisher promises control the estimator norms, the invertibility of $`M`$, and the objective change.

For the current equal-weight theorem, define


```math
b_L=2+(1+\sqrt2)L,
\qquad
A=\frac{768p^{3/2}b_L}{\sqrt\nu\,j^{3/2}}
+\frac{1024p}{\nu j},
```



```math
\gamma_0=
\min\left\{\frac{\nu}{4p},
\frac{\sqrt{\nu j}}{16\sqrt p\,b_L}\right\}.
```


A valid fitted model within factor radius $`\gamma\le\gamma_0`$ satisfies


```math
\widehat J\succeq\frac j8 I_p,
\qquad
|C_H(\widehat K)-C_H(K)|\le A\gamma.
```


Choose $`\gamma\le\min\{\gamma_0,\varepsilon/(2A)\}`$ and label-state tolerance $`\zeta\le\gamma/(8\alpha)`$. The fit guarantees the required factor radius on the acquisition success event.

The classical Holevo program then supplies rational primal and dual witnesses giving $`\ell\le C_H(\widehat K)\le u`$, with $`u-\ell\le\varepsilon`$. The reported interval is


```math
[\ell-A\gamma,\ u+A\gamma].
```


Its midpoint has error at most $`\varepsilon`$. Acquisition supplies the probability statement; the fit and optimization witnesses have deterministic exact checks. The algorithm also specifies failure exits for invalid acquired data.

With $`d_{\mathrm{lab}}=m+1`$, a conservative acquisition count is


```math
Q=
\widetilde O\!\left[
d_{\mathrm{lab}}^2\alpha
\max\left\{\frac{2A}{\varepsilon},\frac1{\gamma_0}\right\}
\right].
```


The [full resource account](../science/circuit_access.md#5-complete-resource-statement) additionally charges preparation gates, workspace operations, small-register tomography gates, finite synthesis, classical fitting, and certificate checking. Inverse-$`\varepsilon`$ precision scaling holds at fixed model size and conditioning. Large $`r,p,\alpha`$, small $`\nu,j`$, or an expensive preparation circuit remain expensive.

## How the sensing example fits the assumptions

The [four-qubit application](../science/applications_and_limits.md#2-a-circuit-family-with-certified-rank-and-explicitly-constructed-access) begins with a parity-preserving circuit $`V`$ and $`|\psi\rangle=V|0000\rangle`$. A parity-flipping Pauli $`P`$ produces an orthogonal partner. The probe


```math
\rho=\frac34|\psi\rangle\langle\psi|
+\frac14P|\psi\rangle\langle\psi|P
```


has known rank two and occupied weights $`3/4,1/4`$. A one-qubit environment explicitly supplies its full-column-rank factor representation. Two local Pauli generators give $`p=2`$, $`L^2=1/2`$, and $`\alpha^2=3/2`$.

The exact rational certificate establishes


```math
J\succeq\frac45 I_2,\qquad
2.250762525\le C_H(I_2)\le2.250762539,
```


while


```math
2.233360231\le\mathrm{Tr}(J^{-1})
\le2.233360232.
```


The tangent span has dimension six inside the sixteen-dimensional physical space. Its numerical compression check agrees with the separate exact full-model enclosure. These are finite algebraic checks of the local bound, not a measured sensing performance or a benchmark of the coherent learner's runtime.

A larger member of the circuit family pays the cost of preparing $`V`$ on every query. If its correlators are efficiently contractible by classical methods, that is the appropriate comparator. A generic large Hilbert space does not establish quantum advantage.

### Three promises with different roles

| Property | Role in the argument |
|---|---|
| Occupied-state floor $`\nu`$ | Bounds estimator action on weakly populated directions. Requires the supplied factor representation to be accounted for. |
| Fisher floor $`j`$ | Keeps the real parameter directions locally distinguishable and permits constraint repair. |
| Full tangent-Gram spectrum | May contain zero or arbitrarily small eigenvalues. Bures acquisition and rational compression accommodate this. |

The [rank-three obstruction](../science/conditioning_limits.md) keeps bounded generators and good Fisher information while an occupied state weight vanishes. Its opaque global preparation oracles imply an $`\Omega(\varepsilon^{-3/2})`$ lower bound. This rules out dropping $`\nu`$ from the uniform inverse-$`\varepsilon`$ claim. It does not give the same lower bound when internal rare-branch rotations or the defining classical parameters are separately available.

Approximate low rank is also a statement about the encoding. For a fixed mixture before the same bounded unitary sensing evolution, the [preparation-tail certificate](../science/applications_and_limits.md#5-a-rigorous-preparation-tail-certificate) bounds the change in $`C_H`$. A tail of weight $`q`$ can change the value at order $`\sqrt q`$, even with $`J=I`$. A small trace-distance error without derivative control does not certify a small precision-value error.

## Research dependencies beyond the anchor

The following are research dependencies for the algorithm, not additional educational anchors.

| Ingredient | Primary dependency | Project-specific connection |
|---|---|---|
| Holevo SDP | [Albarelli, Friel and Datta (2019)](https://arxiv.org/abs/1906.05724) | The certificate preserves the established convex objective and supplies finite rational witnesses. |
| Two-parameter calculation | [Sidhu, Ouyang, Campbell and Kok (2021)](https://arxiv.org/html/1912.09218v4) | The four-qubit verifier uses the established quadratic dual. |
| Purification and nuisance parameters | [Zhou (2026)](https://arxiv.org/html/2605.03975v2) | A close predecessor for evaluation from a small collection of overlaps. |
| Coherent Frobenius tomography | [van Apeldoorn et al., SODA 2023, Theorem 8.2](https://real.mtak.hu/191396/1/1.9781611977554.ch47.pdf) | The residual-state learner invokes this primitive with its preparation cost included. |
| Adaptive geometry and weighted error | [Flammia and O'Donnell (2024)](https://arxiv.org/html/2305.18519v2) | Related spectral peeling and chi-squared control; their copy-access theorem is not silently changed into a query theorem. |
| Mass estimation and residual preparation | [Brassard et al.](https://arxiv.org/abs/quant-ph/0005055) and [Yoder, Low and Chuang](https://arxiv.org/abs/1409.3305) | Amplitude estimation and fixed-point amplification, including every preparation and inverse call. |
| Strong copy comparator | [Pelecanos, Spilecki and Wright](https://arxiv.org/abs/2510.07788); [Pelecanos, Spilecki, Tang and Wright](https://arxiv.org/html/2511.15806v1) | Optimal Bures-copy scaling and a random-purification reduction do not themselves supply a reusable preparation/inverse oracle. |

The [related-work account](../references/RELATED_WORK.md) also compares the current QFI-family characterization and mixed-state geometry algorithms. The [bibliography](../references/references.bib) includes the fidelity-block, weak-SDP, gate-synthesis, and classical-simulation sources used by the full proofs.

## Connecting the route to the complete record

[M01–M12](MILESTONES.md) assemble the estimation, metric, access, and conditioning foundations. The narrative milestones then follow exact compression, label-state preparation, adaptive Bures acquisition, physical fitting, stability, and primal/dual certification. The final examples and limits identify the intended capability and the assumptions under which it is useful.

The [proof files](../science/model.md#reading-route) contain the general statements. The [verification guide](../verification/README.md) identifies which consequences have executable checks. The complete qubit calculation above explains the objective by hand; the four-qubit certificate verifies one less symmetric instance; neither replaces the general tomography or stability proofs.

