# Applications, exact examples, and limits

Certified Holevo-bound evaluation supplies an attainable precision benchmark for a specified multiparameter quantum model. This document develops the useful task, input circuits, matched baselines, and robustness limits. The Bures-tomography theorem and the general factor-stability theorem are proved separately. Results derived here are distinguished from the established metrology and tomography tools cited below.

## 1. The output is a precision benchmark for a known circuit model

For a smooth local model $`\rho_\theta`$, a positive weight matrix $`W`$, and an operating point $`\theta_0`$, the target is a certified interval for the Holevo Cramér–Rao value $`C_H(W)`$. Under the usual regular local asymptotic conditions, this is the attainable weighted covariance benchmark when collective measurements on many identically prepared probes are allowed. It generally differs from the symmetric-logarithmic-derivative benchmark $`\operatorname{Tr}(WJ^{-1})`$. The operational meaning and restrictions are reviewed in [1], and finite-dimensional classical optimization is established in [2,3].

The algorithm evaluates this benchmark for a **specified digital model**. It can compare proposed probe preparations, known preparation errors, control settings, or sensing geometries. It does not by itself produce the optimal measurement, establish a finite-copy error bar for an experiment, or supply coherent access to an uncontrolled laboratory environment.

The input contains a short circuit preparing a purification of a rank-r probe and a coherent implementation of its p local factor derivatives. For unitary sensing, these derivatives are $`\Phi_i=-iG_i\Psi`$, where $`\rho=\Psi\Psi^\dagger`$. The useful label dimension is $`m=r(p+1)`$. All circuit gates, inverses, derivative normalizations, controlled operations, and the small classical optimization are charged. No amplitude loading from an exponentially long classical vector is assumed.

The occupied probe eigenvalues are bounded below by $`\nu>0`$, and $`J\succeq jI`$ with $`j>0`$. These are physical conditioning conditions. They are distinct from the smallest positive eigenvalue of the factor-and-derivative Gram matrix, on which the algorithm imposes no lower bound. In particular, $`r\nu\leq1`$; a thermal state with exponentially many appreciable eigenvalues is not covered by treating $`\nu`$ as constant.

## 2. A circuit family with certified rank and explicitly constructed access

Consider N qubits initially in $`|0\rangle^{\otimes N}`$, with global parity


```math
S=Z_0Z_1\cdots Z_{N-1}.
```


Let $`V`$ be a circuit of local Z rotations and two-qubit XX and ZZ rotations. Every gate commutes with S, so $`|\psi\rangle=V|0^N\rangle`$ has even parity. Choose $`P=X_0`$, which flips parity. The state


```math
\sigma=\tfrac34|\psi\rangle\langle\psi|
       +\tfrac14P|\psi\rangle\langle\psi|P
\tag{1}
```


has exactly two orthogonal occupied eigenvectors and eigenvalues $`3/4,1/4`$, regardless of the circuit depth. The rank and the bound $`\nu=1/4`$ follow from the symmetry; they do not require tomography.

Encode two local field quadratures on qubit 1:


```math
G_1=X_1/2,\qquad G_2=Y_1/2,\qquad
\sigma_\theta=e^{-i\theta_2G_2}e^{-i\theta_1G_1}
\sigma e^{i\theta_1G_1}e^{i\theta_2G_2}.
\tag{2}
```


The remaining qubits are a quantum reference prepared by V. Equation (1) describes a known binary preparation error on that reference before sensing. It is a controlled model of imperfect preparation, not a model of independent noise throughout every gate or every sensing qubit.

### Purification and derivative circuit

Prepare one environment qubit in $`\sqrt{3/4}|0\rangle+\sqrt{1/4}|1\rangle`$, prepare V on the physical register, and apply controlled-P. This gives


```math
|\Psi\rangle=\sqrt{3/4}|0\rangle_E|\psi\rangle
             +\sqrt{1/4}|1\rangle_E P|\psi\rangle.
```


The six unnormalized factor/derivative columns are


```math
v_{0a}=\sqrt{w_a}P^a|\psi\rangle,
\quad v_{1a}=-\tfrac i2 X_1v_{0a},
\quad v_{2a}=-\tfrac i2 Y_1v_{0a},
\qquad (w_0,w_1)=(3/4,1/4).
```


Prepare a three-valued derivative label with amplitudes proportional to $`(1,1/2,1/2)`$, then apply the corresponding controlled Pauli and relative phase. The resulting normalized joint state is


```math
|\Omega\rangle=\frac1\alpha\sum_{i=0}^2\sum_{a=0}^1
 |i,a\rangle v_{ia},\qquad \alpha=\sqrt{3/2}.
\tag{3}
```


There is no postselection in this Pauli-generator example. One call to its preparation unitary costs one V circuit plus a constant number of gates for the environment, derivative label, phases, and controlled Paulis. Its inverse is the reversed circuit. For a general N-qubit, depth-$`\ell`$ local preparation, the gate count is $`G_V=O(N\ell)`$ before fault-tolerant synthesis overhead. Synthesizing the known rotations to the globally allocated error tolerance adds the usual logarithmic accuracy cost; rotations are not treated as exact hardware primitives for a gate-count claim.

Tracing out the N-qubit physical register leaves the six-dimensional label density matrix $`K^T/\alpha^2`$, where $`K_{ia,jb}=\langle v_{ia}|v_{jb}\rangle`$. The classical optimization therefore uses at most a six-dimensional local model.

This is precisely where arbitrary-purifier Bures tomography matters. Although the original probe has a one-qubit environment, the **small label state** in (3) is purified by the entire N-qubit physical register. Pure-state tomography of that supplied purification would learn exponentially many physical amplitudes. The new reduced-state learner reconstructs the six-dimensional label state directly and does not learn or compress those physical amplitudes first.

### Concrete four-qubit instance

The exact arithmetic verifier uses four qubits and the following nine gates, applied in the listed order. Each row denotes $`cI-isP`$, so the rational c and s specify the gate without an angle-rounding convention.

| Order | Pauli P | c | s |
|---|---|---|---|
| 1 | $`X_0X_1`$ | 3/5 | 4/5 |
| 2 | $`Z_1`$ | 5/13 | 12/13 |
| 3 | $`X_1X_2`$ | 8/17 | 15/17 |
| 4 | $`Z_0Z_1`$ | 7/25 | 24/25 |
| 5 | $`X_2X_3`$ | 20/29 | 21/29 |
| 6 | $`Z_2`$ | 12/37 | 35/37 |
| 7 | $`Z_1Z_2`$ | 9/41 | 40/41 |
| 8 | $`X_0X_3`$ | 11/61 | 60/61 |
| 9 | $`Z_0`$ | 28/53 | 45/53 |

For the state and sensing generators in (1)–(2), exact rational primal and dual checks enclose


```math
2.250762525\le C_H(I_2)\le2.250762539,
```


whereas


```math
2.233360231\le\operatorname{Tr}(J^{-1})\le2.233360232.
```


Thus the incompatibility penalty is at least 0.017402294 in these units. The local QFIM is approximately


```math
J=\begin{pmatrix}
0.970712718&-0.068436906\\
-0.068436906&0.840080410
\end{pmatrix},
```


with numerical minimum eigenvalue 0.810793128. Exact rational principal-minor checks certify $`J\succeq(4/5)I_2`$. The exact probe weights remain $`3/4,1/4`$, and both generator norms are $`1/2`$. The derivative array has $`\sum_i\|\Phi_i\|_F^2=1/2`$.

The reproducible verifier and its recorded output are [verify_holevo_circuit.py](../verification/verify_holevo_circuit.py) and [its reference report](../verification/expected/holevo_circuit.json) in the project. Floating-point optimization proposes witnesses; rational arithmetic checks the affine constraints, positive dual quadratic forms, and final objective interval. This is a classical verification of one nontrivial model. It is not a quantum hardware experiment, a simulation of the tomography algorithm, or a demonstrated runtime advantage.

A separate numerical compression check finds a six-dimensional tangent span within the sixteen-dimensional physical space, with projection residual $`8.06\times10^{-16}`$. Optimizing the compressed model gives lower and upper values 2.250762525538578 and 2.2507625593344343, respectively, consistent with the exact full-model enclosure to $`3.4\times10^{-8}`$. These floating-point checks test the dimension reduction; the certified interval above comes from the exact rational full-model witnesses.

### Why the certificate is valid

In the eigenbasis of the rank-two state, $`D_j=-i[G_j,\sigma]`$ has zero diagonal. Removing diagonal entries from a Hermitian estimator preserves its derivative constraints, makes it centered, and subtracts a real positive covariance Gram. Kernel–kernel entries are irrelevant. The remaining off-diagonal Hermitian basis pairs have real covariance weight $`w=\lambda_a+\lambda_b`$ and imaginary antisymmetric weight $`z=\lambda_a-\lambda_b`$.

For two parameters and real coordinate vectors x and y, the Holevo objective is


```math
x^TRx+y^TRy+2|x^TIy|
=\max_{|\beta|\le1}
\begin{pmatrix}x\\y\end{pmatrix}^{T}
\underbrace{\begin{pmatrix}R&\beta I\\-\beta I&R\end{pmatrix}}_{H_\beta}
\begin{pmatrix}x\\y\end{pmatrix}.
```


Here I is real antisymmetric, so $`H_\beta`$ is real symmetric. For $`|\beta|<1`$, each pair block has determinant $`w^2-\beta^2z^2>0`$. If the exact derivative constraints are $`Ax=b`$ on the stacked coordinates, any rational dual multiplier u gives the valid lower bound


```math
2b^Tu-u^TAH_\beta^{-1}A^Tu.
```


Any exactly feasible rational estimator tuple gives the upper objective above. This is the established two-parameter quadratic-dual method [3], implemented as an independent exact arithmetic check. Optimality of the floating-point proposed beta is unnecessary for the validity of the enclosing interval.

## 3. The strongest matched baselines

Let M be the label-state dimension and R its rank; here $`M=6`$, $`R\le6`$. Let A denote the explicit local factor-stability coefficient from the main theorem, so aligned factor error $`\gamma`$ changes $`C_H`$ by at most $`A\gamma`$ within its stated conditioning neighborhood. A label Bures tolerance $`\zeta`$ gives $`\gamma\le\alpha\zeta`$, up to the fixed fitting/normalization factor used by that theorem. Logarithmic confidence factors are suppressed in the comparison below.

| Method | Information acquired and charged cost | Meaning of comparison |
|---|---|---|
| Coherent Bures learner plus small classical optimization | $`\widetilde O(MR\alpha A/\epsilon)`$ preparation/inverse queries; multiply by the joint preparation gate cost and include the learner's small-register gates and classical optimization | Uniform accuracy guarantee without a label-Gram eigenvalue floor |
| Optimal copy-access Bures tomography of the same label state | $`O(MR\alpha^2A^2/\epsilon^2)`$ independently prepared label copies, allowing collective measurements [5] | Strong sample-access comparator, using the same forward preparation circuit |
| Published coherent trace-norm tomography, followed by the generic trace-to-Bures inequality | $`\widetilde O(MR\alpha^2A^2/\epsilon^2)`$ preparation/inverse queries [4] | Strong generic coherent baseline before the Bures upgrade; the square root in metric conversion causes the extra accuracy cost |
| Direct circuit overlaps plus Zhou's reduction and the established classical optimization | Estimate the O(m²) required complex overlaps with Hadamard tests or amplitude estimation; charge each branch/derivative circuit and its precision | Already removes ambient-state reconstruction [6]; no claim that dimension reduction alone is new |
| Classical state-vector, tensor-network, stabilizer, or Gaussian computation of the same overlaps | Use the best available structure of V; the subsequent optimization has the same small dimension | Essential comparator when the circuit is classically tractable [7–10] |

For direct overlap estimation, an entrywise Gram guarantee alone typically certifies factor distance only through a square root near singular Gram matrices. This explains why simply applying amplitude estimation independently to fixed overlap entries does not automatically prove a uniform inverse-$`\epsilon`$ value bound. If the Gram has a known positive eigenvalue floor, or if an analytic formula exposes favorable variance, a direct method can already achieve inverse-$`\epsilon`$ scaling. The present claim is a general guarantee without that additional Gram promise, not a lower bound against every model-specific direct algorithm.

The copy-versus-coherent comparison is an access-model comparison. Coherent algorithms use inverses and long coherent sequences that independently prepared copies do not provide. It is not a claim that an ordinary classical computer equipped with no quantum experiment must take $`\epsilon^{-2}`$ time: it may calculate the relevant correlators deterministically from a tractable circuit description.

For the parity circuit, a state-vector method takes $`O(G_V2^N)`$ arithmetic operations and $`O(2^N)`$ storage, after which the small set of Pauli correlators and the classical optimization are inexpensive. This is only one classical method. In a one-dimensional shallow circuit, tensor-network contractions can be much cheaper; bounded bond dimension gives efficient matrix-product-state simulation [7], and low circuit-treewidth gives efficient contraction [8]. Clifford instances have efficient stabilizer overlap algorithms [9]. Gaussian metrology even admits exact moment-based Holevo programs without a Hilbert-dimension dependence [10]. The four-qubit verifier is deliberately classically easy.

Including non-Clifford ZZ rotations avoids automatically restricting the example to a nearest-neighbor free-fermion circuit. It does not prove hardness of this particular nine-gate example or of every depth family. The worst-case circuit-hardness witness in Section 4 makes a separate, explicit statement.

For this Pauli-generator family, $`\alpha=\sqrt{3/2}`$, r and p are fixed, and the preparation overhead outside V is constant. Consequently the coherent query advantage in the guaranteed high-accuracy regime multiplies the actual cost of V. If V is expensive to prepare, that cost is paid on every query. If V is easy to contract classically, the classical contraction remains the preferred baseline. No ground-state preparation, unknown environmental reset, or free classical data-loading cost is hidden in (3).

## 4. Exact fixed-QFI family and the separate hardness statement

Let $`\sigma`$ be any work-register state, O a Hermitian unitary, and $`c=\operatorname{Tr}(\sigma O)`$. Append a sensor qubit in $`|0\rangle`$, and use


```math
\rho=|0\rangle\langle0|\otimes\sigma,
\quad G_1=X\otimes I/2,
\quad G_2=Y\otimes O/2.
```


At the zero operating point, the exact values are


```math
J=I_2,\qquad
C_H(I_2)=\frac{4}{1+\sqrt{1-c^2}}.
\tag{4}
```


Thus a fixed, perfectly conditioned QFIM leaves a nontrivial precision tradeoff undetermined. The general fact that QFI does not determine simultaneous attainability is established metrology [1]; (4) supplies an explicit algorithmic witness and an exact test family.

### Proof of (4)

In sensor blocks,


```math
\rho=\operatorname{diag}(\sigma,0),\quad
D_1=\begin{pmatrix}0&i\sigma/2\\-i\sigma/2&0\end{pmatrix},\quad
D_2=\begin{pmatrix}0&\sigma O/2\\O\sigma/2&0\end{pmatrix}.
```


Writing estimators as $`X_i=\left(\begin{smallmatrix}A_i&B_i^\dagger\\B_i&C_i\end{smallmatrix}\right)`$, remove $`A_i,C_i`$: the constraints do not involve them, and their covariance contribution is a positive Gram. Set $`x_i=B_i\sqrt\sigma`$, $`u_1=-i\sqrt\sigma`$, $`u_2=O\sqrt\sigma`$. The constraints are $`\operatorname{Re}\langle u_j,x_i\rangle=\delta_{ij}`$, while $`\|u_i\|=1`$ and $`\langle u_1,u_2\rangle=ic`$. The SLD information is therefore $`J=I_2`$.

Projection onto the complex span of $`u_1,u_2`$ only removes another positive covariance Gram. For $`|c|<1`$, write $`u_1=e_1`$, $`u_2=ice_1+se_2`$, $`s=\sqrt{1-c^2}`$, and $`y=c/(1+s)`$. The feasible pair


```math
x_1=e_1+iye_2,\qquad x_2=iye_1+e_2
```


has $`\langle x_1,x_2\rangle=0`$ and objective $`2(1+y^2)=4/(1+s)`$. For a matching lower bound, use the positive quadratic form


```math
q_y=\|x_1\|^2+\|x_2\|^2+2y\operatorname{Im}\langle x_1,x_2\rangle.
```


It is bounded above by the Holevo objective. The displayed feasible pair obeys the constrained stationarity equations $`x_1-iyx_2=(1+y^2)u_1`$, $`x_2+iyx_1=(1+y^2)u_2`$. Convexity proves optimality. The vectors are physically realizable because they are linear combinations of $`\sqrt\sigma`$ and $`O\sqrt\sigma`$; for example $`B_1=i(yO-I)/s`$, $`B_2=(O-yI)/s`$. The endpoints have value 4 by direct evaluation of the constraints.

### Short circuit descriptions can still be hard

Take


```math
\sigma=\tfrac12\bigl(|\psi_U\rangle\langle\psi_U|
+O|\psi_U\rangle\langle\psi_U|O\bigr),
```


where $`|\psi_U\rangle`$ is a short universal-circuit output. The occupied eigenvalues are $`(1\pm|c|)/2`$. Given a bounded-error decision circuit with acceptance probability a, append an output qubit whose Z expectation is 1/3 on the reject branch and 2/3 on the accept branch. Then $`O=Z`$ on that qubit gives $`c=(1+a)/3\in[1/3,2/3]`$, so $`\nu\ge1/6`$, the generator norms are 1/2, and J remains the identity. A constant-accuracy evaluation of the strictly increasing expression (4) distinguishes $`a\le1/3`$ from $`a\ge2/3`$. A general classical polynomial-time circuit-description algorithm for this task would therefore imply $`\mathrm{BQP}\subseteq\mathrm{BPP}`$.

This is a worst-case computational witness. It is not evidence that the four-qubit demonstration is hard. The special family itself has a direct quantum solution: estimate c and use (4). Its closed form must be included as a stronger comparator than the general tomography-and-optimization pipeline.

### Precision lower bounds with conditioning fixed

Choose $`O=Z`$, $`\sigma_t=\operatorname{diag}(\cos^2t,\sin^2t)`$, with purification $`|\Psi_t\rangle=\cos t|00\rangle+\sin t|11\rangle`$. The sensing generators remain fixed. For $`t\in(0,\pi/4)`$,


```math
C_H(t)=4/[1+\sin(2t)].
```


At $`t_0=\pi/6`$, the occupied weights are $`3/4,1/4`$, and the derivative of this value is nonzero. The unitaries $`U_t=\mathrm{CNOT}(R_y(2t)\otimes I)`$ differ in operator norm by $`\Theta(|\Delta t|)`$. The hybrid argument therefore gives $`\Omega(1/\epsilon)`$ preparation/inverse queries for additive-$`\epsilon`$ estimation. Purification copies have overlap $`\cos\Delta t`$; constant-success discrimination requires $`\Omega(1/\epsilon^2)`$ copies. These are black-box accuracy lower bounds. If the only unknown angle t is explicitly supplied as a classical number, this subfamily is easy to evaluate classically.

## 5. A rigorous preparation-tail certificate

Assume a fixed mixture before the same local unitary sensing evolution:


```math
\rho_\theta=(1-\delta)\sigma_\theta+\delta\tau_\theta,
\quad \sigma_\theta=U_\theta\sigma U_\theta^\dagger,
\quad \tau_\theta=U_\theta\tau U_\theta^\dagger,
```


with $`0\le\delta<1`$, $`\|G_j\|\le g`$, retained state $`\lambda_{\min}^+(\sigma)\ge\kappa>0`$, and $`C=C_H(\sigma,I_p)<\infty`$. No bound is imposed on the rank of $`\tau`$. Define


```math
e_u=\delta\left(1+2g\sqrt{2pC/\kappa}\right),\qquad
\overline C=\frac{C(1-\delta+4\delta/\kappa)}{(1-e_u)^2},
```


and


```math
e_l=\frac{\delta+2g\sqrt{p\delta\overline C}}{1-\delta}.
```


Whenever $`e_u<1`$ and $`e_l<1`$,


```math
(1-\delta)(1-e_l)^2 C
\ \le C_H(\rho,I_p)\ \le\overline C.
\tag{5}
```


This is a dimension-independent, two-sided certificate for the full noisy preparation. It retains the actual generator norms, parameter count, retained-state conditioning, and precision value.

### Proof

For a positive covariance Gram Z, let $`h(Z)=\operatorname{Tr}\operatorname{Re}Z+\|\operatorname{Im}Z\|_1`$. It is monotone under positive Gram additions, subadditive, satisfies $`\operatorname{Tr}Z\le h(Z)\le2\operatorname{Tr}Z`$, and obeys $`h(RZR^T)\le\|R\|^2h(Z)`$ for real R.

For the upper bound, take a retained-model optimal estimator X and set its kernel–kernel blocks to zero. Since the retained state is at least $`\kappa`$ on its support,


```math
\sum_i\|X_i\|_{\mathrm{op}}^2
\le\sum_i\|X_i\|_F^2
\le2C/\kappa.
```


Write $`A(\omega)_{ij}=\operatorname{Tr}[X_i(-i[G_j,\omega])]`$. Then $`A(\sigma)=I`$ and $`\|A(\tau)\|\le2g\sqrt{2pC/\kappa}`$. Therefore $`\|A(\rho)-I\|\le e_u`$, so $`A(\rho)^{-1}X`$ repairs local unbiasedness. Centering afterward only removes a positive Gram. Subadditivity gives the pre-repair objective at most $`C(1-\delta+4\delta/\kappa)`$; the repair multiplies it by at most $`(1-e_u)^{-2}`$.

For the lower bound, take a full-model optimal estimator with value $`C_\rho`$. Weighted Cauchy–Schwarz gives


```math
|\operatorname{Tr}[\tau(-i[X_i,G_j])]|\le2g\sqrt{\operatorname{Tr}(\tau X_i^2)},
```


and hence $`\|\delta A(\tau)\|_F\le2g\sqrt{p\delta C_\rho}`$. The identity $`I=(1-\delta)A(\sigma)+\delta A(\tau)`$, together with $`C_\rho\le\overline C`$, bounds $`\|A(\sigma)-I\|\le e_l`$. Repairing the retained-model constraints and using positive-Gram monotonicity yields $`C\le C_\rho/[(1-\delta)(1-e_l)^2]`$. Minimizing sequences suffice if an optimum is not initially assumed attained.

The square-root dependence is necessary. In (4), choose $`O=Z`$ and $`\sigma_\delta=(1-\delta)|0\rangle\langle0|+\delta|1\rangle\langle1|`$. Then


```math
C_H(\delta)=\frac4{1+2\sqrt{\delta(1-\delta)}}
=4-8\sqrt\delta+O(\delta),
```


while the trace distance from the retained pure model is $`\delta`$, $`J=I_2`$, and $`g=1/2`$. A small tail can therefore materially change incompatibility despite perfect Fisher conditioning.

### Physical limits of (5)

The common encoding assumption is essential. A trace-close surrogate state with unrelated derivatives, parameter-dependent mixture weights, or noise acting during sensing does not automatically satisfy it. The mixture is before the encoding; its effect does not contradict monotonicity under a parameter-independent channel applied after encoding.

For fixed other parameters, additive-$`\epsilon`$ truncation generally requires $`\delta=O(\epsilon^2)`$. The constants can be conservative. For the parity family, $`p=2`$, $`g=1/2`$, $`\kappa=1/4`$, giving $`e_u=\delta(1+4\sqrt C)`$ and $`\overline C=C(1+15\delta)/(1-e_u)^2`$. A realistic application must evaluate this interval, rather than describe every small noise probability as negligible.

Independent local errors generally produce exponential exact rank. Truncation to at most k events has up to $`O(N^k)`$ branches for fixed k, and k can increase as the requested tail tolerance decreases. A small error rate per qubit does not by itself bound the total tail. Extensive unnormalized generators also make g grow with N; such growth must be retained in (5) and in the derivative-access normalization.

## References

1. R. Demkowicz-Dobrzański, W. Górecki, M. Guţă, *Multi-parameter estimation beyond Quantum Fisher Information*, J. Phys. A 53, 363001 (2020). https://arxiv.org/abs/2001.11742
2. F. Albarelli, J. F. Friel, A. Datta, *Evaluating the Holevo Cramér–Rao Bound for Multiparameter Quantum Metrology*, PRL 123, 200503 (2019). https://arxiv.org/abs/1906.05724
3. J. S. Sidhu, Y. Ouyang, E. T. Campbell, P. Kok, *Tight bounds on the simultaneous estimation of incompatible parameters*, PRX 11, 011028 (2021). https://arxiv.org/abs/1912.09218
4. J. van Apeldoorn, A. Cornelissen, A. Gilyén, G. Nannicini, *Quantum tomography using state-preparation unitaries*, SODA 2023. Published Theorems 8.2 and 9.3: https://real.mtak.hu/191396/1/1.9781611977554.ch47.pdf
5. A. Pelecanos, J. Spilecki, J. Wright, *The debiased Keyl's algorithm: a new unbiased estimator for full state tomography*, STOC 2026, arXiv:2510.07788. https://arxiv.org/abs/2510.07788
6. S. Zhou, *Quantum metrology of mixed states via purification*, arXiv:2605.03975v2 (2026). https://arxiv.org/html/2605.03975v2
7. G. Vidal, *Efficient classical simulation of slightly entangled quantum computations*, PRL 91, 147902 (2003). https://arxiv.org/abs/quant-ph/0301063
8. I. L. Markov, Y. Shi, *Simulating quantum computation by contracting tensor networks*, SIAM J. Comput. 38, 963–981 (2008). https://arxiv.org/abs/quant-ph/0511069
9. S. Aaronson, D. Gottesman, *Improved Simulation of Stabilizer Circuits*, PRA 70, 052328 (2004). https://arxiv.org/abs/quant-ph/0406196
10. S. Chang, M. G. Genoni, F. Albarelli, *Multiparameter quantum estimation with Gaussian states: efficiently evaluating Holevo, RLD and SLD Cramér–Rao bounds*, Communications Physics 9, 126 (2026), arXiv:2504.17873. https://arxiv.org/abs/2504.17873
