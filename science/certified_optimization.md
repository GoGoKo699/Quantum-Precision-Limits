# Certified finite classical optimization

## What is proved and what is implemented

This note completes the finite classical part of the algorithm in [holevo_compression_and_stability.md](holevo_compression_and_stability.md). It specifies two finite rational semidefinite programs, explicit interior points, bounded search regions, and independently checkable certificates. It does not assume that a floating-point solver's convergence message is a proof.

The generic optimization algorithm is a theorem-level procedure using rational weak SDP optimization with known inner and outer bounds. The accompanying `verify_certified_optimization.py` implements exact rational certificate checks and analytic witness construction for two nontrivial test models, including a singular tangent Gram. It is not a generic SDP solver. The repository's separate circuit verifier provides a further certified two-parameter circuit demonstration.

The Holevo SDP framework is established by [Albarelli, Friel and Datta (2019)](https://arxiv.org/html/1906.05724v1). The fidelity block SDP is standard; see [Watrous (2012), Section 2.1](https://arxiv.org/html/1207.5726v2). Polynomial-time weak optimization is used only after supplying finite rational inner and outer bounds in the relevant affine space. An explicit statement of this requirement and the resulting rational-output theorem appears in [Ciardo and Živný, Theorem 8.2](https://www.cs.ox.ac.uk/standa.zivny/publications/cz25sicomp-semidefinite-preprint.pdf), citing the [Grötschel–Lovász–Schrijver ellipsoid theorem](https://doi.org/10.1007/BF02579273). An arbitrary boundary SDP is not automatically covered by that result.

For exact checks, a complex rational Hermitian matrix can be represented by its real and imaginary rational parts. Positive semidefiniteness is decidable by exact pivoted LDL factorization, or by the standard realification. No tolerance is used to declare a negative eigenvalue nonnegative.

## 1. Data and the common accuracy budget

Let p>=1, r>=1, m=r(p+1). The input promises are occupied state weights at least nu>0, Fisher matrix J>=jI_p, and sum_i||Phi_i||F²<=L². Assume the factor has r independent columns; necessarily r nu<=1. The circuit normalization is alpha. Classical bounds may be supplied as conservative rationals; square roots in accuracy formulas can be rounded outward.

Use


```math
b_L=2+(1+\sqrt2)L,
\quad
A=\frac{768p^{3/2}b_L}{\sqrt\nu\,j^{3/2}}+
\frac{1024p}{\nu j},
```



```math
\gamma_0=\min\left\{\frac\nu{4p},
\frac{\sqrt{\nu j}}{16\sqrt p\,b_L}\right\},
\quad
0<\gamma\le\min\{\gamma_0,\epsilon/(2A)\},
\quad \zeta\le\gamma/(8\alpha).
\tag{1}
```


Choose gamma and zeta dyadic with only a constant-factor reduction from these bounds. The tomography subroutine outputs a finite rational PSD label-state estimate with its stated Bures confidence guarantee; its rounding and normalization error is part of zeta, not an extra omitted operation. Pinching and rescaling its good block produces a rational PSD preliminary Gram Kpre satisfying, on that event,


```math
d_B(K_{\rm true},K_{\rm pre})\le\alpha\zeta\le\gamma/8.
\tag{2}
```


If alpha² is not initially rational, its circuit normalization and classical scaling can be implemented using a conservative finite rational normalization, with the associated circuit approximation charged in the preparation budget. The certification statements below take their finite classical inputs as exact rationals.

## 2. Physical fitting without a hidden boundary assumption

### 2.1 Regularize the reference, not the unknown physical state

Set


```math
\lambda=\frac{\gamma^2}{64m},\qquad
S=K_{\rm pre}+\lambda I_m.
\tag{3}
```


Then S>=lambda I_m. Since Kpre commutes with its added identity,


```math
d_B(K_{\rm pre},S)^2
=\sum_a(\sqrt{s_a+\lambda}-\sqrt{s_a})^2\le m\lambda,
```


and therefore


```math
d_B(K_{\rm true},S)\le\gamma/4.
\tag{4}
```


This is a classical regularization of the fit's reference matrix. It does not change the input sensing state or assume that its tangent Gram is full rank.

### 2.2 The fit SDP and its certificate

Minimize


```math
F(K,Z)=\operatorname{Tr}K+\operatorname{Tr}S
-2\operatorname{Re}\operatorname{Tr}Z
\tag{5}
```


over complex rational matrices K=K† and Z, subject to


```math
\begin{pmatrix}K&Z\\Z^\dagger&S\end{pmatrix}\succeq0,
\quad \operatorname{Tr}K_{00}=1,
\quad K_{00}\succeq\frac\nu2I_r,
\quad\operatorname{Re}\operatorname{Tr}K_{0i}=0,
\quad\sum_i\operatorname{Tr}K_{ii}\le2L^2.
\tag{6}
```


The block PSD condition already implies K>=0. The fidelity SDP gives


```math
d_B(K,S)^2=\min_Z F(K,Z)
```


under the block constraint. Consequently **any exactly feasible rational pair** satisfying


```math
F(K,Z)\le\gamma^2/4
\tag{7}
```


is an acceptable fit certificate. It directly implies d_B(K,S)<=gamma/2. Combined with (4), the accepted fitted Gram obeys


```math
d_B(K_{\rm true},K)\le3\gamma/4\le\gamma.
\tag{8}
```


A fit optimality certificate is unnecessary. The relevant output is a feasible pair with the explicit objective bound (7). The true Gram is feasible and its optimum in (5) is at most gamma²/16, leaving strictly positive objective slack for finite optimization and certificate repair.

### 2.3 Known strict feasibility

Let


```math
K_c=\operatorname{diag}\left(\frac1rI_r,
\frac{L^2}{pr}I_{pr}\right),\qquad Z_c=0.
\tag{9}
```


This satisfies all affine equations. The block matrix in (6) is positive definite, the base-floor slack is at least 1/(2r), and the derivative-budget slack equals L². There is no difficulty at nu=1/r because the fitted floor was relaxed to nu/2.

The case L=0 is excluded by identifiability: for every direction, the SLD identity and Cauchy–Schwarz give J_ii<=4||Phi_i||F², hence L²>=pj/4>0. A common positive matrix/scalar slack at (9) is at least


```math
\mu_{\rm fit}=\min\left\{\lambda,\frac1{2r},\frac{L^2}{pr}\right\}.
\tag{10}
```


Converting this to a Euclidean inner radius only adds explicit polynomial factors from the affine maps and their dimensions. For example, a perturbation of K changes the derivative trace by at most sqrt(pr)||Delta K||F, so the scalar budget can be protected by a radius at most L²/(2sqrt(pr)). Redundant affine equations are removed by exact Gaussian elimination.

The feasible set is bounded: Tr K<=T:=1+2L², while block positivity gives


```math
\|Z\|_F\le\sqrt{\operatorname{Tr}K\operatorname{Tr}S}
\le\sqrt{T\operatorname{Tr}S}.
\tag{11}
```


Thus a rational outer radius is obtained from T and Tr S. Equations (9)–(11) provide the required relative Slater point and inner/outer bounds for rational weak SDP optimization. Solving to an objective tolerance smaller than the gap between gamma²/16 and gamma²/4, and checking the resulting rational certificate, is a finite algorithm polynomial in matrix dimensions, input bit length, log(1/gamma), and the logarithms of the displayed bounds.

### 2.4 How to certify a numerical proposal

Numerical optimization may propose K,Z. Round to complex rationals, Hermitize K, and enforce the affine equations exactly. The trace normalization can be corrected by adding a scalar multiple of I_r to K00; each Re Tr K0i can be corrected by a scalar multiple of I_r in the corresponding cross block and its adjoint. Alternatively, use an exact rational affine projection after removing redundant equations.

Compute certified lower bounds on the PSD slacks and the scalar budget. If a current matrix slack is at least -eI and the corresponding slack at (9) is at least mu I, mixing with (9) by theta>=e/(mu+e) repairs it. Apply the maximum required theta across the constraints. The affine equations remain exact, and the objective changes linearly, so the cost of repair is explicitly known. Accept only if (7) still holds.

With the known positive slack, proposals can be requested sufficiently accurately that this repair has arbitrarily small objective cost. The required precision has polynomial bit length. A failed check triggers refinement, not an unsupported success statement.

## 3. Rational compression without numerical square roots

Let Kfit be an accepted rational Gram. Exact pivoted LDL factorization gives


```math
K_{\rm fit}=R^\dagger D R,
\tag{12}
```


where q=rank Kfit, R is a rational q-by-m matrix, and D is a positive rational diagonal q-by-q matrix. Zero pivots are handled exactly. For PSD matrices a zero residual diagonal forces its entire residual row and column to vanish; otherwise a positive diagonal pivot is available. This gives a constructive polynomial-time exact factorization over the rationals.

Split R into q-by-r blocks R0,R1,...,Rp. A conceptual orthonormal-space factor is


```math
\Psi=\sqrt D R_0,\qquad\Phi_i=\sqrt D R_i.
\tag{13}
```


The optimization will never need to evaluate sqrt(D). All its coefficients are rational expressions in R and D. K00>=nu I_r/2 implies that R0 has independent columns, although the full Kfit can be singular.

The exact rank q is a property of a finite rational fitted matrix, not a claim to recover the unknown input's tangent rank. No small input eigenvalue is measured to determine it. Positive LDL pivots can be very small, but their bit lengths, and the bit lengths of their inverses, are polynomially bounded by the finite rational input length through determinant bounds. This is a bit-complexity statement, not a condition-number-independent claim about ordinary floating-point factorization.

## 4. Action variables and the exact small Holevo SDP

### 4.1 Hermitian realizability

Replace an estimator by its action x_i=X_i Psi. Because Psi has full column rank, a q-by-r matrix x is the action of a Hermitian X if and only if


```math
\Psi^\dagger x=x^\dagger\Psi.
\tag{14}
```


Necessity is immediate. For sufficiency, with Psi^+=(Psi†Psi)^-1Psi†, set


```math
X=x\Psi^++(x\Psi^+)^\dagger
-(\Psi^+)^\dagger(\Psi^\dagger x)\Psi^+.
\tag{15}
```


The middle matrix in the final term is Hermitian by (14); direct multiplication gives X Psi=x. Thus no estimator directions are lost by using action variables.

Write x_i=sqrt(D)B_i, where B_i is a complex q-by-r variable. The exact affine constraints become


```math
R_0^\dagger D B_i\text{ is Hermitian},\qquad
\operatorname{Re}\operatorname{Tr}(R_0^\dagger D B_i)=0,
```



```math
2\operatorname{Re}\operatorname{Tr}(R_j^\dagger D B_i)=\delta_{ij}.
\tag{16}
```


Covariance is


```math
Z_{ij}=\operatorname{Tr}B_i^\dagger D B_j.
\tag{17}
```


Let N=qr, y_i=vec(B_i) using column-major vectorization, Y=[y_1,...,y_p], and


```math
Q=I_r\otimes D^{-1}\succ0.
\tag{18}
```


Then Z=Y†Q^-1Y. Stack real and imaginary parts of Y into a real vector y. Equations (16) are a rational real system A_eq y=b. Remove linearly redundant rows exactly, so A_eq has full row rank. On the tomography success event, consistency follows from the transferred Fisher condition. Exact elimination must first check the full augmented system for consistency; if it is inconsistent, return a failure marker. This handles an unsuccessful acquisition without applying an undefined inverse.

### 4.2 Primal program

The fitted Holevo value is exactly


```math
\begin{array}{ll}
\text{minimize}&\operatorname{Tr}C\\
\text{subject to}&A_{\rm eq}y=b,\quad C\in\mathbb S^p,\\
&M(C,Y):=\begin{pmatrix}C&Y^\dagger\\Y&Q\end{pmatrix}\succeq0.
\end{array}
\tag{19}
```


The Schur complement is valid because Q>0, and imposes C>=Z. The epigraph identity in the companion note completes the proof. The block order is p+qr<=p+mr. Before affine elimination, there are at most p(r²+p+1) real equations and 2qrp action variables. These dimensions depend polynomially on r,p, not on the original physical Hilbert dimension.

### 4.3 Explicit primal strict feasibility and boundedness

Solve the rational equality-constrained positive quadratic problem


```math
\min_{A_{\rm eq}y=b}\operatorname{Tr}Y^\dagger Q^{-1}Y.
\tag{20}
```


After removing dependent constraints, its solution is an exact rational linear solve. Equivalently, for the real quadratic metric H>0,


```math
y_0=H^{-1}A_{\rm eq}^T
(A_{\rm eq}H^{-1}A_{\rm eq}^T)^{-1}b.
```


Its value equals Tr Jhat^-1, and is at most 8p/j on the tomography success event. If the exactly computed value exceeds 8p/j, return a failure marker. This abort cannot occur on that event and keeps the stated bounds valid on every executed branch. This is the standard SLD-dual minimum real covariance; it follows as well from Hilbert-space projection onto the span of the SLDs. Let Z0=Y0†Q^-1Y0 and take


```math
C_0=(\operatorname{Tr}Z_0+1)I_p.
\tag{21}
```


Then C0-Z0>=I_p, so M(C0,Y0)>0. An explicit positive lower bound is


```math
\mu_H\ge
\frac{\min\{1,\lambda_{\min}(Q)\}}
{(1+\|Q^{-1}Y_0\|_{\rm op})^2}.
\tag{22}
```


Use a rational lower bound for this expression if needed. It follows from the block triangular congruence with diagonal blocks C0-Z0 and Q. Frobenius upper bounds may replace operator norms.

Add the harmless objective sublevel Tr C<=U, with U=2Tr C0+2. It contains an optimum and leaves strict slack at (21). On this set, C>=Z>=0 implies


```math
\|C\|_F\le U,\qquad
\|Y\|_F^2\le U\lambda_{\max}(Q).
\tag{23}
```


This supplies a finite rational outer bound. Relative affine coordinates and their norm-conversion factors can be computed by rational Gaussian elimination. Even if some LDL pivots are tiny, the logarithms of the inner and outer bounds have polynomial bit length in the rational fitted data. Thus rational weak SDP optimization applies without an additional promised tangent-Gram spectral gap.

## 5. An explicit dual and an independently checkable interval

Use real vectorization vec_R(Y) that stacks all real entries followed by all imaginary entries. A dual witness consists of a real vector mu and a complex rational Hermitian block matrix


```math
T=\begin{pmatrix}T_{00}&B^\dagger\\B&T_{11}\end{pmatrix}\succeq0
\tag{24}
```


satisfying


```math
\operatorname{Re}T_{00}=I_p,
\qquad A_{\rm eq}^T\mu+2\operatorname{vec}_{\mathbb R}(B)=0.
\tag{25}
```


Its lower value is


```math
\ell=\mu^Tb-\operatorname{Tr}(T_{11}Q).
\tag{26}
```


For any primal feasible C,Y,


```math
\operatorname{Tr}C-\ell
=\operatorname{Tr}[T M(C,Y)]\ge0.
\tag{27}
```


Indeed the off-diagonal trace is 2 Re Tr B†Y=-mu^T A_eq y=-mu^T b, and Re T00=I pairs with real C. Equation (27) proves weak duality directly, without reliance on a solver or on numerical complementarity.

A certified fitted-model interval is obtained by accepting rational primal and dual witnesses satisfying (19), (24), and (25), and returning


```math
[\ell,u],\qquad u=\operatorname{Tr}C.
\tag{28}
```


The certificate checker verifies every equality exactly, both PSD conditions exactly, and the rational inequality u-ell<=epsilon_H. For a general positive definite rational W, replace Tr C by Tr WC and Re T00=I by Re T00=W; the same argument applies.

### 5.1 Existence and finite computation of arbitrarily small gaps

The primal is strictly feasible by (21). The dual has the explicit strict point mu=0 and T=diag(I_p,I_N), with value ell0=-Tr Q. Standard SDP strong duality therefore gives equal attained primal and dual optimal values.

Finite bounded search regions can also be supplied to the dual. For any dual point of value at least ell0, (27) evaluated at the primal strict point implies


```math
\operatorname{Tr}T\le(\operatorname{Tr}C_0-\ell_0)/\mu_H.
\tag{29}
```


Since A_eq has full row rank, (25) bounds mu in terms of B and the explicitly computable norm of (A_eq^T)^+. These bounds have polynomial bit length. For an actual bounded dual program impose the slightly weaker cutoff ell>=ell0-1, replace the numerator in (29) by Tr C0-ell0+1, and choose strictly larger rational boxes. This cutoff leaves the known strict dual point in the relative interior and preserves an optimum. Applying rational weak optimization to the bounded primal and dual produces feasible witnesses with an arbitrarily small requested gap. This proves a finite classical certification algorithm, not merely equality of real-valued infima.

### 5.2 Repair of numerical proposals

For a proposed primal Y, rationally project vec_R(Y) to A_eq y=b. Compute the exact rational Schur complement Z=Y†Q^-1Y. Increase a real symmetric proposed C by tI_p until C-Z>=0 is certified. A rational bisection with exact PSD tests finds a sufficient t to any requested additive slack, and increases the upper value by exactly pt.

For a proposed dual, set Re T00=I exactly, choose B from (25) using the proposed rational mu, and Hermitize the remaining free blocks. If its PSD violation is at most eI, mix it with diag(I,I) by theta>=e/(1+e), scaling mu by 1-theta at the same time. All affine dual equations remain exact. The lower value changes by a known convex combination with -Tr Q. Refine the proposal until the repaired primal–dual gap meets the target.

These repairs are optional practical routes. The rational weak-optimization theorem with the explicit interior points already guarantees that acceptable finite witnesses exist and can be found with polynomial bit complexity. Neither route accepts a numerical tolerance as a mathematical PSD certificate.

## 6. Combined certificate and computational costs

The classical certificate package contains:

1. The rational preliminary reference S and the accepted physical fit Kfit,Zfit, including the exact inequalities (6)–(7).
2. The exact factorization Kfit=R†DR, with positive rational D and rational R.
3. The rational independent affine system A_eq y=b and its construction from (16).
4. A primal witness C,Y and dual witness T,mu, with the exact gap check.
5. The promised bounds nu,j,L, the normalization alpha, the chosen gamma,zeta, and the tomography confidence statement.

On the tomography success event, the stability theorem enlarges (28) to [ell-A gamma,u+A gamma]. With gamma<=epsilon/(2A) and epsilon_H<=epsilon, its midpoint estimates the true value to additive epsilon. An implementation may choose a smaller classical gap to narrow the reported interval further.

Let B_in denote the finite rational input length, including the tomography output and conservative parameter bounds. The fit uses an LMI of order 2m, a base-floor block of order r, and a scalar derivative bound. The Holevo program uses an LMI of order p+qr and O(p(r²+p+1)) affine constraints. Exact elimination, LDL factorization, certificate checking, and rational weak optimization run in time polynomial in these dimensions, B_in, and log(1/epsilon_H), with the explicitly displayed inner/outer bounds. Dependence on nu,j,L,alpha and gamma enters through their finite descriptions and the logarithms of those bounds. This statement charges the polynomial classical work; it does not assert that a dense generic SDP is fast for large r or p.

Quantum acquisition remains


```math
\widetilde O\left((m+1)^2\alpha
\max\{A/\epsilon,1/\gamma_0\}\right)
```


joint-preparation calls. Classical fitting to squared-distance accuracy of order gamma² does not double the quantum precision exponent. It increases the number of classical precision bits by O(log(1/gamma)), with polynomial arithmetic costs.

## 7. Implemented checks and their limits

`verify_certified_optimization.py` verifies all constraints using exact complex rationals for:

- A pure two-parameter model with a singular 3-by-3 tangent Gram of rank two. Its exact primal upper and dual lower are both 4, while the SLD cost is 2.
- The rational compatibility example c=3/5, sqrt(1-c²)=4/5. Its exact primal upper and dual lower are both 20/9, again with SLD cost 2.

Both cases include a feasible rational physical-fit witness against the regularized reference. With gamma=1/100, the fit objective is exactly 1/640000, below the threshold 1/40000. The verifier's Schur-complement recursion checks PSD even when the certificate matrices or tangent Gram are singular.

These tests validate the certificate formulas and the treatment of a boundary Gram. They do not simulate quantum tomography, establish an end-to-end physical resource advantage, or implement a generic SDP optimizer. The independent circuit demonstration and the analytic general theorem have separate roles and should be described separately.
