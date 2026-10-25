# Exact local compression and quantitative Holevo stability

## Scope and conventions

This note proves the geometric part of the certified algorithm. The numerical realization, including singular tangent Grams and rational certificates, is specified in [certified_optimization.md](certified_optimization.md). The final theorem uses a slightly relaxed physical fit to give that classical computation an explicit interior point.

The output is the local Holevo Cramér–Rao value for a specified statistical model. It is a precision benchmark; this result does not construct a measurement attaining that value. The standard Holevo SDP is inherited from [Albarelli, Friel and Datta, PRL 123, 200503 (2019)](https://arxiv.org/html/1906.05724v1). The relation between mixed-state estimation and purification with environment nuisance parameters is established in [Zhou, arXiv:2605.03975v2](https://arxiv.org/html/2605.03975v2). Exact compression alone is not claimed as an independent new metrological principle.

All matrix norms below are Schatten norms unless an operator norm is specified. For positive semidefinite, possibly subnormalized matrices, use the generalized Bures distance


```math
d_B(A,B)^2=\mathrm{Tr}A+\mathrm{Tr}B
-2\mathrm{Tr}\sqrt{\sqrt A B\sqrt A}.
\tag{1}
```


Equivalently, it is the minimum Frobenius distance between aligned factors of A and B, allowing a common enlarged row space. This convention is the square root of twice one minus root fidelity for normalized states. Infidelity conventions differ by constant factors and must not be substituted without conversion.

At a fixed parameter point, let


```math
\rho=\Psi\Psi^\dagger,\qquad
D_i=\Phi_i\Psi^\dagger+\Psi\Phi_i^\dagger,
\quad i=1,\ldots,p,
\tag{2}
```


where Psi and Phi_i have D rows and r columns, with D the physical Hilbert-space dimension. Assume


```math
\mathrm{Tr}\Psi^\dagger\Psi=1,\qquad
\mathrm{Re}\mathrm{Tr}\Psi^\dagger\Phi_i=0.
\tag{3}
```


These conditions make rho normalized and D_i traceless. Any such tangent is locally realizable by normalizing the factor Psi+sum_i theta_i Phi_i. Let


```math
V=[\Psi,\Phi_1,\ldots,\Phi_p],\quad
m=r(p+1),\quad K=V^\dagger V,
\quad\sum_i\|\Phi_i\|_F^2\le L^2.
\tag{4}
```


For the final algorithm, Psi has r linearly independent columns, the occupied state eigenvalues are at least nu>0, and the symmetric-logarithmic-derivative Fisher matrix obeys J>=j I_p with j>0. Equivalently, K_00>=nu I_r. A redundant environment with more columns than rank rho does not satisfy this last inequality merely because rho has a positive occupied-state gap. The theorem directly covers a full-column-rank factor representation; a redundant representation requires a separately charged reduction of its small base Gram.

For Hermitian estimators X_i, impose


```math
\mathrm{Tr}\rho X_i=0,\qquad
\mathrm{Tr}D_jX_i=\delta_{ij}.
\tag{5}
```


Their complex covariance is Z_ij=Tr rho X_i X_j. With a real positive semidefinite cost matrix W, define


```math
h_W(Z)=\mathrm{Tr}(W\mathrm{Re}Z)
+\|\sqrt W\mathrm{Im}Z\sqrt W\|_1.
\tag{6}
```


The Holevo value is the infimum of h_W over (5). The main quantitative theorem below takes W=I_p; a positive definite general W is covered explicitly at the end.

## 1. Objective facts

If Q>=0, then


```math
\|\mathrm{Im}(\sqrt WQ\sqrt W)\|_1
\le\mathrm{Tr}(W\mathrm{Re}Q).
\tag{7}
```


Indeed, for C=sqrt(W)Qsqrt(W)>=0, Im C=(C-overline C)/(2i), and the trace-norm triangle inequality bounds its norm by Tr C. Consequently:

1. h_W(Z+Q)>=h_W(Z) whenever Q>=0.
2. For W=I, Tr Z<=h_I(Z)<=2Tr Z when Z>=0.
3. On Hermitian matrices, |h_I(Z+E)-h_I(Z)|<=2||E||_1.
4. For real R and Z>=0, h_I(RZR^T)<=||R||_op² h_I(Z).

The fourth assertion bounds the real trace using Re Z>=0 and uses ||RAR^T||_1<=||R||_op²||A||_1 for the imaginary part.

For W=I, one also has the exact epigraph identity


```math
h_I(Z)=\min\{\mathrm{Tr}C:C\in\mathbb S^p,\ C\succeq Z\}.
\tag{8}
```


For the upper direction take C=Re Z+sqrt(-(Im Z)²), which is real symmetric. Its difference from Z is positive semidefinite. For the lower direction apply (7) to C-Z. This identity is used in the certified rational SDP.

## 2. Exact compression

**Theorem 1.** Let S be the column span of V and P its projector. Restricting the statistical model to S preserves both its Holevo value and its entire feasible real covariance region. In particular, an effective Hilbert space of dimension at most m=r(p+1) suffices.

**Proof.** Both rho and every D_i are supported in S. Thus Y_i=PX_iP preserves (5). Inserting P+(I-P) between X_i and X_j gives


```math
Z(X)=Z(Y)+Q,
\quad Q_{ij}=\mathrm{Tr}\rho X_i(I-P)X_j.
\tag{9}
```


Q is the positive semidefinite Gram matrix of the Hilbert–Schmidt vectors (I-P)X_i sqrt(rho). Equation (7) shows that compression cannot increase h_W. Conversely, a compressed tuple can be embedded into the original space. More strongly, any real covariance matrix dominating Z(X) also dominates Z(Y), and embedding proves the reverse inclusion of feasible regions. ∎

Any factor of K gives an equivalent compressed model. In particular, taking sqrt(K) as a column matrix works even when K is singular. No lower bound on the smallest positive tangent-Gram eigenvalue is involved. The finite classical implementation need not compute this square root: the rational weighted-factor construction in the companion note avoids it entirely.

## 3. Model-space stability

Let L_i be the canonical SLDs satisfying (rho L_i+L_i rho)/2=D_i, with their kernel–kernel blocks set to zero. Their Fisher matrix is J_ij=Re Tr rho L_iL_j. The SLD-dual tuple T_i=sum_k(J^-1)_ik L_k is feasible and has real covariance J^-1. Therefore


```math
C_H(I_p)\le2\mathrm{Tr}J^{-1}\le2p/j.
\tag{10}
```


An arbitrary feasible estimator may also have its state-kernel–kernel block deleted without changing its constraints or objective. Writing it as


```math
X_i=\begin{bmatrix}A_i&B_i^\dagger\\ B_i&0\end{bmatrix}
```


relative to supp(rho) and its complement yields


```math
\|X_i\|_F^2\le\frac2\nu\mathrm{Tr}\rho X_i^2.
\tag{11}
```


This supplies a bounded optimal representative. It also proves attainment: a minimizing sequence in this quotient can be restricted to a compact Frobenius ball.

**Theorem 2.** Suppose two normalized local models have occupied state weights at least nu_*>0 and Fisher matrices at least j_* I_p. Define


```math
\delta_\rho=\|\rho-\widehat\rho\|_1,
\qquad
\delta_D=\left(\sum_i\|D_i-\widehat D_i\|_F^2\right)^{1/2},
\quad C_{\max}=\frac{2p}{j_*},
\quad B=2\sqrt{\frac p{\nu_*j_*}}.
\tag{12}
```


If B delta_D<1, then


```math
\widehat C_H\le
\frac{C_H+2B^2\delta_\rho}{(1-B\delta_D)^2}.
\tag{13}
```


The same statement holds with the models exchanged. In particular, if B delta_D<=1/2,


```math
|\widehat C_H-C_H|
\le6C_{\max}B\delta_D+8B^2\delta_\rho.
\tag{14}
```


**Proof.** Choose an optimal original tuple with kernel–kernel blocks deleted. Equations (10)–(11) give sum_i||X_i||F²<=B². Its perturbed unbiasedness matrix M_ij=Tr Dhat_j X_i obeys


```math
\|M-I\|_{\mathrm{op}}\le\|M-I\|_F\le B\delta_D.
```


Let mu_i=Tr rhohat X_i. Centering Y_i=X_i-mu_i I preserves derivative constraints and subtracts the PSD covariance mu mu^T. Applying M^-1 in parameter space then repairs unbiasedness exactly.

For covariance perturbations define T|v>=sum_i X_i|v> tensor |i>. Its squared norm is at most B². Up to a harmless transpose, covariance is the partial trace of T rho T†. Trace-norm contraction on Hermitian inputs therefore gives


```math
\|\widehat Z(X)-Z(X)\|_1\le B^2\delta_\rho.
```


Use objective fact 3, the centering decrease, and objective fact 4 with ||M^-1||<=(1-B delta_D)^-1. This proves (13). Exchange the models and use (1-e)^-2<=1+6e for 0<=e<=1/2 to obtain (14). ∎

The two states need not have identical supports. They need only be embedded into a common space, as provided by factor alignment.

## 4. Factor distance and Fisher-condition transfer

For aligned normalized factor tuples, put


```math
\delta_\Psi=\|\Psi-\widehat\Psi\|_F,
\qquad
\delta_\Phi=\left(\sum_i\|\Phi_i-\widehat\Phi_i\|_F^2\right)^{1/2}.
```


If their derivative-factor norms are bounded by L and Lhat, respectively, product expansion gives


```math
\delta_\rho\le2\delta_\Psi,
\qquad
\delta_D\le2\delta_\Phi+(L+\widehat L)\delta_\Psi.
\tag{15}
```


These estimates use Frobenius product inequalities, not ambient-dimension factors. The optimal aligned joint distance sqrt(delta_Psi²+delta_Phi²) is exactly d_B(K,Khat).

**Lemma 3.** Suppose only the true model is initially promised to satisfy the occupied-state lower bound nu and J>=jI. If a valid reconstructed model satisfies


```math
\delta_\rho\le\nu/(2p),\qquad
\delta_D\le\sqrt{\nu j/(8p)},
\tag{16}
```


then its Fisher matrix is at least jI/8.

**Proof.** The true SLD-dual tuple has real covariance J^-1<=I/j and, after kernel deletion, squared Frobenius norm at most B_0²=2p/(nu j). Its perturbed constraint matrix lies within B_0 delta_D<=1/2 of I, so it can be repaired using an inverse of norm at most two. Its real covariance changes in operator norm by at most B_0² delta_rho<=1/j. Centering decreases covariance. The repaired tuple therefore has real covariance at most 8I/j. Its existence proves local identifiability of the reconstructed model; the SLD matrix Cramér–Rao inequality now gives Jhat^-1<=8I/j. ∎

This argument transfers a true-model promise. It does not treat a nominal floating-point Fisher eigenvalue as a certificate about an unknown input.

## 5. A finite-precision version with an interior physical fit

For computation, fit the learned Gram to the relaxed physical set


```math
\mathcal F=\left\{K\succeq0:
\mathrm{Tr}K_{00}=1,\quad K_{00}\succeq\frac\nu2 I_r,
\quad\mathrm{Re}\mathrm{Tr}K_{0i}=0,
\quad\sum_i\mathrm{Tr}K_{ii}\le2L^2\right\}.
\tag{17}
```


The true Gram is feasible. The reconstruction has occupied weights at least nu/2 and derivative-factor norm at most sqrt(2)L. This relaxation has a known strict interior point relative to the displayed affine equations. It avoids an unjustified generic-SDP invocation on a potentially thin feasible face. The fidelity-block regularization and exact fit certificate are given in the companion note.

Set


```math
b_L=2+(1+\sqrt2)L,
\qquad
A=\frac{768p^{3/2}b_L}{\sqrt\nu\,j^{3/2}}
+\frac{1024p}{\nu j},
\tag{18}
```


and


```math
\gamma_0=
\min\left\{\frac\nu{4p},
\frac{\sqrt{\nu j}}{16\sqrt p\,b_L}\right\}.
\tag{19}
```


**Theorem 4.** If the fitted Khat belongs to (17) and d_B(K,Khat)<=gamma<=gamma_0, then


```math
|C_H(\widehat K,I_p)-C_H(K,I_p)|\le A\gamma.
\tag{20}
```


**Proof.** Equation (15) gives delta_rho<=2gamma and delta_D<=b_L gamma. The margins (19) imply (16), hence Jhat>=jI/8. Apply Theorem 2 to both models with nu_*=nu/2 and j_*=j/8. Then Cmax=16p/j, B²=64p/(nu j), and B=8sqrt(p/(nu j)). The second margin in (19) ensures B delta_D<=1/2. Substituting into (14) gives exactly (18). ∎

The full tangent Gram may have arbitrarily small or zero eigenvalues. The required positive spectral promise concerns only the occupied state weights in K_00.

## 6. Acquisition, confidence interval, and resources

Assume a controlled circuit and its inverse prepare a coherent joint factor array with known normalization alpha, charging the circuit and derivative construction costs. The success/failure label state is


```math
\omega=(K^T/\alpha^2)\oplus
\left(1-\mathrm{Tr}K/\alpha^2\right),
\tag{21}
```


on d_label=m+1 dimensions. Failure is retained as a label, so there is no uncharged postselection. The relevant tomography rank can be m+1, even if the physical state has rank r.

For a target error epsilon, choose a positive dyadic gamma, within a constant factor of


```math
\min\{\gamma_0,\epsilon/(2A)\},
\qquad\zeta\le\gamma/(8\alpha).
\tag{22}
```


Learn omega to Bures distance zeta with confidence 1-delta. Pinch its estimate into good/failure blocks and define the preliminary Gram Kpre by transposing and multiplying the good block by alpha². Then d_B(K,Kpre)<=alpha zeta. The finite fit in the companion note returns a rational Khat in (17) with d_B(K,Khat)<=3gamma/4, using no exact projection assumption.

Let [ell,u] be a certified classical interval for C_H(Khat,I_p), with u-ell<=epsilon. Then


```math
[\ell-A\gamma,\ u+A\gamma]
\tag{23}
```


contains the true value with probability at least 1-delta. Its midpoint has error at most epsilon/2+A gamma<=epsilon. One may intersect the interval with nonnegative costs without weakening validity. The interval is conditional on the tomography confidence guarantee and the stated physical promises; the optimization certificates themselves are deterministic rational checks.

Using the coherent Bures learner, the number of joint-preparation or inverse calls is


```math
\widetilde O\!\left((m+1)^2\alpha
\max\left\{\frac A\epsilon,\frac1{\gamma_0}\right\}\right).
\tag{24}
```


This includes the conditioning thresholds rather than only the high-accuracy asymptotic. At fixed r,p,alpha,nu,j,L the precision dependence is inverse epsilon. Gate complexity additionally charges the joint circuit, its workspace, the learner's polynomial small-register overhead, and synthesis accuracy. The classical computation is polynomial in m,p and the finite input/certificate bit length; no physical D-dimensional state is classically reconstructed. Exact resource statements for that computation are in the companion note.

## 7. General weights and limits

For a positive definite W with eigenvalues w_min,w_max, Theorem 2 becomes


```math
C_{\max}=2\mathrm{Tr}W/j_*,\quad
B^2=2C_{\max}/(\nu_* w_{\min}),\quad
a=\sqrt{w_{\max}/w_{\min}}\,B\delta_D,
```



```math
\widehat C_H(W)\le
\frac{C_H(W)+2w_{\max}B^2\delta_\rho}{(1-a)^2}.
\tag{25}
```


For a<=1/2 the symmetric error is at most 6Cmax a+8wmax B²delta_rho. This follows by bounding the repair in the weighted coordinates sqrt(W)M^-1W^-1/2. The rational SDP simply replaces Tr C by Tr WC. A singular weight requires separate nuisance/estimand treatment or a charged regularization; it is not silently included in the quantitative positive-definite statement.

Two limits matter. First, ordinary entrywise Gram accuracy does not imply a linear Holevo-value error: the pure two-parameter family C_H(c)=4/(1+sqrt(1-c²)) has a square-root cusp in Gram entries at c=1 while J=I. Second, occupied-state conditioning cannot simply be discarded. The explicit bounded-generator rank-three obstruction, including the opaque-oracle access qualification, is recorded in [conditioning_limits.md](conditioning_limits.md). No unproved state-gap-free improvement is used in (24).
