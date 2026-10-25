# Bures tomography from a purification circuit of unrestricted width

## Scope and provenance

This note proves a coherent-query tomography theorem by combining three established ingredients: Frobenius tomography from a state-preparation unitary [1], amplitude estimation [2], and fixed-point amplitude amplification [3]. Its adaptive spectral decomposition and weighted error accounting are closely related to the copy-model construction of Flammia and O'Donnell [4]. The proof below is self-contained apart from the explicitly stated primitives. It does not treat the copy-model theorem as a black box whose access assumptions can be silently changed.

The substantive access issue is the purifying register. When a preparation circuit already uses a known purification space of dimension $`r`$, pure-state tomography in dimension $`dr`$ gives the same query scaling. Here the supplied purification may use a much larger register, and the algorithm never reconstructs that register or assumes a free conversion to a minimal purification. The result is a composition theorem; the ingredients and closely related adaptive geometry are inherited. A bounded literature audit is not a proof that no earlier equivalent theorem exists.

## 1. Model, conventions, and theorem

Let $`\mathcal H_S\cong\mathbb C^d`$ be a designated system register, with $`d\ge2`$. A unitary $`U`$ on $`a`$ qubits prepares


```math
U|0^a\rangle=|\Psi\rangle_{SE},\qquad
\rho=\mathrm{Tr}_E|\Psi\rangle\langle\Psi|,
\qquad \mathrm{rank}\rho\le r\le d.
\qquad\text{(1)}
```


The known integer $`r`$ is an upper bound, rather than necessarily the exact rank. The environment $`E`$ can have arbitrary dimension compatible with the circuit width $`a`$. A query is one application of $`U`$ or $`U^\dagger`$. The algorithm also uses ordinary gates, intermediate measurements, and classical computation. It does not require controlled $`U`$, access to matrix entries of $`\rho`$, a Schmidt basis, or QRAM. A non-power-of-two $`d`$ can be padded to the next power of two at a factor smaller than two in dimension.

For density matrices, use the root-fidelity convention


```math
F(\rho,\sigma)=\|\sqrt\rho\sqrt\sigma\|_1,\qquad
D_B(\rho,\sigma)=\sqrt{2-2F(\rho,\sigma)}.
\qquad\text{(2)}
```


For positive semidefinite, not necessarily normalized matrices, replace $`2`$ in (2) by $`\mathrm{Tr}\rho+\mathrm{Tr}\sigma`$.

### Theorem 1 — Arbitrary-purifier Bures tomography

For $`0<\epsilon\le1/4`$ and $`0<\delta\le1/3`$, there is an algorithm that, with probability at least $`1-\delta`$, returns a classical density matrix $`\widehat\rho`$ satisfying


```math
D_B(\rho,\widehat\rho)\le\epsilon
\qquad\text{(3)}
```


using


```math
Q=\widetilde O\!\left(\frac{dr}{\epsilon}\right)
\qquad\text{(4)}
```


queries to $`U,U^\dagger`$. The suppressed factors are polylogarithmic in $`d,r,1/\epsilon,1/\delta`$. The guarantee is uniform over the positive eigenvalues of $`\rho`$.

The ideal spectral output has rank at most


```math
\min\{d,\,O(r\log(1/\epsilon))\}.
\qquad\text{(5)}
```


A finite classical output can instead be an explicitly specified rational factor $`B`$, representing the exactly positive and normalized matrix $`BB^\dagger/\mathrm{Tr}(BB^\dagger)`$. It obeys the same rank bound. A dense matrix can be formed from that factor in polynomial time.

If $`G_U`$ denotes the gate cost of implementing a query, a conservative gate bound is


```math
\widetilde O\!\left[
\frac{dr}{\epsilon}\left(G_U+a+d^{5/2}\right)
\right].
\qquad\text{(6)}
```


This uses $`O(a)`$ clean work qubits for zero-state phases and includes explicit synthesis of small-system basis changes. An implementation of the invoked tomography primitive uses


```math
O\!\left(a+d^2\mathrm{polylog}
\frac{dr}{\epsilon\delta}\right)
\qquad\text{(7)}
```


qubits. Classical arithmetic, circuit generation, and output processing are polynomial in $`d,a,1/\epsilon,\log(1/\delta)`$. More precise inherited bounds appear below. The power of $`d`$ in (6) is a conservative implementation bound, not an asserted gate optimum.

For sufficiently small constant-bounded $`\epsilon`$, the query dependence in (4) is optimal up to logarithms: [1, Theorem 9.3] implies $`\Omega(dr/\epsilon)`$. This is an oracle statement. Reading a transparent circuit description can make particular instances easier; the theorem does not claim a classical computational lower bound for every circuit family.

## 2. Established primitives and their precise role

### 2.1 Frobenius tomography

The published SODA 2023 version of [1, Theorem 8.2], with Schatten parameter $`q=2`$, supplies a positive semidefinite estimate $`Z`$ of a rank-at-most-$`r`$ state $`\sigma`$, with


```math
\|Z-\sigma\|_F\le\eta
\qquad\text{(8)}
```


except with probability at most $`\beta`$. Its input is an $`a_V`$-qubit preparation unitary $`V`$ for any purification of $`\sigma`$, together with $`V^\dagger`$. The stated query bound is


```math
Q_F(d,r,\eta,\beta)
=O\!\left[
\frac{d\sqrt r}{\eta}
\log^{3/2}\!\left(\frac d{\beta\eta}\right)
\sqrt{\left\lceil\frac{\log(1/\beta)}d\right\rceil}
\right].
\qquad\text{(9)}
```


Replacing the last factor by $`O(\sqrt{1+\log(1/\beta)/d})`$ avoids any harmless convention concerning logarithms. We use the rank-dependent branch; a sometimes smaller rank-independent bound is unnecessary.

The additional two-qubit gate bound is


```math
G_F=O\!\left[
\frac{(da_V+d^{7/2})\sqrt r}{\eta}
\mathrm{polylog}\!\left(\frac d{\beta\eta}\right)
\right],
\qquad\text{(10)}
```


excluding the cost of its calls to $`V,V^\dagger`$. The input unitaries need not be controlled [1, Section 8]. The simultaneous-observable construction uses $`O(d^2)`$ coordinate registers of logarithmic precision, together with the preparation workspace; this gives the polynomial-width accounting in (7). The explicit matrix-unit observables have efficient hardcoded block encodings [1, Theorem 8.1], so this invocation contains no hidden QRAM.

We also permit a finite Hermitian representation of $`Z`$ with additional error $`O(\eta)`$. Positivity and the outer algorithm's rank bound are imposed by a separately analyzed classical truncation. The *input state* of the primitive, however, really has rank at most $`r`$; Lemma 2 is needed for that point.

### 2.2 Amplitude estimation

For a known system projector $`P`$, amplitude estimation [2, Theorem 12] estimates


```math
\tau=\langle\Psi|(P\otimes I_E)|\Psi\rangle
\qquad\text{(11)}
```


with error at most


```math
\frac{2\pi\sqrt{\tau(1-\tau)}}M+\frac{\pi^2}{M^2}
\qquad\text{(12)}
```


with probability at least $`8/\pi^2`$, using $`O(M)`$ queries. Repetition and a median give failure probability $`\beta`$ at an $`O(\log(1/\beta))`$ multiplicative cost. Clip the estimate to $`[0,1]`$, which cannot increase its error.

No controlled preparation is hidden here. A controlled phase about $`|\Psi\rangle`$ is


```math
(I_C\otimes U)
\bigl(|0\rangle\langle0|_C\otimes I+
|1\rangle\langle1|_C\otimes S_0\bigr)
(I_C\otimes U^\dagger),
\qquad\text{(13)}
```


where $`S_0`$ is a known phase on $`|0^a\rangle`$. Controlled Grover iterates therefore use uncontrolled $`U,U^\dagger`$, together with controlled known phases.

### 2.3 Fixed-point amplification

For a valid lower bound $`\tau\ge\tau_0>0`$, fixed-point amplitude amplification [3, Eqs. (1)–(3), (6)–(10)] constructs a unitary whose output is within vector distance $`\zeta`$, up to global phase, of


```math
|\phi_P\rangle=
\frac{(P\otimes I_E)|\Psi\rangle}{\sqrt\tau}
\qquad\text{(14)}
```


using


```math
O\!\left(\tau_0^{-1/2}\log(1/\zeta)\right)
\qquad\text{(15)}
```


queries. It is a fixed reversible circuit after $`\tau_0,\zeta`$ are chosen; reversing its gate sequence supplies its inverse.

For an explicit conversion, target failure amplitude at most $`\zeta/4`$, and choose the odd sequence length at least


```math
\frac{\log(8/\zeta)}{\sqrt{\tau_0}}.
\qquad\text{(16)}
```


The successful component is a scalar multiple of (14), so the failure-amplitude bound implies small vector error after global-phase alignment. Analytically prescribed phase angles may be rounded to total circuit error $`O(\zeta)`$. Rounding preserves the form of the phase operations, needs only logarithmic angle precision, and is absorbed by reducing the target failure amplitude by a constant factor.

### Lemma 2 — Exact rank preservation inside amplification

Any sequence of arbitrary phase rotations about $`|\Psi\rangle`$ and range$`(P\otimes I_E)`$, applied to $`|\Psi\rangle`$, has output


```math
|\psi_{\rm amp}\rangle
=\bigl((uP+v(I-P))\otimes I_E\bigr)|\Psi\rangle
\qquad\text{(17)}
```


for complex scalars $`u,v`$. Its reduced system state has rank at most $`r`$.

**Proof.** Both phase operations preserve the span of $`(P\otimes I_E)|\Psi\rangle`$ and $`((I-P)\otimes I_E)|\Psi\rangle`$. If $`\Psi`$ is the matrix of purification amplitudes, the output amplitude matrix is $`(uP+v(I-P))\Psi`$, whose rank is at most $`\mathrm{rank}\Psi=\mathrm{rank}\rho`$. This remains true for rounded phase angles. Clean phase-implementation ancillas are uncomputed. $`\square`$

This lemma concerns the ideal phase circuit, including finitely specified rounded phase angles. Arbitrary gate-synthesis errors can break exact rank. We do not apply a rank-$`r`$ theorem directly to that imperfect hardware output; Section 8 compares the complete compiled algorithm with the ideal algorithm by a global hybrid bound.

## 3. A subnormalized residual-estimation lemma

Fix $`c=1/100`$, $`0<h<1/100`$, and


```math
R=P\rho P,\qquad \tau=\mathrm{Tr}R.
\qquad\text{(18)}
```


### Lemma 3 — Residual estimator

A procedure using $`\widetilde O(dr/\sqrt h)`$ original preparation queries either returns **STOP**, certifying $`\tau\le3h`$, or returns $`\widehat\tau`$ and $`S`$ satisfying


```math
\begin{split}
&\tau\ge h,\qquad
|\widehat\tau-\tau|\le\frac c4\sqrt{\frac{\tau h}r},\\
&S=PSP,\quad S\succeq0,\quad \mathrm{rank}S\le r,\quad
\|S-R\|_F\le c\sqrt{\frac{\tau h}r},
\end{split}
\qquad\text{(19)}
```


except with its assigned failure probability. The construction includes finite-accuracy classical diagonalization; its cost is bounded even on unsuccessful statistical outcomes.

**Proof.** Estimate $`\tau`$ using (12) with


```math
M=\left\lceil C_A\sqrt{r/h}\right\rceil
\qquad\text{(20)}
```


for sufficiently large universal $`C_A`$, and amplify confidence. On success,


```math
|\widehat\tau-\tau|
\le\frac c{16}\left(\sqrt{\frac{\tau h}r}+\frac hr\right).
\qquad\text{(21)}
```


Stop when $`\widehat\tau\le2h`$. Equation (21) implies a stop has $`\tau\le3h`$. Conversely $`\tau\le h`$ implies $`\widehat\tau<2h`$, so a continuing round has $`\tau\ge h`$. Then $`h/r\le\sqrt{\tau h/r}`$, giving (19) and


```math
0.9975\tau\le\widehat\tau\le1.0025\tau.
\qquad\text{(22)}
```


Set


```math
\tau_0=\widehat\tau/2,\qquad
\eta=\frac c{100}\sqrt{\frac h{r\widehat\tau}},\qquad
\zeta=\eta/100.
\qquad\text{(23)}
```


The value $`\tau_0`$ is a valid constant-factor lower bound. Construct $`V_P`$ by fixed-point amplification to approximate (14) to vector error $`\zeta`$. Its actual ideal-circuit marginal $`\sigma_P`$ has rank at most $`r`$, by Lemma 2, and


```math
\|\sigma_P-R/\tau\|_F
\le\|\sigma_P-R/\tau\|_1\le2\zeta.
\qquad\text{(24)}
```


Apply the Frobenius primitive to $`\sigma_P`$ with accuracy $`\eta`$, and compress the Hermitian estimate to $`P`$. Compute a spectral representation $`W\mathrm{diag}(\lambda)W^\dagger`$ with additional Frobenius residual at most $`\eta`$. Here $`W`$ is an exactly specified product of system rotations acting inside range $`P`$; Section 8 supplies a gap-free finite-precision construction. Replace negative diagonal values by zero and retain at most the largest $`r`$ positive values. Call the resulting matrix $`T`$.

For a Hermitian matrix $`C`$, this truncation is a nearest PSD rank-at-most-$`r`$ matrix in Frobenius norm. Since $`R/\tau`$ belongs to that set,


```math
\|T-R/\tau\|_F
\le2(\eta+2\zeta+\eta)\le5\eta.
\qquad\text{(25)}
```


Set $`S=\widehat\tau T`$. Since $`\|R/\tau\|_F\le1`$,


```math
\|S-R\|_F
\le5\widehat\tau\eta+|\widehat\tau-\tau|
\le c\sqrt{\tau h/r}.
\qquad\text{(26)}
```


The final inequality has substantial slack with (22)–(23). Further deterministic numerical error can be included by a constant reduction in $`\eta`$.

The Frobenius primitive uses


```math
\widetilde O(d\sqrt r/\eta)
=\widetilde O\!\left(dr\sqrt{\widehat\tau/h}\right)
\qquad\text{(27)}
```


calls to $`V_P,V_P^\dagger`$, each costing $`\widetilde O(1/\sqrt{\widehat\tau})`$ original queries. Their product is $`\widetilde O(dr/\sqrt h)`$. The mass-estimation cost $`\widetilde O(\sqrt{r/h})`$ is smaller.

For worst-case accounting, clip every mass estimate to $`[0,1]`$ and invoke the conditional procedure only when $`\widehat\tau>2h`$. The costs remain bounded on every branch, whether or not (21) was correct. No assumption about the true $`\tau`$ is needed to bound executed gates. $`\square`$

## 4. The adaptive algorithm

Set


```math
L=\left\lceil8\log_2(8/\epsilon)\right\rceil,\qquad
h=\frac{\epsilon^2}{1000L}.
\qquad\text{(28)}
```


Allow at most $`L+1`$ residual-estimation calls. Assign failure probability


```math
\beta=\frac{\delta}{4(L+1)}
\qquad\text{(29)}
```


to each mass estimate and each Frobenius estimate. A known zero-dimensional residual terminates without an estimate.

Start with $`P_1=I`$. In round $`t`$, apply Lemma 3 to $`P_t\rho P_t`$. If it returns STOP, end the loop. Otherwise, use its spectral representation of $`S_t`$ and retain the eigenvectors whose recorded eigenvalues obey


```math
s_i\ge\theta_t:=\frac{\widehat\tau_t}{4r}.
\qquad\text{(30)}
```


Let $`H_t`$ be their projector. Record the retained values and vectors, and set $`P_{t+1}=P_t-H_t`$. All subsequent basis rotations act only inside range $`P_{t+1}`$; earlier retained vectors remain fixed. If the hard round cap is reached, or a malformed numerical result occurs, return an arbitrary fixed density matrix. Such an abort cannot occur on the simultaneous successful events proved below.

At termination define


```math
H=\sum_tH_t,\qquad
S_{\rm out}=\sum_{i\text{ retained}}s_i|v_i\rangle\langle v_i|,
\qquad \sigma=S_{\rm out}/\mathrm{Tr}S_{\rm out}.
\qquad\text{(31)}
```


If no vector was retained, return an arbitrary fixed state. The initial mass is one, so this case does not occur on successful events.

## 5. Progress and weighted error accumulation

### Lemma 4 — Each continuing round removes at least half the mass

On (19), $`1\le\mathrm{rank}H_t\le r`$,


```math
\mathrm{Tr}(H_t\rho)\ge\tau_t/2,\qquad
\tau_{t+1}\le\tau_t/2,
\qquad\text{(32)}
```


and each retained eigenvalue satisfies


```math
s_i\ge\tau_t/(5r).
\qquad\text{(33)}
```


**Proof.** Write $`R=R_t,S=S_t,e=\|R-S\|_F,\tau=\tau_t`$. Both matrices have rank at most $`r`$, so


```math
|\mathrm{Tr}(S-R)|\le\sqrt{2r}\,e.
\qquad\text{(34)}
```


There are at most $`r`$ positive eigenvalues of $`S`$. The total weight excluded by (30) is at most $`r\theta_t=\widehat\tau/4`$. Since $`H_t`$ has rank at most $`r`$,


```math
\begin{split}
\mathrm{Tr}(H_tR)
&\ge\mathrm{Tr}(H_tS)-\sqrt r\,e\\
&\ge\tau-(\sqrt2+1)\sqrt r\,e-\widehat\tau/4\\
&\ge\tau/2.
\end{split}
\qquad\text{(35)}
```


The last step uses $`\tau\ge h`$, $`e\le.01\sqrt{\tau h/r}`$, and (22). The lower bound actually exceeds $`0.72\tau`$. It proves nonemptiness and progress. Equation (33) follows from (22) and the threshold. $`\square`$

Mass halves at every continuing round, so there are at most $`\lceil\log_2(1/h)\rceil+1\le L`$ continuing rounds. The last inequality follows from (28) for $`0<\epsilon\le1/4`$. Thus the cap is safe. At STOP, the true discarded mass is at most $`3h`$.

For a positive diagonal $`S=\mathrm{diag}(s_i)`$, define


```math
g_S(D)=\mathrm{Tr}[D\mathcal J_S^{-1}(D)]
=\sum_{i,j}\frac{2|D_{ij}|^2}{s_i+s_j},\qquad
\mathcal J_S(D)=\frac{SD+DS}{2}.
\qquad\text{(36)}
```


For normalized density matrices, $`g_S(X-S)`$ is the Bures $`\chi^2`$-divergence and dominates $`D_B(X,S)^2`$ [4, Definition 2.26 and Proposition 2.31]. Only the strictly positive retained space is used in (36).

### Lemma 5 — The retained weighted error

In the recorded retained basis, put $`X=H\rho H`$ and $`S=S_{\rm out}`$. Then


```math
g_S(X-S)\le10c^2hL.
\qquad\text{(37)}
```


**Proof.** Assign an ordered matrix entry $`(i,j)`$ to the earliest round $`t`$ at which either coordinate is retained. One of its recorded eigenvalues satisfies (33), giving


```math
\frac2{s_i+s_j}\le\frac{10r}{\tau_t}.
\qquad\text{(38)}
```


At round $`t`$, $`S_t`$ is diagonal in the selected basis. The squared errors of all entries assigned to the round sum to at most $`\|R_t-S_t\|_F^2`$. To justify this after later basis changes, partition the current matrix into retained and residual blocks. Subsequent rotations are the identity on the retained block and unitary on the residual block. They preserve the Frobenius norm of the retained-to-residual cross block and the retained block itself. Omitting coordinates eventually discarded can only decrease the sum. Ordered pairs are counted once, by their earliest round; both orientations of a cross entry are already included in the Frobenius norm.

Consequently,


```math
g_S(X-S)
\le\sum_t\frac{10r}{\tau_t}\|R_t-S_t\|_F^2
\le10c^2hL.
\qquad\text{(39)}
```


No comparison of exact and approximate spectral projectors occurs, so no eigengap enters. No coordinatewise relative estimate of the true diagonal probabilities is needed: the denominators are recorded positive estimated eigenvalues. $`\square`$

## 6. Normalization and the final Bures guarantee

### Lemma 6 — Exact normalization identity

For PSD $`X`$, positive $`S`$ on the same space, and $`p=\mathrm{Tr}X>0,q=\mathrm{Tr}S>0`$,


```math
\chi_B^2(X/p\|S/q)
=\frac q{p^2}
\left[g_S(X-S)-\frac{(p-q)^2}{q}\right].
\qquad\text{(40)}
```


Also,


```math
(p-q)^2\le q\,g_S(X-S).
\qquad\text{(41)}
```


**Proof.** The inner product in (36) satisfies


```math
\langle D,S\rangle_S=\mathrm{Tr}D,\qquad
g_S(S)=q,\qquad
\mathcal J_{S/q}^{-1}=q\mathcal J_S^{-1}.
\qquad\text{(42)}
```


Expand the square of $`X-(p/q)S=(X-S)-((p-q)/q)S`$, then multiply by $`q/p^2`$. This proves (40). Cauchy–Schwarz gives (41). $`\square`$

Here $`p\ge1-3h\ge1/2`$. Put $`D=g_S(X-S)`$. Equations (37) and (28) give $`D\le1/16`$. From (41),


```math
\sqrt q\le\frac{\sqrt D+\sqrt{D+4p}}2,
\qquad\text{(43)}
```


so $`q\le2`$, and


```math
D_B(X/p,S/q)^2\le\chi_B^2(X/p\|S/q)
\le8D\le80c^2hL.
\qquad\text{(44)}
```


The normalized projected purification $`(H\otimes I)|\Psi\rangle/\sqrt p`$ has overlap $`\sqrt p`$ with $`|\Psi\rangle`$. Therefore


```math
D_B(\rho,X/p)^2
\le2(1-\sqrt p)\le2(1-p)\le6h.
\qquad\text{(45)}
```


The triangle inequality, (28), and $`c=1/100`$ give


```math
D_B(\rho,\sigma)
\le\sqrt{6h}+\sqrt{80c^2hL}
<\epsilon/8.
\qquad\text{(46)}
```


The slack covers finite classical output representation. Condition on the preceding classical history at each round. The probability of any failed mass or Frobenius estimate is at most $`2(L+1)\beta\le\delta/2`$. The projectors are adaptive, but every primitive bound holds for every fixed projector, so this conditional union bound is valid.

## 7. Query, gate, width, and classical costs

### 7.1 Queries, with inherited logarithms visible

Let


```math
\Lambda=\log\!\left(
\frac{Cdr(L+1)^2}{\epsilon\delta}
\right)
\qquad\text{(47)}
```


for sufficiently large universal $`C>1`$. Since each continuing $`\widehat\tau\in(2h,1]`$, the logarithms in (9), (15), and the minimum $`\eta`$ in (23) are $`O(\Lambda)`$. A conservative explicit bound is


```math
\begin{split}
Q={}&O\!\left[
\frac{dr L^{3/2}}{\epsilon}\,
\Lambda^{5/2}
\sqrt{1+\frac{\log((L+1)/\delta)}d}
\right]\\
&+O\!\left[
\frac{\sqrt r L^{3/2}}{\epsilon}
\log\frac{L+1}{\delta}
\right].
\end{split}
\qquad\text{(48)}
```


The second line is mass estimation. This implies (4). Formula (48) is not intended as a log-optimal result. The largest amplification length is $`\widetilde O(1/\sqrt h)=\widetilde O(1/\epsilon)`$, and is charged in (48).

### 7.2 Gates and no-QRAM construction

Represent a known system basis by an explicit sequence of two-level rotations and phases, using $`\widetilde O(d^2)`$ elementary gates [5]. The running product can have $`O(L)`$ basis changes; this affects only logarithmic factors. Coordinate-interval projectors and zero-state phases use explicit reversible control logic. A phase about the $`a`$-qubit all-zero state costs $`O(a)`$ gates with $`O(a)`$ clean ancillary qubits [3, Fig. 2].

Work in the current coordinate frame: apply the running basis inverse after $`U`$, so the residual projector is a coordinate interval. This preparation unitary has width $`O(a)`$ and gate cost $`\widetilde O(G_U+a+d^2)`$. Conditional preparation costs


```math
\widetilde O\!\left[
\widehat\tau^{-1/2}(G_U+a+d^2)
\right].
\qquad\text{(49)}
```


Multiply (49) by (27), add (10), and sum over the rounds. This yields (6). The environment's Hilbert dimension never appears; its qubits are touched through $`U,U^\dagger`$ and the charged zero-state phase.

The width bound (7) follows from the primitive's coordinate registers, $`O(a)`$ preparation and phase workspace, and logarithmically many precision qubits. Repetitions and rounds reuse workspace. No classical $`d\times d`$ matrix is stored in QRAM: measured estimates are processed classically and the resulting basis rotations compiled into a circuit.

### 7.3 Classical computation and output size

The outer algorithm performs $`O(L)`$ matrix diagonalizations, truncations, and basis updates. Dense implementations need $`\widetilde O(d^3L)`$ arithmetic operations, in addition to the invoked primitive's classical postprocessing and the classical control work for executing the quantum circuits. For the maximal-pivot Jacobi construction in Section 8, maintain the off-diagonal squared magnitudes in a max-heap: one rotation changes $`O(d)`$ entries, which can be updated in $`\widetilde O(d)`$ amortized time. Lazy versioned entries and periodic rebuilding keep the heap size $`O(d^2)`$. A naive full pivot search in every rotation instead gives a larger $`O(d^4)`$-type bound and is not the implementation used for this arithmetic estimate. Section 8 gives logarithmic precision without an unknown eigengap.

A factor output uses $`O(d\min\{d,rL\})`$ complex numbers. Forming a dense matrix takes polynomial time and writes $`d^2`$ entries. These costs remain real when the query count is smaller than a particular classical implementation's work. The optimality statement concerns queries, not all resources simultaneously.

## 8. Finite precision and implementation robustness

Three different approximation issues must be kept separate.

### 8.1 Classical spectral computations without a gap

The algorithm needs small residual $`\|C-W\mathrm{diag}(\lambda)W^\dagger\|_F`$, rather than individually accurate eigenvectors. Nearly degenerate eigenvalues therefore create no singular accuracy requirement.

A constructive procedure is maximal-pivot complex Jacobi diagonalization. Diagonalize a selected $`2\times2`$ Hermitian principal block by a two-level unitary. If $`E_{\rm off}`$ is the squared Frobenius norm of the off-diagonal entries, an exact step removes $`2|C_{ij}|^2`$. A maximal pivot in a $`k\times k`$ block gives


```math
E_{\rm off,new}
\le\left(1-\frac2{k(k-1)}\right)E_{\rm off}.
\qquad\text{(50)}
```


Thus $`N=O(k^2\log(\|C\|_F/\xi))`$ rotations suffice for residual $`\xi`$. A finite-precision implementation follows from an elementary error recurrence. Put $`q_k=\sqrt{1-2/[k(k-1)]}`$, for $`k>1`$. Choose the sum of the rotation and matrix-update errors in each step to be at most


```math
\gamma_{\rm step}\le
\min\left\{\frac{\xi}{20N},
\frac{(1-q_k)\xi}{20}\right\}.
```


Then the computed off-diagonal Frobenius norm satisfies
$`e_{j+1}\le q_ke_j+\gamma_{\rm step}`$, while the backward error between the computed matrix and the exactly represented product of rounded rotations is at most $`N\gamma_{\rm step}`$. The final spectral residual is therefore $`O(\xi)`$. A $`1\times1`$ block requires no rotation. A $`2\times2`$ block whose distance from a scalar matrix is below $`\gamma_{\rm step}`$ may be skipped; otherwise its rotation angles can be computed to the required absolute accuracy using a number of bits polynomial in $`\log(1/\gamma_{\rm step})`$. This avoids an unpromised eigenvector-condition assumption. The overall bit precision is $`O(\log(k\|C\|_F/\xi))`$, up to logarithms for the rotation count.

For bounded resource use on unsuccessful branches, clip each real and imaginary component of the primitive's reported matrix to $`[-2,2]`$, then Hermitianize. This gives Frobenius norm at most $`3d`$. Every density-matrix entry lies in the clipping box; coordinate clipping and Hermitianization cannot increase Frobenius error relative to that density matrix. In the current coordinate frame the residual compression is a principal submatrix, so this preprocessing is also valid for the residual calculation. Numerical routines may use fixed iteration budgets from the preceding residual bound.

In Lemma 3, use $`\xi=O(\eta)`$ and reserve this error in (25). Store the basis as an exact product of rotations with finitely specified angles. Choose nonnegative diagonal weights as rational numbers. Truncation acts on that represented diagonal matrix and produces an exactly PSD matrix of rank at most $`r`$, even when a computed value is near zero or near the selection threshold.

### 8.2 Exact nesting of ideal residual spaces

Keep a running product of the rotations. Every new factor acts only on the current residual coordinate interval and is the identity on previously retained coordinates. This defines exactly nested ideal projectors, even with rounded angles. Independently rounding a dense projector matrix and treating it as an exact projector would invalidate the proof and is not the prescribed implementation.

Classical arithmetic evaluates this specified product to the precision needed for subsequent computations; its error is included in the decomposition residual. The physical circuit is allowed a separate synthesis error, handled in Section 8.4. Resynthesis of a running basis is allowed if its error is charged there rather than silently redefining the ideal retained subspaces.

### 8.3 Exactly positive, normalized finite output

Let $`B_0`$ have columns $`\sqrt{s_i/q}\,v_i`$, so $`\sigma=B_0B_0^\dagger`$ and $`\|B_0\|_F=1`$. Compute a complex rational matrix $`\widetilde B`$ with


```math
\|\widetilde B-B_0\|_F\le\kappa,\qquad
\kappa=\epsilon/8.
\qquad\text{(51)}
```


Entrywise precision $`O(\kappa/\sqrt{dk})`$, where $`k`$ is the number of retained columns, suffices. It uses $`O(\log(dk/\epsilon))`$ bits per component. Return


```math
\widehat\rho=
\frac{\widetilde B\widetilde B^\dagger}
{\mathrm{Tr}(\widetilde B\widetilde B^\dagger)}.
\qquad\text{(52)}
```


This is exactly PSD and trace one when represented by its rational numerator and denominator. Since


```math
\left\|\frac{\widetilde B}{\|\widetilde B\|_F}-B_0\right\|_F
\le2\kappa,
\qquad\text{(53)}
```


the common-factor form of Uhlmann's bound gives $`D_B(\widehat\rho,\sigma)\le2\kappa=\epsilon/4`$. Together with (46), the finite output obeys (3), with room to spare. Approximating factors avoids reintroducing a square-root precision loss through a generic trace-to-Bures conversion.

### 8.4 Compiling the adaptive quantum circuit

First analyze the ideal adaptive algorithm with finitely specified phase angles and nested rotations. Every classical branch has a hard round cap, clipped mass estimates, and bounded subroutine lengths. There is therefore a deterministic worst-case bound $`N_{\rm ideal}`$ on its non-oracle elementary gates and oracle calls.

Compile each elementary gate to operator error at most $`\delta/(16N_{\rm ideal})`$. Finite universal gate synthesis incurs only a polylogarithmic overhead [6]. The elementary bound between unitary operator error and channel diamond error, followed by a telescoping hybrid argument, makes the total variation distance between the final classical output distributions at most $`\delta/2`$. The argument applies to adaptive measurements: retain classical records as registers, or condition on histories and use the uniform worst-case bound. This adds at most $`\delta/2`$ to the ideal failure probability $`\delta/2`$.

Compiled intermediate states may have rank larger than $`r`$. The rank-dependent analysis belongs to the ideal algorithm, and the complete compiled experiment is compared to that experiment. Applying the rank-dependent primitive directly to each noisy intermediate state would be an unjustified alternative argument.

If the supplied preparation is itself synthesized rather than treated as an exact oracle, include its errors in the same global budget and evaluate $`G_U`$ at the required synthesis precision. Alternatively, if the scientific target is the state of a specified finite circuit, that circuit defines the exact input. The theorem does not cover uncontrolled constant-strength noise on every query without fault-tolerant suppression.

## 9. Lower bound and inherited baselines

Published [1, Theorem 9.3] states that, for $`0<\xi\le1/128`$, trace-norm accuracy $`\xi/128`$, success probability at least $`2/3`$, and rank bound $`r`$, purified-unitary tomography requires $`\Omega(dr/\xi)`$ queries. Because


```math
\|\rho-\widehat\rho\|_1\le2D_B(\rho,\widehat\rho),
\qquad\text{(54)}
```


the Bures task has the same lower scaling: for example set $`\xi=256\epsilon`$ when $`\epsilon\le1/32768`$. Constants do not affect the asymptotic statement. The trivial $`d=1`$ case is excluded. The published proof's finitely many small-dimensional edge cases are covered by the ordinary two-state hybrid lower bound.

Three comparisons delimit the result:

1. Trace tomography at accuracy $`O(\epsilon^2)`$ gives Bures accuracy $`\epsilon`$, but uses $`\widetilde O(dr/\epsilon^2)`$ queries. The adaptive variance-sensitive construction removes this precision loss.
2. A known $`r`$-dimensional purification space permits pure-state tomography in dimension $`dr`$, already giving $`\widetilde O(dr/\epsilon)`$. That is the inherited baseline when this stronger representation is supplied.
3. The copy-model spectral-peeling and Bures/$`\chi^2`$ geometry of [4] is a close mathematical predecessor. Here conditional coherent Frobenius tomography replaces its statistical estimates, conditional preparation and inversion are charged, and (40) avoids separate diagonal refinement. Spectral peeling and Bures error accounting are not claimed as new ingredients.

## 10. Reproducibility scope

The accompanying verifier, [bures_verify.py](../verification/bures_verify.py), tests deterministic inequalities on finite matrices: rank preservation under generalized phase sequences, mass removal, weighted error, normalization identity (40), and factor rounding (52). It includes rank deficiency and eigenvalues spanning many scales. These checks can reveal algebraic or implementation mistakes; they do not prove the query theorem or simulate a fault-tolerant implementation of the Frobenius primitive. Those claims rest on the analytic argument.

## References

[1] J. van Apeldoorn, A. Cornelissen, A. Gilyén, and G. Nannicini, **Quantum tomography using state-preparation unitaries**, Proceedings of SODA 2023, pp. 1265–1318. Theorem 8.2 supplies Frobenius tomography; Theorem 9.3 supplies the trace-norm lower bound. Use the published version: the arXiv v1 lower bound is weaker.

https://doi.org/10.1137/1.9781611977554.ch47  
https://real.mtak.hu/191396/1/1.9781611977554.ch47.pdf

[2] G. Brassard, P. Høyer, M. Mosca, and A. Tapp, **Quantum Amplitude Amplification and Estimation**, Contemporary Mathematics 305, 53–74 (2002), arXiv:quant-ph/0005055. Theorem 12.

https://arxiv.org/abs/quant-ph/0005055

[3] T. J. Yoder, G. H. Low, and I. L. Chuang, **Fixed-Point Quantum Search with an Optimal Number of Queries**, Physical Review Letters 113, 210501 (2014).

https://doi.org/10.1103/PhysRevLett.113.210501  
https://arxiv.org/abs/1409.3305

[4] S. T. Flammia and R. O'Donnell, **Quantum chi-squared tomography and mutual information testing**, Quantum 8, 1381 (2024). Definition 2.26, Proposition 2.31, Theorem 3.18, and Corollary 3.20.

https://quantum-journal.org/papers/q-2024-06-20-1381/  
https://arxiv.org/html/2305.18519v2

[5] A. Barenco et al., **Elementary gates for quantum computation**, Physical Review A 52, 3457 (1995).

https://doi.org/10.1103/PhysRevA.52.3457  
https://arxiv.org/abs/quant-ph/9503016

[6] C. M. Dawson and M. A. Nielsen, **The Solovay–Kitaev algorithm**, Quantum Information & Computation 6, 81–95 (2006).

https://arxiv.org/abs/quant-ph/0505030
