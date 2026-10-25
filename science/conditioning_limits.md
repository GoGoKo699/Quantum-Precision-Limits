# State weights and a rare-branch obstruction

The [certified algorithm](circuit_access.md) assumes a lower bound on occupied
state weights. This page shows why that assumption cannot simply be removed,
even for fixed-rank states, bounded known generators, and well-conditioned
Fisher information. The example also distinguishes the information in a
global preparation oracle from access to individual rare-branch rotations.

The [compression and stability account](holevo_compression_and_stability.md)
contains the separate square-root cusp in ordinary Gram-entry error. Here the
input factors themselves become close, while the Holevo value changes too fast
for a uniform Lipschitz estimate.

## 1. A rank-three, two-parameter model

Use a flag qubit and a sensor qubit. For fixed $`0<a<1`$ and small $`0<q<1`$,
put $`t=1-2a\sqrt q`$ and

```math
\rho_a(q)=(1-q)|00\rangle\langle00|
 +q(1-a\sqrt q)|10\rangle\langle10|
 +aq^{3/2}|11\rangle\langle11|.
```

Its occupied eigenvalues are $`1-q`$, $`q-aq^{3/2}`$, and $`aq^{3/2}`$.
Take the fixed Hermitian generators

```math
G_1=I_f\otimes X_s/2,
\qquad G_2=Z_f\otimes Y_s/2.
```

Both have operator norm $`1/2`$. For example, encode the parameters with
$`e^{-i\theta_2G_2}e^{-i\theta_1G_1}`$ and evaluate at the origin. Both rotate
the sensor, with the orientation of the second rotation reversed in the rare
flagged branch.

The SLD Fisher matrix is

```math
J_a(q)=[(1-q)+qt^2]I_2
 =[1-4aq^{3/2}+4a^2q^2]I_2.
```

For $`a\in[1/4,3/4]`$ and sufficiently small $`q`$, it is bounded below by
$`I_2/2`$. Identifiability and the generator norms remain well behaved while
the smallest occupied state weight vanishes.

## 2. Exact Holevo value

For $`W=I_2`$, the value is

```math
C_H(a,q)=\frac{2}{\displaystyle\min_{-1\le\beta\le1}F_{a,q}(\beta)},
\qquad
F_{a,q}(\beta)=\frac{1-q}{1+\beta}
 +\frac{qt^2}{1-\beta t}.
```

At a zero denominator, use the corresponding extended-value limit. This
formula can be established directly, without assuming a symmetric form for
the optimal estimators.

The state and derivatives are block diagonal in the flag. Removing cross-flag
estimator blocks preserves the affine constraints and subtracts a positive
semidefinite covariance Gram. Within a sensor block, the state is diagonal and
its derivatives are off diagonal. Removing the estimator's diagonal entries
again subtracts a real positive semidefinite Gram. These reductions preserve
the optimum by the covariance monotonicity proved in
[compression](holevo_compression_and_stability.md).

Write the remaining estimator in branch $`f`$ as
$`X_i=a_{if}X+b_{if}Y`$. The branch weights are
$`w_0=1-q,w_1=q`$, Bloch lengths are $`r_0=1,r_1=t`$, and orientation signs are
$`s_0=1,s_1=-1`$. Local unbiasedness is

```math
\sum_f w_fr_f(-b_{1f})=1,
\quad \sum_f w_fs_fr_f a_{1f}=0,
\quad \sum_f w_fr_f(-b_{2f})=0,
\quad \sum_f w_fs_fr_f a_{2f}=1.
```

The Holevo objective is the maximum over $`\beta\in[-1,1]`$ of

```math
q_\beta(X)=\sum_f w_f\left[
 a_{1f}^2+b_{1f}^2+a_{2f}^2+b_{2f}^2
 +2\beta r_f(a_{1f}b_{2f}-b_{1f}a_{2f})\right].
```

Each quadratic form is positive semidefinite. The pair $`(a_1,b_2)`$ has
homogeneous constraints and can be set to zero. For the remaining pair, put
$`u_f=-b_{1f}`$ and $`v_f=s_fa_{2f}`$. Both obey
$`\sum_f w_fr_fu_f=\sum_f w_fr_fv_f=1`$. Their difference has a homogeneous
constraint and nonnegative cost. The sum reduces to weighted least squares,
giving

```math
\min_Xq_\beta(X)
 =\frac{2}{\displaystyle\sum_f w_fr_f^2/(1+\beta s_fr_f)}.
```

Convex minimax, with the compact scalar interval and coercivity of the maximum
objective on these coordinates, gives the formula for $`C_H`$. One can also
check it constructively: at an interior optimal $`\beta`$, the least-squares
tuple has zero imaginary covariance and attains its quadratic lower bound.
At an endpoint the imaginary covariance has the maximizing sign.

## 3. The singular scaling

Set $`A=\sqrt{1-q}`$ and $`B=\sqrt q\,t^{3/2}`$. The unconstrained minimizer is

```math
\beta_* =\frac{A-B}{At+B}.
```

It satisfies $`\beta_*\le1`$ exactly when $`A(1-t)\le2B`$. For any fixed
$`a<1`$, this holds for sufficiently small $`q`$. In this regime,

```math
C_H(a,q)=\frac{2(1+t)}{t[\sqrt{1-q}+\sqrt{qt}]^2}
 =4-(8-4a)\sqrt q+O(q).
```

In particular, for $`a_-=1/4`$ and $`a_+=3/4`$,

```math
C_H(a_+,q)-C_H(a_-,q)=2\sqrt q+O(q).
```

A common-basis purification is

```math
|\psi_a\rangle=\sqrt{1-q}|00\rangle|0\rangle_E
 +\sqrt{q(1-a\sqrt q)}|10\rangle|1\rangle_E
 +\sqrt a\,q^{3/4}|11\rangle|2\rangle_E.
```

For distinct fixed positive $`a,b`$,

```math
\|\psi_a-\psi_b\|_2^2
 =q[\sqrt{1-a\sqrt q}-\sqrt{1-b\sqrt q}]^2
  +q^{3/2}(\sqrt a-\sqrt b)^2
 =\Theta(q^{3/2}).
```

The derivative factors are $`\Phi_i=-iG_i\Psi`$. Since $`G_i^2=I/4`$, their
pairwise factor differences have half the norm of the corresponding state
factor difference. Thus these displayed joint factors have distance $`\Theta(q^{3/4})`$,
whereas the Holevo-value difference is $`\Theta(q^{1/2})`$. They give an
$`O(q^{3/4})`$ upper bound on the minimized Gram Bures distance; the obstruction
needs only this upper bound.

Consequently there is no uniform state-weight-independent Lipschitz bound in
factor/Bures distance under only fixed rank, bounded generators, and a Fisher
lower bound. Any uniform Hölder exponent for this class is at most $`2/3`$.
This is an obstruction; a matching general Hölder upper bound is not claimed.

## 4. A valid opaque-oracle lower bound

Closeness of output states alone does not make two preparation **unitaries**
hard to distinguish. A separately callable rotation inside a rare branch can
reveal its angle on a deliberately chosen input. The query lower bound must
therefore specify close global oracle extensions.

Let $`|e\rangle=|00\rangle_S|0\rangle_E`$ and define

```math
|u_a\rangle=
 \frac{|e\rangle+|\psi_a\rangle}{\sqrt{2+2\sqrt{1-q}}},
\qquad U_a=2|u_a\rangle\langle u_a|-I.
```

This Hermitian unitary prepares $`U_a|e\rangle=|\psi_a\rangle`$, and its
normalization is independent of $`a`$. Therefore

```math
\|U_a-U_b\|_{\mathrm{op}}
 \le4\|u_a-u_b\|_2=O(q^{3/4}).
```

Controlled and inverse calls obey the same difference bound. Known generator
operations and normalized derivative preparation built from
$`(-2iG_i)U_a`$ preserve it up to constants. The usual query hybrid argument
bounds the output distinguishability of a $`T`$-query procedure by
$`O(Tq^{3/4})`$. Constant-success discrimination needs
$`\Omega(q^{-3/4})`$ queries.

Choose $`q`$ proportional to $`\varepsilon^2`$, with a fixed constant that
separates the two values by more than $`2\varepsilon`$. Any constant-success
additive-$`\varepsilon`$ estimator of this Holevo value then needs

```math
T=\Omega(\varepsilon^{-3/2}).
```

This holds for rank three, two parameters, generator norms $`1/2`$, and
uniformly conditioned Fisher matrices. It rules out a uniform
$`O(1/\varepsilon)`$ global preparation-query theorem over the whole
state-weight-free class.

The lower bound is for opaque global preparation oracles. If $`a,q`$, and the
mixture are supplied as classical data, the displayed elementary formula
evaluates the toy model directly. Access to rare-branch rotations as separate
oracles is also stronger than the model used in this proof.

## 5. Consequence for useful precision computation

Three distinct regularity questions must be kept separate: occupied state
weights, local Fisher identifiability, and the spectrum of the joint tangent
Gram. The algorithm controls the first two and permits the third to be
singular. The example above keeps Fisher identifiability and generator norms
fixed while closing the occupied-state gap.

Approximate low rank can still be useful when a discarded component shares the
same bounded parameter encoding as the retained component. The
[two-sided truncation certificate](applications_and_limits.md) quantifies that
case and its sharp square-root sensitivity. A small discarded trace weight
without control of the parameter derivatives does not establish an accurate
precision benchmark.
