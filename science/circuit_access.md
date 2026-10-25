# Circuit access and the certified algorithm

The algorithm learns a small label register that contains the overlaps needed
for the local precision problem. The physical sensor is traced out. Its large
Hilbert space therefore contributes through the cost of its preparation circuit,
without becoming the dimension of the tomography output.

Use the [local model and notation](model.md). The
[compression theorem](holevo_compression_and_stability.md) shows that the joint
Gram matrix $`K=V^\dagger V`$, with $`V=[\Psi,\Phi_1,\ldots,\Phi_p]`$, determines
the Holevo value. This section states the access needed to acquire it and gives
the complete error allocation.

## 1. The coherent input

A supplied circuit $`U_V`$ and its inverse prepare a pure state whose successful
component is

```math
\frac1\alpha\sum_{i=0}^{p}\sum_{a=1}^{r}
 |i,a\rangle_L\,|V_{i,a}\rangle_S\,|\mathrm{good}\rangle_F.
```

Here $`V_0=\Psi`$, $`V_i=\Phi_i`$ for $`i>0`$, and the vectors in the last
display are **unnormalized** columns. The known normalization satisfies
$`\alpha^2\ge\operatorname{Tr}K`$. The remainder of the prepared state is
orthogonal in the success flag. All garbage used by the circuit remains part of
the coherent workspace until it is explicitly traced out for analysis.

The factor index $`a`$ has exactly $`r`$ values and $`\Psi`$ has full column rank.
A purification with redundant environment columns is a different input:
$`K_{00}\succeq\nu I_r`$ would then be false. The present application does not
quietly identify a physical rank promise with access to a minimal purification.
In the circuit families below the small environment is explicitly constructed.

### Constructing the input for unitary sensing

Suppose a known circuit prepares
$`\sum_a|a\rangle|\Psi_a\rangle`$, and a short, hardcoded linear combination of
known unitaries gives each Hermitian generator $`G_i`$ with normalization
$`\beta_i\ge\|G_i\|`$. Its block encoding has successful block $`G_i/\beta_i`$.
For $`\Phi_i=-iG_i\Psi`$:

1. Prepare a branch register with amplitudes proportional to
   $`1,\beta_1,\ldots,\beta_p`$.
2. Prepare the purification once, keeping its environment index.
3. In branch $`i>0`$, apply the corresponding controlled generator block encoding
   and the phase $`-i`$. Branch zero acts as the identity and is successful.

The resulting successful component is the one above with

```math
\alpha=\sqrt{1+\sum_{i=1}^{p}\beta_i^2}.
```

The branch preparation, SELECT operations, generator ancillas, controls, and
purification circuit are all charged to $`U_V`$. For an explicit list, a
conservative gate count is the preparation cost plus the sum of the controlled
generator costs and polynomial branch bookkeeping. No QRAM lookup or free
loading of a classical state vector is required. If the generator list is long,
or its coefficients have large absolute sum, that cost remains in the bound.

For a parameterized preparation circuit, derivative factors can instead be
assembled from its differentiated gates. The number and normalization of those
terms must be included explicitly; the theorem does not promise that every
circuit derivative is cheap.

## 2. A failure symbol prevents free postselection

Let $`m=r(p+1)`$. On a failure branch, reversibly move the old label into work
storage and replace it by one new label $`\perp`$. Copy the success flag into a
traced workspace register. Successful and failed components then have
orthogonal discarded flags. The reduced label state is exactly

```math
\omega=\frac{K^{\mathsf T}}{\alpha^2}
 \oplus\left(1-\frac{\operatorname{Tr}K}{\alpha^2}\right),
\qquad \dim L=m+1.
```

The transpose follows directly from taking the partial trace of
$`\sum_j|j\rangle|V_j\rangle`$: the label entry $`(j,k)`$ is
$`\langle V_k|V_j\rangle`$. It does not affect Bures distance.

The label state can have rank $`m+1`$. Applying the tomography theorem with the
physical rank $`r`$ here would undercount the resource. The label cleanup and its
inverse are part of the preparation circuit. There is no conditional-success
oracle supplied at unit cost.

## 3. From a confidence ball to a physical model

Apply [coherent Bures tomography](bures_tomography.md) to $`\omega`$ with Bures
error at most $`\zeta`$ and failure probability at most $`\delta`$. Pinch the
classical estimate into its success/failure blocks. Fidelity monotonicity under
this channel preserves the error bound. Let $`H`$ be $`\alpha^2`$ times the
transpose of its successful block. The direct-sum formula for generalized
Bures distance gives

```math
d_B(K,H)\le\alpha\zeta.
```

The acquired $`H`$ need not have the exact normalization and tangent constraints
of a physical local model. Use a strictly feasible convex fit:

```math
K_*\succeq0,\qquad \operatorname{Tr}(K_*)_{00}=1,
\qquad (K_*)_{00}\succeq\frac\nu2 I_r,
```

```math
\operatorname{Re}\operatorname{Tr}(K_*)_{0i}=0,
\qquad \sum_{i=1}^{p}\operatorname{Tr}(K_*)_{ii}\le2L^2.
```

The relaxed constants give room to certify finite arithmetic at boundary
points. They do not remove the original promises $`K_{00}\succeq\nu I_r`$,
$`\sum_i\operatorname{Tr}K_{ii}\le L^2`$ on the input.

For a target factor radius $`\gamma`$, choose

```math
\zeta\le\frac{\gamma}{8\alpha},\qquad
\lambda=\frac{\gamma^2}{64m},\qquad H_+=H+\lambda I_m.
```

Since $`d_B(H,H_+)\le\sqrt{m\lambda}`$, the regularized preliminary matrix
satisfies $`d_B(K,H_+)\le\gamma/4`$. Find a rational feasible pair
$`(K_*,Z)`$ obeying

```math
\begin{pmatrix}K_*&Z\\Z^\dagger&H_+\end{pmatrix}\succeq0,
\qquad
\operatorname{Tr}K_*+\operatorname{Tr}H_+
 -2\operatorname{Re}\operatorname{Tr}Z\le\frac{\gamma^2}{4}.
```

This **primal feasibility certificate** proves
$`d_B(K_*,H_+)\le\gamma/2`$, hence
$`d_B(K_*,K)\le3\gamma/4\le\gamma`$. No unverified solver accuracy is used to
infer proximity. The true $`K`$ has squared distance at most $`\gamma^2/16`$,
so the certificate target has positive slack. The regularized lower block and
the relaxed physical constraints supply the finite-precision construction in
[certified optimization](certified_optimization.md).

## 4. Stability constants and confidence interval

The physical fit bounds its derivative norm by $`\sqrt2L`$. Let

```math
b=2+(1+\sqrt2)L,
\qquad
\gamma_0=\min\left\{\frac\nu{4p},
 \frac{\sqrt{\nu j}}{16\sqrt p\,b}\right\},
```

```math
A=\frac{768p^{3/2}b}{\sqrt\nu\,j^{3/2}}
 +\frac{1024p}{\nu j}.
```

Choose $`\gamma\le\min\{\gamma_0,\varepsilon/(2A)\}`$. Factor alignment gives
$`\|\rho-\rho_*\|_1\le2\gamma`$ and combined derivative error at most
$`b\gamma`$. The [Fisher-transfer and stability theorems](holevo_compression_and_stability.md)
then give

```math
J_*\succeq\frac j8 I_p,
\qquad |C_H(K_*)-C_H(K)|\le A\gamma\le\frac\varepsilon2.
```

No numerical eigenspectrum of the unknown true Fisher matrix is being used as
a certificate: $`J\succeq jI_p`$ is an input promise, and its transfer is proved.
The fitted state weights are at least $`\nu/2`$ by its certified constraints.

Compute certified lower and upper values $`\ell\le C_H(K_*)\le u`$ with
$`u-\ell\le\varepsilon`$. The interval

```math
[\ell-A\gamma,\ u+A\gamma]
```

contains the true value whenever acquisition succeeds. Its midpoint has
absolute error at most $`\varepsilon`$. The lower endpoint may also be clipped
at zero. The [classical certificate formulation](certified_optimization.md)
works directly with rational Gram data and avoids an unverified numerical
square root. Optimizer runtime and certificate checking are charged separately
from quantum acquisition. If exact affine elimination finds an inconsistent fitted
model, or its minimum real-covariance cost exceeds the proved $`8p/j`$ bound,
return a failure marker. These aborts cannot occur on the acquisition success
event and therefore remain within its failure allowance.

## 5. Complete resource statement

Let $`n=m+1`$, let $`G_V`$ be the gate cost of the full label-state preparation,
and let $`a_V`$ be its qubit width. Including the conditioning threshold, the
number of calls to that preparation and its inverse is

```math
Q=\widetilde O\!\left[
 n^2\alpha\max\left\{\frac{2A}{\varepsilon},\frac1{\gamma_0}\right\}
 \right].
```

A conservative total quantum gate bound is

```math
\widetilde O\!\left[
 n^2\alpha\max\left\{\frac{2A}{\varepsilon},\frac1{\gamma_0}\right\}
 (G_V+a_V+n^{5/2})\right].
```

The tomography primitive's classical postprocessing, the small eigendecompositions,
and the finite convex programs remain polynomial classical work with their
accuracy and conditioning parameters included. The certificate construction
uses neither a lower eigenvalue bound on the tangent Gram nor tomography of a
$`D\times D`$ state. It also does not claim a small total cost when $`r,p`$,
$`\alpha`$, or the preparation circuit are large.

For fixed model size and conditioning, the precision dependence is
$`\widetilde O(1/\varepsilon)`$. A fixed-dimensional opaque-oracle subfamily
requires $`\Omega(1/\varepsilon)`$ coherent calls and
$`\Omega(1/\varepsilon^2)`$ purified copies; see
[applications](applications_and_limits.md). Those lower bounds do not apply
unchanged to a white-box toy circuit whose classical parameters reveal its
answer. Removing the state-weight promise permits a stronger lower bound in
the [rare-branch example](conditioning_limits.md).

## 6. What this gives a circuit designer

The output can compare candidate probe circuits at a specified operating point,
including the cost of incompatibility that the scalar SLD bound misses.
Classically contractible circuits, direct overlap estimation, and special
closed-form models remain necessary baselines. The [application account](applications_and_limits.md)
makes those comparisons and separates the proven query advantage from an
unestablished hardware crossover.
