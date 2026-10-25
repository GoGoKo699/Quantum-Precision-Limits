# Related work and the precise contribution

The project combines quantum tomography with a classical evaluation of the Holevo Cramér–Rao bound. Its two relevant outputs are a classical Bures approximation of a reduced state supplied through an arbitrary purification circuit, and a certified interval for the actual Holevo value of a conditioned, low-rank local model. Several ingredients already have close precedents. This account distinguishes those ingredients from the quantitative guarantees proved here. The literature comparison includes primary sources available through 25 October 2026.

## 1. Established metrology reductions and optimization

The Holevo bound is an established operational target. Its role in multiparameter estimation, relation to the SLD bound, local asymptotic interpretation, and attainability qualifications are explained in the review by [Demkowicz-Dobrzański, Górecki and Guţă](https://arxiv.org/html/2001.11742v2). Computing a QFIM alone generally does not settle simultaneous attainable precision. That distinction is background for the present task, not a new metrological principle.

[Albarelli, Friel and Datta (2019)](https://arxiv.org/html/1906.05724v1) express finite-dimensional Holevo evaluation as a semidefinite program and reduce its size for rank-deficient states. In particular, kernel–kernel estimator blocks can be discarded, leaving a real quotient of dimension $`2Dr-r^2`$ for a rank-r state in dimension $`D`$. The present classical optimization uses this established convex structure. It does not propose a new general-purpose SDP algorithm.

[Sidhu, Ouyang, Campbell and Kok (2021)](https://arxiv.org/html/1912.09218v4) give an especially strong comparator for two parameters: a scalar dual optimization and associated linear equations. Their treatment also identifies finite active blocks in noisy bosonic examples. The project's exact two-parameter verifier uses this quadratic-dual structure. Its contribution is a reproducible rational primal/dual certificate for the chosen circuit, including an explicit discrepancy from the SLD benchmark; the optimization principle is inherited.

### The closest dimensional-reduction precedent

[Zhou (2026), Theorem 1 and Appendix B](https://arxiv.org/html/2605.03975v2), relates mixed-state estimation to estimation of a purification with environment nuisance parameters. The pure-state problem is evaluated through tangent-vector inner products. This already yields a route from polynomially many accessible circuit overlaps to a small metrology problem. It is therefore a direct precedent for avoiding full ambient-state tomography. Zhou's Appendix D also discusses fixed-rank continuity, with the positive spectrum locally separated from zero.

The present projection onto the span of factor columns and derivatives makes the compression particularly direct: its dimension is at most $`r(p+1)`$, and positive covariance-Gram monotonicity proves preservation of the Holevo value. This compression is useful for the algorithm and its proof, but exact compression alone is not the novelty claim. The additional work is quantitative: learning the required Gram in a metric appropriate to factor perturbations, repairing local unbiasedness, charging conditioning and access costs, and certifying the resulting numerical interval.

[Yamaguchi and Tajima (September 2026)](https://arxiv.org/html/2609.31601v1) give another close characterization. The Holevo feasible region equals the intersection of matrix Cramér–Rao regions associated with an entire interpolating QFI family, including for rank-deficient locally identifiable models. This is stronger than access to the SLD QFIM. It supplies a possible alternative algorithmic route, but its theorem is not itself a finite-query circuit algorithm. In general, its matrix-intersection characterization must not be replaced by maximizing one scalar weighted inverse-QFI expression. An implementation based on it would have to charge the family evaluation, endpoint conditioning, and outer optimization accuracy.

## 2. The tomography improvement and its inherited ingredients

Fix the Bures convention


```math
B(\rho,\sigma)=\sqrt{2-2\|\sqrt\rho\sqrt\sigma\|_1}.
```


If fidelity is instead defined as the squared quantity $`F=\|\sqrt\rho\sqrt\sigma\|_1^2`$, small infidelity $`1-F`$ is of order $`B^2`$. This conversion is essential when comparing inverse-error formulas.

[van Apeldoorn, Cornelissen, Gilyén and Nannicini (SODA 2023)](https://real.mtak.hu/191396/1/1.9781611977554.ch47.pdf) study precisely the coherent input model needed here: a preparation unitary and its inverse, with arbitrary purification workspace. Their published Theorem 8.2 provides Schatten-norm tomography; the relevant Frobenius cost is approximately $`d\sqrt r/\eta`$. Their trace-norm algorithm uses approximately $`dr/\eta`$ queries, and Theorem 9.3 provides the matching trace-norm query lower bound. The project's Bures learner invokes the Frobenius routine as a subroutine. Applying the generic trace-to-Bures inequality to the published trace algorithm would instead require approximately $`dr/\epsilon_B^2`$ queries for Bures error $`\epsilon_B`$.

[Flammia and O'Donnell (2024)](https://arxiv.org/html/2305.18519v2) supply the closest geometric precedent. Their copy-model tomography construction adaptively resolves spectral mass and controls stronger errors through weighted chi-squared quantities. The present algorithm uses related spectral peeling and weighted error accounting. Those ideas are credited to that work; a copy theorem cannot simply be reused with its resource model silently replaced. The coherent construction must explicitly prepare residual states, estimate their probability mass, account for amplification, and propagate conditional errors.

Amplitude estimation and fixed-point amplification are established primitives, respectively from [Brassard, Høyer, Mosca and Tapp](https://arxiv.org/abs/quant-ph/0005055) and [Yoder, Low and Chuang](https://arxiv.org/abs/1409.3305). Standard explicit basis synthesis and finite-precision gate compilation are charged separately. No array of learned matrix entries is assumed to reside in QRAM.

### What the Bures theorem adds

For a rank-r state on a d-dimensional system, supplied as the reduction of an arbitrary-size coherent purification, the theorem gives a classical approximation with Bures error $`\epsilon_B`$ using


```math
\widetilde O(dr/\epsilon_B)
```


preparation/inverse queries. It imposes no eigenvalue lower bound. The trace-norm lower bound of van Apeldoorn et al. implies the same-order Bures query lower bound, because a Bures approximation is also a trace-distance approximation. This is a query-optimal statement up to logarithms, with gate, workspace, classical processing, and output costs accounted for separately.

If a preparation unitary for a purification in dimension $`d\times r`$ is already supplied, pure-state query tomography immediately provides this scaling. The additional capability is handling the supplied purification when its environment is much larger and an efficient minimal-purification unitary is unavailable. In the Holevo application, the reduced state is a small label Gram, while its purification contains the large physical probe. That is why this distinction matters even for an original probe with a short low-rank purification.

The mechanism combines established tools. Its mathematical content is the uniform error and resource accounting for that combination: the residual preparation cost cancels the scale dependence of the required Frobenius tolerance, without an unknown spectral floor or an uncharged environment-compression step. The related-work distinction should be stated at this level, rather than as a claim to have invented adaptive tomography or low-rank quantum estimation.

## 3. Strong copy-access and fidelity comparators

[Pelecanos, Spilecki and Wright (STOC 2026)](https://arxiv.org/abs/2510.07788) achieve optimal $`O(dr/\epsilon_B^2)`$ copy complexity for Bures tomography. This is a stronger copy baseline than applying elementary entrywise tomography. The coherent improvement uses a stronger interface—reusable preparation and inverse circuits—and must be described as such.

[Pelecanos, Spilecki, Tang and Wright, *Mixed state tomography reduces to pure state tomography*](https://arxiv.org/html/2511.15806v1), strengthen the copy comparison further. Their method uses $`O((dr+\log(1/\delta))/(1-F))`$ copies with polynomial gate complexity. Theorem 1.1 maps n copies of a mixed state into a mixture over n consistent copies of a random minimal purification. This does not supply a reusable preparation unitary and inverse for one fixed chosen purification. Consequently pure-state **query** tomography cannot be substituted into that copy reduction without additional work. The paper is nevertheless a major conceptual predecessor for separating reduced-state dimension from the size of a supplied environment.

The same authors' [July 2026 spectrum-estimation paper](https://arxiv.org/abs/2607.27117) develops refined tomography controls in the copy model. These reinforce the need to compare metrics and access assumptions explicitly; a bound stated in inverse infidelity is an inverse-square Bures-distance bound.

[Utsumi, Nakata, Wang and Takagi (2025)](https://arxiv.org/html/2509.03619v3) study Uhlmann transformations and fidelity estimation using purified access. Their results are close in oracle model and geometry, but are transformation and property-estimation results with their own rank or singular-value dependence. They do not directly supply a classical Bures tomography guarantee with arbitrary supplied environment and the stated optimal query scaling. A minimal-purification preparation oracle is not a free consequence of an Uhlmann transformation.

[Wang (June 2026)](https://arxiv.org/abs/2606.26034) studies fidelity to a known reference state. Its sample upper bound for a rank-r reference is $`O(r^2/\epsilon_F^2)`$ for additive fidelity error, with a different lower bound and task. Neither that scalar output nor its error convention should be identified with full Bures tomography.

The checked primary literature establishes these close constructions, but no source above states the same arbitrary-purifier Bures-query theorem. This is a bounded literature comparison, not a guarantee that no other equivalent theorem exists.

## 4. Current quantum-geometry algorithms and classical structure

[Wang, Zhao, Yuan and Cai (17 October 2026)](https://arxiv.org/html/2610.24686) construct SLD operators through quadrature, obtain near-Heisenberg diagonal-QFI estimation under a spectral promise, and provide coherent QFIM and QGT block encodings. This is the closest recent algorithmic comparison for mixed-state geometry. It explicitly recognizes that the QFIM need not determine joint attainability. Its reported outputs differ from a certified classical Holevo value. Conversely, the present small-model reconstruction is not an improved block encoding for the full QGT. Both approaches must charge derivative access, normalization, and conditioning.

That paper also cites a 2026 Peking University thesis by Yukun Zhang on spectral functionals and information-geometric quantities. An openly accessible full thesis was not located in the checked sources. Its contents beyond the cited description remain unverified; the comparison does not assert comprehensive priority over that work.

For Gaussian models, [Chang, Genoni and Albarelli (2026)](https://arxiv.org/abs/2504.17873) already give efficient moment-based evaluation of Holevo, RLD, and SLD bounds. Infinite or large Hilbert dimension therefore does not by itself establish a useful quantum advantage. More generally, the few overlaps required by the compressed model can be computed by classical stabilizer methods, matrix-product-state simulation, or tensor-network contraction whenever the preparation circuit has the relevant structure. Those are the matched classical baselines, not full density-matrix tomography by default.

The application note gives two deliberately different checks. A four-qubit parity-preserving circuit has certified rank, conditioning, an actual incompatibility penalty, and an exact rational Holevo enclosure; it is classically easy. A separate universal-circuit family shows that fixed QFI and low rank do not force the problem to be classically easy. That family's analytic solution is also its strongest direct quantum comparator. Neither example establishes a practical speed advantage for a proposed device.

## 5. Contribution boundary

The scientific endpoint consists of the coherent Bures theorem, its conversion into a finite-precision Holevo certificate under explicit occupied-spectrum and identifiability promises, and exact application/lower-bound examples that expose the costs those promises control. The low-rank preparation-tail result concerns a fixed mixture before a common bounded unitary encoding; trace closeness by itself is insufficient. The separate small-eigenvalue obstruction prevents promoting the conditioned inverse-error claim to an unconditional bounded-generator statement.

The work does not claim a new Holevo SDP, a new relation between QFI and incompatibility, an optimal-measurement construction, unrestricted noise robustness, or a general separation from all classical circuit simulations. Its useful capability is the certified evaluation of a specified mixed-probe precision limit from explicitly paid coherent circuits, including singular derivative Grams. This boundary is necessary for the algorithmic comparison to be meaningful.
