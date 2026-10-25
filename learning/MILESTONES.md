# Project understanding map: 30 milestones

This map organizes **Quantum Precision Limits** into exactly **30 milestones**, M01–M30. It describes the project's understanding targets for any reader, including an LLM. The targets and their order are determined by the science.

The allocation is **12 background/preliminaries milestones** and **18 narrative/innovation milestones**. The background joins three foundations: local quantum estimation, the geometry of low-rank factors, and the computational meaning of coherent access and certification. The narrative needs more space because the result connects several separate arguments: exact compression, Bures acquisition, physical fitting and stability, finite optimization certificates, and the examples that delimit the claim. Each phase ends by assembling its components into a complete account.

The single educational anchor is Rafał Demkowicz-Dobrzański, Wojciech Górecki and Mădălin Guţă, *Multi-parameter estimation beyond Quantum Fisher Information*, **arXiv:2001.11742v2** (1 September 2020). All anchor section and equation numbers below refer to that version. The [estimation bridge](ESTIMATION_BRIDGE.md) supplies a targeted route, translates notation, works through a complete qubit example, and identifies the steps supplied by this project. The review explains the estimation problem; it does not supply the coherent tomography algorithm or its finite computational certificates. Those mechanisms and their primary research dependencies are identified in the canonical science files.

Use this map alongside the [project overview](../README.md) and [model](../science/model.md). Source pointers distinguish an anchor connection from the project statement that actually establishes a result. The map gives a global account; full technical proofs remain in the linked scientific sources.

## Route through the project

| Phase | Milestones | Connection established |
|---|---|---|
| Background/preliminaries | [M01](#m01-parameterized-probes-and-measurement-records)–[M06](#m06-attainability-and-the-meaning-of-a-local-bound) | Physical measurement, local calibration, Fisher information, incompatibility, and the meaning of the Holevo value |
| Background/preliminaries | [M07](#m07-purifications-tangents-and-occupied-support)–[M11](#m11-rank-and-conditioning-as-separate-promises) | Factor geometry, access, certificates, and distinct conditioning assumptions |
| Background/preliminaries | [M12](#m12-the-complete-computational-question) | The complete input, output, and resource contracts |
| Narrative/innovation | [M13](#m13-exact-tangent-compression)–[M15](#m15-the-bures-acquisition-objective) | The large physical problem becomes a small state-learning task with the right metric |
| Narrative/innovation | [M16](#m16-residual-estimation-with-paid-conditioning)–[M19](#m19-finite-precision-and-full-work-costs) | Paid residual preparation, spectral peeling, normalization, and implementation establish the tomography theorem |
| Narrative/innovation | [M20](#m20-physical-fitting-with-certified-slack)–[M24](#m24-the-end-to-end-confidence-interval) | Acquired data become a physical local model and a certified Holevo interval |
| Narrative/innovation | [M25](#m25-a-four-qubit-sensing-benchmark)–[M28](#m28-approximate-low-rank-with-a-physical-certificate) | A concrete circuit, matched baselines, conditioning obstruction, and preparation-tail guarantee establish usefulness and limits |
| Narrative/innovation | [M29](#m29-attribution-and-the-strength-of-the-evidence)–[M30](#m30-the-complete-scientific-story) | Attribution and evidence support the complete contribution and its boundaries |

## Phase I: background/preliminaries

### M01 Parameterized probes and measurement records

**Central question.** What physical object carries the unknown parameters, and what does a measurement reveal?

**Understanding target.** Interpret a smooth family of density operators $`\rho_\theta`$, an operating point $`\theta_0`$, and a POVM whose outcome probabilities are $`p_\theta(m)=\mathrm{Tr}(\rho_\theta M_m)`$. Distinguish the quantum state, the random outcome, and the estimator constructed from that outcome. The preparation and encoding specify the statistical model whose precision is being evaluated.

**Dependencies and story role.** None within this map. This fixes the physical objects that all later geometry and algorithms describe.

**Conceptual connection.** A known circuit model can have unknown sensing parameters. Computing its precision benchmark is a different task from estimating those parameters in an experiment.

**Source pointers.** [Anchor §1.2, Eqs. (1)–(3)][anchor]; [local sensing model](../science/model.md#local-sensing-model); [the specified circuit-model task](../science/applications_and_limits.md#1-the-output-is-a-precision-benchmark-for-a-known-circuit-model).

### M02 Local unbiasedness and estimation risk

**Central question.** What does it mean for an estimator to track small parameter changes correctly?

**Understanding target.** Read weighted mean-square error as $`\mathrm{Tr}(W\Sigma)`$ and local unbiasedness as calibration of the mean and its first derivatives at one point. Translate these requirements to centered Hermitian estimator operators $`X_i`$ with $`\mathrm{Tr}(\rho X_i)=0`$ and $`\mathrm{Tr}(D_jX_i)=\delta_{ij}`$. Explain why the derivative condition excludes a constant answer that is accidentally correct at the operating point.

**Dependencies and story role.** [M01](#m01-parameterized-probes-and-measurement-records). This supplies the affine constraints that will survive compression, perturbation, and optimization.

**Conceptual connection.** Local unbiasedness is a first-order calibration condition. It does not assert global unbiasedness or finite-sample optimality.

**Source pointers.** [Anchor §1.2, Eqs. (3)–(5), and §2.3, Eqs. (15)–(16)][anchor]; [local calibration in the bridge](ESTIMATION_BRIDGE.md#local-unbiasedness-as-a-calibration-condition); [the quantity computed](../science/model.md#the-quantity-computed).

### M03 Classical information and repeated experiments

**Central question.** How does local sensitivity of an outcome distribution limit estimation error?

**Understanding target.** Interpret the classical Fisher matrix as a weighted quadratic measure of changes in outcome probabilities. Connect local unbiasedness to the matrix Cramér–Rao inequality and connect independent repetitions to additive information. Keep the coefficient of the asymptotic $`1/n`$ estimation cost distinct from a finite-sample performance guarantee.

**Dependencies and story role.** [M02](#m02-local-unbiasedness-and-estimation-risk). This establishes the statistical comparison that quantum measurement incompatibility will modify.

**Conceptual connection.** The number $`n`$ of probes used to estimate a physical parameter is different from the number of circuit queries used to compute a precision coefficient.

**Source pointers.** [Anchor §2.1, Eqs. (7)–(10)][anchor]; [meaning of the computed value](../science/model.md#the-quantity-computed); [computation versus experimental precision](../science/applications_and_limits.md#1-the-output-is-a-precision-benchmark-for-a-known-circuit-model).

### M04 SLD information and local identifiability

**Central question.** What information is available before a measurement is chosen?

**Understanding target.** Use the symmetric logarithmic derivative relation $`D_i=(\rho L_i^{\mathrm{SLD}}+L_i^{\mathrm{SLD}}\rho)/2`$ and its real covariance matrix $`J`$. Interpret $`J\succeq jI_p`$ as a quantitative promise that all parameter directions are locally distinguishable. Understand the SLD-dual estimator tuple and why $`\mathrm{Tr}(WJ^{-1})`$ is a benchmark.

**Dependencies and story role.** [M02](#m02-local-unbiasedness-and-estimation-risk), [M03](#m03-classical-information-and-repeated-experiments). This introduces both the familiar comparator and the identifiability assumption used in the algorithm.

**Conceptual connection.** Identifiability concerns parameter directions. It is not a lower bound on occupied state probabilities, and it does not guarantee simultaneous attainability of the SLD benchmark.

**Source pointers.** [Anchor §2.2, Eqs. (11)–(13), and §2.6, Eqs. (34)–(38)][anchor]; [SLDs and the Fisher promise](../science/model.md#local-sensing-model); [the SLD-dual construction](../science/holevo_compression_and_stability.md#3-model-space-stability).

### M05 Incompatibility and the Holevo objective

**Central question.** Why can two informative parameters still compete for measurement precision?

**Understanding target.** Form the complex covariance $`Z_{ij}=\mathrm{Tr}(\rho X_iX_j)`$ and interpret the two terms in $`h_W(Z)=\mathrm{Tr}(W\mathrm{Re}Z)+\|\sqrt W\,\mathrm{Im}Z\,\sqrt W\|_1`$. Understand why retaining only the real covariance loses incompatibility information. Relate positive covariance additions to the monotonicity of this objective.

**Dependencies and story role.** [M02](#m02-local-unbiasedness-and-estimation-risk), [M04](#m04-sld-information-and-local-identifiability). This identifies the actual scalar target and the order property behind later reductions.

**Conceptual connection.** The bridge's pure-qubit model has $`J=I_2`$, SLD cost $`2`$, and Holevo cost $`4`$. A well-conditioned QFIM therefore does not determine the simultaneous precision cost.

**Source pointers.** [Anchor §§2.3, 2.5–2.6, 2.8, and §3.1.1][anchor]; [complete qubit calculation](ESTIMATION_BRIDGE.md#a-complete-two-parameter-qubit-calculation); [objective facts](../science/holevo_compression_and_stability.md#1-objective-facts).

### M06 Attainability and the meaning of a local bound

**Central question.** When does the local Holevo value describe a precision that measurements can attain?

**Understanding target.** Separate the optimization defining $`C_H`$, its lower-bound meaning, its multiple-copy scaling, and an attainment theorem. Understand the role of localization and collective measurements in regular asymptotic models. Recognize the pure-state attainment result and the full-rank conditions in the anchor's general QLAN treatment; the latter does not automatically establish attainment for every rank-deficient boundary model used by this project.

**Dependencies and story role.** [M03](#m03-classical-information-and-repeated-experiments), [M05](#m05-incompatibility-and-the-holevo-objective). This closes the estimation-theory part of the background before introducing the computation.

**Conceptual connection.** An accurate value certificate can be valid even when the algorithm supplies neither an attaining measurement nor a finite-copy estimator. Its operational interpretation needs an applicable estimation-theory result.

**Source pointers.** [Anchor §§2.10–2.11, 4.5, and 4.6.1–4.6.4][anchor]; [operational scope](../science/model.md#the-quantity-computed); [application qualifications](../science/applications_and_limits.md#1-the-output-is-a-precision-benchmark-for-a-known-circuit-model).

### M07 Purifications, tangents, and occupied support

**Central question.** How can a low-rank local model be represented without listing a large density matrix?

**Understanding target.** Interpret $`\rho=\Psi\Psi^\dagger`$ and $`D_i=\Phi_i\Psi^\dagger+\Psi\Phi_i^\dagger`$, with $`r`$ factor columns and $`p`$ derivative directions. Connect normalization to $`\|\Psi\|_F=1`$ and traceless tangents to $`\mathrm{Re}\,\mathrm{Tr}(\Psi^\dagger\Phi_i)=0`$. For unitary encoding, identify $`\Phi_i=-iG_i\Psi`$.

**Dependencies and story role.** [M01](#m01-parameterized-probes-and-measurement-records), [M02](#m02-local-unbiasedness-and-estimation-risk). This introduces the local data from which the small computational representation is built.

**Conceptual connection.** State support and tangent support need not coincide: derivatives may connect occupied directions to the state kernel. Discarding every kernel-related estimator block would remove relevant information.

**Source pointers.** [Anchor §2.6, Eqs. (35)–(36), for support and kernel blocks][anchor]; [factor model](../science/model.md#local-sensing-model); [bridge to tangent factors](ESTIMATION_BRIDGE.md#from-a-large-probe-to-a-small-tangent-gram).

### M08 Gram data, factor alignment, and Bures distance

**Central question.** Which notion of approximation controls the state and its derivatives together?

**Understanding target.** Form $`V=[\Psi,\Phi_1,\ldots,\Phi_p]`$ and $`K=V^\dagger V`$. Interpret a Gram matrix as all pairwise inner products of its columns. Understand generalized Bures distance as the smallest Frobenius distance between compatible aligned factors, including singular Grams. Keep root fidelity, squared fidelity, Bures distance, and squared Bures distance distinct.

**Dependencies and story role.** [M07](#m07-purifications-tangents-and-occupied-support). This supplies the metric that connects data acquisition to perturbations of a local model.

**Conceptual connection.** The factor comparison is optimized over a common isometric alignment; it is not the distance between arbitrary unaligned purifications. Small entrywise errors can lose a square root when converted to factor accuracy near a rank boundary.

**Source pointers.** [Anchor §2.5 for the covariance geometry motivating the connection][anchor]; [distance conventions](../science/model.md#distances-and-their-operational-role); [Bures acquisition bridge](ESTIMATION_BRIDGE.md#why-bures-distance-is-the-acquired-accuracy); [factor perturbation bounds](../science/holevo_compression_and_stability.md#4-factor-distance-and-fisher-condition-transfer). The factor metric and acquisition mechanism are supplied by the project and its cited research dependencies.

### M09 Coherent access and charged resources

**Central question.** What additional interface makes a coherent quantum algorithm different from a collection of samples?

**Understanding target.** Distinguish independently prepared copies from reusable preparation and inverse circuits. Track queries, gates, workspace, classical arithmetic, and output size separately. Understand why controlled generator operations, derivative normalizations, and conditional preparation must be paid for. A short circuit is a substantive input representation, not permission to load an exponentially long classical vector for free.

**Dependencies and story role.** [M01](#m01-parameterized-probes-and-measurement-records), [M03](#m03-classical-information-and-repeated-experiments), [M07](#m07-purifications-tangents-and-occupied-support). This fixes the resource contract before any speed comparison.

**Conceptual connection.** The general tomography theorem needs $`U,U^\dagger`$ but no controlled unknown $`U`$. Constructing the sensing input can still require controlled known generator circuits, whose costs remain explicit.

**Source pointers.** [Anchor §1.2 and §2.10 for the statistical input and copy setting][anchor]; [tomography access model](../science/bures_tomography.md#1-model-conventions-and-theorem); [coherent circuit input](../science/circuit_access.md#1-the-coherent-input); [paid-conditioning bridge](ESTIMATION_BRIDGE.md#coherent-access-and-paid-conditioning). The coherent oracle model is an additional project assumption.

### M10 Positivity, semidefinite programs, and certificates

**Central question.** How can an optimization result carry evidence that its value is enclosed correctly?

**Understanding target.** Interpret positive semidefinite order, Gram matrices, and a Schur-complement block constraint. Connect the Holevo epigraph to a semidefinite program. Distinguish a primal feasible upper bound, a dual feasible lower bound, and a numerical proposal that has not yet been checked. Understand why strict feasibility and bounded search regions matter for a finite rational algorithm.

**Dependencies and story role.** [M05](#m05-incompatibility-and-the-holevo-objective), [M08](#m08-gram-data-factor-alignment-and-bures-distance). This supplies the language for physical fitting and the final classical certificate.

**Conceptual connection.** A small primal–dual gap is useful only when feasibility is valid. Near a singular Gram, an unchecked floating-point sign is not a proof of positivity.

**Source pointers.** [Anchor §2.4, Eqs. (23)–(27)][anchor]; [objective epigraph](../science/holevo_compression_and_stability.md#1-objective-facts); [finite optimization scope](../science/certified_optimization.md#what-is-proved-and-what-is-implemented); [explicit dual certificate](../science/certified_optimization.md#5-an-explicit-dual-and-an-independently-checkable-interval).

### M11 Rank and conditioning as separate promises

**Central question.** Which small dimensions and positive lower bounds are actually assumed?

**Understanding target.** Distinguish physical dimension $`D`$, tomography dimension $`d`$, occupied state rank, factor-column count, and label-state rank. Separate the occupied-weight floor $`K_{00}\succeq\nu I_r`$, the Fisher floor $`J\succeq jI_p`$, and the unrestricted spectrum of the full tangent Gram. The sensing theorem requires $`\Psi`$ to have $`r`$ independent columns; arbitrary redundant purification columns do not satisfy that promise.

**Dependencies and story role.** [M04](#m04-sld-information-and-local-identifiability), [M07](#m07-purifications-tangents-and-occupied-support), [M08](#m08-gram-data-factor-alignment-and-bures-distance), [M09](#m09-coherent-access-and-charged-resources). This assembles the assumptions that later resource claims must retain.

**Conceptual connection.** The general Bures learner permits an arbitrarily large purifying register and imposes no occupied-eigenvalue floor. Its Holevo application has additional state and Fisher conditioning promises; these are compatible but different theorems.

**Source pointers.** [Anchor §2.6 for state-kernel distinctions][anchor]; [dimension conventions](../science/model.md#two-different-dimensions); [local promises](../science/model.md#local-sensing-model); [factor-input qualification](../science/circuit_access.md#1-the-coherent-input).

### M12 The complete computational question

**Central question.** What precisely must the two algorithms output, and under which contracts?

**Understanding target.** Assemble the background into two tasks: learn a reduced density matrix to Bures accuracy from a purification circuit, and evaluate a specified local Holevo value to absolute accuracy from coherent factor-and-derivative access. Identify the classical output, failure probability, metric or value tolerance, and charged resources for each. The second task retains occupied-state and Fisher promises and allows a singular tangent Gram.

**Dependencies and story role.** [M06](#m06-attainability-and-the-meaning-of-a-local-bound), [M08](#m08-gram-data-factor-alignment-and-bures-distance), [M10](#m10-positivity-semidefinite-programs-and-certificates), [M11](#m11-rank-and-conditioning-as-separate-promises). This integrates Phase I into the full scientific question.

**Conceptual connection.** The computational confidence interval concerns the value of a known local model. It is not an experimental confidence region for the unknown parameters or a construction of their optimal measurement.

**Source pointers.** [Anchor §1.2 and §§2.3–2.4 for the underlying estimation target][anchor]; [general tomography theorem](../science/model.md#general-tomography-theorem); [algorithmic scope](../science/model.md#scope-of-the-algorithmic-claim); [complete resource statement](../science/circuit_access.md#5-complete-resource-statement).

## Phase II: narrative/innovation

### M13 Exact tangent compression

**Central question.** Why can a large physical Hilbert space be replaced by the span of a few local vectors?

**Understanding target.** Let $`S`$ be the column span of $`V`$. Follow the projection $`X_i\mapsto P_SX_iP_S`$: it preserves local unbiasedness and removes a positive covariance Gram. Objective monotonicity then preserves the Holevo value and the feasible real covariance region. The effective dimension is at most $`m=r(p+1)`$, and any factor of $`K`$ represents an equivalent compressed model.

**Dependencies and story role.** [M05](#m05-incompatibility-and-the-holevo-objective), [M07](#m07-purifications-tangents-and-occupied-support), [M08](#m08-gram-data-factor-alignment-and-bures-distance), [M12](#m12-the-complete-computational-question). This identifies the information the algorithm needs.

**Conceptual connection.** Compression preserves a local optimization even when $`S`$ includes state-kernel directions. It is not simply restriction to the occupied state support, and dimension reduction alone is not the project's novelty claim.

**Source pointers.** [Anchor §§2.5–2.6 for covariance and kernel reductions][anchor]; [exact compression theorem](../science/holevo_compression_and_stability.md#2-exact-compression); [closest dimensional-reduction precedent](../references/RELATED_WORK.md#the-closest-dimensional-reduction-precedent).

### M14 The label state and its failure symbol

**Central question.** How do paid circuits place all required overlaps into a small state?

**Understanding target.** Follow coherent preparation of unnormalized factor and derivative columns with known normalization $`\alpha`$. Tracing out the physical register yields the successful block $`K^{\mathsf T}/\alpha^2`$. Retaining one failure label gives a normalized state of dimension $`m+1`$ without free postselection. Its rank can also be $`m+1`$, irrespective of the physical probe's rank $`r`$.

**Dependencies and story role.** [M09](#m09-coherent-access-and-charged-resources), [M11](#m11-rank-and-conditioning-as-separate-promises), [M13](#m13-exact-tangent-compression). This turns the compressed data requirement into an explicit acquisition interface.

**Conceptual connection.** A low-rank probe may have a small preparation environment while the label state is purified by the entire large physical probe. This is the application of arbitrary-purifier tomography.

**Source pointers.** [Anchor §2.4 for the data needed by Holevo optimization][anchor]; [coherent input construction](../science/circuit_access.md#1-the-coherent-input); [failure symbol and transpose](../science/circuit_access.md#2-a-failure-symbol-prevents-free-postselection); [explicit sensing preparation](../science/applications_and_limits.md#purification-and-derivative-circuit). The label construction is project-supplied.

### M15 The Bures acquisition objective

**Central question.** Why learn the label state in Bures distance rather than only estimate its entries?

**Understanding target.** Connect Bures accuracy of the label state to aligned factor accuracy after rescaling by $`\alpha`$. Understand the square-root loss of a generic entrywise or trace-to-Bures conversion near singular Grams. State the separate tomography goal: a $`d`$-dimensional state of rank at most $`r`$, with unrestricted supplied purification, learned using $`\widetilde O(dr/\varepsilon_B)`$ queries without an eigenvalue floor.

**Dependencies and story role.** [M08](#m08-gram-data-factor-alignment-and-bures-distance), [M11](#m11-rank-and-conditioning-as-separate-promises), [M14](#m14-the-label-state-and-its-failure-symbol). This fixes the accuracy contract for the acquisition proof. In M16–M19, $`d,r`$ are the tomography target dimension and rank bound.

**Conceptual connection.** If a minimal $`r`$-dimensional purification is already accessible, pure-state tomography supplies the same scaling. The additional capability concerns a large supplied environment whose compression is not free.

**Source pointers.** [Anchor §2.5 for the target covariance geometry][anchor]; [Bures acquisition bridge](ESTIMATION_BRIDGE.md#why-bures-distance-is-the-acquired-accuracy); [tomography theorem](../science/bures_tomography.md#1-model-conventions-and-theorem); [what the theorem adds](../references/RELATED_WORK.md#what-the-bures-theorem-adds).

### M16 Residual estimation with paid conditioning

**Central question.** How can small remaining spectral mass be examined without an uncharged conditional oracle?

**Understanding target.** For $`R=P\rho P`$ of mass $`\tau`$, combine amplitude estimation, fixed-point amplification, and coherent Frobenius tomography. Explain why required subnormalized error of order $`\sqrt{\tau h/r}`$ makes the conditional tomography cost proportional to $`\sqrt{\tau}`$, cancelling the amplification cost proportional to $`1/\sqrt{\tau}`$. The resulting original-query cost is $`\widetilde O(dr/\sqrt h)`$. The ideal phase sequence acts as $`uP+v(I-P)`$ on the purification and preserves the rank bound.

**Dependencies and story role.** [M09](#m09-coherent-access-and-charged-resources), [M15](#m15-the-bures-acquisition-objective). This establishes the local resource-and-error mechanism of the tomography algorithm.

**Conceptual connection.** Approximate conditional preparation does not justify applying a rank promise to an arbitrary nearby noisy state. Exact rank preservation belongs to the specified ideal phase circuit; compilation error is handled separately.

**Source pointers.** [Anchor §2.4 for the downstream numerical-evaluation objective][anchor]; [established primitives](../science/bures_tomography.md#2-established-primitives-and-their-precise-role); [residual-estimation lemma](../science/bures_tomography.md#3-a-subnormalized-residual-estimation-lemma). These algorithmic ingredients are research dependencies beyond the anchor.

### M17 Spectral peeling and weighted error

**Central question.** How do local residual estimates combine into a guarantee that is uniform over tiny eigenvalues?

**Understanding target.** Follow the adaptive selection of estimated eigenvectors above a mass-dependent threshold. Every continuing round removes at least half the residual mass and retains at most $`r`$ directions. Assign each retained matrix entry to the first round that retains one of its coordinates. The recorded eigenvalue supplies its weighted-error denominator, while later rotations preserve the relevant Frobenius cross-block error.

**Dependencies and story role.** [M08](#m08-gram-data-factor-alignment-and-bures-distance), [M16](#m16-residual-estimation-with-paid-conditioning). This turns one residual guarantee into a controlled sum over logarithmically many scales.

**Conceptual connection.** The proof needs small matrix residuals rather than individually accurate eigenvectors or an unknown eigengap. Its weighted quadratic error is the bridge from Frobenius estimates to Bures accuracy.

**Source pointers.** [Anchor §2.4 for why finite acquired data must support an optimization][anchor]; [adaptive algorithm](../science/bures_tomography.md#4-the-adaptive-algorithm); [mass progress and weighted accumulation](../science/bures_tomography.md#5-progress-and-weighted-error-accumulation); [adaptive-geometry attribution](../references/RELATED_WORK.md#2-the-tomography-improvement-and-its-inherited-ingredients).

### M18 Normalization and query optimality

**Central question.** How does the retained subnormalized estimate become a complete density-matrix approximation?

**Understanding target.** Combine the exact normalization identity for the weighted error with the Bures cost of discarded mass. Understand why both the retained estimate's mass error and the omitted coherences are controlled. Choosing the tail scale of order $`\varepsilon_B^2`$ up to logarithms gives the stated Bures guarantee and query count. Relate the trace-distance inequality to the inherited $`\Omega(dr/\varepsilon_B)`$ query lower bound for sufficiently small accuracy and $`d\ge2`$.

**Dependencies and story role.** [M15](#m15-the-bures-acquisition-objective), [M16](#m16-residual-estimation-with-paid-conditioning), [M17](#m17-spectral-peeling-and-weighted-error). This completes the ideal tomography argument.

**Conceptual connection.** Query optimality is a worst-case statement in the specified coherent oracle model. It does not establish optimal gate cost or classical hardness of a transparent circuit description.

**Source pointers.** [Anchor §2.4 for the motivating numerical task][anchor]; [normalization and Bures guarantee](../science/bures_tomography.md#6-normalization-and-the-final-bures-guarantee); [query costs](../science/bures_tomography.md#71-queries-with-inherited-logarithms-visible); [lower bound and baselines](../science/bures_tomography.md#9-lower-bound-and-inherited-baselines).

### M19 Finite precision and full work costs

**Central question.** Which implementation details are necessary for the query theorem to describe a finite algorithm?

**Understanding target.** Track finitely specified nested basis rotations, residual-based diagonalization, rational factor output, and a global compilation-error budget. Explain why factor rounding preserves positivity and normalized output without another square-root accuracy loss. Include preparation gates, large-workspace phases, small-register synthesis, classical matrix processing, and output size. Fixed branch limits keep work finite even when a statistical estimate fails.

**Dependencies and story role.** [M09](#m09-coherent-access-and-charged-resources), [M10](#m10-positivity-semidefinite-programs-and-certificates), [M18](#m18-normalization-and-query-optimality). This completes the implementation contract of the general learner.

**Conceptual connection.** A compiled intermediate state may have higher rank than the ideal one. The proof compares complete output distributions; it does not reapply a rank-dependent theorem to each noisy intermediate state.

**Source pointers.** [Anchor §2.4 for numerical evaluation as the destination][anchor]; [gate, width, and classical costs](../science/bures_tomography.md#7-query-gate-width-and-classical-costs); [finite precision](../science/bures_tomography.md#8-finite-precision-and-implementation-robustness); [scope of the Bures checks](../verification/README.md#bures-geometry).

### M20 Physical fitting with certified slack

**Central question.** How can acquired Gram data be turned into a valid local sensing model?

**Understanding target.** Pinch the learned label state, rescale its successful block, and fit it to normalization and tangent constraints. Explain why the fit uses an occupied-state floor $`\nu/2`$ and derivative budget $`2L^2`$: these relaxed constraints provide an explicit interior point while preserving the input promises. A positive regularization is added to the classical reference matrix, and an exactly feasible fidelity-block witness certifies proximity.

**Dependencies and story role.** [M10](#m10-positivity-semidefinite-programs-and-certificates), [M11](#m11-rank-and-conditioning-as-separate-promises), [M14](#m14-the-label-state-and-its-failure-symbol), [M19](#m19-finite-precision-and-full-work-costs). This begins the conversion from acquisition to a valid value calculation.

**Conceptual connection.** Regularizing the fit's reference does not add noise to the unknown physical state or impose a floor on its full tangent Gram. A feasible distance certificate suffices; exact fit optimality is unnecessary.

**Source pointers.** [Anchor §§1.2 and 2.3 for normalization and local constraints][anchor]; [physical fit construction](../science/circuit_access.md#3-from-a-confidence-ball-to-a-physical-model); [strict feasibility and certificate](../science/certified_optimization.md#2-physical-fitting-without-a-hidden-boundary-assumption).

### M21 Stability and Fisher-promise transfer

**Central question.** Why does a nearby valid model have a nearby Holevo value?

**Understanding target.** Use occupied-state conditioning to choose bounded estimator representatives. Center a perturbed tuple and invert its almost-identity derivative-response matrix to repair local unbiasedness. Connect factor accuracy to state and derivative errors, transfer the true-model Fisher promise to $`J_*\succeq(j/8)I_p`$, and obtain a linear value bound $`|C_H(K_*)-C_H(K)|\le A\gamma`$ within a specified radius $`\gamma_0`$.

**Dependencies and story role.** [M04](#m04-sld-information-and-local-identifiability), [M05](#m05-incompatibility-and-the-holevo-objective), [M08](#m08-gram-data-factor-alignment-and-bures-distance), [M20](#m20-physical-fitting-with-certified-slack). This justifies converting acquisition accuracy into value accuracy.

**Conceptual connection.** The true Fisher floor is an input promise transferred by proof. A floating-point Fisher eigenvalue of a fitted model is not evidence that the unknown input satisfies that promise.

**Source pointers.** [Anchor §2.6 for SLD-dual geometry and §2.8 for the factor-two comparison][anchor]; [model-space stability](../science/holevo_compression_and_stability.md#3-model-space-stability); [factor and Fisher transfer](../science/holevo_compression_and_stability.md#4-factor-distance-and-fisher-condition-transfer); [explicit finite-fit constants](../science/holevo_compression_and_stability.md#5-a-finite-precision-version-with-an-interior-physical-fit).

### M22 Rational action variables for the Holevo program

**Central question.** How can the compressed optimization remain exact when its Gram matrix is singular?

**Understanding target.** Factor the rational fitted Gram with exact pivoted LDL arithmetic, retaining only its nonzero pivots. Replace each estimator by its action $`x_i=X_i\Psi`$ and characterize Hermitian realizability through $`\Psi^\dagger x_i=x_i^\dagger\Psi`$. Use the rational weighted coordinates to express calibration and covariance without numerically taking a matrix square root. A Schur-complement constraint then gives the small Holevo SDP.

**Dependencies and story role.** [M10](#m10-positivity-semidefinite-programs-and-certificates), [M13](#m13-exact-tangent-compression), [M20](#m20-physical-fitting-with-certified-slack). This supplies a finite optimization representation that tolerates singular tangent data.

**Conceptual connection.** The exact rank of a finite fitted matrix is not an estimate of the unknown input's tangent rank. Small rational pivots are accounted for through bit complexity, without a new physical Gram-gap promise.

**Source pointers.** [Anchor §2.4 for the inherited SDP structure and §2.6 for kernel quotients][anchor]; [rational compression](../science/certified_optimization.md#3-rational-compression-without-numerical-square-roots); [Hermitian actions and the primal program](../science/certified_optimization.md#4-action-variables-and-the-exact-small-holevo-sdp).

### M23 Primal and dual optimization certificates

**Central question.** What finite evidence proves that the optimized value lies between two reported endpoints?

**Understanding target.** Identify the affine equations and PSD blocks checked in a rational primal witness and a rational dual witness. Follow the exact weak-duality identity that makes the difference between upper and lower values nonnegative. Explain how explicit strict points, bounded search regions, and finite precision ensure that arbitrarily small certified gaps can be obtained. Keep proposal generation separate from certificate checking.

**Dependencies and story role.** [M10](#m10-positivity-semidefinite-programs-and-certificates), [M21](#m21-stability-and-fisher-promise-transfer), [M22](#m22-rational-action-variables-for-the-holevo-program). This completes the classical numerical guarantee.

**Conceptual connection.** The theorem supplies a finite generic optimization procedure. The repository's exact witness generators verify selected models; they do not implement that generic optimizer.

**Source pointers.** [Anchor §2.4 for numerical Holevo evaluation][anchor]; [explicit dual and interval](../science/certified_optimization.md#5-an-explicit-dual-and-an-independently-checkable-interval); [combined certificate and work](../science/certified_optimization.md#6-combined-certificate-and-computational-costs); [implemented checks](../science/certified_optimization.md#7-implemented-checks-and-their-limits).

### M24 The end-to-end confidence interval

**Central question.** How do the quantum and classical errors combine into one useful output?

**Understanding target.** Choose the factor radius using both $`\gamma\le\gamma_0`$ and $`A\gamma\le\varepsilon/2`$. A fitted-model interval $`[\ell,u]`$ of width at most $`\varepsilon`$ becomes $`[\ell-A\gamma,u+A\gamma]`$ for the true model; its midpoint has absolute error at most $`\varepsilon`$ on the acquisition success event. Account for label dimension, normalization, conditioning, and inverse calls in

```math
Q=\widetilde O\!\left((m+1)^2\alpha
\max\left\{\frac{A}{\varepsilon},\frac1{\gamma_0}\right\}\right).
```

**Dependencies and story role.** [M14](#m14-the-label-state-and-its-failure-symbol), [M19](#m19-finite-precision-and-full-work-costs), [M20](#m20-physical-fitting-with-certified-slack), [M21](#m21-stability-and-fisher-promise-transfer), [M23](#m23-primal-and-dual-optimization-certificates). This integrates the algorithm from input circuit to certified value.

**Conceptual connection.** The inverse-accuracy scaling holds with model size, normalization, and conditioning fixed. The quantum acquisition confidence is probabilistic; accepted rational optimization certificates are deterministic.

**Source pointers.** [Anchor §§2.3–2.4 for the evaluated target][anchor]; [bridge to a certified value](ESTIMATION_BRIDGE.md#turning-acquired-data-into-a-certified-value); [error allocation](../science/circuit_access.md#4-stability-constants-and-confidence-interval); [complete resources](../science/circuit_access.md#5-complete-resource-statement).

### M25 A four-qubit sensing benchmark

**Central question.** What does the algorithmic task look like for an explicitly specified imperfect probe?

**Understanding target.** Follow the parity-preserving preparation, orthogonal error branch, weights $`3/4,1/4`$, and two local Pauli generators. Derive the small factor-label circuit and its constant normalization overhead. Interpret the exact Holevo enclosure near $`2.25076`$, its separation from the SLD cost near $`2.23336`$, and the certified Fisher floor. Distinguish exact full-model witnesses from the numerical sixteen-to-six-dimensional compression comparison.

**Dependencies and story role.** [M05](#m05-incompatibility-and-the-holevo-objective), [M07](#m07-purifications-tangents-and-occupied-support), [M14](#m14-the-label-state-and-its-failure-symbol), [M24](#m24-the-end-to-end-confidence-interval). This grounds the complete procedure in a reproducible physical model.

**Conceptual connection.** Symmetry certifies the probe rank before acquisition. The example establishes a nonzero incompatibility penalty and valid certificates, while remaining classically easy and representing preparation error before sensing.

**Source pointers.** [Anchor §3.1.1 for the simpler incompatibility example][anchor]; [assumptions in the bridge](ESTIMATION_BRIDGE.md#how-the-sensing-example-fits-the-assumptions); [circuit family and access](../science/applications_and_limits.md#2-a-circuit-family-with-certified-rank-and-explicitly-constructed-access); [four-qubit evidence](../verification/README.md#four-qubit-circuit).

### M26 Baselines and oracle precision limits

**Central question.** Against which alternatives is the algorithm's advantage meaningful?

**Understanding target.** Compare the same required overlaps and label state under coherent access, optimal copy access, direct amplitude estimation, and classical circuit contraction. Include analytic special cases, tensor-network, stabilizer, and Gaussian methods when applicable. Distinguish the fixed-conditioning $`\Omega(1/\varepsilon)`$ coherent-query lower bound from the $`\Omega(1/\varepsilon^2)`$ purified-copy lower bound and from the separate worst-case circuit-description hardness witness.

**Dependencies and story role.** [M09](#m09-coherent-access-and-charged-resources), [M18](#m18-normalization-and-query-optimality), [M24](#m24-the-end-to-end-confidence-interval), [M25](#m25-a-four-qubit-sensing-benchmark). This locates the result among matched alternatives.

**Conceptual connection.** An opaque oracle hides information that a short transparent formula may reveal immediately. A query lower bound does not make every circuit hard, and a classically hard witness can still have a direct quantum method stronger than the general pipeline.

**Source pointers.** [Anchor §§2.4 and 3.1 for numerical and analytic evaluation routes][anchor]; [matched baselines](../science/applications_and_limits.md#3-the-strongest-matched-baselines); [circuit-description hardness](../science/applications_and_limits.md#short-circuit-descriptions-can-still-be-hard); [fixed-conditioning precision lower bounds](../science/applications_and_limits.md#precision-lower-bounds-with-conditioning-fixed).

### M27 Why occupied-state conditioning matters

**Central question.** Can fixed rank and a well-conditioned Fisher matrix replace a lower bound on occupied state weights?

**Understanding target.** Follow the rank-three rare-branch family with bounded fixed generators and uniformly positive Fisher information. Its joint factor changes scale as $`q^{3/4}`$ while its Holevo-value separation scales as $`q^{1/2}`$. Understand why this excludes a uniform state-weight-independent Lipschitz bound. Close global preparation-unitary extensions turn the same family into an $`\Omega(\varepsilon^{-3/2})`$ oracle-query obstruction.

**Dependencies and story role.** [M11](#m11-rank-and-conditioning-as-separate-promises), [M21](#m21-stability-and-fisher-promise-transfer), [M26](#m26-baselines-and-oracle-precision-limits). This identifies a substantive limit of the conditioned algorithm.

**Conceptual connection.** Close prepared states alone do not prove close preparation oracles. Separate access to a rare-branch rotation would strengthen the interface and invalidate that lower-bound argument.

**Source pointers.** [Anchor §2.6 for the separation between Fisher and incompatibility information][anchor]; [rare-branch model](../science/conditioning_limits.md#1-a-rank-three-two-parameter-model); [singular scaling](../science/conditioning_limits.md#3-the-singular-scaling); [valid global-oracle lower bound](../science/conditioning_limits.md#4-a-valid-opaque-oracle-lower-bound). No matching general state-gap-free upper theorem is asserted.

### M28 Approximate low rank with a physical certificate

**Central question.** When may a small preparation tail be discarded while retaining a valid precision statement?

**Understanding target.** Use the two-sided certificate for a mixture prepared before a common bounded unitary encoding. Track the discarded mass, retained occupied-weight floor, generator norms, parameter count, and retained Holevo value. Understand the repair of local unbiasedness in both directions and the sharp square-root sensitivity to tail mass. For fixed other parameters, additive-$`\varepsilon`$ control generally requires tail mass of order $`\varepsilon^2`$.

**Dependencies and story role.** [M02](#m02-local-unbiasedness-and-estimation-risk), [M07](#m07-purifications-tangents-and-occupied-support), [M21](#m21-stability-and-fisher-promise-transfer), [M27](#m27-why-occupied-state-conditioning-matters). This gives a qualified extension to imperfect preparations.

**Conceptual connection.** Trace closeness without derivative control is insufficient. Noise during sensing, parameter-dependent mixture weights, and extensive generators require additional accounting; many independent errors can also make the retained branch count large.

**Source pointers.** [Anchor §1.2 for the complete statistical family and §2.3 for its derivative constraints][anchor]; [preparation-tail theorem and proof](../science/applications_and_limits.md#5-a-rigorous-preparation-tail-certificate); [physical limits](../science/applications_and_limits.md#physical-limits-of-5); [conditioning implications](../science/conditioning_limits.md#5-consequence-for-useful-precision-computation).

### M29 Attribution and the strength of the evidence

**Central question.** Which parts are established ingredients, which are proved here, and what do the finite checks establish?

**Understanding target.** Attribute the Holevo framework, numerical SDP, purification-based reduction, coherent Frobenius tomography, amplitude primitives, and adaptive spectral geometry to their prior sources. Identify the project's composition theorem, explicit stability-and-certificate chain, and application/obstruction statements. Separate analytic proof, imported theorem, exact rational verification, floating-point consistency, and an unperformed quantum implementation.

**Dependencies and story role.** [M18](#m18-normalization-and-query-optimality), [M23](#m23-primal-and-dual-optimization-certificates), [M25](#m25-a-four-qubit-sensing-benchmark), [M26](#m26-baselines-and-oracle-precision-limits), [M27](#m27-why-occupied-state-conditioning-matters), [M28](#m28-approximate-low-rank-with-a-physical-certificate). This assembles the evidential basis of the contribution.

**Conceptual connection.** A bounded primary-literature audit supports a careful comparison; it cannot certify universal priority. Finite examples can expose mistakes without proving a general theorem or demonstrating a hardware crossover.

**Source pointers.** [Anchor §2.4 for inherited optimization and §§2.10–2.11 for estimation background][anchor]; [research dependencies beyond the anchor](ESTIMATION_BRIDGE.md#research-dependencies-beyond-the-anchor); [related work and contribution boundary](../references/RELATED_WORK.md); [verification scope](../verification/README.md#interpret-the-evidence-at-its-stated-scope).

### M30 The complete scientific story

**Central question.** What useful capability has been established, and exactly what makes it possible?

**Understanding target.** Join the full chain: explicitly prepared local factors give a small label state; arbitrary-purifier Bures tomography learns it with inverse-accuracy query scaling; physical fitting, conditioned stability, and rational optimization yield a certified local Holevo interval. Explain why the tomography theorem has no eigenvalue floor while the value theorem retains occupied-state and Fisher floors, and why singular tangent Grams remain allowed. Relate the example, strongest baselines, and obstructions to the claimed capability.

**Dependencies and story role.** [M12](#m12-the-complete-computational-question), [M24](#m24-the-end-to-end-confidence-interval), [M25](#m25-a-four-qubit-sensing-benchmark), [M26](#m26-baselines-and-oracle-precision-limits), [M27](#m27-why-occupied-state-conditioning-matters), [M28](#m28-approximate-low-rank-with-a-physical-certificate), [M29](#m29-attribution-and-the-strength-of-the-evidence). This integrates Phase II and the map as a whole.

**Conceptual connection.** The result evaluates the local multiparameter precision cost of a specified coherent circuit model, including incompatibility. Its usefulness depends on paid preparation and derivative access, controlled conditioning, and comparison with the best model-specific alternative. Attainment, optimal-measurement construction, generic noise robustness, and practical runtime advantage remain separate claims.

**Source pointers.** [Anchor §§1.2, 2.3–2.4, 2.10–2.11, and 4.5–4.6 for the estimation meaning and qualifications][anchor]; [project model](../science/model.md); [complete circuit algorithm](../science/circuit_access.md); [applications and limits](../science/applications_and_limits.md); [precise contribution boundary](../references/RELATED_WORK.md#5-contribution-boundary).

[anchor]: https://arxiv.org/html/2001.11742v2
