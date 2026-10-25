# Quantum Precision Limits

**How can coherent access to a quantum circuit determine a multiparameter
Holevo precision bound without reconstructing its large physical state?**

This project studies Bures tomography from a state-preparation unitary and uses
it to compute certified Holevo bounds for low-rank sensing models. The central
tomography guarantee permits an arbitrarily large purification workspace. The
application learns a small register containing the state-and-derivative Gram
matrix, then certifies a classical optimization over that model.

| Result | Guarantee | Essential conditions |
|---|---|---|
| Coherent Bures tomography | $`\widetilde\Theta(dr/\varepsilon_B)`$ worst-case preparation/inverse queries | A rank-$`r`$, $`d`$-dimensional reduced state; no occupied-eigenvalue floor; arbitrary supplied purification |
| Certified Holevo computation | $`\widetilde O(1/\varepsilon)`$ precision dependence at fixed model size and conditioning | Explicit coherent factor/derivative access; occupied-state and Fisher lower bounds; no tangent-Gram eigenvalue floor |
| State-weight obstruction | $`\Omega(\varepsilon^{-3/2})`$ queries on a rank-three family without an occupied-state floor | Opaque global preparation oracles; bounded fixed generators and uniformly conditioned Fisher information |

The matching tomography lower bound is a worst-case oracle statement for
sufficiently small constant-bounded accuracy and $`d\ge2`$.

The tomography upper bound combines published coherent Frobenius tomography
with adaptive spectral peeling. Its lower bound follows from published
trace-norm tomography. The arbitrary-purifier Bures guarantee and its certified
precision application are developed here; the [related-work account](references/RELATED_WORK.md)
identifies the inherited ingredients and closest competing methods.

## Why Bures accuracy matters

A local model of rank $`r`$ with $`p`$ parameters has factors
$`V=[\Psi,\Phi_1,\ldots,\Phi_p]`$. Its Holevo value depends only on
$`K=V^\dagger V`$, a matrix of order $`m=r(p+1)`$. A controlled preparation
encodes $`K^{\mathsf T}`$ into a label register of at most $`m+1`$ dimensions.
The physical sensor becomes its purification.

Near a singular tangent Gram, entrywise accuracy gives only square-root
control of the factors. Learning directly in Bures distance avoids that generic
accuracy loss. The algorithm returns a classical confidence interval for the
local Holevo value. Constructing an attaining measurement is a separate task.

## Read the science

Start with the [question, model, and exact scope](science/model.md).

| To inspect | Read |
|---|---|
| Optimal Bures query scaling, conditioning cost, and finite gates | [Coherent Bures tomography](science/bures_tomography.md) |
| Small tangent-support model and quantitative perturbation control | [Compression and stability](science/holevo_compression_and_stability.md) |
| Physical Gram fitting and checkable classical value certificates | [Certified optimization](science/certified_optimization.md) |
| Input circuits, normalization, confidence budget, and full resources | [Circuit access and algorithm](science/circuit_access.md) |
| Circuit examples, matched baselines, and approximate-low-rank bounds | [Applications and limits](science/applications_and_limits.md) |
| Why Fisher conditioning alone is insufficient | [Rare-branch obstruction](science/conditioning_limits.md) |
| Primary sources and novelty boundaries | [Related work](references/RELATED_WORK.md) and [bibliography](references/references.bib) |

For specialists, the proof route is tomography → compression/stability →
certified optimization → circuit access. The applications explain the physical
meaning and comparisons. The [LLM reading guide](llms.txt) maps questions to
their authoritative files.

## Study the project

The single educational anchor is Rafał Demkowicz-Dobrzański, Wojciech Górecki,
and Mădălin Guţă's [*Multi-parameter estimation beyond Quantum Fisher
Information*](https://arxiv.org/html/2001.11742v2) (2020). It explains the
statistical meaning of simultaneous precision and the Holevo bound.

Start with the [source-to-project bridge](learning/ESTIMATION_BRIDGE.md), which
selects the relevant sections, translates notation and assumptions, and works
through how the estimation problem becomes a certified quantum computation.
The [30-milestone map](learning/MILESTONES.md) organizes the full project into
background/preliminaries and narrative/innovation, with explicit understanding
targets, dependencies, and links to the proofs. The quantum acquisition and
certification arguments are developed in this repository and retain their
primary research dependencies.

## A checkable circuit example

A four-qubit parity-preserving circuit prepares a rank-two probe with occupied
weights $`3/4`$ and $`1/4`$. It senses two local Pauli quadratures. Exact rational
primal and dual witnesses give

```math
2.250762525\le C_H(I_2)\le2.250762539,
```

while the SLD benchmark lies between $`2.233360231`$ and $`2.233360232`$.
The Fisher matrix is certified to satisfy $`J\succeq(4/5)I_2`$. Thus the
example has a nonzero incompatibility penalty despite good conditioning.
A numerical comparison also agrees after compressing the physical model from
16 dimensions to 6. The [full circuit and certificate](science/applications_and_limits.md)
are reproducible.

## Reproduce the finite checks

With the dependencies in [verification/requirements.txt](verification/requirements.txt), run:

```bash
python3 verification/run_all.py
```

The suite checks spectral-peeling geometry, exact physical-fit and Holevo
certificates, and the four-qubit example. Its [verification guide](verification/README.md)
distinguishes exact rational assertions from floating-point comparisons and
explains the committed expected outputs. These checks do not simulate the
coherent tomography circuit or establish hardware runtime advantage.

## Resource and physical boundaries

Preparation circuits, inverses, controls, derivative normalization, small-register
gates, and classical optimization are charged. No QRAM or free loading of
exponentially many amplitudes is assumed. The Holevo application requires a
full-column-rank factor with an explicitly small environment, even though the
general Bures learner allows an arbitrarily large purifier.

Classically contractible circuits and direct overlap methods remain strong
baselines. The exact closed-form hardness witness has its own simpler quantum
algorithm. The general result establishes a query guarantee for the stated
input model; it does not establish a practical hardware crossover.

Approximate low rank is controlled only with the stated common-encoding and
bounded-generator assumptions. A small discarded state weight can change the
precision benchmark at square-root order.

## Purpose and contact

This repository serves as a record of the work and a guide for the author’s self-directed learning. For discussion or potential collaboration, please contact Ruge Lin at [gogoko699@gmail.com](mailto:gogoko699@gmail.com).
