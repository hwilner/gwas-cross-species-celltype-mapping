# GWAS Cross-Species Cell-Type Mapping

This independent research repository plans and tracks a cross-species (mouse↔human) heritability mapping analysis. It provides data-free alignment utilities for transparent review and extension.

## Introduction

Which disease-associated cell types are conserved between mouse and human, and which are human-specific? This repository addresses that question by mapping GWAS heritability enrichments onto mouse brain cell types and aligning the resulting gene programs with the human atlas — a prerequisite for interpreting mouse models of human disease. Full text: [docs/INTRODUCTION.md](docs/INTRODUCTION.md).

## Extended introduction

Background and related work in depth: [docs/EXTENDED_INTRODUCTION.md](docs/EXTENDED_INTRODUCTION.md). This repository **builds on the companion repository** `gwas-brain-celltype-heritability`, reusing its harmonized GWAS inputs and human cell-type gene programs; the dependence is stated explicitly wherever results rely on its outputs.

## Materials and methods

Mouse cell-type gene programs (ABC Atlas mouse, ~4M cells, 1,201 supertypes) are built with the frozen marker-selection rules of `gwas-brain-celltype-heritability` (no re-tuning); mouse–human orthology uses a versioned HCOP/Ensembl table. Programs are aligned with the homology-confidence-weighted statistic implemented in `src/alignment.py` (SAMap-style refinement; documented `high`/`medium`/`low` confidence tiers), and concordance is tested against a pre-registered null model. All code is data-free and validated on synthetic cases before any real run. Full protocol: [docs/METHODS.md](docs/METHODS.md); run plans and pre-registrations are tracked as issues (#6, #14).

## Results

**Current status:** planning stage; no empirical analysis has been run. The results package is pre-registered as captioned paper-style tables in [docs/RESULTS.md](docs/RESULTS.md) — Table 1 (mouse-side enrichment atlas), Table 2 (alignment summary with confidence tiers), Table 3 (concordance vs. null), Table 4 (conserved vs. human-specific catalogue), Table 5 (mouse AD-model state mapping) — each linked to the issue that authorizes its computation. The alignment utilities themselves are implemented and validated: 35/35 synthetic tests pass (`python -m pytest -q`; closed issue #5).

## Research plan

| Planned work | Expected outcome |
|---|---|
| Mouse cell-type gene programs (ABC Atlas mouse, ~4M cells, 1,201 supertypes) | Versioned mouse marker sets. |
| Cross-species ortholog/homology mapping of gene programs | Aligned mouse–human program pairs with confidence scores. |
| Heritability enrichment in mouse programs | Mouse-side heritability atlas. |
| Cross-species concordance analysis | Map of conserved vs. human-specific disease-cell associations. |

## What is included

| Path | Contents |
|---|---|
| `src/` | In-memory ortholog-mapping and concordance utilities. |
| `tests/` | Synthetic tests for alignment statistics. |
| `docs/` | Research status, methods scope, results skeleton, and contribution guidance. |

## Use and validation

```bash
python -m pytest -q
```

## Keywords

GWAS, cross-species, Allen Brain Atlas, single-cell, orthology, Alzheimer's disease, Parkinson's disease, reproducible research.

## Documentation

- [Introduction for new readers](docs/INTRODUCTION.md)
- [Extended introduction](docs/EXTENDED_INTRODUCTION.md)
- [Materials and methods](docs/METHODS.md)
- [Results (pre-registered tables)](docs/RESULTS.md)
