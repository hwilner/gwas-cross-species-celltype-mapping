# GWAS Cross-Species Cell-Type Mapping (Paper 2)

This independent research repository plans and tracks a cross-species (mouse↔human) heritability mapping analysis. It provides data-free alignment utilities for transparent review and extension.

## Series position

This is **Paper 2** of the GWAS × brain cell-type series (4 papers). **It builds on Paper 1** (`gwas-brain-celltype-heritability`), reusing its harmonized GWAS inputs and human cell-type gene programs.

## Research plan

| Planned work | Expected outcome |
|---|---|
| Mouse cell-type gene programs (ABC Atlas mouse, ~4M cells, 1,201 supertypes) | Versioned mouse marker sets. |
| Cross-species ortholog/homology mapping of gene programs | Aligned mouse–human program pairs with confidence scores. |
| Heritability enrichment in mouse programs | Mouse-side heritability atlas. |
| Cross-species concordance analysis | Map of conserved vs. human-specific disease-cell associations. |

**Current status:** planning stage; no analysis has been run.

## What is included

| Path | Contents |
|---|---|
| `src/` | In-memory ortholog-mapping and concordance utilities. |
| `tests/` | Synthetic tests for alignment statistics. |
| `docs/` | Research status, methods scope, and contribution guidance. |

## Use and validation

```bash
python -m pytest -q
```

## Keywords

GWAS, cross-species, Allen Brain Atlas, single-cell, orthology, Alzheimer's disease, Parkinson's disease, reproducible research.

## Documentation

- [Introduction for new readers](docs/INTRODUCTION.md)
