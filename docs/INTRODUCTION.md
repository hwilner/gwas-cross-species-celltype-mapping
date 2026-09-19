# Introduction — Paper 2: Cross-Species (Mouse↔Human) Heritability Mapping

**Series note:** This is **Paper 2 of 4** in the GWAS × brain cell-type series. **It builds directly on Paper 1** (`gwas-brain-celltype-heritability`): the harmonized GWAS summary statistics, human cell-type gene programs, and enrichment pipeline are inputs here. Papers 3 and 4 reuse this paper's concordance framework.

## Background

Mouse models are central to neurodegeneration research, but whether mouse cell types implicated by models correspond to the human cell types carrying GWAS heritability is largely untested at modern atlas resolution. The paired Allen mouse (~4M cells) and human (~3M nuclei) whole-brain taxonomies make a direct test possible.

## Research questions

1. Do AD/PD heritability enrichments replicate in homologous mouse cell types?
2. Which disease-cell associations are human-specific (caution flags for mouse models)?
3. Where do mouse AD-model cell states map onto human disease enrichments?

## Data

| Resource | Scale | Access |
|---|---|---|
| Paper 1 GWAS inputs + human programs | internal | This series |
| Allen ABC Atlas mouse whole brain | ~4M cells | Open (AWS) |
| Mouse AD-model scRNA-seq | various | Open (GEO) |

## Methods

Ortholog mapping (HCOP/ensembl), gene-program alignment (SAMap-style), parallel S-LDSC runs, concordance statistics with homology-confidence weighting.

## Expected contributions

- The first whole-hierarchy cross-species heritability concordance map for AD/PD.
- Practical guidance on which mouse cell types model human genetic risk.

## Scope and boundary

Planning, software, and synthetic tests live here. Concordance is statistical; translational claims require functional validation outside this repo's scope.
