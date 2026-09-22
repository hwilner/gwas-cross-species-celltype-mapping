# Introduction — Paper 2: Cross-Species Conservation of Brain Cell-Type Heritability

**Series note:** This is **Paper 2 of 4** in the GWAS × brain cell-type series. **It builds directly on Paper 1** (`gwas-brain-celltype-heritability`): the harmonized GWAS summary statistics, human cell-type gene programs, and enrichment pipeline are inputs here. Papers 3 and 4 reuse this paper's concordance framework.

![Concept figure: human and mouse cell-type gene programs are aligned through ortholog mapping and a SAMap homology graph, then the identical heritability-enrichment test is run on both sides to separate concordant cell types from human-specific signals.](figures/concept_figure.png)

*Concept figure: the same enrichment engine from Paper 1 is applied to human cell-type programs and to ortholog-projected mouse programs, and per-species results are compared through a homology-confidence matrix. (Mermaid-renderable version: [figures/concept_figure.md](figures/concept_figure.md).)*

## Background

Most preclinical neuroscience runs on model organisms — above all the mouse. Yet the translational record of neurodegenerative and neuropsychiatric drug development is poor, and a plausible contributor is cellular mismatch: the cell types that carry human genetic risk may not correspond one-to-one to the cell types manipulated in animal models. Human genetics now tells us *where* risk concentrates in the human brain. Partitioned heritability with S-LDSC [1], gene-level enrichment with MAGMA [2], and enhancer-linked disease programs with sc-linker [3] have each shown that common-variant risk for brain traits is not uniformly distributed but concentrated in specific cell populations [4,5]. Paper 1 of this series maps that concentration across the Allen Brain Cell (ABC) Atlas human whole-brain taxonomy — ~3 million nuclei and more than 3,000 cell types spanning the telencephalon to the spinal cord [6].

The missing step is cross-species. Whether the human heritability signal is *conserved* in the homologous mouse (and non-human primate) cell types is almost entirely untested at modern atlas resolution. The question matters in both directions. Where signal is conserved, mouse models of those cell types are validated as genetically faithful systems for mechanistic and therapeutic work. Where signal is human-specific — as might be expected for traits involving human-expanded cortical populations or human-specific regulatory evolution — mouse experiments on the nominally "same" cell type may be systematically misleading.

Two developments make a rigorous test feasible now. First, paired whole-brain atlases exist in matched taxonomic frameworks: the ABC Atlas mouse whole-brain dataset (~4 million cells, 1,201 supertypes) [7] and the human counterpart [6], with additional comparative resources for primate cortex, including human/marmoset/mouse motor cortex comparisons [8]. Second, cross-species integration methods have matured: orthology mapping via Ensembl/HCOP, and alignment of cell types across species through meta-analysis and self-assembling manifold approaches such as SAMap, which maps single-cell atlases onto each other using gene-graph homology [9]. Cell types in these atlases are defined transcriptionally, so a mouse-human homology graph derived from ortholog expression can be built *before* any genetics, and genetic enrichment can then be scored on both sides of each homology edge.

The statistical machinery transfers directly. For each human GWAS, Paper 1 supplies harmonized summary statistics and an S-LDSC/MAGMA/sc-linker pipeline driven by cell-type-specific gene sets [1,2,3]. Because GWAS variants are human variants, "mouse enrichment" cannot be computed from mouse sequence; instead, conservation is tested at the *gene program* level: are the human orthologs of mouse-cell-type-specific genes enriched for human trait heritability? And conversely, do human-cell-type programs score against mouse-cell-type specificity? This ortholog-projected design has precedent in prior cross-species enrichment work [4,5], but it has never been applied systematically across a whole-brain, multi-thousand-type hierarchy in two species with matched taxonomies.

## Prior work and gap

Skene et al. showed that schizophrenia risk maps to pyramidal cells, medium spiny neurons, and interneurons using human and mouse single-cell references, implicitly assuming substantial cross-species concordance [5]. Bryois et al. mapped brain-trait heritability to human cell types and noted broad agreement with mouse-based mappings [4]. Comparative cellular analysis of motor cortex in human, marmoset, and mouse demonstrated both deep conservation of major cell classes and pronounced species divergence at finer granularity — especially in non-neuronal cells and human-specific transcriptomic features [8]. Whole-brain mouse and human taxonomies now exist at near-matched resolution [6,7].

The gap: no study has quantified, cell type by cell type across a whole-brain hierarchy, whether the cell types carrying *human genetic risk* for AD, PD, and related traits are the same types — by transcriptomic homology — in mouse. Prior enrichment studies used coarse references and single species; comparative atlases were not connected to GWAS heritability. The translational consequences (which mouse cell types are valid models of human genetic risk) remain anecdotal.

## Research questions

1. For each trait, do the human orthologs of mouse cell-type-specific genes replicate the heritability enrichments found in the homologous human cell types?
2. At which taxonomy level does human-mouse concordance break down — do classes and subclasses agree while fine types diverge?
3. Which disease-relevant cell types (e.g., AD-associated microglial states, midbrain dopaminergic neurons) are human-specific in their heritability signal?
4. Does marmoset bridge human-mouse divergence for cortical cell types, and does cross-species concordance correlate with sequence/regulatory conservation?

## Data

| Resource | Content | Scale | Access |
|---|---|---|---|
| Paper 1 outputs | Harmonized GWAS, human programs, enrichment pipeline | This series | Internal |
| ABC Atlas human whole brain [6] | snRNA-seq taxonomy | ~3M nuclei, >3,000 types | Open |
| ABC Atlas mouse whole brain [7] | snRNA-seq/spatial taxonomy | ~4M cells, 1,201 supertypes | Open |
| Human/marmoset/mouse motor cortex [8] | Comparative snRNA-seq/snATAC-seq | 3 species | Open |
| Orthology resources | Ensembl Compara / HCOP one-to-one orthologs | Genome-wide | Open |
| Human GWAS summary statistics (AD, PD, SCZ, MDD) | As in Paper 1 [10–13] | Large meta-analyses | Open (GWAS Catalog / PGC) |

## Methods

One-to-one human-mouse(-marmoset) ortholog tables are built from Ensembl Compara and HCOP. For each mouse cell type at each taxonomy level, specificity scores are computed exactly as in Paper 1 and converted to human-gene programs via orthologs. The identical S-LDSC [1] (specifically-expressed-gene mode), MAGMA gene-set [2], and sc-linker-style [3] tests are run against the same GWAS inputs. In parallel, SAMap-style graph alignment [9] produces a homology-confidence matrix between human and mouse cell types. Concordance is quantified per trait as (i) rank correlation of enrichment z-scores across homologous pairs, weighted by homology confidence, and (ii) replication rate of significant human enrichments in mouse ortholog programs (FDR-controlled). Human-specific signals are flagged where a strong human enrichment has no significant mouse-homolog counterpart despite adequate statistical power. Analyses are run across taxonomy levels to localize the divergence point.

## Expected contributions

1. The first whole-hierarchy human-mouse (and marmoset-informed) concordance map of disease heritability across brain cell types.
2. A principled list of mouse cell types that are genetically valid vs. potentially misleading models of human risk for AD/PD and psychiatric traits.
3. A reusable concordance framework and ortholog-projected enrichment pipeline consumed by Papers 3 and 4.
4. Quantification of where in the taxonomy cross-species signal degrades, informing translational study design.

## Scope and boundary

Conservation is assessed statistically at the gene-program level; functional validation in animals is out of scope. GWAS are human and European-ancestry. Ortholog projection loses species-specific genes and regulatory divergence information, which bounds sensitivity; results should be read as tests of conserved *transcriptional* programs, not conserved enhancer sequence.

## References

1. Finucane HK, Bulik-Sullivan B, Gusev A, et al. Partitioning heritability by functional annotation using genome-wide association summary statistics. *Nature Genetics* 47:1228–1235 (2015). doi:10.1038/ng.3404
2. de Leeuw CA, Mooij JM, Heskes T, Posthuma D. MAGMA: generalized gene-set analysis of GWAS data. *PLoS Computational Biology* 11:e1004219 (2015). doi:10.1371/journal.pcbi.1004219
3. Jagadeesh KA, Dey KK, Montoro DT, et al. Identifying disease-critical cell types and cellular processes by integrating single-cell RNA-sequencing and human genetics. *Nature Genetics* 54:1479–1492 (2022). doi:10.1038/s41588-022-01187-9
4. Bryois J, Skene NG, Hansen TF, et al. Genetic identification of cell types underlying brain complex traits yields insights into the etiology of Parkinson's disease. *Nature Genetics* 52:482–493 (2020). doi:10.1038/s41588-020-0610-9
5. Skene NG, Bryois J, Bakken TE, et al. Genetic identification of brain cell types underlying schizophrenia. *Nature Genetics* 50:825–833 (2018). doi:10.1038/s41588-018-0129-5
6. Siletti K, Hodge R, Mossi Albiach A, et al. Transcriptomic diversity of cell types across the adult human brain. *Science* 382:eadd7046 (2023). doi:10.1126/science.add7046
7. Yao Z, van Velthoven CTJ, Kunst M, et al. A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain. *Nature* 624:317–332 (2023). doi:10.1038/s41586-023-06812-z
8. Bakken TE, Jorstad NL, Hu Q, et al. Comparative cellular analysis of motor cortex in human, marmoset and mouse. *Nature* 598:111–119 (2021). doi:10.1038/s41586-021-03465-8
9. Tarashansky AJ, Musser JM, Khariton M, et al. Mapping single-cell atlases throughout Metazoa unravels cell type evolution. *eLife* 10:e66747 (2021). doi:10.7554/eLife.66747
10. Bellenguez C, Küçükali F, Jansen IE, et al. New insights into the genetic etiology of Alzheimer's disease and related dementias. *Nature Genetics* 54:412–436 (2022). doi:10.1038/s41588-022-01024-9
11. Nalls MA, Blauwendraat C, Vallerga CL, et al. Identification of novel risk loci, causal insights, and heritable risk for Parkinson's disease: a meta-analysis of genome-wide association studies. *The Lancet Neurology* 18:1091–1102 (2019). doi:10.1016/S1474-4422(19)30320-5
12. Trubetskoy V, Pardiñas AF, Qi T, et al. Mapping genomic loci implicates genes and synaptic biology in schizophrenia. *Nature* 604:502–508 (2022). doi:10.1038/s41586-022-04434-5
13. Wray NR, Ripke S, Mattheisen M, et al. Genome-wide association analyses identify 44 risk variants and refine the genetic architecture of major depression. *Nature Genetics* 50:668–681 (2018). doi:10.1038/s41588-018-0090-3
14. Finucane HK, Reshef YA, Anttila V, et al. Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types. *Nature Genetics* 50:621–629 (2018). doi:10.1038/s41588-018-0081-4
