# Concept Figure — Cross-Species Conservation of Brain Cell-Type Heritability

> **Note.** The rendered concept figure is [`01-concept-schematic.png`](01-concept-schematic.png) in this folder — a generated scientific illustration. This file keeps the Mermaid source of the same diagram so it remains editable and re-renderable.

**Caption:** Concept figure — the design of this repository: human GWAS summary statistics meet human and mouse cell-type gene programs (aligned via orthologs and SAMap homology); the identical enrichment engine runs on both species' programs, and concordance analysis separates conserved cell types (valid mouse models) from human-specific signals (translation risk).

```mermaid
flowchart LR
    G["Human GWAS summary statistics<br/>AD / PD / SCZ / MDD (from gwas-brain-celltype-heritability)"]
    subgraph Programs
        H["Human cell-type gene programs<br/>ABC human: ~3M nuclei, >3,000 types"]
        M["Mouse cell-type gene programs<br/>ABC mouse: ~4M cells, 1,201 supertypes"]
    end
    O["Ortholog mapping (Ensembl/HCOP)<br/>+ SAMap homology graph"]
    E["Identical enrichment engine<br/>(S-LDSC / MAGMA / sc-linker)<br/>run on both species' programs"]
    C1["Concordant cell types:<br/>mouse model validated"]
    C2["Human-specific signals:<br/>translation risk"]
    G --> E
    H --> E
    M --> O
    H <--> O
    O --> E
    E --> C1
    E --> C2
```
