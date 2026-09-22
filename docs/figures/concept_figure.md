# Concept Figure — Paper 2: Cross-Species Conservation of Brain Cell-Type Heritability

> **Note:** A raster concept figure (`concept_figure.png`, 1536×1024, FigForge/NeurIPS flat-design style) was generated for this repository. Because the available GitHub tooling can only commit text content (verified: binary payloads are stored corrupted), the faithful figure is provided here as a Mermaid diagram that renders directly on GitHub. If `concept_figure.png` is later uploaded to this directory, it can be embedded as:
>
> `![Concept figure: human and mouse cell-type gene programs are aligned through ortholog mapping and a homology graph, then the identical heritability-enrichment test is run on both sides to flag concordant versus human-specific disease cell types.](figures/concept_figure.png)`

**Caption:** Concept figure — the Paper 2 design: human GWAS summary statistics meet human and mouse cell-type gene programs (aligned via orthologs and SAMap homology); the identical enrichment engine runs on both species' programs, and concordance analysis separates conserved cell types (valid mouse models) from human-specific signals (translation risk).

```mermaid
flowchart LR
    G["Human GWAS summary statistics<br/>AD / PD / SCZ / MDD (from Paper 1)"]
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
