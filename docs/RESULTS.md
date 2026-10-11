# Results

**Status: pre-registered skeleton.** This repository is at the planning stage;
no empirical analysis has been run yet (see the README "Current status" line
and issues #1–#4). The tables below fix, in advance, the exact structure of
the results package so that table generation (issue #15) is a fill-in
exercise rather than a post-hoc design decision. Every table lists the issue
that authorizes its computation.

All statistics will be produced by committed code from versioned inputs only;
tables marked *pending* contain no numbers by design.

---

**Table 1. Mouse-side heritability enrichment atlas (planned).** One row per
mouse supertype × GWAS trait pair tested. Enrichment statistics and FDR
procedure are identical by construction to the frozen
`gwas-brain-celltype-heritability` pipeline (no re-tuning; see issue #6).
*Status: pending — blocked by issues #6, #9 (mouse enrichment run) and #10
(mouse supertype gene programs).*

| Supertype | GWAS trait | Enrichment | SE | p | q (FDR) |
|---|---|---|---|---|---|
| — | — | — | — | — | — |

**Table 2. Cross-species program alignment summary (planned).** One row per
aligned mouse–human program pair, reporting the mutual-best-rescaled
alignment score and the contribution-weighted homology confidence with its
tier (`high` / `medium` / `low`; thresholds pre-registered in
`src/alignment.py` as `CONFIDENCE_TIER_THRESHOLDS`). Alignment code is
implemented and synthetically validated (closed issue #5); real-data scores
await versioned inputs. *Status: pending — blocked by issues #7 (ortholog
table) and #10 (mouse programs).*

| Human program | Mouse program | Alignment score | Confidence | Confidence tier |
|---|---|---|---|---|
| — | — | — | — | — |

**Table 3. Concordance metrics against the pre-registered null (planned).**
Observed concordance between aligned mouse and human programs versus the
null distribution for chance-level concordance. Metrics and null family are
fixed before any numbers are computed (issue #14) and implemented with unit
tests (issue #12). *Status: pending — blocked by issues #14, #12, and the
mouse enrichment run (#2).*

| Metric | Observed | Null mean | Null SD | Null exceedance (p) |
|---|---|---|---|---|
| — | — | — | — | — |

**Table 4. Conserved vs. human-specific disease–cell associations
(planned).** Catalogue of associations present in the human atlas that are
conserved, absent, or below null expectation in the aligned mouse programs.
The classification rule is fixed by the pre-registered null (issue #14).
*Status: pending — blocked by issues #11, #12.*

| GWAS trait | Human cell type | Aligned mouse program | Classification (conserved / absent / below null) |
|---|---|---|---|
| — | — | — | — |

**Table 5. Mouse AD-model state mapping onto human enrichments (planned).**
Overlap between mouse AD-model cell-state programs and the human enrichment
results, projected through the alignment. Statistical mapping only; no
translational efficacy claims (issue boundary). *Status: pending — blocked
by issues #12, #13.*

| Mouse AD-model state | Best-aligned human program | Alignment score | Overlap statistic | Null exceedance (p) |
|---|---|---|---|---|
| — | — | — | — | — |

---

## Provenance and versioning

- Every committed table must record its generating code version and input
  manifest, per the acceptance criteria of issue #15.
- Figures are tracked separately (issue #16) and will be regenerated
  deterministically from the versioned tables above.
- The clean-room reproducibility re-run (issue #17) verifies that this file
  regenerates identically from the frozen manifest.
