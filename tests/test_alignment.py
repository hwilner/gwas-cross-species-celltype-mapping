"""Synthetic tests for the cross-species program alignment utilities.

All tests are fully data-free: human/mouse "programs" are planted synthetic
score vectors and the homology graphs are built by hand, so the expected
alignment outcome is known exactly and can be asserted.
"""

import numpy as np
import pytest

from src.alignment import (
    CONFIDENCE_TIER_THRESHOLDS,
    HomologyGraph,
    compute_alignment,
    project_mouse_programs,
)


def _planted_case(
    n_pairs=4,
    n_driver=8,
    n_extra=8,
    noise=0.05,
    homology_type="ortholog_one2one",
):
    """Build a synthetic case with known one-to-one program partners.

    Each of ``n_pairs`` human programs has a private block of ``n_driver``
    high-scoring genes; the mouse programs mirror the same blocks on the
    corresponding mouse genes, and every gene pair is linked by a homology
    edge of the requested type. Program ``p`` in human must therefore align
    to program ``p`` in mouse. ``n_extra`` additional "housekeeping" gene
    pairs carry the same low ``noise`` score in every program, diluting the
    driver signal the way ubiquitously expressed genes would.

    Args:
        n_pairs: Number of planted program pairs.
        n_driver: Number of private driver genes per program.
        n_extra: Number of shared background gene pairs scored ``noise``.
        noise: Background score placed on the shared extra genes.
        homology_type: Edge type used for every planted edge.

    Returns:
        Tuple ``(human_scores, mouse_scores, edges)`` with shapes
        ``(n_pairs, n_pairs * n_driver + n_extra)`` and the edge list.
    """
    n_block = n_pairs * n_driver
    n_genes = n_block + n_extra
    human = np.zeros((n_pairs, n_genes))
    mouse = np.zeros((n_pairs, n_genes))
    for program in range(n_pairs):
        block = slice(program * n_driver, (program + 1) * n_driver)
        human[program, block] = 1.0
        mouse[program, block] = 1.0
    human[:, n_block:] = noise
    mouse[:, n_block:] = noise
    edges = [
        (f"HG{k:02d}", f"MG{k:02d}", homology_type) for k in range(n_genes)
    ]
    return human, mouse, edges


def test_recovers_planted_homologous_programs():
    """Headline criterion: planted homologous programs are recovered.

    For every human program the best-scoring mouse program (and vice versa)
    must be its planted partner, with a wide margin over all off-diagonal
    pairs.
    """
    human, mouse, edges = _planted_case()
    graph = HomologyGraph(edges)
    result = compute_alignment(human, mouse, graph)

    best_mouse = result.scores.argmax(axis=1)
    best_human = result.scores.argmax(axis=0)
    assert list(best_mouse) == [0, 1, 2, 3]
    assert list(best_human) == [0, 1, 2, 3]

    diagonal = np.diag(result.scores)
    off_diagonal = result.scores[~np.eye(4, dtype=bool)].reshape(4, 3)
    assert diagonal.min() >= 0.99  # mutual-best pairs rescale to ~1.0
    assert off_diagonal.max() <= 0.1


def test_alignment_is_deterministic():
    """Identical inputs must produce bit-identical outputs (pure algebra)."""
    human, mouse, edges = _planted_case()
    graph = HomologyGraph(edges)
    first = compute_alignment(human, mouse, graph)
    second = compute_alignment(human, mouse, graph)

    assert np.array_equal(first.scores, second.scores)
    assert np.array_equal(first.confidence, second.confidence)
    assert np.array_equal(first.similarity, second.similarity)
    assert first.confidence_tiers.tolist() == second.confidence_tiers.tolist()


def test_scores_bounded_in_unit_interval():
    """Mutual-best rescaling must keep every score inside [0, 1]."""
    human, mouse, edges = _planted_case(noise=0.2)
    graph = HomologyGraph(edges)
    result = compute_alignment(human, mouse, graph, n_refinements=3)
    assert result.scores.min() >= 0.0
    assert result.scores.max() <= 1.0 + 1e-12


def test_one2one_support_yields_high_confidence_tier():
    """Pairs supported purely by one-to-one edges are tiered ``high``."""
    human, mouse, edges = _planted_case(homology_type="ortholog_one2one")
    graph = HomologyGraph(edges)
    result = compute_alignment(human, mouse, graph)

    diagonal_confidence = np.diag(result.confidence)
    assert diagonal_confidence == pytest.approx(1.0)
    assert all(tier == "high" for tier in np.diag(result.confidence_tiers))


def test_one2many_support_yields_medium_confidence_tier():
    """Pairs supported purely by one-to-many edges are tiered ``medium``."""
    human, mouse, edges = _planted_case(homology_type="ortholog_one2many")
    graph = HomologyGraph(edges)
    result = compute_alignment(human, mouse, graph)

    diagonal_confidence = np.diag(result.confidence)
    assert diagonal_confidence == pytest.approx(0.5)
    assert all(tier == "medium" for tier in np.diag(result.confidence_tiers))


def test_many2many_support_yields_low_confidence_tier():
    """Pairs supported purely by many-to-many edges are tiered ``low``."""
    human, mouse, edges = _planted_case(homology_type="ortholog_many2many")
    graph = HomologyGraph(edges)
    result = compute_alignment(human, mouse, graph)

    diagonal_confidence = np.diag(result.confidence)
    assert diagonal_confidence == pytest.approx(0.25)
    assert all(tier == "low" for tier in np.diag(result.confidence_tiers))


def test_projection_matrix_is_row_stochastic():
    """Each mouse gene distributes exactly one unit of mass to its homologs."""
    _, _, edges = _planted_case()
    graph = HomologyGraph(edges + [("HG00", "MG00", "ortholog_one2many")])
    projection = graph.projection_matrix()
    assert projection.shape == (40, 40)  # 32 block genes + 8 housekeeping extras
    assert projection.sum(axis=1) == pytest.approx(np.ones(40))
    # The duplicated edge keeps only the max weight (1.0), so MG00's row is
    # a single unit on HG00 rather than a split.
    mouse_index = graph.mouse_genes.index("MG00")
    human_index = graph.human_genes.index("HG00")
    assert projection[mouse_index, human_index] == pytest.approx(1.0)


def test_duplicate_edges_collapse_to_max_weight():
    """A duplicated gene pair must keep its most confident weight.

    With edges HG_A--MG_X (one2many, 0.5), HG_B--MG_X (one2one, 1.0) and a
    redundant HG_A--MG_X (one2one, 1.0), max-collapse makes the two weights
    1.0 and 1.0, so MG_X's stochastic row splits mass evenly; keeping the
    first-seen 0.5 instead would split 1/3 vs 2/3.
    """
    graph = HomologyGraph(
        [
            ("HG_A", "MG_X", "ortholog_one2many"),
            ("HG_B", "MG_X", "ortholog_one2one"),
            ("HG_A", "MG_X", "ortholog_one2one"),
        ]
    )
    projection = graph.projection_matrix()
    assert projection[0, 0] == pytest.approx(0.5)
    assert projection[0, 1] == pytest.approx(0.5)


def test_refinement_sharpens_diluted_alignment():
    """Refinement must improve contrast when driver signal is diluted.

    With many shared housekeeping genes (40 extras at 0.3), the plain
    projection yields the correct argmax structure but weak contrast; the
    SAMap-style gene reweighting must strictly increase the diagonal-minus-
    off-diagonal contrast while keeping the planted partners on top.
    """
    human, mouse, edges = _planted_case(n_extra=40, noise=0.3)
    graph = HomologyGraph(edges)

    plain = compute_alignment(human, mouse, graph, n_refinements=0)
    refined = compute_alignment(human, mouse, graph, n_refinements=2)

    def contrast(similarity):
        diag = np.diag(similarity).mean()
        off = similarity[~np.eye(4, dtype=bool)].mean()
        return diag - off

    # Contrast is compared on the cosine similarity (pre-rescaling), which is
    # the quantity the refinement iterations actually reweight.
    assert contrast(refined.similarity) > contrast(plain.similarity)
    assert list(refined.scores.argmax(axis=1)) == [0, 1, 2, 3]
    assert np.diag(refined.scores).min() >= 0.99


def test_refinement_with_zero_signal_is_safe():
    """With no positive alignment signal, refinement is a no-op (no NaN)."""
    human = np.zeros((2, 4))
    mouse = np.ones((2, 4))
    edges = [(f"HG{k}", f"MG{k}", "ortholog_one2one") for k in range(4)]
    graph = HomologyGraph(edges)
    result = compute_alignment(human, mouse, graph, n_refinements=2)
    assert np.isfinite(result.scores).all()
    assert (result.scores == 0.0).all()


def test_no_shared_genes_align_to_zero():
    """Mouse mass projected onto unused human genes yields zero alignment."""
    edges = [(f"HX{k}", f"MX{k}", "ortholog_one2one") for k in range(8)]
    graph = HomologyGraph(edges)
    # Human programs use only HX4..HX7; mouse programs score only MX0..MX3,
    # whose mass projects onto HX0..HX3 -- disjoint from the human usage, so
    # no program pair shares any gene mass.
    human_block = np.zeros((4, 8))
    human_block[:, 4:] = 1.0
    mouse = np.zeros((4, 8))
    mouse[:, 0] = 1.0
    result = compute_alignment(human_block, mouse, graph)
    assert (result.scores == 0.0).all()
    assert (result.confidence == 0.0).all()


def test_zero_program_row_aligns_to_nothing():
    """An all-zero human program gets zero scores and zero confidence."""
    human, mouse, edges = _planted_case()
    human[2, :] = 0.0
    graph = HomologyGraph(edges)
    result = compute_alignment(human, mouse, graph)
    assert (result.scores[2, :] == 0.0).all()
    assert (result.confidence[2, :] == 0.0).all()
    assert np.isfinite(result.scores).all()
    assert np.isfinite(result.similarity).all()


def test_confidence_matches_hand_computed_value():
    """Confidence equals the contribution-weighted mean edge weight.

    Hand computation for edges HG_A--MG_X (weight 1.0) and HG_B--MG_X
    (weight 0.5), one mouse program [1.0]: projection P = [2/3, 1/3];
    confidence for human programs [1,0], [0,1], [1,1] is 1.0, 0.5, 5/6.
    """
    graph = HomologyGraph(
        [("HG_A", "MG_X", "ortholog_one2one"), ("HG_B", "MG_X", "ortholog_one2many")]
    )
    human = np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    mouse = np.array([[1.0]])
    result = compute_alignment(human, mouse, graph, n_refinements=0)

    expected = [1.0, 0.5, 5.0 / 6.0]
    assert result.confidence[:, 0] == pytest.approx(expected)
    assert result.confidence_tiers[:, 0].tolist() == ["high", "medium", "medium"]


def test_project_mouse_programs_matches_manual_projection():
    """Projection through the graph must place mass per the edge weights."""
    graph = HomologyGraph(
        [("HG_A", "MG_X", "ortholog_one2one"), ("HG_B", "MG_X", "ortholog_one2many")]
    )
    projected = project_mouse_programs(np.array([[3.0]]), graph)
    assert projected.shape == (1, 2)
    assert projected[0] == pytest.approx([2.0, 1.0])  # 3.0 * [2/3, 1/3]


def test_gene_order_is_sorted_and_deterministic():
    """Graph gene axes follow sorted label order regardless of edge order."""
    shuffled = HomologyGraph(
        [("HG_B", "MG_Y", "ortholog_one2one"), ("HG_A", "MG_X", "ortholog_one2one")]
    )
    assert shuffled.human_genes == ("HG_A", "HG_B")
    assert shuffled.mouse_genes == ("MG_X", "MG_Y")


def test_labels_default_and_passthrough():
    """Default labels are generated; caller labels pass through unchanged."""
    human, mouse, edges = _planted_case()
    graph = HomologyGraph(edges)

    defaulted = compute_alignment(human, mouse, graph)
    assert defaulted.human_labels == ("human_0", "human_1", "human_2", "human_3")
    assert defaulted.mouse_labels == ("mouse_0", "mouse_1", "mouse_2", "mouse_3")

    named = compute_alignment(
        human,
        mouse,
        graph,
        human_labels=["OPC", "Micro", "Astro", "Vasc"],
        mouse_labels=["mOPC", "mMicro", "mAstro", "mVasc"],
    )
    assert named.human_labels == ("OPC", "Micro", "Astro", "Vasc")
    assert named.mouse_labels == ("mOPC", "mMicro", "mAstro", "mVasc")
    assert named.n_refinements == 2


def test_threshold_constants_are_consistent():
    """The documented tier thresholds must bracket the default tier weights."""
    high = CONFIDENCE_TIER_THRESHOLDS["high"]
    medium = CONFIDENCE_TIER_THRESHOLDS["medium"]
    assert 0.0 < medium < high <= 1.0
    # one2one support (1.0) is high, one2many (0.5) is medium, many2many
    # (0.25) is low under the documented thresholds.
    assert 1.0 >= high
    assert medium <= 0.5 < high
    assert 0.25 < medium


# ---------------------------------------------------------------------------
# Input validation: every malformed input must raise a clear ValueError.
# ---------------------------------------------------------------------------


def test_empty_edges_rejected():
    """A homology graph without edges cannot define a projection."""
    with pytest.raises(ValueError, match="at least one homology edge"):
        HomologyGraph([])


def test_unknown_homology_type_rejected():
    """Edges must use a homology type present in the tier-weight mapping."""
    with pytest.raises(ValueError, match="unknown homology type"):
        HomologyGraph([("HG_A", "MG_X", "paralog")])


def test_malformed_edge_rejected():
    """Edges that are not 3-tuples must be rejected with position info."""
    with pytest.raises(ValueError, match="edge #0"):
        HomologyGraph([("HG_A", "MG_X")])


def test_empty_gene_label_rejected():
    """Blank gene labels would make matrix axes uninterpretable."""
    with pytest.raises(ValueError, match="invalid human gene label"):
        HomologyGraph([("", "MG_X", "ortholog_one2one")])
    with pytest.raises(ValueError, match="invalid mouse gene label"):
        HomologyGraph([("HG_A", "  ", "ortholog_one2one")])


@pytest.mark.parametrize("bad_weight", [0.0, -0.1, 1.5, float("nan"), float("inf")])
def test_invalid_custom_tier_weight_rejected(bad_weight):
    """Custom tier weights must be finite and inside (0, 1]."""
    with pytest.raises(ValueError, match="tier weight"):
        HomologyGraph(
            [("HG_A", "MG_X", "custom")], tier_weights={"custom": bad_weight}
        )


def test_empty_custom_tier_weights_rejected():
    """An empty custom tier-weight mapping types no edge at all."""
    with pytest.raises(ValueError, match="at least one homology type"):
        HomologyGraph([("HG_A", "MG_X", "ortholog_one2one")], tier_weights={})


def test_matrix_width_mismatch_rejected():
    """Score matrices must match the graph's gene-axis sizes."""
    graph = HomologyGraph([("HG_A", "MG_X", "ortholog_one2one")])
    with pytest.raises(ValueError, match="columns"):
        compute_alignment(np.ones((2, 2)), np.ones((2, 1)), graph)
    with pytest.raises(ValueError, match="columns"):
        compute_alignment(np.ones((2, 1)), np.ones((2, 3)), graph)


def test_negative_scores_rejected():
    """Signed statistics must be transformed before alignment."""
    graph = HomologyGraph([("HG_A", "MG_X", "ortholog_one2one")])
    with pytest.raises(ValueError, match="non-negative"):
        compute_alignment(np.array([[-0.5]]), np.ones((1, 1)), graph)


def test_non_finite_scores_rejected():
    """NaN or infinite scores would poison the similarity computation."""
    graph = HomologyGraph([("HG_A", "MG_X", "ortholog_one2one")])
    with pytest.raises(ValueError, match="finite"):
        compute_alignment(np.array([[np.nan]]), np.ones((1, 1)), graph)
    with pytest.raises(ValueError, match="finite"):
        compute_alignment(np.ones((1, 1)), np.array([[np.inf]]), graph)


@pytest.mark.parametrize("bad_refinements", [-1, 2.5, True, "2"])
def test_n_refinements_validated(bad_refinements):
    """n_refinements must be a non-negative integer (bools excluded)."""
    graph = HomologyGraph([("HG_A", "MG_X", "ortholog_one2one")])
    with pytest.raises(ValueError, match="n_refinements"):
        compute_alignment(np.ones((1, 1)), np.ones((1, 1)), graph, n_refinements=bad_refinements)


def test_label_errors_rejected():
    """Wrong-length or duplicated labels would make results unaddressable."""
    human, mouse, edges = _planted_case()
    graph = HomologyGraph(edges)
    with pytest.raises(ValueError, match="entries"):
        compute_alignment(human, mouse, graph, human_labels=["only_one"])
    with pytest.raises(ValueError, match="unique"):
        compute_alignment(
            human, mouse, graph, mouse_labels=["dup", "dup", "dup", "dup"]
        )
