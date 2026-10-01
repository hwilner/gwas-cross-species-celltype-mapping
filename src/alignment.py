"""Data-free cross-species gene-program alignment utilities.

This module implements a SAMap-style, homology-confidence-weighted alignment
between human and mouse cell-type *gene programs* (per-program, per-gene
marker or specificity scores supplied by the caller). Everything operates on
in-memory tables: no data files are read, no network access is performed, and
no real atlas or GWAS data is required. The intended use at the current
planning stage is transparent review and synthetic validation of the
alignment statistic before any real run.

Pipeline (every step is deterministic pure linear algebra):

1. Mouse programs are projected into human gene space through a validated
   homology graph whose edges carry confidence weights derived from the
   homology type (see ``TIER_WEIGHTS``).
2. A cosine-similarity matrix between human programs and projected mouse
   programs is computed.
3. Optional SAMap-style refinement iterations re-weight human genes by how
   strongly they drive the currently aligned program pairs (the program
   graph refines the gene graph, mirroring the alternating structure of
   SAMap), sharpening alignments that are driven by consistently used genes.
4. Final alignment scores are mutual-best rescaled into the unit interval:
   a pair scores 1.0 exactly when the two programs are each other's best
   match.

Every pair also receives a *confidence* value: the contribution-weighted
mean homology-edge weight supporting that pair, tiered as ``high`` /
``medium`` / ``low`` via ``CONFIDENCE_TIER_THRESHOLDS``.

Documented confidence tiers:

* ``high``   -- confidence >= 0.9: pair support comes essentially entirely
  from one-to-one ortholog edges.
* ``medium`` -- 0.4 <= confidence < 0.9: mixed support, e.g. dominated by
  one-to-many ortholog edges.
* ``low``    -- confidence < 0.4: support dominated by many-to-many edges;
  pairs in this tier should not drive conclusions.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Mapping, Optional, Sequence, Tuple

import numpy as np

# Confidence weight assigned to each homology edge type. The defaults mirror
# the usual Ensembl Compara / HCOP-style homology classes: one-to-one
# orthologs carry full confidence, one-to-many half, many-to-many a quarter.
# Rationale: projection mass should flow preferentially through unambiguous
# gene correspondences, and the same weights later define the per-pair
# confidence, so ambiguous support is visible downstream instead of hidden.
TIER_WEIGHTS = {
    "ortholog_one2one": 1.0,
    "ortholog_one2many": 0.5,
    "ortholog_many2many": 0.25,
}

# Thresholds (inclusive lower bounds) mapping a pair's confidence value onto
# the documented ``high`` / ``medium`` / ``low`` tiers. 0.9 requires almost
# pure one-to-one support; 0.4 separates one-to-many-dominated support from
# many-to-many-dominated support.
CONFIDENCE_TIER_THRESHOLDS = {"high": 0.9, "medium": 0.4}

# Tier label assigned below the ``medium`` threshold.
_LOW_TIER_LABEL = "low"


class HomologyGraph:
    """Validated weighted gene--gene homology graph between two species.

    The graph stores homology edges ``(human_gene, mouse_gene, type)`` and
    exposes a row-stochastic projection matrix that moves mouse program
    scores into human gene space, distributing each mouse gene's mass across
    its human homologs in proportion to the edge confidence weights.

    Attributes:
        human_genes: Sorted tuple of unique human gene labels (deterministic
            column order of the projection matrix).
        mouse_genes: Sorted tuple of unique mouse gene labels (deterministic
            row order of the projection matrix).
        tier_weights: Mapping from homology type to its confidence weight.
    """

    def __init__(
        self,
        edges: Iterable[Tuple[str, str, str]],
        tier_weights: Optional[Mapping[str, float]] = None,
    ) -> None:
        """Build and validate the homology graph.

        Args:
            edges: Iterable of ``(human_gene, mouse_gene, homology_type)``
                triples. Gene labels must be non-empty strings and
                ``homology_type`` must be a key of ``tier_weights``. Duplicate
                ``(human_gene, mouse_gene)`` pairs are collapsed, keeping the
                maximum edge weight (the most confident call wins; this makes
                the graph robust to redundant annotation releases).
            tier_weights: Optional custom mapping from homology type to a
                finite confidence weight in ``(0, 1]``. Defaults to
                ``TIER_WEIGHTS``.

        Raises:
            ValueError: If ``tier_weights`` is empty or contains an empty
                type label or a weight outside ``(0, 1]``; if ``edges`` is
                empty; if any edge is malformed (not a 3-tuple, empty gene
                label, or unknown homology type).
        """
        weights = dict(TIER_WEIGHTS if tier_weights is None else tier_weights)
        if not weights:
            raise ValueError("tier_weights must map at least one homology type to a weight.")
        for homology_type, weight in weights.items():
            if not isinstance(homology_type, str) or not homology_type.strip():
                raise ValueError("tier_weights keys must be non-empty strings.")
            if not np.isfinite(weight) or weight <= 0.0 or weight > 1.0:
                raise ValueError(
                    f"tier weight for {homology_type!r} must be finite and in (0, 1]; got {weight!r}."
                )
        self.tier_weights = weights

        edge_list = list(edges)
        if not edge_list:
            raise ValueError("edges must contain at least one homology edge.")

        # Collapse duplicate (human, mouse) pairs to their maximum weight so
        # downstream linear algebra sees one weight per gene pair.
        collapsed = {}
        human_genes = set()
        mouse_genes = set()
        for position, edge in enumerate(edge_list):
            try:
                human_gene, mouse_gene, homology_type = edge
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"edge #{position} must be a (human_gene, mouse_gene, homology_type) triple; "
                    f"got {edge!r}."
                ) from exc
            for label, role in ((human_gene, "human"), (mouse_gene, "mouse")):
                if not isinstance(label, str) or not label.strip():
                    raise ValueError(
                        f"edge #{position} has an invalid {role} gene label: {label!r}."
                    )
            if homology_type not in weights:
                raise ValueError(
                    f"edge #{position} uses unknown homology type {homology_type!r}; "
                    f"known types: {sorted(weights)}."
                )
            weight = float(weights[homology_type])
            key = (human_gene, mouse_gene)
            if key not in collapsed or weight > collapsed[key]:
                collapsed[key] = weight
            human_genes.add(human_gene)
            mouse_genes.add(mouse_gene)

        # Sorted orders make every downstream matrix deterministic.
        self._human_genes = tuple(sorted(human_genes))
        self._mouse_genes = tuple(sorted(mouse_genes))
        self._human_index = {g: i for i, g in enumerate(self._human_genes)}
        self._mouse_index = {g: i for i, g in enumerate(self._mouse_genes)}
        # Raw (unnormalized) weights, keyed by (mouse_index, human_index) so
        # both the projection matrix and the confidence matrix can be built.
        self._raw_weights = {
            (self._mouse_index[m], self._human_index[h]): w
            for (h, m), w in collapsed.items()
        }

    @property
    def human_genes(self) -> Tuple[str, ...]:
        """Sorted tuple of unique human gene labels."""
        return self._human_genes

    @property
    def mouse_genes(self) -> Tuple[str, ...]:
        """Sorted tuple of unique mouse gene labels."""
        return self._mouse_genes

    def projection_matrix(self) -> np.ndarray:
        """Row-stochastic mouse-to-human projection matrix.

        Returns:
            Array of shape ``(n_mouse_genes, n_human_genes)`` where row ``j``
            distributes one unit of mouse-gene mass across that gene's human
            homologs in proportion to the edge confidence weights (each row
            sums to exactly 1.0, since every mouse gene has at least one edge
            by construction).
        """
        matrix = np.zeros((len(self._mouse_genes), len(self._human_genes)))
        for (mouse_i, human_i), weight in self._raw_weights.items():
            matrix[mouse_i, human_i] = weight
        row_sums = matrix.sum(axis=1, keepdims=True)
        # Every row has at least one positive entry, so no division guard is
        # needed; an explicit assert documents the invariant.
        assert (row_sums > 0).all(), "each mouse gene must carry at least one edge"
        return matrix / row_sums

    def _weighted_tier_matrix(self) -> np.ndarray:
        """Unnormalized projection matrix multiplied by its own edge weights.

        Returns:
            Array of shape ``(n_mouse_genes, n_human_genes)`` with entries
            ``P[j, i] * w(j, i)`` where ``P`` is the row-stochastic
            projection and ``w`` the edge confidence weight. Used to compute
            contribution-weighted mean tier weights per aligned pair.
        """
        projection = self.projection_matrix()
        weighted = np.zeros_like(projection)
        for (mouse_i, human_i), weight in self._raw_weights.items():
            weighted[mouse_i, human_i] = projection[mouse_i, human_i] * weight
        return weighted


@dataclass
class AlignmentResult:
    """Output of :func:`compute_alignment`.

    Attributes:
        human_labels: Labels of the human programs (rows of ``scores``).
        mouse_labels: Labels of the mouse programs (columns of ``scores``).
        scores: Mutual-best-rescaled alignment scores of shape
            ``(n_human_programs, n_mouse_programs)`` with values in [0, 1];
            1.0 exactly when two programs are each other's best match.
        confidence: Contribution-weighted mean homology-edge weight per
            pair, same shape as ``scores``; values in (0, 1].
        confidence_tiers: Per-pair tier label (``high`` / ``medium`` /
            ``low``) derived from ``confidence`` via
            ``CONFIDENCE_TIER_THRESHOLDS``; same shape as ``scores``.
        similarity: Cosine-similarity matrix after refinement but *before*
            mutual-best rescaling; kept for diagnostics.
        n_refinements: Number of SAMap-style refinement iterations applied.
    """

    human_labels: Tuple[str, ...]
    mouse_labels: Tuple[str, ...]
    scores: np.ndarray
    confidence: np.ndarray
    confidence_tiers: np.ndarray
    similarity: np.ndarray
    n_refinements: int


def project_mouse_programs(mouse_scores: np.ndarray, graph: HomologyGraph) -> np.ndarray:
    """Project mouse program scores into human gene space.

    Args:
        mouse_scores: Array of shape ``(n_mouse_programs, n_mouse_genes)``
            with finite, non-negative program scores.
        graph: Homology graph whose mouse-gene axis defines the columns of
            ``mouse_scores``.

    Returns:
        Array of shape ``(n_mouse_programs, n_human_genes)``: each mouse
        program re-expressed in human gene space, with mass flowing through
        confidence-weighted homology edges.

    Raises:
        ValueError: If ``mouse_scores`` fails validation (see
            ``_validate_score_matrix``).
    """
    matrix = _validate_score_matrix(mouse_scores, len(graph.mouse_genes), "mouse_scores")
    return matrix @ graph.projection_matrix()


def compute_alignment(
    human_scores: np.ndarray,
    mouse_scores: np.ndarray,
    graph: HomologyGraph,
    human_labels: Optional[Sequence[str]] = None,
    mouse_labels: Optional[Sequence[str]] = None,
    n_refinements: int = 2,
) -> AlignmentResult:
    """Align human and mouse gene programs with confidence weighting.

    Args:
        human_scores: Array of shape ``(n_human_programs, n_human_genes)``
            with finite, non-negative program scores (e.g. marker
            specificity weights). Zero rows are allowed and align to nothing.
        mouse_scores: Array of shape ``(n_mouse_programs, n_mouse_genes)``,
            same constraints.
        graph: Homology graph relating the two gene axes.
        human_labels: Optional unique labels for the human programs;
            defaults to ``human_0, human_1, ...``.
        mouse_labels: Optional unique labels for the mouse programs;
            defaults to ``mouse_0, mouse_1, ...``.
        n_refinements: Number of SAMap-style gene-reweighting iterations
            (``0`` = plain projection + cosine similarity). Each iteration
            up-weights human genes that drive the currently aligned pairs,
            so small values (1-3) sharpen the alignment; large values are
            unnecessary because the reweighting converges quickly.

    Returns:
        An :class:`AlignmentResult` with deterministic scores, per-pair
        confidence values, and confidence tiers.

    Raises:
        ValueError: If either score matrix is malformed (wrong width,
            negative or non-finite entries), if labels are missing,
            duplicated, or of the wrong length, or if ``n_refinements`` is
            negative or not an integer.
    """
    if isinstance(n_refinements, bool) or not isinstance(n_refinements, int):
        raise ValueError(f"n_refinements must be an integer; got {n_refinements!r}.")
    if n_refinements < 0:
        raise ValueError(f"n_refinements must be >= 0; got {n_refinements}.")

    human = _validate_score_matrix(human_scores, len(graph.human_genes), "human_scores")
    mouse = _validate_score_matrix(mouse_scores, len(graph.mouse_genes), "mouse_scores")
    human_labels = _validate_labels(human_labels, human.shape[0], "human")
    mouse_labels = _validate_labels(mouse_labels, mouse.shape[0], "mouse")

    projected = project_mouse_programs(mouse, graph)

    # Step 2: initial cross-species similarity in the shared gene space.
    similarity = _cosine_similarity(human, projected)

    # Step 3: SAMap-style refinement -- the current program alignment
    # re-weights genes, and the re-weighted genes sharpen the alignment.
    refined_human = human
    refined_projected = projected
    for _ in range(n_refinements):
        gene_weights = _driving_gene_weights(refined_human, refined_projected, similarity)
        if gene_weights is None:
            # No positive alignment signal exists yet (e.g. all-zero inputs);
            # further refinement cannot change anything, so stop early.
            break
        refined_human = refined_human * gene_weights
        refined_projected = refined_projected * gene_weights
        similarity = _cosine_similarity(refined_human, refined_projected)

    # Step 4: mutual-best rescaling into [0, 1]. Dividing by the geometric
    # mean of the row and column maxima yields exactly 1.0 for mutual-best
    # pairs and dampens one-sided matches.
    row_max = similarity.max(axis=1, keepdims=True)
    col_max = similarity.max(axis=0, keepdims=True)
    scale = np.sqrt(row_max * col_max)
    scores = np.divide(
        similarity,
        scale,
        out=np.zeros_like(similarity),
        where=scale > 0,
    )

    confidence = _pair_confidence(human, mouse, projected, graph)
    tiers = _assign_tiers(confidence)

    return AlignmentResult(
        human_labels=tuple(human_labels),
        mouse_labels=tuple(mouse_labels),
        scores=scores,
        confidence=confidence,
        confidence_tiers=tiers,
        similarity=similarity,
        n_refinements=n_refinements,
    )


def _validate_score_matrix(matrix: np.ndarray, n_genes: int, name: str) -> np.ndarray:
    """Check a program-by-gene score matrix and return it as float64.

    Args:
        matrix: Candidate 2-D array of program scores.
        n_genes: Required number of columns (gene axis of the graph).
        name: Parameter name used in error messages.

    Returns:
        The matrix as a ``float64`` NumPy array.

    Raises:
        ValueError: If the matrix is not 2-D, has the wrong number of
            columns, or contains negative or non-finite entries. Scores are
            required to be non-negative because they represent marker /
            specificity weights; signed statistics (e.g. log fold-changes)
            must be transformed by the caller first.
    """
    array = np.asarray(matrix, dtype=float)
    if array.ndim != 2:
        raise ValueError(f"{name} must be a 2-D array; got {array.ndim} dimensions.")
    if array.shape[1] != n_genes:
        raise ValueError(
            f"{name} has {array.shape[1]} columns but the homology graph defines "
            f"{n_genes} genes on that axis."
        )
    if not np.isfinite(array).all():
        raise ValueError(f"{name} must contain only finite values (no NaN or infinities).")
    if (array < 0).any():
        raise ValueError(f"{name} must be non-negative; transform signed statistics first.")
    return array


def _validate_labels(
    labels: Optional[Sequence[str]], n_programs: int, species: str
) -> Tuple[str, ...]:
    """Validate or default the program labels for one species.

    Args:
        labels: Caller-supplied labels, or ``None`` for defaults.
        n_programs: Number of programs (rows of the score matrix).
        species: ``"human"`` or ``"mouse"``; used for defaults and errors.

    Returns:
        Tuple of ``n_programs`` unique string labels.

    Raises:
        ValueError: If labels are provided but have the wrong length, contain
            non-string or empty entries, or contain duplicates (ambiguous
            labels would make result rows/columns unaddressable).
    """
    if labels is None:
        return tuple(f"{species}_{i}" for i in range(n_programs))
    result = tuple(labels)
    if len(result) != n_programs:
        raise ValueError(
            f"{species}_labels has {len(result)} entries but {n_programs} programs were given."
        )
    for label in result:
        if not isinstance(label, str) or not label.strip():
            raise ValueError(f"{species}_labels entries must be non-empty strings; got {label!r}.")
    if len(set(result)) != len(result):
        raise ValueError(f"{species}_labels must be unique; got duplicates in {result!r}.")
    return result


def _cosine_similarity(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Zero-safe pairwise cosine similarity between rows of ``a`` and ``b``.

    Args:
        a: Array of shape ``(n, g)`` with non-negative entries.
        b: Array of shape ``(m, g)`` with non-negative entries.

    Returns:
        Array of shape ``(n, m)`` with values in [0, 1]; rows or columns
        involving a zero vector are exactly 0 (never NaN).
    """
    a_norm = np.linalg.norm(a, axis=1, keepdims=True)
    b_norm = np.linalg.norm(b, axis=1, keepdims=True)
    denominator = a_norm * b_norm.T
    return np.divide(
        a @ b.T,
        denominator,
        out=np.zeros((a.shape[0], b.shape[0])),
        where=denominator > 0,
    )


def _driving_gene_weights(
    human: np.ndarray, projected: np.ndarray, similarity: np.ndarray
) -> Optional[np.ndarray]:
    """Compute SAMap-style gene weights from the current alignment.

    Genes receive weight in proportion to how much they contribute to the
    dot products of currently aligned program pairs; genes used mainly by
    unaligned pairs are down-weighted. This mirrors the SAMap alternation in
    which the cell-graph alignment refines the gene graph.

    Args:
        human: Human program matrix, shape ``(n_h, g)``.
        projected: Projected mouse program matrix, shape ``(n_m, g)``.
        similarity: Current alignment matrix, shape ``(n_h, n_m)``.

    Returns:
        Array of shape ``(g,)`` normalized to mean 1.0 (so refinement does
        not change the overall scale), or ``None`` when no positive
        alignment signal exists (in which case refinement is a no-op and the
        caller stops iterating).
    """
    # Mass each human program drives onto each gene through its aligned
    # partners: (n_h, n_m) @ (n_m, g) -> (n_h, g).
    driven = similarity @ projected
    weights = (human * driven).sum(axis=0)
    total = weights.sum()
    if total <= 0:
        return None
    # Normalize to mean 1 so repeated iterations stay scale-stable.
    return weights * (weights.size / total)


def _pair_confidence(
    human: np.ndarray,
    mouse: np.ndarray,
    projected: np.ndarray,
    graph: HomologyGraph,
) -> np.ndarray:
    """Contribution-weighted mean homology-edge weight per program pair.

    For a pair ``(h, m)``, each homology edge contributes
    ``H[h, i] * M[m, j] * P[j, i]`` units of alignment mass flowing from
    mouse gene ``j`` into human gene ``i``; the confidence is the weighted
    average of the edge confidence weights under that flow. A pair whose
    support flows entirely through one-to-one edges scores 1.0.

    Args:
        human: Validated human program matrix, shape ``(n_h, g_h)``.
        mouse: Validated mouse program matrix, shape ``(n_m, g_m)``.
        projected: Mouse programs in human gene space, shape ``(n_m, g_h)``.
        graph: The homology graph used for the projection.

    Returns:
        Array of shape ``(n_h, n_m)`` with values in ``(0, 1]``; pairs with
        no shared mass (zero denominator) receive exactly 0.0.
    """
    confidence_space = mouse @ graph._weighted_tier_matrix()
    numerator = np.einsum("hi,mi->hm", human, confidence_space)
    denominator = np.einsum("hi,mi->hm", human, projected)
    return np.divide(
        numerator,
        denominator,
        out=np.zeros_like(numerator),
        where=denominator > 0,
    )


def _assign_tiers(confidence: np.ndarray) -> np.ndarray:
    """Map confidence values onto the documented tier labels.

    Args:
        confidence: Array of per-pair confidence values in [0, 1].

    Returns:
        Array of strings (``"high"`` / ``"medium"`` / ``"low"``) with the
        same shape, using ``CONFIDENCE_TIER_THRESHOLDS``.
    """
    high = CONFIDENCE_TIER_THRESHOLDS["high"]
    medium = CONFIDENCE_TIER_THRESHOLDS["medium"]
    tiers = np.full(confidence.shape, _LOW_TIER_LABEL, dtype=object)
    tiers[confidence >= medium] = "medium"
    tiers[confidence >= high] = "high"
    return tiers
