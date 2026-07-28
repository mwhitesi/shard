"""Tests for panel definition, spaced seeds, and fragment size estimation."""

import json
from pathlib import Path

import numpy as np
import pytest

from shard.data.panel import (
    AnchorHit,
    AnchorIndex,
    AnchorKmer,
    PanelDefinition,
    SpacedSeed,
    SpacedSeedSet,
    TargetRegion,
    _parse_ales_output,
    _parse_iedera_output,
    build_fragment_distribution,
    estimate_fragment_sizes,
)


class TestSpacedSeed:
    def test_basic_properties(self):
        seed = SpacedSeed(pattern="11101110111")
        assert seed.weight == 9
        assert seed.span == 11
        assert seed.care_positions == [0, 1, 2, 4, 5, 6, 8, 9, 10]
        assert seed.wildcard_positions == [3, 7]

    def test_exact_seed(self):
        seed = SpacedSeed(pattern="11111")
        assert seed.weight == 5
        assert seed.span == 5
        assert seed.wildcard_positions == []

    def test_invalid_pattern(self):
        with pytest.raises(ValueError, match="only '0' and '1'"):
            SpacedSeed(pattern="111-000")

    def test_extract_key(self):
        seed = SpacedSeed(pattern="11011")
        assert seed.extract_key("ACGTG") == "ACTG"

    def test_extract_key_wrong_length(self):
        seed = SpacedSeed(pattern="11011")
        with pytest.raises(ValueError, match="K-mer length"):
            seed.extract_key("ACG")


class TestSpacedSeedSet:
    def test_exact_match_factory(self):
        ss = SpacedSeedSet.exact_match(31)
        assert ss.n_seeds == 1
        assert ss.weight == 31
        assert ss.kmer_size == 31
        assert ss.fraction_positions_covered == 0.0  # no wildcards

    def test_wildcard_coverage(self):
        seeds = SpacedSeedSet(
            seeds=[
                SpacedSeed(pattern="11011"),
                SpacedSeed(pattern="10111"),
            ],
            kmer_size=5,
        )
        coverage = seeds.wildcard_coverage()
        # Seed 1 wildcards: pos 2. Seed 2 wildcards: pos 1.
        np.testing.assert_array_equal(coverage, [False, True, True, False, False])
        assert seeds.fraction_positions_covered == pytest.approx(0.4)

    def test_span_mismatch_rejected(self):
        with pytest.raises(ValueError, match="kmer_size"):
            SpacedSeedSet(
                seeds=[SpacedSeed(pattern="11011")],
                kmer_size=10,
            )

    def test_serialisation_roundtrip(self):
        original = SpacedSeedSet(
            seeds=[
                SpacedSeed(pattern="11011"),
                SpacedSeed(pattern="10111"),
            ],
            kmer_size=5,
        )
        d = original.to_dict()
        restored = SpacedSeedSet.from_dict(d)
        assert restored.kmer_size == original.kmer_size
        assert [s.pattern for s in restored.seeds] == [s.pattern for s in original.seeds]


class TestParseOutput:
    def test_parse_ales_output(self):
        stdout = (
            "Some header line\n"
            "1110111011101110111011101110111\n"
            "1011101110111011101110111011101\n"
            "Sensitivity: 0.9876\n"
        )
        seeds = _parse_ales_output(stdout, expected_span=31)
        assert len(seeds) == 2
        assert seeds[0].pattern == "1110111011101110111011101110111"
        assert seeds[0].weight == 24
        assert seeds[0].span == 31

    def test_parse_iedera_output(self):
        stdout = (
            "###-#--#-#--##-###\t0.999\t0.467\t0.533\n"
            "#-###--##-#--#-###\t0.999\t0.450\t0.550\n"
        )
        seeds = _parse_iedera_output(stdout)
        assert len(seeds) == 2
        assert seeds[0].pattern == "111010010100110111"
        assert seeds[0].weight == 11

    def test_parse_ales_ignores_wrong_length(self):
        stdout = "11011\n1110111011101110111011101110111\n"
        seeds = _parse_ales_output(stdout, expected_span=31)
        assert len(seeds) == 1


class TestTargetRegion:
    def test_basic(self):
        t = TargetRegion(
            name="BRCA1_exon10",
            chrom="chr17",
            start=100,
            end=300,
            gc_content=0.45,
            anchor_offsets=list(range(0, 200, 4)),
        )
        assert t.length == 200
        assert t.n_anchor_kmers == 50
        assert t.anchor_density == pytest.approx(0.25)

    def test_default_empty_anchors(self):
        t = TargetRegion(name="t0", chrom="chr1", start=0, end=100)
        assert t.n_anchor_kmers == 0
        assert t.anchor_density == 0.0


class TestPanelDefinition:
    def _make_panel(self) -> PanelDefinition:
        return PanelDefinition(
            panel_id="test_panel",
            kmer_size=5,
            targets=[
                TargetRegion("t0", "chr1", 100, 300, 0.45, [0, 10, 20, 50, 100, 150]),
                TargetRegion("t1", "chr2", 500, 800, 0.55, [0, 30, 60, 90, 120, 200, 250]),
            ],
            seed_set=SpacedSeedSet.exact_match(5),
        )

    def test_properties(self):
        panel = self._make_panel()
        assert panel.n_targets == 2
        assert panel.target_names == ["t0", "t1"]

    def test_save_load_roundtrip(self, tmp_path: Path):
        panel = self._make_panel()
        path = tmp_path / "panel.json"
        panel.save(path)

        loaded = PanelDefinition.load(path)
        assert loaded.panel_id == panel.panel_id
        assert loaded.kmer_size == panel.kmer_size
        assert loaded.n_targets == panel.n_targets
        assert loaded.target_names == panel.target_names
        assert loaded.seed_set.n_seeds == 1
        assert loaded.seed_set.seeds[0].pattern == "1" * 5
        assert loaded.targets[0].gc_content == pytest.approx(0.45)
        assert loaded.targets[0].anchor_offsets == [0, 10, 20, 50, 100, 150]
        assert loaded.targets[1].n_anchor_kmers == 7


# ---------------------------------------------------------------------------
# Anchor index tests
# ---------------------------------------------------------------------------


def _make_test_panel_and_anchors():
    """Create a small panel with k=5 and known anchors for testing."""
    panel = PanelDefinition(
        panel_id="test",
        kmer_size=5,
        targets=[
            TargetRegion("t0", "chr1", 1000, 1200, 0.5, [0, 50, 100, 150]),
            TargetRegion("t1", "chr2", 2000, 2300, 0.5, [0, 80, 160, 240]),
        ],
        seed_set=SpacedSeedSet.exact_match(5),
    )
    anchors = [
        # Target 0 anchors at offsets 0, 50, 100, 150
        AnchorKmer("ACGTG", target_idx=0, offset=0),
        AnchorKmer("CGTGA", target_idx=0, offset=50),
        AnchorKmer("GTGAC", target_idx=0, offset=100),
        AnchorKmer("TGACC", target_idx=0, offset=150),
        # Target 1 anchors at offsets 0, 80, 160, 240
        AnchorKmer("TTAAG", target_idx=1, offset=0),
        AnchorKmer("AAGCC", target_idx=1, offset=80),
        AnchorKmer("GCCTT", target_idx=1, offset=160),
        AnchorKmer("CTTAA", target_idx=1, offset=240),
    ]
    return panel, anchors


class TestAnchorIndex:
    def test_exact_lookup(self):
        panel, anchors = _make_test_panel_and_anchors()
        index = AnchorIndex(panel, anchors)

        hits = index.lookup("ACGTG")
        assert hits == [(0, 0)]

        hits = index.lookup("AAGCC")
        assert hits == [(1, 80)]

    def test_no_match(self):
        panel, anchors = _make_test_panel_and_anchors()
        index = AnchorIndex(panel, anchors)
        assert index.lookup("NNNNN") == []

    def test_spaced_seed_tolerates_snp(self):
        """A SNP at a wildcard position should still match."""
        panel, anchors = _make_test_panel_and_anchors()
        # Replace exact seed with a spaced seed: care positions 0,1,3,4 (skip pos 2)
        panel.seed_set = SpacedSeedSet(
            seeds=[SpacedSeed(pattern="11011")],
            kmer_size=5,
        )
        index = AnchorIndex(panel, anchors)

        # Original anchor: ACGTG → spaced key = ACTG (positions 0,1,3,4)
        # Query with SNP at position 2: ACATG → spaced key = ACTG (same!)
        hits = index.lookup("ACATG")
        assert (0, 0) in hits

        # SNP at care position should NOT match
        # XCGTG → spaced key = XCTG ≠ ACTG
        hits = index.lookup("TCGTG")
        assert (0, 0) not in hits

    def test_n_entries(self):
        panel, anchors = _make_test_panel_and_anchors()
        index = AnchorIndex(panel, anchors)
        assert index.n_entries == 8  # 8 anchors × 1 seed


# ---------------------------------------------------------------------------
# Fragment size estimation tests
# ---------------------------------------------------------------------------


class TestEstimateFragmentSizes:
    def test_same_target_pair(self):
        """R1 and R2 anchors in the same target give a fragment estimate."""
        panel, _ = _make_test_panel_and_anchors()
        r1_hits = [AnchorHit(target_idx=0, offset=0, read_position=10)]
        r2_hits = [AnchorHit(target_idx=0, offset=150, read_position=5)]

        estimates = estimate_fragment_sizes(r1_hits, r2_hits, panel, read_length=150)
        assert len(estimates) == 1
        target_idx, frag_size = estimates[0]
        assert target_idx == 0
        assert frag_size == 150 + 150  # |0 - 150| + read_length = 300

    def test_different_targets_no_estimate(self):
        """R1 in target 0, R2 in target 1 → no fragment estimate."""
        panel, _ = _make_test_panel_and_anchors()
        r1_hits = [AnchorHit(target_idx=0, offset=50, read_position=10)]
        r2_hits = [AnchorHit(target_idx=1, offset=80, read_position=5)]

        estimates = estimate_fragment_sizes(r1_hits, r2_hits, panel)
        assert estimates == []

    def test_multiple_anchors_multiple_estimates(self):
        """Multiple anchor hits per read produce multiple estimates."""
        panel, _ = _make_test_panel_and_anchors()
        r1_hits = [
            AnchorHit(target_idx=0, offset=0, read_position=5),
            AnchorHit(target_idx=0, offset=50, read_position=55),
        ]
        r2_hits = [
            AnchorHit(target_idx=0, offset=150, read_position=10),
        ]

        estimates = estimate_fragment_sizes(r1_hits, r2_hits, panel, read_length=150)
        assert len(estimates) == 2
        sizes = sorted(frag for _, frag in estimates)
        assert sizes == [250, 300]  # |50-150|+150=250, |0-150|+150=300

    def test_typical_fragment_size(self):
        """A pair 200bp apart with 150bp reads gives ~350bp fragment."""
        panel, _ = _make_test_panel_and_anchors()
        r1_hits = [AnchorHit(target_idx=0, offset=0, read_position=0)]
        r2_hits = [AnchorHit(target_idx=0, offset=200, read_position=0)]

        estimates = estimate_fragment_sizes(r1_hits, r2_hits, panel, read_length=150)
        assert estimates[0] == (0, 350)


class TestBuildFragmentDistribution:
    def test_basic_distribution(self):
        estimates = [
            (0, 200),
            (0, 250),
            (0, 200),
            (1, 300),
        ]
        per_target, global_hist = build_fragment_distribution(
            estimates, n_targets=2, frag_range=(50, 500)
        )
        assert per_target.shape == (2, 450)
        assert global_hist.shape == (450,)

        # Target 0: two counts at bin 150 (200-50), one at bin 200 (250-50)
        assert per_target[0, 150] == 2.0
        assert per_target[0, 200] == 1.0
        # Target 1: one count at bin 250 (300-50)
        assert per_target[1, 250] == 1.0
        # Global: sum
        assert global_hist.sum() == 4.0

    def test_out_of_range_discarded(self):
        estimates = [
            (0, 30),   # below frag_range
            (0, 200),  # in range
            (0, 600),  # above frag_range
        ]
        per_target, global_hist = build_fragment_distribution(
            estimates, n_targets=1, frag_range=(50, 500)
        )
        assert global_hist.sum() == 1.0  # only the 200bp estimate counted

    def test_empty_estimates(self):
        per_target, global_hist = build_fragment_distribution(
            [], n_targets=3, frag_range=(50, 500)
        )
        assert per_target.sum() == 0.0
        assert global_hist.sum() == 0.0
