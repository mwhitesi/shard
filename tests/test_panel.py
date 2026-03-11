"""Tests for panel definition and spaced seed design."""

import json
import tempfile
from pathlib import Path

import numpy as np
import pytest

from shard.data.panel import (
    PanelDefinition,
    SpacedSeed,
    SpacedSeedSet,
    TargetRegion,
    _parse_ales_output,
    _parse_iedera_output,
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
            n_anchor_kmers=50,
        )
        assert t.length == 200
        assert t.anchor_density == pytest.approx(0.25)


class TestPanelDefinition:
    def _make_panel(self) -> PanelDefinition:
        return PanelDefinition(
            panel_id="test_panel",
            kmer_size=31,
            targets=[
                TargetRegion("t0", "chr1", 100, 300, 0.45, 50),
                TargetRegion("t1", "chr2", 500, 800, 0.55, 80),
            ],
            seed_set=SpacedSeedSet.exact_match(31),
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
        assert loaded.seed_set.seeds[0].pattern == "1" * 31
        assert loaded.targets[0].gc_content == pytest.approx(0.45)
        assert loaded.targets[1].n_anchor_kmers == 80
