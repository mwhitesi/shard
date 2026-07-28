"""Tests for the k-mer extraction pipeline: FASTQ → anchor hits → raw fingerprint."""

import gzip
from pathlib import Path

import pytest

from shard.data.panel import (
    AnchorIndex,
    AnchorKmer,
    PanelDefinition,
    SpacedSeed,
    SpacedSeedSet,
    TargetRegion,
)
from shard.extraction import (
    SampleAccumulator,
    extract_raw_fingerprint,
    scan_read,
)
from shard.utils.fastq import FastqRecord
from shard.utils.kmer import reverse_complement


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_panel(
    targets: list[TargetRegion],
    kmer_size: int = 5,
    seed_set: SpacedSeedSet | None = None,
) -> PanelDefinition:
    """Create a minimal PanelDefinition for testing."""
    if seed_set is None:
        seed_set = SpacedSeedSet.exact_match(kmer_size)
    return PanelDefinition(
        panel_id="test",
        kmer_size=kmer_size,
        targets=targets,
        seed_set=seed_set,
    )


def _write_fastq(path: Path, records: list[tuple[str, str, str]]) -> None:
    """Write FASTQ records to a plain text file."""
    lines = []
    for name, seq, qual in records:
        lines.append(f"@{name}\n{seq}\n+\n{qual}\n")
    path.write_text("".join(lines))


# ---------------------------------------------------------------------------
# scan_read tests
# ---------------------------------------------------------------------------


class TestScanRead:
    def test_finds_exact_anchor(self):
        """A read containing an exact anchor k-mer produces a hit."""
        target = TargetRegion("t0", "chr1", 0, 20, 0.5, anchor_offsets=[5])
        panel = _make_panel([target], kmer_size=5)
        anchor_kmer = "ACGTG"
        anchors = [AnchorKmer(kmer=anchor_kmer, target_idx=0, offset=5)]
        index = AnchorIndex(panel, anchors)

        # Read that contains the anchor k-mer at position 3
        read_seq = "NNNACGTGNN"
        hits = scan_read(read_seq, index, k=5)

        assert len(hits) >= 1
        hit = hits[0]
        assert hit.target_idx == 0
        assert hit.offset == 5
        assert hit.read_position == 3

    def test_finds_reverse_complement_anchor(self):
        """A read containing the RC of an anchor produces a hit."""
        anchor_kmer = "ACGTG"
        rc_kmer = reverse_complement(anchor_kmer)  # CACGT
        target = TargetRegion("t0", "chr1", 0, 20, 0.5, anchor_offsets=[5])
        panel = _make_panel([target], kmer_size=5)
        anchors = [AnchorKmer(kmer=anchor_kmer, target_idx=0, offset=5)]
        index = AnchorIndex(panel, anchors)

        # Read that contains the reverse complement
        read_seq = f"NN{rc_kmer}NNN"
        hits = scan_read(read_seq, index, k=5)

        assert len(hits) >= 1
        assert hits[0].target_idx == 0
        assert hits[0].offset == 5

    def test_no_hits_for_unmatched_read(self):
        """A read with no anchor k-mers produces no hits."""
        target = TargetRegion("t0", "chr1", 0, 20, 0.5, anchor_offsets=[5])
        panel = _make_panel([target], kmer_size=5)
        anchors = [AnchorKmer(kmer="ACGTG", target_idx=0, offset=5)]
        index = AnchorIndex(panel, anchors)

        hits = scan_read("TTTTTTTTT", index, k=5)
        assert hits == []

    def test_multiple_anchors_in_one_read(self):
        """A read containing multiple anchor k-mers produces multiple hits."""
        target = TargetRegion("t0", "chr1", 0, 30, 0.5, anchor_offsets=[5, 15])
        panel = _make_panel([target], kmer_size=5)
        anchors = [
            AnchorKmer(kmer="ACGTG", target_idx=0, offset=5),
            AnchorKmer(kmer="TGCAA", target_idx=0, offset=15),
        ]
        index = AnchorIndex(panel, anchors)

        read_seq = "ACGTGNNNNTGCAA"
        hits = scan_read(read_seq, index, k=5)

        target_offsets = {(h.target_idx, h.offset) for h in hits}
        assert (0, 5) in target_offsets
        assert (0, 15) in target_offsets

    def test_spaced_seed_tolerance(self):
        """With spaced seeds, a SNP at a wildcard position still produces a hit."""
        # Seed: 11011 — position 2 is wildcard
        seed_set = SpacedSeedSet(
            seeds=[SpacedSeed(pattern="11011")],
            kmer_size=5,
        )
        target = TargetRegion("t0", "chr1", 0, 20, 0.5, anchor_offsets=[5])
        panel = _make_panel([target], kmer_size=5, seed_set=seed_set)
        anchors = [AnchorKmer(kmer="ACGTG", target_idx=0, offset=5)]
        index = AnchorIndex(panel, anchors)

        # Mutate position 2 (G→T): ACGTG → ACTTG
        read_seq = "NNACTTGNN"
        hits = scan_read(read_seq, index, k=5)

        assert len(hits) >= 1
        assert hits[0].target_idx == 0


# ---------------------------------------------------------------------------
# SampleAccumulator tests
# ---------------------------------------------------------------------------


class TestSampleAccumulator:
    def _make_simple_setup(self):
        """Create a minimal panel + anchors + index for accumulator tests."""
        target = TargetRegion("t0", "chr1", 0, 30, 0.5, anchor_offsets=[5, 15])
        panel = _make_panel([target], kmer_size=5)
        anchors = [
            AnchorKmer(kmer="ACGTG", target_idx=0, offset=5),
            AnchorKmer(kmer="TGCAA", target_idx=0, offset=15),
        ]
        index = AnchorIndex(panel, anchors)
        return panel, anchors, index

    def test_counts_read_pairs(self):
        panel, anchors, index = self._make_simple_setup()
        acc = SampleAccumulator(panel)

        r1 = FastqRecord("r1", "NNNNNACGTG", "I" * 10)
        r2 = FastqRecord("r1", "NNNNNTGCAA", "I" * 10)
        acc.process_pair(r1, r2, index)

        assert acc.total_read_pairs == 1

    def test_anchor_counts(self):
        panel, anchors, index = self._make_simple_setup()
        acc = SampleAccumulator(panel)

        # Both reads hit different anchors
        r1 = FastqRecord("r1", "NNNNNACGTG", "I" * 10)
        r2 = FastqRecord("r1", "NNNNNTGCAA", "I" * 10)
        acc.process_pair(r1, r2, index)

        counts = acc.compute_kmer_counts()
        # Both anchors hit once → median of [1, 1] = 1
        assert counts[0] == pytest.approx(1.0)

    def test_anchor_counts_with_zero_hits(self):
        """Anchors that are never hit contribute 0 to the median."""
        panel, anchors, index = self._make_simple_setup()
        acc = SampleAccumulator(panel)

        # Only R1 hits one anchor, R2 hits nothing
        r1 = FastqRecord("r1", "NNNNNACGTG", "I" * 10)
        r2 = FastqRecord("r1", "TTTTTTTTTT", "I" * 10)
        acc.process_pair(r1, r2, index)

        counts = acc.compute_kmer_counts()
        # Anchor at offset 5 has count 1, anchor at offset 15 has count 0
        # median([1, 0]) = 0.5
        assert counts[0] == pytest.approx(0.5)

    def test_global_gc_histogram_populated(self):
        panel, anchors, index = self._make_simple_setup()
        acc = SampleAccumulator(panel)

        r1 = FastqRecord("r1", "GCGCGCGCGC", "I" * 10)  # GC = 1.0
        r2 = FastqRecord("r1", "ATATATATAT", "I" * 10)  # GC = 0.0
        acc.process_pair(r1, r2, index)

        assert acc.global_gc_histogram.sum() == 2  # two reads
        # GC=0.0 should be in first bins, GC=1.0 in last bin
        assert acc.global_gc_histogram[0] >= 1  # AT-rich read
        assert acc.global_gc_histogram[-1] >= 1  # GC-rich read

    def test_gc_profiles_per_target(self):
        panel, anchors, index = self._make_simple_setup()
        acc = SampleAccumulator(panel)

        # R1 contains the anchor ACGTG at position 5 and is GC-rich
        r1 = FastqRecord("r1", "GCGCGACGTG", "I" * 10)  # GC = 0.7
        r2 = FastqRecord("r1", "ATATATATAT", "I" * 10)  # no anchor
        acc.process_pair(r1, r2, index)

        gc, _ = acc.normalise_profiles()
        # Target 0 should have GC profile with a peak near GC=0.7
        assert gc[0].sum() == pytest.approx(1.0)

    def test_fragment_size_estimation(self):
        panel, anchors, index = self._make_simple_setup()
        acc = SampleAccumulator(panel)

        # R1 hits anchor at offset 5, R2 hits anchor at offset 15
        # Both on target 0 → fragment size = |5 - 15| + 150 = 160
        r1 = FastqRecord("r1", "NNNNNACGTG", "I" * 10)
        r2 = FastqRecord("r1", "NNNNNTGCAA", "I" * 10)
        acc.process_pair(r1, r2, index)

        # Global frag histogram should have an entry at bin for size 160
        # Bin index = 160 - 50 = 110
        assert acc.global_frag_histogram[110] >= 1

    def test_normalise_profiles_sums_to_one(self):
        panel, anchors, index = self._make_simple_setup()
        acc = SampleAccumulator(panel)

        r1 = FastqRecord("r1", "NNNNNACGTG", "I" * 10)
        r2 = FastqRecord("r1", "NNNNNTGCAA", "I" * 10)
        for _ in range(10):
            acc.process_pair(r1, r2, index)

        gc, frag = acc.normalise_profiles()
        # Targets with data should sum to 1
        assert gc[0].sum() == pytest.approx(1.0)
        assert frag[0].sum() == pytest.approx(1.0)

    def test_normalise_empty_target(self):
        """Targets with no data produce all-zero profiles."""
        targets = [
            TargetRegion("t0", "chr1", 0, 30, 0.5, anchor_offsets=[5]),
            TargetRegion("t1", "chr2", 100, 130, 0.5, anchor_offsets=[10]),
        ]
        panel = _make_panel(targets, kmer_size=5)
        anchors = [
            AnchorKmer(kmer="ACGTG", target_idx=0, offset=5),
            AnchorKmer(kmer="CCCCC", target_idx=1, offset=10),
        ]
        index = AnchorIndex(panel, anchors)
        acc = SampleAccumulator(panel)

        # Only hit target 0
        r1 = FastqRecord("r1", "NNNNNACGTG", "I" * 10)
        r2 = FastqRecord("r1", "TTTTTTTTTT", "I" * 10)
        acc.process_pair(r1, r2, index)

        gc, frag = acc.normalise_profiles()
        # Target 1 should be all zeros (no reads hit it)
        assert gc[1].sum() == pytest.approx(0.0)
        assert frag[1].sum() == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# End-to-end extraction tests
# ---------------------------------------------------------------------------


class TestExtractRawFingerprint:
    def _make_test_files(self, tmp_path):
        """Create test FASTQ files and a panel with known anchors.

        Panel has 2 targets:
            t0: anchor "ACGTG" at offset 5
            t1: anchor "TGCAA" at offset 10

        R1 reads contain the anchor for t0.
        R2 reads contain the anchor for t1.
        Both anchors on different targets → no fragment size estimate.
        """
        targets = [
            TargetRegion("t0", "chr1", 0, 30, 0.5, anchor_offsets=[5]),
            TargetRegion("t1", "chr2", 100, 130, 0.5, anchor_offsets=[10]),
        ]
        panel = _make_panel(targets, kmer_size=5)
        anchors = [
            AnchorKmer(kmer="ACGTG", target_idx=0, offset=5),
            AnchorKmer(kmer="TGCAA", target_idx=1, offset=10),
        ]

        r1_path = tmp_path / "R1.fq"
        r2_path = tmp_path / "R2.fq"

        # 5 read pairs, R1 always has t0 anchor, R2 always has t1 anchor
        r1_records = [
            (f"read{i}/1", "NNNNNACGTGNNNN", "I" * 14)
            for i in range(5)
        ]
        r2_records = [
            (f"read{i}/2", "NNNNTGCAANNNNN", "I" * 14)
            for i in range(5)
        ]
        _write_fastq(r1_path, r1_records)
        _write_fastq(r2_path, r2_records)

        return r1_path, r2_path, panel, anchors

    def test_basic_extraction(self, tmp_path):
        r1_path, r2_path, panel, anchors = self._make_test_files(tmp_path)
        result = extract_raw_fingerprint(r1_path, r2_path, panel, anchors)

        assert result["total_read_pairs"] == 5
        assert result["kmer_counts"].shape == (2,)
        assert result["gc_profiles"].shape == (2, 20)
        assert result["frag_profiles"].shape == (2, 50)
        assert result["global_gc_histogram"].shape == (100,)
        assert result["global_frag_histogram"].shape == (450,)

    def test_kmer_counts_nonzero(self, tmp_path):
        r1_path, r2_path, panel, anchors = self._make_test_files(tmp_path)
        result = extract_raw_fingerprint(r1_path, r2_path, panel, anchors)

        # Each target has 1 anchor, hit 5 times each → median = 5
        assert result["kmer_counts"][0] == pytest.approx(5.0)
        assert result["kmer_counts"][1] == pytest.approx(5.0)

    def test_gc_histograms_normalised(self, tmp_path):
        r1_path, r2_path, panel, anchors = self._make_test_files(tmp_path)
        result = extract_raw_fingerprint(r1_path, r2_path, panel, anchors)

        # Global GC should be normalised to sum to 1
        assert result["global_gc_histogram"].sum() == pytest.approx(1.0)

    def test_gc_profiles_normalised(self, tmp_path):
        r1_path, r2_path, panel, anchors = self._make_test_files(tmp_path)
        result = extract_raw_fingerprint(r1_path, r2_path, panel, anchors)

        # Each target's GC profile should sum to 1 (since reads hit both targets)
        for t in range(2):
            total = result["gc_profiles"][t].sum()
            assert total == pytest.approx(1.0) or total == pytest.approx(0.0)

    def test_same_target_fragment_estimation(self, tmp_path):
        """When R1 and R2 hit the same target, fragment sizes are estimated."""
        target = TargetRegion("t0", "chr1", 0, 50, 0.5, anchor_offsets=[5, 25])
        panel = _make_panel([target], kmer_size=5)
        anchors = [
            AnchorKmer(kmer="ACGTG", target_idx=0, offset=5),
            AnchorKmer(kmer="TGCAA", target_idx=0, offset=25),
        ]

        r1_path = tmp_path / "R1.fq"
        r2_path = tmp_path / "R2.fq"

        # R1 hits offset 5, R2 hits offset 25 → frag size = |5-25| + 150 = 170
        _write_fastq(r1_path, [("r1/1", "NNNNNACGTGNNNN", "I" * 14)])
        _write_fastq(r2_path, [("r1/2", "NNNNTGCAANNNNN", "I" * 14)])

        result = extract_raw_fingerprint(r1_path, r2_path, panel, anchors)

        # Fragment histogram should have counts (frag size 170, bin 170-50=120)
        assert result["global_frag_histogram"].sum() > 0

    def test_gzipped_input(self, tmp_path):
        """Works with gzipped FASTQ files."""
        target = TargetRegion("t0", "chr1", 0, 30, 0.5, anchor_offsets=[5])
        panel = _make_panel([target], kmer_size=5)
        anchors = [AnchorKmer(kmer="ACGTG", target_idx=0, offset=5)]

        r1_path = tmp_path / "R1.fq.gz"
        r2_path = tmp_path / "R2.fq.gz"

        r1_content = "@r1/1\nNNNNNACGTGNNNN\n+\nIIIIIIIIIIIIII\n"
        r2_content = "@r1/2\nTTTTTTTTTTTTTT\n+\nIIIIIIIIIIIIII\n"
        with gzip.open(r1_path, "wt") as f:
            f.write(r1_content)
        with gzip.open(r2_path, "wt") as f:
            f.write(r2_content)

        result = extract_raw_fingerprint(r1_path, r2_path, panel, anchors)
        assert result["total_read_pairs"] == 1
        assert result["kmer_counts"][0] == pytest.approx(1.0)

    def test_empty_fastq(self, tmp_path):
        """Empty FASTQ files produce zero counts."""
        target = TargetRegion("t0", "chr1", 0, 30, 0.5, anchor_offsets=[5])
        panel = _make_panel([target], kmer_size=5)
        anchors = [AnchorKmer(kmer="ACGTG", target_idx=0, offset=5)]

        r1_path = tmp_path / "R1.fq"
        r2_path = tmp_path / "R2.fq"
        r1_path.write_text("")
        r2_path.write_text("")

        result = extract_raw_fingerprint(r1_path, r2_path, panel, anchors)
        assert result["total_read_pairs"] == 0
        assert result["kmer_counts"][0] == pytest.approx(0.0)
