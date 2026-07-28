"""Tests for k-mer utilities and anchor selection pipeline."""

from pathlib import Path

import numpy as np
import pytest

from shard.data.panel import AnchorKmer, PanelDefinition, SpacedSeedSet, TargetRegion
from shard.utils.kmer import (
    canonical_kmer,
    extract_kmers,
    filter_variant_sites,
    find_unique_target_kmers,
    gc_content,
    kmer_gc_content,
    kmer_overlaps_variant,
    reverse_complement,
    select_anchors_for_panel,
    select_spaced_anchors,
)


class TestReverseComplement:
    def test_basic(self):
        assert reverse_complement("ACGT") == "ACGT"  # palindrome
        assert reverse_complement("AAAA") == "TTTT"
        assert reverse_complement("ACGTA") == "TACGT"

    def test_preserves_length(self):
        seq = "ACGTACGTACGT"
        assert len(reverse_complement(seq)) == len(seq)

    def test_double_rc_is_identity(self):
        seq = "ACGTAACCGGTT"
        assert reverse_complement(reverse_complement(seq)) == seq


class TestCanonicalKmer:
    def test_returns_lexicographically_smaller(self):
        assert canonical_kmer("ACGTA") == "ACGTA"
        rc = reverse_complement("ACGTA")
        assert canonical_kmer(rc) == "ACGTA"

    def test_palindrome_unchanged(self):
        assert canonical_kmer("ACGT") == "ACGT"

    def test_strand_consistency(self):
        """A k-mer and its RC should produce the same canonical form."""
        kmer = "GATTACA"
        assert canonical_kmer(kmer) == canonical_kmer(reverse_complement(kmer))


class TestExtractKmers:
    def test_basic(self):
        kmers = extract_kmers("ACGTG", 3)
        assert kmers == [("ACG", 0), ("CGT", 1), ("GTG", 2)]

    def test_skips_n(self):
        kmers = extract_kmers("ACNGT", 3)
        # ACN and CNG contain N, only NGT... wait, NGT also has N at pos 0 of the kmer
        # Actually: ACN (skip), CNG (skip), NGT (skip)
        assert kmers == []

    def test_n_in_middle(self):
        kmers = extract_kmers("ACGTNACGT", 3)
        # ACG(0), CGT(1), GTN(skip), TNa(skip), NAC(skip), ACG(5), CGT(6)
        assert ("ACG", 0) in kmers
        assert ("CGT", 1) in kmers
        assert ("ACG", 5) in kmers
        assert ("CGT", 6) in kmers
        assert len(kmers) == 4

    def test_sequence_shorter_than_k(self):
        assert extract_kmers("AC", 3) == []

    def test_case_insensitive(self):
        kmers = extract_kmers("acgt", 3)
        assert kmers == [("ACG", 0), ("CGT", 1)]


class TestFindUniqueTargetKmers:
    def test_filters_non_unique(self):
        target = TargetRegion("t0", "chr1", 0, 10, 0.5)
        seq = "ACGTACGTAC"  # k=3: ACG, CGT, GTA, TAC, ACG, CGT, GTA, TAC
        # Note: GTA and TAC are reverse complements, so they share a canonical form.
        # canonical("ACG") = "ACG", canonical("CGT") = "ACG" (RC of CGT is ACG)
        # Actually let's verify: RC("CGT") = "ACG", so canonical("CGT") = "ACG"
        # This means ACG and CGT share a canonical form!
        # Use a sequence with clearly distinct canonical k-mers instead.
        seq = "AACCGGTT"  # k=3: AAC, ACC, CCG, CGG, GGT, GTT
        # canonical forms: AAC→AAC(vs GTT→AAC, same!), ACC→ACC(vs GGT→ACC, same!)
        # CCG→CCG(vs CGG→CCG, same!)
        # So unique canonicals: AAC, ACC, CCG — each appears twice in the seq
        # via forward + RC equivalence. Let's use a simpler approach:
        genome_counts = {
            canonical_kmer("AAC"): 2,   # not unique
            canonical_kmer("ACC"): 1,   # unique
            canonical_kmer("CCG"): 1,   # unique
        }
        unique = find_unique_target_kmers(target, 0, seq, genome_counts, k=3)
        kmer_seqs = {a.kmer for a in unique}
        # ACC (offset 1) and GGT (offset 4) both have canonical "ACC" with count 1
        # CCG (offset 2) and CGG (offset 3) both have canonical "CCG" with count 1
        assert "ACC" in kmer_seqs
        assert "CCG" in kmer_seqs
        # AAC (offset 0) and GTT (offset 5) have canonical "AAC" with count 2
        assert "AAC" not in kmer_seqs

    def test_all_unique(self):
        target = TargetRegion("t0", "chr1", 0, 7, 0.5)
        seq = "ACGTGCA"
        genome_counts = {canonical_kmer(k): 1 for k, _ in extract_kmers(seq, 3)}
        unique = find_unique_target_kmers(target, 0, seq, genome_counts, k=3)
        assert len(unique) == 5  # 7 - 3 + 1 = 5 k-mers

    def test_missing_from_counts_treated_as_zero(self):
        """K-mers not in genome_counts are treated as count=0 (not unique)."""
        target = TargetRegion("t0", "chr1", 0, 5, 0.5)
        seq = "ACGTG"
        genome_counts = {}  # empty — nothing is unique
        unique = find_unique_target_kmers(target, 0, seq, genome_counts, k=3)
        assert unique == []


class TestFilterVariantSites:
    def test_removes_overlapping_variants(self):
        target = TargetRegion("t0", "chr1", 100, 200, 0.5)
        candidates = [
            AnchorKmer("ACG", target_idx=0, offset=10),  # genomic 110-112
            AnchorKmer("CGT", target_idx=0, offset=50),  # genomic 150-152
        ]
        # Variant at position 111 overlaps the first k-mer
        variant_sites = {"chr1": {111}}
        filtered = filter_variant_sites(candidates, target, variant_sites, k=3)
        assert len(filtered) == 1
        assert filtered[0].offset == 50

    def test_no_variants_passes_all(self):
        target = TargetRegion("t0", "chr1", 0, 100, 0.5)
        candidates = [AnchorKmer("ACG", 0, 10), AnchorKmer("CGT", 0, 50)]
        filtered = filter_variant_sites(candidates, target, {}, k=3)
        assert len(filtered) == 2

    def test_variant_on_different_chrom(self):
        target = TargetRegion("t0", "chr1", 0, 100, 0.5)
        candidates = [AnchorKmer("ACG", 0, 10)]
        variant_sites = {"chr2": {10}}  # different chromosome
        filtered = filter_variant_sites(candidates, target, variant_sites, k=3)
        assert len(filtered) == 1


class TestKmerOverlapsVariant:
    def test_overlap_at_start(self):
        assert kmer_overlaps_variant("chr1", 100, 5, {"chr1": {100}})

    def test_overlap_at_end(self):
        assert kmer_overlaps_variant("chr1", 100, 5, {"chr1": {104}})

    def test_no_overlap_just_outside(self):
        assert not kmer_overlaps_variant("chr1", 100, 5, {"chr1": {105}})
        assert not kmer_overlaps_variant("chr1", 100, 5, {"chr1": {99}})


class TestSelectSpacedAnchors:
    def test_respects_min_spacing(self):
        candidates = [
            AnchorKmer("A" * 3, 0, offset=0),
            AnchorKmer("A" * 3, 0, offset=3),
            AnchorKmer("A" * 3, 0, offset=5),
            AnchorKmer("A" * 3, 0, offset=10),
            AnchorKmer("A" * 3, 0, offset=12),
        ]
        selected = select_spaced_anchors(candidates, target_length=20, min_spacing=5)
        offsets = [a.offset for a in selected]
        assert offsets == [0, 5, 10]

    def test_max_anchors(self):
        candidates = [AnchorKmer("A" * 3, 0, offset=i * 10) for i in range(20)]
        selected = select_spaced_anchors(
            candidates, target_length=200, min_spacing=5, max_anchors=5
        )
        assert len(selected) == 5

    def test_empty_candidates(self):
        assert select_spaced_anchors([], target_length=100) == []

    def test_single_candidate(self):
        candidates = [AnchorKmer("ACG", 0, offset=50)]
        selected = select_spaced_anchors(candidates, target_length=100)
        assert len(selected) == 1

    def test_unsorted_input(self):
        """Candidates don't need to be pre-sorted."""
        candidates = [
            AnchorKmer("A" * 3, 0, offset=50),
            AnchorKmer("A" * 3, 0, offset=0),
            AnchorKmer("A" * 3, 0, offset=100),
        ]
        selected = select_spaced_anchors(candidates, target_length=200, min_spacing=10)
        offsets = [a.offset for a in selected]
        assert offsets == [0, 50, 100]


class TestSelectAnchorsForPanel:
    def test_full_pipeline(self):
        panel = PanelDefinition(
            panel_id="test",
            kmer_size=3,
            targets=[
                TargetRegion("t0", "chr1", 0, 10, 0.5),
                TargetRegion("t1", "chr2", 100, 112, 0.5),
            ],
            seed_set=SpacedSeedSet.exact_match(3),
        )
        seqs = ["ACGTGCATGC", "TTAAGCCTTAAG"]

        # Build genome counts where all k-mers are unique
        genome_counts: dict[str, int] = {}
        for seq in seqs:
            for kmer, _ in extract_kmers(seq, 3):
                canon = canonical_kmer(kmer)
                genome_counts[canon] = 1  # pretend all are unique

        anchors, low_density = select_anchors_for_panel(
            panel, seqs, genome_counts, min_spacing=1, min_anchors_per_target=1
        )

        assert len(anchors) > 0
        assert all(isinstance(a, AnchorKmer) for a in anchors)
        # Each target should have updated anchor_offsets
        assert len(panel.targets[0].anchor_offsets) > 0
        assert len(panel.targets[1].anchor_offsets) > 0

    def test_low_density_warning(self):
        panel = PanelDefinition(
            panel_id="test",
            kmer_size=3,
            targets=[TargetRegion("t0", "chr1", 0, 5, 0.5)],
            seed_set=SpacedSeedSet.exact_match(3),
        )
        # All k-mers non-unique
        genome_counts = {canonical_kmer("ACG"): 5, canonical_kmer("CGT"): 5, canonical_kmer("GTG"): 5}
        anchors, low_density = select_anchors_for_panel(
            panel, ["ACGTG"], genome_counts, min_anchors_per_target=3
        )
        assert 0 in low_density  # target 0 flagged as low density
        assert len(anchors) == 0

    def test_variant_filtering_integrated(self):
        panel = PanelDefinition(
            panel_id="test",
            kmer_size=3,
            targets=[TargetRegion("t0", "chr1", 100, 110, 0.5)],
            seed_set=SpacedSeedSet.exact_match(3),
        )
        seq = "ACGTGCATGC"
        genome_counts = {}
        for kmer, _ in extract_kmers(seq, 3):
            genome_counts[canonical_kmer(kmer)] = 1

        # Place a variant that overlaps every k-mer position (pos 103 = offset 3)
        # This affects k-mers at offsets 1, 2, 3 (any k-mer spanning position 103)
        variant_sites = {"chr1": {103}}

        anchors_with, _ = select_anchors_for_panel(
            panel, [seq], genome_counts, variant_sites=variant_sites, min_spacing=1
        )
        anchors_without, _ = select_anchors_for_panel(
            panel, [seq], genome_counts, variant_sites=None, min_spacing=1
        )
        assert len(anchors_with) < len(anchors_without)


class TestGCContent:
    def test_basic(self):
        assert gc_content("GCGC") == pytest.approx(1.0)
        assert gc_content("ATAT") == pytest.approx(0.0)
        assert gc_content("ACGT") == pytest.approx(0.5)

    def test_ignores_n(self):
        assert gc_content("GCNA") == pytest.approx(2 / 3)

    def test_empty(self):
        assert gc_content("NNN") == 0.0

    def test_kmer_gc(self):
        assert kmer_gc_content("GCGC") == pytest.approx(1.0)
