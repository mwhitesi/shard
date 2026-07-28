"""K-mer extraction pipeline: FASTQ reads → anchor hits → raw fingerprint data.

This module connects raw sequencing data to the fingerprint representation.
It streams paired-end FASTQ files, extracts k-mers from each read, looks
them up against the panel's AnchorIndex (with spaced-seed SNP tolerance),
and accumulates per-target statistics:

    - Per-anchor k-mer counts (for median-based depth estimation)
    - Per-target GC profiles (from anchor-containing reads)
    - Per-target fragment size distributions (from R1/R2 anchor offsets)
    - Global GC histogram and fragment size histogram

The output is a dict of raw numpy arrays ready for SampleFingerprint assembly
(after batch normalisation via RunContext).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from shard.data.fingerprint import (
    FRAG_BINS,
    GC_BINS,
    GLOBAL_FRAG_BINS,
    GLOBAL_FRAG_RANGE,
    GLOBAL_GC_BINS,
)
from shard.data.panel import (
    AnchorHit,
    AnchorIndex,
    AnchorKmer,
    PanelDefinition,
)
from shard.utils.fastq import FastqRecord, read_paired_fastq
from shard.utils.kmer import canonical_kmer, extract_kmers, gc_content


# ---------------------------------------------------------------------------
# Per-read anchor scanning
# ---------------------------------------------------------------------------


def scan_read(
    sequence: str,
    index: AnchorIndex,
    k: int,
) -> list[AnchorHit]:
    """Extract k-mers from a read and look each up in the anchor index.

    Args:
        sequence: DNA sequence from a FASTQ read (uppercase).
        index: Panel anchor index (supports spaced-seed lookup).
        k: K-mer length (must match panel kmer_size).

    Returns:
        List of AnchorHit for every k-mer that matches a panel anchor.
    """
    hits: list[AnchorHit] = []
    for kmer, read_pos in extract_kmers(sequence, k):
        # Try the k-mer as-is first, then its canonical form.
        # The index stores keys from the original anchor k-mer sequence,
        # but reads can come from either strand. We need to try both
        # the forward k-mer and its reverse complement.
        matches = index.lookup(kmer)
        if not matches:
            rc = canonical_kmer(kmer)
            # Only retry if canonical differs (i.e., we haven't tried the RC yet)
            if rc != kmer:
                from shard.utils.kmer import reverse_complement

                matches = index.lookup(reverse_complement(kmer))
        for target_idx, offset in matches:
            hits.append(
                AnchorHit(
                    target_idx=target_idx,
                    offset=offset,
                    read_position=read_pos,
                )
            )
    return hits


# ---------------------------------------------------------------------------
# Sample-level accumulator
# ---------------------------------------------------------------------------


class SampleAccumulator:
    """Accumulates per-target statistics from all read pairs in a sample.

    Tracks:
        - Per-anchor hit counts (individual k-mer level, for median computation)
        - Per-target GC profiles (GC content of reads with anchor hits)
        - Per-target fragment size estimates (from R1/R2 anchor offsets)
        - Global GC histogram (all reads)
        - Global fragment size histogram (all pairs with estimates)
        - Total read pair count
    """

    def __init__(self, panel: PanelDefinition) -> None:
        self.panel = panel
        n_targets = panel.n_targets

        # Per-anchor k-mer hit counts: for each target, for each anchor offset,
        # count how many reads contained that anchor.
        # Store as dict: (target_idx, offset) → count
        self._anchor_counts: dict[tuple[int, int], int] = {}

        # Per-target GC profile: binned GC content of reads hitting each target.
        # Shape: (n_targets, GC_BINS)
        self.gc_profiles = np.zeros((n_targets, GC_BINS), dtype=np.float64)

        # Per-target fragment size distribution.
        # Shape: (n_targets, FRAG_BINS)
        # Bins cover a subsampled range for the per-target profile (50 bins).
        self._frag_min = GLOBAL_FRAG_RANGE[0]
        self._frag_max = GLOBAL_FRAG_RANGE[1]
        self.frag_profiles = np.zeros((n_targets, FRAG_BINS), dtype=np.float64)

        # Global GC histogram (all reads, 100 bins).
        self.global_gc_histogram = np.zeros(GLOBAL_GC_BINS, dtype=np.float64)

        # Global fragment size histogram (450 bins).
        self.global_frag_histogram = np.zeros(GLOBAL_FRAG_BINS, dtype=np.float64)

        self.total_read_pairs: int = 0

    def process_pair(
        self,
        r1: FastqRecord,
        r2: FastqRecord,
        index: AnchorIndex,
    ) -> None:
        """Process a single read pair, updating all accumulators.

        Args:
            r1: Forward read (R1).
            r2: Reverse read (R2).
            index: Panel anchor index.
        """
        k = self.panel.kmer_size
        self.total_read_pairs += 1

        # --- Global GC ---
        r1_gc = gc_content(r1.sequence)
        r2_gc = gc_content(r2.sequence)
        self._add_to_gc_histogram(r1_gc, self.global_gc_histogram, GLOBAL_GC_BINS)
        self._add_to_gc_histogram(r2_gc, self.global_gc_histogram, GLOBAL_GC_BINS)

        # --- Anchor scanning ---
        r1_hits = scan_read(r1.sequence, index, k)
        r2_hits = scan_read(r2.sequence, index, k)

        # --- Per-anchor counts ---
        for hit in r1_hits:
            key = (hit.target_idx, hit.offset)
            self._anchor_counts[key] = self._anchor_counts.get(key, 0) + 1
        for hit in r2_hits:
            key = (hit.target_idx, hit.offset)
            self._anchor_counts[key] = self._anchor_counts.get(key, 0) + 1

        # --- Per-target GC profiles ---
        # Record the read's GC for each target it hits
        r1_targets = {h.target_idx for h in r1_hits}
        r2_targets = {h.target_idx for h in r2_hits}
        for t in r1_targets:
            self._add_to_gc_profile(r1_gc, t)
        for t in r2_targets:
            self._add_to_gc_profile(r2_gc, t)

        # --- Fragment size estimation ---
        self._estimate_fragments(r1_hits, r2_hits)

    def _add_to_gc_histogram(
        self, gc_frac: float, histogram: NDArray, n_bins: int
    ) -> None:
        """Add a GC fraction to a histogram."""
        bin_idx = int(gc_frac * n_bins)
        if bin_idx >= n_bins:
            bin_idx = n_bins - 1
        if bin_idx >= 0:
            histogram[bin_idx] += 1

    def _add_to_gc_profile(self, gc_frac: float, target_idx: int) -> None:
        """Add a GC fraction to a per-target GC profile."""
        bin_idx = int(gc_frac * GC_BINS)
        if bin_idx >= GC_BINS:
            bin_idx = GC_BINS - 1
        if bin_idx >= 0:
            self.gc_profiles[target_idx, bin_idx] += 1

    def _estimate_fragments(
        self,
        r1_hits: list[AnchorHit],
        r2_hits: list[AnchorHit],
    ) -> None:
        """Estimate fragment sizes from anchor hit pairs and update distributions."""
        # Group by target
        r1_by_target: dict[int, list[AnchorHit]] = {}
        for h in r1_hits:
            r1_by_target.setdefault(h.target_idx, []).append(h)

        r2_by_target: dict[int, list[AnchorHit]] = {}
        for h in r2_hits:
            r2_by_target.setdefault(h.target_idx, []).append(h)

        for target_idx in r1_by_target:
            if target_idx not in r2_by_target:
                continue
            for h1 in r1_by_target[target_idx]:
                for h2 in r2_by_target[target_idx]:
                    frag_size = abs(h1.offset - h2.offset) + 150  # read_length default
                    self._add_fragment_estimate(target_idx, frag_size)

    def _add_fragment_estimate(self, target_idx: int, frag_size: int) -> None:
        """Add a fragment size estimate to both per-target and global distributions."""
        # Global histogram (1-bp bins)
        global_bin = frag_size - self._frag_min
        if 0 <= global_bin < GLOBAL_FRAG_BINS:
            self.global_frag_histogram[global_bin] += 1

        # Per-target profile (coarser bins: FRAG_BINS bins over the range)
        bin_width = (self._frag_max - self._frag_min) / FRAG_BINS
        target_bin = int((frag_size - self._frag_min) / bin_width)
        if 0 <= target_bin < FRAG_BINS:
            self.frag_profiles[target_idx, target_bin] += 1

    def compute_kmer_counts(self) -> NDArray[np.float32]:
        """Compute per-target median anchor k-mer counts.

        For each target, groups anchor counts by the anchors belonging to
        that target and takes the median. Targets with no anchor hits get 0.

        Returns:
            Array of shape (n_targets,) with median k-mer count per target.
        """
        n_targets = self.panel.n_targets
        kmer_counts = np.zeros(n_targets, dtype=np.float32)

        # Group counts by target
        target_counts: dict[int, list[int]] = {}
        for (target_idx, offset), count in self._anchor_counts.items():
            target_counts.setdefault(target_idx, []).append(count)

        # Also include zeros for anchors that were never hit
        for idx, target in enumerate(self.panel.targets):
            hit_offsets = {
                off for (t, off) in self._anchor_counts if t == idx
            }
            n_missing = max(0, target.n_anchor_kmers - len(hit_offsets))
            counts = target_counts.get(idx, [])
            counts.extend([0] * n_missing)
            if counts:
                kmer_counts[idx] = float(np.median(counts))

        return kmer_counts

    def normalise_profiles(self) -> tuple[NDArray[np.float32], NDArray[np.float32]]:
        """Normalise GC and fragment profiles to probability distributions.

        Returns:
            Tuple of (gc_profiles, frag_profiles) both as float32,
            each row summing to 1 (or all zeros if no data).
        """
        gc = self.gc_profiles.copy()
        frag = self.frag_profiles.copy()

        gc_sums = gc.sum(axis=1, keepdims=True)
        gc_sums[gc_sums == 0] = 1.0
        gc = (gc / gc_sums).astype(np.float32)

        frag_sums = frag.sum(axis=1, keepdims=True)
        frag_sums[frag_sums == 0] = 1.0
        frag = (frag / frag_sums).astype(np.float32)

        return gc, frag

    def normalise_global_histograms(
        self,
    ) -> tuple[NDArray[np.float32], NDArray[np.float32]]:
        """Normalise global GC and fragment histograms to probability distributions.

        Returns:
            Tuple of (gc_histogram, frag_histogram) as float32.
        """
        gc = self.global_gc_histogram.copy()
        gc_total = gc.sum()
        if gc_total > 0:
            gc /= gc_total

        frag = self.global_frag_histogram.copy()
        frag_total = frag.sum()
        if frag_total > 0:
            frag /= frag_total

        return gc.astype(np.float32), frag.astype(np.float32)


# ---------------------------------------------------------------------------
# Top-level extraction
# ---------------------------------------------------------------------------


def extract_raw_fingerprint(
    r1_path: Path,
    r2_path: Path,
    panel: PanelDefinition,
    anchors: list[AnchorKmer],
) -> dict:
    """Extract raw fingerprint data from paired FASTQ files.

    Streams reads, scans for anchor k-mers, and accumulates all per-target
    statistics needed to build a SampleFingerprint (before batch normalisation).

    Args:
        r1_path: Path to R1 FASTQ file (plain or gzipped).
        r2_path: Path to R2 FASTQ file (plain or gzipped).
        panel: Panel definition with targets and spaced seeds.
        anchors: All anchor k-mers for the panel.

    Returns:
        Dictionary with keys:
            - "kmer_counts": NDArray[float32], shape (n_targets,)
            - "gc_profiles": NDArray[float32], shape (n_targets, 20)
            - "frag_profiles": NDArray[float32], shape (n_targets, 50)
            - "global_gc_histogram": NDArray[float32], shape (100,)
            - "global_frag_histogram": NDArray[float32], shape (450,)
            - "total_read_pairs": int
    """
    index = AnchorIndex(panel, anchors)
    acc = SampleAccumulator(panel)

    for r1, r2 in read_paired_fastq(r1_path, r2_path):
        acc.process_pair(r1, r2, index)

    kmer_counts = acc.compute_kmer_counts()
    gc_profiles, frag_profiles = acc.normalise_profiles()
    global_gc, global_frag = acc.normalise_global_histograms()

    return {
        "kmer_counts": kmer_counts,
        "gc_profiles": gc_profiles,
        "frag_profiles": frag_profiles,
        "global_gc_histogram": global_gc,
        "global_frag_histogram": global_frag,
        "total_read_pairs": acc.total_read_pairs,
    }
