"""Fingerprint generation from FASTQ files.

Computes per-sample fingerprints: anchor k-mer counts, GC distributions,
fragment size distributions, and per-target joint profiles.

The primary representation is a 2D matrix (N_targets, N_features) where each
row is a panel target in genomic order and columns are:

    [0]      kmer_count         - median anchor k-mer count (raw, unnormalised)
    [1]      kmer_count_ratio   - sample kmer count / run median (batch-normalised RD signal)
    [2]      kmer_count_zscore  - z-score within run (per target)
    [3:23]   gc_profile         - per-target GC distribution (20 bins)
    [23:73]  frag_profile       - per-target fragment size distribution (50 bins)
    [73]     gc_deviation       - per-target GC profile distance from run median
    [74]     frag_deviation     - per-target fragment profile distance from run median
    [75]     target_gc          - target region GC content (from panel design)
    [76]     target_length      - target region length in bp
    [77]     anchor_density     - unique anchor k-mers / target length

See research/fingerprint_storage_formats.md for format rationale.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray


# Column layout constants
GC_BINS: int = 20
FRAG_BINS: int = 50
N_BATCH_FEATURES: int = 4  # kmer_count_ratio, kmer_count_zscore, gc_deviation, frag_deviation
N_SCALARS: int = 3
N_FEATURES: int = 1 + N_BATCH_FEATURES + GC_BINS + FRAG_BINS + N_SCALARS  # 78

# Column slices for indexing into the feature axis
KMER_COUNT_COL: int = 0
KMER_RATIO_COL: int = 1
KMER_ZSCORE_COL: int = 2
GC_PROFILE_SLICE: slice = slice(3, 3 + GC_BINS)
FRAG_PROFILE_SLICE: slice = slice(3 + GC_BINS, 3 + GC_BINS + FRAG_BINS)
GC_DEVIATION_COL: int = 3 + GC_BINS + FRAG_BINS
FRAG_DEVIATION_COL: int = GC_DEVIATION_COL + 1
TARGET_GC_COL: int = GC_DEVIATION_COL + 2
TARGET_LENGTH_COL: int = TARGET_GC_COL + 1
ANCHOR_DENSITY_COL: int = TARGET_GC_COL + 2

# Global distribution sizes (stored separately, not per-target)
GLOBAL_GC_BINS: int = 100
GLOBAL_FRAG_RANGE: tuple[int, int] = (50, 500)
GLOBAL_FRAG_BINS: int = GLOBAL_FRAG_RANGE[1] - GLOBAL_FRAG_RANGE[0]  # 450


@dataclass
class RunContext:
    """Summary statistics from the sequencing run a sample belongs to.

    Captures run-level distributions so that each sample can be
    batch-normalised independently at training/inference time without
    requiring all batch-mates to be present.

    Attributes:
        run_id: Unique identifier for the sequencing run.
        sequencer_id: Instrument identifier (e.g. "NB501234").
        run_date: ISO-format date string (YYYY-MM-DD).
        n_samples_in_run: Number of samples in the run.
        median_total_reads: Median total read count across samples in the run.
        median_kmer_counts: Per-target median k-mer counts across the run,
            shape (n_targets,). The primary reference for RD ratio calculation.
        std_kmer_counts: Per-target std dev of k-mer counts across the run,
            shape (n_targets,). Used for z-score calculation.
        median_gc_histogram: Run-wide median per-read GC distribution,
            shape (100,).
        median_frag_histogram: Run-wide median fragment size distribution,
            shape (450,).
        median_gc_profiles: Run-wide median per-target GC profiles,
            shape (n_targets, 20).
        median_frag_profiles: Run-wide median per-target fragment profiles,
            shape (n_targets, 50).
    """

    run_id: str
    sequencer_id: str
    run_date: str
    n_samples_in_run: int
    median_total_reads: float
    median_kmer_counts: NDArray[np.float32]
    std_kmer_counts: NDArray[np.float32]
    median_gc_histogram: NDArray[np.float32]
    median_frag_histogram: NDArray[np.float32]
    median_gc_profiles: NDArray[np.float32]
    median_frag_profiles: NDArray[np.float32]


@dataclass
class SampleFingerprint:
    """Alignment-free fingerprint for a single gene panel sample.

    The core representation is a 2D feature matrix of shape (n_targets, 78)
    containing per-target k-mer counts, batch-normalised ratios, distribution
    profiles, and metadata.  Global (sample-level) distributions are stored
    alongside for QC and bias correction.

    Raw k-mer counts and distributions are stored unnormalised. Batch-relative
    features (ratio, z-score, deviations) are precomputed from RunContext at
    fingerprinting time because the run median is a sufficient statistic for
    the dominant batch effects (depth, GC bias, fragment size shift).

    Attributes:
        sample_id: Unique sample identifier.
        panel_id: Panel design identifier (ties to anchor k-mer set).
        target_names: Ordered target names, length n_targets.
        features: Per-target feature matrix, shape (n_targets, 78).
            Column layout:
                [0]      median anchor k-mer count (raw)
                [1]      kmer count ratio (sample / run median)
                [2]      kmer count z-score within run
                [3:23]   GC profile (20 bins)
                [23:73]  fragment size profile (50 bins)
                [73]     GC deviation (L2 dist from run median GC profile)
                [74]     fragment deviation (L2 dist from run median frag profile)
                [75]     target GC content
                [76]     target length (bp)
                [77]     anchor density (unique k-mers / target length)
        global_gc_histogram: Sample-wide per-read GC distribution, shape (100,).
        global_frag_histogram: Sample-wide fragment size distribution, shape (450,).
        total_reads: Total read count in the sample (for depth normalisation).
        run_context: Summary statistics from the sequencing run this sample
            belongs to, used for batch normalisation.
        cn_labels: Optional ground-truth copy number per target, shape (n_targets,).
            Values: 0, 1, 2, 3, 4, ... or -1 for unknown.
    """

    sample_id: str
    panel_id: str
    target_names: list[str]
    features: NDArray[np.float32]
    global_gc_histogram: NDArray[np.float32]
    global_frag_histogram: NDArray[np.float32]
    total_reads: int
    run_context: RunContext
    cn_labels: NDArray[np.int8] | None = None

    def __post_init__(self) -> None:
        n_targets = len(self.target_names)
        if self.features.shape != (n_targets, N_FEATURES):
            raise ValueError(
                f"features shape must be ({n_targets}, {N_FEATURES}), "
                f"got {self.features.shape}"
            )
        if self.global_gc_histogram.shape != (GLOBAL_GC_BINS,):
            raise ValueError(
                f"global_gc_histogram shape must be ({GLOBAL_GC_BINS},), "
                f"got {self.global_gc_histogram.shape}"
            )
        if self.global_frag_histogram.shape != (GLOBAL_FRAG_BINS,):
            raise ValueError(
                f"global_frag_histogram shape must be ({GLOBAL_FRAG_BINS},), "
                f"got {self.global_frag_histogram.shape}"
            )
        if self.cn_labels is not None and self.cn_labels.shape != (n_targets,):
            raise ValueError(
                f"cn_labels shape must be ({n_targets},), "
                f"got {self.cn_labels.shape}"
            )

    @property
    def n_targets(self) -> int:
        return len(self.target_names)

    @property
    def kmer_counts(self) -> NDArray[np.float32]:
        """Per-target median anchor k-mer counts, shape (n_targets,)."""
        return self.features[:, KMER_COUNT_COL]

    @property
    def kmer_ratios(self) -> NDArray[np.float32]:
        """Per-target kmer count / run median ratio (batch-normalised RD signal), shape (n_targets,)."""
        return self.features[:, KMER_RATIO_COL]

    @property
    def kmer_zscores(self) -> NDArray[np.float32]:
        """Per-target kmer count z-score within run, shape (n_targets,)."""
        return self.features[:, KMER_ZSCORE_COL]

    @property
    def gc_profiles(self) -> NDArray[np.float32]:
        """Per-target GC distributions, shape (n_targets, 20)."""
        return self.features[:, GC_PROFILE_SLICE]

    @property
    def frag_profiles(self) -> NDArray[np.float32]:
        """Per-target fragment size distributions, shape (n_targets, 50)."""
        return self.features[:, FRAG_PROFILE_SLICE]

    @property
    def gc_deviations(self) -> NDArray[np.float32]:
        """Per-target GC profile L2 distance from run median, shape (n_targets,)."""
        return self.features[:, GC_DEVIATION_COL]

    @property
    def frag_deviations(self) -> NDArray[np.float32]:
        """Per-target fragment profile L2 distance from run median, shape (n_targets,)."""
        return self.features[:, FRAG_DEVIATION_COL]

    @property
    def target_gc(self) -> NDArray[np.float32]:
        """Per-target GC content from panel design, shape (n_targets,)."""
        return self.features[:, TARGET_GC_COL]

    @property
    def target_lengths(self) -> NDArray[np.float32]:
        """Per-target region lengths in bp, shape (n_targets,)."""
        return self.features[:, TARGET_LENGTH_COL]

    @property
    def anchor_densities(self) -> NDArray[np.float32]:
        """Per-target anchor density (unique k-mers / length), shape (n_targets,)."""
        return self.features[:, ANCHOR_DENSITY_COL]


@dataclass
class CohortFingerprints:
    """Collection of fingerprints for a cohort, pre-stacked for ML training.

    All samples must share the same panel (same targets in the same order).
    Samples may come from different sequencing runs; ``run_ids`` tracks which
    run each sample belongs to, enabling run-stratified train/test splits
    (critical to avoid leaking batch structure into evaluation).

    Attributes:
        panel_id: Panel design identifier.
        target_names: Ordered target names, length n_targets.
        sample_ids: Sample identifiers, length n_samples.
        run_ids: Per-sample run identifier, length n_samples.  Used for
            stratified splitting — never split samples from the same run
            across train and test.
        run_contexts: Mapping from run_id to its RunContext.
        features: Stacked feature matrices, shape (n_samples, n_targets, 78).
        global_gc_histograms: shape (n_samples, 100).
        global_frag_histograms: shape (n_samples, 450).
        total_reads: Per-sample read counts, shape (n_samples,).
        cn_labels: Optional ground-truth CN states, shape (n_samples, n_targets).
            -1 indicates unknown.
    """

    panel_id: str
    target_names: list[str]
    sample_ids: list[str]
    run_ids: list[str]
    run_contexts: dict[str, RunContext]
    features: NDArray[np.float32]
    global_gc_histograms: NDArray[np.float32]
    global_frag_histograms: NDArray[np.float32]
    total_reads: NDArray[np.int64]
    cn_labels: NDArray[np.int8] | None = None

    def __post_init__(self) -> None:
        n_samples = len(self.sample_ids)
        n_targets = len(self.target_names)
        if len(self.run_ids) != n_samples:
            raise ValueError(
                f"run_ids length must be {n_samples}, got {len(self.run_ids)}"
            )
        expected_feat = (n_samples, n_targets, N_FEATURES)
        if self.features.shape != expected_feat:
            raise ValueError(
                f"features shape must be {expected_feat}, got {self.features.shape}"
            )
        if self.global_gc_histograms.shape != (n_samples, GLOBAL_GC_BINS):
            raise ValueError(
                f"global_gc_histograms shape must be ({n_samples}, {GLOBAL_GC_BINS}), "
                f"got {self.global_gc_histograms.shape}"
            )
        if self.global_frag_histograms.shape != (n_samples, GLOBAL_FRAG_BINS):
            raise ValueError(
                f"global_frag_histograms shape must be ({n_samples}, {GLOBAL_FRAG_BINS}), "
                f"got {self.global_frag_histograms.shape}"
            )
        if self.cn_labels is not None and self.cn_labels.shape != (n_samples, n_targets):
            raise ValueError(
                f"cn_labels shape must be ({n_samples}, {n_targets}), "
                f"got {self.cn_labels.shape}"
            )

    @property
    def n_samples(self) -> int:
        return len(self.sample_ids)

    @property
    def n_targets(self) -> int:
        return len(self.target_names)

    @property
    def unique_run_ids(self) -> list[str]:
        """Distinct run IDs in order of first appearance."""
        seen: set[str] = set()
        result: list[str] = []
        for rid in self.run_ids:
            if rid not in seen:
                seen.add(rid)
                result.append(rid)
        return result

    def get_sample(self, idx: int) -> SampleFingerprint:
        """Extract a single SampleFingerprint by index."""
        run_id = self.run_ids[idx]
        return SampleFingerprint(
            sample_id=self.sample_ids[idx],
            panel_id=self.panel_id,
            target_names=self.target_names,
            features=self.features[idx],
            global_gc_histogram=self.global_gc_histograms[idx],
            global_frag_histogram=self.global_frag_histograms[idx],
            total_reads=int(self.total_reads[idx]),
            run_context=self.run_contexts[run_id],
            cn_labels=self.cn_labels[idx] if self.cn_labels is not None else None,
        )

    @classmethod
    def from_samples(cls, samples: list[SampleFingerprint]) -> CohortFingerprints:
        """Stack a list of SampleFingerprints into a CohortFingerprints."""
        if not samples:
            raise ValueError("Cannot create CohortFingerprints from empty list")
        panel_id = samples[0].panel_id
        target_names = samples[0].target_names
        for s in samples[1:]:
            if s.panel_id != panel_id:
                raise ValueError(
                    f"All samples must share panel_id: expected {panel_id!r}, "
                    f"got {s.panel_id!r} for sample {s.sample_id!r}"
                )
            if s.target_names != target_names:
                raise ValueError(
                    f"Target name mismatch for sample {s.sample_id!r}"
                )
        has_labels = all(s.cn_labels is not None for s in samples)
        run_contexts: dict[str, RunContext] = {}
        for s in samples:
            rid = s.run_context.run_id
            if rid not in run_contexts:
                run_contexts[rid] = s.run_context
        return cls(
            panel_id=panel_id,
            target_names=target_names,
            sample_ids=[s.sample_id for s in samples],
            run_ids=[s.run_context.run_id for s in samples],
            run_contexts=run_contexts,
            features=np.stack([s.features for s in samples]),
            global_gc_histograms=np.stack(
                [s.global_gc_histogram for s in samples]
            ),
            global_frag_histograms=np.stack(
                [s.global_frag_histogram for s in samples]
            ),
            total_reads=np.array([s.total_reads for s in samples], dtype=np.int64),
            cn_labels=(
                np.stack([s.cn_labels for s in samples]) if has_labels else None
            ),
        )


def build_run_context(
    run_id: str,
    sequencer_id: str,
    run_date: str,
    sample_kmer_counts: NDArray[np.float32],
    sample_total_reads: NDArray[np.int64],
    sample_gc_histograms: NDArray[np.float32],
    sample_frag_histograms: NDArray[np.float32],
    sample_gc_profiles: NDArray[np.float32],
    sample_frag_profiles: NDArray[np.float32],
) -> RunContext:
    """Build a RunContext from all samples in a sequencing run.

    Args:
        run_id: Unique run identifier.
        sequencer_id: Instrument identifier.
        run_date: ISO-format date string.
        sample_kmer_counts: shape (n_samples, n_targets).
        sample_total_reads: shape (n_samples,).
        sample_gc_histograms: shape (n_samples, 100).
        sample_frag_histograms: shape (n_samples, 450).
        sample_gc_profiles: shape (n_samples, n_targets, 20).
        sample_frag_profiles: shape (n_samples, n_targets, 50).
    """
    return RunContext(
        run_id=run_id,
        sequencer_id=sequencer_id,
        run_date=run_date,
        n_samples_in_run=sample_kmer_counts.shape[0],
        median_total_reads=float(np.median(sample_total_reads)),
        median_kmer_counts=np.median(sample_kmer_counts, axis=0),
        std_kmer_counts=np.std(sample_kmer_counts, axis=0).astype(np.float32),
        median_gc_histogram=np.median(sample_gc_histograms, axis=0).astype(np.float32),
        median_frag_histogram=np.median(sample_frag_histograms, axis=0).astype(np.float32),
        median_gc_profiles=np.median(sample_gc_profiles, axis=0).astype(np.float32),
        median_frag_profiles=np.median(sample_frag_profiles, axis=0).astype(np.float32),
    )


def compute_batch_features(
    kmer_counts: NDArray[np.float32],
    gc_profiles: NDArray[np.float32],
    frag_profiles: NDArray[np.float32],
    run_ctx: RunContext,
) -> tuple[NDArray[np.float32], NDArray[np.float32], NDArray[np.float32], NDArray[np.float32]]:
    """Compute batch-relative features for a single sample.

    Args:
        kmer_counts: Raw per-target k-mer counts, shape (n_targets,).
        gc_profiles: Per-target GC distributions, shape (n_targets, 20).
        frag_profiles: Per-target fragment distributions, shape (n_targets, 50).
        run_ctx: RunContext with run-level summary statistics.

    Returns:
        Tuple of (kmer_ratios, kmer_zscores, gc_deviations, frag_deviations),
        each shape (n_targets,).
    """
    # Ratio: sample / run_median (the primary batch-normalised RD signal)
    # Add small epsilon to avoid division by zero for targets with no counts
    eps = np.float32(1e-6)
    kmer_ratios = kmer_counts / (run_ctx.median_kmer_counts + eps)

    # Z-score: (sample - median) / std within run
    std = run_ctx.std_kmer_counts + eps
    kmer_zscores = (kmer_counts - run_ctx.median_kmer_counts) / std

    # Per-target L2 distance of GC profile from run median
    gc_diff = gc_profiles - run_ctx.median_gc_profiles
    gc_deviations = np.sqrt(np.sum(gc_diff ** 2, axis=1))

    # Per-target L2 distance of fragment profile from run median
    frag_diff = frag_profiles - run_ctx.median_frag_profiles
    frag_deviations = np.sqrt(np.sum(frag_diff ** 2, axis=1))

    return (
        kmer_ratios.astype(np.float32),
        kmer_zscores.astype(np.float32),
        gc_deviations.astype(np.float32),
        frag_deviations.astype(np.float32),
    )
