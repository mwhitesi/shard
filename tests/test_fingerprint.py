"""Tests for fingerprint data model with batch-effect handling."""

import numpy as np
import pytest

from shard.data.fingerprint import (
    ANCHOR_DENSITY_COL,
    FRAG_DEVIATION_COL,
    FRAG_PROFILE_SLICE,
    GC_DEVIATION_COL,
    GC_PROFILE_SLICE,
    GLOBAL_FRAG_BINS,
    GLOBAL_GC_BINS,
    KMER_COUNT_COL,
    KMER_RATIO_COL,
    KMER_ZSCORE_COL,
    N_FEATURES,
    TARGET_GC_COL,
    TARGET_LENGTH_COL,
    CohortFingerprints,
    RunContext,
    SampleFingerprint,
    build_run_context,
    compute_batch_features,
)

N_TARGETS = 5


def _make_run_context(n_targets: int = N_TARGETS) -> RunContext:
    return RunContext(
        run_id="run_001",
        sequencer_id="NB501234",
        run_date="2026-01-15",
        n_samples_in_run=8,
        median_total_reads=1_000_000.0,
        median_kmer_counts=np.full(n_targets, 500.0, dtype=np.float32),
        std_kmer_counts=np.full(n_targets, 50.0, dtype=np.float32),
        median_gc_histogram=np.ones(GLOBAL_GC_BINS, dtype=np.float32) / GLOBAL_GC_BINS,
        median_frag_histogram=np.ones(GLOBAL_FRAG_BINS, dtype=np.float32) / GLOBAL_FRAG_BINS,
        median_gc_profiles=np.ones((n_targets, 20), dtype=np.float32) / 20,
        median_frag_profiles=np.ones((n_targets, 50), dtype=np.float32) / 50,
    )


def _make_sample(
    sample_id: str = "sample_A",
    run_context: RunContext | None = None,
    cn_labels: np.ndarray | None = None,
) -> SampleFingerprint:
    run_ctx = run_context or _make_run_context()
    return SampleFingerprint(
        sample_id=sample_id,
        panel_id="panel_v1",
        target_names=[f"target_{i}" for i in range(N_TARGETS)],
        features=np.zeros((N_TARGETS, N_FEATURES), dtype=np.float32),
        global_gc_histogram=np.zeros(GLOBAL_GC_BINS, dtype=np.float32),
        global_frag_histogram=np.zeros(GLOBAL_FRAG_BINS, dtype=np.float32),
        total_reads=500_000,
        run_context=run_ctx,
        cn_labels=cn_labels,
    )


class TestColumnLayout:
    def test_n_features_is_78(self):
        assert N_FEATURES == 78

    def test_column_indices_are_consistent(self):
        assert KMER_COUNT_COL == 0
        assert KMER_RATIO_COL == 1
        assert KMER_ZSCORE_COL == 2
        assert GC_PROFILE_SLICE == slice(3, 23)
        assert FRAG_PROFILE_SLICE == slice(23, 73)
        assert GC_DEVIATION_COL == 73
        assert FRAG_DEVIATION_COL == 74
        assert TARGET_GC_COL == 75
        assert TARGET_LENGTH_COL == 76
        assert ANCHOR_DENSITY_COL == 77

    def test_slices_dont_overlap(self):
        used = set()
        cols = [KMER_COUNT_COL, KMER_RATIO_COL, KMER_ZSCORE_COL,
                GC_DEVIATION_COL, FRAG_DEVIATION_COL,
                TARGET_GC_COL, TARGET_LENGTH_COL, ANCHOR_DENSITY_COL]
        for c in cols:
            assert c not in used, f"Column {c} used twice"
            used.add(c)
        for i in range(GC_PROFILE_SLICE.start, GC_PROFILE_SLICE.stop):
            assert i not in used
            used.add(i)
        for i in range(FRAG_PROFILE_SLICE.start, FRAG_PROFILE_SLICE.stop):
            assert i not in used
            used.add(i)
        assert len(used) == N_FEATURES


class TestSampleFingerprint:
    def test_construct_valid(self):
        s = _make_sample()
        assert s.n_targets == N_TARGETS
        assert s.run_context.run_id == "run_001"

    def test_wrong_features_shape(self):
        with pytest.raises(ValueError, match="features shape"):
            SampleFingerprint(
                sample_id="bad",
                panel_id="panel_v1",
                target_names=[f"target_{i}" for i in range(N_TARGETS)],
                features=np.zeros((N_TARGETS, 74), dtype=np.float32),  # wrong
                global_gc_histogram=np.zeros(GLOBAL_GC_BINS, dtype=np.float32),
                global_frag_histogram=np.zeros(GLOBAL_FRAG_BINS, dtype=np.float32),
                total_reads=100,
                run_context=_make_run_context(),
            )

    def test_property_accessors(self):
        s = _make_sample()
        assert s.kmer_counts.shape == (N_TARGETS,)
        assert s.kmer_ratios.shape == (N_TARGETS,)
        assert s.kmer_zscores.shape == (N_TARGETS,)
        assert s.gc_profiles.shape == (N_TARGETS, 20)
        assert s.frag_profiles.shape == (N_TARGETS, 50)
        assert s.gc_deviations.shape == (N_TARGETS,)
        assert s.frag_deviations.shape == (N_TARGETS,)
        assert s.target_gc.shape == (N_TARGETS,)
        assert s.target_lengths.shape == (N_TARGETS,)
        assert s.anchor_densities.shape == (N_TARGETS,)


class TestCohortFingerprints:
    def test_from_samples_roundtrip(self):
        ctx1 = _make_run_context()
        ctx2 = RunContext(
            run_id="run_002",
            sequencer_id="NB505678",
            run_date="2026-02-10",
            n_samples_in_run=4,
            median_total_reads=800_000.0,
            median_kmer_counts=np.full(N_TARGETS, 400.0, dtype=np.float32),
            std_kmer_counts=np.full(N_TARGETS, 40.0, dtype=np.float32),
            median_gc_histogram=np.ones(GLOBAL_GC_BINS, dtype=np.float32) / GLOBAL_GC_BINS,
            median_frag_histogram=np.ones(GLOBAL_FRAG_BINS, dtype=np.float32) / GLOBAL_FRAG_BINS,
            median_gc_profiles=np.ones((N_TARGETS, 20), dtype=np.float32) / 20,
            median_frag_profiles=np.ones((N_TARGETS, 50), dtype=np.float32) / 50,
        )
        samples = [
            _make_sample("s1", run_context=ctx1),
            _make_sample("s2", run_context=ctx1),
            _make_sample("s3", run_context=ctx2),
        ]
        cohort = CohortFingerprints.from_samples(samples)
        assert cohort.n_samples == 3
        assert cohort.n_targets == N_TARGETS
        assert cohort.run_ids == ["run_001", "run_001", "run_002"]
        assert cohort.unique_run_ids == ["run_001", "run_002"]
        assert set(cohort.run_contexts.keys()) == {"run_001", "run_002"}

        # Roundtrip through get_sample
        s0 = cohort.get_sample(0)
        assert s0.sample_id == "s1"
        assert s0.run_context.run_id == "run_001"
        s2 = cohort.get_sample(2)
        assert s2.run_context.run_id == "run_002"

    def test_run_ids_length_validation(self):
        with pytest.raises(ValueError, match="run_ids length"):
            CohortFingerprints(
                panel_id="panel_v1",
                target_names=[f"target_{i}" for i in range(N_TARGETS)],
                sample_ids=["s1", "s2"],
                run_ids=["run_001"],  # wrong length
                run_contexts={"run_001": _make_run_context()},
                features=np.zeros((2, N_TARGETS, N_FEATURES), dtype=np.float32),
                global_gc_histograms=np.zeros((2, GLOBAL_GC_BINS), dtype=np.float32),
                global_frag_histograms=np.zeros((2, GLOBAL_FRAG_BINS), dtype=np.float32),
                total_reads=np.zeros(2, dtype=np.int64),
            )


class TestBuildRunContext:
    def test_basic(self):
        rng = np.random.default_rng(42)
        n_samples, n_targets = 8, N_TARGETS
        kmer = rng.poisson(500, (n_samples, n_targets)).astype(np.float32)
        reads = rng.poisson(1_000_000, n_samples).astype(np.int64)
        gc_hist = rng.dirichlet(np.ones(GLOBAL_GC_BINS), n_samples).astype(np.float32)
        frag_hist = rng.dirichlet(np.ones(GLOBAL_FRAG_BINS), n_samples).astype(np.float32)
        gc_prof = rng.dirichlet(np.ones(20), (n_samples, n_targets)).astype(np.float32)
        frag_prof = rng.dirichlet(np.ones(50), (n_samples, n_targets)).astype(np.float32)

        ctx = build_run_context(
            run_id="run_test",
            sequencer_id="NB999",
            run_date="2026-03-01",
            sample_kmer_counts=kmer,
            sample_total_reads=reads,
            sample_gc_histograms=gc_hist,
            sample_frag_histograms=frag_hist,
            sample_gc_profiles=gc_prof,
            sample_frag_profiles=frag_prof,
        )
        assert ctx.n_samples_in_run == 8
        assert ctx.median_kmer_counts.shape == (n_targets,)
        assert ctx.std_kmer_counts.shape == (n_targets,)
        assert ctx.median_gc_profiles.shape == (n_targets, 20)
        assert ctx.median_frag_profiles.shape == (n_targets, 50)


class TestComputeBatchFeatures:
    def test_ratio_for_diploid(self):
        """A sample with counts equal to run median should have ratio ~1.0."""
        ctx = _make_run_context()
        kmer = ctx.median_kmer_counts.copy()
        gc = ctx.median_gc_profiles.copy()
        frag = ctx.median_frag_profiles.copy()

        ratios, zscores, gc_dev, frag_dev = compute_batch_features(kmer, gc, frag, ctx)
        np.testing.assert_allclose(ratios, 1.0, atol=1e-4)
        np.testing.assert_allclose(zscores, 0.0, atol=1e-4)
        np.testing.assert_allclose(gc_dev, 0.0, atol=1e-6)
        np.testing.assert_allclose(frag_dev, 0.0, atol=1e-6)

    def test_deletion_signal(self):
        """A heterozygous deletion (half the counts) should give ratio ~0.5."""
        ctx = _make_run_context()
        kmer = ctx.median_kmer_counts.copy()
        kmer[2] = kmer[2] * 0.5  # simulate het-del at target 2

        ratios, zscores, _, _ = compute_batch_features(
            kmer, ctx.median_gc_profiles.copy(), ctx.median_frag_profiles.copy(), ctx
        )
        np.testing.assert_allclose(ratios[2], 0.5, atol=1e-3)
        assert zscores[2] < -1.0  # should be significantly negative

    def test_duplication_signal(self):
        """A duplication (1.5x counts) should give ratio ~1.5."""
        ctx = _make_run_context()
        kmer = ctx.median_kmer_counts.copy()
        kmer[0] = kmer[0] * 1.5

        ratios, _, _, _ = compute_batch_features(
            kmer, ctx.median_gc_profiles.copy(), ctx.median_frag_profiles.copy(), ctx
        )
        np.testing.assert_allclose(ratios[0], 1.5, atol=1e-3)
