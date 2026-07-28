# Alignment-based fragment tensor over alignment-free fingerprinting

**Status:** accepted (2026-07-27, supersedes the original alignment-free design)

The project began as *alignment-free* panel CNV detection: characterise a
sample's whole read population as a fixed-dimension k-mer "fingerprint" and
detect CNVs as perturbations of it. We are abandoning that framing. For panels
the alignment step is not the bottleneck (minutes, not the WGS-scale
30+ minutes the speed argument assumed), so "alignment-free" bought little,
while it discarded exactly the machinery that makes panel CNV calling work
(per-coordinate GC/PCA/tangent normalisation against a panel of normals), and
the global-fingerprint signal for a single-target CNV was buried under
technical noise. We are pivoting to an **alignment-based** design: align with
BWA-MEM, build a per-tile **fragment tensor** (position × insert size, plus
off-target-mate channels), and detect CNVs from the **residual** after
**tangent normalisation** against a panel of normals.

## Considered options

- **Anchor-guided placement** — reuse the built `AnchorIndex` to place
  fragments to tiles without a full aligner. Kept the "no aligner" identity and
  most existing code, but placement is sparse and the speed payoff is small for
  panels.
- **Full alignment + rich features (chosen)** — accept BWA-MEM and spend all
  novelty on the fragment-tensor representation and its normalisation. Simplest
  and most robust; fastest route to testing the core idea.
- **Stay strictly alignment-free** — keep the fingerprint direction and its
  difficulty. Rejected: weak single-target signal and no source of training
  labels.

## Consequences

- The alignment-free implementation (k-mer/anchor selection, spaced seeds,
  `AnchorIndex`, FASTQ extraction, the 78-feature fingerprint matrix) is moved
  to `archive/alignment_free/` as a coherent, non-runnable snapshot rather than
  deleted — the "relaxed placement to sidestep alignment bias" idea is a
  deliberate future tangent, not a dead end.
- New runtime dependencies: **BWA-MEM** (external) and **pysam** (BAM reading).
  The tangent-normalisation maths reuses scikit-learn (already a dependency).
- Tangent normalisation is **unsupervised** — it needs only normal samples, not
  labelled CNVs. This removes the "no source of training labels" problem that
  had no answer under the old design and lets the core thesis be validated on
  synthetic data before any real cohort or supervised model exists.
- The `TargetRegion` / `PanelDefinition` types and the
  `RunContext` → panel-of-normals → cohort-stacking *pattern* carry forward
  conceptually, but are re-introduced fresh (trimmed of k-mer/seed fields) in
  the new slices rather than kept in place.
