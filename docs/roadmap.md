# shard — roadmap

## Goal

Call copy-number variants in targeted gene panels more sensitively than
scalar-depth callers by using **fragment geometry** — the joint pattern of
read start position and insert size at each locus — instead of a single depth
number per target. The pattern is normalised against a panel of normals with
**tangent normalisation**, and CNVs are called from the residual.

The bet: whole-gene deletions are already easy for depth-only tools; the extra
signal in fragment geometry (insert-size shifts, off-target-mate fraction,
sub-exon coverage shape) should help on the *hard* cases — partial/sub-exon
events, low-coverage or noisy targets, and duplications — where scalar depth is
weakest.

## How it works (target pipeline)

```
paired-end FASTQ
   │  BWA-MEM
   ▼
aligned BAM
   │  tiling: targets + flanks → fixed-width tiles
   ▼
per-tile fragment records (start, insert size, orientation, off-target mate)
   │  histogram per tile
   ▼
fragment tensor per sample   (n_tiles × position × insert-size, + channels)
   │  panel of normals → SVD → tangent normalisation
   ▼
residual per tile   (shared technical structure removed)
   │  segmentation (HMM / changepoint) over tiles in genomic order
   ▼
copy-number segments → CNV calls
```

Everything upstream of the residual is deterministic feature extraction;
tangent normalisation is unsupervised (needs only normals); segmentation is the
only place a learned/tuned model is required.

## Status

| Area | Module | State |
|---|---|---|
| Panel / tiling | `data/panel.py`, tiling | **to build** (Slice 1; lift trimmed `TargetRegion`/`PanelDefinition` from archive) |
| Fragment extraction | `extraction` from BAM | **to build** (Slice 2; needs pysam) |
| Fragment tensor | `data/` tensor builder | **to build** (Slice 3) |
| Storage | `utils/io.py` | stub → Slice 4 |
| Normalisation | `data/normalisation.py` | stub → **Slice 5 (core thesis)** |
| Segmentation / calling | `models/`, `inference/` | stub → Slice 6 |
| Training (optional) | `training/*` | stub → only if a supervised head is added later |
| Scripts | `scripts/*` | stub → Slice 7 |
| Alignment-free approach | `archive/alignment_free/` | archived snapshot (see ADR 0001) |

## Slices

Each slice is one independently-testable feature. Slices 1–6 are testable with
synthetic / hand-built inputs — **no real cohort required**. Recommended order:
**0 → 1 → 2 → 3 → 5 → 4/6**, so the make-or-break normalisation test (Slice 5)
is reached as early as possible.

- **Slice 0 — Archive & reset.** Move the alignment-free modules + their tests
  to `archive/alignment_free/`. *Done as part of the pivot.*
- **Slice 1 — Tiling.** `(targets BED, tile_width, flank) → ordered tiles`.
  Pure function. *Test:* toy panel tiles deterministically; edge cases (target
  shorter than a tile, flank off chromosome end, adjacent targets).
  Design decisions and activity breakdown: [slice-1-tiling.md](slice-1-tiling.md).
- **Slice 2 — Fragment extraction from BAM.** `(BAM, tiles) → per-tile fragment
  records`. Uses pysam. *Test:* hand-written SAM → known fragments → known tile
  assignment.
- **Slice 3 — Fragment tensor.** `(fragments in a tile) → (position × insert
  size) histogram + channels`. Pure function. *Test:* known fragments → known
  tensor cells.
- **Slice 4 — Sample assembly + storage.** Stack tiles into a per-sample
  tensor; HDF5 write/read. *Test:* roundtrip equality.
- **Slice 5 — Panel of normals + tangent normalisation.** `(normal tensors) →
  PoN (SVD)`; `(sample, PoN) → residual`. *Test (fully synthetic, the crux):*
  build normals from a known low-rank subspace + noise, inject a scaled-down
  region, assert tangent normalisation suppresses shared structure and surfaces
  the injected signal in the residual. If this fails, the pivot is refuted
  cheaply.
- **Slice 6 — Segmentation / calling.** `(residual over ordered tiles) → CN
  segments`. *Test:* synthetic residual with a planted deletion → recovered
  boundaries.
- **Slice 7 — End-to-end on real data.** Wire `build_fingerprints` / `predict`
  scripts; validate against a truth set. **Only slice needing a real cohort.**

## Open questions (deferred to Slice 7)

- **Panel scale** — tens of genes vs full exome. Sets tile count and PoN
  dimensions.
- **Truth set** — existing ExomeDepth/DECoN calls, MLPA/ddPCR-confirmed CNVs,
  or simulated spike-ins. Decides how calling accuracy is validated.
- **Supervised head?** — whether Slice 6 stays unsupervised (segment the
  residual) or a learned model is added on top. Depends on whether labelled
  data materialises.
