# Fingerprint Storage Formats for ML Training

> Comparing storage and tensor layout options for the alignment-free sample
> fingerprint described in alignment_free_fingerprinting.md.

---

## Fingerprint Dimensions Recap

For a 200-target panel, the full fingerprint per sample comprises:

| Layer | Content | Dimensions |
|---|---|---|
| 1 | Per-target median anchor k-mer counts | 200 |
| 2 | Per-read GC histogram | 100 bins |
| 3 | Fragment size histogram | 450 bins |
| 4 | Per-target GC profiles | 200 x 20 |
| 5 | Per-target fragment size profiles | 200 x 50 |
| **Total** | | **~14,750** |

---

## File Format Options

### 1. HDF5

Hierarchical, typed, memory-mappable. Native support in PyTorch/TensorFlow via
`h5py` or `torch.utils.data.Dataset`.

```
cohort.h5
├── /metadata
│   ├── sample_ids         (string[N_samples])
│   ├── panel_id           (string)
│   └── fingerprint_version (string)
├── /targets
│   ├── names              (string[N_targets])
│   ├── gc_content         (float32[N_targets])
│   └── anchor_count       (uint16[N_targets])
├── /fingerprint
│   ├── kmer_counts        (float32[N_samples, N_targets])
│   ├── gc_histogram       (float32[N_samples, 100])
│   ├── fragment_histogram (float32[N_samples, 450])
│   ├── target_gc_profiles (float32[N_samples, N_targets, 20])
│   └── target_frag_profiles (float32[N_samples, N_targets, 50])
└── /labels
    ├── cn_states          (int8[N_samples, N_targets])
    └── sample_class       (string[N_samples])
```

| Strengths | Weaknesses |
|---|---|
| Random access by layer and sample index | Single-writer limitation (no concurrent writes) |
| Memory-mapped reads for large cohorts | File corruption risk if process crashes mid-write |
| Hierarchical organisation keeps metadata with data | Slightly more complex API than NumPy |
| Chunked compression (gzip/lz4) built in | Not human-readable |
| Widely supported in scientific Python | |

**Best for**: Production ML training pipelines, large cohorts (>1K samples).

### 2. Apache Parquet

Row-per-sample, columnar compression. Native in pandas, polars, Spark.

```
schema:
  sample_id:            string
  panel_id:             string
  total_reads:          uint64
  kmer_counts:          list<float32>     # [N_targets]
  gc_histogram:         list<float32>     # [100]
  fragment_histogram:   list<float32>     # [450]
  target_gc_profiles:   list<float32>     # [N_targets * 20], reshape on load
  target_frag_profiles: list<float32>     # [N_targets * 50], reshape on load
  cn_labels:            list<int8>        # optional
```

| Strengths | Weaknesses |
|---|---|
| Columnar layout: load only the layers you need | Nested arrays require flatten/reshape on load |
| Excellent compression ratios | No native multi-dimensional array support |
| Integrates with dataframe ecosystems (pandas, polars, Spark) | Less natural for tensor-oriented workflows |
| Append-friendly: easy to add samples | Higher per-row overhead than dense binary |
| Human-inspectable schema | |

**Best for**: Cohort-level analytics, feature selection, exploratory work.

### 3. NumPy `.npz`

Compressed archive of named arrays. Zero-dependency in Python.

```python
# Per-sample
np.savez_compressed("sample.npz",
    kmer_counts=kmer_counts,            # (N_targets,)
    gc_histogram=gc_histogram,          # (100,)
    fragment_histogram=frag_hist,       # (450,)
    target_gc_profiles=tgt_gc,          # (N_targets, 20)
    target_frag_profiles=tgt_frag,      # (N_targets, 50)
)

# Pre-stacked cohort
np.savez_compressed("cohort.npz",
    kmer_counts=all_kmer,               # (N_samples, N_targets)
    labels=all_labels,                  # (N_samples, N_targets)
)
```

| Strengths | Weaknesses |
|---|---|
| Simplest possible API (`np.load`) | Loads entire array into memory (no memory-mapping for compressed) |
| No dependencies beyond NumPy | No metadata/schema support |
| Fast for small-to-medium datasets | Doesn't scale well beyond ~10K samples |
| Easy to inspect and debug | No random access to individual samples |

**Best for**: Prototyping, small cohorts (<1K samples), quick experiments.

### 4. Flat Binary Vector + Manifest

Concatenate all layers into a single float32 vector per sample. Store a JSON
manifest mapping index ranges to layers. Store vectors as a memory-mapped raw
binary matrix.

```json
{
  "fingerprint_version": "0.1.0",
  "panel_id": "panel_v3",
  "n_targets": 200,
  "vector_length": 14750,
  "layout": [
    {"layer": "kmer_counts",          "start": 0,     "end": 200,   "shape": [200]},
    {"layer": "gc_histogram",         "start": 200,   "end": 300,   "shape": [100]},
    {"layer": "fragment_histogram",   "start": 300,   "end": 750,   "shape": [450]},
    {"layer": "target_gc_profiles",   "start": 750,   "end": 4750,  "shape": [200, 20]},
    {"layer": "target_frag_profiles", "start": 4750,  "end": 14750, "shape": [200, 50]}
  ]
}
```

```
cohort.bin  →  raw float32 matrix (N_samples x 14750), memory-mapped
```

| Strengths | Weaknesses |
|---|---|
| Maximum I/O throughput (mmap, no decompression) | No self-describing metadata in the data file |
| Zero-copy random access to any sample | Manifest must be kept in sync with data |
| Minimal storage overhead | Not human-inspectable |
| Trivial `DataLoader` integration | Fragile to schema changes (layout is positional) |
| Scales to millions of samples | Requires custom reader code |

**Best for**: Maximum throughput training at scale, when the schema is stable.

---

## Tensor Layout Options

Orthogonal to file format, the fingerprint can be arranged into different tensor
shapes for model consumption.

### 2D Matrix: Target x Feature

Each row is a target (in genomic order), each column is a feature.

```
Shape: (N_targets, N_features_per_target)
     = (200, 74)

where N_features_per_target = 1 (kmer count) + 20 (GC bins) + 50 (frag bins) + K (scalars)

         kmer   GC profile (20)        Frag profile (50)       Scalars
Target 0 [ 487 | 0.01 0.03 ... 0.02 | 0.00 0.01 ... 0.03  | gc=0.45 len=200 ]
Target 1 [ 251 | 0.02 0.05 ... 0.01 | 0.01 0.02 ... 0.02  | gc=0.62 len=180 ]
...
```

| Strengths | Weaknesses |
|---|---|
| Simple, no padding waste | Loses distinction between feature types |
| Genomic row ordering enables 1D CNN segmentation | Global distributions (Layers 2-3) must be stored separately or broadcast |
| Works with any architecture (CNN, transformer, MLP) | All features share the same scale (needs careful normalisation) |
| Transformer can treat each target as a token (74-dim) | |

**Best for**: General-purpose starting point. Most flexible layout.

### 3D Tensor: Target x Channel x Bin

Separate feature types into channels, pad to uniform bin width.

```
Shape: (N_targets, N_channels, max_bins)
     = (200, 4, 50)

Channel 0: K-mer count (1 value, broadcast to 50 or zero-padded)
Channel 1: GC profile (20 bins, zero-padded to 50)
Channel 2: Fragment size profile (50 bins)
Channel 3: Target metadata (gc, length, anchor_density, ..., padded)
```

Treat as a small image: height = targets (genomic axis), width = bins, depth =
channels. A deletion dims an entire row across all channels. A GC batch effect
shifts channel 1 across all rows.

| Strengths | Weaknesses |
|---|---|
| Preserves feature-type structure for 2D CNNs | Padding wastes space for unequal bin counts |
| Channel separation allows per-feature-type processing | More complex data loading code |
| Natural for architectures with per-channel convolutions | Harder to incorporate scalar features |
| Visually interpretable (plot as heatmap per channel) | |

**Best for**: 2D CNN architectures, when feature-type separation matters.

### 3D Tensor: Chromosome x Position x Feature

Organise by chromosome, with targets placed at intra-chromosomal positions.

```
Shape: (N_chromosomes, max_targets_per_chrom, N_features)
     = (24, ~20, 74)
```

| Strengths | Weaknesses |
|---|---|
| Encodes chromosomal structure explicitly | Very sparse for most panels (many empty slots) |
| Enables chromosome-level pattern learning (whole-arm events) | Variable targets per chromosome requires padding |
| Natural for models with per-chromosome branches | Adds complexity with little gain unless panel is dense |

**Best for**: Dense panels with known chromosome-level biology. Not recommended
for typical gene panels.

---

## Summary Comparison

```
                    File Formats
                    ────────────
Format          Scalability   Simplicity   ML Integration   Metadata
─────────────────────────────────────────────────────────────────────
HDF5            Excellent     Moderate     Excellent        Excellent
Parquet         Excellent     Good         Moderate         Good
NumPy .npz      Poor-Mod      Excellent    Good             Poor
Flat binary     Excellent     Poor         Excellent        Poor

                    Tensor Layouts
                    ──────────────
Layout              Shape              Architectures       Padding
─────────────────────────────────────────────────────────────────────
2D (target x feat)  (200, 74)          CNN/Transformer/MLP None
3D (target x ch x   (200, 4, 50)      2D CNN              Moderate
    bin)
3D (chrom x pos x   (24, ~20, 74)     Chrom-aware models  Heavy
    feat)
```

---

## Recommendations

1. **Start with HDF5 + 2D matrix layout.** This combination offers the best
   balance of simplicity, performance, and flexibility. The 2D layout works
   with all model architectures and the HDF5 container handles metadata,
   compression, and memory-mapped access.

2. **Store raw (unnormalised) values.** Keep raw k-mer counts and raw
   histograms in storage. Apply normalisation (depth correction, GC-LOESS,
   sum-to-1 for distributions) at training time. This preserves the option to
   experiment with normalisation strategies without regenerating fingerprints.

3. **Graduate to flat binary + manifest** if training throughput becomes a
   bottleneck at scale, once the schema has stabilised.

4. **Use Parquet as a secondary export** for cohort-level analytics, QC
   dashboards, and feature importance analysis where dataframe tools are more
   natural than tensor tools.

---

*Companion to: alignment_free_fingerprinting.md (fingerprint definition),
cnv_calling_mind_map.md (software landscape context).*
