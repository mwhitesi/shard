# Read-Pair Encoder: Transformer Architecture for K-mer-Level CNV Evidence

> A transformer-based approach that processes raw k-mers from paired-end reads,
> implicitly capturing fragment size, GC bias, and mapping confidence without
> explicit feature engineering.

---

## 1. Motivation

The current fingerprint model pre-aggregates read evidence into fixed per-target
features (median k-mer count, binned GC/fragment distributions). This discards
per-read structure:

- **Fragment size** requires alignment or explicit estimation from adapter
  read-through or anchor co-occurrence.
- **GC bias** is summarised as a histogram, losing the correlation between a
  specific fragment's GC content and its depth contribution.
- **Mapping ambiguity** (paralog-derived reads, off-target fragments) is
  invisible after median aggregation.

A read-pair encoder processes k-mers at the individual read level and lets
attention discover these patterns end-to-end.

---

## 2. Input Representation

### 2.1 From FASTQ to Token Sequence

For each read pair (R1, R2), extract all k-mers (k=31) from both reads and
look up each against the panel anchor set (using spaced seeds for SNP
tolerance).

```
R1 (150bp) → ~120 k-mers → match against anchor set → hits + misses
R2 (150bp) → ~120 k-mers → match against anchor set → hits + misses
```

Each anchor hit produces a token with structured information:

```python
@dataclass
class AnchorHit:
    target_id: int        # which panel target this anchor belongs to
    offset: int           # position of the anchor within the target region
    read_id: int          # 0 = R1, 1 = R2
    read_position: int    # position within the read (0-indexed from 5')
    gc_content: float     # GC fraction of the k-mer itself
```

Non-anchor k-mers (no match in the panel) are dropped. A read pair with no
anchor hits in either read is discarded entirely.

### 2.2 Token Sequence Construction

Tokens from R1 and R2 are concatenated into a single sequence with special
delimiter tokens:

```
[CLS] [R1] hit hit hit ... [R2] hit hit hit ... [SEP]

Example (target 5, offsets shown):
[CLS] [R1] (t5,+12) (t5,+43) (t5,+74) [R2] (t5,+180) (t5,+149) (t5,+118) [SEP]
              ↑ ascending (R1 forward)          ↑ descending (R2 reverse complement)
```

Typical sequence length: 5-30 anchor hits per read pair for a well-designed
panel. Some pairs will have very few hits (off-target reads) or hits
spanning multiple targets (chimeric fragments).

### 2.3 Token Embedding

Each token is embedded as the sum of several component embeddings:

```
token_embedding = target_emb + offset_emb + read_emb + gc_emb + position_emb

target_emb:    learned embedding, dim d_t (~16)
               index: target_id ∈ {0, ..., n_targets-1, UNK_TARGET}

offset_emb:    continuous positional encoding of the anchor's offset within
               the target region. Use sinusoidal encoding or a small MLP
               on the normalised offset (offset / target_length).
               dim d_o (~16)

read_emb:      learned embedding, dim d_r (~4)
               index: {R1=0, R2=1, CLS=2, R1_DELIM=3, R2_DELIM=4, SEP=5}

gc_emb:        linear projection of k-mer GC fraction (scalar → dim d_g (~8))

position_emb:  position within the token sequence (standard transformer
               positional encoding), dim d_model
               This encodes ordering of hits within a read.
```

Total embedding dimension: d_model = d_t + d_o + d_r + d_g (or use a
projection layer to map the concatenation to d_model).

**Recommended d_model**: 64. This is a very short sequence (<32 tokens
typically) with simple structure — a small model suffices.

---

## 3. Encoder Architecture

### 3.1 Read-Pair Encoder (Level 1)

A small transformer encoder processes one read pair at a time.

```
ReadPairEncoder:
    input:  token sequence, shape (seq_len, d_model)
    layers: 2 transformer encoder layers
            - d_model = 64
            - n_heads = 4 (head_dim = 16)
            - d_ff = 128
            - dropout = 0.1
    output: [CLS] token embedding → pair_embedding, shape (d_model,)
```

**What attention learns at this level:**

- **Fragment size proxy**: Attention between the last R1 anchor and first R2
  anchor captures the offset gap — directly proportional to insert size.
  No explicit TLEN needed.

- **Read coherence**: Pairs where R1 and R2 anchors hit the same target with
  consistent offset ordering (ascending in R1, descending in R2) produce a
  coherent attention pattern. Chimeric or mismapped pairs produce scattered
  patterns that attention can learn to downweight.

- **Per-fragment GC**: The GC embeddings of all k-mers in a pair collectively
  represent the fragment's GC content. Attention can learn the relationship
  between fragment GC and depth contribution (the core of GC bias).

### 3.2 Per-Target Aggregation (Level 2)

All read pairs are assigned to targets based on their anchor hits. A pair
that hits only target T is assigned to T. Pairs hitting multiple targets
are assigned to each hit target (counted fractionally or duplicated).

For each target, pool the pair embeddings:

```
TargetAggregator:
    input:  all pair embeddings assigned to target t
            shape (n_pairs_t, d_model)
    method: one of:
            (a) mean pooling (simplest)
            (b) attention-weighted pooling (learned importance per pair)
            (c) small transformer (2 layers) if n_pairs_t is bounded

    output: target_embedding, shape (d_model,)
```

Option (b) is recommended — it lets the model learn to downweight noisy
pairs (few anchor hits, multi-target, high GC deviation from target mean)
without explicit filtering rules.

```python
# Attention-weighted pooling
class AttentionPool(nn.Module):
    def __init__(self, d_model: int):
        super().__init__()
        self.attn = nn.Linear(d_model, 1)

    def forward(self, pair_embeddings: Tensor, mask: Tensor) -> Tensor:
        # pair_embeddings: (n_pairs, d_model)
        # mask: (n_pairs,) boolean, True for valid pairs
        scores = self.attn(pair_embeddings).squeeze(-1)  # (n_pairs,)
        scores = scores.masked_fill(~mask, float('-inf'))
        weights = F.softmax(scores, dim=0)               # (n_pairs,)
        return (weights.unsqueeze(-1) * pair_embeddings).sum(dim=0)
```

After aggregation, we have one embedding per target:

```
target_representations: shape (n_targets, d_model)
```

### 3.3 Auxiliary Count Features

The target embeddings capture nuanced per-read evidence, but the simple
count ratio is still the strongest single CNV feature. Inject the current
fingerprint features as a parallel branch:

```
fingerprint_features: shape (n_targets, 78)    # current feature matrix
target_embeddings:    shape (n_targets, d_model) # from read-pair encoder

fused: shape (n_targets, 78 + d_model)          # concatenated
     → linear projection → (n_targets, d_fused)
```

This lets the model use the reliable count-ratio signal as a backbone and
learn corrections from the read-level evidence.

### 3.4 CNV Segmentation Head (Level 3)

Operates on the fused per-target representations in genomic order.

```
CNVSegmentationHead:
    input:  (n_targets, d_fused)  — targets ordered by genomic position
    architecture: one of:
        (a) 1D CNN (causal convolutions, kernel_size=5-11)
        (b) Transformer encoder (targets as tokens, positional = genomic order)
        (c) Bidirectional LSTM

    output: (n_targets, n_cn_states)  — logits per target per CN state
            n_cn_states = 5 (CN0, CN1, CN2, CN3, CN4+)
```

The segmentation head captures spatial dependencies — a deletion spanning
3 consecutive targets should be called as a single event, not 3 independent
calls. The 1D CNN is simplest and sufficient for typical panel sizes
(~200 targets).

---

## 4. Full Architecture Diagram

```
FASTQ Read Pairs
       │
       ▼
  K-mer extraction + anchor lookup (spaced seeds)
       │
       ▼
  Per-pair token sequences: [CLS] [R1] hits... [R2] hits... [SEP]
       │
       ▼
┌──────────────────────┐
│  Read-Pair Encoder   │  Level 1: pair → embedding
│  (2-layer transformer│  Processes each read pair independently
│   d_model=64)        │  Output: [CLS] embedding per pair
└──────┬───────────────┘
       │
       ▼
  Group pair embeddings by target assignment
       │
       ▼
┌──────────────────────┐
│  Target Aggregator   │  Level 2: pairs → target embedding
│  (attention-weighted │  Pools all pair embeddings per target
│   pooling)           │  Output: one vector per target
└──────┬───────────────┘
       │
       ├──── concatenate ◄── Fingerprint features (n_targets, 78)
       │                     (count ratios, batch features, metadata)
       ▼
┌──────────────────────┐
│  Fusion + Projection │  Linear: (78 + d_model) → d_fused
└──────┬───────────────┘
       │
       ▼
┌──────────────────────┐
│  CNV Segmentation    │  Level 3: target sequence → CN states
│  (1D CNN or          │  Captures multi-target events (deletions,
│   transformer)       │  duplications spanning consecutive targets)
└──────┬───────────────┘
       │
       ▼
  Per-target CN state predictions: (n_targets, 5)
```

---

## 5. Training Considerations

### 5.1 Data Pipeline

The read-pair encoder requires access to raw reads during training, not just
pre-computed fingerprints. Two options:

**Option A: Online processing**
Stream FASTQ pairs, extract anchor hits on-the-fly, feed to model. Memory
efficient but I/O-bound and slow for training.

**Option B: Pre-extracted anchor hits (recommended)**
At fingerprinting time, save the per-read anchor hit records to a compact
format alongside the aggregated fingerprint. This is a one-time cost and
keeps training fast.

```
sample_hits.h5
├── /pairs
│   ├── target_ids       int16[n_pairs, max_hits]
│   ├── offsets          int16[n_pairs, max_hits]
│   ├── read_ids         int8[n_pairs, max_hits]
│   ├── read_positions   int8[n_pairs, max_hits]
│   ├── gc_contents      float16[n_pairs, max_hits]
│   ├── n_hits           int8[n_pairs]
│   └── pair_targets     int16[n_pairs]  (primary target assignment)
└── /metadata
    ├── sample_id        string
    └── total_pairs      int64
```

Estimated storage: ~20 bytes per pair × 200K pairs/sample ≈ 4 MB per sample
(compact, much smaller than the FASTQ).

### 5.2 Subsampling Strategy

Not all read pairs are equally informative. For training efficiency:

- **Cap at N pairs per target** (e.g., 200-500). CNV calling needs depth
  ratios, not every individual read. Subsampling also acts as regularisation.
- **Stratified sampling**: Ensure each target gets a representative share of
  pairs. Targets with very high coverage are downsampled; targets with low
  coverage keep all pairs.
- **Hard example mining**: After initial training, oversample pairs from
  targets where the model makes errors (high-loss targets). These are likely
  the hard targets with paralog contamination or low anchor density.

With 200 targets × 300 pairs/target = 60K pairs per sample, and a max
sequence length of ~30 tokens per pair, the Level 1 encoder processes
60K × 30 ≈ 1.8M tokens per sample. At d_model=64, this is feasible on a
single GPU.

### 5.3 Loss Function

```
loss = cross_entropy(cn_predictions, cn_labels)
     + λ_seg * segmentation_smoothness_penalty
     + λ_count * count_ratio_reconstruction_loss  (auxiliary)
```

- **Primary**: Per-target cross-entropy over CN states (0-4+).
- **Segmentation smoothness**: Penalise adjacent targets with different CN
  predictions unless there's strong evidence for a breakpoint. Encourages
  contiguous segments.
- **Auxiliary count reconstruction**: Predict the log count ratio from the
  pair encoder output alone (without the fingerprint features). This
  encourages the pair encoder to learn depth-related representations and
  prevents it from ignoring the read-level signal in favour of the
  already-strong fingerprint features.

### 5.4 Curriculum

1. **Phase 1**: Train the segmentation head on fingerprint features only
   (no read-pair encoder). Establishes a strong count-ratio baseline.
2. **Phase 2**: Freeze the segmentation head, train the read-pair encoder
   + aggregator with the auxiliary count reconstruction loss. Forces the
   encoder to learn useful read-level representations.
3. **Phase 3**: Unfreeze everything, fine-tune end-to-end with reduced
   learning rate. The read-pair encoder learns to correct the fingerprint
   baseline on hard targets.

---

## 6. What the Model Implicitly Learns

| Explicit feature (current) | Implicit equivalent (learned by attention) |
|---|---|
| Fragment size histogram | Offset gap between R1 and R2 anchors |
| Per-target GC profile | GC content of k-mers grouped by target |
| Median k-mer count | Number of high-confidence pair embeddings per target |
| Batch normalisation (RunContext ratio) | Still injected via fingerprint features |
| Mapping quality filtering | Attention weights downweight incoherent pairs |

The key advantage is that these aren't independent pre-computed features —
the model sees their *joint* distribution per read pair and can learn
interactions (e.g., "short fragments with high GC are overrepresented at
this target, consistent with GC-biased PCR duplication, not a real
duplication event").

---

## 7. Complexity and Scaling

| Component | Time per sample | Memory |
|---|---|---|
| Anchor hit extraction | O(n_reads × k) | Streaming, ~constant |
| Level 1 (pair encoder) | O(n_pairs × seq_len² × d_model) | O(batch × seq_len × d_model) |
| Level 2 (aggregation) | O(n_pairs × d_model) | O(n_targets × d_model) |
| Level 3 (segmentation) | O(n_targets × d_fused) | O(n_targets × d_fused) |

With n_pairs=60K, seq_len=30, d_model=64: Level 1 dominates at ~3.5G FLOPs
per sample. This is comparable to a single forward pass of a small ResNet
on a 224×224 image — feasible for training on a single GPU.

**Inference optimisation**: The Level 1 encoder processes pairs independently,
so it parallelises trivially across pairs within a sample. On a modern GPU
with batched attention, 60K short sequences is fast.

---

## 8. Incremental Adoption Path

The architecture is designed for incremental adoption:

1. **Current state**: Fingerprint features only → segmentation head.
   No read-pair encoder. This is the baseline.

2. **Add pre-extracted hits**: Modify fingerprinting pipeline to save anchor
   hit records alongside aggregated features. No model change yet, just
   data pipeline preparation.

3. **Train pair encoder offline**: Use the auxiliary count reconstruction
   loss to pre-train the pair encoder. Evaluate whether the learned target
   embeddings improve over the fixed features on held-out samples.

4. **Fuse and fine-tune**: Concatenate pair encoder output with fingerprint
   features, fine-tune the full model. Compare against fingerprint-only
   baseline.

5. **Ablation**: Measure contribution of pair encoder on easy targets (high
   anchor density, no known paralogs) vs. hard targets. If the gain is
   concentrated on hard targets, consider using the pair encoder selectively.

---

*Companion to: alignment_free_fingerprinting.md (fingerprint definition),
fingerprint_storage_formats.md (storage), kmer_encoding_schemes_for_ml.md
(encoding background).*
