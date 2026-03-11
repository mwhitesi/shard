# Alignment-Free Sample Fingerprinting for Panel CNV Detection

> Can we characterise the total read population of a gene panel sample —
> without alignment — in a way that creates a unique, informative fingerprint
> sensitive to copy number perturbations?

---

## The Core Idea

Alignment maps each read to a coordinate, then counts reads per region.
An alignment-free approach inverts this: instead of asking "where does this
read come from?", ask "what is the composition of the entire read population,
and how does that composition shift when copy number changes?"

A gene panel sample is a *mixture* of reads drawn from a finite set of target
regions. Each target contributes a sub-population of reads with characteristic
properties:

```
Gene Panel Sample (FASTQ)
│
├── Reads from Target A  ──→  characteristic k-mer set, GC%, fragment sizes
├── Reads from Target B  ──→  characteristic k-mer set, GC%, fragment sizes
├── Reads from Target C  ──→  ...
├── Off-target reads     ──→  background k-mer noise
└── Adapter/artifact     ──→  filterable
```

If Target B is deleted (CN=1 instead of CN=2), the read population loses ~half
of Target B's reads. This shifts:
- The abundance of B-specific k-mers
- The overall GC distribution (by the GC signature of B's reads)
- The fragment size distribution (by B's fragment size signature)
- The relative proportions of all other targets (they become slightly enriched)

The fingerprint is the multi-dimensional descriptor of the total read
population. CNVs perturb specific, predictable components of this fingerprint.

---

## 1. K-mer Target Signatures

### Anchor K-mers

The most direct alignment-free proxy for per-target read depth is **anchor
k-mer counting**. For each panel target, identify k-mers that are:

1. **Unique to that target** — appear nowhere else in the genome
2. **Interior to the target** — not at edges where capture drops off
3. **Mappability-equivalent to unique alignment** — a k-mer that appears
   exactly once in the reference is equivalent to a read mapping uniquely

```
Target Region Design
│
├── Reference sequence for target (e.g., exon of BRCA1)
│   ACGTACGT...ACGTACGT  (200bp captured region)
│
├── Enumerate all k-mers (k=31 is standard for uniqueness in human genome)
│   Position 1: ACGTACGTACGTACGTACGTACGTACGTACG
│   Position 2: CGTACGTACGTACGTACGTACGTACGTACGT
│   ...
│
├── Filter to ANCHOR k-mers:
│   ├── Appears exactly once in hg38 (unique)
│   ├── Not in first/last 50bp of target (edge effects)
│   ├── GC content 30-70% (avoid extreme-GC k-mers with biased counting)
│   └── No homopolymer runs >5bp (sequencing error hotspots)
│
└── Result: N anchor k-mers per target (typically 10-100 for a 200bp exon)
```

### K-mer Counting as Depth Proxy

Count anchor k-mers directly in the FASTQ using exact matching (Jellyfish,
KMC, or a custom hash table). The count of a target's anchor k-mers is
directly proportional to the number of reads from that target — i.e., to
depth, and therefore to copy number.

```
Per-Target K-mer Count Vector
│
├── Target A: median anchor k-mer count = 487
├── Target B: median anchor k-mer count = 251   ← half expected → CN=1?
├── Target C: median anchor k-mer count = 502
├── Target D: median anchor k-mer count = 493
│   ...
└── Target N: median anchor k-mer count = 478

→ This vector IS the alignment-free equivalent of a depth profile
```

**Why median?** Individual k-mers within a target will have variable counts
due to position-specific sequencing biases. The median across a target's
anchor k-mers is robust to outliers — analogous to how alignment-based tools
use median depth per bin.

### Advantages Over Alignment-Based Depth

| Property | Alignment-based depth | Anchor k-mer counts |
|---|---|---|
| Multi-mapping ambiguity | Major problem (MAPQ heuristic) | Eliminated by design (anchors are unique) |
| Speed | Minutes to hours (full alignment) | Seconds to minutes (k-mer counting) |
| Reference bias | Yes (aligner penalises variants) | Reduced (can include variant k-mers) |
| Segmental duplications | Depth inflated/ambiguous | Anchors exclude non-unique k-mers |
| Sensitivity to aligner | MAPQ differs between aligners | Deterministic (exact match) |

### Limitations

- **Target regions in repetitive sequence**: If a target has few or no unique
  k-mers (e.g., a pseudogene), anchor k-mer counting cannot assess it. This
  is the same fundamental limitation as mappability, just made explicit.
- **SNPs/variants in anchor k-mers**: A heterozygous SNP in an anchor k-mer
  destroys the match for ~half the reads. Mitigation: use enough anchors per
  target that losing a few to variants doesn't matter, or include known
  common variant k-mers in the anchor set.
- **Sequencing errors**: An error in an anchor k-mer position prevents the
  match. At ~0.1% per-base error rate and k=31, ~3% of k-mer instances are
  corrupted. This is a small, uniform loss — it reduces sensitivity slightly
  but doesn't create bias.

---

## 2. GC Content Distribution as a Population Descriptor

### Per-Read GC Distribution

Every read has a GC fraction. The distribution of per-read GC across the
entire sample is a fingerprint of which genomic regions contributed reads.

```
GC Distribution Histogram (per sample)
│
│  Bin the GC fraction of every read into (e.g.) 100 bins from 0.0 to 1.0
│
│        ╭───╮
│       ╭┤   ├╮
│      ╭┤│   ││╮
│     ╭┤││   │││╮
│    ╭┤│││   ││││╮
│   ╭┤││││   │││││╮
│  ╭┤│││││   ││││││╮
│ ─┴┴┴┴┴┴┴───┴┴┴┴┴┴┴──
│ 0.0              1.0   GC fraction
│
│ Shape reflects the GC composition of all targets weighted by their
│ read counts (which are proportional to copy number × capture efficiency)
```

### How CNVs Perturb the GC Distribution

Each panel target has a characteristic GC content. A CNV affecting a target
shifts the GC distribution by removing or adding reads at that target's GC%.

For a small panel (50-500 targets), individual targets contribute measurable
fractions of the total read population. A deletion of a high-GC target
produces a detectable dip in the high-GC portion of the distribution.

```
Example: 100-target panel, target B has GC=65%, all others ~40-55%

Normal:     GC distribution has a bump at 65% from target B
CN=1 at B:  Bump at 65% shrinks by ~50%
CN=0 at B:  Bump at 65% disappears

→ The GC distribution is a PROJECTION of the copy number profile
  onto the GC axis
```

**Key insight**: The GC distribution is a *lossy projection* — multiple
targets with similar GC content are conflated. It cannot uniquely identify
which target is affected, but it provides a complementary signal to k-mer
counts and can detect global perturbations (like whole-arm events or batch
effects) that shift the overall GC profile.

### Differential GC Analysis

Compare a test sample's GC distribution to a reference distribution
(from a PoN or cohort median). The residual reveals:

- **Localised bumps/dips**: Suggest CNVs at targets with that GC content
- **Broad shifts**: Suggest batch effects or GC bias differences
- **Asymmetric tails**: Suggest contamination or adapter artifacts

This is useful as a **QC metric** and as a **supporting signal** for k-mer-
based calls, not as a primary CNV detector.

---

## 3. Fragment Size Distribution

### Population-Level Fragment Signature

The insert size distribution of the library is shaped by:
- Fragmentation method (sonication, enzymatic, tagmentation)
- Size selection parameters
- Local sequence properties at each target
- For cfDNA: nucleosome protection patterns (biological, not technical)

```
Fragment Size Distribution (from paired-end reads, no alignment needed)
│
│  Infer insert size from:
│  ├── Read overlap in paired-end data (if reads overlap, insert < 2×read_length)
│  ├── Or: k-mer-based insert size estimation (find shared k-mers between R1/R2
│  │       and the target reference to estimate fragment length)
│  └── Or: for SE data, use read-level features only (less informative)
│
│     ╭──────╮
│    ╭┤      ├──╮
│   ╭┤│      │  ├──╮
│  ─┴┴┴──────┴──┴──┴──
│  100   200   300   400   Fragment size (bp)
```

### How CNVs Perturb Fragment Size

Each target produces fragments with a characteristic size distribution
(influenced by local GC, chromatin, and distance to capture probe edges).
A CNV removes or adds a target's contribution to the overall fragment
size profile.

The effect is subtle and highly degenerate (many targets produce similar
fragment sizes), so fragment size alone is a weak CNV signal. Its value is
as a **joint feature** with k-mer counts and GC:

```
For each target's anchor k-mers, you can ask:
  "What is the fragment size distribution of reads containing these k-mers?"

This gives a PER-TARGET fragment size profile without alignment:
  Target A: fragments 180-220bp (tight, GC-moderate target)
  Target B: fragments 150-280bp (broad, GC-extreme target)

A deletion at target A would specifically deplete 180-220bp fragments
that also contain target A's anchor k-mers.
```

### cfDNA-Specific: Fragmentomic Fingerprinting

For cfDNA panels, fragment size is biologically informative. Tumor-derived
cfDNA fragments are shorter (~140bp, sub-nucleosomal) than normal cfDNA
(~167bp, mono-nucleosomal). A target-specific shift toward shorter fragments
suggests tumor-derived reads — which correlates with somatic CNV status.

---

## 4. The Composite Fingerprint

### Multi-Dimensional Sample Descriptor

Combine all three axes into a single structured fingerprint per sample:

```
Sample Fingerprint
│
├── Layer 1: K-mer Count Vector (primary CNV signal)
│   ├── Per-target median anchor k-mer count
│   ├── Dimension: N_targets (e.g., 200 targets → 200-dimensional vector)
│   └── Analogous to: alignment-based depth profile
│
├── Layer 2: GC Distribution (QC + supporting signal)
│   ├── Histogram of per-read GC fractions (100 bins)
│   ├── Dimension: 100
│   └── Analogous to: GC bias curve
│
├── Layer 3: Fragment Size Distribution (supporting signal)
│   ├── Histogram of insert sizes (e.g., 50-500bp in 1bp bins → 450 bins)
│   ├── Dimension: 450
│   └── Analogous to: library QC metrics
│
├── Layer 4: Per-Target K-mer × GC Joint Distribution (cross-signal)
│   ├── For each target: GC distribution of reads containing that target's anchors
│   ├── Detects shifts in which sub-population of reads contributes to a target
│   └── Sensitive to contamination, mismapping equivalents, and allelic bias
│
└── Layer 5: Per-Target K-mer × Fragment Size Joint (cross-signal)
    ├── For each target: fragment size distribution of anchor-containing reads
    ├── Sensitive to cfDNA tumor fraction per target
    └── Detects technical artifacts (e.g., preferential loss of long fragments)
```

### Dimensionality

```
Layer                       Dimensions    Information Content
──────────────────────────────────────────────────────────────
K-mer count vector          N_targets     HIGH (direct CN proxy)
GC distribution             ~100          MEDIUM (global QC)
Fragment size distribution  ~450          LOW-MEDIUM (library QC)
Per-target GC profiles      N_targets×20  MEDIUM (per-target QC)
Per-target frag profiles    N_targets×50  LOW (subtle, noisy)
──────────────────────────────────────────────────────────────
Total                       ~N_targets×71 + 550
For 200 targets:            ~14,750 dimensions
```

Most of the discriminative power for CNV detection lives in Layer 1 (the
k-mer count vector). Layers 2-5 serve as:
- **Bias correction features**: GC and fragment distributions help normalise
  the k-mer counts without explicit GC correction
- **QC signals**: Global distribution shifts flag batch effects before they
  corrupt CNV calls
- **Complementary evidence**: Joint distributions provide per-target context
  that pure counts lack

---

## 5. Detecting CNV Perturbations

### Reference Model

Build a reference fingerprint from a cohort of known-normal samples:

```
Reference Construction
│
├── Process N normal samples (N ≥ 20) through the fingerprinting pipeline
│
├── For Layer 1 (k-mer counts):
│   ├── Compute per-target median and variance across normals
│   ├── GC-normalise: regress k-mer counts against per-target GC
│   │   (the alignment-free equivalent of LOESS GC correction)
│   └── Result: expected count ± variance per target
│
├── For Layers 2-3 (GC + fragment distributions):
│   ├── Compute median histogram across normals
│   └── Result: expected distribution ± per-bin variance
│
└── For Layers 4-5 (joint distributions):
    ├── Per-target expected GC and fragment profiles
    └── Used for per-target bias correction
```

### Test Sample Scoring

```
CNV Detection from Fingerprint
│
├── 1. Compute test sample fingerprint (all layers)
│
├── 2. Normalise k-mer count vector:
│   ├── Total-count normalisation (divide by sum of all k-mer counts)
│   │   → accounts for overall sequencing depth
│   ├── GC normalisation using Layer 2 or per-target GC profiles
│   │   → removes GC bias without alignment
│   └── Reference ratio: log2(test_count / reference_count) per target
│
├── 3. Detect outlier targets:
│   ├── Z-score against reference distribution
│   ├── Or: likelihood ratio under Poisson/NB model
│   ├── Threshold: |z| > 3 for single-target, HMM for multi-target
│   └── Expected log2 ratios:
│       CN=0: -∞ (no reads)
│       CN=1: -1.0
│       CN=2:  0.0 (baseline)
│       CN=3: +0.58
│       CN=4: +1.0
│
├── 4. Segment adjacent targets (if panel is dense enough):
│   ├── CBS or HMM on the ordered k-mer ratio profile
│   └── Panel gene order provides genomic ordering without alignment
│
└── 5. Cross-validate with GC and fragment distribution shifts:
    ├── Does the GC distribution shift match what's expected if
    │   the flagged target(s) changed copy number?
    └── Does the fragment profile shift corroborate?
```

### What Makes This a "Fingerprint" vs Just K-mer Counting

Pure k-mer counting (Layer 1 alone) is effectively alignment-free depth
profiling — a known technique. The fingerprint concept adds value by:

1. **Encoding bias structure alongside signal**: Layers 2-5 capture the
   sample's bias state (GC curve shape, fragment distribution, per-target
   technical profiles). This enables self-normalisation — the bias correction
   features travel with the data, not with an external reference.

2. **Enabling distance-based sample comparison**: Two fingerprints can be
   compared directly (cosine similarity, Euclidean distance, Mahalanobis
   distance) to detect whether samples are from the same batch, same
   individual, same capture kit — or whether something has changed.

3. **Providing a compact, fixed-dimension representation**: Regardless of read
   count, every sample maps to the same-dimensional vector. This makes the
   fingerprint suitable as input to ML models, clustering, and anomaly
   detection without per-sample preprocessing.

---

## 6. GC Normalisation Without Alignment

A key challenge: GC bias must still be corrected, but without alignment we
don't have per-bin genomic GC content to regress against.

### Approach: Target-Intrinsic GC

Each panel target has a known GC content (from the reference sequence used to
design the panel). Each target's anchor k-mers have a known GC. So the GC
correction becomes:

```
Alignment-Free GC Correction
│
├── For each target t with known GC content gc_t:
│   ├── Observed: k-mer count c_t
│   ├── Fit LOESS: c ~ gc across all targets (using the reference set)
│   ├── Expected from LOESS: ĉ_t = f(gc_t)
│   └── Corrected: c_t / ĉ_t  (or log2 ratio)
│
└── This is exactly the standard GC-LOESS correction,
    but using target-level GC instead of bin-level GC,
    and k-mer counts instead of aligned read counts
```

The per-read GC distribution (Layer 2) provides an additional handle:
if the sample's GC distribution is shifted relative to the reference, this
indicates a global GC bias difference that should be corrected before
target-level analysis.

---

## 7. Practical Considerations

### K-mer Size Selection

```
k       Uniqueness in hg38    Error Sensitivity    Memory
─────────────────────────────────────────────────────────────
21      ~85% unique            ~2% corrupted        ~4^21 = 4T (too large
                                                     for naive table)
25      ~95% unique            ~2.5% corrupted      Compressed OK
31      ~99% unique            ~3% corrupted        Standard choice ★
35      ~99.5% unique          ~3.5% corrupted      Marginal gain
─────────────────────────────────────────────────────────────

k=31 is the standard tradeoff: nearly all 31-mers in the human genome are
unique, error sensitivity is manageable, and efficient k-mer counting tools
(Jellyfish, KMC3) handle this natively.
```

### Panel Design Implications

Not all panel targets are equally amenable to anchor k-mer fingerprinting:

```
Target Amenability
│
├── Excellent: unique-sequence exons, moderate GC, no known SVs nearby
│   → Many anchor k-mers, stable counts, clean signal
│
├── Good: unique sequence but extreme GC or short target
│   → Fewer usable anchors, higher variance, but workable
│
├── Poor: targets in segmental duplications or recent paralogs
│   → Few or no unique k-mers; these targets are fundamentally
│     unresolvable without long reads (same as with alignment)
│
└── Informative metric: "anchor density" = unique k-mers / target length
    Targets with anchor density < 0.1 should be flagged as unreliable
```

### Computational Cost

```
Step                        Time (200-target panel, 20M read pairs)
────────────────────────────────────────────────────────────────────
K-mer counting (KMC3)      ~10-30 seconds
Anchor k-mer lookup         ~5-10 seconds (hash table query)
GC/fragment distributions   ~10-20 seconds (single pass over FASTQ)
Normalisation + calling     ~1 second
────────────────────────────────────────────────────────────────────
Total                       ~30-60 seconds per sample

vs alignment-based:
BWA-MEM alignment           ~5-20 minutes
Sort + duplicate marking    ~5-10 minutes
Coverage computation        ~2-5 minutes
────────────────────────────────────────────────────────────────────
Total                       ~15-35 minutes per sample
```

The alignment-free approach is 20-50x faster. For high-throughput clinical
panels where turnaround time matters, this is significant.

---

## 8. Relationship to Existing Work

### Precedents

| Tool/Method | What It Does | Relevance |
|---|---|---|
| Mash/MinHash | Genome distance via k-mer sketches | Shows k-mer distributions characterise genomes |
| Kraken/Bracken | Metagenomic classification via k-mer matching | Target-specific k-mer counting for abundance estimation |
| novoBreak | Tumor-normal k-mer subtraction for SV detection | K-mer differences between samples detect structural changes |
| Jellyfish/KMC | Fast k-mer counting | Core computational primitive for this approach |
| FastQ Screen | Contamination detection via k-mer matching | Sample-level QC from raw reads |
| WisecondorX | NIPT CNV calling | Uses binned profiles + reference set (aligned, but similar normalisation logic) |
| FASTCN (Pendleton et al.) | Alignment-free CN estimation for segmental dups | Uses unique k-mer depth as CN proxy in duplicated regions |

### What's Novel in the Fingerprint Framing

The individual components (k-mer counting, GC distributions, fragment sizes)
are established. The contribution of the fingerprint concept is:

1. **Structured composition** of multiple alignment-free signals into a single
   per-sample descriptor optimised for CNV detection in panels
2. **Self-normalising** design where bias-descriptive layers (GC, fragment)
   enable correction without external reference alignment
3. **Fixed-dimensional** representation that enables direct sample-to-sample
   comparison, cohort-level analysis, and ML-based anomaly detection

---

## 9. Open Questions

1. **How many anchor k-mers per target are needed for robust CN estimation?**
   Theory suggests the variance of the median scales as ~1/(0.5 × N_anchors).
   Empirically, 20-50 anchors per target may suffice for single-copy
   resolution at typical panel coverage (500-1000x).

2. **Can variant-aware anchor sets improve robustness?** Including k-mers
   for common SNP alleles (from gnomAD) would prevent anchor dropout at
   polymorphic sites. The anchor set becomes population-aware.

3. **How well does the GC distribution shift predict single-target CNVs?**
   For small panels, each target contributes 0.5-2% of reads. A het deletion
   removes ~0.25-1% of reads at a specific GC. Is this detectable above GC
   distribution noise?

4. **Can the fingerprint detect mosaicism/subclonal CNVs?** The k-mer count
   sensitivity is fundamentally limited by the same Poisson statistics as
   alignment-based depth. High panel coverage (500-1000x) might enable
   detection of mosaic events at 10-20% cell fraction.

5. **What is the minimum panel size for useful GC/fragment distribution
   signals?** Very small panels (10-20 targets) may have too few targets
   to build meaningful distributions. The k-mer count vector remains
   useful regardless of panel size.

---

*Companion to: sources_of_rd_variance.md (what this approach must handle),
rd_normalization_approaches.md (alignment-based counterparts),
cnv_calling_mind_map.md (software landscape context).*
