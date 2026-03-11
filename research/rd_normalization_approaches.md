# Read-Depth CNV Calling: Normalization Approaches

> A detailed summary of normalization strategies used by read-depth CNV callers.
> Covers what each method corrects, which genomic regions are used as the
> normalization baseline, and key assumptions and limitations.

---

## The Core Problem

Raw read depth is a poor proxy for copy number. Systematic biases -- GC content,
mappability, capture efficiency, batch effects, replication timing -- often produce
larger signals than actual CNVs. Normalization must remove these biases while
preserving true copy number signal.

```
                    RAW DEPTH SIGNAL
                          │
        ┌────────┬────────┼────────┬────────┬────────┐
        ▼        ▼        ▼        ▼        ▼        ▼
    GC Bias   Mappab-  Capture   Batch    Replic.  True CN
              ility    Efficiency Effects  Timing   Signal
        │        │        │        │        │        │
     ~20-40%  ~5-15%   ~10-30%   ~5-20%   ~5-10%   What we
     of var.  of var.  of var.   of var.  of var.   actually
                       (WES)                         want
```

**Key assumption shared by all methods**: Most of the genome is diploid (CN=2), so
systematic patterns across bins/targets reflect bias, not biology. Methods differ in
how they estimate and remove that bias.

---

## 1. GC Bias Correction (LOESS/LOWESS)

```
GC Correction
│
├── What: Fit smooth curve of coverage vs GC content per sample,
│         then divide observed depth by expected depth from the curve
│
├── Regions used for normalization:
│   All autosomal bins/targets genome-wide
│   ├── WGS: typically 100bp-10kb bins, genome-wide
│   ├── WES: all captured exonic targets
│   └── Assumes bulk of regions are CN=2 → trend = bias, not biology
│
├── Tools: CNVnator, QDNAseq, HMMcopy, CNVkit, Control-FREEC
│
├── Strengths:
│   ├── Simple, fast, well-understood
│   ├── Works per-sample (no cohort needed)
│   └── Corrects the single largest source of RD variance
│
└── Limitations:
    ├── GC-coverage relationship isn't always smooth
    ├── Fragment-level GC effects differ from bin-level
    └── Residual "waves" from replication timing remain
```

---

## 2. Wave / Replication-Timing Correction

```
Wave Correction
│
├── What: Post-GC correction for long-range (~100kb-1Mb) systematic
│         coverage fluctuations correlated with replication timing
│
├── Regions used for normalization:
│   Genome-wide smoothed residuals after GC correction
│   ├── Rolling median/mean across large windows (hundreds of kb)
│   ├── Or pre-computed replication-timing reference profiles
│   └── Operates on autosomal regions; sex chromosomes handled separately
│
├── Tools: QDNAseq, CNVkit
│
├── Strengths:
│   ├── Removes periodic wave artifacts that GC correction misses
│   └── Particularly important for shallow WGS and FFPE samples
│
└── Limitations:
    ├── Can overcorrect real large CNVs if smoothing window too small
    └── Replication-timing profiles are cell-type dependent
```

---

## 3. PCA-Based Normalization

```
PCA / SVD Normalization
│
├── What: Build samples × targets matrix, run PCA/SVD, remove top
│         principal components that capture systematic noise
│         (GC, batch effects, capture probe efficiency)
│
├── Regions used for normalization:
│   All targets across entire cohort simultaneously
│   ├── Typically all autosomal exonic targets (WES) or bins (WGS)
│   ├── Cohort size: 30-100+ samples recommended
│   └── Key insight: systematic biases are correlated across samples;
│       true CNVs are rare and sample-specific → land in lower PCs
│
├── Tools: XHMM, CODEX2, CLAMMS, GATK gCNV (PCA denoising step)
│
├── How many PCs to remove?
│   ├── XHMM: user-specified (typically 10-30)
│   ├── CODEX2: BIC-based automatic selection
│   └── GATK gCNV: learned during model training
│
├── Strengths:
│   ├── Implicitly corrects GC, batch, and probe effects simultaneously
│   ├── No explicit bias model needed
│   └── Very effective for WES where per-target biases are complex
│
└── Limitations:
    ├── Requires a cohort (single samples cannot be PCA-normalized)
    ├── Common/polymorphic CNVs can load onto top PCs → removed as "noise"
    ├── Risk of overcorrection if too many PCs removed
    └── Cohort composition matters: mixed capture kits degrade performance
```

---

## 4. Panel of Normals (PoN)

```
Panel of Normals
│
├── What: Build a reference model from a set of known-normal samples;
│         normalize new test samples against this baseline
│
├── Regions used for normalization:
│   All targets/bins across the PoN samples
│   ├── PoN captures position-specific biases implicitly:
│   │   ├── Capture probe efficiency (per-target)
│   │   ├── Local GC content
│   │   ├── Mappability
│   │   └── Systematic artifacts specific to the platform
│   ├── Typical PoN size: 20-50 normals (more is better, diminishing after ~40)
│   └── Must match test sample: same capture kit, library prep, sequencer
│
├── Tools: GATK (CreateReadCountPanelOfNormals), CNVkit (reference),
│          ichorCNA
│
├── Strengths:
│   ├── Implicitly handles multiple bias types at once
│   ├── New samples normalized individually (no cohort reprocessing)
│   └── Industry standard for somatic pipelines (tumor vs PoN)
│
└── Limitations:
    ├── PoN must closely match test sample platform/protocol
    ├── Stale PoNs (from old batches) degrade performance over time
    └── Building a good PoN requires curating truly normal samples
```

---

## 5. Reference Sample Optimization

```
Automatic Reference Selection
│
├── What: Instead of a fixed PoN, automatically select the best-matching
│         subset of reference samples per test sample based on correlation
│
├── Regions used for normalization:
│   All targets used to compute pairwise correlation
│   ├── Step 1: compute correlation of depth profiles across all targets
│   │           between test sample and every candidate reference
│   ├── Step 2: select top-N most correlated samples as the reference set
│   ├── Step 3: build sample-specific reference from the selected normals
│   └── Typically autosomal targets only; sex chromosomes separate
│
├── Tools: ExomeDepth, DECoN
│
├── Strengths:
│   ├── Adapts to per-sample variation (batch, GC, capture differences)
│   ├── More robust than fixed PoN when samples come from mixed batches
│   └── No manual PoN curation needed
│
└── Limitations:
    ├── Requires a pool of candidate references at analysis time
    ├── Selection can be unstable with small reference pools (<10)
    └── Computationally heavier than fixed PoN (pairwise correlations)
```

---

## 6. Poisson Latent Factor Model

```
Poisson Latent Factor
│
├── What: Model the samples × targets count matrix as Poisson-distributed
│         with latent factors representing systematic biases. More principled
│         than PCA for count data (which is non-Gaussian).
│
├── Regions used for normalization:
│   All targets across cohort (same matrix as PCA)
│   ├── But accounts for mean-variance relationship of count data
│   ├── Low-coverage targets appropriately down-weighted
│   └── Avoids PCA's implicit Gaussian assumption
│
├── Tools: CODEX2 (offers both PCA and Poisson factor modes)
│
├── Strengths:
│   ├── Statistically principled for sequencing count data
│   ├── Handles heteroscedasticity (variance scales with mean)
│   └── BIC-based automatic selection of number of latent factors
│
└── Limitations:
    ├── Computationally heavier than PCA/SVD
    ├── Still requires cohort
    └── Minimal practical improvement over PCA in most benchmarks
```

---

## 7. On-Target + Off-Target Normalization (CNVkit)

```
On-Target + Off-Target (CNVkit-specific)
│
├── What: Use scattered off-target reads (covering much of the genome)
│         as an additional normalization signal alongside on-target counts
│
├── Regions used for normalization:
│   ├── ON-TARGET: captured exonic regions (high coverage, sparse genomic)
│   ├── OFF-TARGET: genomic regions between capture targets
│   │   ├── Reads that land outside targets during hybrid capture
│   │   ├── Tiled into large bins (avg ~5kb) to accumulate sufficient reads
│   │   └── Provide genome-wide baseline unavailable from targets alone
│   └── Both combined into a unified log2 ratio profile
│
├── Tools: CNVkit (unique to this tool)
│
├── Strengths:
│   ├── Genome-wide coverage context from a WES experiment
│   ├── Off-target bins help correct biases that sparse targets can't reveal
│   └── Enables detection of large off-target CNVs from WES data
│
└── Limitations:
    ├── Off-target coverage is low and noisy (~0.5-2x)
    ├── Only works with hybrid capture (not amplicon-based panels)
    └── Off-target bin size must be tuned per experiment
```

---

## 8. Mappability Filtering / Correction

```
Mappability Correction
│
├── What: Exclude or adjust depth in regions where short reads cannot
│         map uniquely, producing unreliable depth estimates
│
├── Regions used:
│   Pre-computed genome-wide mappability tracks
│   ├── Calculated as: fraction of unique k-mers per bin
│   │   (k = read length, typically 100-150bp)
│   ├── Bins below threshold (e.g., mappability < 0.5) excluded
│   ├── Or: depth adjusted by dividing by mappability score
│   └── Problematic regions: segmental duplications, centromeric
│       repeats, recent paralogs, collapsed reference sequences
│
├── Tools:
│   ├── Control-FREEC: explicit mappability correction (adjusts counts)
│   ├── GATK: uses mappability in PoN construction
│   └── Most tools: implicit via blacklist filtering (see below)
│
├── Strengths:
│   ├── Prevents false calls in regions with inflated/ambiguous depth
│   └── Essential for WGS where unmappable regions are a major problem
│
└── Limitations:
    ├── Aggressive filtering loses real CNVs in segmental duplications
    ├── Mappability tracks are reference-build specific
    └── T2T/CHM13 reference reduces (but doesn't eliminate) the problem
```

---

## 9. Blacklist / Region Exclusion

```
Blacklist Filtering
│
├── What: Exclude known problematic regions before any normalization,
│         preventing extreme outliers from distorting bias estimates
│
├── Regions excluded:
│   ├── ENCODE blacklist (DAC + Duke merged lists)
│   │   └── Regions with anomalous signal across many cell types/assays
│   ├── Centromeres and pericentromeric regions
│   ├── Telomeres and subtelomeric regions
│   ├── Satellite repeat arrays
│   ├── Decoy/alt sequences (if using GRCh38 with alts)
│   └── Sex chromosomes (handled separately with sex-matched references)
│
├── Tools: Nearly all tools (standard preprocessing step)
│
├── Strengths:
│   ├── Prevents gross normalization distortion from outlier regions
│   ├── Community-maintained lists (ENCODE, Heng Li's blacklist)
│   └── Fast, simple, uncontroversial
│
└── Limitations:
    ├── Overly aggressive blacklists remove real biology
    ├── Lists are reference-build specific (GRCh37 vs GRCh38 vs T2T)
    └── Custom capture kits may need custom blacklists
```

---

## Comparison: Which Tools Use Which Methods

```
                    GC     Wave   PCA    PoN    Ref-   Poisson  On+Off  Mapp.  Black-
                    LOESS  Corr.  /SVD          Select  Latent  Target  Corr.  list
                    ─────  ─────  ─────  ─────  ─────  ──────  ──────  ─────  ─────
CNVnator            ✓                                                          ✓
QDNAseq             ✓      ✓                                                   ✓
HMMcopy             ✓                                                   ✓      ✓
Control-FREEC       ✓                                                   ✓      ✓
CNVkit              ✓      ✓             ✓                       ✓             ✓
XHMM                              ✓                                            ✓
CODEX2                            ✓                    ✓                       ✓
CLAMMS                            ✓                                            ✓
ExomeDepth                                      ✓                              ✓
DECoN                                           ✓                              ✓
GATK gCNV                        ✓      ✓                                     ✓
panelcn.MOPS                                                                   ✓
ichorCNA            ✓                    ✓                                     ✓
FACETS              ✓                                                          ✓
BIC-seq2            ✓                                                          ✓
```

---

## Order of Operations

Most tools apply normalization in a specific sequence. The general pattern:

```
Typical Normalization Pipeline
│
├── 1. Blacklist filtering
│      Remove known-bad regions before anything else
│
├── 2. Mappability filtering
│      Exclude or flag low-mappability bins
│
├── 3. GC bias correction (LOESS)
│      Remove the dominant single source of variance
│
├── 4. Wave / replication-timing correction
│      Remove long-range residual patterns post-GC
│
├── 5. Systematic bias removal (choose one or combine):
│      ├── PCA across cohort (WES, large cohorts)
│      ├── PoN subtraction (WGS somatic, clinical)
│      ├── Reference selection (WES, clinical panels)
│      └── Poisson latent factor (WES, principled alternative to PCA)
│
├── 6. Log2 ratio transformation
│      Convert normalized counts to log2(sample/reference) ratios
│
└── 7. Segmentation
       CBS, HMM, or BIC-based partitioning of the normalized profile
```

**Key insight**: Steps 1-4 are largely universal and uncontroversial. Step 5 is where
tools diverge most and where tool selection matters most. The choice depends primarily
on whether you have a cohort (PCA, Poisson factor), a matched PoN (PoN subtraction),
or neither (reference selection, per-sample GC-only).

---

*Companion to: cnv_calling_mind_map.md. Covers normalization in depth; see the main
mind map for tool selection, ensemble strategies, and benchmarking.*
