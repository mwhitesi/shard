# Short-Read CNV Calling: Software Landscape Mind Map

> A high-level survey of software approaches for copy number variation detection
> from short-read sequencing data. Organised by signal type, application context,
> and common challenges.

---

## 1. Signal Types Used for CNV Detection

```
                        ┌─────────────────────────────────┐
                        │   SHORT-READ CNV CALLING SIGNALS │
                        └──────────────┬──────────────────┘
               ┌───────────┬───────────┼───────────┬───────────────┐
               ▼           ▼           ▼           ▼               ▼
          Read Depth   Paired-End   Split-Read   Assembly    B-Allele Freq
            (RD)         (PE)         (SR)                     (BAF)
               │           │           │           │               │
          Coverage     Discordant   Soft-clipped  Local/global  Het SNP
          in bins      insert size  read halves   de novo       allele
                       or orient.   map to diff.  assembly      ratios
                                    locations     of contigs
```

### 1.1 Read Depth (RD)
- **Principle**: CNVs alter local coverage -- deletions reduce depth, duplications increase it
- **Resolution**: Bin-dependent (100 bp -- 10 kb); cannot pinpoint breakpoints
- **Strengths**: Detects all sizes >1 kb; works for WGS, WES, and panels
- **Weaknesses**: Confounded by GC bias, mappability, batch effects; poor breakpoint resolution; duplications harder to detect than deletions (1.5x vs 0.5x signal for hets)

### 1.2 Paired-End / Read-Pair (PE)
- **Principle**: Discordant pairs (abnormal insert size or orientation) indicate SVs
- **Resolution**: ~50-150 bp (insert size SD)
- **Strengths**: Detects inversions and translocations (not just CNV); orientation-aware
- **Weaknesses**: Insensitive to small events; limited in WES/panels (breakpoints usually off-target)

### 1.3 Split-Read (SR)
- **Principle**: Soft-clipped reads spanning breakpoints reveal exact junction sequences
- **Resolution**: 1-10 bp
- **Strengths**: Best breakpoint resolution of alignment-based methods
- **Weaknesses**: Requires breakpoint in unique sequence; fails in repeats

### 1.4 Local Assembly
- **Principle**: De novo assemble reads around candidate breakpoints into contigs
- **Resolution**: Base-pair (if assembly succeeds)
- **Strengths**: Resolves complex events; finds novel insertions
- **Weaknesses**: Computationally expensive; assembly can fail in low-complexity regions

### 1.5 B-Allele Frequency (BAF)
- **Principle**: Heterozygous SNP allele ratios shift with copy number changes (e.g., LOH shifts BAF from 0.5 to 0 or 1)
- **Strengths**: Detects CN-neutral LOH; resolves purity/ploidy degeneracy in somatic calling; confirms RD-based calls
- **Weaknesses**: Requires heterozygous SNPs (density varies); resolution limited by SNP spacing

---

## 2. Software by Primary Approach

### 2.1 Read-Depth-Only Tools

```
RD-Only Tools
├── WGS-focused
│   ├── CNVnator ──── mean-shift segmentation, GC-corrected bins (2011)
│   ├── Control-FREEC ── LASSO seg + GC/mappability correction (2012)
│   ├── BIC-seq2 ──── BIC-based segmentation, base-pair resolution (2016)
│   └── QDNAseq ───── CBS segmentation, excellent for shallow WGS (2014)
│
├── WES-focused
│   ├── XHMM ──────── PCA normalization + Gaussian HMM, 3 states (2012)
│   ├── cn.MOPS ───── Bayesian Poisson mixture, cohort-based (2012)
│   ├── CODEX2 ────── Poisson latent-factor norm + BIC seg (2018)
│   ├── ExomeDepth ── beta-binomial + HMM, auto reference selection (2012)
│   ├── CLAMMS ────── lattice-aligned mixture + HMM, large cohorts (2016)
│   └── GATK gCNV ── hierarchical HMM + NB, PCA denoising (2019-2023)
│
├── Panel-focused
│   ├── panelcn.MOPS ── Poisson mixture adapted for small panels (2017)
│   ├── DECoN ─────── ExomeDepth-based + correlation QC (2016)
│   └── SavvySuite ── off-target reads + HMM (2023)
│
└── Somatic-specific
    ├── CNVkit ────── on+off-target log2 ratio, CBS/HMM (2016) ★ top WES somatic
    ├── FACETS ────── joint purity/ploidy + allele-specific CN (2016)
    ├── PureCN ────── Bayesian, coverage+BAF+mutations (2016)
    ├── ASCAT ─────── BAF+LRR joint purity/ploidy estimation (2010)
    ├── Sequenza ─── BAF+depth grid search (2015)
    ├── ABSOLUTE ─── probabilistic purity/ploidy model (2012)
    ├── Battenberg ── phased BAF, subclonal CN (2014)
    └── ichorCNA ─── ultra-low-pass WGS / cfDNA (2017)
```

### 2.2 Split-Read Tools
| Tool | Year | Notes |
|------|------|-------|
| Pindel | 2009 | Pattern-growth algorithm; one of the earliest SR tools |
| CREST | 2011 | Assembles soft-clips into contigs; cancer genomes |
| SplazerS | 2012 | Split-read alignment with gapped aligner |
| Gustaf | 2014 | Multi-breakpoint split-read alignment |

### 2.3 Paired-End Tools
| Tool | Year | Notes |
|------|------|-------|
| BreakDancer | 2009 | Foundational PE caller; Poisson model on discordant pairs |
| PEMer | 2009 | Cluster-based discordant pair analysis |
| GASV/GASVPro | 2009/2012 | Geometric breakpoint-space intersection; GASVPro adds Bayesian model |
| HYDRA | 2010 | Handles mapping ambiguity with multiple alignments |
| CLEVER | 2012 | Max-clique enumeration on insert-size interval graph |

### 2.4 Assembly-Based Tools
| Tool | Year | Scope | Notes |
|------|------|-------|-------|
| fermikit | 2015 | WGS germline | Whole-genome de novo assembly; high specificity |
| Manta | 2016 | WGS/WES both | Breakpoint graph + local assembly; fast and accurate ★ |
| TIGRA | 2014 | WGS | Targeted de Bruijn graph assembly for refinement |
| SvABA | 2018 | WGS somatic | String-graph local assembly; bp-level resolution |
| novoBreak | 2017 | WGS somatic | K-mer subtraction between tumor/normal + assembly |

### 2.5 Hybrid / Multi-Signal Tools

```
Multi-Signal Integration
├── PE + SR
│   ├── LUMPY ────── breakpoint probability framework (2014)
│   └── Wham ──────── anomalous read distribution testing (2015)
│
├── PE + SR + Assembly
│   ├── Manta ────── breakpoint graph + local assembly (2016) ★
│   ├── GRIDSS2 ──── positional de Bruijn graph + quality scoring (2017/2022)
│   └── SvABA ────── string-graph assembly at breakpoints (2018)
│
├── PE + SR + RD
│   ├── DELLY2 ───── Bayesian genotyping across 3 signals (2012/2023) ★
│   └── TIDDIT ───── clustering + coverage analysis (2017)
│
├── RD + BAF
│   ├── ERDS ──────── HMM with depth + SNV zygosity (2012)
│   ├── CNVpytor ── successor to CNVnator + BAF (2021)
│   └── PURPLE ──── purity/ploidy-aware somatic CN (Hartwig) ★
│
└── Full Ensemble Pipelines
    ├── GATK-SV ──── Manta+DELLY+MELT+Wham+gCNV → RF filtering (2020/2023) ★★
    ├── Parliament2 ── Manta+DELLY+LUMPY+BreakDancer+CNVnator → SURVIVOR (2020)
    └── GRIDSS2/PURPLE/LINX ── full somatic SV+CN+interpretation (2022) ★★
```

### 2.6 Machine Learning / Deep Learning Approaches

| Tool | Model | Approach | Year | Status |
|------|-------|---------|------|--------|
| DeepSV | CNN | Encodes alignment signals as image matrices | 2019 | Research |
| SVcnn | CNN | Pileup images for SV classification | 2020 | Research |
| CNV-Net | CNN+RNN | Coverage profiles → boundary detection | 2020 | Research |
| CNV-espresso | CNN | Post-hoc validation/refinement of CNV calls | 2023 | Emerging |
| Cue | Multi-channel CNN | Multiple signals as image channels | 2024 | Promising ★ |
| SVFormer | Transformer | Attention on alignment features | 2023-24 | Emerging |
| GATK-SV (RF) | Random Forest | Filters merged ensemble calls | 2020 | Production ★ |

**Key insight**: ML/DL has had more success in *filtering* CNV calls than in *primary calling*. Cue (2024) is the most competitive standalone DL approach.

### 2.7 Ensemble / Consensus Merging Tools

| Tool | Strategy | Year |
|------|----------|------|
| SURVIVOR | Rule-based merging by proximity, type, size | 2017 ★ most widely used |
| Jasmine | Improved distance-based merging | 2023 |
| MetaSV | Re-evaluates candidates via local assembly | 2015 |
| FusorSV | Bayesian weighting of caller accuracy | 2018 |
| Truvari | Benchmarking + merge/collapse | 2022 |
| smoove | LUMPY + SVTyper + duphold pipeline | ongoing |

---

## 3. Comparison by Application Context

### 3.1 Recommended Tools by Assay + Application

```
                    ┌─────────── GERMLINE ───────────┐  ┌────────── SOMATIC ──────────┐
                    │                                 │  │                              │
WGS (30-60x)       │ Manta + CNVnator/DELLY          │  │ GRIDSS2/PURPLE/LINX          │
                    │ GATK-SV (gold standard)         │  │ Battenberg, ASCAT            │
                    │ smoove (population-scale)        │  │ Control-FREEC, CNVkit        │
                    │                                 │  │                              │
WES (100-300x)     │ GATK gCNV ★ (large cohort)      │  │ CNVkit ★ (top performer)     │
                    │ ExomeDepth ★ (rare CNV)          │  │ FACETS, PureCN               │
                    │ XHMM, CODEX2 (cohort)           │  │ Sequenza, GATK ModelSeg      │
                    │                                 │  │                              │
Targeted Panel     │ ExomeDepth, DECoN                │  │ CNVkit, PureCN               │
                    │ panelcn.MOPS, SavvySuite         │  │ Custom approaches            │
                    │                                 │  │                              │
Shallow WGS        │ QDNAseq                          │  │ ichorCNA ★ (cfDNA/tumor)     │
(0.1-1x)           │ WisecondorX (NIPT)              │  │ QDNAseq                      │
                    └─────────────────────────────────┘  └──────────────────────────────┘
```

### 3.2 Somatic-Specific: Purity/Ploidy Estimation

```
Purity/Ploidy Challenge
│
├── The Problem: purity and ploidy are confounded in RD alone
│   (50% pure tetraploid ≈ 100% pure diploid)
│
├── Solution: integrate BAF from heterozygous SNPs
│
└── Tools ranked by approach:
    ├── ASCAT ──── BAF + LRR → joint estimation (widely used)
    ├── FACETS ── two-pass EM + segmentation (popular for WES)
    ├── ABSOLUTE ── probabilistic model + mutation AFs
    ├── PureCN ── Bayesian, coverage + BAF + somatic mutations
    ├── Battenberg ── phased BAF → subclonal detection
    ├── Sequenza ── grid search over purity/ploidy space
    └── PURPLE ── integrated with GRIDSS2 pipeline (WGS gold standard)
```

### 3.3 Special Applications

```
cfDNA / Liquid Biopsy
├── ichorCNA ★ ── ULP-WGS, estimates tumor fraction + arm-level CN
├── DELFI ──────── fragmentomic approach (size + coverage patterns)
├── WisecondorX ── adapted from NIPT for cancer monitoring
└── Challenge: ctDNA fraction often 0.1-10%; CHIP confounds results

NIPT / Prenatal
├── WisecondorX ★ ── reference-set based, trisomies + sub-chromosomal
├── NIPTeR ────────── R package for NIPT analysis
└── Challenge: confined placental mosaicism → false positives;
    fetal fraction estimation critical; sub-chromosomal CNVs harder

Clinical Diagnostics (panels)
├── ExomeDepth ★ ── auto reference selection, beta-binomial model
├── DECoN ─────── ExomeDepth-based + correlation QC for panels
├── panelcn.MOPS ── optimised for small gene panels
└── Challenge: single-exon calls need orthogonal confirmation (MLPA/ddPCR)
```

---

## 4. Common Challenges

### 4.1 Technical Challenges Mind Map

```
                         CHALLENGES
        ┌────────┬────────┬────────┬────────┬────────┐
        ▼        ▼        ▼        ▼        ▼        ▼
    GC Bias   Mappab-  Segmental  Batch   Reference  Breakpoint
              ility    Duplicat.  Effects  Bias       Resolution
        │        │        │        │        │        │
   LOESS,PCA  Masks,   SUN/PSN   Matched  T2T-     SR: 1-10bp
   PoN,frag-  MAPQ     methods,  PoN,PCA  CHM13    PE: 50-150bp
   level GC   filters  long-read  residual resolves RD: bin-level
   correction          needed    ization  many     Assembly: bp
```

### 4.2 GC Bias Correction Approaches

| Method | Used By | Best For |
|--------|---------|----------|
| LOESS/LOWESS regression | CNVnator, QDNAseq, HMMcopy | WGS; simple first-pass |
| PCA-based removal | XHMM, CODEX2, CLAMMS | WES/panels with cohorts |
| Panel of Normals | GATK, ichorCNA | Any platform; implicit GC correction |
| Poisson latent factor | CODEX2 | WES; principled for count data |
| HMM with GC covariates | HMMcopy, ichorCNA | Direct modeling |
| Wave correction (post-GC) | QDNAseq, CNVkit | Replication-timing artifacts |

### 4.3 Size-Dependent Performance

```
CNV Size Spectrum
│
│  < 300 bp        SR + assembly only; confused with indels
│  300 bp - 1 kb   PE + SR; depth signal too weak
│  1 kb - 100 kb   ★ SWEET SPOT: PE + SR + emerging RD signal
│  100 kb - 1 Mb   RD excels; PE/SR provide breakpoints
│  > 1 Mb          RD + BAF; arm-level events
│  Whole chrom.     Aneuploidy detection; explicit ploidy modeling
│
└── KEY: No single method covers all sizes → ensemble needed
```

### 4.4 Deletions vs Duplications

- **Deletions** are easier: het deletion = 0.5x depth (50% signal change)
- **Duplications** are harder: het duplication = 1.5x depth (50% signal change, but from a higher baseline → lower relative change in noise)
- BAF integration substantially helps duplication detection
- Benchmarks consistently show lower recall for duplications across all tools

---

## 5. Normalization & Preprocessing

```
Preprocessing Pipeline
│
├── 1. Read Alignment (BWA-MEM/BWA-MEM2)
│
├── 2. Duplicate Marking (Picard/samblaster)
│
├── 3. Coverage Computation
│   ├── Bin reads into windows (WGS: 100bp-10kb)
│   ├── Count reads per target (WES/panels)
│   └── Separate on-target vs off-target (CNVkit)
│
├── 4. GC Correction
│   ├── LOESS regression (most tools)
│   └── + Wave/replication-timing correction (QDNAseq, CNVkit)
│
├── 5. Systematic Bias Removal
│   ├── Panel of Normals (GATK, CNVkit)
│   ├── PCA residualization (XHMM, CODEX2)
│   ├── Mappability filtering (exclude low-mappability bins)
│   └── Blacklist filtering (centromeres, telomeres, ENCODE blacklist)
│
├── 6. Segmentation
│   ├── CBS (Circular Binary Segmentation) ── CNVkit, QDNAseq
│   ├── HMM ── GATK gCNV, XHMM, ExomeDepth, DELLY
│   ├── Mean-shift ── CNVnator
│   ├── BIC-based ── BIC-seq2, CODEX2
│   └── LASSO ── Control-FREEC
│
└── 7. Copy Number Calling / Classification
    ├── Integer CN states (DEL/NEUT/DUP)
    ├── Continuous log2 ratios
    └── Allele-specific CN (major/minor allele)
```

---

## 6. Benchmarking & Gold Standards

### 6.1 Reference Resources

| Resource | Description | Use |
|----------|-------------|-----|
| GIAB HG002 SV benchmark v0.6 | ~7,200 high-confidence SVs (Zook et al., 2022) | Primary germline benchmark |
| gnomAD-SV v2/v4 | ~433K SVs from 14,891+ WGS samples (Collins et al., 2020) | Population frequency filtering |
| 1000 Genomes Phase 3 SVs | ~68K SVs across 2,504 individuals (Sudmant et al., 2015) | Population reference |
| BAMsurgeon | Spike SVs into real BAMs (Ewing et al., 2015) | Somatic benchmarking |
| DREAM Challenge | Community somatic calling benchmark | Standardised comparison |

### 6.2 Key Benchmarking Papers

| Study | Year | Key Finding |
|-------|------|-------------|
| Kosugi et al., *Genome Biology* | 2019 | Compared 69 SV tools; no single tool dominates; ensemble recommended |
| Zare et al., *BMC Bioinformatics* | 2017 | CNVkit + Control-FREEC best for somatic WES |
| Roca et al., *Nature Biotechnology* | 2019 | ExomeDepth + GATK gCNV best for clinical germline WES |
| Moreno-Cabrera et al., *Human Mutation* | 2020 | ExomeDepth + DECoN best for diagnostic panels |

### 6.3 Concordance Metrics

- **Precision** (PPV), **Recall** (Sensitivity), **F1 score**
- **Reciprocal overlap** (standard: >=50%)
- **Breakpoint distance** (for SR/PE methods)
- **Size-stratified** and **region-stratified** metrics essential
- Tools: **Truvari** (standard), SVanalyzer, Wittyer

---

## 7. Key Takeaways

### 7.1 Consensus from the Field

1. **No single tool is sufficient** -- all benchmarks converge on this. Ensemble approaches with 2-3 callers consistently outperform any individual tool.

2. **Tool selection is assay-dependent:**
   - WGS germline → GATK-SV or Manta+DELLY+CNVnator
   - WGS somatic → GRIDSS2/PURPLE/LINX or Battenberg
   - WES germline → GATK gCNV + ExomeDepth
   - WES somatic → CNVkit (consistently top) + FACETS
   - Panels → ExomeDepth or panelcn.MOPS
   - Shallow WGS → QDNAseq (germline) or ichorCNA (somatic)

3. **Multi-signal integration is the future** -- the field has moved from single-signal tools (2009-2014) to multi-signal (2015-2020) to full ensemble pipelines (2020+).

4. **ML/DL is filtering, not yet calling** -- random forests (GATK-SV) succeed at filtering; standalone DL callers (Cue, DeepSV) are emerging but not yet production-ready.

5. **Duplications remain the Achilles' heel** -- consistently lower recall than deletions across all tools and signal types.

6. **Cohort size matters enormously** for WES tools -- GATK gCNV, XHMM, cn.MOPS gain substantial power at 30-100+ samples.

7. **BAF is essential for somatic calling** -- resolves purity/ploidy degeneracy, detects CN-neutral LOH, confirms RD calls.

8. **Short reads have fundamental blind spots** in segmental duplications and centromeric regions -- long-read or targeted methods needed for these regions.

### 7.2 Emerging Trends (2023-2025)

- **GATK-SV** as the reference germline WGS pipeline (gnomAD, All of Us)
- **GRIDSS2/PURPLE/LINX** as the reference somatic WGS pipeline (Hartwig)
- **ClinCNV** for clinical all-assay CNV calling
- **T2T reference genomes** reducing reference bias
- **Transformer architectures** (SVFormer) beginning to appear
- **Fragmentomics** (DELFI) adding new signal types for cfDNA
- Move toward **pangenome references** reducing CNV calling artifacts

---

## 8. Tool Quick-Reference Table

| Tool | Signal(s) | Model | Assay | Germ/Som | Year | GitHub/Source |
|------|-----------|-------|-------|----------|------|---------------|
| CNVnator | RD | Mean-shift | WGS | G | 2011 | abyzovlab/CNVnator |
| Control-FREEC | RD+mapp | LASSO | WGS/WES | G+S | 2012 | BoevaLab/FREEC |
| cn.MOPS | RD cohort | Bayesian Poisson | WGS/WES | G | 2012 | Bioconductor |
| QDNAseq | RD | CBS | lpWGS | S | 2014 | Bioconductor |
| CNVkit | RD on+off | CBS/HMM | WES/panels | G+S | 2016 | etal/cnvkit |
| GATK gCNV | RD PCA | Hierarchical HMM | WES/WGS | G | 2019+ | broadinstitute/gatk |
| CODEX2 | RD latent | Poisson factor + BIC | WES | G+S | 2018 | Bioconductor |
| ExomeDepth | RD ref-opt | Beta-binomial+HMM | WES/panels | G | 2012 | CRAN |
| XHMM | RD PCA | Gaussian HMM | WES | G | 2012 | atgu.mgh.harvard.edu |
| CLAMMS | RD matched | Mixture+HMM | WES | G | 2016 | rgcgithub/clamms |
| panelcn.MOPS | RD cohort | Poisson mixture | Panels | G | 2017 | Bioconductor |
| DECoN | RD | ExomeDepth-based | Panels | G | 2016 | RahmanTeam/DECoN |
| BIC-seq2 | RD | BIC segmentation | WGS | G+S | 2016 | Harvard |
| Pindel | SR | Pattern-growth | WGS/WES | G+S | 2009 | genome/pindel |
| BreakDancer | PE | Poisson on discordant | WGS | G+S | 2009 | genome/breakdancer |
| Manta | PE+SR+asm | Breakpoint graph | WGS/WES | G+S | 2016 | Illumina/manta |
| DELLY2 | PE+SR+RD | Bayesian | WGS | G+S | 2012/23 | dellytools/delly |
| LUMPY | PE+SR | Breakpoint probability | WGS | G+S | 2014 | arq5x/lumpy-sv |
| GRIDSS2 | PE+SR+asm | Quality scoring | WGS | G+S | 2017/22 | PapenfussLab/gridss |
| SvABA | SR+asm | String-graph+LRT | WGS | G+S | 2018 | walaj/svaba |
| GATK-SV | All signals | RF ensemble | WGS | G | 2020/23 | broadinstitute/gatk-sv |
| Parliament2 | Multi-caller | SURVIVOR merge | WGS | G | 2020 | dnanexus/parliament2 |
| FACETS | RD+BAF | EM+segmentation | WES | S | 2016 | mskcc/facets |
| ASCAT | RD+BAF | Joint purity/ploidy | WGS/WES | S | 2010 | VanLoo-lab/ascat |
| PureCN | RD+BAF+mut | Bayesian | WES/panels | S | 2016 | lima1/PureCN |
| PURPLE | RD+BAF | Purity/ploidy fit | WGS | S | 2019+ | hartwigmedical/hmftools |
| Battenberg | Phased BAF | Subclonal CN | WGS | S | 2014 | Wedge-lab/battenberg |
| ichorCNA | RD | HMM (low-pass) | lpWGS/cfDNA | S | 2017 | broadinstitute/ichorCNA |
| SURVIVOR | Merging | Rule-based consensus | Any | Any | 2017 | fritzsedlazeck/SURVIVOR |
| Truvari | Benchmark | Matching+metrics | Any | Any | 2022 | ACEnglish/truvari |
| ClinCNV | RD+BAF | HMM multi-sample | All | G+S | 2022+ | imgag/ClinCNV |
| CNVpytor | RD+BAF | Mean-shift | WGS | G | 2021 | abyzovlab/CNVpytor |
| Cue | Multi-signal | CNN | WGS | G | 2024 | Emerging |

---

*Document compiled from literature review spanning 2009-2025. Key references include Kosugi et al. (2019) Genome Biology, Collins et al. (2020) Nature, Babadi et al. (2023) Nature Genetics, Roca et al. (2019) Nature Biotechnology, and Zare et al. (2017) BMC Bioinformatics. For the latest tools, search PubMed for "CNV detection short read" and GitHub for "CNV calling" sorted by recently updated.*
