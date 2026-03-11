# Sources of Read-Depth Variance in Short-Read CNV Calling

> A systematic catalog of the factors that obscure the true copy number signal
> in read depth data. Organised from the physics of library prep through to
> computational artifacts. Each source is characterised by its magnitude,
> spatial scale, and how (or whether) it can be corrected.

---

## The Signal We Want vs The Signal We Get

A read-depth CNV caller treats binned read counts as a noisy proxy for integer
copy number. For a heterozygous deletion at diploid baseline, the expected
signal is a 50% drop in coverage. For a heterozygous duplication, a 50%
increase. In practice, these signals are buried under multiple layers of
systematic and stochastic variation that can easily produce 2-3x fold-changes
in coverage across the genome *at constant copy number*.

```
                         OBSERVED READ DEPTH
                                │
     ┌──────────┬───────────┬───┴───┬───────────┬──────────┬──────────┐
     ▼          ▼           ▼       ▼           ▼          ▼          ▼
  Library    Sequencing   Capture  Alignment  Reference  Stochastic  TRUE
  Prep       Chemistry    (WES)    & Mapping  Genome     Sampling    COPY
  Biases     Biases       Biases   Biases     Biases     Noise       NUMBER
     │          │           │       │           │          │
  GC bias,   Flow cell   Probe    Multi-map,  Collapsed  Poisson
  PCR dup,   position,   hybrid.  MAPQ, soft  SDs, alt   counting
  fragment   cycle err,  kinetics clip, read  haplotype  noise,
  length,    index hop    & decay  length     gaps       library
  FFPE deg                                               complexity
```

---

## 1. Library Preparation Biases

### 1.1 GC Content Bias

- **What**: Non-uniform amplification and cluster formation as a function of
  fragment GC content. Both extremes (AT-rich and GC-rich) are under-represented.
- **Magnitude**: Typically the single largest source of RD variance. Can produce
  2-4x coverage differences between GC-extreme and GC-moderate bins.
  Accounts for ~20-40% of total depth variance in a typical WGS sample.
- **Spatial scale**: Correlated over ~100bp-10kb (the scale of local GC content
  variation). Creates a smooth, genome-wide modulation of coverage.
- **Root cause**: Multiple compounding mechanisms:
  - **PCR amplification**: Denaturation efficiency varies with GC%; high-GC
    templates re-anneal faster and amplify less efficiently. AT-rich templates
    may denature too readily, creating secondary structures.
  - **Bridge amplification** (Illumina): Cluster generation on the flow cell
    surface recapitulates PCR GC bias.
  - **Tagmentation** (Nextera): Tn5 transposase has slight insertion site
    preferences correlated with local sequence context.
- **Key subtlety**: The GC bias is a property of the *fragment*, not the
  *genomic bin*. A bin's GC content is a proxy for the GC of the fragments
  overlapping it, but fragments span variable lengths and the relationship is
  imperfect. Fragment-level GC correction (used by some newer tools) is more
  accurate than bin-level.
- **Correction**: LOESS/LOWESS regression of coverage vs GC per sample.
  Effective but leaves residuals, especially at GC extremes where data is sparse.

### 1.2 PCR Duplication

- **What**: Identical fragments amplified during library PCR produce duplicate
  reads that inflate apparent depth at some loci more than others.
- **Magnitude**: Duplicate rates typically 5-30% of total reads. After marking
  and removal (Picard/samblaster), the remaining unique reads reflect the
  *library complexity*, which sets a ceiling on usable depth.
- **Spatial scale**: Affects individual loci stochastically. The propensity for
  duplication correlates with fragment GC and length, so it compounds GC bias.
- **Root cause**: PCR cycles amplify some templates more than others based on
  GC content, fragment length, and secondary structure.
- **Correction**: Duplicate marking (coordinate-based or UMI-based). Standard
  practice. Optical/clustering duplicates (from patterned flow cells) are
  distinct from PCR duplicates and handled separately.
- **Residual risk**: After duplicate removal, the effective coverage may vary
  across the genome in a GC-correlated pattern. Low-complexity libraries
  (low input DNA) are especially affected -- duplicate removal leaves
  coverage deserts.

### 1.3 Fragment Length Distribution

- **What**: The insert size distribution of the sequencing library is not
  uniform across the genome. Fragment length correlates with local sequence
  features (GC, nucleosome positioning, chromatin accessibility).
- **Magnitude**: Modest direct effect on RD, but fragments at the tails of the
  size distribution may be lost during size selection, creating subtle
  coverage variation (~2-5%).
- **Spatial scale**: ~150bp-1kb (nucleosome scale). Particularly relevant for
  cfDNA, where fragment sizes are biologically determined by nucleosome
  protection.
- **Root cause**: Sonication/enzymatic fragmentation is not perfectly random.
  AT-rich sequences shear more easily. Size selection during library prep
  removes fragments outside the target range, and the retained fraction varies
  by locus.
- **Correction**: Not directly corrected by most CNV callers. Fragment-level GC
  correction partially accounts for this. cfDNA-specific tools (DELFI,
  ichorCNA) explicitly model fragmentomics.

### 1.4 DNA Quality / FFPE Degradation

- **What**: Degraded input DNA (especially formalin-fixed paraffin-embedded
  tissue) produces shorter, damaged fragments with sequence artifacts.
- **Magnitude**: Severe. FFPE samples can show 2-5x higher coverage variance
  compared to fresh-frozen tissue. Deamination (C→T) artifacts further
  corrupt the signal.
- **Spatial scale**: Genome-wide increase in noise floor, plus locus-specific
  dropout where DNA is too damaged to amplify.
- **Root cause**: Formalin cross-links fragment DNA and cause base damage.
  Low-input samples (needle biopsies, cfDNA) have inherently low library
  complexity.
- **Correction**: Aggressive QC filtering. Some tools (QDNAseq, CNVkit) have
  FFPE-specific modes. UMI-based duplicate tracking helps with low-input
  samples. Generally, RD-based CNV calling in FFPE requires larger bin sizes
  and detects only larger events.

---

## 2. Sequencing Chemistry Biases

### 2.1 Flow Cell Position Effects

- **What**: Coverage varies systematically by position on the flow cell surface.
  Tiles, lanes, and flow cell edges can show distinct coverage patterns.
- **Magnitude**: Usually small (~1-3% of variance) in modern instruments, but
  can be significant in older or poorly maintained sequencers.
- **Spatial scale**: Affects all reads from a given tile/lane uniformly, so the
  genomic signature depends on which reads land where -- effectively random
  at the genomic coordinate level but creates inter-sample variance.
- **Correction**: Not corrected explicitly. Averaging across many tiles/lanes
  dilutes the effect. Multi-lane pooling helps.

### 2.2 Index Hopping (Multiplexing Crosstalk)

- **What**: On patterned flow cells (NovaSeq, NovaSeqX), index sequences can
  swap between clusters during ExAmp chemistry, assigning reads to the wrong
  sample.
- **Magnitude**: 0.1-2% of reads may be mis-assigned. Creates low-level
  contamination between multiplexed samples.
- **Spatial scale**: Random genomic positions. Creates a uniform background
  noise floor.
- **Correction**: Unique dual indexing (UDI) reduces hopping to <0.01%.
  Computational detection via unexpected alleles at known SNP sites.

### 2.3 Sequencing Cycle / Base-Calling Errors

- **What**: Error rates increase toward the 3' end of reads and in
  homopolymer runs, potentially affecting mapping and coverage uniformity.
- **Magnitude**: Minor direct impact on RD (~0.5-1% of variance). Mainly
  relevant through its effect on mapping quality.
- **Spatial scale**: Affects mapping at homopolymeric and low-complexity loci.
- **Correction**: Base quality recalibration (BQSR) addresses systematic errors
  in quality scores. Most RD callers are robust to base-level errors since they
  count reads, not bases.

---

## 3. Target Capture Biases (WES / Panel-Specific)

### 3.1 Probe Hybridization Efficiency

- **What**: Each capture probe has a different hybridization efficiency
  determined by its sequence composition, length, Tm, and secondary structure.
  This creates a per-target "capture fingerprint."
- **Magnitude**: The dominant source of RD variance in WES/panel data. Coverage
  across targets routinely varies 10-100x even at constant copy number.
  Accounts for ~10-30% of total variance.
- **Spatial scale**: Per-target (per-exon). Adjacent exons of the same gene can
  have dramatically different coverage.
- **Root cause**: Hybridization kinetics (Tm, GC, secondary structure, probe
  length, target accessibility). Some exons are simply hard to capture.
- **Key property**: Highly reproducible across samples processed with the same
  capture kit. This reproducibility is what makes PoN and PCA normalization
  so effective for WES -- the per-target bias is systematic and removable.
- **Correction**: PoN subtraction, PCA, or reference-sample normalization.
  All exploit the cross-sample reproducibility of per-target efficiency.

### 3.2 Capture Bait Decay / Edge Effects

- **What**: Coverage drops sharply at the edges of capture targets (probe
  tiling gaps) and at targets with degraded or old baits.
- **Magnitude**: ~5-15% of targets may have significantly reduced capture
  efficiency, reducing effective coverage.
- **Spatial scale**: Per-target, concentrated at target boundaries.
- **Root cause**: Probe tiling design, bait synthesis quality, probe
  degradation during storage.
- **Correction**: Target interval padding (extend intervals by 50-100bp).
  PoN implicitly captures per-target efficiency. Exclude targets with
  consistently low coverage across the PoN.

### 3.3 Capture Kit Version Differences

- **What**: Different versions of the same capture kit (or different kits
  entirely) have different probe sets, creating batch-like effects.
- **Magnitude**: Can be severe. Mixing samples from different capture kits
  without correction introduces large systematic differences that mimic
  CNVs.
- **Spatial scale**: Per-target. Targets present in one kit but absent from
  another will show complete dropout.
- **Correction**: Never mix capture kits without explicit cross-kit
  normalization. PoN must be built from the same kit. PCA can partially
  compensate if both kits are represented in the cohort, but this is fragile.

---

## 4. Alignment and Mapping Biases

### 4.1 Multi-Mapping / Mappability

- **What**: Reads originating from repetitive or duplicated sequences cannot
  be uniquely placed. Aligners either assign them to one location with low
  MAPQ, distribute them, or discard them. This makes depth unreliable in
  low-mappability regions.
- **Magnitude**: ~5-15% of the human genome has mappability <0.5 for 150bp
  reads. In these regions, RD is systematically biased (often inflated for
  segmental duplications, deflated for collapsed repeats).
- **Spatial scale**: Kilobases to megabases, concentrated in segmental
  duplications, pericentromeric regions, subtelomeric regions, and gene
  families with recent paralogs.
- **Root cause**: The reference genome contains sequences that are so similar
  that short reads cannot distinguish between them. Segmental duplications
  (>1kb, >90% identity) cover ~5% of the human genome.
- **Correction**: Mappability filtering (exclude bins with mappability below
  threshold), MAPQ filtering (require MAPQ >=20), or mappability-adjusted
  depth (Control-FREEC divides depth by mappability score). T2T-CHM13
  reference resolves many previously collapsed regions.

### 4.2 Reference Genome Bias

- **What**: The reference genome is an incomplete, haploid representation of
  human sequence diversity. Regions where the reference is structurally
  different from the sample create systematic mapping artifacts.
- **Magnitude**: Variable. For common structural variants (polymorphic
  inversions, indels, copy number polymorphisms), the reference allele
  produces baseline depth artifacts in all samples carrying the non-reference
  allele.
- **Specific manifestations**:
  - **Collapsed duplications**: Paralogous sequences merged in the reference
    attract reads from all copies, inflating depth.
  - **Reference gaps/Ns**: Unmappable by definition.
  - **ALT haplotypes** (GRCh38): Reads from ALT haplotypes may map to the
    primary assembly with reduced MAPQ, deflating depth at the true locus.
  - **Population-specific insertions**: Sequences absent from the reference
    produce soft-clipped reads that reduce effective depth at flanking regions.
- **Correction**: Use the most complete reference available (T2T-CHM13 for
  analyses where it's appropriate). GRCh38+ALTs with alt-aware alignment
  (bwa-mem ALT pipeline) helps for known haplotypes. Pangenome references
  (HPRC) are an emerging solution.

### 4.3 MAPQ Assignment and Read Filtering

- **What**: Aligner-assigned mapping quality (MAPQ) is used by most tools to
  filter unreliably placed reads. The MAPQ model varies between aligners and
  is imperfect.
- **Magnitude**: MAPQ filtering removes ~3-8% of aligned reads in typical WGS.
  The threshold choice (0, 10, 20, 30) meaningfully affects depth profiles
  in multi-copy regions.
- **Spatial scale**: Concentrated in repetitive regions (same as mappability).
- **Root cause**: MAPQ is an estimate of placement confidence. BWA-MEM assigns
  MAPQ=0 to reads with equally good alternative alignments. Different aligners
  (BWA-MEM2, minimap2, Bowtie2) assign different MAPQs to the same reads.
- **Correction**: Use a consistent aligner and MAPQ threshold. Some tools are
  sensitive to the choice of aligner. Aligner version changes between batches
  can create subtle batch effects in depth profiles.

### 4.4 Soft-Clipping at Structural Variant Breakpoints

- **What**: Reads spanning SV breakpoints are partially aligned (soft-clipped),
  reducing the effective depth at and near breakpoints.
- **Magnitude**: Localised to within ~read-length of breakpoints. Reduces
  depth by the proportion of reads spanning the junction.
- **Spatial scale**: ~150-300bp around each breakpoint.
- **Correction**: Not corrected; this is actually useful signal for SR-based
  callers. For RD-only callers, the effect is absorbed into bin-level noise.

---

## 5. Genomic Context Biases

### 5.1 Replication Timing / Genomic Waves

- **What**: Long-range (100kb-1Mb) oscillations in coverage that persist after
  GC correction. Correlated with DNA replication timing -- early-replicating
  regions tend to have slightly higher coverage in dividing cell populations.
- **Magnitude**: ~5-10% of RD variance. Creates smooth "waves" of coverage
  variation across chromosomes.
- **Spatial scale**: Megabase-scale. Visible as broad undulations in log2-ratio
  profiles.
- **Root cause**: Cells in S-phase have already replicated early-replicating
  regions, contributing extra template copies. The effect scales with the
  fraction of S-phase cells in the DNA source. More pronounced in cell lines
  and tumors (high proliferative index) than in blood.
- **Correction**: Post-GC wave correction (QDNAseq, CNVkit). Smoothing-based
  approaches subtract the megabase-scale trend. Risk of overcorrecting real
  large CNVs.

### 5.2 Segmental Duplications and Paralogs

- **What**: Regions with recent segmental duplications (>1kb, >90% identity)
  cause multi-mapping of reads. Depth in these regions reflects the *total
  copy number across all paralogous copies*, not the copy number at a single
  locus.
- **Magnitude**: ~5% of the genome is in segmental duplications. These regions
  are enriched for CNV polymorphism and disease-relevant genes (e.g., CYP2D6,
  SMN1/SMN2, FCGR, STRC). RD in SDs can be inflated 2-10x above the
  expected diploid baseline.
- **Spatial scale**: 1kb to several Mb per duplicated block.
- **Root cause**: Evolutionary: recent duplications produce nearly identical
  sequences. Short reads (~150bp) cannot resolve copies that differ by fewer
  than a few SNPs per read.
- **Correction**: Specialised paralog-aware tools (e.g., Gauchian for CYP2D6,
  SMNCopyNumberCaller for SMN). WGS-depth CNV callers typically filter or
  mask SDs. Long-read sequencing is the definitive solution.

### 5.3 Centromeric and Pericentromeric Regions

- **What**: Centromeres consist of highly repetitive alpha-satellite arrays
  that are largely unmappable with short reads. Pericentromeric regions have
  reduced mappability and enrichment for segmental duplications.
- **Magnitude**: ~5-10% of the genome is effectively invisible to short-read
  mapping. Pericentromeric regions produce unreliable depth estimates.
- **Spatial scale**: Megabases around each centromere.
- **Correction**: Blacklist exclusion. T2T-CHM13 includes centromeric sequence
  but even with T2T, short-read mapping in satellite arrays is unreliable.

### 5.4 Telomeric and Subtelomeric Regions

- **What**: Telomeric repeats (TTAGGG)n are unmappable. Subtelomeric regions
  are enriched for segmental duplications and structural polymorphisms.
- **Magnitude**: Modest in total genome fraction but clinically relevant
  (subtelomeric CNVs are a common cause of intellectual disability).
- **Correction**: Blacklist exclusion for telomeric repeats. Subtelomeric
  regions require careful per-region validation.

### 5.5 Sex Chromosomes

- **What**: Males have one X and one Y; females have two X and no Y.
  Pseudoautosomal regions (PAR1/PAR2) recombine and have diploid depth in
  both sexes. The Y has massive ampliconic and palindromic regions.
- **Magnitude**: Baseline copy number differs by sex, and the X inactivation
  pattern can affect coverage in some assays.
- **Correction**: Sex-matched reference or explicit sex-chromosome models.
  Most tools normalise autosomes and sex chromosomes separately. Sample sex
  must be correctly inferred (or provided).

---

## 6. Batch and Process Effects

### 6.1 Sequencing Batch Effects

- **What**: Samples processed at different times, on different instruments,
  or with different reagent lots show systematic coverage differences even
  when all other protocols are identical.
- **Magnitude**: ~5-20% of RD variance. Can dominate over GC bias if
  cross-batch comparisons are made without correction.
- **Spatial scale**: Affects entire depth profile -- not localised.
  Manifests as shifts in the GC-coverage curve shape, overall coverage
  level, and noise characteristics.
- **Root cause**: Reagent lot variation, instrument calibration drift,
  operator differences, environmental conditions. Even within the same
  sequencing center, month-to-month variation is measurable.
- **Correction**: PoN constructed from the same batch. PCA across a
  cohort that includes the batch. Reference sample selection (ExomeDepth)
  naturally adapts. Process controls (run the same reference sample in
  each batch) allow explicit batch correction.

### 6.2 Library Prep Protocol Variation

- **What**: Differences in library preparation protocol (PCR cycles,
  enzyme lot, fragmentation method, adapter ligation efficiency) create
  systematic coverage variation.
- **Magnitude**: Can be large, especially when switching between PCR-free
  and PCR-based protocols (PCR-free libraries have dramatically reduced
  GC bias).
- **Spatial scale**: Genome-wide, modulated by sequence content.
- **Correction**: Never mix PCR-free and PCR-based libraries without
  explicit correction. PoN and PCA can compensate if the cohort is
  internally consistent.

### 6.3 Aligner and Software Version Changes

- **What**: Updating the aligner (BWA-MEM → BWA-MEM2), the reference genome
  build, duplicate marking tools, or even software versions can shift depth
  profiles enough to create false CNV calls in cross-batch comparisons.
- **Magnitude**: Usually small (~1-3%) for minor version updates. Can be
  substantial for aligner changes or reference build switches.
- **Correction**: Reprocess all samples with the same pipeline version when
  performing cohort-level CNV calling. Version-lock pipeline components.

---

## 7. Stochastic / Fundamental Noise

### 7.1 Poisson Sampling Noise

- **What**: Read counts in a bin follow (approximately) a Poisson distribution.
  The coefficient of variation (CV) scales as 1/√n, where n is the expected
  count.
- **Magnitude**: At 30x WGS with 1kb bins, expected ~200 reads/bin →
  CV ≈ 7%. At 1x shallow WGS with 10kb bins → ~67 reads/bin → CV ≈ 12%.
  For WES at 100x on a 200bp exon → ~133 reads → CV ≈ 8.7%.
- **Spatial scale**: Independent per bin (no spatial correlation).
- **Key property**: This is the irreducible noise floor. All other biases
  sit on top of it. Increasing coverage is the only way to reduce it.
- **Correction**: Cannot be corrected, only modeled. Tools use statistical
  models (Poisson, negative binomial, Gaussian on log-ratios) that account
  for expected sampling variance. Larger bins reduce variance but sacrifice
  resolution -- this is the fundamental resolution-sensitivity tradeoff.

### 7.2 Overdispersion (Negative Binomial)

- **What**: Real sequencing data is more variable than Poisson predicts.
  The variance-to-mean ratio is >1 (overdispersed), often 1.5-5x.
- **Magnitude**: The "extra" variance beyond Poisson is substantial and
  represents the aggregate of all the small, uncorrectable biases listed
  above.
- **Root cause**: Residual GC effects, fragment sampling, local chromatin
  structure, and other sources that are too fine-grained to model explicitly
  contribute additional variance.
- **Correction**: Use negative binomial or beta-binomial models instead of
  Poisson. ExomeDepth uses beta-binomial specifically to handle
  overdispersion. GATK gCNV uses negative binomial.

### 7.3 Library Complexity

- **What**: The total number of unique molecules in the library sets an
  upper bound on useful depth. Beyond this point, additional sequencing
  produces only duplicates, not new coverage.
- **Magnitude**: Depends on input DNA amount and library prep efficiency.
  Clinical samples with low input (<10ng) may have complexity of only
  100-500M unique fragments, limiting usable depth to ~10-15x even if
  sequenced to 60x raw.
- **Correction**: UMI-based deduplication gives exact counts. Without UMIs,
  duplicate marking is heuristic (coordinate-based). Monitor duplication
  rate as a QC metric -- rates >30% indicate low complexity.

---

## 8. Sample-Level Biological Variation

### 8.1 Tumor Purity (Somatic Only)

- **What**: Tumor samples are admixed with normal cells. A deletion in 100%
  of tumor cells at 50% purity produces only a 25% depth drop, not 50%.
- **Magnitude**: Dominates somatic CNV detection difficulty. Purity ranges
  from <10% (liquid biopsies) to >90% (cell lines). At 30% purity, a
  heterozygous deletion produces only a 15% signal -- well within the noise
  of many methods.
- **Correction**: Joint purity/ploidy estimation (FACETS, ASCAT, PureCN,
  PURPLE, ABSOLUTE). Requires BAF from heterozygous SNPs to break the
  purity-ploidy degeneracy.

### 8.2 Tumor Ploidy (Somatic Only)

- **What**: Many tumors are not diploid. Whole-genome duplication (tetraploidy)
  shifts the baseline copy number, so the "normal" depth represents CN=4, not
  CN=2. A "deletion" back to CN=2 looks like a loss from the baseline.
- **Magnitude**: ~30-40% of solid tumors have undergone whole-genome
  duplication. Without ploidy estimation, all copy number calls are
  systematically wrong.
- **Correction**: Same as purity -- joint estimation using RD + BAF. Grid
  search (Sequenza) or probabilistic models (ABSOLUTE, PureCN) over
  purity-ploidy space.

### 8.3 Subclonal Copy Number / Mosaicism

- **What**: Not all cells in the sample carry the same CNVs. Subclonal events
  produce fractional depth changes that may fall below detection thresholds.
  In germline, somatic mosaicism produces the same effect.
- **Magnitude**: A subclonal deletion present in 20% of cells at 100% purity
  produces only a 10% depth drop. Most callers require >20-30% cell fraction
  for reliable detection.
- **Correction**: Battenberg explicitly models subclonal CN. Most other tools
  call only clonal events. Increasing depth and using larger bin sizes
  improves sensitivity to subclonal events.

### 8.4 Sample Contamination

- **What**: Cross-contamination between samples (during library prep,
  sequencing, or demultiplexing) adds a background of foreign reads that
  dilute the true depth signal.
- **Magnitude**: Even 1-3% contamination can create false heterozygous SNP
  calls and subtly shift allele frequencies, confounding BAF-based methods.
  Higher contamination directly affects RD accuracy.
- **Correction**: QC tools (VerifyBamID, Conpair, ContEst) detect
  contamination. Contaminated samples should be excluded or, if the
  contamination source is known, computationally decontaminated.

---

## Summary: Variance Sources Ranked by Magnitude

```
Source                      Typical % of RD    Spatial       Correctable?
                            Variance           Scale
────────────────────────────────────────────────────────────────────────
Capture probe efficiency    10-30% (WES)       Per-target    Yes (PoN/PCA)
GC content bias             20-40% (WGS)       100bp-10kb    Mostly (LOESS)
Batch effects               5-20%              Genome-wide   Yes (PoN/PCA)
Poisson sampling noise      5-15%              Per-bin       No (fundamental)
Mappability / multi-map     5-15%              1kb-1Mb       Partially (filter)
Replication timing waves    5-10%              100kb-1Mb     Mostly (wave corr.)
Tumor purity (somatic)      0-50%              Genome-wide   Yes (purity est.)
PCR duplication residual    2-10%              GC-correlated Mostly (dedup)
Fragment length effects     2-5%               ~150bp-1kb    Partially
Reference genome bias       1-5%               Variable      Partially (T2T)
DNA quality (FFPE)          0-30%              Genome-wide   Partially (QC)
Aligner/software version    1-3%               Genome-wide   Yes (reprocess)
Flow cell position          1-3%               Random        Averaged out
Sequencing cycle errors     0.5-1%             Homopolymers  Minimal impact
Index hopping               0.1-2%             Random        Yes (UDI)
```

---

## Implications for Caller Design

### The Normalization Hierarchy

Sources of variance are not independent -- they interact and compound.
Effective normalization follows a hierarchy from largest to smallest effect:

```
1. Remove catastrophic outliers     (blacklist, mappability filter)
2. Remove dominant systematic bias  (GC correction via LOESS)
3. Remove residual systematic bias  (wave correction, PoN, PCA)
4. Model remaining overdispersion   (negative binomial / beta-binomial)
5. Segment the cleaned signal       (CBS, HMM, BIC)
6. Classify segments                (integer CN states, with or without
                                     purity/ploidy adjustment)
```

### Why WES Is Harder Than WGS

WES adds a massive additional variance source (probe hybridization efficiency)
that doesn't exist in WGS. This is why WES tools universally require either a
cohort (PCA) or a matched reference set (PoN / reference selection) -- single-
sample WES CNV calling without a reference baseline is essentially impossible.

### Why Duplications Are Harder Than Deletions

A heterozygous deletion drops depth by 50% from baseline -- a large, clean
signal. A heterozygous duplication increases depth by 50% -- the same absolute
magnitude, but from a higher baseline, so the *relative* signal-to-noise ratio
is slightly worse. More importantly, many of the bias sources listed above
(GC, capture efficiency, mappability) can produce *increases* in coverage that
mimic duplications, while coverage *decreases* are harder to produce
artifactually. The false positive rate for duplications is therefore inherently
higher.

### The Resolution-Sensitivity Tradeoff

Poisson sampling noise is the irreducible floor. Larger bins average over more
reads, reducing noise, but sacrificing spatial resolution. The optimal bin size
depends on coverage depth and the minimum CNV size of interest:

```
Coverage    Bin Size    Reads/Bin    CV (Poisson)    Min Detectable CNV
────────────────────────────────────────────────────────────────────────
0.5x        50kb        ~167         ~7.7%           ~500kb
1x          10kb        ~67          ~12.2%          ~100kb
5x          1kb         ~33          ~17.4%          ~10kb
30x         1kb         ~200         ~7.1%           ~3kb
30x         100bp       ~20          ~22.4%          ~1kb (noisy)
100x        100bp       ~67          ~12.2%          ~500bp
```

---

*Companion to: cnv_calling_mind_map.md and rd_normalization_approaches.md.
This document catalogs the sources of variance; the normalization doc covers
how each is addressed; the mind map covers tool selection and ensemble
strategies.*
