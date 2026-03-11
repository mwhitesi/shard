# K-mer Encoding Schemes for ML: Locality-Preserving Representations

## 1. Locality-Sensitive Orderings of K-mers

### 1.1 De Bruijn Graph Traversals

A de Bruijn graph of order k is a directed graph where nodes are k-mers and edges connect k-mers that overlap by k-1 characters. An Eulerian or Hamiltonian path through this graph visits k-mers in an order where consecutive entries differ by a single character shift --- inherently preserving sequence-context similarity.

- **Similarity preserved**: Sequence-context overlap (k-1 shared bases between consecutive k-mers).
- **Used in genomics ML**: Kmer-Node2Vec (Ren et al., 2023) learns embeddings by performing random walks on the k-mer co-occurrence graph (a weighted de Bruijn-like structure), then feeds walks to a Word2Vec-style skip-gram model. This is 29x faster than dna2vec for training on multi-GB datasets while achieving comparable downstream accuracy.
- **Key tools**: Genome assemblers (SPAdes, MEGAHIT) use de Bruijn graphs as their core data structure, though primarily for assembly rather than ML encoding.

### 1.2 Locality-Preserving Minimal Perfect Hashing (LPHash)

LPHash (Pibiri, 2023) constructs a minimal perfect hash function for a k-mer set that maps consecutive k-mers (sharing k-1 overlap) to consecutive integer identifiers. It exploits the fact that successive k-mers in a sequence share k-1 bases.

- **Similarity preserved**: Sequence contiguity --- k-mers adjacent in the genome get adjacent hash values.
- **Benefits**: (1) Satellite data associated with k-mers (e.g., counts) compresses much better because correlated values cluster together. (2) Query throughput improves via cache locality when streaming k-mers from a sequence.
- **Used in genomics ML**: Primarily a data structure contribution, but directly applicable to any ML pipeline that needs to index into k-mer feature vectors efficiently.
- **Key paper**: Pibiri (2023), *Bioinformatics* 39(Suppl 1):i534. Tool: [github.com/jermp/lphash](https://github.com/jermp/lphash).

### 1.3 Hilbert Curves on K-mer Space

Hilbert curves are space-filling curves that map a 1D index to a 2D (or higher-dimensional) grid while preserving locality much better than row-major or Z-order mappings. In genomic ML, DNA sequences are mapped onto 2D images via a Hilbert curve, then processed by CNNs.

- **Similarity preserved**: Sequential proximity in the 1D sequence translates to spatial proximity in the 2D image. Gray codes are used internally to encode nucleotides such that successive integer values differ by a single bit, minimising Hamming distance between adjacent Hilbert curve positions.
- **Used in genomics ML**: Yes --- CNN models using Hilbert curve representations for enhancer prediction (Gao & Ren, 2019, bioRxiv:552141) and general DNA classification. The approach converts each k-mer to a one-hot vector, then lays out k-mers along a Hilbert curve to form an image, enabling 2D convolutions to capture both local and long-range sequence patterns.
- **Key paper**: Gao & Ren (2019), "CNN Model With Hilbert Curve Representation of DNA Sequence For Enhancer Prediction". Also: Basu et al. (2023), "Hilbert Curve Based Molecular Sequence Analysis".

### 1.4 Gray Codes over Nucleotide Alphabets

A Gray code orders binary strings so that consecutive entries differ in exactly one bit. Applied to the 4-letter nucleotide alphabet (encoded as 2-bit per base), a Gray code ordering of all k-mers would ensure consecutive k-mers in the ordering differ by a single nucleotide substitution at one position.

- **Similarity preserved**: Hamming distance / single-nucleotide edit distance.
- **Used in genomics ML**: Not widely as a standalone k-mer ordering scheme, but Gray code encodings of nucleotides are used as a component within Hilbert curve DNA image construction (see 1.3) to maintain smoothness in the mapping. The concept also underlies some genetic algorithm representations of DNA parameters.
- **Practical note**: A pure Gray code ordering of k-mers does not preserve *sequence context* (unlike de Bruijn traversals), only *point-mutation proximity*.

### 1.5 Chaos Game Representation (CGR)

CGR (Jeffrey, 1990) is an iterative mapping that places each nucleotide of a DNA sequence as a point in a unit square, where each corner represents one of {A, C, G, T}. Each successive nucleotide is plotted at the midpoint between the current position and the corner of the new nucleotide.

- **Similarity preserved**: The resulting image is a fractal that is *exactly equivalent* to the full k-mer frequency spectrum for all k. Subdividing the CGR image into a 2^k x 2^k grid yields a Frequency CGR (FCGR) matrix that is the count of each k-mer. Crucially, k-mers sharing a common prefix or suffix map to nearby regions of the image.
- **Used in genomics ML**: Extensively --- FCGR images are used as input to CNNs for taxonomic classification, phylogenetic analysis, and genome comparison. FCGR inherently preserves a notion of k-mer similarity where k-mers with shared subsequences are spatially proximate.
- **Key paper**: Jeffrey (1990); Deschavanne et al. (1999), "Genomic signature"; review by Loesch & Hoinka (2022), *Computational and Structural Biotechnology Journal*.

---

## 2. K-mer Embedding Methods

### 2.1 dna2vec

dna2vec (Ng, 2017) applies Word2Vec (skip-gram) to DNA, treating overlapping k-mers from genomic sequences as "words" in a "sentence". The model learns dense vector representations for variable-length k-mers.

- **Similarity preserved**: Sequence-context similarity. K-mers that appear in similar genomic contexts get similar vectors. A key property: the cosine similarity between dna2vec vectors correlates with the Needleman-Wunsch alignment score between the k-mer sequences, providing a consistent notion of sequence similarity.
- **Used in genomics ML**: Yes --- used for downstream tasks including sequence classification, promoter/enhancer prediction, and splice site identification.
- **Key paper**: Ng (2017), arXiv:1701.06279.

### 2.2 kmer2vec

kmer2vec (Ren et al., 2022) also uses Word2Vec but focuses on alignment-free sequence comparison. K-mers from entire genomes are embedded, and a sequence is represented as an aggregation (e.g., average) of its constituent k-mer vectors. Pairwise sequence distances are then computed in the embedding space.

- **Similarity preserved**: Contextual co-occurrence in genomes --- captures evolutionary/functional relationships.
- **Used in genomics ML**: Phylogenetic tree construction and species clustering. Demonstrated to be faster than multiple sequence alignment methods while producing more accurate phylogenetic relationships than conventional k-mer alignment-free approaches.
- **Key paper**: Ren et al. (2022), *Journal of Computational Biology*.

### 2.3 Kmer-Node2Vec

Kmer-Node2Vec (2023) constructs a k-mer co-occurrence graph and applies the Node2Vec algorithm (biased random walks + skip-gram) to learn k-mer embeddings from the graph structure.

- **Similarity preserved**: Graph-topological similarity of k-mer co-occurrence patterns, which captures both sequence context and frequency-based relationships.
- **Used in genomics ML**: Yes --- benchmarked on DNA sequence classification tasks. 29x faster training than dna2vec on a 4 GB dataset with comparable accuracy.
- **Key paper**: bioRxiv:2022.08.30.505832.

### 2.4 BPE-style Tokenisation (DNABERT-2)

Byte Pair Encoding (BPE) iteratively merges the most frequently co-occurring adjacent tokens to build a vocabulary of variable-length subword units. Applied to DNA, BPE discovers recurrent motifs of varying length.

- **Similarity preserved**: Frequency co-occurrence at the corpus level --- frequent genomic motifs become single tokens. This is distinct from edit-distance or context similarity; it captures *statistical recurrence patterns*.
- **Used in genomics ML**: Central to DNABERT-2 (Zhou et al., 2023, ICLR 2024), which replaced DNABERT's fixed k-mer tokenisation with BPE. DNABERT-2 achieves comparable performance to state-of-the-art models with 21x fewer parameters and ~92x less GPU pre-training time. BPE reduces sequence length ~5x compared to single-nucleotide tokenisation.
- **Limitations**: BPE does not efficiently learn known biological motifs (motifs are important but not necessarily the most frequent substrings). The optimal tokeniser for DNA remains an open question.
- **Key papers**: DNABERT (Ji et al., 2021, *Bioinformatics*); DNABERT-2 (Zhou et al., 2023, arXiv:2306.15006). Also: Nucleotide Transformer (Dalla-Torre et al., 2023) scales single-nucleotide tokenisation to very large models.

### 2.5 Hybrid Tokenisation

Recent work (2025) proposes combining 6-mer tokenisation with BPE to capture both short motifs (via k-mers) and longer recurrent patterns (via BPE), creating a balanced and context-aware vocabulary.

- **Key paper**: arXiv:2507.18570.

---

## 3. Frequency-Domain and Spectral Encodings

### 3.1 K-mer Spectrum as Feature Vector

The simplest spectral encoding: count the frequency of every possible k-mer in a sequence, yielding a vector in R^{4^k}. This is the basis of many alignment-free comparison methods.

- **Similarity preserved**: Composition similarity (Jaccard, cosine, or Euclidean distance on frequency vectors correlates with overall sequence similarity).
- **Used in genomics ML**: Ubiquitous --- k-mer frequency vectors are used for LTR retrotransposon classification (F1=95%), metagenomic binning, species identification, and more.
- **Key paper**: Chor et al. (2009), "Genomic DNA k-mer spectra: models and modalities", *Genome Biology*.

### 3.2 KPop: Correspondence Analysis on K-mer Spectra

KPop (Utro et al., 2025) takes the full k-mer frequency spectrum and applies Correspondence Analysis (CA) --- a spectral/SVD-based dimensionality reduction technique --- to produce a "twisted" low-dimensional embedding optimised for a given dataset.

- **Similarity preserved**: The CA transformation maximises the chi-squared distance between samples, which captures compositional differences in k-mer usage. It preserves both species-level and sub-species-level separation even when genomic diversity is low.
- **Advantages over MinHash**: KPop maps sequences to coordinates in a continuous space (enabling ML directly), whereas MinHash methods only produce pairwise distances with lower resolution.
- **Used in genomics ML**: Yes --- microbial genome classification and comparison at scale.
- **Key paper**: Utro et al. (2025), *Genome Biology*.

### 3.3 The Folded K-Spectrum Kernel

The folded k-spectrum kernel (Hooghe et al., 2017) extends the basic k-mer spectrum by "folding" in gapped k-mer features. For each k-mer, it also computes features for all possible gapped variants (where subsets of positions are masked), effectively incorporating variable-length gaps from 0 to k-1.

- **Similarity preserved**: Both contiguous and gapped subsequence composition, capturing complex nucleotide dependencies including distal base interactions within binding sites.
- **Used in genomics ML**: Transcription factor binding site prediction.
- **Key paper**: Hooghe et al. (2017), *PLOS ONE*.

### 3.4 Numerical / Signal-Processing Encodings

DNA sequences can be converted to numerical signals (e.g., via Voss, Z-curve, electron-ion interaction potential, or CGR mappings) and then analysed with Fourier transforms, wavelet transforms, or other spectral methods. These approaches treat k-mer composition as a signal in the frequency domain.

- **Similarity preserved**: Periodicity and spectral characteristics of the sequence.
- **Used in genomics ML**: Gene prediction (detection of period-3 signal in coding regions), repeat detection, and genome-wide structural analysis.

---

## 4. Minimizer and Syncmer Schemes: Locality Properties

### 4.1 Minimizers

A minimizer of a window of w consecutive k-mers is the k-mer with the smallest hash value in that window. Minimizers subsample k-mers such that overlapping reads sharing a region will select the same representative k-mer.

- **Locality properties**: Minimizers are *context-dependent* --- a k-mer's selection status depends on its flanking sequence. This means a mutation in flanking sequence can change which k-mer is selected, even if the k-mer itself is unchanged. Consecutive k-mers in a genome tend to share the same minimizer (this is the whole point --- it partitions the sequence into runs), which provides a form of locality. The ordering within the minimizer hash function (lexicographic, random, or designed) determines how well locality is preserved in the hash space.
- **Used in genomics ML**: Kraken2 uses minimizers for taxonomic classification; minimap2 uses them for read mapping and seeding.
- **Key paper**: Roberts et al. (2004); Schleimer et al. (2003), WINNOWING.

### 4.2 Syncmers

Syncmers (Edgar, 2021) select a k-mer based solely on the position of the smallest s-mer *within* the k-mer itself (e.g., open syncmer: smallest s-mer at the start; closed syncmer: at start or end).

- **Locality properties**: A syncmer's selection is determined by its own sequence alone --- it is *context-independent*. If a k-mer is selected in one sequence, it will also be selected in any other sequence containing it. This "synchronisation" property means syncmers are more robust to mutations than minimizers. Experiments show syncmers simultaneously achieve lower density (fewer selected k-mers) and higher conservation (more selected k-mers are shared between related sequences) than minimizers.
- **Similarity preserved**: Sequence identity / conservation --- syncmers are more sensitive for detecting conserved k-mers across sequences.
- **Key paper**: Edgar (2021), *PeerJ*.

### 4.3 Mod-Minimizers

Mod-minimizers (Groot Koerkamp & Pibiri, 2024) combine minimizer and modular sampling ideas: select a k-mer if its minimizer hash modulo t equals a target value. This achieves lower density than standard minimizers while maintaining locality (consecutive k-mers sharing a minimizer are grouped).

- **Key paper**: Groot Koerkamp & Pibiri (2024), WABI.

### 4.4 Do These Schemes Preserve Locality for ML?

Minimizers and syncmers are primarily **subsampling** schemes rather than **encoding** schemes. They reduce the number of k-mers that need to be stored or compared. For ML purposes:
- They preserve **sequence-level locality** (nearby positions share representatives).
- They do **not** directly provide an embedding or vector representation.
- They are typically combined with a downstream encoding (hash table lookup, MinHash sketch, or frequency counting) for use in ML pipelines.

---

## 5. Kernel Methods and Distance Metrics on K-mer Space

### 5.1 Spectrum Kernel

The spectrum kernel (Leslie et al., 2002) defines the similarity between two sequences as the inner product of their k-mer frequency vectors. Equivalent to counting shared k-mers.

- **Similarity preserved**: Exact k-mer composition overlap.
- **Used in genomics ML**: Protein family classification (original application); later adapted for DNA.
- **Key paper**: Leslie et al. (2002), *Pacific Symposium on Biocomputing*.

### 5.2 Mismatch Kernel

The (k,m)-mismatch kernel (Leslie et al., 2004) extends the spectrum kernel by counting not just exact k-mer matches but also k-mers within Hamming distance m.

- **Similarity preserved**: Approximate sequence similarity (tolerates point mutations).
- **Used in genomics ML**: Protein classification, TF binding site prediction.
- **Key paper**: Leslie et al. (2004), *JMLR*.

### 5.3 Gapped K-mer Kernel (gkm-SVM / LS-GKM)

The gapped k-mer kernel counts the number of shared *gapped* k-mers (subsequences with allowed gaps/wildcards at certain positions) between two sequences. This generalises both the spectrum and mismatch kernels.

- **Similarity preserved**: Subsequence similarity with tolerance for insertions/deletions at specific positions. Captures distal nucleotide dependencies within regulatory elements.
- **Used in genomics ML**: Extensively in regulatory genomics. gkm-SVM and its scalable successor LS-GKM (Lee, 2016) are used for TF binding prediction, enhancer prediction, chromatin accessibility prediction, and variant impact scoring. GkmExplain (Shrikumar et al., 2019) provides interpretability for these models.
- **Key tools**: [gkmSVM](https://github.com/Dongwon-Lee/lsgkm), gkmQC, GkmExplain.
- **Key papers**: Ghandi et al. (2014), *PLOS Computational Biology*; Lee (2016), *Bioinformatics*; Shrikumar et al. (2019), *Bioinformatics*.

### 5.4 GaKCo and FastSK

These are computational accelerations of gapped k-mer kernels:
- **GaKCo** (Singh et al., 2017): Uses associative arrays and cumulative counting to avoid the O(|Sigma|^M) trie-based bottleneck.
- **FastSK** (Blakely et al., 2020): Decomposes the kernel into independent counting operations over mismatch positions, enabling parallelisation.

- **Key papers**: Singh et al. (2017), arXiv:1704.07468; Blakely et al. (2020), *Bioinformatics*.

### 5.5 Edit-Distance and Alignment-Based Kernels

Some kernels directly use edit distance or local alignment scores (e.g., the local alignment kernel of Vert et al., 2004). These are more expressive but computationally expensive.

- **Similarity preserved**: Full edit distance (substitutions, insertions, deletions).

---

## 6. Hashing and Ordering in Kraken, Mash, and Sourmash

### 6.1 Mash (MinHash)

Mash (Ondov et al., 2016) applies the MinHash locality-sensitive hashing scheme to k-mer sets. It hashes every k-mer in a genome using a single hash function (MurmurHash3), retains the s smallest hash values (the "sketch"), and estimates Jaccard similarity between two genomes by comparing their sketches.

- **Locality properties**: MinHash is a form of *locality-sensitive hashing* for Jaccard similarity --- two k-mer sets with high Jaccard similarity are likely to share minimum hash values. However, this operates at the *set level*, not at the *individual k-mer level*. There is no meaningful locality in how individual k-mers map to hash values (MurmurHash is designed to scatter uniformly).
- **Similarity preserved**: Jaccard similarity of k-mer sets, which Mash converts to an estimated mutation distance (Mash distance approximates 1 - ANI).
- **Used in genomics ML**: Primarily for fast approximate distance computation; Mash distances are used as input to clustering/classification.
- **Key paper**: Ondov et al. (2016), *Genome Biology*.

### 6.2 Sourmash (MinHash / FracMinHash)

Sourmash (Brown & Irber, 2016) uses a similar approach to Mash but supports FracMinHash (also called "modulo hash" or "scaled MinHash"), where a k-mer is retained if hash(kmer) < H/s for some fraction s. This produces sketches whose size scales with genome size, enabling containment estimation.

- **Locality properties**: Same as Mash --- LSH for set-level Jaccard/containment similarity. No individual k-mer locality in hash space.
- **Hash function**: MurmurHash (same as Mash).
- **Used in genomics ML**: Metagenome search, contamination screening, taxonomic classification.
- **Key tool**: [sourmash.readthedocs.io](https://sourmash.readthedocs.io).

### 6.3 Dashing 2 (SetSketch + LSH)

Dashing 2 (Baker & Langmead, 2023) replaces HyperLogLog sketches (used in Dashing 1) with SetSketch, a data structure that supports multiplicity-aware sketching (via ProbMinHash). It integrates locality-sensitive hashing to scale all-pairs comparisons to millions of sequences.

- **Locality properties**: The LSH component groups similar sketches (i.e., similar genomes) into the same hash buckets, enabling fast nearest-neighbour search. This is *sketch-level* LSH, not individual-k-mer-level locality.
- **Similarity preserved**: Jaccard coefficient and average nucleotide identity (ANI).
- **Key paper**: Baker & Langmead (2023), *Genome Research*.

### 6.4 Kraken 2 (Minimizer + Compact Hash Table)

Kraken 2 (Wood et al., 2019) does not store full k-mers. Instead, it extracts the minimizer (smallest l-mer) from each k-mer and stores the minimizer in a compact hash table (CHT) that uses linear probing.

- **Locality properties**: (1) The minimizer scheme provides *sequence-level locality* --- consecutive k-mers in a read tend to share the same minimizer, so they map to the same or nearby hash table entries. (2) The CHT with linear probing provides *memory-access locality* --- probing is sequential in memory, yielding good cache performance. (3) However, there is no *k-mer-similarity locality* in the hash space --- the CHT is a lookup structure, not an embedding.
- **Spaced seeds**: Kraken 2 uses spaced seeds (masks with gaps) when computing minimizers, which improves sensitivity to approximate matches.
- **Key paper**: Wood et al. (2019), *Genome Biology*.

---

## Summary Table

| Approach | Similarity Preserved | Genomics ML Use | Key Reference |
|---|---|---|---|
| De Bruijn traversal | Sequence context (k-1 overlap) | Kmer-Node2Vec embeddings | Ren et al. (2023) |
| LPHash | Sequence contiguity | Efficient k-mer indexing | Pibiri (2023) |
| Hilbert curve encoding | Sequential proximity | CNN on DNA images | Gao & Ren (2019) |
| CGR / FCGR | Full k-mer spectrum (prefix/suffix locality) | CNN classification, phylogenetics | Jeffrey (1990) |
| dna2vec | Sequence context (Word2Vec) | General sequence classification | Ng (2017) |
| kmer2vec | Contextual co-occurrence | Phylogenetics, species clustering | Ren et al. (2022) |
| BPE (DNABERT-2) | Statistical co-occurrence | Foundation model tokenisation | Zhou et al. (2023) |
| K-mer spectrum + CA (KPop) | Compositional chi-squared distance | Microbial genome comparison | Utro et al. (2025) |
| Spectrum kernel | Exact k-mer composition | Protein/DNA classification | Leslie et al. (2002) |
| Mismatch kernel | Hamming distance tolerance | TF binding, protein families | Leslie et al. (2004) |
| Gapped k-mer kernel (gkm-SVM) | Subsequence with gaps | Regulatory genomics (TF, chromatin) | Ghandi et al. (2014) |
| Minimizers | Sequence-level subsampling | Kraken2, minimap2 | Roberts et al. (2004) |
| Syncmers | Context-independent conservation | Read mapping, indexing | Edgar (2021) |
| MinHash (Mash/sourmash) | Set-level Jaccard similarity | Genome distance estimation | Ondov et al. (2016) |
| Dashing 2 (SetSketch + LSH) | Jaccard + ANI with multiplicities | All-pairs genome comparison | Baker & Langmead (2023) |

---

## Key Takeaways for ML Pipeline Design

1. **For preserving edit-distance locality**: Gray codes, mismatch kernels, and gapped k-mer kernels directly operate in edit-distance space. CGR/FCGR also places edit-distance-similar k-mers in spatial proximity.

2. **For preserving sequence-context locality**: dna2vec, kmer2vec, Kmer-Node2Vec, and de Bruijn traversals capture which k-mers co-occur in genomic sequences. BPE tokenisation captures statistical co-occurrence at a coarser granularity.

3. **For preserving compositional similarity**: K-mer spectra, KPop (with CA), and MinHash-based sketches capture overall genome composition.

4. **For efficient indexing with locality**: LPHash provides contiguity-preserving perfect hashing; Hilbert curves provide 2D spatial locality for CNN input.

5. **For kernel-based ML (SVMs)**: The gkm-SVM family (LS-GKM, GaKCo, FastSK) is the most mature and widely validated approach in regulatory genomics, with interpretability tools (GkmExplain).

6. **For foundation models / deep learning**: BPE tokenisation (DNABERT-2) or single-nucleotide tokenisation (Nucleotide Transformer) are the current standards. The choice of tokeniser significantly impacts performance and efficiency, and the optimal approach remains an open research question.
