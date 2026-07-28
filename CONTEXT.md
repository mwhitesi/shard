# shard

Copy-number variation (CNV) calling for targeted gene panels from aligned
short reads, using a per-locus **fragment tensor** normalised against a panel
of normals. This file is the project glossary: it fixes the language, not the
implementation.

## Language

**Panel**:
The fixed set of genomic regions a sequencing assay targets. One panel design
is shared by every sample in a cohort.
_Avoid_: kit, assay (use for the wet-lab product, not the region set)

**Target**:
A single region the panel is designed to capture (typically an exon).
_Avoid_: region, interval, amplicon

**Tile**:
A fixed-width sub-division of a target (and its flanks) — the unit along the
position axis of a fragment tensor. Tiles give a target internal spatial
resolution that a single per-target count discards.
_Avoid_: bin (reserve "bin" for histogram buckets within a tensor), window

**Flank**:
The margin on either side of a target that is tiled alongside it, to capture
the off-target-mate signal from fragments that straddle a capture boundary.

**Fragment**:
The single sequenced molecule behind a read pair, characterised by its start
position and insert size. The fragment — not the individual read — is the unit
of evidence.
_Avoid_: read (a read is one end; a fragment is the whole molecule), template

**Insert size**:
The full length of a fragment, inferred from the mapped positions of its two
read ends.
_Avoid_: TLEN, fragment length (use "insert size" everywhere)

**Off-target mate**:
The read end of an on-target fragment whose mate maps outside any target —
the signal SavvyCNV/CNVkit exploit, here carried as a tensor channel rather
than a separate count.

**Fragment tensor**:
The per-tile representation of fragment geometry: a histogram over
(position × insert size) plus auxiliary channels (e.g. off-target-mate
fraction). The central object of the project — it replaces the scalar read
depth used by conventional callers.
_Avoid_: fingerprint (retired term from the alignment-free design), profile

**Panel of Normals (PoN)**:
A cohort of samples believed copy-number-normal, whose fragment tensors define
the expected signal against which a test sample is compared.
_Avoid_: reference set, control cohort, baseline

**Tangent normalisation**:
Removing the shared technical structure from a sample by expressing its signal
relative to the subspace spanned by the panel of normals, leaving a residual.
_Avoid_: denoising, batch correction, PoN correction

**Residual**:
What remains of a sample's signal after tangent normalisation — the
copy-number evidence with shared technical structure removed. Segmentation and
calling operate on the residual.
_Avoid_: log-ratio, deviation, score

**Segment**:
A run of consecutive tiles sharing one copy-number state, produced by
segmenting the residual profile. A CNV call is a segment whose state is not
diploid.
_Avoid_: event, call region

**Copy-number state**:
The integer copy count assigned to a tile or segment (0, 1, 2, 3, 4+), where 2
is the diploid baseline.
_Avoid_: CN level, ploidy (ploidy is genome-wide, not per-segment)

**Run**:
The sequencing batch a sample was produced in. Recorded per sample so that
train/test splits never place same-run samples on both sides.
_Avoid_: batch, lane, flowcell
