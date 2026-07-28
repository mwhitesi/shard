# Slice 1 — Tiling

Expanded plan for the first slice of the alignment-based design (see
[roadmap](roadmap.md), [ADR 0001](adr/0001-alignment-based-fragment-tensor.md)).

**Scope:** turn a targets BED into the ordered set of tiles that every
downstream slice is written against, and persist it as a panel artifact.

**In scope:** BED reading, tile geometry, tile↔target annotation, chromosome
clipping, tile ordering, panel save/load with integrity checking.
**Out of scope:** anything needing a BAM (Slice 2), a reference genome, or a
mappability/blacklist track (a later panel-annotation step).

## Decisions

| # | Decision | Choice | Why |
|---|---|---|---|
| 1 | Tile anchoring | **Absolute genomic grid** — cells at fixed multiples of `tile_width` | Uniform width makes tiles comparable; adjacent/overlapping targets dedupe to a set of cells for free; stable under BED edits |
| 2 | Flank unit | **Base pairs in config, rounded up to whole tiles** at build time | Knob stays in biological units; every target gets the same flank tile count regardless of grid phase |
| 3 | Flank default | **500bp** (5 tiles at `tile_width=100`) | Matches `insert_size_range` max — flank covers one full fragment length, so the two config knobs agree |
| 4 | Tile annotation | **Names of all overlapping targets + in-target bp** | A call must be reportable as "BRCA1 exon 10"; handles two exons in one cell honestly; overlap bp is free QC |
| 5 | Chromosome lengths | **Optional, caller-supplied** mapping | Keeps tiling pure and synthetically testable; no mandatory reference file in Slice 1 |
| 6 | Tile off chromosome end | **Truncate and keep short**; tiles wholly past the end are dropped | Never loses coverage, including on short contigs like chrM |
| 7 | Tile order | **Reference order** (from `chrom_lengths` if given, else BED first-appearance), **frozen into the panel file** | Tile order is the row order of every tensor; freezing removes the risk of an optional input silently changing it |
| 8 | Panel file contents | **Inputs + tile checksum**; tiles rebuilt on load and verified | Small at any panel scale; a tiling change fails loudly instead of misaligning a cohort |
| 9 | Overlaps and names in BED | **Accept overlaps**; same name at different coordinates is an **error** | Overlapping transcripts are normal input; ambiguous names make calls unreportable |
| 10 | In-memory layout | **Columnar numpy arrays** + on-demand `Tile` view | Scales to exome (~2.4M tiles); keeps per-fragment lookup allocation-free |
| 11 | Module split | **`tiling.py`** (pure geometry) + **`panel.py`** (BED reader, artifact) | Separates pure logic from I/O; matches the roadmap's status table |
| 12 | Tile annotations | **Geometry only** — no GC, mappability, or blacklist | No new dependencies; gives Slice 5 a clean test of whether tangent normalisation alone suffices |
| 13 | Test input | **Inline for geometry**, `tmp_path` for loader, **one committed messy BED** | Input sits next to expectation; the fixture doubles as a worked example of real input |
| 14 | `tile_width` default | **100bp** | Most exons become two tiles; ~333 fragments/tile at 1000x. Sub-tile resolution comes from the tensor's position axis |

### Smaller decisions recorded without debate

- 3-column BED gets auto-generated names (`chr1:1030-1250`); name comes from
  column 4 when present; further columns ignored.
- `#`, `track`, and `browser` lines are skipped.
- Zero-length or reversed intervals are an error.
- Exactly duplicated entries (same name, same coordinates) collapse silently.
- A target on a chromosome absent from a supplied `chrom_lengths` is an error,
  not a skip.
- Library code takes explicit arguments and never reads `configs/default.yaml`;
  config wiring belongs to the Slice 7 scripts.
- `TargetRegion.gc_content` from the archived code is **dropped** (decision 12).
- "Checksum", never "fingerprint" — that term is retired in
  [CONTEXT.md](../CONTEXT.md).

## How tiling works

```
tile_width = 100, flank_bp = 500  ->  n_flank_tiles = ceil(500/100) = 5

ruler   | 1000 | 1100 | 1200 | 1300 |     fixed cells, never move
exon              [== 1030-1250 ==]
                  ^               ^       lands mid-cell at both ends

target cells : 1000-1100, 1100-1200, 1200-1300
flank cells  : 5 more each side (500-1000 and 1300-1800)
total        : 13 tiles for this target
```

Algorithm:

1. Validate `tile_width > 0`, `flank_bp >= 0`.
2. `n_flank_tiles = ceil(flank_bp / tile_width)`.
3. Resolve chromosome order: keys of `chrom_lengths` if supplied, else order of
   first appearance in the targets.
4. Per target, compute the cells it touches:
   `first = start // tile_width`, `last = (end - 1) // tile_width`.
   Extend by `n_flank_tiles` each side; clamp at cell 0.
5. Accumulate cells into a map `(chrom, cell) → {target names, in-target bp}`.
   Cells reached only as flank contribute no name and no bp — so a shared cell
   between two nearby targets appears **once**, carrying both names.
6. In-target bp is the **union** of target coverage within the cell, not the
   sum, so overlapping transcripts do not double-count.
7. If `chrom_lengths` is supplied: drop cells starting at or past the
   chromosome length; truncate a straddling cell's `end` to the length.
8. Sort by `(chrom rank, cell index)` and assign tile indices.
9. Compute the tile checksum over the ordered `(chrom, start, end)` tuples plus
   `tile_width` and `flank_bp`.

Two consequences worth stating outright:

- **Fragment→tile assignment is arithmetic, not a search.** `pos // tile_width`
  gives the cell, then one dict lookup. That matters in Slice 2, which does this
  hundreds of millions of times.
- **Flank-only tiles carry no target name.** Reporting a call therefore derives
  its name from the target-overlapping tiles inside the segment, not from every
  tile in it.

## Data model

```python
# shard/data/panel.py
@dataclass(frozen=True)
class TargetRegion:
    name: str
    chrom: str
    start: int          # 0-based
    end: int            # exclusive

def load_targets_bed(path: Path) -> list[TargetRegion]: ...

@dataclass
class PanelDefinition:
    panel_id: str
    tile_width: int
    flank_bp: int
    chroms: tuple[str, ...]                  # frozen order (decision 7)
    chrom_lengths: dict[str, int] | None
    targets: list[TargetRegion]
    source_bed: str | None                   # provenance
    source_bed_sha256: str | None
    tile_checksum: str | None

    @property
    def tiles(self) -> TileSet: ...          # built on first access, cached
    def save(self, path: Path) -> None: ...
    @classmethod
    def load(cls, path: Path) -> PanelDefinition: ...   # verifies checksum
```

```python
# shard/data/tiling.py
def tile_targets(
    targets: Sequence[TargetRegion],
    tile_width: int,
    flank_bp: int,
    chrom_order: Sequence[str] | None = None,
    chrom_lengths: Mapping[str, int] | None = None,
) -> TileSet: ...

class TileSet:
    chroms: tuple[str, ...]
    chrom_index: NDArray[np.int32]           # parallel arrays, one per tile
    start: NDArray[np.int64]
    end: NDArray[np.int64]
    target_bp: NDArray[np.int32]
    target_names: list[tuple[str, ...]]
    tile_width: int

    def __len__(self) -> int: ...
    def __getitem__(self, i: int) -> Tile: ...           # view, built on demand
    def lookup(self, chrom: str, pos: int) -> int | None: ...
    def checksum(self) -> str: ...

@dataclass(frozen=True)
class Tile:
    index: int
    chrom: str
    start: int
    end: int
    target_names: tuple[str, ...]
    target_bp: int
    # width = end - start; NEVER assume it equals tile_width (decision 6)
    # is_flank = target_bp == 0
```

## Activities

Each is independently testable and lands as its own commit.

1. **`TargetRegion` + BED reader** — parsing, auto-naming, comment/track
   skipping, duplicate collapsing, ambiguous-name and malformed-interval errors.
2. **Grid geometry + `TileSet`** — `tile_targets` for the no-clipping case,
   columnar arrays, `Tile` view, union in-target bp, target-name accumulation.
3. **Chromosome clipping** — clamp at 0, drop past the end, truncate the
   straddling tile.
4. **Ordering, `lookup`, checksum** — reference order resolution, sort, O(1)
   position lookup, stable checksum.
5. **`PanelDefinition` save/load** — JSON artifact, provenance fields, checksum
   verification on load.
6. **Config + fixture + docs** — `flank: 200 → 500` in `configs/default.yaml`,
   commit the messy BED fixture, flip the roadmap status row to built.

## Test plan

Roadmap-mandated cases are marked ★.

**Geometry**
- Grid-aligned target → exact expected tile list
- ★ Off-grid target (1030–1250) → 3 target tiles + 5 flank each side
- ★ Target shorter than a tile, wholly inside one cell → 1 target tile
- Target shorter than a tile but straddling a cell boundary → 2 target tiles
- ★ Adjacent targets sharing cells → no duplicate tiles; shared cell lists both
- Overlapping transcripts → both names present, in-target bp not double-counted
- Determinism: shuffled target input → byte-identical tiles and checksum

**Chromosome bounds**
- Flank running below coordinate 0 → clamped, no negative starts
- ★ Flank past the chromosome end with lengths given → tiles beyond dropped,
  straddling tile truncated to a short width
- Same panel without `chrom_lengths` → no end clipping, still clamped at 0
- Target on a chromosome missing from `chrom_lengths` → error

**Ordering and lookup**
- Multi-chromosome panel follows reference order; within a chromosome, by start
- `lookup` returns the right index inside a tile and `None` in an untiled gap
- `lookup` still resolves correctly for a truncated end-of-chromosome tile

**BED reader**
- 3-column BED → auto-generated names
- 6-column BED → name from column 4, extras ignored
- Comment, `track`, `browser` lines skipped
- Exact duplicate entries collapse
- Same name, different coordinates → error
- Reversed / zero-length interval → error
- Committed messy fixture loads and tiles deterministically

**Panel artifact**
- Save/load roundtrip equality
- Checksum mismatch after tampering with saved parameters → loud error
- Frozen chromosome order honoured on load with no `chrom_lengths` supplied

**Validation**
- `tile_width <= 0`, `flank_bp < 0` → error

## Definition of done

- All tests above pass; `uv run ruff check shard/ tests/` clean.
- A panel builds from the fixture BED, saves, reloads, and verifies.
- `configs/default.yaml` reflects `flank: 500`.
- Roadmap status row for Panel/tiling flipped from "to build".

## Deferred, with the slice that owns it

- **Slice 3** — what a position bin means inside a truncated tile.
- **Slice 3** — histogram sparsity: ~333 fragments per tile at 1000x spread over
  10 × 45 = 450 tensor cells is under one count per cell. Fix by coarsening bin
  counts, which does not require rebuilding the panel.
- **Slice 3/5** — whether in-target bp becomes a tensor channel or a covariate.
- **Slice 7** — sweep `tile_width` and `flank_bp` against a truth set; resolve
  the panel-scale question (tens of genes vs full exome).
- **Later** — a panel-annotation step adding GC, mappability, and blacklist per
  tile (research notes' universal steps 1–3).
