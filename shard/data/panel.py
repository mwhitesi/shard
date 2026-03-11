"""Panel definition: target regions, anchor k-mer sets, and spaced seed design.

A panel defines the genomic targets for CNV calling. Each target has:
- Genomic coordinates (chrom, start, end)
- A set of anchor k-mers used for alignment-free depth estimation
- Spaced seed patterns that allow anchor matching to tolerate SNPs and
  sequencing errors

Spaced seeds are binary masks applied to k-mers so that only "care" positions
must match exactly. A weight-28, span-31 seed has 3 don't-care positions —
a SNP at any of those positions won't break the match. Multiple complementary
seeds cover different SNP positions for broad robustness.

Seed patterns are designed at panel build time using external optimisers
(ALeS or Iedera) and stored as part of the panel definition.
"""

from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from numpy.typing import NDArray


@dataclass
class SpacedSeed:
    """A single spaced seed pattern.

    Attributes:
        pattern: Binary string, e.g. "1111011110111101111011110111011".
            '1' = care position (must match), '0' = don't-care (wildcard).
        weight: Number of care positions (number of '1's).
        span: Total length of the pattern.
    """

    pattern: str

    def __post_init__(self) -> None:
        if not all(c in "01" for c in self.pattern):
            raise ValueError(
                f"Seed pattern must contain only '0' and '1', got {self.pattern!r}"
            )

    @property
    def weight(self) -> int:
        return self.pattern.count("1")

    @property
    def span(self) -> int:
        return len(self.pattern)

    @property
    def care_positions(self) -> list[int]:
        """Indices of care ('1') positions."""
        return [i for i, c in enumerate(self.pattern) if c == "1"]

    @property
    def wildcard_positions(self) -> list[int]:
        """Indices of don't-care ('0') positions."""
        return [i for i, c in enumerate(self.pattern) if c == "0"]

    def extract_key(self, kmer: str) -> str:
        """Extract the spaced key from a k-mer (bases at care positions only)."""
        if len(kmer) != self.span:
            raise ValueError(
                f"K-mer length {len(kmer)} != seed span {self.span}"
            )
        return "".join(kmer[i] for i in self.care_positions)


@dataclass
class SpacedSeedSet:
    """A set of complementary spaced seeds for SNP-tolerant k-mer matching.

    A k-mer matches an anchor if it matches via *any* seed in the set.
    Multiple seeds with wildcards at different positions provide broad
    coverage against SNPs at any location.

    Attributes:
        seeds: List of SpacedSeed patterns.
        kmer_size: The k-mer span these seeds are designed for.
    """

    seeds: list[SpacedSeed]
    kmer_size: int

    def __post_init__(self) -> None:
        for seed in self.seeds:
            if seed.span != self.kmer_size:
                raise ValueError(
                    f"Seed span {seed.span} != kmer_size {self.kmer_size}"
                )

    @property
    def n_seeds(self) -> int:
        return len(self.seeds)

    @property
    def weight(self) -> int:
        """Weight of the seeds (assumed uniform across the set)."""
        return self.seeds[0].weight

    def wildcard_coverage(self) -> NDArray[np.bool_]:
        """Boolean mask of positions covered by at least one wildcard.

        Shape (kmer_size,). True at positions where at least one seed has
        a don't-care, meaning a SNP at that position is tolerated by at
        least one seed.
        """
        covered = np.zeros(self.kmer_size, dtype=np.bool_)
        for seed in self.seeds:
            for pos in seed.wildcard_positions:
                covered[pos] = True
        return covered

    @property
    def fraction_positions_covered(self) -> float:
        """Fraction of positions where a SNP would be tolerated by at least one seed."""
        return float(self.wildcard_coverage().mean())

    @classmethod
    def exact_match(cls, kmer_size: int) -> SpacedSeedSet:
        """Create a seed set equivalent to exact k-mer matching (all care positions)."""
        return cls(
            seeds=[SpacedSeed(pattern="1" * kmer_size)],
            kmer_size=kmer_size,
        )

    def to_dict(self) -> dict:
        return {
            "kmer_size": self.kmer_size,
            "seeds": [s.pattern for s in self.seeds],
        }

    @classmethod
    def from_dict(cls, d: dict) -> SpacedSeedSet:
        kmer_size = d["kmer_size"]
        return cls(
            seeds=[SpacedSeed(pattern=p) for p in d["seeds"]],
            kmer_size=kmer_size,
        )


# ---------------------------------------------------------------------------
# Seed generation via external tools
# ---------------------------------------------------------------------------


def generate_seeds_ales(
    weight: int,
    n_seeds: int,
    similarity: float,
    alignment_length: int,
    span: int,
    *,
    ales_binary: str = "ALeS",
) -> SpacedSeedSet:
    """Generate optimised spaced seeds using ALeS.

    ALeS (Adaptive Length Seeds) supports high-weight seeds (weight > 16)
    that Iedera cannot handle, making it suitable for k=31 anchors.

    Args:
        weight: Number of care positions per seed (e.g., 28).
        n_seeds: Number of complementary seeds to generate.
        similarity: Expected sequence similarity (e.g., 0.85 for 15% mismatch).
        alignment_length: Alignment length for sensitivity evaluation.
        span: Maximum seed span (typically = kmer_size).
        ales_binary: Path to the ALeS executable.

    Returns:
        SpacedSeedSet with the optimised seeds.

    Raises:
        FileNotFoundError: If ALeS binary is not found.
        RuntimeError: If ALeS returns no seeds or exits with an error.
    """
    cmd = [ales_binary, str(weight), str(n_seeds), str(similarity),
           str(alignment_length), str(span)]
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True, timeout=300,
        )
    except FileNotFoundError:
        raise FileNotFoundError(
            f"ALeS binary not found at {ales_binary!r}. "
            "Install from https://github.com/lucian-ilie/ALeS"
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"ALeS failed (exit {e.returncode}): {e.stderr}")

    seeds = _parse_ales_output(result.stdout, span)
    if not seeds:
        raise RuntimeError(f"ALeS produced no seeds. stdout:\n{result.stdout}")
    return SpacedSeedSet(seeds=seeds, kmer_size=span)


def generate_seeds_iedera(
    weight: int,
    n_seeds: int,
    span_min: int,
    span_max: int,
    *,
    match_prob: float = 0.85,
    alignment_length: int = 64,
    n_random: int = 100_000,
    iedera_binary: str = "iedera",
) -> SpacedSeedSet:
    """Generate optimised spaced seeds using Iedera.

    Note: Iedera's documented weight limit is 16. For weight > 16,
    use generate_seeds_ales() instead.

    Args:
        weight: Number of care positions per seed.
        n_seeds: Number of complementary seeds.
        span_min: Minimum seed span.
        span_max: Maximum seed span.
        match_prob: Per-base match probability in foreground model.
        alignment_length: Alignment length for sensitivity evaluation.
        n_random: Number of random seeds to sample before hill-climbing.
        iedera_binary: Path to the Iedera executable.

    Returns:
        SpacedSeedSet with the optimised seeds.
    """
    mismatch_prob = 1.0 - match_prob
    cmd = [
        iedera_binary, "-spaced",
        "-w", f"{weight},{weight}",
        "-s", f"{span_min},{span_max}",
        "-n", str(n_seeds),
        "-l", str(alignment_length),
        "-r", str(n_random),
        "-k",
        "-f", f"{mismatch_prob:.4f},{match_prob:.4f}",
    ]
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True, timeout=600,
        )
    except FileNotFoundError:
        raise FileNotFoundError(
            f"Iedera binary not found at {iedera_binary!r}. "
            "Install from https://github.com/laurentnoe/iedera"
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Iedera failed (exit {e.returncode}): {e.stderr}")

    seeds = _parse_iedera_output(result.stdout)
    if not seeds:
        raise RuntimeError(f"Iedera produced no seeds. stdout:\n{result.stdout}")
    span = seeds[0].span
    return SpacedSeedSet(seeds=seeds, kmer_size=span)


def _parse_ales_output(stdout: str, expected_span: int) -> list[SpacedSeed]:
    """Parse ALeS stdout into SpacedSeed objects.

    ALeS outputs one binary seed per line (e.g., "1110111011101110111011101110111").
    Lines not matching the expected pattern are skipped.
    """
    seeds = []
    for line in stdout.strip().splitlines():
        line = line.strip()
        if line and all(c in "01" for c in line) and len(line) == expected_span:
            seeds.append(SpacedSeed(pattern=line))
    return seeds


def _parse_iedera_output(stdout: str) -> list[SpacedSeed]:
    """Parse Iedera stdout into SpacedSeed objects.

    Iedera outputs tab-separated lines. The first column is the seed pattern
    using '#' (care) and '-' (don't-care).
    """
    seeds = []
    for line in stdout.strip().splitlines():
        parts = line.split()
        if not parts:
            continue
        pattern_str = parts[0]
        if all(c in "#-" for c in pattern_str):
            binary = pattern_str.replace("#", "1").replace("-", "0")
            seeds.append(SpacedSeed(pattern=binary))
    return seeds


# ---------------------------------------------------------------------------
# Panel target and panel definition
# ---------------------------------------------------------------------------


@dataclass
class AnchorKmer:
    """A single anchor k-mer with its known genomic position.

    Attributes:
        kmer: The k-mer sequence string.
        target_idx: Index of the target this anchor belongs to.
        offset: 0-based offset of the k-mer's start within the target region.
            The genomic position is target.start + offset.
    """

    kmer: str
    target_idx: int
    offset: int


@dataclass
class TargetRegion:
    """A single panel target region.

    Attributes:
        name: Target name (e.g., "BRCA1_exon10").
        chrom: Chromosome (e.g., "chr17").
        start: 0-based start coordinate.
        end: 0-based exclusive end coordinate.
        gc_content: GC fraction of the target sequence.
        anchor_offsets: Sorted offsets (within target) of each anchor k-mer.
            Length equals the number of unique anchors for this target.
    """

    name: str
    chrom: str
    start: int
    end: int
    gc_content: float = 0.0
    anchor_offsets: list[int] = field(default_factory=list)

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def n_anchor_kmers(self) -> int:
        return len(self.anchor_offsets)

    @property
    def anchor_density(self) -> float:
        """Unique anchor k-mers per base pair."""
        return self.n_anchor_kmers / self.length if self.length > 0 else 0.0


@dataclass
class PanelDefinition:
    """Complete panel definition with targets and spaced seed configuration.

    Attributes:
        panel_id: Unique panel identifier (e.g., "panel_v3").
        kmer_size: K-mer length used for anchor design.
        targets: Ordered list of target regions.
        seed_set: Spaced seed patterns for SNP-tolerant anchor matching.
    """

    panel_id: str
    kmer_size: int
    targets: list[TargetRegion]
    seed_set: SpacedSeedSet

    @property
    def n_targets(self) -> int:
        return len(self.targets)

    @property
    def target_names(self) -> list[str]:
        return [t.name for t in self.targets]

    def save(self, path: Path) -> None:
        """Serialise panel definition to JSON."""
        data = {
            "panel_id": self.panel_id,
            "kmer_size": self.kmer_size,
            "seed_set": self.seed_set.to_dict(),
            "targets": [
                {
                    "name": t.name,
                    "chrom": t.chrom,
                    "start": t.start,
                    "end": t.end,
                    "gc_content": t.gc_content,
                    "anchor_offsets": t.anchor_offsets,
                }
                for t in self.targets
            ],
        }
        path.write_text(json.dumps(data, indent=2))

    @classmethod
    def load(cls, path: Path) -> PanelDefinition:
        """Load panel definition from JSON."""
        data = json.loads(path.read_text())
        return cls(
            panel_id=data["panel_id"],
            kmer_size=data["kmer_size"],
            targets=[
                TargetRegion(
                    name=t["name"],
                    chrom=t["chrom"],
                    start=t["start"],
                    end=t["end"],
                    gc_content=t.get("gc_content", 0.0),
                    anchor_offsets=t.get("anchor_offsets", []),
                )
                for t in data["targets"]
            ],
            seed_set=SpacedSeedSet.from_dict(data["seed_set"]),
        )


# ---------------------------------------------------------------------------
# Anchor index and fragment size estimation
# ---------------------------------------------------------------------------


@dataclass
class AnchorHit:
    """A single anchor k-mer match from a read.

    Attributes:
        target_idx: Panel target this anchor belongs to.
        offset: Anchor position within the target (0-based).
        read_position: Position of the k-mer within the read (0-based from 5').
    """

    target_idx: int
    offset: int
    read_position: int


class AnchorIndex:
    """Lookup structure mapping k-mer keys to anchor positions.

    Supports both exact and spaced-seed matching. For each anchor k-mer
    in the panel, stores the (target_idx, offset) so that a k-mer match
    during FASTQ scanning can be resolved to a known genomic position.

    The index maps spaced keys (care-position bases) rather than full
    k-mers, so a single anchor produces one entry per seed pattern.
    """

    def __init__(self, panel: PanelDefinition, anchors: list[AnchorKmer]) -> None:
        """Build the anchor index.

        Args:
            panel: Panel definition (provides seed set and targets).
            anchors: All anchor k-mers with their target assignments and offsets.
        """
        self.panel = panel
        self.seed_set = panel.seed_set
        # Map: spaced_key → list of (target_idx, offset)
        # Multiple seeds means each anchor is indexed multiple times under
        # different keys.
        self._index: dict[str, list[tuple[int, int]]] = {}
        for anchor in anchors:
            for seed in self.seed_set.seeds:
                key = seed.extract_key(anchor.kmer)
                if key not in self._index:
                    self._index[key] = []
                self._index[key].append((anchor.target_idx, anchor.offset))

    def lookup(self, kmer: str) -> list[tuple[int, int]]:
        """Look up a k-mer against all seed patterns.

        Returns a list of (target_idx, offset) for each matching anchor.
        Tries each seed pattern and returns the union of matches.
        De-duplicates by (target_idx, offset).
        """
        seen: set[tuple[int, int]] = set()
        results: list[tuple[int, int]] = []
        for seed in self.seed_set.seeds:
            key = seed.extract_key(kmer)
            for hit in self._index.get(key, []):
                if hit not in seen:
                    seen.add(hit)
                    results.append(hit)
        return results

    @property
    def n_entries(self) -> int:
        """Total number of (key → anchor) entries in the index."""
        return sum(len(v) for v in self._index.values())


def estimate_fragment_sizes(
    r1_hits: list[AnchorHit],
    r2_hits: list[AnchorHit],
    panel: PanelDefinition,
    read_length: int = 150,
) -> list[tuple[int, int]]:
    """Estimate fragment sizes from anchor hits in a read pair.

    For each pair of anchor hits (one from R1, one from R2) that map to the
    same target, compute the implied fragment size from their known genomic
    distance.

    Fragment size is estimated as:
        |offset_r1 - offset_r2| + read_length

    This works because:
    - offset_r1 is where R1's anchor sits in the target (forward strand)
    - offset_r2 is where R2's anchor sits (reverse complement, so it
      represents the other end of the fragment)
    - The distance between them plus one read length approximates the
      full fragment span

    When multiple anchor pairs are available (common — each read may hit
    several anchors), we get multiple fragment size estimates. The median
    of consistent estimates is the best single estimate.

    Args:
        r1_hits: Anchor hits from R1.
        r2_hits: Anchor hits from R2.
        panel: Panel definition (for target coordinates).
        read_length: Sequencing read length in bp.

    Returns:
        List of (target_idx, estimated_fragment_size) tuples, one per
        concordant anchor pair. Empty if no same-target pairs found.
    """
    # Group hits by target for efficient pairing
    r1_by_target: dict[int, list[AnchorHit]] = {}
    for hit in r1_hits:
        r1_by_target.setdefault(hit.target_idx, []).append(hit)

    r2_by_target: dict[int, list[AnchorHit]] = {}
    for hit in r2_hits:
        r2_by_target.setdefault(hit.target_idx, []).append(hit)

    estimates: list[tuple[int, int]] = []
    for target_idx in r1_by_target:
        if target_idx not in r2_by_target:
            continue
        for h1 in r1_by_target[target_idx]:
            for h2 in r2_by_target[target_idx]:
                genomic_dist = abs(h1.offset - h2.offset)
                frag_size = genomic_dist + read_length
                estimates.append((target_idx, frag_size))

    return estimates


def build_fragment_distribution(
    all_estimates: list[tuple[int, int]],
    n_targets: int,
    frag_range: tuple[int, int] = (50, 500),
) -> tuple[NDArray[np.float32], NDArray[np.float32]]:
    """Build per-target and global fragment size distributions from estimates.

    Args:
        all_estimates: List of (target_idx, fragment_size) from all read pairs
            in a sample.
        n_targets: Number of targets in the panel.
        frag_range: (min_size, max_size) for binning.

    Returns:
        Tuple of:
        - per_target_hist: shape (n_targets, n_bins), raw counts per target.
        - global_hist: shape (n_bins,), raw counts across all targets.
    """
    frag_min, frag_max = frag_range
    n_bins = frag_max - frag_min

    per_target = np.zeros((n_targets, n_bins), dtype=np.float32)
    global_hist = np.zeros(n_bins, dtype=np.float32)

    for target_idx, frag_size in all_estimates:
        bin_idx = frag_size - frag_min
        if 0 <= bin_idx < n_bins:
            per_target[target_idx, bin_idx] += 1
            global_hist[bin_idx] += 1

    return per_target, global_hist
