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
class TargetRegion:
    """A single panel target region.

    Attributes:
        name: Target name (e.g., "BRCA1_exon10").
        chrom: Chromosome (e.g., "chr17").
        start: 0-based start coordinate.
        end: 0-based exclusive end coordinate.
        gc_content: GC fraction of the target sequence.
        n_anchor_kmers: Number of unique anchor k-mers for this target.
    """

    name: str
    chrom: str
    start: int
    end: int
    gc_content: float = 0.0
    n_anchor_kmers: int = 0

    @property
    def length(self) -> int:
        return self.end - self.start

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
                    "n_anchor_kmers": t.n_anchor_kmers,
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
                    n_anchor_kmers=t.get("n_anchor_kmers", 0),
                )
                for t in data["targets"]
            ],
            seed_set=SpacedSeedSet.from_dict(data["seed_set"]),
        )
