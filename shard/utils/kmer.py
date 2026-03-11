"""K-mer utilities: anchor selection, canonical forms, and external tool wrappers.

Anchor k-mer selection pipeline:

    1. Extract all k-mers from each panel target's reference sequence.
    2. Filter to genome-unique k-mers (appear exactly once in the reference).
    3. Optionally exclude k-mers overlapping known common variants.
    4. Select a well-spaced subset to ensure even coverage across each target.

Genome-wide k-mer counting is delegated to external tools (Jellyfish or KMC)
which are called via subprocess. The selection logic operates on the resulting
count dictionaries.
"""

from __future__ import annotations

import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from shard.data.panel import AnchorKmer, PanelDefinition, TargetRegion


# ---------------------------------------------------------------------------
# Canonical k-mer representation
# ---------------------------------------------------------------------------

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(_COMPLEMENT)[::-1]


def canonical_kmer(kmer: str) -> str:
    """Return the lexicographically smaller of a k-mer and its reverse complement.

    Using canonical k-mers ensures that a sequence and its reverse complement
    produce the same key, which is necessary because paired-end reads can
    come from either strand.
    """
    rc = reverse_complement(kmer)
    return kmer if kmer <= rc else rc


# ---------------------------------------------------------------------------
# K-mer extraction from sequences
# ---------------------------------------------------------------------------


def extract_kmers(sequence: str, k: int) -> list[tuple[str, int]]:
    """Extract all k-mers from a sequence with their offsets.

    Args:
        sequence: DNA sequence string (uppercase).
        k: K-mer length.

    Returns:
        List of (kmer, offset) tuples. K-mers containing 'N' are skipped.
    """
    kmers = []
    seq_upper = sequence.upper()
    for i in range(len(seq_upper) - k + 1):
        kmer = seq_upper[i : i + k]
        if "N" not in kmer:
            kmers.append((kmer, i))
    return kmers


# ---------------------------------------------------------------------------
# Genome-wide k-mer counting via external tools
# ---------------------------------------------------------------------------


def count_kmers_jellyfish(
    fasta_path: Path,
    k: int,
    *,
    jellyfish_binary: str = "jellyfish",
    threads: int = 4,
    hash_size: str = "1G",
) -> dict[str, int]:
    """Count k-mers in a FASTA file using Jellyfish.

    Args:
        fasta_path: Path to the reference genome FASTA.
        k: K-mer length.
        jellyfish_binary: Path to the jellyfish executable.
        threads: Number of threads for counting.
        hash_size: Initial hash size (e.g., "1G" for human genome).

    Returns:
        Dictionary mapping canonical k-mer → count.

    Raises:
        FileNotFoundError: If Jellyfish is not installed.
        RuntimeError: If Jellyfish exits with an error.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        jf_file = Path(tmpdir) / "counts.jf"
        dump_file = Path(tmpdir) / "counts.tsv"

        # Count
        count_cmd = [
            jellyfish_binary, "count",
            "-m", str(k),
            "-s", hash_size,
            "-t", str(threads),
            "-C",  # canonical (count both strands together)
            "-o", str(jf_file),
            str(fasta_path),
        ]
        try:
            subprocess.run(count_cmd, check=True, capture_output=True, text=True, timeout=3600)
        except FileNotFoundError:
            raise FileNotFoundError(
                f"Jellyfish not found at {jellyfish_binary!r}. "
                "Install with: conda install -c bioconda jellyfish"
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Jellyfish count failed: {e.stderr}")

        # Dump to tab-separated
        dump_cmd = [
            jellyfish_binary, "dump",
            "-c",  # column format: kmer\tcount
            "-t",  # tab separator
            str(jf_file),
        ]
        try:
            result = subprocess.run(
                dump_cmd, check=True, capture_output=True, text=True, timeout=3600
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Jellyfish dump failed: {e.stderr}")

        counts: dict[str, int] = {}
        for line in result.stdout.splitlines():
            parts = line.split("\t")
            if len(parts) == 2:
                counts[parts[0]] = int(parts[1])

        return counts


def count_kmers_kmc(
    fasta_path: Path,
    k: int,
    *,
    kmc_binary: str = "kmc",
    kmc_dump_binary: str = "kmc_dump",
    threads: int = 4,
    max_ram_gb: int = 12,
) -> dict[str, int]:
    """Count k-mers in a FASTA file using KMC.

    Args:
        fasta_path: Path to the reference genome FASTA.
        k: K-mer length.
        kmc_binary: Path to the kmc executable.
        kmc_dump_binary: Path to the kmc_dump executable.
        threads: Number of threads.
        max_ram_gb: Maximum RAM in GB.

    Returns:
        Dictionary mapping canonical k-mer → count.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        db_prefix = Path(tmpdir) / "kmc_db"
        dump_file = Path(tmpdir) / "counts.tsv"

        # Count
        count_cmd = [
            kmc_binary,
            "-k" + str(k),
            "-t" + str(threads),
            "-m" + str(max_ram_gb),
            "-fm",  # FASTA/multi-FASTA input
            "-ci1",  # min count = 1
            "-cs65535",  # max count stored
            "-b",  # canonical form
            str(fasta_path),
            str(db_prefix),
            tmpdir,
        ]
        try:
            subprocess.run(count_cmd, check=True, capture_output=True, text=True, timeout=3600)
        except FileNotFoundError:
            raise FileNotFoundError(
                f"KMC not found at {kmc_binary!r}. "
                "Install with: conda install -c bioconda kmc"
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"KMC count failed: {e.stderr}")

        # Dump
        dump_cmd = [kmc_dump_binary, str(db_prefix), str(dump_file)]
        try:
            subprocess.run(dump_cmd, check=True, capture_output=True, text=True, timeout=3600)
        except FileNotFoundError:
            raise FileNotFoundError(
                f"kmc_dump not found at {kmc_dump_binary!r}. "
                "Install with: conda install -c bioconda kmc"
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"kmc_dump failed: {e.stderr}")

        counts: dict[str, int] = {}
        for line in dump_file.read_text().splitlines():
            parts = line.split("\t")
            if len(parts) == 2:
                counts[parts[0]] = int(parts[1])

        return counts


# ---------------------------------------------------------------------------
# Variant-aware filtering
# ---------------------------------------------------------------------------


@dataclass
class VariantSite:
    """A known variant position to avoid when selecting anchors.

    Attributes:
        chrom: Chromosome.
        position: 0-based genomic position of the variant.
        allele_frequency: Population allele frequency (0-1).
    """

    chrom: str
    position: int
    allele_frequency: float = 0.0


def load_variant_sites(
    bed_path: Path,
    min_af: float = 0.01,
) -> dict[str, set[int]]:
    """Load variant sites from a BED-like file.

    Expects tab-separated: chrom, position, end, allele_frequency.
    Returns a dict mapping chrom → set of positions with AF >= min_af.

    This is intentionally simple — for production use, parse a VCF
    (e.g., gnomAD sites VCF) with cyvcf2 or pysam.
    """
    sites: dict[str, set[int]] = {}
    for line in bed_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        chrom = parts[0]
        pos = int(parts[1])
        af = float(parts[3])
        if af >= min_af:
            if chrom not in sites:
                sites[chrom] = set()
            sites[chrom].add(pos)
    return sites


def kmer_overlaps_variant(
    chrom: str,
    genomic_start: int,
    k: int,
    variant_sites: dict[str, set[int]],
) -> bool:
    """Check if a k-mer at the given genomic position overlaps any variant site."""
    chrom_sites = variant_sites.get(chrom, set())
    if not chrom_sites:
        return False
    for pos in range(genomic_start, genomic_start + k):
        if pos in chrom_sites:
            return True
    return False


# ---------------------------------------------------------------------------
# Anchor selection pipeline
# ---------------------------------------------------------------------------


def find_unique_target_kmers(
    target: TargetRegion,
    target_idx: int,
    reference_sequence: str,
    genome_counts: dict[str, int],
    k: int,
) -> list[AnchorKmer]:
    """Extract k-mers from a target that appear exactly once in the genome.

    Args:
        target: The target region.
        target_idx: Index of the target in the panel.
        reference_sequence: The target's reference DNA sequence.
        genome_counts: Genome-wide canonical k-mer counts.
        k: K-mer length.

    Returns:
        List of AnchorKmer for k-mers with genome count == 1.
    """
    candidates = []
    for kmer, offset in extract_kmers(reference_sequence, k):
        canon = canonical_kmer(kmer)
        if genome_counts.get(canon, 0) == 1:
            candidates.append(AnchorKmer(kmer=kmer, target_idx=target_idx, offset=offset))
    return candidates


def filter_variant_sites(
    candidates: list[AnchorKmer],
    target: TargetRegion,
    variant_sites: dict[str, set[int]],
    k: int,
) -> list[AnchorKmer]:
    """Remove anchor candidates that overlap known common variants.

    Args:
        candidates: Candidate anchors from find_unique_target_kmers.
        target: The target region (provides chrom and start for genomic coords).
        variant_sites: Chrom → set of variant positions.
        k: K-mer length.

    Returns:
        Filtered list with variant-overlapping k-mers removed.
    """
    filtered = []
    for anchor in candidates:
        genomic_start = target.start + anchor.offset
        if not kmer_overlaps_variant(target.chrom, genomic_start, k, variant_sites):
            filtered.append(anchor)
    return filtered


def select_spaced_anchors(
    candidates: list[AnchorKmer],
    target_length: int,
    min_spacing: int = 5,
    max_anchors: int | None = None,
) -> list[AnchorKmer]:
    """Select a well-spaced subset of anchor k-mers across a target.

    Greedily selects anchors to maximise coverage across the target,
    ensuring a minimum spacing between consecutive anchors.

    Args:
        candidates: Candidate anchors sorted by offset (output of
            find_unique_target_kmers or filter_variant_sites).
        target_length: Length of the target region in bp.
        min_spacing: Minimum distance (in bp) between selected anchors.
        max_anchors: Optional cap on anchors per target. None = no cap.

    Returns:
        Selected anchors, sorted by offset.
    """
    if not candidates:
        return []

    # Sort by offset
    sorted_cands = sorted(candidates, key=lambda a: a.offset)

    # Greedy selection with minimum spacing
    selected = [sorted_cands[0]]
    for cand in sorted_cands[1:]:
        if cand.offset - selected[-1].offset >= min_spacing:
            selected.append(cand)
            if max_anchors is not None and len(selected) >= max_anchors:
                break

    return selected


def select_anchors_for_panel(
    panel: PanelDefinition,
    target_sequences: list[str],
    genome_counts: dict[str, int],
    variant_sites: dict[str, set[int]] | None = None,
    min_spacing: int = 5,
    max_anchors_per_target: int | None = None,
    min_anchors_per_target: int = 3,
) -> tuple[list[AnchorKmer], list[int]]:
    """Run the full anchor selection pipeline for all panel targets.

    Args:
        panel: Panel definition (provides targets and k-mer size).
        target_sequences: Reference DNA sequence for each target,
            in the same order as panel.targets.
        genome_counts: Genome-wide canonical k-mer counts (from Jellyfish/KMC).
        variant_sites: Optional chrom → variant positions to avoid.
        min_spacing: Minimum bp between consecutive anchors.
        max_anchors_per_target: Optional cap per target.
        min_anchors_per_target: Warn if a target has fewer than this many anchors.

    Returns:
        Tuple of:
        - all_anchors: List of all selected AnchorKmer across all targets.
        - low_density_targets: Indices of targets with fewer than
          min_anchors_per_target anchors (may need manual review).
    """
    k = panel.kmer_size
    all_anchors: list[AnchorKmer] = []
    low_density_targets: list[int] = []

    for idx, (target, seq) in enumerate(zip(panel.targets, target_sequences)):
        # Step 1: find genome-unique k-mers
        candidates = find_unique_target_kmers(target, idx, seq, genome_counts, k)

        # Step 2: filter variant sites
        if variant_sites:
            candidates = filter_variant_sites(candidates, target, variant_sites, k)

        # Step 3: select well-spaced subset
        selected = select_spaced_anchors(
            candidates, target.length, min_spacing, max_anchors_per_target
        )

        # Track low-density targets
        if len(selected) < min_anchors_per_target:
            low_density_targets.append(idx)

        # Update target's anchor offsets
        target.anchor_offsets = [a.offset for a in selected]

        all_anchors.extend(selected)

    return all_anchors, low_density_targets


# ---------------------------------------------------------------------------
# GC content utilities
# ---------------------------------------------------------------------------


def gc_content(sequence: str) -> float:
    """Compute GC fraction of a DNA sequence."""
    seq = sequence.upper()
    gc = seq.count("G") + seq.count("C")
    total = len(seq) - seq.count("N")
    return gc / total if total > 0 else 0.0


def kmer_gc_content(kmer: str) -> float:
    """Compute GC fraction of a k-mer."""
    return gc_content(kmer)
