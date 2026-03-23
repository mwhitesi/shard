"""FASTQ reading utilities for paired-end sequencing data.

Provides streaming parsers for FASTQ files (plain text or gzip-compressed)
and a paired-end reader that yields matched (R1, R2) records.

FASTQ format (4 lines per record):
    @read_name
    SEQUENCE
    +
    QUALITY
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Iterator


@dataclass(slots=True)
class FastqRecord:
    """A single FASTQ record.

    Attributes:
        name: Read name (without leading '@').
        sequence: DNA sequence string (uppercase).
        quality: Phred quality string (same length as sequence).
    """

    name: str
    sequence: str
    quality: str


def read_fastq(path: Path) -> Iterator[FastqRecord]:
    """Stream FASTQ records from a file.

    Handles both plain text (.fastq, .fq) and gzip-compressed
    (.fastq.gz, .fq.gz) files.

    Args:
        path: Path to the FASTQ file.

    Yields:
        FastqRecord for each read in the file.

    Raises:
        ValueError: If the file format is malformed.
    """
    opener = gzip.open if _is_gzipped(path) else open
    with opener(path, "rt") as fh:
        yield from _parse_fastq(fh, str(path))


def read_paired_fastq(
    r1_path: Path,
    r2_path: Path,
) -> Iterator[tuple[FastqRecord, FastqRecord]]:
    """Stream matched read pairs from paired-end FASTQ files.

    R1 and R2 files must have the same number of records in the same order.
    Read names are checked for consistency (ignoring /1, /2 or space-delimited
    suffixes).

    Args:
        r1_path: Path to the R1 (forward) FASTQ file.
        r2_path: Path to the R2 (reverse) FASTQ file.

    Yields:
        (R1 record, R2 record) tuples.

    Raises:
        ValueError: If read names don't match or files have different lengths.
    """
    r1_opener = gzip.open if _is_gzipped(r1_path) else open
    r2_opener = gzip.open if _is_gzipped(r2_path) else open

    with r1_opener(r1_path, "rt") as r1_fh, r2_opener(r2_path, "rt") as r2_fh:
        r1_iter = _parse_fastq(r1_fh, str(r1_path))
        r2_iter = _parse_fastq(r2_fh, str(r2_path))

        for r1_rec in r1_iter:
            try:
                r2_rec = next(r2_iter)
            except StopIteration:
                raise ValueError(
                    f"R1 file has more records than R2 file: "
                    f"R1={r1_path}, R2={r2_path}"
                )
            r1_base = _base_read_name(r1_rec.name)
            r2_base = _base_read_name(r2_rec.name)
            if r1_base != r2_base:
                raise ValueError(
                    f"Read name mismatch: R1={r1_rec.name!r}, R2={r2_rec.name!r}"
                )
            yield r1_rec, r2_rec

        # Check R2 doesn't have extra records
        try:
            next(r2_iter)
            raise ValueError(
                f"R2 file has more records than R1 file: "
                f"R1={r1_path}, R2={r2_path}"
            )
        except StopIteration:
            pass


def _parse_fastq(fh: IO[str], source: str) -> Iterator[FastqRecord]:
    """Parse FASTQ records from an open file handle.

    Args:
        fh: Open text-mode file handle.
        source: File path string for error messages.

    Yields:
        FastqRecord for each 4-line block.
    """
    line_num = 0
    while True:
        header = fh.readline()
        if not header:
            break  # EOF
        line_num += 1
        header = header.rstrip("\n\r")
        if not header.startswith("@"):
            raise ValueError(
                f"{source}:{line_num}: expected '@' header, got {header[:40]!r}"
            )

        sequence = fh.readline()
        if not sequence:
            raise ValueError(f"{source}:{line_num}: truncated record (missing sequence)")
        line_num += 1
        sequence = sequence.rstrip("\n\r").upper()

        plus = fh.readline()
        if not plus:
            raise ValueError(f"{source}:{line_num}: truncated record (missing '+')")
        line_num += 1

        quality = fh.readline()
        if not quality:
            raise ValueError(f"{source}:{line_num}: truncated record (missing quality)")
        line_num += 1
        quality = quality.rstrip("\n\r")

        if len(quality) != len(sequence):
            raise ValueError(
                f"{source}:{line_num}: quality length ({len(quality)}) != "
                f"sequence length ({len(sequence)})"
            )

        yield FastqRecord(
            name=header[1:],  # strip '@'
            sequence=sequence,
            quality=quality,
        )


def _base_read_name(name: str) -> str:
    """Extract base read name, stripping /1, /2 or Illumina-style suffixes.

    Illumina names: "INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y 1:N:0:INDEX"
    Old-style names: "read_name/1" or "read_name/2"
    """
    # Strip /1 or /2 suffix
    if name.endswith(("/1", "/2")):
        return name[:-2]
    # Strip space-delimited Illumina suffix (e.g., " 1:N:0:ATCACG")
    space_idx = name.find(" ")
    if space_idx != -1:
        return name[:space_idx]
    return name


def _is_gzipped(path: Path) -> bool:
    """Check if a file is gzip-compressed based on extension."""
    suffixes = path.suffixes
    return ".gz" in suffixes
