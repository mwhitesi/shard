"""Tests for FASTQ reading utilities."""

import gzip
from pathlib import Path

import pytest

from shard.utils.fastq import (
    _base_read_name,
    read_fastq,
    read_paired_fastq,
)


def _write_fastq(path: Path, records: list[tuple[str, str, str]]) -> None:
    """Write FASTQ records to a plain text file."""
    lines = []
    for name, seq, qual in records:
        lines.append(f"@{name}\n{seq}\n+\n{qual}\n")
    path.write_text("".join(lines))


def _write_fastq_gz(path: Path, records: list[tuple[str, str, str]]) -> None:
    """Write FASTQ records to a gzipped file."""
    lines = []
    for name, seq, qual in records:
        lines.append(f"@{name}\n{seq}\n+\n{qual}\n")
    with gzip.open(path, "wt") as fh:
        fh.write("".join(lines))


class TestReadFastq:
    def test_basic(self, tmp_path):
        fq = tmp_path / "reads.fq"
        _write_fastq(fq, [
            ("read1", "ACGTACGT", "IIIIIIII"),
            ("read2", "GGCCTTAA", "IIIIIIII"),
        ])
        records = list(read_fastq(fq))
        assert len(records) == 2
        assert records[0].name == "read1"
        assert records[0].sequence == "ACGTACGT"
        assert records[0].quality == "IIIIIIII"
        assert records[1].name == "read2"

    def test_gzipped(self, tmp_path):
        fq = tmp_path / "reads.fq.gz"
        _write_fastq_gz(fq, [
            ("read1", "ACGT", "IIII"),
        ])
        records = list(read_fastq(fq))
        assert len(records) == 1
        assert records[0].sequence == "ACGT"

    def test_uppercases_sequence(self, tmp_path):
        fq = tmp_path / "reads.fq"
        _write_fastq(fq, [("r1", "acgt", "IIII")])
        records = list(read_fastq(fq))
        assert records[0].sequence == "ACGT"

    def test_empty_file(self, tmp_path):
        fq = tmp_path / "empty.fq"
        fq.write_text("")
        assert list(read_fastq(fq)) == []

    def test_malformed_header(self, tmp_path):
        fq = tmp_path / "bad.fq"
        fq.write_text("BAD_HEADER\nACGT\n+\nIIII\n")
        with pytest.raises(ValueError, match="expected '@' header"):
            list(read_fastq(fq))

    def test_truncated_record(self, tmp_path):
        fq = tmp_path / "trunc.fq"
        fq.write_text("@read1\nACGT\n")
        with pytest.raises(ValueError, match="truncated"):
            list(read_fastq(fq))

    def test_quality_length_mismatch(self, tmp_path):
        fq = tmp_path / "bad_qual.fq"
        fq.write_text("@read1\nACGT\n+\nII\n")
        with pytest.raises(ValueError, match="quality length"):
            list(read_fastq(fq))


class TestReadPairedFastq:
    def test_basic_pairing(self, tmp_path):
        r1 = tmp_path / "R1.fq"
        r2 = tmp_path / "R2.fq"
        _write_fastq(r1, [
            ("read1/1", "ACGTACGT", "IIIIIIII"),
            ("read2/1", "GGCCTTAA", "IIIIIIII"),
        ])
        _write_fastq(r2, [
            ("read1/2", "TTAACCGG", "IIIIIIII"),
            ("read2/2", "AACCGGTT", "IIIIIIII"),
        ])
        pairs = list(read_paired_fastq(r1, r2))
        assert len(pairs) == 2
        assert pairs[0][0].name == "read1/1"
        assert pairs[0][1].name == "read1/2"

    def test_illumina_naming(self, tmp_path):
        r1 = tmp_path / "R1.fq"
        r2 = tmp_path / "R2.fq"
        _write_fastq(r1, [
            ("INST:1:FC:1:1:100:200 1:N:0:ATCACG", "ACGT", "IIII"),
        ])
        _write_fastq(r2, [
            ("INST:1:FC:1:1:100:200 2:N:0:ATCACG", "TGCA", "IIII"),
        ])
        pairs = list(read_paired_fastq(r1, r2))
        assert len(pairs) == 1

    def test_name_mismatch(self, tmp_path):
        r1 = tmp_path / "R1.fq"
        r2 = tmp_path / "R2.fq"
        _write_fastq(r1, [("readA/1", "ACGT", "IIII")])
        _write_fastq(r2, [("readB/2", "ACGT", "IIII")])
        with pytest.raises(ValueError, match="name mismatch"):
            list(read_paired_fastq(r1, r2))

    def test_r1_longer_than_r2(self, tmp_path):
        r1 = tmp_path / "R1.fq"
        r2 = tmp_path / "R2.fq"
        _write_fastq(r1, [
            ("r1/1", "ACGT", "IIII"),
            ("r2/1", "ACGT", "IIII"),
        ])
        _write_fastq(r2, [("r1/2", "ACGT", "IIII")])
        with pytest.raises(ValueError, match="R1 file has more"):
            list(read_paired_fastq(r1, r2))

    def test_r2_longer_than_r1(self, tmp_path):
        r1 = tmp_path / "R1.fq"
        r2 = tmp_path / "R2.fq"
        _write_fastq(r1, [("r1/1", "ACGT", "IIII")])
        _write_fastq(r2, [
            ("r1/2", "ACGT", "IIII"),
            ("r2/2", "ACGT", "IIII"),
        ])
        with pytest.raises(ValueError, match="R2 file has more"):
            list(read_paired_fastq(r1, r2))

    def test_mixed_gz_and_plain(self, tmp_path):
        r1 = tmp_path / "R1.fq"
        r2 = tmp_path / "R2.fq.gz"
        _write_fastq(r1, [("read1/1", "ACGT", "IIII")])
        _write_fastq_gz(r2, [("read1/2", "TGCA", "IIII")])
        pairs = list(read_paired_fastq(r1, r2))
        assert len(pairs) == 1


class TestBaseReadName:
    def test_strip_slash_suffix(self):
        assert _base_read_name("read1/1") == "read1"
        assert _base_read_name("read1/2") == "read1"

    def test_strip_illumina_suffix(self):
        assert _base_read_name("INST:1:FC:1:1:100:200 1:N:0:IDX") == "INST:1:FC:1:1:100:200"

    def test_no_suffix(self):
        assert _base_read_name("read1") == "read1"
