# shard

Copy-number variation (CNV) calling for targeted gene panels from aligned
short reads.

Conventional depth-based panel callers reduce each target to a single number —
its read depth. `shard` instead builds a per-locus **fragment tensor** that
keeps the *geometry* of the fragments (where reads start, how long the inserts
are, where mates land), normalises it against a **panel of normals** with
**tangent normalisation**, and calls CNVs from the residual. The aim is better
sensitivity on the cases scalar depth handles worst: partial/sub-exon events,
noisy targets, and duplications.

## Pipeline

```
FASTQ → BWA-MEM → BAM → tiling → fragment tensor → tangent normalisation
      → residual → segmentation → CNV calls
```

## Docs

- [`docs/roadmap.md`](docs/roadmap.md) — goal, pipeline, current status, and the
  testable build slices
- [`CONTEXT.md`](CONTEXT.md) — glossary (fragment tensor, tile, panel of
  normals, tangent normalisation, residual, …)
- [`docs/adr/`](docs/adr/) — architecture decision records
  ([0001](docs/adr/0001-alignment-based-fragment-tensor.md): the pivot away from
  alignment-free fingerprinting)
- [`research/`](research/) — background surveys (read-depth normalisation,
  sources of read-depth variance, the CNV-caller landscape). The alignment-free
  fingerprinting and read-pair encoder docs there are **historical** — see ADR
  0001.
- [`archive/alignment_free/`](archive/alignment_free/) — the original
  alignment-free implementation, frozen for reference

## Development

```
uv sync                              # create .venv, install deps + dev group
uv run pytest                        # tests
uv run ruff check shard/ tests/      # lint
uv run ruff format shard/ tests/     # format
```

Status: pivot in progress. See the roadmap for what's built vs. planned.
