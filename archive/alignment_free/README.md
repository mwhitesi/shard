# archive/alignment_free

A frozen snapshot of `shard`'s original **alignment-free** approach: characterise
a panel sample's whole read population as a fixed-dimension k-mer *fingerprint*
and detect CNVs as perturbations of it.

Retired in favour of the alignment-based fragment-tensor design — see
[ADR 0001](../../docs/adr/0001-alignment-based-fragment-tensor.md) for why.

## Status

**Not runnable in place and not part of the test suite or package.** The tests
here still `import shard.…`, which now resolves to the active (pivoted) package,
not this snapshot. To resurrect, move the `shard/` subtree back to the repo root
(and restore the archived dependencies).

## Contents

```
shard/
  extraction.py          FASTQ → raw fingerprint accumulation
  utils/fastq.py         paired-end FASTQ streaming
  utils/kmer.py          canonical k-mers, anchor selection, Jellyfish/KMC wrappers
  data/panel.py          spaced seeds (ALeS/Iedera), AnchorIndex, fragment-size estimation
  data/fingerprint.py    the 78-feature per-target matrix, RunContext, cohort stacking
tests/                   the original 95 tests for the above
```

## What carried forward

The `TargetRegion` / `PanelDefinition` types and the
`RunContext` → panel-of-normals → cohort-stacking *pattern* are re-introduced
(trimmed) in the new design. The idea of *relaxed / alignment-free placement to
sidestep alignment bias* remains a deliberate future tangent, not a dead end.

Background: `research/alignment_free_fingerprinting.md` and
`research/read_pair_encoder.md` describe this approach in full.
