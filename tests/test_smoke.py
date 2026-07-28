"""Smoke tests: the pivoted package imports cleanly.

The alignment-free implementation and its tests now live under
``archive/alignment_free/`` (see ADR 0001). Feature-level tests return as the
new slices land (tiling, fragment extraction, tensor, tangent normalisation).
"""

import importlib


def test_package_imports():
    assert importlib.import_module("shard") is not None


def test_forward_stubs_import():
    for name in [
        "shard.data.dataset",
        "shard.data.normalisation",
        "shard.utils.io",
        "shard.models.cnv_caller",
        "shard.training.trainer",
        "shard.inference.predict",
    ]:
        assert importlib.import_module(name) is not None
