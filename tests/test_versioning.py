"""Tests for dynamic package version exposure and CI release version resolver."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import oilprops

RESOLVER_PATH = (
    Path(__file__).resolve().parents[1] / ".github" / "scripts" / "resolve_version.py"
)
_SPEC = importlib.util.spec_from_file_location("resolve_version", RESOLVER_PATH)
assert _SPEC is not None and _SPEC.loader is not None
resolve_version = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(resolve_version)


def test_package_exposes_version_string() -> None:
    assert isinstance(oilprops.__version__, str)
    assert len(oilprops.__version__.strip()) > 0


def test_initial_release_starts_at_0_1_0() -> None:
    res = resolve_version.resolve_release(
        event_name="push",
        ref_type="branch",
        ref_name="main",
        published_versions=[],
        repo_tags=[],
        head_tag=None,
    )
    assert res["version"] == "0.1.0"
    assert res["tag"] == "v0.1.0"
    assert res["needs_tag"] == "true"
    assert res["should_publish"] == "true"


def test_push_to_main_auto_bumps_patch_above_pypi_and_tags() -> None:
    res = resolve_version.resolve_release(
        event_name="push",
        ref_type="branch",
        ref_name="main",
        published_versions=[(0, 1, 0), (0, 1, 1)],
        repo_tags=[(0, 1, 0), (0, 1, 1)],
        head_tag=None,
    )
    assert res["version"] == "0.1.2"
    assert res["tag"] == "v0.1.2"
    assert res["needs_tag"] == "true"
    assert res["should_publish"] == "true"


def test_push_with_tagged_head_does_not_double_bump() -> None:
    res = resolve_version.resolve_release(
        event_name="push",
        ref_type="branch",
        ref_name="main",
        published_versions=[(0, 1, 0)],
        repo_tags=[(0, 1, 0), (0, 2, 0)],
        head_tag=(0, 2, 0),
    )
    assert res["version"] == "0.2.0"
    assert res["needs_tag"] == "false"
    assert res["should_publish"] == "true"


def test_already_published_tag_skips_cleanly() -> None:
    res = resolve_version.resolve_release(
        event_name="push",
        ref_type="tag",
        ref_name="v0.1.0",
        published_versions=[(0, 1, 0)],
        repo_tags=[(0, 1, 0)],
        head_tag=(0, 1, 0),
    )
    assert res["version"] == "0.1.0"
    assert res["should_publish"] == "false"
