"""Resolve the next semantic release version for GitHub Actions PyPI publishing."""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import urllib.error
import urllib.request
from typing import Optional, Sequence, Tuple

SEMVER_RE = re.compile(r"^v?(\d+)\.(\d+)\.(\d+)$")
PUBLISH_DIRECTIVE_RE = re.compile(r"\[publish(?::(patch|minor|major))?\]", re.IGNORECASE)
INITIAL_VERSION: Tuple[int, int, int] = (0, 1, 0)


def parse_publish_directive(commit_message: str) -> Optional[str]:
    """Return 'patch', 'minor', or 'major' if commit_message requests publishing, else None."""
    match = PUBLISH_DIRECTIVE_RE.search(commit_message or "")
    if not match:
        return None
    return (match.group(1) or "patch").lower()


def parse_semver(raw: str) -> Optional[Tuple[int, int, int]]:
    """Parse a strict X.Y.Z or vX.Y.Z version string."""
    match = SEMVER_RE.match(raw.strip())
    if not match:
        return None
    return int(match.group(1)), int(match.group(2)), int(match.group(3))


def format_semver(parts: Tuple[int, int, int]) -> str:
    """Format (major, minor, patch) tuple as X.Y.Z."""
    return f"{parts[0]}.{parts[1]}.{parts[2]}"


def bump_semver(parts: Tuple[int, int, int], bump: str) -> Tuple[int, int, int]:
    """Increment major, minor, or patch component."""
    major, minor, patch = parts
    if bump == "major":
        return (major + 1, 0, 0)
    if bump == "minor":
        return (major, minor + 1, 0)
    if bump == "patch":
        return (major, minor, patch + 1)
    raise ValueError(f"Unsupported bump level: {bump!r}")


def fetch_pypi_versions(package_name: str = "oilprops") -> list[Tuple[int, int, int]]:
    """Query PyPI JSON API for already-published X.Y.Z versions."""
    url = f"https://pypi.org/pypi/{package_name}/json"
    try:
        with urllib.request.urlopen(url, timeout=15) as resp:
            data = json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return []
        raise
    versions: list[Tuple[int, int, int]] = []
    for key in data.get("releases", {}):
        parsed = parse_semver(key)
        if parsed is not None:
            versions.append(parsed)
    return sorted(versions)


def git_tags() -> list[Tuple[int, int, int]]:
    """Return all semver git tags present in the repository."""
    out = subprocess.check_output(["git", "tag", "-l"], text=True)
    tags: list[Tuple[int, int, int]] = []
    for line in out.splitlines():
        parsed = parse_semver(line)
        if parsed is not None:
            tags.append(parsed)
    return sorted(tags)


def head_exact_tag() -> Optional[Tuple[int, int, int]]:
    """Return the semver tag pointing exactly at HEAD, if any."""
    out = subprocess.check_output(["git", "tag", "--points-at", "HEAD"], text=True)
    candidates = [parse_semver(line) for line in out.splitlines()]
    valid = [c for c in candidates if c is not None]
    return max(valid) if valid else None


def resolve_release(
    *,
    event_name: str,
    ref_type: str,
    ref_name: str,
    commit_message: str = "",
    bump_input: str = "patch",
    explicit_version: str = "",
    published_versions: Sequence[Tuple[int, int, int]],
    repo_tags: Sequence[Tuple[int, int, int]],
    head_tag: Optional[Tuple[int, int, int]],
) -> dict[str, str]:
    """
    Determine version, tag, whether a new tag must be created, and whether to publish.
    Branch pushes require '[publish]' (or '[publish:patch|minor|major]') in commit_message.
    """
    published_set = set(published_versions)

    if event_name == "workflow_dispatch" and explicit_version.strip():
        target = parse_semver(explicit_version)
        if target is None:
            raise ValueError(
                f"Invalid explicit version {explicit_version!r}; expected X.Y.Z"
            )
        needs_tag = target not in set(repo_tags)
    elif ref_type == "tag":
        target = parse_semver(ref_name)
        if target is None:
            return {
                "version": "",
                "tag": "",
                "needs_tag": "false",
                "should_publish": "false",
                "reason": f"Tag {ref_name!r} is not strict vX.Y.Z semver",
            }
        needs_tag = False
    elif event_name == "push":
        directive = parse_publish_directive(commit_message)
        if directive is None:
            return {
                "version": "",
                "tag": "",
                "needs_tag": "false",
                "should_publish": "false",
                "reason": "Commit message does not contain [publish]; skipping release",
            }
        if head_tag is not None:
            target = head_tag
            needs_tag = False
        else:
            all_known = list(published_versions) + list(repo_tags)
            target = INITIAL_VERSION if not all_known else bump_semver(max(all_known), directive)
            needs_tag = target not in set(repo_tags)
    else:
        all_known = list(published_versions) + list(repo_tags)
        if not all_known:
            target = INITIAL_VERSION
        else:
            target = bump_semver(max(all_known), bump_input or "patch")
        needs_tag = target not in set(repo_tags)

    version_str = format_semver(target)
    if target in published_set:
        return {
            "version": version_str,
            "tag": f"v{version_str}",
            "needs_tag": "false",
            "should_publish": "false",
            "reason": f"Version {version_str} is already published on PyPI",
        }

    return {
        "version": version_str,
        "tag": f"v{version_str}",
        "needs_tag": "true" if needs_tag else "false",
        "should_publish": "true",
        "reason": f"Ready to publish {version_str}",
    }


def main() -> int:
    commit_msg = os.environ.get("COMMIT_MESSAGE", "")
    if not commit_msg:
        try:
            commit_msg = subprocess.check_output(
                ["git", "log", "-1", "--pretty=%B"], text=True
            )
        except Exception:
            commit_msg = ""

    result = resolve_release(
        event_name=os.environ.get("GITHUB_EVENT_NAME", "push"),
        ref_type=os.environ.get("GITHUB_REF_TYPE", "branch"),
        ref_name=os.environ.get("GITHUB_REF_NAME", "main"),
        commit_message=commit_msg,
        bump_input=os.environ.get("INPUT_BUMP", "patch"),
        explicit_version=os.environ.get("INPUT_VERSION", ""),
        published_versions=fetch_pypi_versions("oilprops"),
        repo_tags=git_tags(),
        head_tag=head_exact_tag(),
    )
    print(json.dumps(result, indent=2))
    output_path = os.environ.get("GITHUB_OUTPUT")
    if output_path:
        with open(output_path, "a", encoding="utf-8") as fh:
            for k, v in result.items():
                fh.write(f"{k}={v}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
