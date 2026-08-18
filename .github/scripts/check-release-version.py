#!/usr/bin/env python3
"""Require a v-prefixed release tag to match project.version exactly."""

from __future__ import annotations

import sys
import tomllib
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: check-release-version.py v<project.version>", file=sys.stderr)
        return 2

    tag = sys.argv[1]
    if not tag.startswith("v"):
        print(f"release tag must start with 'v': {tag}", file=sys.stderr)
        return 1

    with Path("pyproject.toml").open("rb") as stream:
        version = tomllib.load(stream)["project"]["version"]

    expected = f"v{version}"
    if tag != expected:
        print(f"release tag {tag!r} does not match project version {expected!r}", file=sys.stderr)
        return 1

    print(f"release tag {tag} matches project version {version}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
