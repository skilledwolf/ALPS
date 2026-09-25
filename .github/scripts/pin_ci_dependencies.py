"""Import a complete published SHA256SUMS file into the dependency manifest."""

import argparse
import json
from pathlib import Path
import re

from ci_matrix import select_matrix


ROOT = Path(__file__).resolve().parents[2]


def parse_checksums(text, required):
    checksums = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        match = re.fullmatch(r"([a-fA-F0-9]{64})\s+\*?(?:\./)?([A-Za-z0-9_.-]+)\.tar\.gz", line)
        if not match:
            raise ValueError(f"Invalid checksum line: {line!r}")
        checksum, package = match.groups()
        if package in checksums:
            raise ValueError(f"Duplicate archive: {package}")
        checksums[package] = checksum.lower()
    if set(checksums) != set(required):
        raise ValueError(f"Incomplete dependency set: missing {set(required) - checksums.keys()}, unexpected {checksums.keys() - set(required)}")
    return dict(sorted(checksums.items()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checksums", type=Path)
    parser.add_argument("--repository", required=True)
    parser.add_argument("--release", required=True)
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", args.repository):
        parser.error("Expected owner/repository")
    if not re.fullmatch(r"ci-dependencies-[A-Za-z0-9_.-]+", args.release):
        parser.error("Expected a ci-dependencies- release tag")
    source = json.loads((ROOT / ".github/ci-matrix.json").read_text())
    required = {build["dependency"] for build in select_matrix(source, "full")["include"]}
    required.update({"windows-x64", "windows-arm64", f"boost-{source['boost_default']}-manylinux-x64"})
    checksums = parse_checksums(args.checksums.read_text(), required)
    path = ROOT / ".github/dependencies.json"
    manifest = json.loads(path.read_text())
    manifest.update(repository=args.repository, release=args.release, sha256=checksums)
    path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Pinned {len(checksums)} archives from {args.repository}/{args.release}. Run full CI before merging.")


if __name__ == "__main__":
    main()
