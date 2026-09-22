"""Download checksum-pinned CI dependencies; never build a missing package."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import tarfile
import tempfile
import urllib.request


ROOT = Path(__file__).resolve().parents[2]


def download(package, destination, repository, manifest):
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", package):
        raise ValueError("Invalid dependency package name")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", repository):
        raise ValueError("Expected a GitHub owner/repository")
    release = manifest["release"]
    if not re.fullmatch(r"ci-dependencies-[A-Za-z0-9_.-]+", release):
        raise ValueError("Invalid dependency release tag")
    checksum = manifest["sha256"].get(package)
    if not checksum or not re.fullmatch(r"[a-f0-9]{64}", checksum):
        raise ValueError(f"No pinned binary for {package}. Publish dependencies.yml and update .github/dependencies.json first.")
    if destination.exists():
        raise ValueError(f"Dependency destination must be new: {destination}")
    url = f"https://github.com/{repository}/releases/download/{release}/{package}.tar.gz"
    with tempfile.TemporaryDirectory(prefix="alps-dependencies-") as temporary:
        archive = Path(temporary) / "package.tar.gz"
        digest = hashlib.sha256()
        print(f"Downloading {url}", flush=True)
        with urllib.request.urlopen(url, timeout=120) as response, archive.open("wb") as stream:
            while chunk := response.read(1024 * 1024):
                digest.update(chunk)
                stream.write(chunk)
        if digest.hexdigest() != checksum:
            raise ValueError(f"Checksum mismatch for {package}")
        extracted = Path(temporary) / "extracted"
        with tarfile.open(archive) as package_archive:
            package_archive.extractall(extracted, filter="data")
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(extracted), destination)
    print(f"Installed {package} in {destination}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("package")
    parser.add_argument("destination", type=Path)
    parser.add_argument("--repository", default=os.environ.get("GITHUB_REPOSITORY", ""))
    args = parser.parse_args()
    manifest = json.loads((ROOT / ".github/dependencies.json").read_text())
    download(args.package, args.destination, args.repository, manifest)


if __name__ == "__main__":
    main()
