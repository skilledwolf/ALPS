"""Check release versions before building or publishing (Python >= 3.11)."""

import argparse
from email.parser import BytesParser
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import tarfile
import tomllib
import zipfile

from packaging.utils import (
    canonicalize_name,
    parse_sdist_filename,
    parse_wheel_filename,
)
from packaging.version import Version


CORE_PATTERN = r"(?:0|[1-9][0-9]*)\.(?:0|[1-9][0-9]*)\.(?:0|[1-9][0-9]*)"
TAG_PATTERN = rf"v{CORE_PATTERN}(?:-(?:alpha|beta|rc|dev)\.[0-9]+)?"


def check_version(root: Path, ref: str) -> Version:
    core = (root / "cmake/ALPS_VERSION.txt").read_text().strip()
    if not re.fullmatch(CORE_PATTERN, core):
        raise ValueError("ALPS_VERSION.txt must contain MAJOR.MINOR.PATCH")

    project_dir = root / "python/pyalps"
    with (project_dir / "pyproject.toml").open("rb") as stream:
        project = tomllib.load(stream)["project"]
    if "version" in project.get("dynamic", []):
        provider_path = Path(__file__).resolve().parents[2] / "python/pyalps/_build_support/alps_version.py"
        spec = importlib.util.spec_from_file_location("alps_version", provider_path)
        provider = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(provider)
        version = Version(provider.version(project_dir, ref=ref))
    else:
        version = Version(project["version"])
    if version.release != Version(core).release:
        raise ValueError(
            f"pyproject.toml version {version} disagrees with ALPS_VERSION.txt ({core})"
        )

    if ref.startswith("refs/tags/"):
        tag = ref.removeprefix("refs/tags/")
        if not re.fullmatch(TAG_PATTERN, tag):
            raise ValueError(
                f"Invalid release tag {tag!r}; expected vMAJOR.MINOR.PATCH"
                "[-{alpha,beta,rc,dev}.N]"
            )
        if Version(tag) != version:
            raise ValueError(
                f"Release tag {tag} disagrees with pyproject.toml ({version})"
            )

    return version


def load_matrix_policy(root: Path, profile: str) -> dict:
    policy = json.loads((root / ".github/wheel-matrix.json").read_text())
    families = policy["profiles"][profile]
    platforms = [platform for platform in policy["platforms"] if platform["family"] in families]
    if {platform["family"] for platform in platforms} != set(families):
        raise ValueError(f"Incomplete platform policy for {profile}")
    return {"platforms": platforms, "abis": policy["abis"]}


def wheel_matrix(policy: dict) -> list[dict]:
    return [
        {key: platform[key] for key in ("os", "target", "arch", "family")}
        for platform in policy["platforms"]
    ]


def check_distributions(
    directory: Path, expected: Version, policy: dict | None = None,
) -> list[Path]:
    # Keep downloads in their artifact directories until this check completes:
    # flattening first can silently overwrite a conflicting distribution.
    artifacts = sorted(path for path in directory.rglob("*") if path.is_file())
    if not artifacts:
        raise ValueError(f"No distributions found in {directory}")

    filenames = {}
    wheels = {}
    sdists = []
    expected_wheels = {
        (platform["family"], platform["arch"], interpreter, abi)
        for platform in policy["platforms"]
        for interpreter, abi in policy["abis"]
    } if policy is not None else set()

    for artifact in artifacts:
        if artifact.name in filenames:
            raise ValueError(
                f"Duplicate distribution filename {artifact.name}: "
                f"{filenames[artifact.name]} and {artifact}"
            )
        filenames[artifact.name] = artifact
        if artifact.suffix == ".whl":
            name, version, _, tags = parse_wheel_filename(artifact.name)
            if policy is not None:
                variants = set()
                for tag in tags:
                    matches = {
                        (platform["family"], platform["arch"], tag.interpreter, tag.abi)
                        for platform in policy["platforms"]
                        if re.fullmatch(platform["tag_pattern"], tag.platform)
                    }
                    if len(matches) != 1 or not matches <= expected_wheels:
                        raise ValueError(f"{artifact.name}: unexpected wheel tag {tag}")
                    variants.update(matches)
                # Multiple repair platform tags are fine when they describe the
                # same release variant; a wheel must not replace two variants.
                if len(variants) != 1:
                    raise ValueError(f"{artifact.name}: overlapping release variants")
                variant = variants.pop()
                if variant in wheels:
                    raise ValueError(
                        f"Conflicting wheels for {'/'.join(variant)}: "
                        f"{wheels[variant].name} and {artifact.name}"
                    )
                wheels[variant] = artifact
            with zipfile.ZipFile(artifact) as archive:
                members = [
                    member for member in archive.namelist()
                    if member.count("/") == 1 and member.endswith(".dist-info/METADATA")
                ]
                if len(members) != 1:
                    raise ValueError(f"{artifact.name}: expected one wheel METADATA file")
                metadata = archive.read(members[0])
        elif artifact.name.endswith(".tar.gz"):
            name, version = parse_sdist_filename(artifact.name)
            sdists.append(artifact)
            with tarfile.open(artifact) as archive:
                members = [
                    member for member in archive.getmembers()
                    if member.isfile() and member.name.count("/") == 1
                    and member.name.endswith("/PKG-INFO")
                ]
                if len(members) != 1:
                    raise ValueError(f"{artifact.name}: expected one sdist PKG-INFO file")
                with archive.extractfile(members[0]) as stream:
                    metadata = stream.read()
        else:
            raise ValueError(f"Unexpected distribution: {artifact.name}")

        if name != "pyalps" or version != expected:
            raise ValueError(f"{artifact.name}: expected a pyalps {expected} distribution")
        headers = BytesParser().parsebytes(metadata)
        if (
            canonicalize_name(headers.get("Name", "")) != "pyalps"
            or Version(headers.get("Version", "")) != expected
        ):
            raise ValueError(f"{artifact.name}: metadata does not describe pyalps {expected}")
        print(f"Verified {artifact.name}")

    if policy is not None:
        if len(sdists) != 1:
            raise ValueError(f"Expected exactly one sdist, found {len(sdists)}")
        missing = expected_wheels - wheels.keys()
        if missing:
            raise ValueError("Missing release wheels: " + ", ".join(
                "/".join(variant) for variant in sorted(missing)
            ))
    return artifacts


def stage_distributions(artifacts: list[Path], destination: Path) -> None:
    # Refuse an existing directory so stale files cannot enter publication.
    destination.mkdir(parents=True, exist_ok=False)
    for artifact in artifacts:
        shutil.copyfile(artifact, destination / artifact.name)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--ref", default=os.environ.get("GITHUB_REF", ""))
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--dist", type=Path)
    mode.add_argument("--matrix", choices=("minimal", "full"), help="Print the CI platform matrix as JSON")
    parser.add_argument("--matrix-policy", choices=("minimal", "full"), help="Require the complete distribution matrix")
    parser.add_argument("--stage-dist", type=Path, help="Copy verified distributions into a new flat directory")
    args = parser.parse_args()
    if args.matrix_policy and args.dist is None:
        parser.error("--matrix-policy requires --dist")
    if args.stage_dist and (args.dist is None or args.matrix_policy is None):
        parser.error("--stage-dist requires --dist and --matrix-policy")
    try:
        if args.matrix:
            print(json.dumps(wheel_matrix(load_matrix_policy(args.root, args.matrix)), separators=(",", ":")))
            return
        version = check_version(args.root, args.ref)
        print(f"Release version: {version}")
        if args.dist is not None:
            policy = load_matrix_policy(args.root, args.matrix_policy) if args.matrix_policy else None
            artifacts = check_distributions(args.dist, version, policy)
            if args.stage_dist:
                stage_distributions(artifacts, args.stage_dist)
    except (OSError, ValueError, RuntimeError, KeyError, tarfile.TarError, zipfile.BadZipFile) as error:
        parser.exit(1, f"Release version check failed: {error}\n")


if __name__ == "__main__":
    main()
