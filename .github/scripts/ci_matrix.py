"""Select explicitly budgeted PR coverage and complete scheduled/release coverage.

Unknown paths and unavailable diffs fail open to the complete matrix. Only prose
and the independently tested optical-lattice tutorial can avoid native builds.
"""

import json
import os
from pathlib import Path, PurePosixPath
import re
import subprocess

ROOT = Path(__file__).resolve().parents[2]


def select_matrix(manifest, profile="full", *, dependencies=None):
    if profile not in {"full", "representative", "fast"}:
        raise ValueError(f"Unknown coverage profile: {profile}")
    if dependencies is None:
        dependencies = json.loads((ROOT / ".github/dependencies.json").read_text())
    defaults = {
        "standard": 17,
        "packages": "", "repository": "", "extras": False,
        "extensive": False, "python": "3.14", "mpi": "ON", "pr": False,
    }
    builds = []
    identifiers = set()
    for entry in manifest["builds"]:
        build = defaults | entry
        identifier = build["id"]
        if not re.fullmatch(r"[a-z0-9-]+", identifier) or identifier in identifiers:
            raise ValueError(f"Invalid or duplicate build id: {identifier}")
        identifiers.add(identifier)
        checksum = dependencies["boost"][build["boost"]]
        if not re.fullmatch(r"[a-f0-9]{64}", checksum):
            raise ValueError(f"Invalid Boost checksum for {identifier}")
        if profile == "full" or build["pr"]:
            builds.append(build)
    if not builds:
        raise ValueError("Empty source matrix")
    return {"include": builds}


def select_profile(event, paths):
    if event != "pull_request" or not paths:
        return "full"

    def independent(path):
        name = PurePosixPath(path)
        return (
            (name.suffix in {".md", ".rst"} and (
                path.startswith("docs/") or len(name.parts) == 1))
            or (name.suffix in {".py", ".md", ".rst"} and (
                path.startswith("tutorials/12-optical-lattice/01-bandstructure/")
                or path.startswith("tests/tutorials/")))
        )

    if all(independent(path) for path in paths):
        return "fast"
    # Build machinery, shared headers and dependency policy deserve every
    # supported configuration. Unrecognized locations do, too.
    for path in paths:
        name = PurePosixPath(path)
        if independent(path):
            continue
        if (
            path.startswith(("cmake/", ".github/", "third_party/"))
            or name.name in {"CMakeLists.txt", "pyproject.toml", "pytest.ini", "conftest.py"}
            or name.suffix in {".cmake", ".h", ".hpp", ".ipp", ".tcc", ".hxx"}
            or not path.startswith(("src/", "python/", "tests/", "tutorials/"))
        ):
            return "full"
    return "representative"


def changed_paths(environment):
    """Read GitHub's trusted event metadata; never interpolate paths in a shell."""
    if environment.get("GITHUB_EVENT_NAME") != "pull_request":
        return None
    try:
        event = json.loads(Path(environment["GITHUB_EVENT_PATH"]).read_text())
        base = event["pull_request"]["base"]["sha"]
        head = event["pull_request"]["head"]["sha"]
        if not all(re.fullmatch(r"[a-f0-9]{40}", ref) for ref in (base, head)):
            return None
        result = subprocess.run(
            ["git", "diff", "--name-only", "--no-renames", "-z", f"{base}...{head}", "--"],
            cwd=ROOT, check=True, capture_output=True, timeout=30,
        )
        return result.stdout.decode().rstrip("\0").split("\0") if result.stdout else []
    except (KeyError, ValueError, OSError, subprocess.SubprocessError):
        return None


def packaging_matrices(profile):
    versions = "11,12,13,14" if profile == "full" else "11,14"
    platforms = [
        {"os": "ubuntu-24.04", "target": "", "arch": "x86_64", "family": "manylinux"},
        {"os": "ubuntu-24.04", "target": "", "arch": "x86_64", "family": "musllinux"},
        {"os": "macos-15", "target": "15.0", "arch": "arm64", "family": "macos"},
    ]
    for platform in platforms:
        suffix = "macosx" if platform["family"] == "macos" else platform["family"]
        platform["build"] = f"cp3{{{versions}}}-{suffix}*"
    runners = [
        {"os": "ubuntu-24.04", "architecture": "x64", "artifact": "cibw-wheels-manylinux"},
        {"os": "macos-15", "architecture": "arm64", "artifact": "cibw-wheels-macos"},
    ]
    if profile == "full":
        runners.append({"os": "macos-26", "architecture": "arm64", "artifact": "cibw-wheels-macos"})
    interpreters = ["3.11", "3.12", "3.13", "3.14"] if profile == "full" else ["3.11", "3.14"]
    return {"plat": platforms}, {"plat": runners, "python": interpreters}


def main():
    manifest = json.loads((ROOT / ".github/ci-matrix.json").read_text())
    paths = changed_paths(os.environ)
    profile = select_profile(os.environ.get("GITHUB_EVENT_NAME"), paths)
    matrix = select_matrix(manifest, profile)
    wheels, smoke = packaging_matrices(profile)
    values = {
        "profile": profile,
        "native": str(profile != "fast").lower(),
        "matrix": matrix,
        "wheel_matrix": wheels,
        "smoke_matrix": smoke,
        "developer_matrix": {"os": ["ubuntu-24.04", "macos-15"] if profile == "full" else ["ubuntu-24.04"]},
    }
    output = "".join(f"{key}={json.dumps(value, separators=(',', ':')) if isinstance(value, dict) else value}\n"
                     for key, value in values.items())
    print(output, end="")
    if path := os.environ.get("GITHUB_OUTPUT"):
        with Path(path).open("a", encoding="utf-8") as stream:
            stream.write(output)
    if summary := os.environ.get("GITHUB_STEP_SUMMARY"):
        with Path(summary).open("a", encoding="utf-8") as stream:
            stream.write(f"Coverage profile: **{profile}**. " + (
                "Native and wheel builds are unnecessary for this prose/standalone-tutorial change.\n"
                if profile == "fast" else f"Source matrix: {len(matrix['include'])} builds; native, MPI, sanitizer and packaging checks required.\n"))
            if paths is None:
                stream.write("Non-PR event or unavailable diff: selected the complete matrix.\n")


if __name__ == "__main__":
    main()
