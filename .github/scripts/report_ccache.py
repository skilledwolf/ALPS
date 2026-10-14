"""Retain cumulative wheel compiler-cache snapshots without masking build failures."""

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import subprocess


ROOT = Path(__file__).resolve().parents[2]


def report(phase, cache_dir, output_dir):
    cache_dir = Path(cache_dir).resolve()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    environment = dict(os.environ, CCACHE_DIR=str(cache_dir))
    errors = []

    def capture(option):
        try:
            result = subprocess.run(
                ["ccache", option], env=environment, text=True,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30,
            )
        except (OSError, subprocess.TimeoutExpired) as error:
            errors.append(f"{option}: {error}")
            return ""
        if result.returncode:
            errors.append(f"{option}: exit {result.returncode}: {result.stderr.strip()}")
        return result.stdout

    version = capture("--version").splitlines()
    summary = capture("--show-stats")
    raw_counters = capture("--print-stats")
    # Both ccache 3.7 and 4.x provide TSV counters. Preserve their native names
    # and units rather than parsing their different human-readable summaries.
    counters = {}
    for line in raw_counters.splitlines():
        try:
            name, value = line.split()
            counters[name] = int(value)
        except ValueError:
            errors.append(f"Unrecognized statistics counter: {line!r}")
    snapshot = {
        "phase": phase,
        "recorded_at": datetime.now(timezone.utc).isoformat(),
        "cache_dir": str(cache_dir),
        "ccache_version": version[0] if version else None,
        "counters": counters,
        "errors": errors,
    }
    text = f"Compiler cache snapshot: {phase}\nCache directory: {cache_dir}\n{summary}"
    if errors:
        text += "\nStatistics reporting errors:\n" + "\n".join(errors) + "\n"
    (output_dir / f"cache-{phase}.txt").write_text(text, encoding="utf-8")
    (output_dir / f"cache-{phase}.tsv").write_text(raw_counters, encoding="utf-8")
    (output_dir / f"cache-{phase}.json").write_text(json.dumps(snapshot, indent=2) + "\n", encoding="utf-8")
    print(text, end="")
    return snapshot


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", required=True)
    parser.add_argument("--cache-dir", type=Path, default=ROOT / "_build/ccache")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "_build/wheel-reports")
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", args.phase):
        parser.error("phase must contain only letters, digits, dots, underscores, or hyphens")
    # Missing or incompatible ccache only affects observability. The build or
    # test command remains responsible for the job's success status.
    report(args.phase, args.cache_dir, args.output_dir)


if __name__ == "__main__":
    main()
