"""Check compiler-cache reports across versions and failed builds."""

import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("report_ccache", ROOT / ".github/scripts/report_ccache.py")
reporter = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(reporter)


@pytest.mark.parametrize("version,summary,counters", [
    ("3.7.12", "cache hit (direct) 24\ncache miss 205\nmax cache size 256.0 MB\n",
     "cache_hit_direct\t24\ncache_miss\t205\ncache_size_kib\t230000\n"),
    ("4.13.5", "Cacheable calls: 229\n  Hits: 24\n  Misses: 205\n",
     "direct_cache_hit\t24\ncache_miss\t205\nlocal_storage_size_kibibyte\t230000\n"),
])
def test_retains_native_counters_and_uses_requested_cache(tmp_path, monkeypatch, capsys, version, summary, counters):
    cache = tmp_path / "wheel cache"
    output = tmp_path / "reports"
    monkeypatch.setenv("CCACHE_DIR", str(tmp_path / "wrong-cache"))
    responses = {"--version": f"ccache version {version}\n", "--show-stats": summary, "--print-stats": counters}

    def run(command, **kwargs):
        assert kwargs["env"]["CCACHE_DIR"] == str(cache)
        return subprocess.CompletedProcess(command, 0, responses[command[1]], "")

    monkeypatch.setattr(reporter.subprocess, "run", run)
    snapshot = reporter.report("sdk", cache, output)
    assert snapshot["counters"]["cache_miss"] == 205
    assert snapshot["ccache_version"] == f"ccache version {version}"
    assert snapshot["errors"] == []
    assert json.loads((output / "cache-sdk.json").read_text()) == snapshot
    assert (output / "cache-sdk.tsv").read_text() == counters
    assert summary in (output / "cache-sdk.txt").read_text()
    assert summary in capsys.readouterr().out


@pytest.mark.parametrize("failure", [
    FileNotFoundError("ccache is unavailable"),
    subprocess.TimeoutExpired("ccache", 30),
    subprocess.CompletedProcess(["ccache"], 1, "", "Cannot read cache"),
])
def test_reporting_failure_is_retained_without_failing_build(tmp_path, monkeypatch, failure):
    def run(*args, **kwargs):
        if isinstance(failure, Exception):
            raise failure
        return failure

    monkeypatch.setattr(reporter.subprocess, "run", run)
    snapshot = reporter.report("full", tmp_path / "cache", tmp_path / "reports")
    assert len(snapshot["errors"]) == 3
    assert snapshot["counters"] == {}
    assert "Statistics reporting errors" in (tmp_path / "reports/cache-full.txt").read_text()


def test_default_paths_are_checkout_relative_and_preserve_dotted_phase(tmp_path, monkeypatch):
    root = tmp_path / "checkout"
    monkeypatch.setattr(reporter, "ROOT", root)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("CCACHE_DIR", str(tmp_path / "unrelated-cache"))
    monkeypatch.setattr(sys, "argv", ["report_ccache.py", "--phase", "wheels-3.12.15"])

    def run(command, **kwargs):
        assert kwargs["env"]["CCACHE_DIR"] == str(root / "_build/ccache")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(reporter.subprocess, "run", run)
    reporter.main()
    assert (root / "_build/wheel-reports/cache-wheels-3.12.15.json").exists()
