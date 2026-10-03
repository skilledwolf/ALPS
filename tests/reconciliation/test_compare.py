"""Runner checks remain useful without a reference SDK or a native build."""

import importlib.util
import json
from pathlib import Path
import subprocess

import pytest


spec = importlib.util.spec_from_file_location("reconciliation_compare", Path(__file__).with_name("compare.py"))
compare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(compare)

COMMON = {"load": "ok", "integer": "42", "real": "1.25", "text": "example", "boolean": "1", "vector": "ok"}


@pytest.fixture
def run_report(tmp_path, monkeypatch):
    programs = {}
    for provider in ("ALPS", "ALPSCore"):
        programs[provider] = tmp_path / provider
        programs[provider].write_text("probe fixture\n")

    def run(*, core=False, corrupt=None, revisions=()):
        def probe(command, **kwargs):
            provider, action = Path(command[0]).name, command[1]
            if action == "semantics":
                values = {"native_long_bits": "64"}
            elif action == "archive-semantics":
                values = {"initial_context": "root"}
            elif action.startswith("write-"):
                values = {"save": "ok"}
            elif action == "read-archive":
                values = {"/int": "ok", "bool_type_marker": "1"}
            else:
                values = dict(COMMON)
                if action == "read-extended-params":
                    values.update(unsigned="42", float="1.25", wide="1099511627776")
            if corrupt and (provider, action) == corrupt[:2]:
                values[corrupt[2]] = corrupt[3]
            output = "".join(f"{key}\t{value}\n" for key, value in values.items())
            return subprocess.CompletedProcess(command, 0, output, "")

        monkeypatch.setattr(compare.subprocess, "run", probe)
        output = tmp_path / "report.json"
        args = ["--alps", str(programs["ALPS"]), "--output", str(output), *revisions]
        if core:
            args += ["--alpscore", str(programs["ALPSCore"])]
        compare.main(args)
        return json.loads(output.read_text())

    return run


def test_alps_only_report_does_not_claim_interchange(run_report):
    report = run_report()
    assert report["scope"] == "same-provider"
    assert report["source_revision_labels"] == {"ALPS": "unknown"}
    assert report["probe_sha256"].keys() == {"ALPS"}
    assert "ALPS->ALPS/extended-params" in report["measurements"]
    assert all("ALPSCore" not in name for name in report["measurements"])


@pytest.mark.parametrize("provider,field", [("ALPS", "wide"), ("ALPSCore", "float")])
def test_extended_self_read_rejects_silently_changed_value(run_report, provider, field):
    with pytest.raises(SystemExit, match=f"Extended params self-check failed: {provider}"):
        run_report(core=provider == "ALPSCore", corrupt=(provider, "read-extended-params", field, "0"))


def test_reference_report_keeps_explicit_revision_labels(run_report):
    report = run_report(core=True, revisions=("--alps-revision", "alps-sdk-revision", "--alpscore-revision", "core-sdk-revision"))
    assert report["scope"] == "same-provider-and-cross-provider"
    assert report["source_revision_labels"] == {"ALPS": "alps-sdk-revision", "ALPSCore": "core-sdk-revision"}
    assert "ALPS->ALPSCore/archive" in report["measurements"]
    assert "ALPSCore->ALPS/extended-params" in report["measurements"]


def test_archive_self_read_failure_is_not_labeled_interchange(run_report):
    with pytest.raises(SystemExit, match="Typed archive value check failed: ALPS->ALPS/archive"):
        run_report(corrupt=("ALPS", "read-archive", "/int", "mismatch"))
