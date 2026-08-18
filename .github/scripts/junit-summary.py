#!/usr/bin/env python3
"""Render a ctest JUnit XML report as GitHub step-summary Markdown.

Never exits non-zero: reporting must not fail the job.
"""
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


def main():
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("build/junit.xml")
    if not path.is_file():
        print(f"_No JUnit report found at `{path}`._")
        return
    root = ET.parse(path).getroot()
    suites = [root] if root.tag == "testsuite" else list(root)
    total = failures = skipped = 0
    elapsed = 0.0
    failed = []
    for s in suites:
        total += int(s.get("tests", 0) or 0)
        failures += int(s.get("failures", 0) or 0) + int(s.get("errors", 0) or 0)
        skipped += int(s.get("skipped", 0) or 0) + int(s.get("disabled", 0) or 0)
        elapsed += float(s.get("time", 0) or 0)
        for case in s.iter("testcase"):
            if case.find("failure") is not None or case.find("error") is not None:
                failed.append(case.get("name", "?"))
    icon = "✅" if failures == 0 else "❌"
    passed = total - failures - skipped
    print(f"### {icon} Tests: {passed} passed, {failures} failed, "
          f"{skipped} skipped ({elapsed:.0f}s)")
    if failed:
        print()
        print("| Failed test |")
        print("|---|")
        for name in failed[:50]:
            print(f"| `{name}` |")
        if len(failed) > 50:
            print(f"| … and {len(failed) - 50} more |")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"_JUnit summary failed: {exc}_")
