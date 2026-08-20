#!/usr/bin/env python3
# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT

"""Split Python-independent ALPS payload out of repaired pyalps wheels.

cibuildwheel must build and test a complete wheel because auditwheel and
delocate need to see every native dependency together.  After repair, however,
the application programs, shared libraries, and XML data are byte-identical in
all CPython wheels for one platform tag.  Keeping those files in every ABI
wheel makes release artifacts several times larger than the installed package.

This script verifies that the runtime payloads really are identical, creates
one ``pyalps-runtime`` ``py3-none-<platform>`` wheel for each platform tag, and
rewrites the CPython wheels to depend on that exact runtime version.  The
runtime wheel owns disjoint paths below the existing ``pyalps`` package, so
installing both distributions recreates the original repaired-wheel layout
without changing any RPATH or loader name.
"""

from __future__ import annotations

import argparse
import base64
import copy
import csv
from dataclasses import dataclass
from email.parser import BytesParser
import hashlib
import io
from pathlib import Path
import re
import sys
import zipfile


RUNTIME_PREFIXES = (
    "pyalps/bin/",
    "pyalps/lib/",
    "pyalps/xml/",
    "pyalps/.dylibs/",
    "pyalps.libs/",
)
RUNTIME_REQUIREMENT = "pyalps-runtime"


@dataclass(frozen=True)
class WheelIdentity:
    path: Path
    filename_prefix: str
    python_tag: str
    abi_tag: str
    platform_tag: str


@dataclass
class WheelMember:
    info: zipfile.ZipInfo
    data: bytes


def _parse_wheel_filename(path: Path) -> WheelIdentity:
    if path.suffix != ".whl":
        raise ValueError(f"not a wheel filename: {path.name}")
    try:
        filename_prefix, python_tag, abi_tag, platform_tag = path.stem.rsplit("-", 3)
    except ValueError as exc:
        raise ValueError(f"invalid wheel filename: {path.name}") from exc
    if not re.match(r"^pyalps-[^-]+(?:-[^-]+)?$", filename_prefix):
        raise ValueError(f"expected a pyalps wheel, got: {path.name}")
    return WheelIdentity(
        path=path,
        filename_prefix=filename_prefix,
        python_tag=python_tag,
        abi_tag=abi_tag,
        platform_tag=platform_tag,
    )


def _read_wheel(path: Path) -> list[WheelMember]:
    members: list[WheelMember] = []
    seen: set[str] = set()
    with zipfile.ZipFile(path) as wheel:
        for info in wheel.infolist():
            if info.filename in seen:
                raise ValueError(f"{path.name} contains duplicate member {info.filename}")
            seen.add(info.filename)
            members.append(WheelMember(info=info, data=wheel.read(info)))
    return members


def _is_runtime_member(name: str) -> bool:
    return name.startswith(RUNTIME_PREFIXES)


def _one_member(members: list[WheelMember], suffix: str) -> WheelMember:
    matches = [member for member in members if member.info.filename.endswith(suffix)]
    if len(matches) != 1:
        names = [member.info.filename for member in matches]
        raise ValueError(f"expected exactly one *{suffix} member; found {names}")
    return matches[0]


def _metadata(members: list[WheelMember]) -> tuple[WheelMember, str, str]:
    member = _one_member(members, ".dist-info/METADATA")
    message = BytesParser().parsebytes(member.data, headersonly=True)
    name = message.get("Name", "")
    version = message.get("Version", "")
    if re.sub(r"[-_.]+", "-", name).lower() != "pyalps" or not version:
        raise ValueError(f"unexpected wheel metadata Name={name!r}, Version={version!r}")
    return member, name, version


def _runtime_fingerprint(members: list[WheelMember]) -> dict[str, tuple[str, int]]:
    fingerprint: dict[str, tuple[str, int]] = {}
    for member in members:
        if not _is_runtime_member(member.info.filename):
            continue
        mode = member.info.external_attr >> 16
        digest = hashlib.sha256(member.data).hexdigest()
        fingerprint[member.info.filename] = (digest, mode)
    required = ("pyalps/bin/", "pyalps/lib/", "pyalps/xml/")
    missing = [prefix for prefix in required if not any(n.startswith(prefix) for n in fingerprint)]
    if missing:
        raise ValueError(f"wheel is missing runtime payload under: {missing}")
    return fingerprint


def _assert_same_runtime(
    reference_path: Path,
    reference: dict[str, tuple[str, int]],
    candidate_path: Path,
    candidate: dict[str, tuple[str, int]],
) -> None:
    if reference == candidate:
        return
    missing = sorted(reference.keys() - candidate.keys())
    extra = sorted(candidate.keys() - reference.keys())
    changed = sorted(
        name for name in reference.keys() & candidate.keys() if reference[name] != candidate[name]
    )
    raise ValueError(
        "runtime payload differs between ABI wheels for one platform tag: "
        f"{reference_path.name} vs {candidate_path.name}; "
        f"missing={missing}, extra={extra}, changed={changed}"
    )


def _add_runtime_requirement(metadata: bytes, version: str) -> bytes:
    newline = b"\r\n" if b"\r\n" in metadata else b"\n"
    separator = newline + newline
    if separator in metadata:
        headers, body = metadata.split(separator, 1)
    else:
        headers, body = metadata.rstrip(b"\r\n"), b""
    requirement = f"Requires-Dist: {RUNTIME_REQUIREMENT}=={version}".encode()
    normalized_headers = headers.lower().replace(b"_", b"-")
    if b"requires-dist: pyalps-runtime" not in normalized_headers:
        headers += newline + requirement
    return headers + separator + body


def _hash_record_value(data: bytes) -> str:
    digest = base64.urlsafe_b64encode(hashlib.sha256(data).digest()).rstrip(b"=")
    return "sha256=" + digest.decode("ascii")


def _record_bytes(members: list[WheelMember], record_name: str) -> bytes:
    output = io.StringIO(newline="")
    writer = csv.writer(output, lineterminator="\n")
    for member in members:
        name = member.info.filename
        if name == record_name or member.info.is_dir():
            continue
        writer.writerow((name, _hash_record_value(member.data), len(member.data)))
    writer.writerow((record_name, "", ""))
    return output.getvalue().encode("utf-8")


def _cloned_info(source: zipfile.ZipInfo, filename: str | None = None) -> zipfile.ZipInfo:
    info = copy.copy(source)
    if filename is not None:
        info.filename = filename
        info.orig_filename = filename
    return info


def _new_file_info(filename: str, template: zipfile.ZipInfo) -> zipfile.ZipInfo:
    info = zipfile.ZipInfo(filename=filename, date_time=template.date_time)
    info.compress_type = zipfile.ZIP_DEFLATED
    info.create_system = 3
    info.external_attr = 0o100644 << 16
    return info


def _write_wheel(path: Path, members: list[WheelMember], record_name: str) -> None:
    record_template = _one_member(members, ".dist-info/RECORD")
    without_record = [m for m in members if m.info.filename != record_name]
    record = WheelMember(
        info=_cloned_info(record_template.info, record_name),
        data=_record_bytes(without_record, record_name),
    )
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as wheel:
        for member in [*without_record, record]:
            wheel.writestr(_cloned_info(member.info), member.data)


def _rewrite_main_wheel(
    identity: WheelIdentity,
    members: list[WheelMember],
    output_dir: Path,
) -> Path:
    metadata_member, _, version = _metadata(members)
    record_member = _one_member(members, ".dist-info/RECORD")
    rewritten: list[WheelMember] = []
    for member in members:
        name = member.info.filename
        if _is_runtime_member(name) or name == record_member.info.filename:
            continue
        data = member.data
        if name == metadata_member.info.filename:
            data = _add_runtime_requirement(data, version)
        rewritten.append(WheelMember(info=_cloned_info(member.info), data=data))
    rewritten.append(WheelMember(info=_cloned_info(record_member.info), data=b""))
    output = output_dir / identity.path.name
    _write_wheel(output, rewritten, record_member.info.filename)
    return output


def _runtime_metadata(version: str) -> bytes:
    return (
        "Metadata-Version: 2.1\n"
        "Name: pyalps-runtime\n"
        f"Version: {version}\n"
        "Summary: Platform runtime for the pyalps simulation package\n"
        "License: MIT\n"
        "Requires-Python: >=3.10\n"
        "Description-Content-Type: text/plain\n"
        "\n"
        "Shared ALPS applications, libraries, and XML data used by pyalps.\n"
        "\n"
    ).encode()


def _runtime_wheel_metadata(platform_tag: str) -> bytes:
    tags = "".join(f"Tag: py3-none-{tag}\n" for tag in platform_tag.split("."))
    return (
        "Wheel-Version: 1.0\n"
        "Generator: pyalps split_pyalps_runtime_wheels\n"
        "Root-Is-Purelib: false\n"
        f"{tags}"
        "\n"
    ).encode()


def _create_runtime_wheel(
    identity: WheelIdentity,
    members: list[WheelMember],
    output_dir: Path,
) -> Path:
    metadata_member, _, version = _metadata(members)
    record_member = _one_member(members, ".dist-info/RECORD")
    license_member = _one_member(members, ".dist-info/licenses/LICENSE.txt")
    runtime_dist_info = f"pyalps_runtime-{version}.dist-info"
    runtime_members = [
        WheelMember(info=_cloned_info(member.info), data=member.data)
        for member in members
        if _is_runtime_member(member.info.filename)
    ]
    runtime_members.extend(
        [
            WheelMember(
                info=_new_file_info(f"{runtime_dist_info}/METADATA", metadata_member.info),
                data=_runtime_metadata(version),
            ),
            WheelMember(
                info=_new_file_info(f"{runtime_dist_info}/WHEEL", metadata_member.info),
                data=_runtime_wheel_metadata(identity.platform_tag),
            ),
            WheelMember(
                info=_cloned_info(
                    license_member.info,
                    f"{runtime_dist_info}/licenses/LICENSE.txt",
                ),
                data=license_member.data,
            ),
            WheelMember(
                info=_cloned_info(record_member.info, f"{runtime_dist_info}/RECORD"),
                data=b"",
            ),
        ]
    )
    output = output_dir / (
        f"pyalps_runtime-{version}-py3-none-{identity.platform_tag}.whl"
    )
    _write_wheel(output, runtime_members, f"{runtime_dist_info}/RECORD")
    return output


def verify_record(path: Path) -> None:
    """Raise if a wheel's RECORD does not match its contents."""
    members = _read_wheel(path)
    record_member = _one_member(members, ".dist-info/RECORD")
    expected = _record_bytes(
        [m for m in members if m.info.filename != record_member.info.filename],
        record_member.info.filename,
    )
    if record_member.data != expected:
        raise ValueError(f"invalid RECORD in {path.name}")


def split_wheel_directory(input_dir: Path, output_dir: Path) -> list[Path]:
    input_dir = input_dir.resolve()
    output_dir = output_dir.resolve()
    if input_dir == output_dir:
        raise ValueError("input and output directories must be different")
    inputs = sorted(input_dir.glob("pyalps-*.whl"))
    if not inputs:
        raise ValueError(f"no pyalps wheels found in {input_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    stale = sorted(output_dir.glob("*.whl"))
    if stale:
        raise ValueError(f"output directory already contains wheels: {stale}")

    groups: dict[str, list[tuple[WheelIdentity, list[WheelMember]]]] = {}
    for path in inputs:
        identity = _parse_wheel_filename(path)
        groups.setdefault(identity.platform_tag, []).append((identity, _read_wheel(path)))

    outputs: list[Path] = []
    for platform_tag, wheels in sorted(groups.items()):
        reference_identity, reference_members = wheels[0]
        reference_fingerprint = _runtime_fingerprint(reference_members)
        reference_version = _metadata(reference_members)[2]
        for identity, members in wheels[1:]:
            version = _metadata(members)[2]
            if version != reference_version:
                raise ValueError(
                    f"version mismatch in {platform_tag} wheels: "
                    f"{reference_version} vs {version}"
                )
            _assert_same_runtime(
                reference_identity.path,
                reference_fingerprint,
                identity.path,
                _runtime_fingerprint(members),
            )

        outputs.append(
            _create_runtime_wheel(reference_identity, reference_members, output_dir)
        )
        for identity, members in wheels:
            outputs.append(_rewrite_main_wheel(identity, members, output_dir))

    for path in outputs:
        verify_record(path)
    return sorted(outputs)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path, help="directory of repaired pyalps wheels")
    parser.add_argument("output_dir", type=Path, help="empty output directory")
    args = parser.parse_args(argv)
    try:
        outputs = split_wheel_directory(args.input_dir, args.output_dir)
    except (OSError, ValueError, zipfile.BadZipFile) as exc:
        parser.error(str(exc))
    input_size = sum(path.stat().st_size for path in args.input_dir.glob("pyalps-*.whl"))
    output_size = sum(path.stat().st_size for path in outputs)
    print(f"wrote {len(outputs)} wheels to {args.output_dir}")
    print(f"wheel bytes: {input_size} -> {output_size} ({input_size - output_size} saved)")
    for path in outputs:
        print(f"  {path.name}: {path.stat().st_size}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
