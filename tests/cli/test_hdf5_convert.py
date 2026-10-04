"""Validate offline migration with independent h5py-produced input files."""

import importlib.util
from contextlib import closing
import os
from pathlib import Path
import shutil
import subprocess
import sys

import h5py
import numpy as np
import pytest


SCRIPT = Path(__file__).resolve().parents[2] / "src/tools/hdf5/convert.py"


@pytest.fixture
def converter():
    spec = importlib.util.spec_from_file_location("alps_hdf5_converter", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def marked_complex(parent, name, values, **options):
    dataset = parent.create_dataset(name, data=values, **options)
    dataset.attrs["__complex__"] = np.int8(1)
    return dataset


def compound_values(dataset, component_dtype, shape):
    """Read the actual compound fields, bypassing h5py's complex convenience."""
    dtype = np.dtype([("r", component_dtype), ("i", component_dtype)])
    value = np.empty(shape, dtype=dtype)
    with closing(h5py.h5t.py_create(dtype)) as memory_type:
        dataset.id.read(h5py.h5s.ALL, h5py.h5s.ALL, value, mtype=memory_type)
    with closing(dataset.id.get_type()) as stored:
        assert stored.get_class() == h5py.h5t.COMPOUND
        assert [stored.get_member_name(i) for i in range(stored.get_nmembers())] == [b"r", b"i"]
        for i in range(2):
            with closing(stored.get_member_type(i)) as member:
                assert member.get_size() == np.dtype(component_dtype).itemsize
                if np.dtype(component_dtype).byteorder in ("<", ">"):
                    assert member.get_order() == (
                        h5py.h5t.ORDER_LE if np.dtype(component_dtype).byteorder == "<"
                        else h5py.h5t.ORDER_BE
                    )
    return value


@pytest.mark.parametrize("dtype", ["<f4", ">f4", "<f8", ">f8"])
@pytest.mark.parametrize("scalar", [False, True])
def test_complex_preserves_precision_endian_and_ieee_components(converter, tmp_path, dtype, scalar):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    values = np.array([[0.0, -0.0], [np.inf, -np.inf], [np.nan, 1.25], [-4.5, np.nan]], dtype=dtype)
    if scalar:
        values = values[0]
    with h5py.File(source, "w") as archive:
        marked_complex(archive, "science/value", values)
    before = source.read_bytes()
    report = converter.convert(source, output)
    assert "/science/value: complex" in report
    assert source.read_bytes() == before
    with h5py.File(output, "r") as archive:
        dataset = archive["science/value"]
        assert dataset.shape == values.shape[:-1]
        actual = compound_values(dataset, dtype, values.shape[:-1])
        assert actual["r"].tobytes() == values[..., 0].tobytes()
        assert actual["i"].tobytes() == values[..., 1].tobytes()
        assert "__complex__" not in dataset.attrs


@pytest.mark.parametrize("dtype,values", [
    (">i8", [[-(2**63), 2**63 - 1], [2**53 + 1, -(2**53) - 3]]),
    ("<u8", [[2**64 - 1, 2**63 + 1], [2**53 + 1, 0]]),
])
def test_integral_complex_does_not_round_through_float(converter, tmp_path, dtype, values):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    values = np.array(values, dtype=dtype)
    with h5py.File(source, "w") as archive:
        marked_complex(archive, "value", values)
    converter.convert(source, output)
    with h5py.File(output, "r") as archive:
        actual = compound_values(archive["value"], dtype, values.shape[:-1])
        assert actual["r"].tobytes() == values[:, 0].tobytes()
        assert actual["i"].tobytes() == values[:, 1].tobytes()


def test_boolean_and_int8_markers_and_explicit_declarations(converter, tmp_path):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        for name, kind, values in [("mask", "bool", [0, 1, 0]), ("bytes", "int8", [-7, 0, 42])]:
            dataset = archive.create_dataset(name, data=np.array(values, dtype="i1"))
            dataset.attrs["__alps_type__"] = kind
        archive["legacy"] = np.array([0, 1, 0], dtype="i1")
        archive["unmarked"] = np.array([0, 1], dtype="i1")
        archive.attrs["legacy"] = np.int8(1)
        archive.attrs["byte"] = np.int8(-7)
        archive.attrs["__alps_type__:byte"] = "int8"
    report = converter.convert(source, output, boolean_datasets=["/legacy"],
                               boolean_attributes=[("/", "legacy")])
    assert any("/unmarked:" in line and "preserved as int8" in line for line in report)
    with h5py.File(output, "r") as archive:
        for name in ("mask", "legacy"):
            assert archive[name].dtype == np.dtype(bool)
            with closing(archive[name].id.get_type()) as datatype:
                assert datatype.get_class() == h5py.h5t.ENUM
            np.testing.assert_array_equal(np.arange(3)[archive[name][...]], [1])
        np.testing.assert_array_equal(archive["bytes"], [-7, 0, 42])
        assert archive["bytes"].dtype == np.dtype("i1")
        assert archive["unmarked"].dtype == np.dtype("i1")
        assert isinstance(archive.attrs["legacy"], np.bool_) and archive.attrs["legacy"]
        assert archive.attrs["byte"] == np.int8(-7)
        assert not any(name.startswith("__alps_type__") for name in archive.attrs)
        assert "__alps_type__" not in archive["bytes"].attrs


@pytest.mark.parametrize("parent", ["/", "/group", "/dataset"])
def test_attributes_null_empty_and_raw_strings(converter, tmp_path, parent):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        archive.create_group("group")
        archive["dataset"] = 42
        obj = archive[parent]
        obj.attrs["z"] = np.array([[1.25, -0.0], [np.inf, -np.inf]], dtype="f4")
        obj.attrs["__complex__:z"] = np.int8(1)
        obj.attrs["scalar"] = np.array([1.25, -0.0], dtype=">f4")
        obj.attrs["__complex__:scalar"] = np.int8(1)
        obj.attrs.create("null_z", h5py.Empty("f4"))
        obj.attrs["__complex__:null_z"] = np.int8(1)
        obj.attrs.create("null_flag", h5py.Empty("i1"))
        obj.attrs["__alps_type__:null_flag"] = "bool"
        obj.attrs["flag"] = np.array([0, 1], dtype="i1")
        obj.attrs["__alps_type__:flag"] = "bool"
        obj.attrs.create("null", h5py.Empty("f8"))
        obj.attrs.create("empty", np.empty((2, 0), dtype="f8"))
        obj.attrs.create("fixed", np.bytes_(b"a\x00b"), dtype="S3")
        obj.attrs.create("bytes", b"\xff", dtype=h5py.string_dtype("ascii"))
        archive.create_dataset("null", data=h5py.Empty("f8"))
        archive.create_dataset("null_complex", data=h5py.Empty("f4")).attrs["__complex__"] = np.int8(1)
        archive.create_dataset("null_bool", data=h5py.Empty("i1")).attrs["__alps_type__"] = "bool"
        marked_complex(archive, "empty_complex", np.empty((2, 0, 2), dtype="f4"))
        archive.create_dataset("empty", data=np.empty((2, 0), dtype="u8"))
        archive.create_dataset("strings", data=[b"first", b"", b"\xff"], dtype=h5py.string_dtype("ascii"))
        archive.create_dataset("fixed", data=np.array([b"a\x00b", b"last"], dtype="S4"))
        archive["ragged/0"] = np.array([1, 2], dtype="i8")
        archive["ragged/1"] = np.array([3], dtype="i8")
        archive["map/a&#47;b"] = 9
    converter.convert(source, output)
    with h5py.File(output, "r") as archive:
        obj = archive[parent]
        np.testing.assert_array_equal(obj.attrs["z"], np.array([1.25 - 0j, complex(np.inf, -np.inf)], dtype="c8"))
        assert np.signbit(obj.attrs["z"].imag[0])
        assert isinstance(obj.attrs["scalar"], np.complex64)
        assert obj.attrs["scalar"].real == 1.25 and np.signbit(obj.attrs["scalar"].imag)
        assert isinstance(obj.attrs["null_z"], h5py.Empty) and obj.attrs["null_z"].dtype.kind == "c"
        assert isinstance(obj.attrs["null_flag"], h5py.Empty) and obj.attrs["null_flag"].dtype == np.dtype(bool)
        assert obj.attrs["flag"].dtype == np.dtype(bool)
        assert isinstance(obj.attrs["null"], h5py.Empty)
        assert obj.attrs["empty"].shape == (2, 0)
        assert obj.attrs["fixed"].tobytes() == b"a\x00b"
        with closing(obj.attrs.get_id("bytes")) as attr:
            raw = np.empty((), attr.dtype)
            attr.read(raw)
            assert raw.item() == b"\xff"
        assert archive["null"].shape is None
        assert archive["null_complex"].shape is None
        with closing(archive["null_complex"].id.get_type()) as datatype:
            assert datatype.get_class() == h5py.h5t.COMPOUND
        assert archive["null_bool"].shape is None and archive["null_bool"].dtype == np.dtype(bool)
        assert archive["empty_complex"].shape == (2, 0)
        assert archive["empty"].shape == (2, 0) and archive["empty"].dtype == np.dtype("u8")
        np.testing.assert_array_equal(archive["strings"], [b"first", b"", b"\xff"])
        assert archive["fixed"][...].tobytes() == b"a\x00b\x00last"
        assert set(archive["ragged"]) == {"0", "1"}
        np.testing.assert_array_equal(archive["ragged/0"], [1, 2])
        assert archive["map/a&#47;b"][()] == 9
        assert not any(name.startswith("__complex__") or name.startswith("__alps_type__") for name in obj.attrs)


def test_aliases_cycles_and_soft_links_preserve_graph(converter, tmp_path):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        group = archive.create_group("science")
        values = marked_complex(group, "values", np.array([[1., 2.], [3., 4.]]))
        archive["alias"] = values
        archive["group_alias"] = group
        group["self"] = group
        group["root"] = archive["/"]
        archive["soft"] = h5py.SoftLink("/science/values")
        archive["dangling"] = h5py.SoftLink("/missing")
        archive["soft_cycle"] = h5py.SoftLink("/soft_cycle")
    converter.convert(source, output)
    with h5py.File(output, "r") as archive:
        address = lambda path: h5py.h5o.get_info(archive[path].id).addr
        assert address("alias") == address("science/values")
        assert address("group_alias") == address("science") == address("science/self")
        assert address("/") == address("science/root")
        for name, target in [("soft", "/science/values"), ("dangling", "/missing"), ("soft_cycle", "/soft_cycle")]:
            link = archive.get(name, getlink=True)
            assert isinstance(link, h5py.SoftLink) and link.path == target
        np.testing.assert_array_equal(archive["soft"], [1 + 2j, 3 + 4j])


@pytest.mark.parametrize("parent", ["/", "/group", "/dataset"])
@pytest.mark.parametrize("kind", ["utf8_spacepad", "bitfield"])
def test_unchanged_attribute_preserves_exact_hdf5_type(converter, tmp_path, parent, kind):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        archive.create_group("group")
        archive["dataset"] = 42
        with closing(h5py.h5s.create_simple((2,))) as space:
            with closing(h5py.h5t.C_S1.copy() if kind == "utf8_spacepad" else h5py.h5t.STD_B8LE.copy()) as datatype:
                if kind == "utf8_spacepad":
                    datatype.set_size(8)
                    datatype.set_cset(h5py.h5t.CSET_UTF8)
                    datatype.set_strpad(h5py.h5t.STR_SPACEPAD)
                    values = np.array([b"\xce\xbb      ", b"abc     "], dtype="S8")
                else:
                    values = np.array([0xA5, 0x5A], dtype="u1")
                with closing(h5py.h5a.create(archive[parent].id, b"metadata", datatype, space)) as attribute:
                    attribute.write(values, mtype=datatype)
    converter.convert(source, output)
    with h5py.File(source, "r") as original, h5py.File(output, "r") as migrated:
        with closing(original[parent].attrs.get_id("metadata")) as before:
            with closing(migrated[parent].attrs.get_id("metadata")) as after:
                with closing(before.get_type()) as input_type, closing(after.get_type()) as output_type:
                    assert input_type.equal(output_type)
                    if kind == "utf8_spacepad":
                        assert output_type.get_cset() == h5py.h5t.CSET_UTF8
                        assert output_type.get_strpad() == h5py.h5t.STR_SPACEPAD
                    else:
                        assert output_type.get_class() == h5py.h5t.BITFIELD
                    actual = np.empty_like(values)
                    after.read(actual, mtype=output_type)
                    assert actual.tobytes() == values.tobytes()


@pytest.mark.parametrize("case", [
    "complex_false", "complex_array_marker", "complex_bad_shape", "complex_string",
    "unknown_type", "type_array_marker", "bool_wrong_dtype", "bool_bad_value",
    "conflicting_markers", "orphan_attribute_marker", "group_marker",
    "explicit_bad_value", "explicit_int8_conflict", "complex_extensible_axis",
    "lossy_filter", "missing_declaration", "attribute_marker_declaration",
])
def test_invalid_input_never_publishes_or_changes_source(converter, tmp_path, case):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    arguments = {}
    with h5py.File(source, "w") as archive:
        archive["already_copied"] = np.arange(20)
        if case == "complex_extensible_axis":
            dataset = marked_complex(archive, "bad", np.ones((3, 2)), chunks=(2, 2), maxshape=(None, None))
        elif case == "lossy_filter":
            dataset = marked_complex(archive, "bad", np.ones((3, 2)), scaleoffset=2)
        else:
            data = [b"a", b"b"] if case == "complex_string" else [1, 2, 3] if case == "complex_bad_shape" else [0, 2] if case in ("bool_bad_value", "explicit_bad_value") else [0, 1]
            dtype = "S1" if case == "complex_string" else "f8" if case == "bool_wrong_dtype" else "i1"
            dataset = archive.create_dataset("bad", data=data, dtype=dtype)
        if case.startswith("complex_") or case == "lossy_filter":
            dataset.attrs["__complex__"] = np.int8(0) if case == "complex_false" else np.array([1], dtype="i1") if case == "complex_array_marker" else np.int8(1)
        elif case == "unknown_type":
            dataset.attrs["__alps_type__"] = "unknown"
        elif case == "type_array_marker":
            dataset.attrs["__alps_type__"] = ["bool"]
        elif case.startswith("bool_"):
            dataset.attrs["__alps_type__"] = "bool"
        elif case == "conflicting_markers":
            dataset.attrs["__complex__"] = np.int8(1)
            dataset.attrs["__alps_type__"] = "bool"
        elif case == "orphan_attribute_marker":
            archive.attrs["__complex__:missing"] = np.int8(1)
        elif case == "group_marker":
            archive.attrs["__complex__"] = np.int8(1)
        elif case in ("explicit_bad_value", "explicit_int8_conflict"):
            arguments["boolean_datasets"] = ["/bad"]
            if case == "explicit_int8_conflict":
                dataset.attrs["__alps_type__"] = "int8"
        elif case == "missing_declaration":
            arguments["boolean_datasets"] = ["/missing"]
        elif case == "attribute_marker_declaration":
            dataset.attrs["__alps_type__"] = "bool"
            arguments["boolean_attributes"] = [("/bad", "__alps_type__")]
    before = source.read_bytes()
    with pytest.raises((ValueError, KeyError, TypeError)):
        converter.convert(source, output, **arguments)
    assert source.read_bytes() == before
    assert not output.exists()
    assert sorted(path.name for path in tmp_path.iterdir()) == [source.name]


@pytest.mark.parametrize("case", ["object_reference", "region_reference", "attribute_reference", "nested_reference", "external_link", "external_storage", "virtual"])
def test_references_and_external_storage_are_rejected(converter, tmp_path, case):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        dataset = archive.create_dataset("values", data=np.arange(4))
        if case == "object_reference":
            archive.create_dataset("ref", data=[dataset.ref], dtype=h5py.ref_dtype)
        elif case == "region_reference":
            archive.create_dataset("ref", data=[dataset.regionref[1:3]], dtype=h5py.regionref_dtype)
        elif case == "attribute_reference":
            dataset.attrs["ref"] = dataset.ref
        elif case == "nested_reference":
            value = np.empty(1, dtype=np.dtype([("payload", "i4"), ("ref", h5py.ref_dtype)]))
            value["payload"], value["ref"] = 7, dataset.ref
            archive.create_dataset("nested", data=value)
        elif case == "external_link":
            archive["external"] = h5py.ExternalLink("missing.h5", "/values")
        elif case == "external_storage":
            archive.create_dataset("external", shape=(4,), dtype="i4", external=[(str(tmp_path / "raw.bin"), 0, h5py.h5f.UNLIMITED)])
        elif case == "virtual":
            layout = h5py.VirtualLayout(shape=(4,), dtype="i8")
            layout[:] = h5py.VirtualSource(str(source), "/values", shape=(4,))
            archive.create_virtual_dataset("virtual", layout)
    before = source.read_bytes()
    with pytest.raises(ValueError, match="unsupported"):
        converter.convert(source, output)
    assert source.read_bytes() == before
    assert not output.exists()
    assert not list(tmp_path.glob(".*.tmp"))


@pytest.mark.parametrize("kind", ["bool", "complex"])
def test_compressed_extensible_arrays_are_streamed_in_bounded_tiles(converter, tmp_path, monkeypatch, kind):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    shape = (3, 37)
    values = (np.arange(np.prod(shape)).reshape(shape) % 2).astype("i1") if kind == "bool" else np.arange(np.prod(shape) * 2, dtype="f8").reshape(shape + (2,))
    chunks = (1, 7) if kind == "bool" else (1, 7, 2)
    maxshape = (None, 37) if kind == "bool" else (None, 37, 2)
    with h5py.File(source, "w") as archive:
        dataset = archive.create_dataset("values", data=values, chunks=chunks, maxshape=maxshape,
                                         compression="gzip", compression_opts=4, shuffle=True, fletcher32=True)
        dataset.attrs["__alps_type__" if kind == "bool" else "__complex__"] = "bool" if kind == "bool" else np.int8(1)
    monkeypatch.setattr(converter, "BUFFER_BYTES", 64)
    reads = []
    original = h5py.Dataset.__getitem__
    def recorded(dataset, selection):
        value = original(dataset, selection)
        if Path(dataset.file.filename) == source:
            reads.append(np.asarray(value).nbytes)
        return value
    monkeypatch.setattr(h5py.Dataset, "__getitem__", recorded)
    converter.convert(source, output)
    assert len(reads) > 1 and max(reads) <= 64
    with h5py.File(output, "r") as archive:
        dataset = archive["values"]
        assert dataset.shape == shape and dataset.maxshape == (None, 37)
        assert dataset.chunks == (1, 7)
        assert dataset.compression == "gzip" and dataset.compression_opts == 4
        assert dataset.shuffle and dataset.fletcher32
        if kind == "bool":
            np.testing.assert_array_equal(dataset[...], values.astype(bool))
        else:
            actual = compound_values(dataset, "f8", shape)
            np.testing.assert_array_equal(actual["r"], values[..., 0])
            np.testing.assert_array_equal(actual["i"], values[..., 1])


def run_cli(*arguments):
    return subprocess.run([sys.executable, str(SCRIPT), *map(str, arguments)], capture_output=True, text=True)


def test_cli_explicit_boolean_and_existing_destination_protection(tmp_path):
    source, output = tmp_path / "source with spaces.h5", tmp_path / "output with spaces.h5"
    with h5py.File(source, "w") as archive:
        archive["mask"] = np.array([0, 1], dtype="i1")
        archive.attrs["flag"] = np.int8(1)
    before = source.read_bytes()
    result = run_cli(source, output, "--boolean", "/mask", "--boolean-attribute", "/", "flag")
    assert result.returncode == 0, result.stderr
    assert "Converted" in result.stdout
    with h5py.File(output, "r") as archive:
        assert archive["mask"].dtype == np.dtype(bool)
        assert isinstance(archive.attrs["flag"], np.bool_) and archive.attrs["flag"]
    output_before = output.read_bytes()
    result = run_cli(source, output)
    assert result.returncode == 1 and "already exists" in result.stderr
    assert result.stdout == ""
    assert source.read_bytes() == before and output.read_bytes() == output_before
    assert not list(tmp_path.glob(".*.tmp"))


@pytest.mark.parametrize("destination_kind", ["source", "hard_link", "dangling_symlink"])
def test_existing_alias_destinations_are_never_replaced(converter, tmp_path, destination_kind):
    source = tmp_path / "source.h5"
    with h5py.File(source, "w") as archive:
        archive["value"] = 42
    output = source if destination_kind == "source" else tmp_path / "output.h5"
    if destination_kind == "hard_link":
        output.hardlink_to(source)
    elif destination_kind == "dangling_symlink":
        try:
            output.symlink_to(tmp_path / "missing.h5")
        except OSError:
            pytest.skip("file symlinks unavailable")
    before = source.read_bytes()
    with pytest.raises(FileExistsError):
        converter.convert(source, output)
    assert source.read_bytes() == before
    if destination_kind == "dangling_symlink":
        assert output.is_symlink()


def test_installed_script_runs_after_standalone_relocation(tmp_path):
    sdk = os.environ.get("ALPS_DIR")
    if not sdk:
        pytest.skip("requires an installed SDK in ALPS_DIR")
    installed = Path(sdk).resolve().parents[1] / "bin" / "alps-hdf5-convert"
    assert installed.is_file(), "Install the SDK tools before running this integration test"
    directory = tmp_path / "standalone tools with spaces"
    directory.mkdir()
    script = directory / "alps-hdf5-convert"
    shutil.copy2(installed, script)
    source, output = directory / "input.h5", directory / "output.h5"
    with h5py.File(source, "w") as archive:
        marked_complex(archive, "science/value", np.array([[1.25, -0.0], [3., 4.]]))
        legacy_params(archive)
        pair = archive.create_group("pair")
        pair["first"], pair["second"] = 1, 2
        padded_matrix(archive, "matrix", 2, 3, 4)
    result = subprocess.run([sys.executable, "-I", str(script), str(source), str(output),
                             "--pair", "/pair", "--matrix", "/matrix"],
                            cwd=directory, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    with h5py.File(output, "r") as archive:
        actual = compound_values(archive["science/value"], "f8", (2,))
        np.testing.assert_array_equal(actual["r"], [1.25, 3.])
        np.testing.assert_array_equal(actual["i"], [-0.0, 4.])
        assert np.signbit(actual["i"][0])
        assert archive["parameters/format"].asstr()[()] == "alps.params.v2"
        assert set(archive["pair"]) == {"0", "1"}
        assert archive["matrix"].shape == (3, 2)


def legacy_params(archive, path="parameters"):
    """An independent v1 checkpoint, including historically ambiguous NULLs."""
    group = archive.create_group(path)
    group["format"] = "alps.params.v1"
    entries = group.create_group("entries")
    examples = [
        ("flag", "bool", np.int8(1)),
        ("wide", "int64", np.int64(2**53 + 1)),
        ("unsigned", "uint64", np.uint64(2**64 - 1)),
        ("real", "float64", np.float64(-0.0)),
        ("complex", "complex128", np.array([np.inf, -0.0])),
        ("label/with.dots", "string", "a,b"),
        ("flags", "bool[]", np.array([0, 1], dtype="i1")),
        ("ints", "int64[]", np.array([-7, 2**53 + 1], dtype="i8")),
        ("uints", "uint64[]", np.array([2**64 - 1], dtype="u8")),
        ("reals", "float64[]", np.array([np.nan, -np.inf])),
        ("complexes", "complex128[]", np.array([[1.25, -2.5], [0.0, -0.0]])),
        ("labels", "string[]", np.array(["", "x/y"], dtype=h5py.string_dtype())),
    ]
    for name, logical, values in examples:
        entry = entries.create_group(str(len(entries)))
        entry["name"], entry["type"] = name, logical
        dataset = entry.create_dataset("value", data=values)
        if logical.startswith("complex"):
            dataset.attrs["__complex__"] = np.int8(1)
    for logical, dtype in [("bool", "i1"), ("int64", "i8"), ("uint64", "u8"),
                           ("float64", "f8"), ("complex128", "f8"),
                           ("string", h5py.string_dtype())]:
        entry = entries.create_group(str(len(entries)))
        entry["name"], entry["type"] = "empty " + logical, logical + "[]"
        dataset = entry.create_dataset("value", data=h5py.Empty(dtype))
        if logical == "complex128":
            dataset.attrs["__complex__"] = np.int8(1)
    group.attrs["scientific-note"] = "keep this"
    return group


def test_params_v1_migrates_declared_types_and_recovers_empty_vector_rank(converter, tmp_path):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        group = legacy_params(archive, "custom/checkpoint")
        # A hard link encountered before the checkpoint group must see the same
        # logical-type conversion and preserve object identity.
        archive["alias-mask"] = group["entries/6/value"]
        archive.create_group("empty-dictionary")["format"] = "alps.params.v1"
        archive["empty-dictionary"].create_group("entries")
        archive["generic-null"] = h5py.Empty("f8")
    before = source.read_bytes()
    report = converter.convert(source, output)
    assert source.read_bytes() == before
    assert "/custom/checkpoint: alps.params.v1 -> alps.params.v2" in report
    with h5py.File(output, "r") as archive:
        group = archive["custom/checkpoint"]
        assert group["format"].asstr()[()] == "alps.params.v2"
        assert group.attrs["scientific-note"] == "keep this"
        assert archive["empty-dictionary/format"].asstr()[()] == "alps.params.v2"
        assert len(archive["empty-dictionary/entries"]) == 0
        entries = {entry["name"].asstr()[()]: entry for entry in group["entries"].values()}
        for entry in entries.values():
            logical = entry["type"].asstr()[()]
            value = entry["value"]
            assert len(value.shape) == (1 if logical.endswith("[]") else 0)
            assert not any(name.startswith("__complex__") or name.startswith("__alps_type__")
                           for name in value.attrs)
        assert entries["flag"]["value"].dtype == np.dtype(bool)
        assert entries["flag"]["value"][()] == np.bool_(True)
        assert entries["wide"]["value"][()] == 2**53 + 1
        assert entries["unsigned"]["value"][()] == 2**64 - 1
        assert np.signbit(entries["real"]["value"][()])
        scalar = compound_values(entries["complex"]["value"], "f8", ())
        assert np.isposinf(scalar["r"]) and np.signbit(scalar["i"])
        np.testing.assert_array_equal(entries["flags"]["value"], [False, True])
        assert archive["alias-mask"].id == entries["flags"]["value"].id
        for name, entry in entries.items():
            if name.startswith("empty "):
                assert entry["value"].shape == (0,)
        assert archive["generic-null"].shape is None


@pytest.mark.parametrize("fault", ["unknown-type", "wrong-type", "wrong-rank", "null-scalar",
                                  "invalid-bool", "duplicate-name", "sparse-entries", "extra-field"])
def test_malformed_params_v1_fails_without_publishing(converter, tmp_path, fault):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        group = legacy_params(archive)
        entry = group["entries/0"]
        if fault == "unknown-type":
            del entry["type"]
            entry["type"] = "guess"
        elif fault == "wrong-type":
            del entry["value"]
            entry["value"] = np.float64(1)
        elif fault == "wrong-rank":
            del entry["value"]
            entry["value"] = np.array([0, 1], dtype="i1")
        elif fault == "null-scalar":
            del entry["value"]
            entry["value"] = h5py.Empty("i1")
        elif fault == "invalid-bool":
            entry["value"][()] = np.int8(2)
        elif fault == "duplicate-name":
            del group["entries/1/name"]
            group["entries/1/name"] = "flag"
        elif fault == "sparse-entries":
            group["entries"].move("1", "99")
        else:
            entry["unexpected"] = 1
    before = source.read_bytes()
    with pytest.raises(ValueError):
        converter.convert(source, output)
    assert source.read_bytes() == before
    assert not output.exists()
    assert set(tmp_path.iterdir()) == {source}


def test_params_converter_output_loads_in_native_sdk(converter, tmp_path):
    if not os.environ.get("ALPS_DIR"):
        pytest.skip("installed ALPS SDK is unavailable")
    from pyalps import hdf5, ngs
    source, output = tmp_path / "legacy-params.h5", tmp_path / "params.h5"
    with h5py.File(source, "w") as archive:
        legacy_params(archive)
    converter.convert(source, output)
    with hdf5.archive(str(output), "r") as archive:
        parameters = ngs.params(archive, "/parameters")
    assert parameters["wide"] == 2**53 + 1
    assert parameters["unsigned"] == 2**64 - 1
    assert parameters["flag"] is True
    np.testing.assert_array_equal(parameters["flags"], [False, True])
    for logical in ("bool", "int64", "uint64", "float64", "complex128", "string"):
        assert len(parameters["empty " + logical]) == 0


def test_params_old_ascii_charset_with_utf8_bytes_becomes_utf8_text(converter, tmp_path):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        group = archive.create_group("parameters")
        group["format"] = "alps.params.v1"
        entry = group.create_group("entries/0")
        dtype = h5py.string_dtype("ascii")
        entry.create_dataset("name", data="Δ/name".encode(), dtype=dtype)
        entry.create_dataset("type", data=b"string", dtype=dtype)
        entry.create_dataset("value", data="λ value".encode(), dtype=dtype)
    converter.convert(source, output)
    with h5py.File(output, "r") as archive:
        entry = archive["parameters/entries/0"]
        assert entry["name"].asstr()[()] == "Δ/name"
        assert entry["value"].asstr()[()] == "λ value"
        for name in ("name", "type", "value"):
            assert h5py.check_string_dtype(entry[name].dtype).encoding == "utf-8"


@pytest.mark.parametrize("fault", ["nul-string", "invalid-utf8", "raw-bool-enum", "format-alias"])
def test_params_preflight_rejects_unrepresentable_text_and_boolean_codes(converter, tmp_path, fault):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        group = legacy_params(archive)
        entry = group["entries/0"]
        if fault in ("nul-string", "invalid-utf8"):
            del entry["type"], entry["value"]
            entry["type"] = "string"
            entry["value"] = np.bytes_(b"abc\0def" if fault == "nul-string" else b"\xff")
        elif fault == "raw-bool-enum":
            del entry["value"]
            dataset = entry.create_dataset("value", shape=(), dtype=bool)
            with closing(dataset.id.get_type()) as datatype:
                dataset.id.write(h5py.h5s.ALL, h5py.h5s.ALL, np.array(2, dtype="u1"), mtype=datatype)
        else:
            del entry["type"], entry["value"]
            entry["type"] = "string"
            entry["value"] = group["format"]
    before = source.read_bytes()
    with pytest.raises((ValueError, UnicodeError)):
        converter.convert(source, output)
    assert source.read_bytes() == before
    assert not output.exists()
    assert set(tmp_path.iterdir()) == {source}


def padded_matrix(parent, name, rows, columns, stride, *, complex_values=False, **options):
    group = parent.create_group(name)
    for field, value in [("size1", rows), ("size2", columns), ("reserved_size1", stride)]:
        group[field] = np.int64(value)
    # Padding deliberately differs from the visible matrix values.
    values = np.arange(stride * columns, dtype="f8") + 0.125
    if complex_values:
        values = np.column_stack((values, -values))
    dataset = group.create_dataset("values", data=values, **options)
    if complex_values:
        dataset.attrs["__complex__"] = np.int8(1)
    return group


def test_explicit_pair_migration_preserves_hard_links_cycles_and_attributes(converter, tmp_path):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        pair = archive.create_group("pair")
        pair.attrs["label"] = "keep"
        pair["first"] = np.int64(42)
        pair["first"].attrs["unit"] = "energy"
        pair["second"] = pair
        archive["same-pair"] = pair
        archive["same-value"] = pair["first"]
        archive["ordinary"] = 17
        archive["safe-soft"] = h5py.SoftLink("/same-value")
    converter.convert(source, output, pair_groups=["/pair", "/same-pair"])
    with h5py.File(output, "r") as archive:
        assert set(archive["pair"]) == {"0", "1"}
        assert archive["pair"].id == archive["same-pair"].id
        assert archive["pair/0"].id == archive["same-value"].id == archive["safe-soft"].id
        assert archive["pair/1"].id == archive["pair"].id
        assert archive["pair"].attrs["label"] == "keep"
        assert archive["pair/0"].attrs["unit"] == "energy"


@pytest.mark.parametrize("complex_values", [False, True])
@pytest.mark.parametrize("rows,columns,stride", [(3, 2, 5), (0, 3, 4), (2, 0, 7), (0, 0, 0)])
def test_explicit_matrix_migration_removes_only_padding_and_preserves_orientation(
        converter, tmp_path, rows, columns, stride, complex_values):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        group = padded_matrix(archive, "matrix", rows, columns, stride, complex_values=complex_values)
        group.attrs["label"] = "keep group attribute"
        group["values"].attrs["unit"] = "keep value attribute"
        archive["soft-matrix"] = h5py.SoftLink("/matrix")
        archive["unselected/first"] = 3
        archive["unselected/second"] = 4
    before = source.read_bytes()
    converter.convert(source, output, matrix_groups=["/matrix"])
    assert source.read_bytes() == before
    with h5py.File(output, "r") as archive:
        matrix = archive["matrix"]
        assert isinstance(matrix, h5py.Dataset)
        assert matrix.shape == (columns, rows)
        assert archive["soft-matrix"].id == matrix.id
        assert matrix.attrs["label"] == "keep group attribute"
        assert matrix.attrs["unit"] == "keep value attribute"
        expected = (np.arange(stride * columns, dtype="f8") + .125).reshape(columns, stride)[:, :rows]
        if complex_values:
            values = compound_values(matrix, "f8", (columns, rows))
            np.testing.assert_array_equal(values["r"], expected)
            np.testing.assert_array_equal(values["i"], -expected)
        else:
            np.testing.assert_array_equal(matrix, expected)
        assert set(archive["unselected"]) == {"first", "second"}


def test_genuine_old_matrix_fixture_requires_explicit_selection(converter, tmp_path):
    source = Path(__file__).with_name("fixtures") / "legacy-alps-matrix.h5"
    output = tmp_path / "matrix.h5"
    before = source.read_bytes()
    converter.convert(source, output, matrix_groups=["/matrix_old_hdf5_format"])
    assert source.read_bytes() == before
    with h5py.File(source, "r") as old, h5py.File(output, "r") as new:
        group = old["matrix_old_hdf5_format"]
        rows, columns, stride = (int(group[name][()]) for name in
                                  ("size1", "size2", "reserved_size1"))
        expected = group["values"][...].reshape(columns, stride)[:, :rows]
        assert new["matrix_old_hdf5_format"].shape == (columns, rows)
        np.testing.assert_array_equal(new["matrix_old_hdf5_format"], expected)


def test_matrix_conversion_streams_compressed_padded_storage(converter, tmp_path, monkeypatch):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        padded_matrix(archive, "matrix", 7, 9, 11, compression="gzip", compression_opts=3,
                      shuffle=True, fletcher32=True, chunks=(22,))
    monkeypatch.setattr(converter, "BUFFER_BYTES", 64)
    converter.convert(source, output, matrix_groups=["/matrix"])
    with h5py.File(output, "r") as archive:
        dataset = archive["matrix"]
        assert dataset.compression == "gzip" and dataset.compression_opts == 3
        assert dataset.shuffle and dataset.fletcher32 and dataset.chunks is not None
        expected = (np.arange(99) + .125).reshape(9, 11)[:, :7]
        np.testing.assert_array_equal(dataset, expected)


@pytest.mark.parametrize("fault", ["extra-field", "wrong-rank", "negative", "short-storage", "small-stride",
                                  "group-alias", "value-alias", "metadata-attr", "attribute-conflict",
                                  "soft-field", "extensible-storage", "soft-to-metadata", "not-numeric"])
def test_matrix_schema_ambiguities_fail_without_publishing(converter, tmp_path, fault):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        group = padded_matrix(archive, "matrix", 2, 3, 4)
        if fault == "extra-field":
            group["extra"] = 1
        elif fault == "wrong-rank":
            del group["size1"]
            group["size1"] = [2]
        elif fault == "negative":
            group["size1"][()] = -1
        elif fault == "short-storage":
            del group["values"]
            group["values"] = np.arange(11.)
        elif fault == "small-stride":
            group["reserved_size1"][()] = 1
        elif fault == "group-alias":
            archive["alias"] = group
        elif fault == "value-alias":
            archive["alias"] = group["values"]
        elif fault == "metadata-attr":
            group["size1"].attrs["meaningful"] = "cannot discard"
        elif fault == "attribute-conflict":
            group.attrs["label"] = "one"
            group["values"].attrs["label"] = "different"
        elif fault == "soft-field":
            archive["outside"] = group["values"][...]
            del group["values"]
            group["values"] = h5py.SoftLink("/outside")
        elif fault == "extensible-storage":
            del group["values"]
            group.create_dataset("values", data=np.arange(12.), maxshape=(None,))
        elif fault == "soft-to-metadata":
            archive["linked-size"] = h5py.SoftLink("/matrix/size1")
        else:
            del group["values"]
            group["values"] = np.array([b"text"] * 12)
    before = source.read_bytes()
    with pytest.raises(ValueError):
        converter.convert(source, output, matrix_groups=["/matrix"])
    assert source.read_bytes() == before
    assert not output.exists()
    assert set(tmp_path.iterdir()) == {source}


@pytest.mark.parametrize("fault", ["extra-field", "soft-target", "dangling-becomes-valid", "soft-selection",
                                  "overlap", "conflicting-selection"])
def test_pair_schema_ambiguities_fail_without_publishing(converter, tmp_path, fault):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        group = archive.create_group("pair")
        group["first"], group["second"] = 3, 4
        pairs, matrices = ["/pair"], []
        if fault == "extra-field":
            group["extra"] = 1
        elif fault == "soft-target":
            archive["linked-first"] = h5py.SoftLink("/pair/first")
        elif fault == "dangling-becomes-valid":
            archive["dangling"] = h5py.SoftLink("/pair/0")
        elif fault == "soft-selection":
            archive["alias"] = h5py.SoftLink("/pair")
            pairs = ["/alias"]
        elif fault == "overlap":
            del group["first"]
            child = group.create_group("first")
            child["first"], child["second"] = 1, 2
            pairs.append("/pair/first")
        else:
            matrices.append("/pair")
    before = source.read_bytes()
    with pytest.raises(ValueError):
        converter.convert(source, output, pair_groups=pairs, matrix_groups=matrices)
    assert source.read_bytes() == before
    assert not output.exists()
    assert set(tmp_path.iterdir()) == {source}


def test_cli_repeatable_pair_and_matrix_options(converter, tmp_path):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        for name in ("pair-a", "pair-b"):
            group = archive.create_group(name)
            group["first"], group["second"] = 1, 2
        for name in ("matrix-a", "matrix-b"):
            padded_matrix(archive, name, 2, 3, 4)
    assert converter.main([str(source), str(output), "--pair", "/pair-a", "--pair", "/pair-b",
                           "--matrix", "/matrix-a", "--matrix", "/matrix-b"]) == 0
    with h5py.File(output, "r") as archive:
        for name in ("pair-a", "pair-b"):
            assert set(archive[name]) == {"0", "1"}
        for name in ("matrix-a", "matrix-b"):
            assert archive[name].shape == (3, 2)


def test_params_text_conversion_preserves_layout_and_utf8_fill(converter, tmp_path):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        group = archive.create_group("parameters")
        group["format"] = "alps.params.v1"
        entry = group.create_group("entries/0")
        entry["name"], entry["type"] = "text", "string[]"
        entry.create_dataset("value", data=np.array([b"one"], dtype=object),
                             dtype=h5py.string_dtype("ascii"), chunks=(1,), maxshape=(None,),
                             compression="gzip", compression_opts=4, fillvalue=b"next")
    converter.convert(source, output)
    with h5py.File(output, "r") as archive:
        dataset = archive["parameters/entries/0/value"]
        assert dataset.chunks == (1,) and dataset.maxshape == (None,)
        assert dataset.compression == "gzip" and dataset.compression_opts == 4
        assert dataset.fillvalue == b"next"
        dataset_values = dataset.asstr()[...]
        np.testing.assert_array_equal(dataset_values, ["one"])


@pytest.mark.parametrize("field", ["int64", "uint64"])
def test_params_numeric_enum_is_not_relabelled_as_a_plain_integer(converter, tmp_path, field):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        group = archive.create_group("parameters")
        group["format"] = "alps.params.v1"
        entry = group.create_group("entries/0")
        entry["name"], entry["type"] = "x", field
        entry.create_dataset("value", data=1, dtype=h5py.enum_dtype({"ONE": 1},
                             basetype="i8" if field == "int64" else "u8"))
    before = source.read_bytes()
    with pytest.raises(ValueError, match="disagrees"):
        converter.convert(source, output)
    assert source.read_bytes() == before
    assert set(tmp_path.iterdir()) == {source}


@pytest.mark.parametrize("rows,columns,stride", [(0, 2, 0), (2, 0, 5)])
@pytest.mark.parametrize("complex_values", [False, True])
def test_declared_matrix_shape_recovers_null_empty_storage(converter, tmp_path, rows, columns,
                                                           stride, complex_values):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        group = padded_matrix(archive, "matrix", rows, columns, stride)
        del group["values"]
        values = group.create_dataset("values", data=h5py.Empty("f8"))
        if complex_values:
            values.attrs["__complex__"] = np.int8(1)
        archive["unresolved"] = h5py.Empty("f8")
    report = converter.convert(source, output, matrix_groups=["/matrix"])
    assert any("/unresolved: NULL dataspace preserved" in line for line in report)
    assert not any("/matrix/values: NULL dataspace preserved" in line for line in report)
    with h5py.File(output, "r") as archive:
        assert archive["matrix"].shape == (columns, rows)
        assert archive["unresolved"].shape is None


@pytest.mark.parametrize("fault", ["payload", "fill"])
def test_boolean_matrix_rejects_invalid_raw_enum_codes(converter, tmp_path, fault):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        group = padded_matrix(archive, "matrix", 1, 2, 1)
        del group["values"]
        if fault == "payload":
            dataset = group.create_dataset("values", shape=(2,), dtype=bool)
            with closing(dataset.id.get_type()) as datatype:
                dataset.id.write(h5py.h5s.ALL, h5py.h5s.ALL, np.array([0, 2], dtype="u1"),
                                 mtype=datatype)
        else:
            with closing(h5py.h5t.py_create(np.dtype(bool), logical=True)) as datatype, \
                    closing(h5py.h5s.create_simple((2,))) as space, \
                    closing(h5py.h5p.create(h5py.h5p.DATASET_CREATE)) as props:
                raw_fill = np.array(2, dtype="u1").view(bool)
                props.set_fill_value(raw_fill)
                with closing(h5py.h5d.create(group.id, b"values", datatype, space, dcpl=props)):
                    pass
    before = source.read_bytes()
    with pytest.raises(ValueError, match="other than 0 or 1"):
        converter.convert(source, output, matrix_groups=["/matrix"])
    assert source.read_bytes() == before
    assert set(tmp_path.iterdir()) == {source}
