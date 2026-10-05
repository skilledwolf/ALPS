"""Validate offline migration with independent h5py-produced input files."""

import importlib.util
import hashlib
import json
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


def core_alea_fixture():
    folder = Path(__file__).with_name("fixtures")
    fixture = folder / "alpscore-v2.3.3-alea.h5"
    metadata = json.loads(fixture.with_suffix(".json").read_text())
    assert hashlib.sha256(fixture.read_bytes()).hexdigest() == metadata["fixture_sha256"]
    selections = [(kind, "/results/" + name) for name, kind in metadata["selections"].items()]
    return fixture, selections


def test_released_core_alea_results_preserve_statistics_and_load_natively(converter, tmp_path):
    source, selections = core_alea_fixture()
    output = tmp_path / "converted.h5"
    converter.convert(source, output, core_alea_groups=selections)
    with h5py.File(output, "r") as archive:
        for kind, path in selections:
            assert archive[path].attrs["version"] == 1
            assert archive[path].attrs["kind"] == converter.CORE_ALEA_KINDS[kind]
        np.testing.assert_array_equal(archive["results/mean/mean/value"],
                                      [33., 33., sum(i % 7 for i in range(67)) / 67.])
        assert archive["results/variance/count"][()] == 67
        assert archive["results/batch/batch/count"][...].sum() == 67
        assert archive["results/batch/batch/sum"].shape == (16, 3)
        assert archive["results/circular/cov"].shape == (2, 2)
        assert archive["results/circular/cov"].dtype == np.dtype("complex128")
        assert archive["results/elliptic/cov"].shape == (2, 2, 2, 2)
        assert archive["results/elliptic/cov"].dtype == np.dtype("float64")
        for group in archive["results/autocorr/level"].values():
            assert group.attrs["version"] == 1 and group.attrs["kind"] == 2
    native = os.environ.get("ALPS_ALEA_MIGRATION_READER")
    if native:
        subprocess.run([native, str(output)], check=True)


@pytest.mark.parametrize("fault", ["unknown-kind", "wrong-kind", "wrong-components", "missing-count2",
                                  "negative-count", "already-versioned", "wrong-level-size",
                                  "zero-count2", "negative-count2", "nan-count2",
                                  "narrow-count", "narrow-covariance", "narrow-size", "enum-count"])
def test_core_alea_conversion_rejects_invalid_result_layouts(converter, tmp_path, fault):
    fixture, _ = core_alea_fixture()
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    shutil.copy2(fixture, source)
    kind, path = "covariance", "/results/covariance"
    with h5py.File(source, "a") as archive:
        group = archive[path]
        if fault == "unknown-kind": kind = "accumulator"
        elif fault == "wrong-kind": kind = "batch"
        elif fault == "wrong-components": group.attrs["size"] = np.uint64(4)
        elif fault == "missing-count2": del group["count2"]
        elif fault == "negative-count":
            del group["count"]
            group["count"] = np.int64(-1)
        elif fault == "already-versioned": group.attrs["version"] = np.uint64(1)
        elif fault == "narrow-count":
            del group["count"]
            group["count"] = np.int32(67)
        elif fault == "narrow-covariance":
            values = group["cov"][...].astype(np.float32)
            del group["cov"]
            group["cov"] = values
        elif fault == "narrow-size": group.attrs["size"] = np.uint32(3)
        elif fault == "enum-count":
            del group["count"]
            group.create_dataset("count", data=67, dtype=h5py.enum_dtype({"N": 67}, basetype="u8"))
        elif fault.endswith("count2"):
            group["count2"][()] = {"zero-count2": 0., "negative-count2": -1., "nan-count2": np.nan}[fault]
        else:
            kind, path = "autocorr", "/results/autocorr"
            archive[path + "/level/1"].attrs["size"] = np.uint64(4)
    before = source.read_bytes()
    with pytest.raises((ValueError, KeyError)):
        converter.convert(source, output, core_alea_groups=[(kind, path)])
    assert source.read_bytes() == before
    assert set(tmp_path.iterdir()) == {source}


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
        release_parameters(archive)
        pair = archive.create_group("pair")
        pair["first"], pair["second"] = 1, 2
        padded_matrix(archive, "matrix", 2, 3, 4)
    result = subprocess.run([sys.executable, "-I", str(script), str(source), str(output),
                             "--pair", "/pair", "--matrix", "/matrix",
                             "--parameters", "/parameters", "--boolean", "/parameters/flag"],
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


def release_parameters(archive, path="parameters"):
    """Independent inputs using the official flat native-leaf writer contract."""
    group = archive.create_group(path)
    group["flag"] = np.int8(1)
    group["N"] = np.int32(7)
    group["real"] = np.float64(-0.0)
    group["wide"] = np.int64(2**53 + 1)
    group["unsigned"] = np.uint64(2**64 - 1)
    group["ints"] = np.array([-7, 23], dtype="i4")
    group["empty-int"] = h5py.Empty("i4")
    group["empty-real"] = h5py.Empty("f8")
    group["empty-string"] = h5py.Empty(h5py.string_dtype("ascii"))
    group.create_dataset("EXPRESSION", data=b"sqrt(2) + unknown", dtype=h5py.string_dtype("ascii"))
    marked_complex(group, "complex", np.array([np.inf, -0.0]))
    group.attrs["note"] = "keep this"
    return group


def parameter_values(group):
    assert group["format"].asstr()[()] == "alps.params.v2"
    assert set(group) == {"format", "entries"}
    values = {}
    for entry in group["entries"].values():
        assert set(entry) == {"name", "value"}
        values[entry["name"].asstr()[()]] = entry["value"]
    return values


def test_released_flat_parameters_preserve_names_expressions_aliases_and_empty_types(converter, tmp_path):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        group = release_parameters(archive)
        # The official writers do NOT encode names. Literal escape text is a key.
        group["literal&#47;and&#38;name"] = np.int32(19)
        group["entries"], group["format"] = np.int32(1), "a legitimate parameter"
        archive["alias-value"] = group["N"]
        archive["alias-parameters"] = group
        archive["unresolved"] = h5py.Empty("f8")
    before = source.read_bytes()
    report = converter.convert(source, output, parameter_groups=["/parameters"],
                               boolean_datasets=["/parameters/flag"])
    assert source.read_bytes() == before
    assert any("ALPS 3.0.0 flat parameters" in line for line in report)
    with h5py.File(output, "r") as archive:
        values = parameter_values(archive["parameters"])
        assert values["N"].dtype == np.dtype("i8") and values["N"][()] == 7
        assert archive["alias-value"].id == values["N"].id
        assert archive["alias-parameters"].id == archive["parameters"].id
        assert values["flag"].dtype == np.dtype(bool) and values["flag"][()]
        assert values["wide"][()] == 2**53 + 1
        assert values["unsigned"][()] == 2**64 - 1
        assert np.signbit(values["real"][()])
        complex_value = compound_values(values["complex"], "f8", ())
        assert np.isposinf(complex_value["r"]) and np.signbit(complex_value["i"])
        assert values["literal&#47;and&#38;name"][()] == 19
        assert values["format"].asstr()[()] == "a legitimate parameter"
        assert values["entries"][()] == 1
        assert values["EXPRESSION"].asstr()[()] == "sqrt(2) + unknown"
        assert values["empty-int"].shape == (0,) and values["empty-int"].dtype == np.dtype("i8")
        assert values["empty-real"].shape == (0,)
        assert values["empty-string"].shape == (0,)
        assert h5py.check_string_dtype(values["empty-string"].dtype).encoding == "utf-8"
        assert archive["parameters"].attrs["note"] == "keep this"
        assert archive["unresolved"].shape is None


@pytest.mark.parametrize("fault", ["ambiguous-byte", "rank-two", "nested-key", "numeric-enum",
                                  "nul-string", "invalid-utf8", "raw-bool-enum", "soft-key"])
def test_malformed_released_parameters_leave_source_and_destination_untouched(converter, tmp_path, fault):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        group = release_parameters(archive)
        if fault in ("rank-two", "nested-key", "numeric-enum", "soft-key"):
            del group["N"]
        if fault == "rank-two":
            group["N"] = np.ones((2, 2))
        elif fault == "nested-key":
            group["N/key"] = 1
        elif fault == "numeric-enum":
            group.create_dataset("N", data=1, dtype=h5py.enum_dtype({"ONE": 1}, basetype="i8"))
        elif fault == "soft-key":
            group["N"] = h5py.SoftLink("/parameters/wide")
        elif fault in ("nul-string", "invalid-utf8"):
            group["bad"] = np.bytes_(b"abc\0def" if fault == "nul-string" else b"\xff")
        elif fault == "raw-bool-enum":
            del group["flag"]
            dataset = group.create_dataset("flag", shape=(), dtype=bool)
            with closing(dataset.id.get_type()) as datatype:
                dataset.id.write(h5py.h5s.ALL, h5py.h5s.ALL, np.array(2, dtype="u1"), mtype=datatype)
    before = source.read_bytes()
    with pytest.raises((ValueError, UnicodeError)):
        converter.convert(source, output, parameter_groups=["/parameters"],
                          boolean_datasets=[] if fault == "ambiguous-byte" else ["/parameters/flag"])
    assert source.read_bytes() == before
    assert set(tmp_path.iterdir()) == {source}


def test_released_parameter_text_uses_utf8_and_preserves_layout_fill_and_expression(converter, tmp_path):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        group = archive.create_group("parameters")
        group.create_dataset("Δ key", data=np.array([b"unknown + 1", "λ".encode()], dtype=object),
                             dtype=h5py.string_dtype("ascii"), chunks=(1,), maxshape=(None,),
                             compression="gzip", compression_opts=4, fillvalue="δ".encode())
    converter.convert(source, output, parameter_groups=["/parameters"])
    with h5py.File(output, "r") as archive:
        dataset = parameter_values(archive["parameters"])["Δ key"]
        assert dataset.chunks == (1,) and dataset.maxshape == (None,)
        assert dataset.compression == "gzip" and dataset.compression_opts == 4
        assert dataset.fillvalue == "δ".encode()
        np.testing.assert_array_equal(dataset.asstr()[...], ["unknown + 1", "λ"])
        assert h5py.check_string_dtype(dataset.dtype).encoding == "utf-8"


def test_official_v3_writer_fixture_profiles_are_pinned_and_preserve_science(converter, tmp_path):
    folder = Path(__file__).with_name("fixtures")
    fixture = folder / "alps-v3.0.0-profiles.h5"
    metadata = json.loads((folder / "alps-v3.0.0-profiles.json").read_text())
    assert metadata["source_revision"] == "1950cc6f682d7c4c1deae8b816f283857b1819d1"
    assert metadata["release"] == "ALPS v3.0.0"
    assert "Reconstructed" in metadata["provenance"]
    assert hashlib.sha256(fixture.read_bytes()).hexdigest() == metadata["sha256"]
    source, output = tmp_path / "released.h5", tmp_path / "converted.h5"
    shutil.copy2(fixture, source)
    before = source.read_bytes()
    report = converter.convert(source, output,
                               parameter_groups=["/ngs-parameters", "/legacy-parameters"],
                               boolean_datasets=["/ngs-parameters/ENABLED"],
                               alea_groups=["/scalar-mcdata", "/scalar-evaluator", "/scalar-observable", "/vector-observable"])
    assert source.read_bytes() == before
    assert any("unrecoverable-empty-vector-observable" in line and "NULL" in line for line in report)
    assert not any("/vector-observable/timeseries/data: NULL" in line for line in report)
    with h5py.File(output, "r") as archive:
        values = parameter_values(archive["ngs-parameters"])
        assert values["TITLE"].asstr()[()] == "Δ experiment"
        assert values["EXPRESSION"].asstr()[()] == "sqrt(2) + unknown"
        assert values["N"].dtype == np.dtype("i8") and values["EMPTY_INT"].shape == (0,)
        assert parameter_values(archive["legacy-parameters"])["EXPRESSION"].asstr()[()] == "2 * unresolved"
        scalar = archive["scalar-mcdata"]
        assert scalar["timeseries/data"].shape == scalar["jacknife/data"].shape == (0,)
        assert scalar.attrs["cannotrebin"] == np.bool_(False)
        assert scalar["timeseries/data"].attrs["binsize"] == 0
        evaluator = archive["scalar-evaluator"]
        assert evaluator["count"].dtype == np.dtype("f8") and evaluator["count"][()] == 6
        assert evaluator.attrs["changed"] == np.bool_(False)
        assert evaluator.attrs["nonlinearoperations"] == np.bool_(True)
        assert evaluator["timeseries/data2"].shape == (0,)
        np.testing.assert_array_equal(evaluator["timeseries/data"], [1., 2., 3.])
        labels = archive["scalar-observable/labels"]
        assert labels.shape == () and labels.asstr()[()] == "Δ energy"
        assert h5py.check_string_dtype(labels.dtype).encoding == "utf-8"
        vector = archive["vector-observable"]
        np.testing.assert_array_equal(vector["labels"].asstr()[...], ["x", "λ"])
        for field in ("data", "data2"):
            assert vector["timeseries/" + field].shape == (0, 2)
            assert vector["timeseries/" + field].dtype == np.dtype("f8")
        np.testing.assert_array_equal(vector["timeseries/logbinning"], [[1., 2.]])
        np.testing.assert_array_equal(vector["timeseries/partialbin"], [1., 2.])
        assert vector["timeseries/partialbin"].attrs["count"] == 1
        assert archive["unrecoverable-empty-vector-observable/timeseries/data"].shape is None


@pytest.mark.parametrize("fault", ["no-element-shape", "no-value-type", "conflicting-shape", "flagless",
                                  "invalid-flag", "vector-flag", "invalid-count", "bad-counts",
                                  "scalar-vector-label-mismatch", "rank-two-label", "numeric-label"])
def test_alea_profile_rejects_ambiguous_or_malformed_released_fields(converter, tmp_path, fault):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        group = archive.create_group("result")
        group["count"] = np.uint64(3)
        group.attrs["cannotrebin"] = np.int8(0)
        group["mean/value"] = np.array([1., 2.])
        group["mean/error"] = np.array([.1, .2])
        group["timeseries/data"] = np.array([[.5, 1.]])
        if fault == "no-element-shape":
            del group["mean"], group["timeseries"]
            group["timeseries/logbinning"] = h5py.Empty("i4")
        elif fault == "no-value-type":
            del group["timeseries/data"]
            group["timeseries/data"] = h5py.Empty("i4")
        elif fault == "conflicting-shape":
            del group["mean/error"]
            group["mean/error"] = np.ones(3)
        elif fault == "flagless":
            del group.attrs["cannotrebin"]
        elif fault == "invalid-flag":
            group.attrs["cannotrebin"] = np.int8(2)
        elif fault == "vector-flag":
            group.attrs["cannotrebin"] = np.array([0], dtype="i1")
        elif fault == "invalid-count":
            del group["count"]
            group["count"] = np.nan
        elif fault == "bad-counts":
            group["timeseries/logbinning_counts"] = np.ones((1, 2), dtype="u8")
        elif fault == "scalar-vector-label-mismatch":
            group.create_dataset("labels", data="scalar", dtype=h5py.string_dtype("utf-8"))
        elif fault == "rank-two-label":
            group.create_dataset("labels", data=[[b"x", b"y"]], dtype=h5py.string_dtype("ascii"))
        elif fault == "numeric-label":
            group["labels"] = np.array([1, 2])
    before = source.read_bytes()
    with pytest.raises(ValueError):
        converter.convert(source, output, alea_groups=["/result"])
    assert source.read_bytes() == before
    assert set(tmp_path.iterdir()) == {source}


def test_alea_null_bin_conversion_uses_value_exemplar_and_preserves_aliases(converter, tmp_path):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        group = archive.create_group("result")
        group["count"] = np.float64(2)
        group.attrs["changed"], group.attrs["nonlinearoperations"] = np.int8(0), np.int8(0)
        group["mean/value"] = np.array([1., 2.])
        group["mean/error"] = np.array([.1, .2])
        # Integer measurement bins and floating averages intentionally differ.
        group["timeseries/data"] = np.array([[1, 2], [1, 2]], dtype="i4")
        group["timeseries/data2"] = h5py.Empty("i4")
        group["jacknife/data"] = h5py.Empty("i4")
        archive["alias"] = group["timeseries/data2"]
    converter.convert(source, output, alea_groups=["/result"])
    with h5py.File(output, "r") as archive:
        assert archive["result/timeseries/data2"].shape == (0, 2)
        assert archive["result/timeseries/data2"].dtype == np.dtype("i4")
        assert archive["alias"].id == archive["result/timeseries/data2"].id
        assert archive["result/jacknife/data"].dtype == np.dtype("f8")
        np.testing.assert_array_equal(archive["result/timeseries/data"], [[1, 2], [1, 2]])


def test_release_profiles_are_explicit_and_repeatable_on_cli(converter, tmp_path):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        for name in ("one", "two"):
            group = archive.create_group(name)
            group["N"] = np.int32(7)
        for name in ("a", "b"):
            group = archive.create_group(name)
            group["count"] = np.uint64(1)
            group["mean/value"] = 1.
            group.attrs["cannotrebin"] = np.int8(0)
            group["timeseries/data"] = h5py.Empty("f8")
    assert converter.main([str(source), str(output), "--parameters", "/one", "--parameters", "/two",
                           "--alea", "/a", "--alea", "/b"]) == 0
    with h5py.File(output, "r") as archive:
        assert parameter_values(archive["one"])["N"][()] == 7
        assert parameter_values(archive["two"])["N"][()] == 7
        assert archive["a/timeseries/data"].shape == archive["b/timeseries/data"].shape == (0,)


def test_release_params_profile_output_loads_in_native_sdk(converter, tmp_path):
    if not os.environ.get("ALPS_DIR"):
        pytest.skip("installed ALPS SDK is unavailable")
    from pyalps import hdf5, ngs
    source, output = tmp_path / "released-params.h5", tmp_path / "params.h5"
    with h5py.File(source, "w") as archive:
        release_parameters(archive)
    converter.convert(source, output, parameter_groups=["/parameters"],
                      boolean_datasets=["/parameters/flag"])
    with hdf5.archive(str(output), "r") as archive:
        parameters = ngs.params(archive, "/parameters")
    assert parameters["wide"] == 2**53 + 1 and parameters["unsigned"] == 2**64 - 1
    assert parameters["flag"] is True
    assert parameters["EXPRESSION"] == "sqrt(2) + unknown"
    for name in ("empty-int", "empty-real", "empty-string"):
        assert len(parameters[name]) == 0


def test_alea_null_complex_vector_bins_use_known_partialbin_datatype(converter, tmp_path):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        group = archive.create_group("observable")
        group["count"] = np.uint64(1)
        values = np.array([[np.inf, -0.0], [3., -4.]])
        marked_complex(group, "mean/value", values)
        marked_complex(group, "timeseries/partialbin", values)
        marked_complex(group, "timeseries/logbinning", values[np.newaxis, ...])
        # Unlike the scalar typed-NULL case, vector<valarray<complex>> lost T.
        group["timeseries/data"] = h5py.Empty("i4")
    converter.convert(source, output, alea_groups=["/observable"])
    with h5py.File(output, "r") as archive:
        dataset = archive["observable/timeseries/data"]
        assert dataset.shape == (0, 2)
        compound_values(dataset, "f8", (0, 2))
        values = compound_values(archive["observable/timeseries/partialbin"], "f8", (2,))
        assert np.isposinf(values["r"][0]) and np.signbit(values["i"][0])


def test_null_declared_boolean_parameter_rejects_invalid_fill(converter, tmp_path):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        group = archive.create_group("parameters")
        group.create_dataset("flags", shape=None, dtype="i1", fillvalue=np.int8(2))
    before = source.read_bytes()
    with pytest.raises(ValueError, match="Boolean"):
        converter.convert(source, output, parameter_groups=["/parameters"],
                          boolean_datasets=["/parameters/flags"])
    assert source.read_bytes() == before
    assert set(tmp_path.iterdir()) == {source}


@pytest.mark.parametrize("link", ["/parameters/N", "/parameters/format"])
def test_parameter_schema_changes_cannot_retarget_soft_links(converter, tmp_path, link):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        archive["parameters/N"] = np.int32(7)
        archive["view"] = h5py.SoftLink(link)
    before = source.read_bytes()
    with pytest.raises(ValueError, match="changes a soft-link target"):
        converter.convert(source, output, parameter_groups=["/parameters"])
    assert source.read_bytes() == before
    assert set(tmp_path.iterdir()) == {source}


def test_released_fixture_profiles_load_through_native_params_and_result_readers(converter, tmp_path):
    if not os.environ.get("ALPS_DIR"):
        pytest.skip("installed ALPS SDK is unavailable")
    from pyalps import alea, hdf5, ngs
    source = Path(__file__).with_name("fixtures") / "alps-v3.0.0-profiles.h5"
    output = tmp_path / "converted.h5"
    converter.convert(source, output,
                      parameter_groups=["/ngs-parameters", "/legacy-parameters"],
                      boolean_datasets=["/ngs-parameters/ENABLED"],
                      alea_groups=["/scalar-mcdata", "/scalar-evaluator", "/scalar-observable", "/vector-observable"])
    # Each provider owns the file separately; no h5py/native overlap is needed.
    with hdf5.archive(str(output), "r") as archive:
        parameters = ngs.params(archive, "/ngs-parameters")
    assert parameters["N"] == 7 and parameters["ENABLED"] is True
    assert parameters["TITLE"] == "Δ experiment"
    assert parameters["EXPRESSION"] == "sqrt(2) + unknown"
    np.testing.assert_array_equal(parameters["VECTOR"], [1, 2, 3])
    for name in ("EMPTY_INT", "EMPTY_REAL", "EMPTY_TEXT"):
        assert len(parameters[name]) == 0
    scalar = alea.MCScalarData()
    scalar.load(str(output), "/scalar-mcdata")
    assert scalar.count == 0 and scalar.bins.size == 0
    with pytest.raises(RuntimeError, match="No measurements available"):
        _ = scalar.mean
    evaluator = alea.MCScalarData()
    evaluator.load(str(output), "/scalar-evaluator")
    assert evaluator.count == 6 and evaluator.mean == 1 and evaluator.error == .25
    # SimpleObservableData stores bin sums. The live domain reader supplies the
    # cross-schema sums-to-means rule; the offline converter kept [1,2,3] intact.
    np.testing.assert_array_equal(evaluator.bins, [.5, 1., 1.5])
    vector = alea.MCVectorData()
    vector.load(str(output), "/vector-observable")
    assert vector.count == 1 and vector.bins.size == 0
    np.testing.assert_array_equal(vector.mean, [1., 2.])


def test_parameter_profile_respects_an_explicit_integer_byte_marker(converter, tmp_path):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        dataset = archive.create_dataset("parameters/code", data=np.int8(1))
        dataset.attrs["__alps_type__"] = "int8"
    converter.convert(source, output, parameter_groups=["/parameters"])
    with h5py.File(output, "r") as archive:
        dataset = parameter_values(archive["parameters"])["code"]
        assert dataset.dtype == np.dtype("i8") and dataset[()] == 1
        assert "__alps_type__" not in dataset.attrs


def released_bins(archive, family, values, binsize=3, remainder=2):
    group = archive.create_group("observable")
    values = np.asarray(values)
    partial = np.asarray(values[0] + 1)
    group["count"] = np.uint64(len(values) * binsize + (remainder if family == "observable" else 0))
    data = group.create_dataset("timeseries/data", data=values)
    data.attrs.update(binsize=np.uint64(binsize), binningtype="linear", discard=np.uint32(0))
    if family == "observable":
        group["timeseries/logbinning"] = values[:1]
        group["timeseries/partialbin"] = partial
        group["timeseries/partialbin"].attrs["count"] = np.uint64(remainder)
        group["sum"] = values.sum(axis=0) + (partial if remainder else 0)
    elif family == "mcdata":
        group.attrs["cannotrebin"] = np.int8(0)
    else:
        group.attrs.update(changed=np.int8(0), nonlinearoperations=np.int8(0))
    group["mean/value"] = values.sum(axis=0) / max(1, int(group["count"][()]))
    group["mean/error"] = np.ones(values.shape[1:])
    group["user-note"] = "preserve me"
    return group


@pytest.mark.parametrize("family", ["observable", "evaluator", "mcdata"])
@pytest.mark.parametrize("components", [1, 2])
def test_released_linear_bins_recover_native_joint_analysis(converter, tmp_path, family, components):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    values = np.asarray([3., 9., 6.])
    if components == 2:
        values = np.column_stack((values, 2 * values))
    with h5py.File(source, "w") as archive:
        released_bins(archive, family, values)
    converter.convert(source, output, alea_batch_groups=["/observable"])
    sums = values * (3 if family == "mcdata" else 1)
    counts = np.full(3, 3, dtype="u8")
    if family == "observable":
        sums = np.concatenate((sums, np.asarray(values[:1] + 1)))
        counts = np.append(counts, 2).astype("u8")
    sums = sums.reshape(-1, components)
    average = sums.sum(axis=0) / counts.sum()
    weight = counts.astype(float)
    centered = (weight[:, None] * abs(sums / weight[:, None] - average)**2).sum(axis=0)
    error = np.sqrt(centered * (weight**2).sum() /
                    (weight.sum() * (weight.sum()**2 - (weight**2).sum())))
    with h5py.File(output, "r") as archive:
        group = archive["observable"]
        assert group.attrs["version"] == 1 and group.attrs["kind"] == 5
        assert group.attrs["size"] == components
        assert "timeseries" not in group and "count" not in group
        assert group["user-note"].asstr()[()] == "preserve me"
        np.testing.assert_array_equal(group["batch/count"], counts)
        np.testing.assert_array_equal(group["batch/sum"], sums)
        np.testing.assert_allclose(group["mean/value"], average)
        np.testing.assert_allclose(group["mean/error"], error)
    if os.environ.get("ALPS_DIR"):
        from pyalps import alea, hdf5
        with hdf5.archive(str(output), "r") as archive:
            native = alea.BatchResult.read(archive, "/observable")
        assert native.count == int(counts.sum())
        np.testing.assert_allclose(native.mean, average)
        np.testing.assert_allclose(native.error, error)
        if components == 2:
            np.testing.assert_allclose(native.covariance[0, 1], 2 * native.variance[0])


@pytest.mark.parametrize("fault", ["summary", "missing-tail", "discarded", "nonlinear", "changed",
                                  "nan", "wrong-total", "zero-bin", "alias", "softlink", "unknown-flags"])
def test_native_batch_conversion_rejects_unrecoverable_histories(converter, tmp_path, fault):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        group = released_bins(archive, "observable", [3., 9., 6.])
        if fault == "summary": del group["timeseries"]
        elif fault == "missing-tail": del group["timeseries/partialbin"]
        elif fault == "discarded": group["timeseries/data"].attrs["discard"] = np.uint32(1)
        elif fault in ("changed", "nonlinear"):
            del group["timeseries/logbinning"]
            group.attrs.update(changed=np.int8(fault == "changed"), nonlinearoperations=np.int8(fault == "nonlinear"))
        elif fault == "nan": group["timeseries/data"][0] = np.nan
        elif fault == "wrong-total": group["sum"][()] += 1
        elif fault == "zero-bin": group["timeseries/data"].attrs["binsize"] = np.uint64(0)
        elif fault == "alias": archive["alias"] = group["timeseries/data"]
        elif fault == "softlink": archive["view"] = h5py.SoftLink("/observable/count")
        elif fault == "unknown-flags": del group["timeseries/logbinning"]
    before = source.read_bytes()
    with pytest.raises(ValueError):
        converter.convert(source, output, alea_batch_groups=["/observable"])
    assert source.read_bytes() == before
    assert set(tmp_path.iterdir()) == {source}


def test_released_partial_vector_and_empty_mcdata_migrate_to_native_batches(converter, tmp_path):
    source = Path(__file__).with_name("fixtures") / "alps-v3.0.0-profiles.h5"
    output = tmp_path / "converted.h5"
    converter.convert(source, output, alea_batch_groups=["/vector-observable", "/scalar-mcdata"])
    with h5py.File(output, "r") as archive:
        group = archive["vector-observable"]
        np.testing.assert_array_equal(group["batch/count"], [1, 0])
        np.testing.assert_array_equal(group["batch/sum"], [[1., 2.], [0., 0.]])
        np.testing.assert_array_equal(group["mean/value"], [1., 2.])
        assert np.isinf(group["mean/error"][:]).all()
        group = archive["scalar-mcdata"]
        np.testing.assert_array_equal(group["batch/count"], [0, 0])
        assert np.isnan(group["mean/value"][:]).all()
        assert np.isnan(group["mean/error"][:]).all()


@pytest.mark.parametrize("family", ["evaluator", "mcdata"])
def test_released_result_final_row_retains_its_partial_weight(converter, tmp_path, family):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        group = released_bins(archive, family, [3., 9., 6.])
        group["count"][()] = np.uint64(8)
    converter.convert(source, output, alea_batch_groups=["/observable"])
    with h5py.File(output, "r") as archive:
        group = archive["observable"]
        np.testing.assert_array_equal(group["batch/count"], [3, 3, 2])
        factor = 3 if family == "mcdata" else 1
        np.testing.assert_array_equal(group["batch/sum"], np.array([[3.], [9.], [6.]]) * factor)
        np.testing.assert_array_equal(group["mean/value"], [18. * factor / 8])


def test_native_complex_bin_conversion_is_bounded_and_preserves_covariance(converter, tmp_path, monkeypatch):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    values = np.array([[3+2j, 6+4j], [9-2j, 18-4j], [6+1j, 12+2j]])
    with h5py.File(source, "w") as archive:
        released_bins(archive, "mcdata", values)
    monkeypatch.setattr(converter, "BUFFER_BYTES", 64)
    converter.convert(source, output, alea_batch_groups=["/observable"])
    sums = 3 * values
    with h5py.File(output, "r") as archive:
        group = archive["observable"]
        assert group["batch/sum"].dtype == np.dtype("c16")
        np.testing.assert_array_equal(group["batch/sum"], sums)
        np.testing.assert_allclose(group["mean/value"], sums.sum(axis=0) / 9)
        np.testing.assert_allclose(group["mean/error"], np.sqrt((abs(values-values.mean(axis=0))**2).sum(axis=0) / 6))
    if os.environ.get("ALPS_DIR"):
        from pyalps import alea, hdf5
        with hdf5.archive(str(output), "r") as archive:
            native = alea.ComplexBatchResult.read(archive, "/observable")
        np.testing.assert_array_equal(native.batch_sums, sums)
        np.testing.assert_allclose(native.covariance[0, 1], 2 * native.variance[0])


def test_native_conversion_rejects_nonempty_zero_weight_bins(converter, tmp_path):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    with h5py.File(source, "w") as archive:
        group = released_bins(archive, "evaluator", [1., 2.], binsize=0)
    with pytest.raises(ValueError, match="incomplete/inconsistent bin counts"):
        converter.convert(source, output, alea_batch_groups=["/observable"])
    assert not output.exists()


def test_released_summary_histogram_and_history_preservation(converter, tmp_path):
    # Reconstruct v3.0.0 histogram.h:234-250 and SimpleObservableData's
    # summary fields. A histogram is already ordinary HDF5; no new codec.
    source, output = tmp_path / "released.h5", tmp_path / "converted.h5"
    with h5py.File(source, "w") as archive:
        histogram = archive.create_group("histogram")
        histogram["histogram"] = np.array([2, 4, 1], dtype="u4")
        histogram["count"] = np.uint64(7)
        histogram.attrs.update(min=-1., max=2., stepsize=1.)
        summary = archive.create_group("summary")
        summary["count"] = np.uint64(103)
        summary["mean/value"] = np.array([1., 2.])
        summary["mean/error"] = np.array([.3, .6])
        summary["variance/value"] = np.array([2., 8.])
        summary["mean/error_convergence"] = np.array([0, 1], dtype="i4")
        summary.attrs.update(changed=np.int8(0), nonlinearoperations=np.int8(1))
        history = released_bins(archive, "observable", [3., 9., 6.])
        history["timeseries/logbinning2"] = np.array([18.])
        history["timeseries/logbinning_counts"] = np.array([3], dtype="u8")
        archive["checkpoint/version"] = np.int32(42)
        archive["checkpoint/random"] = "opaque released generator state"
    before = source.read_bytes()
    converter.convert(source, output, alea_groups=["/summary", "/observable"])
    assert source.read_bytes() == before
    with h5py.File(source) as old, h5py.File(output) as new:
        for group in ("histogram", "summary", "observable", "checkpoint"):
            def check(name, obj):
                if isinstance(obj, h5py.Dataset):
                    converted = new[group + "/" + name]
                    np.testing.assert_array_equal(obj[()], converted[()])
                    assert obj.dtype == converted.dtype
            old[group].visititems(check)
        for name in ("min", "max", "stepsize"):
            assert old["histogram"].attrs[name] == new["histogram"].attrs[name]
        assert new["summary"].attrs["nonlinearoperations"] == np.bool_(True)
        assert "kind" not in new["summary"].attrs
        assert "kind" not in new["histogram"].attrs
    # Preserved summary uncertainty cannot manufacture batches or a live cursor.
    rejected = tmp_path / "invented.h5"
    with pytest.raises(ValueError):
        converter.convert(source, rejected, alea_batch_groups=["/summary"])
    assert not rejected.exists()


def test_released_qwl_per_run_conversion(converter, tmp_path):
    source = Path(__file__).with_name("fixtures") / "alps-v3.0.0-qwl.h5"
    output = tmp_path / "native.h5"
    converter.convert(source, output, qwl_sites=4)
    root = "simulation/realizations/0/clones/0"
    with h5py.File(source) as old, h5py.File(output) as native:
        assert native["parameters/format"][()] == b"alps.params.v2"
        assert native["simulation/number_of_sites"][()] == 4
        assert native[root + "/complete"][()]
        for name in ("Coefficients", "Offset", "Histogram", "Fraction", "Uniform Structure Factor Coefficients"):
            new = native[root + "/results/" + name]
            assert new.attrs["version"] == 1 and new.attrs["kind"] == 1
            np.testing.assert_array_equal(new["mean/value"], np.atleast_1d(old[new.name + "/mean/value"][()]))
            assert new["count"].dtype == np.dtype("uint64")
        # Generic conversion preserves the original physical evidence.
        np.testing.assert_array_equal(native[root + "/results/Time Up/timeseries/data"],
                                      old[root + "/results/Time Up/timeseries/data"])
    executable = os.environ.get("ALPS_QWL_EXECUTABLE")
    if executable:
        evaluated = subprocess.run([str(Path(executable).with_name("qwl_evaluate")),
            "--T_MIN", "1", "--T_MAX", "1", str(output)], text=True, capture_output=True)
        assert evaluated.returncode == 0, evaluated.stderr
        assert (tmp_path / "native.plot.energy.xml").exists()


@pytest.mark.parametrize("fault", ["sites", "averaged", "incomplete", "aliases", "window", "shape"])
def test_qwl_conversion_rejects_missing_or_ambiguous_evidence(converter, tmp_path, fault):
    fixture = Path(__file__).with_name("fixtures") / "alps-v3.0.0-qwl.h5"
    source, output = tmp_path / "old.h5", tmp_path / "new.h5"
    shutil.copyfile(fixture, source)
    with h5py.File(source, "a") as ar:
        result = ar["simulation/realizations/0/clones/0/results"]
        if fault in ("averaged", "incomplete"):
            result["Coefficients/count"][()] = 2 if fault == "averaged" else 0
        elif fault == "aliases": ar["alias"] = result
        elif fault == "window": ar["parameters/EXPANSION_ORDER_MINIMUM"] = 3
        elif fault == "shape":
            del result["Offset/mean/value"]
            result["Offset/mean/value"] = [1., 2.]
    before = source.read_bytes()
    with pytest.raises(ValueError):
        converter.convert(source, output, qwl_sites=5 if fault == "sites" else 4)
    assert source.read_bytes() == before
    assert not output.exists()


def spinmc_pair(tmp_path, model):
    """HDF5 framework reconstruction + independently compiled XDR protocol fixture."""
    fixtures = Path(__file__).with_name('fixtures')
    stem = 'O4' if model == 'O(4)' else model
    xdr = fixtures / ('spinmc-' + stem + '.xdr')
    provenance = json.loads((fixtures / 'spinmc-state.json').read_text())
    assert hashlib.sha256(xdr.read_bytes()).hexdigest() == provenance['sha256'][xdr.name]
    source = tmp_path / 'worker.h5'
    with h5py.File(source, 'w') as ar:
        ar['parameters/MODEL'] = model
        if model == 'Potts':
            ar['parameters/q'] = 3
        ar['rng'] = 'opaque released RNG state'
        ar['rng'].attrs['name'] = 'mt19937'
        ar.create_group('simulation/realizations/0/clones/0/results')
    return source, xdr


@pytest.mark.parametrize('model', ['Ising', 'Potts', 'XY', 'Heisenberg', 'O(4)'])
def test_recover_released_spin_state(converter, tmp_path, model):
    source, xdr = spinmc_pair(tmp_path, model)
    original = source.read_bytes(), xdr.read_bytes()
    output = tmp_path / 'converted.h5'
    converter.convert(source, output, spinmc_state=xdr, parameter_groups=['/parameters'])
    with h5py.File(output, 'r') as ar:
        state = ar['migration/spinmc']
        assert state['model'].asstr()[()] == model
        assert state['sweeps_done'][()] == 0x100000005
        assert state['thermalization_sweeps'][()] == 0x100000003
        assert state['thermalization_fraction'][()] == 1.125
        assert bytes(state['source_xdr'][...]) == xdr.read_bytes()
        assert ar['rng'].asstr()[()] == 'opaque released RNG state'
        assert ar['rng'].attrs['name'] == 'mt19937'
        if model == 'Ising':
            expected = [[-1.], [1.]]
        elif model == 'Potts':
            expected = [[0.], [2.]]
        else:
            expected = np.zeros((2, {'XY':2, 'Heisenberg':3, 'O(4)':4}[model]))
            expected[0, 0] = 1.
            expected[1, 1] = -1.
        np.testing.assert_array_equal(state['spins'], expected)
        assert ar['parameters/format'].asstr()[()] == 'alps.params.v2'
    assert (source.read_bytes(), xdr.read_bytes()) == original


@pytest.mark.parametrize('fault', ['version', 'truncated', 'trailing', 'sites', 'boolean',
                                  'fraction', 'counters', 'vector', 'potts', 'rng', 'results'])
def test_spin_state_rejects_invalid_input_atomically(converter, tmp_path, fault):
    import struct
    model = 'XY' if fault == 'vector' else 'Potts' if fault == 'potts' else 'Ising'
    source, fixture = spinmc_pair(tmp_path, model)
    data = bytearray(fixture.read_bytes())
    if fault == 'version': struct.pack_into('>i', data, 8, 310)
    elif fault == 'truncated': del data[-1]
    elif fault == 'trailing': data.append(0)
    elif fault == 'sites': struct.pack_into('>I', data, 36, 0xffffffff)
    elif fault == 'boolean': struct.pack_into('>I', data, 40, 2)
    elif fault == 'fraction': struct.pack_into('>d', data, 20, float('nan'))
    elif fault == 'counters': struct.pack_into('>Q', data, 28, 0xffffffffffffffff)
    elif fault == 'vector': struct.pack_into('>d', data, 40, 2.)
    elif fault == 'potts': struct.pack_into('>I', data, 40, 3)
    else:
        with h5py.File(source, 'a') as ar:
            del ar['rng' if fault == 'rng' else 'simulation/realizations/0/clones/0/results']
    xdr = tmp_path / 'bad.xdr'
    xdr.write_bytes(data)
    before = {path.name:path.read_bytes() for path in tmp_path.iterdir()}
    with pytest.raises((ValueError, KeyError)):
        converter.convert(source, tmp_path / 'output.h5', spinmc_state=xdr)
    assert before == {path.name:path.read_bytes() for path in tmp_path.iterdir()}


def test_spin_state_cli_accepts_released_numeric_text(tmp_path):
    source, xdr = spinmc_pair(tmp_path, 'Potts')
    with h5py.File(source, 'a') as ar:
        del ar['parameters/q']
        ar['parameters/q'] = '3.0'
    output = tmp_path / 'recovered.h5'
    run = subprocess.run([sys.executable, str(SCRIPT), str(source), str(output),
                          '--spinmc-state', str(xdr)], capture_output=True, text=True)
    assert run.returncode == 0, run.stdout + run.stderr
    assert 'native restart conversion remains pending' in run.stdout
    with h5py.File(output, 'r') as ar:
        np.testing.assert_array_equal(ar['migration/spinmc/spins'], [[0.], [2.]])
