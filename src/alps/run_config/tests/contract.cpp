// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/hdf5/archive.hpp>
#include <alps/run_config.hpp>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
const char *schema = R"(
application="test"
schema_version=7
[parameters.count]
type="int64"
required=true
min=1
[parameters.rate]
type="float64"
default=0.5
[parameters.names]
type="string[]"
[parameters.flags]
type="bool[]"
[parameters.z]
type="complex128"
[input.data]
type="path"
[output.results]
type="path"
default="results.h5"
[output.log]
type="path"
[execution]
)";
void require(bool ok) {
    if (!ok)
        throw std::runtime_error("run configuration contract failed");
}
template <class F> void rejects(F f, const std::string &needle) {
    try {
        f();
    } catch (const std::exception &e) {
        require(std::string(e.what()).find(needle) != std::string::npos);
        return;
    }
    throw std::runtime_error("invalid configuration accepted");
}
int main() {
    const std::filesystem::path file = "run-config-contract.toml";
    auto load = [&](const std::string &text, std::string_view rules = schema) {
        {
            std::ofstream out(file);
            out << text;
        }
        return alps::load_run_configuration(file, rules);
    };
    auto run = load("[parameters]\ncount=9007199254740993\nnames=[\"a,b\",\"c\"]\nflags=[]"
                    "\nz={real=1.0,imag=-2.0}\n[input]\ndata=\"data.h5\"\n");
    require(run.parameters["count"].as<std::int64_t>() == 9007199254740993LL);
    require(run.parameters["names"].as<std::vector<std::string>>() ==
            std::vector<std::string>({"a,b", "c"}));
    require(run.parameters["flags"].as<std::vector<bool>>().empty());
    require(run.parameters["z"].as<std::complex<double>>() == std::complex<double>(1, -2));
    require(run.parameters["rate"].as<double>() == 0.5);
    require(run.origins.at("parameters.rate") == "default");
    require(run.origins.at("parameters.count") == "input");
    require(run.input["data"].as<std::string>() ==
            (std::filesystem::current_path() / "data.h5").string());
    rejects([&] { load("[parameters]\ncount=0"); }, "minimum");
    rejects([&] { load("[parameters]\ncount=true"); }, "count");
    rejects([&] { load("[parameters]\ncount=1\nrate=9007199254740993"); }, "exactly");
    rejects([&] { load("[parameters]\ncount=1\nnames=[\"a\",1]"); }, "names");
    rejects([&] { load("[parameters]\ncount=1\ntypo=2"); }, "unknown key");
    rejects([&] { load("[parameters]\ncount=1\n[output]\nresults=\"run-config-contract.toml\""); },
            "replace the TOML run file");
    rejects([&] { load("[parameters]\ncount=1\n[input]\ndata=\"x.h5\"\n[output]\nresults=\"./x.h5\""); },
            "replace input.data");
    rejects([&] { load("[parameters]\ncount=1\n[output]\nlog=\"results.h5\""); },
            "must not replace output.");
    rejects([&] { load("[parameters]\n"); }, "required");
    // Run files cannot override the executable's schema identity or version.
    for (const auto *metadata : {"format_version=1", "schema_version=7", "application=\"test\""})
        rejects([&] { load(metadata); }, "unknown run-file key");
    for (const auto *version : {"0", "-1", "1.5", "true", "\"7\"", "2147483648"})
        rejects([&] {
            load("", std::string("application=\"test\"\nschema_version=") + version);
        }, "schema_version");
    rejects([&] { load("", "application=\"test\""); }, "schema_version");
    rejects([&] { load("", "schema_version=7"); }, "application");
    rejects([&] { load("", "application=\"\"\nschema_version=7"); }, "application");
    require(run.application == "test" && run.schema_version == 7);
    alps::run_configuration supplied_run;
    supplied_run.parameters["count"] = std::int64_t(9007199254740993LL);
    supplied_run.input["data"] = "data.h5";
    auto resolved_run = alps::resolve_run_configuration(supplied_run, schema);
    require(resolved_run.application == "test" && resolved_run.schema_version == 7);
    require(resolved_run.parameters["rate"].as<double>() == 0.5);
    require(resolved_run.input["data"].as<std::string>() == "data.h5");
    auto based_run = alps::resolve_run_configuration(supplied_run, schema,
                                                    std::filesystem::current_path());
    require(based_run.input["data"].as<std::string>() == run.input["data"].as<std::string>());
    require(alps::resolve_run_configuration(run, schema).origins == run.origins);
    resolved_run.application = "another-application";
    rejects([&] { alps::resolve_run_configuration(resolved_run, schema); }, "application schema");
    supplied_run.execution["count"] = 1;
    rejects([&] { alps::resolve_run_configuration(supplied_run, schema); }, "execution.count");
    const std::filesystem::path output = "run-config-contract.h5";
    {
        alps::hdf5::archive ar(output.string(), "w");
        run.save(ar);
    }
    {
        alps::hdf5::archive ar(output.string(), "r");
        std::string application, format;
        int version;
        alps::params parameters;
        ar["application"] >> application;
        ar["schema_version"] >> version;
        ar["format"] >> format;
        ar["parameters"] >> parameters;
        require(application == "test" && version == 7 && format == "alps.run_config.v1");
        require(parameters["count"].as<std::int64_t>() == 9007199254740993LL);
        require(parameters["rate"].as<double>() == 0.5);
        alps::run_configuration restored;
        restored.load(ar);
        require(restored.application == run.application && restored.origins == run.origins);
        require(restored.parameters == run.parameters && restored.input == run.input);
        require(restored.output == run.output && restored.execution == run.execution);
        require(alps::resolve_run_configuration(restored, schema).parameters == run.parameters);
    }
    // A malformed archive must not replace a configuration already held by a caller.
    {
        alps::hdf5::archive ar(output.string(), "a");
        ar["execution/format"] << std::string("unsupported");
    }
    {
        alps::hdf5::archive ar(output.string(), "r");
        auto restored = run;
        rejects([&] { restored.load(ar); }, "format");
        require(restored.parameters == run.parameters && restored.origins == run.origins);
    }
    const auto empty_run = alps::resolve_run_configuration({}, R"(
application="empty-contract"
schema_version=1
[parameters]
[input]
[output]
[execution]
)");
    require(empty_run.origins.empty());
    {
        alps::hdf5::archive ar(output.string(), "a");
        // A fresh empty run and an overwrite must both retain an empty origins
        // group; the overwrite must remove the earlier run's provenance.
        ar["fresh"] << empty_run;
        ar["overwritten"] << run;
        ar["overwritten"] << empty_run;
        for (const auto *path : {"fresh", "overwritten"}) {
            require(ar.is_group(std::string(path) + "/origins"));
            require(ar.list_children(std::string(path) + "/origins").empty());
            auto restored = run;
            ar[path] >> restored;
            require(restored.application == empty_run.application && restored.schema_version == 1);
            require(restored.parameters.empty() && restored.input.empty() && restored.output.empty());
            require(restored.execution.empty() && restored.origins.empty());
        }
    }
    std::filesystem::remove(output);
    alps::params supplied;
    supplied["count"] = 2;
    auto effective = alps::resolve_parameters(supplied, schema);
    require(effective["rate"].as<double>() == 0.5 && !supplied.exists("rate"));
    supplied["parent_only"] = "ignored";
    const auto selected = alps::select_parameters(supplied, schema);
    require(selected.size() == 1 && selected["count"].as<std::int64_t>() == 2);
    require(!selected.exists("rate") && !selected.exists("parent_only"));
    require(alps::select_parameters({}, schema).empty()); // no required checks or defaults
    require(alps::select_parameters(supplied, schema, "execution").empty());
    rejects([&] { alps::resolve_parameters(supplied, schema); }, "unknown key");
    supplied.erase("parent_only");
    supplied["count"] = "2";
    require(alps::select_parameters(supplied, schema)["count"].as<std::string>() == "2");
    rejects([&] { alps::resolve_parameters(supplied, schema); }, "count");
    rejects(
        [&] {
            alps::resolve_parameters({}, R"([parameters.unused]
type="int64"
requred=true
)");
        },
        "unknown schema field");
    rejects(
        [&] {
            alps::resolve_parameters({}, R"([parameters.unused]
type="int64"
min="oops"
)");
        },
        "wrong TOML");
    alps::params exact;
    exact["n"] = 9007199254740992LL;
    rejects(
        [&] {
            alps::resolve_parameters(exact, R"([parameters.n]
type="int64"
min=9007199254740993
)");
        },
        "minimum");
    const char *serialization_schema = R"(
application="serialization-contract"
schema_version=1
[parameters]
"MEASURE[spin \"α\"]\u0001"={type="string"}
message={type="string"}
minimum={type="int64"}
maximum={type="int64"}
exact={type="int64"}
unsigned={type="int64"}
unsigneds={type="int64[]"}
ints={type="int64[]"}
reals={type="float64[]"}
booleans={type="bool[]"}
z={type="complex128"}
zs={type="complex128[]"}
names={type="string[]"}
empty={type="string[]"}
[input]
data={type="string"}
[output]
results={type="string"}
[execution]
seed={type="int64"}
)";
    alps::run_configuration raw;
    const std::string measurement = "MEASURE[spin \"α\"]\x01";
    raw.parameters[measurement] = "Sz";
    raw.parameters["message"] = std::string("Unicode 🎲; \\\n\t\r\"'''\x01\x7f") + '\0';
    raw.parameters["minimum"] = std::numeric_limits<std::int64_t>::min();
    raw.parameters["maximum"] = std::numeric_limits<std::int64_t>::max();
    raw.parameters["exact"] = std::int64_t(9007199254740993LL);
    raw.parameters["unsigned"] = std::uint64_t(9007199254740993ULL);
    raw.parameters["unsigneds"] = std::vector<std::uint64_t>{0, 9007199254740993ULL};
    raw.parameters["ints"] = std::vector<std::int64_t>{-1, 0, 9007199254740993LL};
    raw.parameters["reals"] = std::vector<double>{-0.0, std::nextafter(1., 2.),
        std::numeric_limits<double>::min(), std::numeric_limits<double>::max()};
    raw.parameters["booleans"] = std::vector<bool>{true, false};
    raw.parameters["z"] = std::complex<double>(1., -2.);
    raw.parameters["zs"] = std::vector<std::complex<double>>{{-0.0, 2.}, {3., -4.}};
    raw.parameters["names"] = std::vector<std::string>{"a,b", "c\n\"d"};
    raw.parameters["empty"] = std::vector<std::string>{};
    raw.input["data"] = "data.h5";
    raw.output["results"] = "results.h5";
    raw.execution["seed"] = std::int64_t(17);
    const auto formatted = alps::format_run_configuration(raw);
    const auto roundtrip = load(formatted, serialization_schema);
    raw.parameters["unsigned"] = raw.parameters["unsigned"].as<std::int64_t>();
    raw.parameters["unsigneds"] = raw.parameters["unsigneds"].as<std::vector<std::int64_t>>();
    require(roundtrip.parameters == raw.parameters && roundtrip.input == raw.input);
    require(roundtrip.output == raw.output && roundtrip.execution == raw.execution);
    require(std::signbit(roundtrip.parameters["reals"].as<std::vector<double>>().front()));
    require(std::signbit(roundtrip.parameters["zs"].as<std::vector<std::complex<double>>>().front().real()));
    for (const auto *section : {"parameters", "input", "output", "execution"})
        require(alps::format_run_configuration({}).find(std::string("[") + section + "]") != std::string::npos);
    alps::run_configuration bad;
    bad.parameters["unset"];
    rejects([&] { alps::format_run_configuration(bad); }, "parameters.unset");
    bad.parameters.erase("unset");
    bad.parameters["unsigned"] = std::numeric_limits<std::uint64_t>::max();
    rejects([&] { alps::format_run_configuration(bad); }, "signed 64-bit");
    bad.parameters["unsigned"] = std::vector<std::uint64_t>{0, std::uint64_t(1) << 63};
    rejects([&] { alps::format_run_configuration(bad); }, "signed 64-bit");
    bad.parameters.erase("unsigned");
    bad.parameters["invalid_utf8"] = std::string("\xff");
    rejects([&] { alps::format_run_configuration(bad); }, "roundtrip exactly");
    bad.parameters.erase("invalid_utf8");
    bad.parameters[std::string("\xff")] = "bad key";
    rejects([&] { alps::format_run_configuration(bad); }, "roundtrip exactly");
    bad.parameters.erase(std::string("\xff"));
    bad.parameters["invalid_utf8"] = std::vector<std::string>{"valid", std::string("\xff")};
    rejects([&] { alps::format_run_configuration(bad); }, "roundtrip exactly");
    bad.parameters.erase("invalid_utf8");
    // Some TOML providers use a strtod fallback that rejects subnormal values.
    // Such providers must reject formatting too, rather than emit unreadable runs.
    bool parses_subnormal = true;
    try {
        load("[parameters]\ncount=1\nrate=4.9406564584124654e-324\n");
    } catch (const std::exception &) {
        parses_subnormal = false;
    }
    alps::run_configuration subnormal;
    subnormal.parameters["count"] = 1;
    subnormal = alps::resolve_run_configuration(subnormal, schema);
    subnormal.parameters["rate"] = std::numeric_limits<double>::denorm_min();
    if (parses_subnormal)
        require(load(alps::format_run_configuration(subnormal)).parameters["rate"].as<double>() ==
                std::numeric_limits<double>::denorm_min());
    else
        rejects([&] { alps::format_run_configuration(subnormal); }, "TOML provider");
    for (const auto nonfinite : {std::numeric_limits<double>::infinity(),
                                 std::numeric_limits<double>::quiet_NaN()}) {
        bad.input["bad"] = nonfinite;
        rejects([&] { alps::format_run_configuration(bad); }, "finite");
        bad.input["bad"] = std::vector<double>{1., nonfinite};
        rejects([&] { alps::format_run_configuration(bad); }, "finite");
        bad.input["bad"] = std::complex<double>(1., nonfinite);
        rejects([&] { alps::format_run_configuration(bad); }, "finite");
        bad.input["bad"] = std::vector<std::complex<double>>{{nonfinite, 0.}};
        rejects([&] { alps::format_run_configuration(bad); }, "finite");
    }
    std::filesystem::remove(file);
}
