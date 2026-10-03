// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/hdf5/archive.hpp>
#include <alps/run_config.hpp>
#include <filesystem>
#include <fstream>
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
    std::filesystem::remove(output);
    alps::params supplied;
    supplied["count"] = 2;
    auto effective = alps::resolve_parameters(supplied, schema);
    require(effective["rate"].as<double>() == 0.5 && !supplied.exists("rate"));
    supplied["count"] = "2";
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
    std::filesystem::remove(file);
}
