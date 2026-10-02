// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/run_config.hpp>
#include <filesystem>
#include <fstream>
#include <stdexcept>
const char *schema = R"(
application="test"
schema_version=1
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
    auto load = [&](const std::string &text) {
        {
            std::ofstream out(file);
            out << text;
        }
        return alps::load_run_configuration(file, schema);
    };
    const std::string header = "format_version=1\napplication=\"test\"\nschema_version=1\n";
    auto run = load(header + "[parameters]\ncount=9007199254740993\nnames=[\"a,b\",\"c\"]\nflags=[]"
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
    rejects([&] { load(header + "[parameters]\ncount=0"); }, "minimum");
    rejects([&] { load(header + "[parameters]\ncount=true"); }, "count");
    rejects([&] { load(header + "[parameters]\ncount=1\nrate=9007199254740993"); }, "exactly");
    rejects([&] { load(header + "[parameters]\ncount=1\nnames=[\"a\",1]"); }, "names");
    rejects([&] { load(header + "[parameters]\ncount=1\ntypo=2"); }, "unknown key");
    rejects([&] { load(header + "[parameters]\n"); }, "required");
    rejects([&] { load("format_version=0\n"); }, "format_version");
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
