// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/map.hpp>
#include <alps/run_config.hpp>
#include <map>
#include <set>
#include <sstream>
#include <toml++/toml.hpp>
namespace alps {
namespace {
using value = params::value_type;
[[noreturn]] void fail(const std::string &key, const std::string &reason) {
    throw std::invalid_argument("Configuration '" + key + "': " + reason);
}
template <class T> auto toml_value(const T &x, const std::string &key) {
    if constexpr (params_ns::detail::vector_type<T>::value) {
        toml::array out;
        for (const auto &element : x)
            out.push_back(toml_value(static_cast<typename T::value_type>(element), key));
        return out;
    } else if constexpr (std::is_same_v<T, std::uint64_t>) {
        if (x > std::uint64_t(std::numeric_limits<std::int64_t>::max()))
            fail(key, "integer exceeds TOML's signed 64-bit range");
        return static_cast<std::int64_t>(x);
    } else if constexpr (std::is_same_v<T, std::complex<double>>) {
        if (!std::isfinite(x.real()) || !std::isfinite(x.imag()))
            fail(key, "value must be finite");
        toml::table out{{"real", x.real()}, {"imag", x.imag()}};
        out.is_inline(true);
        return out;
    } else {
        if constexpr (std::is_same_v<T, double>)
            if (!std::isfinite(x))
                fail(key, "value must be finite");
        return x;
    }
}
toml::table format_section(const params &values, const std::string &section) {
    toml::table out;
    for (const auto &[key, item] : values)
        item.apply_visitor([&](const auto &x) {
            using T = std::decay_t<decltype(x)>;
            if constexpr (std::is_same_v<T, params_ns::detail::None>)
                fail(section + "." + key, "cannot serialize an unset value");
            else
                out.insert(key, toml_value(x, section + "." + key));
        });
    return out;
}
std::string text(const toml::node_view<const toml::node> &node, const std::string &key) {
    if (auto s = node.value<std::string>())
        return *s;
    fail(key, "expected a string");
}
std::string location(const toml::node &node) {
    std::ostringstream out;
    out << node.source();
    return out.str();
}
template <class T> T scalar(const toml::node &node, const std::string &key) {
    if constexpr (std::is_same_v<T, double>) {
        if (auto n = node.as_integer())
            return params_ns::detail::convert<double>(n->get(), key);
        if (auto n = node.as_floating_point()) {
            if (!std::isfinite(n->get()))
                fail(key, "expected a finite real");
            return n->get();
        }
    } else if constexpr (std::is_same_v<T, std::complex<double>>) {
        const auto *t = node.as_table();
        if (t && t->size() == 2 && t->contains("real") && t->contains("imag"))
            return {scalar<double>(*t->get("real"), key), scalar<double>(*t->get("imag"), key)};
    } else if constexpr (std::is_same_v<T, std::uint64_t>) {
        if (auto n = node.as_integer())
            return params_ns::detail::convert<T>(n->get(), key);
    } else if (auto n = node.as<T>())
        return n->get();
    fail(key, "wrong TOML value type at " + location(node));
}
template <class T>
void assign(value &out, const toml::node &node, const std::string &key, bool array) {
    if (!array) {
        out = scalar<T>(node, key);
        return;
    }
    const auto *a = node.as_array();
    if (!a)
        fail(key, "expected an array at " + location(node));
    std::vector<T> values;
    values.reserve(a->size());
    for (const auto &x : *a)
        values.push_back(scalar<T>(x, key));
    out = values;
}
value read_value(const toml::node &node, const std::string &type, const std::string &key) {
    value out(key);
    const bool array = type.size() > 2 && type.substr(type.size() - 2) == "[]";
    const auto base = array ? type.substr(0, type.size() - 2) : type;
    if (base == "int64")
        assign<std::int64_t>(out, node, key, array);
    else if (base == "float64")
        assign<double>(out, node, key, array);
    else if (base == "bool")
        assign<bool>(out, node, key, array);
    else if (base == "string" || base == "path")
        assign<std::string>(out, node, key, array);
    else if (base == "complex128")
        assign<std::complex<double>>(out, node, key, array);
    else
        fail(key, "unsupported schema type '" + type + "'");
    return out;
}
template <class T> void coerce(value &target, const value &source, bool array) {
    if (array)
        target = source.as<std::vector<T>>();
    else
        target = source.as<T>();
}
value convert(const value &source, const std::string &type, const std::string &key) {
    value out(key);
    const bool array = type.size() > 2 && type.substr(type.size() - 2) == "[]";
    const auto base = array ? type.substr(0, type.size() - 2) : type;
    if (base == "int64")
        coerce<std::int64_t>(out, source, array);
    else if (base == "float64")
        coerce<double>(out, source, array);
    else if (base == "bool")
        coerce<bool>(out, source, array);
    else if (base == "string" || base == "path")
        coerce<std::string>(out, source, array);
    else if (base == "complex128")
        coerce<std::complex<double>>(out, source, array);
    else
        fail(key, "unsupported schema type '" + type + "'");
    return out;
}
void validate_rule(const toml::table &rule, const std::string &key) {
    const std::set<std::string_view> fields{"type", "required", "default",    "min",
                                            "max",  "choices",  "description"};
    for (const auto &[name, node] : rule)
        if (!fields.count(name.str()))
            fail(key, "unknown schema field '" + std::string(name.str()) + "'");
    const auto type = text(rule["type"], key);
    const std::set<std::string> types{"int64",  "float64",    "bool",        "string",
                                      "path",   "complex128", "int64[]",     "float64[]",
                                      "bool[]", "string[]",   "complex128[]", "path[]"};
    if (!types.count(type))
        fail(key, "unsupported schema type '" + type + "'");
    if (rule.contains("required") && !rule["required"].is_boolean())
        fail(key, "required must be boolean");
    if (rule.contains("description") && !rule["description"].is_string())
        fail(key, "description must be a string");
    const auto base = type.substr(0, type.find('['));
    for (const auto *bound : {"min", "max"})
        if (const auto *node = rule.get(bound)) {
            if (base != "int64" && base != "float64")
                fail(key, "bounds require a numeric type");
            read_value(*node, base, key);
        }
    if (rule.contains("min") && rule.contains("max")) {
        const bool ordered = base == "int64" ? scalar<std::int64_t>(*rule.get("min"), key) <=
                                                   scalar<std::int64_t>(*rule.get("max"), key)
                                             : scalar<double>(*rule.get("min"), key) <=
                                                   scalar<double>(*rule.get("max"), key);
        if (!ordered)
            fail(key, "schema minimum exceeds maximum");
    }
    if (rule.contains("choices")) {
        const auto *choices = rule["choices"].as_array();
        if (!choices || choices->empty())
            fail(key, "choices must be a nonempty array");
        for (const auto &choice : *choices)
            read_value(choice, type, key);
    }
}
template <class T> void check_element(const toml::table &rule, const T &x, const std::string &key) {
    if constexpr (params_ns::detail::vector_type<T>::value) {
        for (const auto &element : x)
            check_element(rule, static_cast<typename T::value_type>(element), key);
    } else if constexpr (std::is_arithmetic_v<T> && !std::is_same_v<T, bool>) {
        if constexpr (std::is_floating_point_v<T>)
            if (!std::isfinite(x))
                fail(key, "value must be finite");
        if (const auto *min = rule.get("min"); min && x < scalar<T>(*min, key))
            fail(key, "value is below minimum");
        if (const auto *max = rule.get("max"); max && x > scalar<T>(*max, key))
            fail(key, "value is above maximum");
    } else if constexpr (std::is_same_v<T, std::complex<double>>) {
        if (!std::isfinite(x.real()) || !std::isfinite(x.imag()))
            fail(key, "value must be finite");
    }
}
void check_rule(const toml::table &rule, const value &v, const std::string &key) {
    v.apply_visitor([&](const auto &x) { check_element(rule, x, key); });
    if (const auto *choices = rule["choices"].as_array()) {
        bool found = false;
        const auto type = text(rule["type"], key);
        for (const auto &choice : *choices)
            found |= read_value(choice, type, key).equals(v);
        if (!found)
            fail(key, "value is not one of the allowed choices");
    }
}
const toml::table &definitions(const toml::table &schema, const std::string &section) {
    const auto *table = schema[section].as_table();
    if (!table)
        fail(section, "application schema has no definition table");
    return *table;
}
params resolve(const params &supplied, const toml::table &rules, const std::filesystem::path &base,
               const std::string &section, std::map<std::string, std::string> *origins = nullptr) {
    for (const auto &entry : supplied)
        if (!rules.contains(entry.first))
            fail(section + "." + entry.first, "unknown key");
    params out;
    for (const auto &[name, node] : rules) {
        const std::string key(name.str()), qualified = section + "." + key;
        const auto *rule = node.as_table();
        if (!rule)
            fail(qualified, "schema definition must be a table");
        validate_rule(*rule, qualified);
        const auto type = text((*rule)["type"], qualified);
        if (const auto *d = rule->get("default"))
            check_rule(*rule, read_value(*d, type, qualified), qualified);
        if (supplied.exists(key))
            out[key] = convert(supplied[key], type, qualified);
        else if (const auto *d = rule->get("default"))
            out[key] = read_value(*d, type, qualified);
        else if ((*rule)["required"].value_or(false))
            fail(qualified, "required value is missing");
        else
            continue;
        check_rule(*rule, out[key], qualified);
        if ((type == "path" || type == "path[]") && !base.empty()) {
            auto paths = run_paths(out[key]);
            for (auto& filename : paths) {
                auto path = std::filesystem::path(filename);
                filename = (path.is_relative() ? base / path : path).lexically_normal().string();
            }
            if (type == "path") out[key] = paths.front();
            else out[key] = paths;
        }
        if (origins)
            (*origins)[qualified] = supplied.exists(key) ? "input" : "default";
    }
    return out;
}
run_configuration identity(const toml::table &schema) {
    run_configuration run;
    run.application = text(schema["application"], "schema.application");
    if (run.application.empty())
        fail("schema.application", "expected a nonempty application name");
    const auto version = schema["schema_version"].value<std::int64_t>();
    if (!schema["schema_version"].is_integer() || !version || *version < 1)
        fail("schema.schema_version", "expected a positive integer");
    run.schema_version = params_ns::detail::convert<int>(*version, "schema.schema_version");
    return run;
}
// Anchored path outputs must not replace the run file, an input or another
// output. Relative paths have no base yet and are checked once resolved.
void check_outputs(const run_configuration &run, const toml::table &schema) {
    std::map<std::filesystem::path, std::string> claimed;
    if (!run.source_file.empty())
        claimed.emplace(run.source_file, "the TOML run file");
    for (const auto &[section, values] :
         {std::pair<std::string, const params *>{"input", &run.input}, {"output", &run.output}})
        for (const auto &[name, node] : definitions(schema, section)) {
            const std::string key(name.str()), qualified = section + "." + key;
            const auto *rule = node.as_table();
            if (!rule || !values->exists(key))
                continue;
            const auto type = (*rule)["type"].value_or(std::string());
            if (type != "path" && type != "path[]") continue;
            for (std::filesystem::path path : run_paths((*values)[key])) {
                if (path.is_relative()) continue;
                const auto [claim, added] = claimed.emplace(std::filesystem::weakly_canonical(path), qualified);
                if (!added && section == "output") fail(qualified, "output must not replace " + claim->second);
            }
        }
}
run_configuration resolve_run(const run_configuration &supplied, const toml::table &schema,
                              const std::filesystem::path &base) {
    auto run = identity(schema);
    if (!supplied.application.empty() &&
        (supplied.application != run.application || supplied.schema_version != run.schema_version))
        fail("application", "configuration does not match the application schema");
    run.source_file = supplied.source_file;
    const auto directory = base.empty() ? base : std::filesystem::absolute(base);
    for (const auto &[section, member] :
         std::vector<std::pair<std::string, params run_configuration::*>>{
             {"parameters", &run_configuration::parameters}, {"input", &run_configuration::input},
             {"output", &run_configuration::output}, {"execution", &run_configuration::execution}}) {
        run.*member = resolve(supplied.*member, definitions(schema, section), directory, section,
                              &run.origins);
    }
    // Revalidation by an application must retain the origin of defaults
    // already inserted by the TOML loader.
    if (!supplied.application.empty())
        for (auto &[key, origin] : run.origins)
            if (const auto it = supplied.origins.find(key); it != supplied.origins.end())
                origin = it->second;
    check_outputs(run, schema);
    return run;
}
} // namespace
run_configuration resolve_run_configuration(const run_configuration &supplied,
                                           std::string_view schema_text,
                                           const std::filesystem::path &base_directory) {
    return resolve_run(supplied, toml::parse(schema_text), base_directory);
}
std::string extend_run_schema(const std::filesystem::path &file, std::string_view base,
                             const std::function<char const*(std::string const&)> &type) {
    auto definitions = toml::parse(base);
    auto* parameters = definitions["parameters"].as_table();
    if (!parameters) fail("parameters", "schema must contain a parameter table");
    if (!file.empty()) {
        const auto document = toml::parse_file(file.string());
        if (const auto* supplied = document["parameters"].as_table())
            for (auto const& [name, node] : *supplied) {
                const auto key = std::string(name.str());
                if (parameters->contains(key)) continue;
                const char* value_type = type ? type(key) : nullptr;
                if (!value_type) value_type = node.is_integer() ? "int64" : node.is_floating_point() ?
                    "float64" : node.is_string() ? "string" : node.is_boolean() ? "bool" : nullptr;
                if (!value_type) fail(key, "model parameter must be a scalar or have an explicit schema type");
                parameters->insert(key, toml::table{{"type", value_type}});
            }
    }
    std::ostringstream output;
    output << definitions;
    return output.str();
}

params resolve_parameters(const params &supplied, std::string_view schema_text,
                          const std::string &section) {
    const auto schema = toml::parse(schema_text);
    return resolve(supplied, definitions(schema, section), {}, section);
}
params select_parameters(const params &supplied, std::string_view schema_text,
                         const std::string &section) {
    const auto schema = toml::parse(schema_text);
    params out;
    for (const auto &[name, node] : definitions(schema, section)) {
        const std::string key(name.str()), qualified = section + "." + key;
        const auto *rule = node.as_table();
        if (!rule)
            fail(qualified, "schema definition must be a table");
        validate_rule(*rule, qualified);
        if (supplied.exists(key))
            out[key] = supplied[key];
    }
    return out;
}
std::string format_run_configuration(const run_configuration &run) {
    toml::table document;
    document.insert("parameters", format_section(run.parameters, "parameters"));
    document.insert("input", format_section(run.input, "input"));
    document.insert("output", format_section(run.output, "output"));
    document.insert("execution", format_section(run.execution, "execution"));
    std::ostringstream out;
    out << toml::toml_formatter(document) << '\n';
    auto encoded = out.str();
    try {
        if (toml::parse(encoded) != document)
            fail("run", "values do not roundtrip exactly through TOML, including UTF-8 text");
    } catch (const toml::parse_error &error) {
        fail("run", "TOML provider cannot read this representation: " + std::string(error.description()));
    }
    return encoded;
}
run_configuration load_run_configuration(const std::filesystem::path &filename,
                                         std::string_view schema_text) {
    const auto schema = toml::parse(schema_text);
    const auto document = toml::parse_file(filename.string());
    const std::set<std::string_view> top{"parameters", "input", "output", "execution"};
    for (const auto &[key, node] : document)
        if (!top.count(key.str()))
            fail(std::string(key.str()), "unknown run-file key at " + location(node));
    auto run = identity(schema);
    const auto base = std::filesystem::absolute(filename).parent_path();
    for (const auto &item :
         std::vector<std::pair<std::string, params *>>{{"parameters", &run.parameters},
                                                       {"input", &run.input},
                                                       {"output", &run.output},
                                                       {"execution", &run.execution}}) {
        const auto &section = item.first;
        const auto &rules = definitions(schema, section);
        params supplied;
        if (const auto *n = document.get(section)) {
            const auto *table = n->as_table();
            if (!table)
                fail(section, "expected a table");
            for (const auto &[key, node] : *table) {
                const auto *rule = rules[key.str()].as_table();
                const std::string qualified = section + "." + std::string(key.str());
                if (!rule)
                    fail(qualified, "unknown key at " + location(node));
                supplied[std::string(key.str())] =
                    read_value(node, text((*rule)["type"], qualified), qualified);
            }
        }
        *item.second = resolve(supplied, rules, base, section, &run.origins);
    }
    run.source_file = std::filesystem::weakly_canonical(filename).string();
    check_outputs(run, schema);
    return run;
}
void run_configuration::save(hdf5::archive &ar) const {
    ar["format"] << std::string("alps.run_config.v1");
    ar["application"] << application;
    ar["schema_version"] << schema_version;
    ar["parameters"] << parameters;
    ar["input"] << input;
    ar["output"] << output;
    ar["execution"] << execution;
    ar["origins"] << origins;
    // Empty maps do not create an archive group; retain an explicit empty
    // origins group so every saved run can be loaded by the same format.
    if (origins.empty())
        ar.create_group("origins");
}
void run_configuration::load(hdf5::archive &ar) {
    run_configuration restored;
    std::string format;
    ar["format"] >> format;
    if (format != "alps.run_config.v1")
        fail("format", "unsupported archived run configuration");
    ar["application"] >> restored.application;
    ar["schema_version"] >> restored.schema_version;
    if (restored.application.empty() || restored.schema_version < 1)
        fail("application", "invalid archived configuration identity");
    ar["parameters"] >> restored.parameters;
    ar["input"] >> restored.input;
    ar["output"] >> restored.output;
    ar["execution"] >> restored.execution;
    ar["origins"] >> restored.origins;
    *this = std::move(restored);
}
} // namespace alps
