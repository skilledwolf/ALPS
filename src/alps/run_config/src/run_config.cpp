// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/map.hpp>
#include <alps/run_config.hpp>
#include <set>
#include <sstream>
#include <toml++/toml.hpp>
namespace alps {
namespace {
using value = params::value_type;
[[noreturn]] void fail(const std::string &key, const std::string &reason) {
    throw std::invalid_argument("Configuration '" + key + "': " + reason);
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
                                      "bool[]", "string[]",   "complex128[]"};
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
        if (type == "path" && !base.empty()) {
            auto path = std::filesystem::path(out[key].as<std::string>());
            if (path.is_relative())
                path = base / path;
            out[key] = path.lexically_normal().string();
        }
        if (origins)
            (*origins)[qualified] = supplied.exists(key) ? "input" : "default";
    }
    return out;
}
} // namespace
params resolve_parameters(const params &supplied, std::string_view schema_text,
                          const std::string &section) {
    const auto schema = toml::parse(schema_text);
    return resolve(supplied, definitions(schema, section), {}, section);
}
run_configuration load_run_configuration(const std::filesystem::path &filename,
                                         std::string_view schema_text) {
    const auto schema = toml::parse(schema_text);
    const auto document = toml::parse_file(filename.string());
    const std::set<std::string_view> top{"parameters", "input", "output", "execution"};
    for (const auto &[key, node] : document)
        if (!top.count(key.str()))
            fail(std::string(key.str()), "unknown run-file key at " + location(node));
    run_configuration run;
    run.application = text(schema["application"], "schema.application");
    if (run.application.empty())
        fail("schema.application", "expected a nonempty application name");
    const auto version = schema["schema_version"].value<std::int64_t>();
    if (!schema["schema_version"].is_integer() || !version || *version < 1)
        fail("schema.schema_version", "expected a positive integer");
    run.schema_version = params_ns::detail::convert<int>(*version, "schema.schema_version");
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
}
} // namespace alps
