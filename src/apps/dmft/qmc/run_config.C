// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include "run_config.h"
#include "dmft_schema.hpp"
#include "U_matrix.h"
#include "hirschfyeaux.h"
#include "interaction_expansion_choice.h"
#include "hybridization/input.hpp"
#include <alps/cthyb.hpp>
#include <alps/ctint.hpp>
#include <toml++/toml.hpp>
#include <charconv>
#include <fstream>
#include <iterator>
#include <limits>
#include <set>
#include <sstream>
namespace alps::dmft {
namespace {
// The driver names the grid and flavor counts N, NMATSUBARA and FLAVORS; CT-HYB
// uses N_TAU, N_MATSUBARA and N_ORBITALS for the same quantities.
constexpr std::pair<std::string_view, std::string_view> renamed[] = {
    {"N_TAU", "N"}, {"N_MATSUBARA", "NMATSUBARA"}, {"N_ORBITALS", "FLAVORS"}};
std::string driver_name(std::string_view solver_key) {
  for (const auto& [solver, driver] : renamed)
    if (solver == solver_key) return std::string(driver);
  return std::string(solver_key);
}
// Scientific inputs the driver passes to a solver that declares them.
constexpr const char* forwarded_inputs[] = {"retarded_interaction", "retarded_interaction_format",
                                            "retarded_interaction_coordinate"};
std::string read_text(const std::filesystem::path& file) {
  std::ifstream stream(file);
  if (!stream) throw std::invalid_argument("Cannot open DMFT input.solver_schema: " + file.string());
  return {std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>()};
}
std::string compose(const run_configuration& selectors, const std::filesystem::path& base) {
  auto root = toml::parse(base_schema);
  const auto child = toml::parse(solver_schema(selectors, base));
  auto& parameters = *root["parameters"].as_table();
  // The driver's rule wins for shared keys; validation also applies the solver's rule.
  if (const auto* source = child["parameters"].as_table())
    for (const auto& [name, rule] : *source)
      if (!parameters.contains(driver_name(name.str()))) parameters.insert(driver_name(name.str()), rule);
  if (external(selected_solver(selectors)))
    if (const auto* source = child["execution"].as_table())
      for (const auto& [name, rule] : *source)
        if (!root["execution"].as_table()->contains(name)) root["execution"].as_table()->insert(name, rule);
  for (const auto* key : forwarded_inputs)
    if (const auto* rule = child["input"][key].node()) root["input"].as_table()->insert(key, *rule);
  const auto flavors = selectors.parameters.value_or<std::int64_t>(
      "FLAVORS", *root["parameters"]["FLAVORS"]["default"].value<std::int64_t>());
  if (flavors < 1 || flavors > 128) throw std::invalid_argument("DMFT FLAVORS must be in [1,128]");
  for (int f = 0; f < flavors; ++f)
    for (const auto* prefix : {"EPS_", "EPSSQ_"})
      parameters.insert_or_assign(std::string(prefix) + std::to_string(f), toml::table{{"type", "float64"}});
  for (int band = 0; band < flavors / 2; ++band)
    parameters.insert_or_assign("t" + std::to_string(band), toml::table{{"type", "float64"}});
  std::ostringstream out; out << root; return out.str();
}
bool input_path(const toml::table& schema, const std::string& name) {
  return schema["input"][name]["type"].value_or(std::string{}) == "path";
}
bool text_name(const std::string& name, const run_configuration& run) {
  const bool omega = run.execution["loop"].as<std::string>() == "omega";
  if (name == "G_tau" || (omega && name == "G_omega")) return true;
  for (const auto* prefix : {"G0_omega", "G0_omegareal", "G0_tau", "G_omega", "G_omegareal", "G_tau", "selfenergy"}) {
    if (!omega && prefix != std::string_view("G_tau")) continue;
    const auto beginning = std::string(prefix) + "_";
    if (name.compare(0, beginning.size(), beginning) != 0) continue;
    unsigned iteration = 0;
    const auto parsed = std::from_chars(name.data() + beginning.size(), name.data() + name.size(), iteration);
    if (parsed.ec == std::errc{} && parsed.ptr == name.data() + name.size() &&
        iteration > 0 && iteration <= run.execution["max_iterations"].as<unsigned>()) return true;
  }
  return false;
}
void validate_text_outputs(const run_configuration& run, const toml::table& schema,
                           std::set<std::filesystem::path> protected_paths) {
  if (!run.output["text"].as<bool>()) return;
  const auto directory = std::filesystem::weakly_canonical(run.output["text_directory"].as<std::string>());
  if (!std::filesystem::is_directory(directory)) throw std::invalid_argument("DMFT output.text_directory must exist");
  if (!run.source_file.empty()) protected_paths.insert(run.source_file);
  for (const auto& [name, value] : run.input)
    if (input_path(schema, name)) protected_paths.insert(std::filesystem::weakly_canonical(value.as<std::string>()));
  for (const auto& path : protected_paths)
    if (path.parent_path() == directory && text_name(path.filename().string(), run))
      throw std::invalid_argument("DMFT text output would overwrite a run, input or explicit output file: " + path.string());
  // Check existing aliases without enumerating every possible iteration.
  for (const auto& entry : std::filesystem::directory_iterator(directory))
    if (text_name(entry.path().filename().string(), run) &&
        (entry.is_directory() || protected_paths.count(std::filesystem::weakly_canonical(entry.path()))))
      throw std::invalid_argument("DMFT text output collides with a protected file or directory: " + entry.path().string());
}
void validate_hirschfye_execution(const params& execution) {
  if (execution["bins"].as<std::uint64_t>() % 2)
    throw std::invalid_argument("Hirsch-Fye execution.bins must be even and at least two");
}
void validate_solver(const run_configuration& run) {
  const auto& p = run.parameters;
  const auto kind = selected_solver(run);
  if (kind == solver_kind::custom) {
    if (!run.execution.exists("solver_input"))
      throw std::invalid_argument("A custom DMFT solver requires execution.solver_input");
  } else if (run.input.exists("solver_schema") || run.execution.exists("solver_input")) {
    throw std::invalid_argument("input.solver_schema and execution.solver_input apply only to custom solvers");
  }
  if (run.execution["loop"].as<std::string>() == "tau" && !receives_delta(run))
    throw std::invalid_argument("The DMFT tau loop requires a solver that receives a hybridization function");
  const auto schema = solver_schema(run);
  switch (kind) {
    case solver_kind::hybridization: {
      const auto parameters = cthyb::prepare_parameters(solver_parameters(p, schema));
      if (run.input.exists("retarded_interaction")) cthyb_input::retarded_kernel(parameters, run.input);
      break;
    }
    case solver_kind::interaction:
      ctint::prepare_parameters(solver_parameters(p, schema));
      break;
    case solver_kind::hirsch_fye:
      prepare_hirschfye_parameters(solver_parameters(p, schema));
      validate_hirschfye_execution(run.execution);
      break;
    case solver_kind::custom:
      resolve_parameters(solver_parameters(p, schema), schema);
      break;
    case solver_kind::interaction_expansion:
      // The scheduler CT-INT supports multiband densities; the standalone schema does not.
      if (select_interaction_expansion(p["FLAVORS"].as<int>(), p["SITES"].as<int>()) ==
          interaction_expansion_choice::unsupported)
        throw std::invalid_argument("Unsupported Interaction Expansion dimensions");
      if (p["THERMALIZATION"].as<std::uint64_t>() > std::numeric_limits<unsigned int>::max())
        throw std::invalid_argument("Interaction Expansion THERMALIZATION exceeds the scheduler counter");
      if (p.exists("NMATSUBARA_MEASUREMENTS") &&
          p["NMATSUBARA_MEASUREMENTS"].as<int>() > p["NMATSUBARA"].as<int>())
        throw std::invalid_argument("NMATSUBARA_MEASUREMENTS must not exceed NMATSUBARA");
      break;
  }
}
void validate(const run_configuration& run, const toml::table& schema) {
  const auto& p = run.parameters;
  if (p["BETA"].as<double>() <= 0 || p["RELAX_RATE"].as<double>() <= 0)
    throw std::invalid_argument("DMFT BETA and RELAX_RATE must be positive");
  if (p["SITES"].as<int>() != 1)
    throw std::invalid_argument("Use the cluster framework for DMFT SITES != 1");
  if ((p["ANTIFERROMAGNET"].as<bool>() || p["SYMMETRIZATION"].as<bool>()) && p["FLAVORS"].as<int>() % 2)
    throw std::invalid_argument("DMFT symmetry operations require paired flavors");
  if (p["ANTIFERROMAGNET"].as<bool>() && p["SYMMETRIZATION"].as<bool>())
    throw std::invalid_argument("ANTIFERROMAGNET and SYMMETRIZATION are incompatible");
  if (run.input.exists("dos") && (p["SEMICIRCLE_HILBERT"].as<bool>() || run.execution["loop"].as<std::string>() == "tau"))
    throw std::invalid_argument("input.dos requires the general omega Hilbert transform");
  if (p["N"].as<int>() == std::numeric_limits<int>::max())
    throw std::invalid_argument("DMFT N + 1 exceeds the supported index range");
  const auto flavors = p["FLAVORS"].as<unsigned>();
  if ((p["N"].as<std::uint64_t>() + 1) * flavors > std::numeric_limits<unsigned>::max() ||
      p["NMATSUBARA"].as<std::uint64_t>() * flavors > std::numeric_limits<unsigned>::max())
    throw std::invalid_argument("DMFT Green-function dimensions exceed the supported storage range");
  for (const auto& [name, value] : run.input)
    if (input_path(schema, name) && !std::filesystem::is_regular_file(value.as<std::string>()))
      throw std::invalid_argument("DMFT input." + name + " must be an existing file");
  std::set<std::filesystem::path> outputs;
  for (const auto* key : {"results", "final_tau", "final_omega"})
    if (run.output.exists(key)) {
      const auto path = std::filesystem::path(run.output[key].as<std::string>());
      outputs.insert(std::filesystem::weakly_canonical(path));
      if (!std::filesystem::is_directory(path.parent_path().empty() ? "." : path.parent_path()))
        throw std::invalid_argument(std::string("DMFT output.") + key + " parent directory does not exist");
      if (std::filesystem::is_directory(path))
        throw std::invalid_argument(std::string("DMFT output.") + key + " is a directory");
    }
  validate_text_outputs(run, schema, outputs);
  static_cast<void>(U_matrix(p, run.input));  // A malformed matrix fails before any output.
  validate_solver(run);
}
}
params prepare_hirschfye_parameters(const params& supplied) {
  auto parameters = resolve_parameters(supplied, hirschfye_schema);
  if (parameters["BETA"].as<double>() <= 0.)
    throw std::invalid_argument("Hirsch-Fye BETA must be positive");
  const auto slices = parameters["N"].as<std::size_t>();
  if (slices > std::numeric_limits<std::size_t>::max() / sizeof(double) / slices)
    throw std::invalid_argument("Hirsch-Fye matrix dimensions exceed the storage range");
  const auto sweeps = parameters["SWEEPS"].as<std::uint64_t>();
  const auto thermalization = parameters["THERMALIZATION"].as<std::uint64_t>();
  if (sweeps > std::numeric_limits<std::uint64_t>::max() - thermalization)
    throw std::invalid_argument("Hirsch-Fye total sweep count exceeds the counter range");
  hirschfye_lambda(parameters["BETA"].as<double>(), parameters["U"].as<double>(), slices);
  return parameters;
}
matsubara_green_function_t prepare_hirschfye_run(run_configuration& run) {
  run = resolve_run_configuration(run, hirschfye_schema, std::filesystem::current_path());
  run.parameters = prepare_hirschfye_parameters(run.parameters);
  validate_hirschfye_execution(run.execution);
  const auto output = std::filesystem::path(run.output["results"].as<std::string>());
  if (std::filesystem::is_directory(output) || !std::filesystem::is_directory(output.parent_path()))
    throw std::invalid_argument("Hirsch-Fye output.results must name a file in an existing directory");
  const auto input = run.input["g0"].as<std::string>();
  if (!std::filesystem::is_regular_file(input))
    throw std::invalid_argument("Hirsch-Fye input.g0 must be an existing file");
  matsubara_green_function_t green(run.parameters["NMATSUBARA"].as<unsigned>(), 1, 2);
  alps::hdf5::archive archive(input, "r");
  read_flavor_vectors(archive, "/G0", green);
  return green;
}
solver_kind selected_solver(const run_configuration& run) {
  static const auto fallback = *toml::parse(base_schema)["execution"]["solver"]["default"].value<std::string>();
  const auto name = run.execution.value_or("solver", fallback);
  if (name == "hybridization") return solver_kind::hybridization;
  if (name == "interaction") return solver_kind::interaction;
  if (name == "hirschfye") return solver_kind::hirsch_fye;
  if (name == "Interaction Expansion") return solver_kind::interaction_expansion;
  return solver_kind::custom;
}
bool external(solver_kind kind) {
  return kind == solver_kind::hybridization || kind == solver_kind::interaction ||
         kind == solver_kind::hirsch_fye || kind == solver_kind::custom;
}
bool receives_delta(const run_configuration& run) {
  const auto kind = selected_solver(run);
  return kind == solver_kind::hybridization ||
         (kind == solver_kind::custom && run.execution.value_or("solver_input", "") == "delta");
}
std::string solver_schema(const run_configuration& run, const std::filesystem::path& base) {
  switch (selected_solver(run)) {
    case solver_kind::hybridization: return std::string(cthyb::schema());
    case solver_kind::interaction:
    case solver_kind::interaction_expansion: return std::string(ctint::schema());
    case solver_kind::hirsch_fye: return std::string(hirschfye_schema);
    case solver_kind::custom: break;
  }
  if (!run.input.exists("solver_schema"))
    throw std::invalid_argument("A custom DMFT solver requires input.solver_schema");
  auto path = std::filesystem::path(run.input["solver_schema"].as<std::string>());
  if (path.is_relative() && !base.empty()) path = base / path;
  return read_text(path);
}
params solver_parameters(const params& dmft, std::string_view schema) {
  auto mapped = dmft;
  for (const auto& [solver, driver] : renamed)
    if (dmft.exists(std::string(driver))) mapped[std::string(solver)] = dmft[std::string(driver)];
  return select_parameters(mapped, schema);
}
params solver_inputs(const params& dmft, std::string_view schema) {
  params forwarded;
  for (const auto* key : forwarded_inputs)
    if (dmft.exists(key)) forwarded[key] = dmft[key];
  return select_parameters(forwarded, schema, "input");
}
std::string schema_for_run(const std::filesystem::path& file) {
  if (file.empty()) return compose({}, {});
  const auto raw = toml::parse_file(file.string());
  run_configuration selectors;
  if (const auto* node = raw["execution"]["solver"].node()) {
    if (auto value = node->value<std::string>()) selectors.execution["solver"] = *value;
    else throw std::invalid_argument("DMFT execution.solver must be a string");
  }
  if (const auto* node = raw["input"]["solver_schema"].node()) {
    if (auto value = node->value<std::string>()) selectors.input["solver_schema"] = *value;
    else throw std::invalid_argument("DMFT input.solver_schema must be a path");
  }
  if (const auto* node = raw["parameters"]["FLAVORS"].node()) {
    if (auto value = node->value<std::int64_t>()) selectors.parameters["FLAVORS"] = *value;
    else throw std::invalid_argument("DMFT FLAVORS must be an integer");
  }
  return compose(selectors, std::filesystem::absolute(file).parent_path());
}
run_configuration load_run(const std::filesystem::path& file) {
  const auto schema = schema_for_run(file);
  auto run = load_run_configuration(file, schema);
  validate(run, toml::parse(schema));
  return run;
}
}
