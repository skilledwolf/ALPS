/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1994-2009 by Matthias Troyer <troyer@comp-phys.org>,
*                            Andreas Honecker <ahoneck@uni-goettingen.de>
* Modifications (C) 2026 ALPS Collaboration
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// Thermodynamics from a fulldiag spectrum. The HDF5 result supplies the
// energies, quantum numbers and eigenstate measurements; the plots are written
// to <prefix>.plot.*.xml and <prefix>.measurements.*.plot.xml.
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/params.hpp>
#include <alps/plot.h>
#include <alps/utility/copyright.hpp>
#include <boost/tuple/tuple.hpp>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <sstream>

namespace {

void print_usage(std::ostream& out, const char* pname)
{
  out << "Usage:\n"
      << pname << " [--T_MIN ...] [--T_MAX ...] [--DELTA_T ...] [--H_MIN ...] [--H_MAX ... ] [--DELTA_H ... ] [--versus h] [--DENSITIES ...] results.h5 ...\n"
      << "or:\n"
      << pname << " --couple mu [--T_MIN ...] [--T_MAX ...] [--DELTA_T ...] [--MU_MIN ...] [--MU_MAX ... ] [--DELTA_MU ...] [--versus mu] [--DENSITIES ...] results.h5 ...\n"
      << "\nOptions override the run's parameters. Without T_MAX and DELTA_T, T or beta\n"
      << "selects one temperature; T_MIN defaults to DELTA_T.\n"
      << "  -h, --help     produce help message\n"
      << "  -l, --license  print license conditions\n";
}

struct sector {
  std::map<std::string, std::string> quantumnumbers;
  std::vector<double> energies;
  std::map<std::string, std::vector<double> > averages;               // [state]
  std::map<std::string, std::vector<std::vector<double> > > profiles; // [state][label]
};

struct spectrum {
  alps::params parameters;
  double sites = 0;
  std::vector<sector> sectors;
  std::map<std::string, std::vector<std::string> > labels;
};

spectrum read(std::filesystem::path const& file)
{
  alps::hdf5::archive ar(file.string());
  if (!ar.is_data("/spectrum/number_of_sites"))
    throw std::invalid_argument(file.string() + " is not a fulldiag result of this ALPS version; rerun fulldiag");
  spectrum s;
  ar["/parameters"] >> s.parameters;
  std::uint64_t sites;
  ar["/spectrum/number_of_sites"] >> sites;
  s.sites = double(sites);
  for (auto const& id : ar.list_children("/spectrum/sectors")) {
    const std::string path = "/spectrum/sectors/" + id;
    sector sec;
    ar[path + "/energies"] >> sec.energies;
    if (ar.is_group(path + "/quantumnumbers"))
      for (auto const& name : ar.list_children(path + "/quantumnumbers"))
        ar[path + "/quantumnumbers/" + name] >> sec.quantumnumbers[name];
    if (ar.is_group(path + "/results"))
      for (auto const& segment : ar.list_children(path + "/results")) {
        const std::string observable = path + "/results/" + segment, name = ar.decode_segment(segment);
        if (ar.is_data(observable + "/labels")) {
          ar[observable + "/mean/value"] >> sec.profiles[name];
          ar[observable + "/labels"] >> s.labels[name];
          for (auto const& values : sec.profiles[name])
            if (values.size() != s.labels[name].size())
              throw std::invalid_argument("Inconsistent labels for " + name);
        }
        else
          ar[observable + "/mean/value"] >> sec.averages[name];
      }
    for (auto const& [name, values] : sec.averages)
      if (values.size() != sec.energies.size())
        throw std::invalid_argument("Inconsistent eigenstate count for " + name);
    for (auto const& [name, values] : sec.profiles)
      if (values.size() != sec.energies.size())
        throw std::invalid_argument("Inconsistent eigenstate count for " + name);
    s.sectors.push_back(std::move(sec));
  }
  if (s.sectors.empty())
    throw std::invalid_argument(file.string() + " contains no spectrum");
  return s;
}

double number(std::string const& text, std::string const& name)
{
  // Quantum numbers are written as integers or halves, e.g. "-3/2".
  std::size_t used = 0, slash = text.find('/');
  double value;
  try {
    value = std::stod(text.substr(0, slash), &used);
    if (used == slash || (slash == std::string::npos && used == text.size()))
      return slash == std::string::npos ? value : value / std::stod(text.substr(slash + 1));
  } catch (std::exception const&) {}
  throw std::invalid_argument("Invalid number for " + name + ": " + text);
}

// Ensemble averages and the eigenstate measurements weighted alike.
struct averages {
  double energy, free_energy, specific_heat, entropy, magnetization, susceptibility;
  std::map<std::string, double> scalars;
  std::map<std::string, std::vector<double> > profiles;
};

class evaluation {
public:
  evaluation(spectrum const& s, std::string const& conserved, double field0)
    : s_(s), field0_(field0), conserved_(true)
  {
    for (auto const& sec : s.sectors) {
      auto found = sec.quantumnumbers.find(conserved);
      conserved_ = conserved_ && found != sec.quantumnumbers.end();
      values_.push_back(found == sec.quantumnumbers.end() ? 0. : number(found->second, conserved));
    }
    if (!conserved_)
      values_.assign(values_.size(), 0.);
  }

  bool conserved() const { return conserved_; }

  averages operator()(double temperature, double field) const
  {
    const double beta = 1. / temperature, shift = field - field0_;
    double ground = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < s_.sectors.size(); ++i)
      for (double e : s_.sectors[i].energies)
        ground = std::min(ground, e - values_[i] * shift);
    double z = 0, en = 0, en2 = 0, m = 0, m2 = 0;
    averages a;
    for (std::size_t i = 0; i < s_.sectors.size(); ++i) {
      auto const& sec = s_.sectors[i];
      for (std::size_t j = 0; j < sec.energies.size(); ++j) {
        const double value = sec.energies[j] - values_[i] * shift, w = std::exp(-beta * (value - ground));
        z += w;
        en += value * w;
        en2 += value * value * w;
        m += values_[i] * w;
        m2 += values_[i] * values_[i] * w;
        for (auto const& [name, values] : sec.averages)
          a.scalars[name] += w * values[j];
        for (auto const& [name, values] : sec.profiles) {
          auto& sum = a.profiles[name];
          sum.resize(values[j].size());
          for (std::size_t k = 0; k < sum.size(); ++k)
            sum[k] += w * values[j][k];
        }
      }
    }
    a.energy = en / z;
    a.free_energy = -std::log(z) / beta + ground;
    a.specific_heat = (en2 / z - a.energy * a.energy) * beta * beta;
    a.entropy = beta * (a.energy - a.free_energy);
    a.magnetization = m / z;
    a.susceptibility = (m2 / z - a.magnetization * a.magnetization) * beta;
    for (auto& entry : a.scalars)
      entry.second /= z;
    for (auto& entry : a.profiles)
      for (double& value : entry.second)
        value /= z;
    return a;
  }

private:
  spectrum const& s_;
  double field0_;
  bool conserved_;
  std::vector<double> values_;
};

std::vector<double> grid(double low, double high, double step, std::string const& name)
{
  if (!std::isfinite(low) || !std::isfinite(high) || !std::isfinite(step) || step <= 0 || high < low ||
      (high - low) / step > 1000000 || (high > low && low + step == low))
    throw std::invalid_argument("Require finite " + name + "_MIN <= " + name + "_MAX, DELTA_" + name +
                                " > 0 and at most 1000001 points");
  std::vector<double> values;
  for (std::size_t i = 0; i <= std::size_t(std::floor((high - low) / step + .5)); ++i)
    values.push_back(low + i * step);
  return values;
}

std::string legend(std::string const& name, double value)
{
  std::ostringstream out;
  out << name << '=' << value;
  return out.str();
}

void write(std::string const& file, alps::plot::Plot<double> const& plot)
{
  alps::oxstream out{boost::filesystem::path(file)};
  out << plot;
}

void evaluate(std::filesystem::path const& file, std::map<std::string, std::string> const& options)
{
  const spectrum s = read(file);
  auto option = [&](std::string const& name) -> std::optional<double> {
    if (auto found = options.find(name); found != options.end())
      return number(found->second, name);
    if (!s.parameters.exists(name))
      return std::nullopt;
    auto const& value = s.parameters[name];
    return value.isType<std::string>() ? number(value.as<std::string>(), name) : value.as<double>();
  };
  auto word = [&](std::string const& name) -> std::string {
    if (auto found = options.find(name); found != options.end())
      return found->second;
    return s.parameters.exists<std::string>(name) ? s.parameters[name].as<std::string>() : std::string();
  };
  auto flag = [&](std::string const& name) {
    const std::string text = word(name);
    if (text == "true" || text == "false")
      return text == "true";
    if (s.parameters.exists<bool>(name) && !options.count(name))
      return s.parameters[name].as<bool>();
    return option(name).value_or(1.) != 0.;
  };

  // The field couples to Sz by default and to N with --couple mu.
  const bool mu = word("couple") == "mu" || word("couple") == "MU";
  const std::string conserved = mu ? "N" : "Sz", field_name = mu ? "mu" : "h", field_upper = mu ? "MU" : "H";
  const std::string field_label = mu ? "Chemical Potential" : "Magnetic Field";
  const std::string conserved_label = mu ? "Particle number" : "Magnetization";
  const std::string moment = mu ? "Compressibility" : "Uniform Susceptibility";
  const bool measure = flag(mu ? "MEASURE_CHARGE_PROPERTIES" : "MEASURE_MAGNETIC_PROPERTIES");
  const bool densities = flag("DENSITIES");
  const bool versus_field = word("versus") == field_name || word("versus") == field_upper;

  std::vector<double> temperatures;
  if (auto high = option("T_MAX"), step = option("DELTA_T"); high && step)
    temperatures = grid(option("T_MIN").value_or(*step), *high, *step, "T");
  else if (auto t = option("T"))
    temperatures = {*t};
  else if (auto beta = option("beta"))
    temperatures = {1. / *beta};
  else
    throw std::invalid_argument("Specify T_MAX and DELTA_T, T or beta");
  for (double t : temperatures)
    if (!(t > 0) || !std::isfinite(t))
      throw std::invalid_argument("Temperatures must be positive and finite");

  // The spectrum was computed at the run's field; other fields shift each
  // sector by its conserved quantum number.
  const double field0 = s.parameters.exists(field_name) ? s.parameters[field_name].as<double>() : 0.;
  std::vector<double> fields;
  if (auto high = option(field_upper + "_MAX"), step = option("DELTA_" + field_upper); high && step)
    fields = grid(option(field_upper + "_MIN").value_or(0.), *high, *step, field_upper);
  else if (auto found = options.count(field_name) ? options.find(field_name) : options.find(field_upper);
           found != options.end())
    fields = {number(found->second, found->first)};
  else
    fields = {field0};

  const evaluation at(s, conserved, field0);
  const bool magnetic = at.conserved() && measure;
  const double per_site = densities ? s.sites : 1.;
  const std::string dstr = densities ? " Density" : "", sstr = densities ? " per Site" : "";
  const std::string xname = versus_field ? field_label : "Temperature";
  const alps::Parameters parameters(s.parameters);

  struct curve { std::string name, ylabel, suffix; double averages::*value; };
  std::vector<curve> curves{
    {"Energy" + dstr, "Energy" + dstr, "energy", &averages::energy},
    {"Free Energy" + dstr, "Free Energy" + dstr, "free_energy", &averages::free_energy},
    {"Entropy" + dstr, "Entropy" + dstr, "entropy", &averages::entropy},
    {"Specific Heat" + sstr, "Specific Heat" + sstr, "specific_heat", &averages::specific_heat}};
  if (magnetic) {
    auto lower = [](std::string text) {
      for (auto& c : text) c = c == ' ' ? '_' : char(std::tolower(static_cast<unsigned char>(c)));
      return text;
    };
    curves.push_back({conserved_label + sstr, conserved_label + sstr, lower(conserved_label), &averages::magnetization});
    curves.push_back({moment + sstr, moment + sstr, lower(moment), &averages::susceptibility});
  }
  std::vector<alps::plot::Plot<double> > plots;
  for (auto const& c : curves) {
    plots.emplace_back(c.name + " versus " + xname, parameters, false);
    plots.back().set_labels(versus_field ? field_name : "T", c.ylabel);
  }
  // Measurements collect one curve per measurement and label over all points.
  std::map<std::string, alps::plot::Set<double> > scalar_sets;
  std::map<std::string, std::vector<alps::plot::Set<double> > > profile_sets;

  const std::vector<double>& outer = versus_field ? temperatures : fields;
  const std::vector<double>& inner = versus_field ? fields : temperatures;
  for (double o : outer) {
    std::vector<alps::plot::Set<double> > sets(curves.size());
    for (auto& set : sets)
      set << legend(versus_field ? "T" : field_name, o);
    for (double x : inner) {
      const averages a = versus_field ? at(o, x) : at(x, o);
      for (std::size_t c = 0; c < curves.size(); ++c)
        sets[c] << boost::make_tuple(x, a.*(curves[c].value) / per_site);
      for (auto const& [name, value] : a.scalars)
        scalar_sets[name] << boost::make_tuple(x, value);
      for (auto const& [name, values] : a.profiles) {
        auto& set = profile_sets[name];
        set.resize(values.size());
        for (std::size_t k = 0; k < values.size(); ++k)
          set[k] << boost::make_tuple(x, values[k]);
      }
    }
    for (std::size_t c = 0; c < curves.size(); ++c)
      plots[c] << sets[c];
  }

  auto prefix = file;
  prefix.replace_extension();
  if (prefix.extension() == ".out")
    prefix.replace_extension();
  for (std::size_t c = 0; c < curves.size(); ++c)
    write(prefix.string() + ".plot." + curves[c].suffix + ".xml", plots[c]);
  for (auto& [name, set] : scalar_sets) {
    alps::plot::Plot<double> plot("", parameters);
    plot.set_labels(xname, name);
    set << name;
    plot << set;
    write(prefix.string() + ".measurements." + name + ".plot.xml", plot);
  }
  for (auto& [name, sets] : profile_sets) {
    alps::plot::Plot<double> plot("", parameters);
    plot.set_labels(xname, name);
    auto const& labels = s.labels.at(name);
    for (std::size_t k = 0; k < sets.size(); ++k) {
      sets[k] << labels[k];
      plot << sets[k];
    }
    write(prefix.string() + ".measurements." + name + ".plot.xml", plot);
  }
}

} // namespace

int main(int argc, char** argv)
{
  try {
    int i = 1;
    std::map<std::string, std::string> options;
    while (i < argc && argv[i][0] == '-') {
      if (!std::strcmp(argv[i], "--help") || !std::strcmp(argv[i], "-h")) {
        print_usage(std::cout, argv[0]);
        alps::print_copyright(std::cout);
        return 0;
      }
      if (!std::strcmp(argv[i], "--license") || !std::strcmp(argv[i], "-l")) {
        alps::print_license(std::cout);
        return 0;
      }
      if (!std::strcmp(argv[i], "--")) {
        ++i;
        break;
      }
      if (argv[i][1] != '-' || argv[i][2] == '\0') {
        std::cerr << "Unknown option: " << argv[i] << "\n";
        print_usage(std::cerr, argv[0]);
        return 1;
      }
      if (i + 1 >= argc) {
        std::cerr << "Missing value for option: " << argv[i] << "\n";
        print_usage(std::cerr, argv[0]);
        return 1;
      }
      options[argv[i] + 2] = argv[i + 1];
      i += 2;
    }
    if (i >= argc) {
      print_usage(std::cerr, argv[0]);
      return 1;
    }
    for (; i < argc; ++i)
      evaluate(argv[i], options);
    return 0;
  } catch (std::exception const& e) {
    std::cerr << "fulldiag_evaluate: " << e.what() << "\n";
    return 1;
  }
}
