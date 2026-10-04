// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/cthyb.hpp>
#include <alps/solvers.hpp>
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include "../hyb.hpp"
#include "../hybevaluate.hpp"
#include <array>
#include <cmath>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#ifdef ALPS_HAVE_MPI
#include <boost/mpi/environment.hpp>
#endif

namespace {
void check(bool ok,const char* message) { if(!ok) throw std::runtime_error(message); }
template<class F> void rejects(F f,const char* expected) {
  try { f(); } catch(const std::exception& e) {
    if(std::string(e.what()).find(expected)!=std::string::npos) return;
    throw std::runtime_error(std::string("Wrong rejection: ")+e.what());
  }
  throw std::runtime_error(std::string("Accepted invalid configuration: ")+expected);
}
alps::params minimal() {
  alps::params p;
  p["BETA"]=2.; p["N_ORBITALS"]=2; p["N_TAU"]=2; p["N_MEAS"]=1;
  p["THERMALIZATION"]=0; p["SWEEPS"]=10; p["U"]=1.;
  return p;
}

class measured_hybridization : public hybridization {
public:
  using hybridization::hybridization;
  void block(double density, std::array<double,3> const& f, double final_sign=-1.) {
    // Three update signs (+,+,-), followed by a final negative configuration.
    sgn=1.;
    for (std::size_t i=0; i<n_orbitals; ++i) {
      densities[i]=density;
      orders[i]=N_meas;
      G[i]={0., (-2*density+.1)*beta*beta/N_t, 0.};
      F[i]={f[0]*beta*beta/N_t, f[1]*beta*beta/N_t, f[2]*beta*beta/N_t};
    }
    accumulate_G();
    accumulate_order();
    sign=final_sign;
    record_measurement("gw_re_0", std::vector<double>{sign*.75, sign*1.25}, sign);
    record_measurement("g2w_re_0_0", std::vector<double>{sign, 2*sign, 3*sign, 4*sign}, sign);
  }
  auto raw(std::string const& name) const {
    return std::get<alps::alea::batch_acc<double>>(measurements.at(name).accumulator).result();
  }
  void density_sample(double sample_sign) {
    record_measurement("density_0", .25*sample_sign, sample_sign);
  }
};

void measurement_contract(alps::run_configuration run, std::filesystem::path const& directory) {
  run.parameters["N_MEAS"]=3;
  run.parameters["MEASURE_time"]=true;
  run.parameters["MEASURE_freq"]=true;
  run.parameters["N_MATSUBARA"]=2;
  run.parameters["MEASURE_g2w"]=true;
  run.parameters["N_w2"]=2;
  run.parameters["N_W"]=1;
  run.execution["bins"]=8;
  alps::cthyb::prepare_run(run);
  measured_hybridization single(run,0);
  auto empty=single.collect_results();
  check(empty.at("g_0").count()==0 && empty.at("g_0").size()==3,
        "Empty physical time result lost its shape");
  check(empty.at("g2w_re_0_0").count()==0 && empty.at("g2w_re_0_0").size()==4,
        "Empty two-particle result lost its shape");
  single.block(.25,{3.,5.,7.});
  const auto one=single.collect_results();
  check(std::abs(one.at("Sign").mean<double>()(0)-1./3.)<1e-12,
        "N_MEAS sign average changed");
  check(std::abs(one.at("density_0").mean<double>()(0)-.25)<1e-12 &&
        std::abs(one.at("gw_re_0").mean<double>()(0)-.75)<1e-12,
        "Averaged and final-configuration sign cadences were mixed");
  check(!one.at("density_0").stderror<double>().array().isFinite().any(),
        "Single-bin density claims an independent error estimate");
  const auto g=one.at("g_0").mean<double>(), f=one.at("f_0").mean<double>();
  check(std::abs(g(0)+.75)<1e-12 && std::abs(g(2)+.25)<1e-12 &&
        std::abs(f(0)-6.)<1e-12 && std::abs(f(1)-5.)<1e-12 && std::abs(f(2)-14.)<1e-12,
        "Time endpoint conventions were not applied to signed raw blocks");
  measured_hybridization cancelled(run,0);
  cancelled.density_sample(1.); cancelled.density_sample(-1.);
  const auto before=cancelled.raw("density_0");
  rejects([&]{cancelled.collect_results();},"zero average sign");
  check(cancelled.raw("density_0")==before,"Failed signed analysis consumed live measurements");

  measured_hybridization sampled(run,0);
  for (std::size_t i=0; i<37; ++i)
    sampled.block(.1+(i%7)*.1, {3.+i%3, 2.+i%5, 7.+i%4}, i%5 ? 1. : -1.);
  const auto raw_g=sampled.raw("g_0"), raw_f=sampled.raw("f_0");
  const auto results=sampled.collect_results();
  check(sampled.raw("g_0")==raw_g && sampled.raw("f_0")==raw_f,
        "Collection consumed live time measurements");
  check(results.at("g_0").count()==37 && results.at("g2w_re_0_0").count()==37,
        "Collection lost partial-bin samples");
  check((results.at("g2w_re_0_0").mean<double>()-alps::alea::column<double>::LinSpaced(4,1.,4.)).norm()<1e-12 &&
        results.at("g2w_re_0_0").stderror<double>().norm()<1e-12,
        "Two-particle numerator/sign covariance lost a constant physical ratio");

  const auto output=(directory/"time-contract.h5").string();
  {
    alps::hdf5::archive archive(output,"w");
    evaluate_time(results,run.parameters,run.output,archive);
  }
  alps::hdf5::archive archive(output,"r");
  for (auto const& name : {"g_0","f_0"}) {
    auto const& raw=std::string(name)=="g_0" ? raw_g : raw_f;
    auto const& result=results.at(name);
    // Constant block-sign makes the independent weighted-bin covariance exact.
    const auto count=raw.count();
    long double count2=0;
    std::array<long double,3> mean{};
    for (std::size_t b=0; b<raw.num_batches(); ++b) {
      const auto weight=raw.store().count()(b);
      count2+=static_cast<long double>(weight)*weight;
      for (std::size_t i=0;i<3;++i) mean[i]+=3.L*raw.store().batch()(i,b)/count;
    }
    std::vector<double> covariance, errors;
    const std::string path=std::string(name)=="g_0" ? "/G_tau/0/mean/" : "/F_tau/0/mean/";
    archive[path+"covariance"]>>covariance;
    archive[path+"error"]>>errors;
    check(covariance.size()==9 && errors.size()==3,"Published time statistics shape changed");
    for (std::size_t i=0;i<3;++i) {
      check(std::abs(result.mean<double>()(i)-mean[i])<1e-12,"Time mean differs from independent signed bins");
      for (std::size_t j=0;j<3;++j) {
        long double squared=0;
        for (std::size_t b=0;b<raw.num_batches();++b) {
          const auto weight=raw.store().count()(b);
          if (!weight) continue;
          const auto x=3.L*raw.store().batch()(i,b)/weight-mean[i];
          const auto y=3.L*raw.store().batch()(j,b)/weight-mean[j];
          squared+=weight*x*y;
        }
        const auto expected=squared/(count-count2/count)*count2/(static_cast<long double>(count)*count);
        check(std::abs(covariance[i*3+j]-expected)<1e-12,
              "Published G/F endpoint covariance differs from independent signed bins");
      }
      check(std::abs(errors[i]*errors[i]-covariance[i*3+i])<1e-12,
            "Published time error disagrees with covariance diagonal");
    }
  }
  // A filled zero-order orbital has exact n(w=0)=beta, with the same sign as every frequency.
  local_configuration local(run.parameters,run.input,0);
  local.set_zero_order_orbital_occupied(0,true);
  std::vector<double> nnw(2,0.);
  local.measure_nnw(0,nnw,-1.);
  check(nnw[0]==-2. && nnw[1]==0.,"Zero-frequency density correlator omitted its sign");
}
}
int main(int argc,char** argv) {
#ifdef ALPS_HAVE_MPI
  boost::mpi::environment environment(argc,argv);
#endif
  const auto directory=std::filesystem::temp_directory_path()/
    ("alps-cthyb-config-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
  std::filesystem::create_directory(directory);
  struct cleanup { std::filesystem::path path; ~cleanup(){std::error_code error;std::filesystem::remove_all(path,error);} } guard{directory};
  try {
    const auto p=alps::cthyb::prepare_parameters(minimal());
    check(p["MEASURE_time"].as<bool>() && !p["MEASURE_freq"].as<bool>(),"Measurement defaults changed");
    auto bad=minimal(); bad["N_MEAS"]=0;
    rejects([&]{alps::cthyb::prepare_parameters(bad);},"N_MEAS");
    bad=minimal(); bad["BETA"]=0.;
    rejects([&]{alps::cthyb::prepare_parameters(bad);},"BETA");
    bad=minimal(); bad["MEASURE_freq"]=true;
    rejects([&]{alps::cthyb::prepare_parameters(bad);},"N_MATSUBARA");
    bad=minimal(); bad["MEASURE_g2w"]=true; bad["N_w2"]=3; bad["N_W"]=1;
    rejects([&]{alps::cthyb::prepare_parameters(bad);},"even");
    bad["N_w2"]=4; bad["COMPUTE_VERTEX"]=true;
    rejects([&]{alps::cthyb::prepare_parameters(bad);},"MEASURE_freq");
    bad["MEASURE_freq"]=true; bad["N_MATSUBARA"]=1;
    rejects([&]{alps::cthyb::prepare_parameters(bad);},"N_MATSUBARA >=");
    bad["N_MATSUBARA"]=2;
    alps::cthyb::prepare_parameters(bad);
    bad=minimal(); bad["MAX_TIME"]=10;
    rejects([&]{alps::cthyb::prepare_parameters(bad);},"unknown key");
    bad=minimal(); bad["N_ORBITALS"]=50000;
    rejects([&]{alps::cthyb::prepare_parameters(bad);},"index range");
    const std::string delta_text="0 -0.5 -0.5\n1 -0.5 -0.5\n2 -0.5 -0.5\n";
    std::ofstream(directory/"delta.dat")<<delta_text;
    {
      struct restore_directory {
        std::filesystem::path previous=std::filesystem::current_path();
        ~restore_directory(){std::error_code error;std::filesystem::current_path(previous,error);}
      } restore;
      std::filesystem::current_path(directory);
      const auto base=std::filesystem::current_path();
      const auto input_unchanged=[&]{
        std::ifstream input(base/"delta.dat");
        check(std::string(std::istreambuf_iterator<char>(input),{})==delta_text,
              "Rejected path collision changed the scientific input");
      };
      alps::run_configuration relative;
      relative.parameters=minimal(); relative.input["delta"]="delta.dat";
      relative.output["results"]="delta.dat";
      rejects([&]{alps::cthyb::prepare_run(relative);},"must not replace input.delta");
      input_unchanged();
      std::error_code symlink_error;
      std::filesystem::create_symlink("delta.dat","delta-alias.dat",symlink_error);
      if(!symlink_error){
        relative.output["results"]="delta-alias.dat";
        rejects([&]{alps::cthyb::prepare_run(relative);},"must not replace input.delta");
        input_unchanged();
      } else {
#ifdef _WIN32
        std::cout<<"Symlink collision check unavailable: "<<symlink_error.message()<<'\n';
#else
        throw std::filesystem::filesystem_error("Create input alias",base/"delta-alias.dat",symlink_error);
#endif
      }
      relative.output["results"]="relative-result.h5";
      alps::cthyb::prepare_run(relative);
      check(relative.input["delta"].as<std::string>()==(base/"delta.dat").string(),
            "Prepared scientific input must be anchored to the working directory");
      check(relative.output["results"].as<std::string>()==(base/"relative-result.h5").string(),
            "Prepared results must be anchored to the working directory");
      std::filesystem::create_directory(base/"other");
      std::filesystem::current_path(base/"other");
      alps::cthyb::prepare_run(relative);
      check(relative.input["delta"].as<std::string>()==(base/"delta.dat").string() &&
            relative.output["results"].as<std::string>()==(base/"relative-result.h5").string(),
            "A working-directory change reinterpreted prepared paths");
      input_unchanged();
    }
    alps::run_configuration run;
    run.parameters=minimal(); run.input["delta"]=(directory/"delta.dat").string();
    run.output["results"]=(directory/"result.h5").string();
    alps::cthyb::prepare_run(run);
    measurement_contract(run,directory);
    for (auto bins : {0,1,3,127}) {
      auto invalid_bins=run;
      invalid_bins.execution["bins"]=bins;
      rejects([&]{alps::cthyb::prepare_run(invalid_bins);},"bins");
    }
    check(run.application=="cthyb" && run.schema_version==1,"Application identity must come from schema");
    check(run.execution["time_limit"].as<int>()==0,"Finite-sweep run must not require a time limit");
    auto invalid=run; invalid.output["results"]=invalid.input["delta"];
    rejects([&]{alps::cthyb::prepare_run(invalid);},"must not replace input.");
    invalid=run; invalid.output["text"]=true;
    invalid.output["text_directory"]=directory.string();
    invalid.output["results"]=(directory/"simulation.dat").string();
    rejects([&]{alps::cthyb::prepare_run(invalid);},"replace HDF5 results");
    std::ofstream(directory/"positive.dat")<<"0 0.5 -0.5\n1 -0.5 -0.5\n2 -0.5 -0.5\n";
    invalid=run; invalid.input["delta"]=(directory/"positive.dat").string();
    rejects([&]{alps::cthyb::prepare_run(invalid);},"nonpositive");
    std::ofstream(directory/"short.dat")<<"0 -0.5 -0.5\n";
    invalid.input["delta"]=(directory/"short.dat").string();
    rejects([&]{alps::cthyb::prepare_run(invalid);},"malformed");
    std::ofstream(directory/"coordinate.dat")<<"0 -0.5 -0.5\n0.5 -0.5 -0.5\n2 -0.5 -0.5\n";
    invalid.input["delta"]=(directory/"coordinate.dat").string();
    rejects([&]{alps::cthyb::prepare_run(invalid);},"coordinate grid");
    {
      alps::hdf5::archive ar((directory/"delta.h5").string(),"w");
      ar["/Delta_0"]<<std::vector<double>{-0.5,-0.5,-0.5};
      ar["/Delta_1"]<<std::vector<double>{-0.5,-0.5};
    }
    invalid=run; invalid.input["delta"]=(directory/"delta.h5").string(); invalid.input["delta_format"]="hdf5";
    rejects([&]{alps::cthyb::prepare_run(invalid);},"length 3");
    {
      alps::hdf5::archive ar((directory/"delta.h5").string(),"a");
      ar["/Delta_1"]<<std::vector<double>{-0.5,std::numeric_limits<double>::infinity(),-0.5};
    }
    rejects([&]{alps::cthyb::prepare_run(invalid);},"finite");
    std::ofstream(directory/"mu.dat")<<"0.2 0.3";
    run.input["chemical_potential"]=(directory/"mu.dat").string();
    alps::cthyb::prepare_run(run);
    run.parameters["MEASURE_time"]=false;
    alps::solvers::cthyb(run);
    {
      alps::hdf5::archive ar(run.output["results"].as<std::string>(),"r");
      std::string application;
      ar["/run_config/application"]>>application;
      check(application=="cthyb","Results lost application identity");
      alps::params execution,parameters;
      ar["/run_config/execution"]>>execution;
      ar["/parameters"]>>parameters;
      check(execution["time_limit"].as<int>()==0,"Finite-sweep execution metadata changed");
      check(!parameters.exists("DELTA") && !parameters.exists("MAX_TIME"),"Orchestration leaked into scientific parameters");
      check(ar.is_group("/simulation/results/Sign"),"Simulation did not produce measurements");
    }
    std::cout<<"CT-HYB run configuration contracts passed\n";
    return 0;
  } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
