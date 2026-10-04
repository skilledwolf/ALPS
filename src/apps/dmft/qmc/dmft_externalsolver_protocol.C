// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
// The executable doubles as a protocol-aware fake solver. This exercises argv
// paths with spaces, scientific conventions, malformed results, and cleanup.
#include "externalsolver.h"
#include <alps/run_config.hpp>
#include <alps/hdf5.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/hdf5/complex.hpp>
#include <boost/filesystem/operations.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <cstdint>
#include <boost/math/constants/constants.hpp>

namespace {
constexpr const char* schema=R"toml(application="dmft-contract"
schema_version=1
[parameters.N]
type="int64"
required=true
[parameters.NMATSUBARA]
type="int64"
required=true
[parameters.FLAVORS]
type="int64"
required=true
[parameters.N_TAU]
type="int64"
[parameters.N_MATSUBARA]
type="int64"
[parameters.N_ORBITALS]
type="int64"
[parameters.MU]
type="float64"
required=true
[parameters.MODE]
type="string"
default="time"
[parameters.AUDIT]
type="path"
[input.g0]
type="path"
[input.delta]
type="path"
[input.delta_format]
type="string"
[input.chemical_potential]
type="path"
[input.chemical_potential_format]
type="string"
[input.retarded_interaction]
type="path"
[input.retarded_interaction_format]
type="string"
[input.retarded_interaction_coordinate]
type="string"
[output.results]
type="path"
required=true
[output.text]
type="bool"
[output.text_directory]
type="path"
[execution.seed]
type="int64"
required=true
)toml";
void require(bool value,const std::string& message) {
  if(!value) throw std::runtime_error(message);
}
void close(double actual,double expected,const std::string& message) {
  require(std::abs(actual-expected)<1e-12,message);
}
int fake_solver(const std::filesystem::path& file) {
  const auto run=alps::load_run_configuration(file,schema);
  if(run.parameters.exists("AUDIT")) std::ofstream(run.parameters["AUDIT"].as<std::string>())<<file.string();
  const auto mode=run.parameters["MODE"].as<std::string>();
  if(mode=="fail") return 7;
  if(mode=="empty") return 0;
  const bool delta=run.input.exists("delta");
  const auto data=run.input[delta?"delta":"g0"].as<std::string>();
  require(std::filesystem::path(data).is_absolute(),"solver scientific input must be absolute");
  require(std::filesystem::path(run.output["results"].as<std::string>()).is_absolute(),"solver results must be absolute");
  require(run.execution["seed"].as<int>()==17,"execution seed lost");
  require(!run.parameters.exists("CONVERGED"),"driver settings leaked into solver parameters");
  require(!run.output.exists("text") && !run.output.exists("text_directory"),
          "driver output settings leaked into the solver run");
  const auto nt=run.parameters["N"].as<unsigned>()+1;
  const auto nw=run.parameters["NMATSUBARA"].as<unsigned>();
  const auto nf=run.parameters["FLAVORS"].as<unsigned>();
  {
    alps::hdf5::archive input(data,"r");
    require(!input.is_group("/parameters"),"scientific input archive contains legacy run parameters");
    if(delta) {
      require(run.input.exists("retarded_interaction"),"retarded interaction input was dropped");
      const std::filesystem::path retarded(run.input["retarded_interaction"].as<std::string>());
      require(retarded.is_absolute() && retarded.filename()=="retarded kernel.dat",
              "retarded interaction path was changed");
      require(run.input["retarded_interaction_format"].as<std::string>()=="text" &&
              run.input["retarded_interaction_coordinate"].as<std::string>()=="tau",
              "retarded interaction format or coordinate was dropped");
      std::ifstream kernel(retarded);
      const std::string contents((std::istreambuf_iterator<char>(kernel)),std::istreambuf_iterator<char>());
      require(contents=="0 0 0\n1 0.01 0\n2 0 0\n","retarded interaction data was changed");
      require(run.parameters["N_TAU"].as<unsigned>()==nt-1 &&
              run.parameters["N_ORBITALS"].as<unsigned>()==nf &&
              run.parameters["N_MATSUBARA"].as<unsigned>()==nw,"solver grid/flavor mapping lost");
      close(run.parameters["MU"].as<double>(),1.25,"interaction chemical-potential shift lost");
      std::vector<double> chemical_potential;
      input["/MUvector"]>>chemical_potential;
      require(chemical_potential.size()==2,"field-vector shape");
      close(chemical_potential[0],0.85,"even flavor field sign");
      close(chemical_potential[1],1.65,"odd flavor field sign");
      std::vector<double> values;
      input["/Delta_0"]>>values;
      close(values[0],-0.15,"Bethe second moment or AFM flavor exchange lost");
      input["/Delta_1"]>>values;
      close(values[0],-0.05,"Bethe second moment or AFM flavor exchange lost");
    } else {
      close(run.parameters["MU"].as<double>(),0.25,"bare solver chemical potential changed");
      std::vector<std::complex<double>> values;
      input["/G0_0"]>>values;
      close(values[1].real(),0.2,"complex bare data changed");
      close(values[1].imag(),-0.3,"complex bare data changed");
    }
  }
  itime_green_function_t time(mode=="badshape"?nt+1:nt,1,nf);
  for(unsigned flavor=0;flavor<nf;++flavor)
    for(unsigned i=0;i<time.ntime();++i) time(i,flavor)=-0.5+0.25*flavor;
  if(mode=="nonfinite") time(1,0)=std::numeric_limits<double>::quiet_NaN();
  alps::hdf5::archive output(run.output["results"].as<std::string>(),"w");
  time.write_hdf5(output,"/G_tau");
  if(mode=="badextent") output["/G_tau/0/mean/value"]<<std::vector<double>(nt+1,-0.5);
  // Keep the payload shape correct: only the dimension metadata is corrupt.
  // These cases must fail before a scalar conversion or raw-buffer read can
  // disguise a malformed dimension as the expected value.
  if(mode=="floatdimension") output["/G_tau/nt"]<<double(nt)+0.5;
  if(mode=="vectordimension") output["/G_tau/nt"]<<std::vector<unsigned>{nt};
  if(mode=="widedimension") output["/G_tau/nt"]<<(std::uint64_t(1)<<32)+nt;
  if(mode=="negativedimension") output["/G_tau/nt"]<<-static_cast<long long>(nt);
  if(mode=="omega") {
    matsubara_green_function_t frequency(nw,1,nf);
    for(unsigned flavor=0;flavor<nf;++flavor)
      for(unsigned i=0;i<nw;++i) frequency(i,flavor)={double(flavor),-1./(i+1)};
    frequency.write_hdf5(output,"/G_omega");
  }
  return 0;
}
// reference: "exact" (U=0 equals G0), "sampled" (analytic within noise) or "retarded".
int verify_dmft(const std::filesystem::path& file,const std::string& reference) {
  alps::hdf5::archive archive(file.string(), "r");
  alps::run_configuration run;
  archive["/run_config"]>>run;
  require(run.application=="dmft", "wrong saved application identity");
  require(run.execution["seed"].as<int>()==42 && run.execution["max_iterations"].as<int>()==2,
          "driver execution provenance lost");
  require(!run.parameters.exists("SOLVER") && !run.parameters.exists("MAX_TIME") &&
          !run.parameters.exists("SEED"), "orchestration settings leaked into scientific parameters");
  if(reference=="retarded")
    require(run.input.exists("retarded_interaction") &&
            run.input["retarded_interaction_format"].as<std::string>()=="text" &&
            run.input["retarded_interaction_coordinate"].as<std::string>()=="tau",
            "retarded interaction provenance lost");
  require(!archive.is_group("/simulation/iteration/3"), "iteration limit was ignored");
  const auto nt=run.parameters["N"].as<unsigned>()+1;
  const auto nw=run.parameters["NMATSUBARA"].as<unsigned>();
  const auto nf=run.parameters["FLAVORS"].as<unsigned>();
  const auto beta=run.parameters["BETA"].as<double>(),t=run.parameters["t"].as<double>();
  const bool omega=run.execution["loop"].as<std::string>()=="omega";
  for(unsigned iteration=1;iteration<=2;++iteration) {
    const auto base="/simulation/iteration/"+std::to_string(iteration)+"/results/";
    require(archive.is_group(base+"G_tau"), "missing scientific DMFT iteration");
    for(unsigned flavor=0;flavor<nf;++flavor) {
      unsigned stored_nt,stored_ns,stored_nf;
      archive[base+"G_tau/nt"]>>stored_nt;
      archive[base+"G_tau/ns"]>>stored_ns;
      archive[base+"G_tau/nf"]>>stored_nf;
      require(stored_nt==nt && stored_ns==1 && stored_nf==nf, "wrong time-grid/flavor metadata");
      std::vector<double> time;
      const auto time_path=base+"G_tau/"+std::to_string(flavor)+"/mean/value";
      require(archive.extent(time_path)==std::vector<std::size_t>{nt}, "wrong real time data extent");
      archive[time_path]>>time;
      require(time.size()==nt, "wrong time data length");
      for(double value:time) require(std::isfinite(value), "nonfinite scientific time output");
      close(time.front()+time.back(),-1., "fermion endpoint discontinuity changed");
      require(time.back()<=0. && time.back()>=-1., "density endpoint outside [0,1]");
      if(omega) {
        std::vector<std::complex<double>> frequency,bare;
        const auto frequency_path=base+"G_omega/"+std::to_string(flavor)+"/mean/value";
        require(archive.is_complex(frequency_path) &&
                archive.extent(frequency_path)==std::vector<std::size_t>{nw}, "wrong complex frequency shape");
        archive[frequency_path]>>frequency;
        archive[base+"G0_omega/"+std::to_string(flavor)+"/mean/value"]>>bare;
        require(frequency.size()==nw && bare.size()==nw, "wrong Matsubara data length");
        for(unsigned i=0;i<nw;++i) {
          require(std::isfinite(frequency[i].real()) && std::isfinite(frequency[i].imag()),
                  "nonfinite scientific frequency output");
          if(reference!="retarded") {
            const auto w=(2.*i+1)*boost::math::constants::pi<double>()/beta;
            const std::complex<double> analytic(0., -2./(w+std::sqrt(w*w+4*t*t)));
            const auto tolerance=reference=="exact"?1e-10:0.15;
            require(std::abs(frequency[i]-analytic)<tolerance, "U=0 Bethe Green function differs from analytic reference");
          }
          if(reference=="exact") require(std::abs(frequency[i]-bare[i])<1e-12,
                                          "noninteracting Green function differs from its bare input");
        }
      }
    }
  }
  require(archive.is_group("/simulation/results/G_tau"), "missing final scientific result");
  return 0;
}
struct work_directory {
  std::filesystem::path path=std::filesystem::temp_directory_path()/
    boost::filesystem::unique_path("alps-dmft-contract-%%%%%%%%%%%%").string();
  work_directory() { std::filesystem::create_directory(path); }
  ~work_directory() { std::error_code error; std::filesystem::remove_all(path,error); }
};
void check_cleanup(const std::filesystem::path& audit) {
  std::ifstream stream(audit);
  std::string file;
  std::getline(stream,file);
  require(!file.empty(),"child was not launched through its TOML protocol");
  require(!std::filesystem::exists(std::filesystem::path(file).parent_path()),"solver invocation files leaked");
}
}
int main(int argc,char** argv) {
  try {
    if(argc==4 && std::string(argv[1])=="--verify") return verify_dmft(argv[2],argv[3]);
    if(argc==2) return fake_solver(argv[1]);
    work_directory work;
    const auto executable=work.path/"solver with spaces.exe";
    std::filesystem::copy_file(std::filesystem::absolute(argv[0]),executable);
    std::filesystem::permissions(executable,std::filesystem::perms::owner_exec,
                                std::filesystem::perm_options::add);
    const auto schema_file=work.path/"solver schema.toml",audit=work.path/"last invocation.txt";
    std::ofstream(schema_file)<<schema;
    alps::run_configuration run;
    run.input["solver_schema"]=schema_file.string();
    run.execution["solver"]=executable.string();
    run.execution["solver_input"]="g0";
    run.execution["seed"]=17;
    run.output["text"]=true;
    run.output["text_directory"]=work.path.string();
    auto& p=run.parameters;
    p["N"]=4; p["NMATSUBARA"]=3; p["SITES"]=1; p["FLAVORS"]=2;
    p["MU"]=0.25; p["U"]=2.0; p["BETA"]=2.0; p["H"]=0.4; p["t"]=0.5;
    p["ANTIFERROMAGNET"]=true; p["CONVERGED"]=0.001; p["AUDIT"]=audit.string();
    p["MODE"]="omega";
    matsubara_green_function_t frequency(3,1,2);
    for(unsigned flavor=0;flavor<2;++flavor)
      for(unsigned i=0;i<3;++i) frequency(i,flavor)={0.1*(i+1),-0.3};
    ExternalSolver bare(run);
    auto pair=bare.solve_omega(frequency,p);
    close(pair.first(1,1).real(),1.,"frequency output flavor ordering");
    close(pair.first(1,1).imag(),-0.5,"frequency output sign");
    close(pair.second(4,1),-0.25,"time output endpoint returned with frequency output");
    check_cleanup(audit);
    p["MODE"]="time";
    run.execution["solver_input"]="delta";
    const auto retarded=work.path/"retarded kernel.dat";
    std::ofstream(retarded)<<"0 0 0\n1 0.01 0\n2 0 0\n";
    run.input["retarded_interaction"]=retarded.string();
    run.input["retarded_interaction_format"]="text";
    run.input["retarded_interaction_coordinate"]="tau";
    ExternalSolver hybridization(run);
    itime_green_function_t time(5,1,2);
    for(unsigned i=0;i<5;++i) { time(i,0)=-0.2-0.01*i; time(i,1)=-0.6+0.01*i; }
    auto result=hybridization.solve(time,p);
    close(result(4,0),-0.5,"time output endpoint convention");
    close(result(4,1),-0.25,"time output flavor ordering");
    check_cleanup(audit);
    for(const auto* mode:{"badshape","badextent","floatdimension","vectordimension",
                          "widedimension","negativedimension","nonfinite","fail","empty"}) {
      p["MODE"]=mode;
      bool rejected=false;
      try { hybridization.solve(time,p); } catch(const std::exception&) { rejected=true; }
      require(rejected,std::string("invalid solver result accepted: ")+mode);
      check_cleanup(audit);
    }
    require(std::filesystem::exists(retarded),"retarded interaction input was removed");
    require(std::filesystem::exists(schema_file),"solver schema was removed");
    std::cout<<"DMFT TOML protocol, physical conventions, result validation, and cleanup passed\n";
    return 0;
  } catch(const std::exception& error) { std::cerr<<error.what()<<'\n'; return 1; }
}
