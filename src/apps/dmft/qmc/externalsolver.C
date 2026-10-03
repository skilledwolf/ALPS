/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2010 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

/* $Id: externalsolver.C 360 2009-06-01 02:32:00Z gullc $ */

/// @file externalsolver.C
/// @brief implements the external solver
/// @sa ExternalSolver
#include "externalsolver.h"
#include "fouriertransform.h"
#include "U_matrix.h"
#include "bandstructure.h"
#include <alps/cthyb.hpp>
#include <alps/ctint.hpp>
#include <alps/utility/os.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <fstream>
#include <system_error>
#include <type_traits>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <spawn.h>
#include <sys/wait.h>
#include <cerrno>
extern char** environ;
#endif

struct ExternalSolver::invocation_files {
  std::filesystem::path directory, input, run, output;
  invocation_files() {
    directory=std::filesystem::temp_directory_path()/boost::filesystem::unique_path("alps-dmft-%%%%%%%%%%%%").string();
    if(!std::filesystem::create_directory(directory)) throw std::runtime_error("Cannot reserve solver work directory");
    input=directory/"input.h5"; run=directory/"run.toml"; output=directory/"output.h5";
  }
  ~invocation_files() { std::error_code error; std::filesystem::remove_all(directory,error); }
};

namespace {
template<class T>
void check_input(const green_function<T>& green, unsigned nt, unsigned nf) {
  if(green.ntime()!=nt || green.nsite()!=1 || green.nflavor()!=nf)
    throw std::invalid_argument("Impurity Green-function dimensions do not match DMFT parameters");
  require_finite(green,"Impurity Green function");
}
template<class T>
green_function<T> read_green(alps::hdf5::archive& archive, const std::string& path, unsigned nt, unsigned nf) {
  for (const auto& [name, expected] : std::vector<std::pair<std::string,unsigned>>{{"nt",nt},{"ns",1},{"nf",nf}}) {
    const auto dimension = path + "/" + name;
    // Reject non-integers before a conversion can disguise them as the expected value.
    if (!archive.is_data(dimension) || !archive.is_scalar(dimension) || archive.is_complex(dimension) ||
        archive.is_datatype<double>(dimension) || archive.is_datatype<float>(dimension))
      throw std::runtime_error(dimension + ": solver output dimension must be a scalar integer");
    long long stored;
    archive[dimension] >> stored;
    if (stored != expected) throw std::runtime_error(dimension + ": solver output dimension does not match DMFT parameters");
  }
  // Validate extents before the Green-function reader writes into its fixed-size buffers.
  std::vector<std::size_t> shape{nt};
  if constexpr(std::is_same_v<T,std::complex<double>>) shape.push_back(2);
  for(unsigned flavor=0;flavor<nf;++flavor) {
    const auto dataset=path+"/"+std::to_string(flavor)+"/mean/value";
    if(!archive.is_data(dataset) || archive.extent(dataset)!=shape ||
       archive.is_complex(dataset)!=std::is_same_v<T,std::complex<double>>)
      throw std::runtime_error(dataset+": invalid solver output dataset shape or type");
  }
  green_function<T> result(nt,1,nf);
  result.read_hdf5(archive,path);
  require_finite(result,"Impurity solver output "+path);
  return result;
}
void write_run(const std::filesystem::path& file, const alps::run_configuration& run) {
  std::ofstream stream;
  stream.exceptions(std::ios::badbit|std::ios::failbit);
  stream.open(file);
  stream<<alps::format_run_configuration(run);
  stream.close();
}
void launch(const std::filesystem::path& executable,const std::filesystem::path& argument) {
#ifdef _WIN32
  const auto quote=[](const std::wstring& text) {
    std::wstring result=L"\"";
    std::size_t slashes=0;
    for(wchar_t character:text) {
      if(character==L'\\') { ++slashes; continue; }
      result.append(character==L'"'?2*slashes+1:slashes,L'\\');
      result+=character; slashes=0;
    }
    result.append(2*slashes,L'\\'); result+=L'"';
    return result;
  };
  auto command=quote(executable.wstring())+L" "+quote(argument.wstring());
  STARTUPINFOW startup{}; startup.cb=sizeof(startup);
  PROCESS_INFORMATION process{};
  if(!CreateProcessW(executable.c_str(),command.data(),nullptr,nullptr,FALSE,0,nullptr,nullptr,&startup,&process))
    throw std::system_error(GetLastError(),std::system_category(),"Start impurity solver");
  CloseHandle(process.hThread);
  const auto waited=WaitForSingleObject(process.hProcess,INFINITE);
  DWORD status=0,error=0;
  if(waited!=WAIT_OBJECT_0 || !GetExitCodeProcess(process.hProcess,&status)) error=GetLastError();
  CloseHandle(process.hProcess);
  if(error) throw std::system_error(error,std::system_category(),"Wait for impurity solver");
  if(status) throw std::runtime_error("Impurity solver exited with status "+std::to_string(status));
#else
  auto program=executable.string(),file=argument.string();
  char* arguments[]{program.data(),file.data(),nullptr};
  pid_t process;
  const auto error=posix_spawn(&process,program.c_str(),nullptr,nullptr,arguments,environ);
  if(error) throw std::system_error(error,std::generic_category(),"Start impurity solver");
  int status;
  while(waitpid(process,&status,0)<0)
    if(errno!=EINTR) throw std::system_error(errno,std::generic_category(),"Wait for impurity solver");
  if(!WIFEXITED(status) || WEXITSTATUS(status))
    throw std::runtime_error("Impurity solver failed"+(WIFEXITED(status)?
      std::string(" with status ")+std::to_string(WEXITSTATUS(status)):std::string(" after a signal")));
#endif
}

}
  
ExternalSolver::ExternalSolver(const alps::run_configuration& configuration)
  : configuration_(configuration), kind_(alps::dmft::selected_solver(configuration)),
    schema_(alps::dmft::solver_schema(configuration)), delta_(alps::dmft::receives_delta(configuration)) {
  std::filesystem::path name(configuration.execution["solver"].as<std::string>());
  if(kind_==alps::dmft::solver_kind::custom) {
    if(name.is_relative()) name=std::filesystem::path(alps::bin_directory().string())/name;
  } else {
#ifdef _WIN32
    name+=".exe";  // CreateProcessW does not supply an extension.
#endif
    name=std::filesystem::path(alps::bin_directory().string())/name;
  }
  executable_=name;
}

alps::run_configuration ExternalSolver::solver_run(const alps::params& parameters,
                                                   const invocation_files& files) const {
  alps::run_configuration run;
  auto scientific=parameters;
  if(delta_) scientific["MEASURE_time"]=true;  // DMFT needs both endpoints of G(tau).
  run.parameters=alps::dmft::solver_parameters(scientific,schema_);
  run.input=alps::dmft::solver_inputs(configuration_.input,schema_);
  run.execution=alps::select_parameters(configuration_.execution,schema_,"execution");
  run.output["results"]=files.output.string();
  if(delta_) {
    run.input["delta"]=files.input.string();
    run.input["delta_format"]="hdf5";
  } else run.input["g0"]=files.input.string();
  return run;
}
      
void ExternalSolver::write_hybridization(alps::run_configuration& run, const alps::params& parameters,
                                         const itime_green_function_t& delta, alps::hdf5::archive& input) const {
  write_flavor_vectors(input,"/Delta",delta);
  U_matrix interaction(parameters,configuration_.input);
  const double mu=parameters["MU"].as<double>()+interaction.mu_shift();
  const double field=parameters["H"].as<double>();
  run.parameters["MU"]=mu;
  if(field!=0.) {
    std::vector<double> potentials(parameters["FLAVORS"].as<std::size_t>());
    for(std::size_t flavor=0;flavor<potentials.size();++flavor)
      potentials[flavor]=mu+(flavor%2?field:-field);
    input["/MUvector"]<<potentials;
    run.input["chemical_potential"]=run.input["delta"];
    run.input["chemical_potential_format"]="hdf5";
  } 
  if(configuration_.input.exists("interaction_matrix")) {
    std::vector<double> matrix(std::size_t(interaction.nf())*interaction.nf());
    for(unsigned i=0;i<interaction.nf();++i)
      for(unsigned j=0;j<interaction.nf();++j) matrix[i*interaction.nf()+j]=interaction(i,j);
    input["/Umatrix"]<<matrix;
    run.input["interaction_matrix"]=run.input["delta"];
    run.input["interaction_format"]="hdf5";
  }
}

void ExternalSolver::call(alps::run_configuration run, const invocation_files& files) const {
  if(kind_==alps::dmft::solver_kind::hybridization) alps::cthyb::prepare_run(run);
  else if(kind_==alps::dmft::solver_kind::interaction) alps::ctint::prepare_run(run);
  else {
    // A custom schema receives the subset of the built-in protocol it declares.
    run.parameters=alps::select_parameters(run.parameters,schema_);
    run.input=alps::select_parameters(run.input,schema_,"input");
    run.output=alps::select_parameters(run.output,schema_,"output");
    run=alps::resolve_run_configuration(run,schema_);
    if(!run.input.exists(delta_?"delta":"g0"))
      throw std::invalid_argument("Custom solver schema must declare its scientific input path");
  }
  write_run(files.run,run);
  launch(executable_,files.run);
  if(!std::filesystem::is_regular_file(files.output))
    throw std::runtime_error("The external impurity solver did not write its configured results");
}
      
ImpuritySolver::result_type ExternalSolver::solve(const itime_green_function_t& green,
                                                 const alps::params& parameters) {
  // The tau loop passes the Bethe-lattice hybridization Delta(tau) = t^2 G(tau).
  if(!delta_) throw std::logic_error("The DMFT tau loop requires a hybridization-function solver");
  const auto nt=parameters["N"].as<unsigned>()+1,nf=parameters["FLAVORS"].as<unsigned>();
  check_input(green,nt,nf);
  SemicircleBandstructure band(parameters);
  itime_green_function_t delta(nt,1,nf);
  const bool afm=parameters["ANTIFERROMAGNET"].as<bool>();
  for(unsigned flavor=0;flavor<nf;++flavor)
    for(unsigned time=0;time<nt;++time)
      delta(time,flavor)=band.second_moment(flavor)*green(time,afm?(flavor%2?flavor-1:flavor+1):flavor);
  invocation_files files;
  auto run=solver_run(parameters,files);
  {
    alps::hdf5::archive input(files.input.string(),"w");
    write_hybridization(run,parameters,delta,input);
  }
  call(run,files);
  alps::hdf5::archive output(files.output.string(),"r");
  return read_green<double>(output,"/G_tau",nt,nf);
}
  
MatsubaraImpuritySolver::result_type ExternalSolver::solve_omega(
    const matsubara_green_function_t& green,const alps::params& parameters) {
  const auto nw=parameters["NMATSUBARA"].as<unsigned>(),nt=parameters["N"].as<unsigned>()+1;
  const auto nf=parameters["FLAVORS"].as<unsigned>();
  check_input(green,nw,nf);
  invocation_files files;
  auto run=solver_run(parameters,files);
  {
    alps::hdf5::archive input(files.input.string(),"w");
    if(delta_) {
      FFunctionFourierTransformer fourier(parameters);
      matsubara_green_function_t delta_omega(nw,1,nf);
      itime_green_function_t delta_tau(nt,1,nf);
      const auto mu=parameters["MU"].as<double>(),beta=parameters["BETA"].as<double>();
      const auto field=parameters["H"].as<double>();
      for(unsigned flavor=0;flavor<nf;++flavor)
        for(unsigned frequency=0;frequency<nw;++frequency)
          delta_omega(frequency,flavor)=-1./green(frequency,flavor)+std::complex<double>(
            mu+(flavor%2?field:-field),(2.*frequency+1)*boost::math::constants::pi<double>()/beta);
      fourier.backward_ft(delta_tau,delta_omega);
      write_hybridization(run,parameters,delta_tau,input);
    } else write_flavor_vectors(input,"/G0",green);
  }
  call(run,files);
  alps::hdf5::archive output(files.output.string(),"r");
  auto time=read_green<double>(output,"/G_tau",nt,nf);
  matsubara_green_function_t frequency(nw,1,nf);
  if(output.is_group("/G_omega")) frequency=read_green<std::complex<double>>(output,"/G_omega",nw,nf);
  else {
    std::vector<double> density(nf);
    for(unsigned flavor=0;flavor<nf;++flavor) density[flavor]=-time(nt-1,0,0,flavor);
    boost::shared_ptr<FourierTransformer> fourier;
    FourierTransformer::generate_transformer_U(parameters,fourier,density);
    fourier->forward_ft(time,frequency);
    require_finite(frequency,"Impurity Green function from the solver's G(tau)");
  }
  return std::make_pair(frequency,time);
}
