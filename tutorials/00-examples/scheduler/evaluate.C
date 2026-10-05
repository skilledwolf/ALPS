/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2002-2006 by Matthias Troyer <troyer@itp.phys.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/alea/batch.hpp>
#include <alps/alea/hdf5.hpp>
#include <alps/alea/transform.hpp>
#include <alps/alea/transformer.hpp>
#include <alps/run_config.hpp>
#include <iostream>

constexpr char schema[] = R"toml(
application = "ising-evaluate"
schema_version = 1
[parameters]
[input.results]
type = "path"
required = true
[output.results]
type = "path"
required = true
[execution]
)toml";

struct binder : alps::alea::transformer<double> {
    explicit binder(size_t size) : size_(size) {}
    size_t in_size() const override { return size_; }
    size_t out_size() const override { return 1; }
    alps::alea::column<double> operator()(alps::alea::column<double> const& x) const override {
        if (!(x(2)>0)) throw std::domain_error("Binder analysis needs positive <m^2> in every jackknife sample");
        return alps::alea::column<double>{BINDER_FACTOR*x(3)/(x(2)*x(2))};
    }
private:
    size_t size_;
};

int main(int argc, char** argv) {
    try {
        bool validate = false;
        std::string file;
        for (int i=1; i<argc; ++i) {
            const std::string argument(argv[i]);
            if (argument == "--schema") { std::cout << schema; return 0; }
            if (argument == "--help" || argument == "-h") {
                std::cout << "Usage: " << argv[0] << " [--validate] run.toml | --schema\n"; return 0;
            }
            if (argument == "--validate") { validate = true; continue; }
            if (argument.empty() || argument.front()=='-' || !file.empty())
                throw std::invalid_argument("Expected one TOML run file");
            file = argument;
        }
        if (file.empty()) throw std::invalid_argument("Expected one TOML run file");
        const auto run = alps::load_run_configuration(file,schema);
        alps::hdf5::archive input(run.input["results"].as<std::string>(),"r");
        alps::run_configuration source;
        alps::params parameters;
        input["/run_config"] >> source; input["/parameters"] >> parameters;
        if (source.application!="scheduler-ising" || source.schema_version!=1
                || !parameters.value_or<bool>("CORRELATIONS",false))
            throw std::invalid_argument("Expected native scheduler Ising chain results with joint moments");
        alps::alea::batch_result<double> joint;
        alps::alea::hdf5_serializer raw(input,"/simulation");
        deserialize(raw,"joint",joint);
        if (joint.size()!=4+parameters["L"].as<size_t>() || !(joint.observations()>1))
            throw std::invalid_argument("Joint chain moments and at least two effective batches are required");
        const auto result = alps::alea::transform(alps::alea::jackknife_prop(),binder(joint.size()),joint);
        Eigen::MatrixXd select = Eigen::MatrixXd::Zero(joint.size()-4,joint.size());
        select.rightCols(joint.size()-4).setIdentity();
        const auto correlations = alps::alea::transform(alps::alea::jackknife_prop(),
            alps::alea::linear_transformer<double>(select),joint);
        if (validate) { std::cout << "Valid Ising analysis\n"; return 0; }
        alps::hdf5::save_checkpoint(run.output["results"].as<std::string>(),[&](auto& ar) {
            ar["/parameters"] << parameters; ar["/run_config"] << run;
            ar["/analysis/binder_factor"] << double(BINDER_FACTOR);
            alps::alea::hdf5_serializer evidence(ar,"/simulation"), output(ar,"/simulation/results");
            serialize(evidence,"joint",joint);
            serialize(output,"Correlations",correlations);
            serialize(output,"Binder cumulant of Magnetization",result);
        });
        std::cout << "Correlations: " << correlations.mean().transpose() << " +/- " << correlations.stderror().transpose()
                  << "\nBinder cumulant of Magnetization: " << result.mean()(0) << " +/- " << result.stderror()(0) << '\n';
        return 0;
    } catch (std::exception const& error) {
        std::cerr << error.what() << '\n'; return 1;
    }
}
