/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2009-2010 by Sergei Isakov <isakov@itp.phys.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef __SSE_RUN_H__
#define __SSE_RUN_H__

#include <alps/lattice.h>

#include "../qmc.h"
#include "lattice.h"
#include "model.h"
#include "measurement.h"
#include "sse_alg.h"

class SSE_run : public QMCRun<> {    
public:
    typedef SSE_run self_type;
    typedef QMCRun<> super_type;
    typedef Lattice<self_type> lattice_type;
    typedef Model<self_type, Lattice<self_type> > model_type;
    typedef SSE_alg<lattice_type, model_type, self_type> algorithm_type;
    
    SSE_run(alps::params const& params,size_t bins=128,size_t chain=0) :
        super_type(params,bins,chain), lattice(parms,*this), model(parms,*this,lattice),
        algorithm(lattice,model,*this,parms), chain_(chain),
        nsweeps(params["SWEEPS"].as<uint64_t>()),
        nthermalization(params["THERMALIZATION"].as<uint64_t>()),
        sweeps_done(0), measurements_skip(params["SKIP"].as<unsigned>()), measurements_done(0) {}

    void save(alps::hdf5::archive& ar) const override {
        super_type::save(ar);
        ar["checkpoint/version"] << uint64_t(1);
        ar["checkpoint/chain"] << uint64_t(chain_);
        ar["checkpoint/sweeps"] << sweeps_done;
        ar["checkpoint/measurement_cursor"] << measurements_done;
        algorithm.save(ar);
    }
    void load(alps::hdf5::archive& ar) override {
        alps::params saved;
        uint64_t version,sweeps,chain;
        unsigned cursor;
        ar["/parameters"] >> saved;
        ar["checkpoint/version"] >> version;
        ar["checkpoint/chain"] >> chain;
        ar["checkpoint/sweeps"] >> sweeps;
        ar["checkpoint/measurement_cursor"] >> cursor;
        auto production=sweeps>nthermalization ? sweeps-nthermalization : 0;
        if (version!=1 || chain!=chain_ || checkpoint_parameters(saved)!=checkpoint_parameters(parameters) ||
            production>nsweeps || cursor!=production%measurements_skip)
            throw std::invalid_argument("Invalid directed-loop checkpoint parameters or progress");
        auto state=algorithm.read_checkpoint(ar);
        validate_measurements(ar,production/measurements_skip);
        auto current=parameters;
        super_type::load(ar);
        parameters=std::move(current);
        algorithm.restore(std::move(state));
        sweeps_done=sweeps; measurements_done=cursor;
    }
    void measure() override {} // Measurements occur at the configured sweep interval.
    void update() override { dostep(); }
    double fraction_completed() const override { return work_done(); }
    uint64_t completed_sweeps() const { return sweeps_done; }

    void dostep()
    {
//        if (sweeps_done >= nthermalization + nsweeps)
//            return;
        
        algorithm.do_step();
        
        if (is_thermalized() && ++measurements_done == measurements_skip) {
            algorithm.do_measurement();
            measurements_done = 0;
        }

        ++sweeps_done;
    }
    
    bool is_thermalized() const
    {
        return sweeps_done >= nthermalization;
    }
    
    double work_done() const
    {
        if (is_thermalized())
            return (sweeps_done - nthermalization) / double(nsweeps);
        return 0.;
    }
    
    bool is_thermalization_done(double percentage) const
    {
        return double(sweeps_done) > percentage * nthermalization;
    }
    
    int mrandom_int(int n)
    {
        return super_type::random_int(n);
    }
    
    double mrandom_real()
    {
        return super_type::random_01();
    }
    
    bool is_signed() const
    {
        return super_type::is_signed_;
    }
    
    bool measure_green_function() const
    {
        return super_type::measure_green_function_;
    }
    
    bool do_measurement_origin() const
    {
        if (super_type::measurement_origin_)
            return true;
        else
            return false;
    }
    
    int measurement_origin() const
    {
        return super_type::measurement_origin_.get();
    }
    
    std::vector<unsigned> const& distance_mult() const
    {
        return super_type::distance_mult;
    }
    
    bool measure_site_compressibility() const
    {
        return super_type::measure_site_compressibility_;
    }
    
    bool do_common_measurements(double sign,
        std::vector<state_type> const& state, std::valarray<double> const& localint)
    {
        return super_type::do_common_measurements(sign, state, localint);
    }
    
    void create_common_observables()
    {
        super_type::create_common_observables();
    }
    
    void initialize_site_states()
    {
        super_type::initialize_site_states();
    }
    
    double beta() const
    {
        return super_type::beta;
    }
    
    std::vector<std::vector<double> > const& phys_states_n()
    {
        return super_type::diagonal_matrix_element["n"];
    }
private:
    lattice_type lattice;
    model_type model;
    algorithm_type algorithm;
    
    size_t chain_;
    boost::uint64_t nsweeps;
    boost::uint64_t nthermalization;
    boost::uint64_t sweeps_done;
    unsigned measurements_skip;
    unsigned measurements_done;
};

#endif
