/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2003-2010 by Sergei Isakov <isakov@itp.phys.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef __MEASUREMENT_H__
#define __MEASUREMENT_H__

#include <vector>
#include <valarray>

#include "sse_alg_def.h"

template<typename L, typename M, typename W>
class Measurement {
public:
    typedef L lattice_type;
    typedef M model_type;
    typedef W worker_type;
    
    typedef typename lattice_type::vector_type vector_type;
    typedef typename lattice_type::lat_unit_type lat_unit_type;
    typedef typename model_type::vertex_type vertex_type;
    
    typedef typename worker_type::state_type state_type;
    
    Measurement(lattice_type const& lattice,
            model_type const& model,
            worker_type& worker,
            alps::SymbolTable const& params,
            std::vector<state_type>& state,
            std::vector<Operator> const& opstring) :
        lattice(lattice),
        model(model),
        worker(worker),
        state(state),
        opstring(opstring)
    {
        nsites = lattice.nsites();
        
        skip=params.value_or_default("SKIP",1);
        escale = 1.0 / worker.beta() / nsites;
        sscale = worker.winding_dimension() ? 1.0 / worker.beta() / worker.winding_dimension() : 0.;
        
        worker.add_measurement("Kinetic Energy");
        worker.add_measurement("Kinetic Energy Density");
        
        worker.add_measurement("n");
        worker.add_measurement("n^2");
        worker.add_measurement("n^3");
            
        worker.initialize_site_states();
        worker.create_common_observables();
        
        if (worker.measure_green_function()) {
            green.resize(worker.do_measurement_origin() ?
                        nsites : lattice.ndistances());
            green = 0.0;
        }
        
        if (worker.measure_site_compressibility()) {
            localint.resize(nsites);
            localint2.resize(nsites);
            lasti.resize(nsites);
        }
    }
    
    void set_nworms(unsigned nworms)
    {
        this->nworms = nworms;
    }
    
    std::valarray<double>& get_green()
    {
        return green;
    }
    
    void do_measurement(unsigned nops, unsigned nnonzero, double cc)
    {
        std::vector<std::vector<double> > const& pstates = worker.phys_states_n();
        
        if (worker.measure_site_compressibility()) {
            localint = 0.0;
            localint2 = 0.0;
            lasti = 0;
        }
        
        unsigned noff_diag = 0;
        
        std::vector<double> wns(worker.winding_dimension(), 0.0);
        
        int sign = 1;
        
        for (op_c_iterator op = opstring.begin(); op != opstring.end(); ++op) {
            if (op->vertex_index == IDENTITY)
                continue;
                
            vertex_type const& vertex = model.vertex(op->vertex_index);
            
            if (!vertex.diagonal) {
                // winding numbers
                vector_type const& v = lattice.bond_vector_relative(op->unit_ref);
                int delta = vertex.state[0] < vertex.state[UNIT_SIZE] ? 1 : -1;
                unsigned imax = std::min(unsigned(v.size()), worker.winding_dimension());
                for (unsigned i = 0; i < imax; ++i)
                    wns[i] += delta * v[i];
                
                ++noff_diag;
                
                if (worker.is_signed() && vertex.me > 0.0)
                    sign = -sign;
                    
                if (worker.measure_site_compressibility()) {
                    lat_unit_type const& lat_unit = lattice.lat_units()[op->unit_ref];

                    for (unsigned i = 0; i < UNIT_SIZE; i++) {
                        unsigned si = lat_unit.sites[i];

                        unsigned ii = op - opstring.begin();
                        double p = pstates[lattice.sitei2alps_type(si)][state[si]];
                        localint[si] += p * (ii - lasti[si]);
                        localint2[si] += p * p * (ii - lasti[si]);
                        lasti[si] = ii;

                        // propagate state
                        state[si] = vertex.state[UNIT_SIZE + i];
                    }
                }
            }
        }
        
        if (worker.measure_site_compressibility()) {
            for (unsigned i = 0; i < nsites; ++i) {
                double p = pstates[lattice.sitei2alps_type(i)][state[i]];
                localint[i] += p * (nops - lasti[i]);
                localint2[i] += p * p * (nops - lasti[i]);
            }

            localint = (localint2 + localint * localint) * worker.beta()
                    / double(nops) / double(nops + 1);
        }

        if (worker.do_common_measurements(double(sign), state, localint)) {            
            double e = (cc - nnonzero) * escale * sign;
            worker.record("Energy",e * nsites,sign);
            worker.record("Energy Density",e,sign);
        
            double ke = -double(noff_diag) * escale * sign;
            worker.record("Kinetic Energy",ke * nsites,sign);
            worker.record("Kinetic Energy Density",ke,sign);
        
            double wn = 0.0;
            for (unsigned i = 0; i < worker.winding_dimension(); ++i)
                wn += wns[i] * wns[i];    
            if (worker.winding_dimension()) worker.record("Stiffness",wn * sscale * sign,sign);
        
            double n = double(nnonzero) * sign;
            worker.record("n",n,sign);
            worker.record("n^2",n * nnonzero,sign);
            worker.record("n^3",n * nnonzero * nnonzero,sign);
        
            if (worker.measure_green_function()) {
                double scale = static_cast<double>(nsites) / nworms / skip;
                if (worker.do_measurement_origin())
                    for (unsigned i = 0; i < green.size(); ++i)
                        green[i] *= scale;
                else
                    for (unsigned i = 0; i < green.size(); ++i)
                        green[i] *= scale / worker.distance_mult()[i];
                    
                worker.record("Green's Function",green,sign);

            }
        }
        green=0.;
    }
private:
    Measurement();
    
    lattice_type const& lattice;
    model_type const& model;
    worker_type& worker;
    
    std::vector<state_type>& state;
    std::vector<Operator> const& opstring;
    
    std::valarray<double> green;
    
    std::valarray<double> localint;
    std::valarray<double> localint2;
    std::valarray<unsigned> lasti;
    
    unsigned skip;
    unsigned nsites;
    unsigned nworms;
    
    double escale;
    double sscale;
};

#endif
