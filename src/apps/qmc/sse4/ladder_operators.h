// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
// Bespoke discovery of a requested Green observable independent of a zero
// Hamiltonian coefficient. Retained nonzero Hamiltonian operator choices are
// preserved; the named ALPS spin/boson basis supplies actual evaluated ladders.
#pragma once
#include <alps/model.h>
#include <set>
#include <string>

namespace alps::dirloop_sse {
template<class Model>
std::set<std::string> ladder_operator_names(const Model &model,const Parameters &parameters,
    unsigned type,const std::set<std::string> &hamiltonian_operators) {
    if(!parameters.value_or_default("MEASURE[Green Function]",false))return hamiltonian_operators;
    const auto count=model.basis().site_basis(type).num_states();
    bool raising=false,lowering=false;
    for(const auto &name:hamiltonian_operators) {
        const auto matrix=get_matrix(double(),SiteOperator(name),model.basis().site_basis(type),parameters);
        for(int i=0;i+1<count;++i) {
            raising|=matrix[i][i+1]!=0;lowering|=matrix[i+1][i]!=0;
        }
    }
    if(raising || lowering)return hamiltonian_operators;
    const auto family=parameters.value_or_default("MODEL",std::string("spin"));
    auto selected=hamiltonian_operators;
    if(family=="spin") {selected.insert("Splus");selected.insert("Sminus");}
    else if(family=="hardcore boson" || family=="boson Hubbard") {selected.insert("bdag");selected.insert("b");}
    return selected;
}
}
