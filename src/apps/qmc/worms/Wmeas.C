/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2005 by Matthias Troyer <troyer@comp-phys.org>,
*                            Simon Trebst <trebst@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include "WRun.h"

void WRun::create_observables()
{
  create_common_observables();

  add_measurement("Statistics",4,false);
  if (winding_dimension()) add_measurement("Winding number^2",winding_dimension());

  if(use_1D_stiffness){   //@#$br: take care of D=1 ---> need a histogram of the winding numbers
     add_measurement("Winding number histogram",3);
  }

  // warn a user if using the D>1 stiffness estimator in 1D
  //if(1==winding_dimension() && !use_1D_stiffness){
  //  std::cout<<" *** Warning: using the <W^2> estimator for the stiffness \n *** Consider switching USE_1D_STIFFNESS=true \n" ;
  //}



  num_chains=0;
  if (chain_kappa) {
    if (!is_charge_model_)
      throw std::invalid_argument("CHAIN_KAPPA requires a charge model");
    std::map<double,int> chains;
    for (auto [it,end]=sites();it!=end;++it) {
      auto const& position=coordinate(*it);
      if (position.size()<2 || !std::isfinite(position[1]))
        throw std::invalid_argument("CHAIN_KAPPA requires finite transverse coordinates on every site");
      chains.emplace(position[1],0);
    }
    for (auto& [position,index]:chains) index=num_chains++;
    chain_size.resize(num_chains);
    for (auto [it,end]=sites();it!=end;++it) {
      chain_number[*it]=chains.at(coordinate(*it)[1]);
      ++chain_size[chain_number[*it]];
    }
    add_measurement("Chain density",num_chains);
    add_measurement("Chain density^2",num_chains);
  }

}

int WRun::get_particle_number() {
  double n=0.;
  for (int i=0;i<num_sites();i++)
    n += QMCRun<>::diagonal_matrix_element["n"][site_type(i)][initial_state(i)];
  return static_cast<int>(n);
}

// Tuning may change the fugacity, but not relative weights within the selected
// particle-number sector. Report energies for the requested Hamiltonian so that
// independent chains with different tuned fugacities remain comparable.
void WRun::set_adjusted_parameter(double value) {
  if (!std::isfinite(value)) throw std::invalid_argument("Nonfinite adjusted parameter");
  auto interactions=hamiltonian_state(false);
  parms[adjust_parameter]=value;
  initialize_hamiltonian();
  if (hamiltonian_state(false)!=interactions)
    throw std::invalid_argument("ADJUST must not change hopping or bond interactions");
  double slope=std::numeric_limits<double>::quiet_NaN(), offset=0.;
  for (auto [it,end]=sites();it!=end && !std::isfinite(slope);++it) {
    auto type=inhomogeneous_site_type(*it);
    auto const& n=QMCRun<>::diagonal_matrix_element.at("n")[site_type(*it)];
    for (size_t state=1;state<n.size();++state) if (n[state]!=n[0]) {
      slope=((site_matrix[type][state]-requested_site_matrix[type][state])-
             (site_matrix[type][0]-requested_site_matrix[type][0]))/(n[state]-n[0]);
      break;
    }
  }
  for (auto [it,end]=sites();it!=end;++it) {
    auto type=inhomogeneous_site_type(*it);
    auto const& n=QMCRun<>::diagonal_matrix_element.at("n")[site_type(*it)];
    auto const& original=requested_site_matrix.at(type);
    auto const& adjusted=site_matrix.at(type);
    double delta0=adjusted[0]-original[0];
    for (size_t state=0;state<n.size();++state) {
      double delta=adjusted[state]-original[state];
      double expected=delta0+(std::isfinite(slope) ? slope*(n[state]-n[0]) : 0.);
      if (!std::isfinite(delta) || std::abs(delta-expected)>1e-10*std::max({1.,std::abs(delta),std::abs(expected)}))
        throw std::invalid_argument("ADJUST must preserve the Hamiltonian within each particle-number sector");
    }
    offset+=delta0-(std::isfinite(slope) ? slope*n[0] : 0.);
  }
  adjustment_energy_shift=offset+(std::isfinite(slope) ? slope*double(parms["NUMBER_OF_PARTICLES"]) : 0.);
}

void WRun::adjustment() {
  if (steps>=thermal_sweeps/4)
    throw std::runtime_error("Particle-number adjustment exceeded 25% of thermalization; increase THERMALIZATION");
  int target=static_cast<int>(parms["NUMBER_OF_PARTICLES"]);
  double correction=static_cast<double>(parms["CORRECTION"]);
  if (!preadjustment_done) {
    if (steps%25) return;
    int n=get_particle_number();
    if (n==target || (corrections_upwards && corrections_downwards)) {
      preadjustment_done=true;
      nob.clear();
    } else {
      if (n<target) ++corrections_upwards;
      else ++corrections_downwards;
      set_adjusted_parameter(double(parms[adjust_parameter])+(n<target ? 10. : -10.)*correction);
    }
    return;
  }
  nob.push_back(get_particle_number());
  if (nob.size()<60) return;
  double n=std::accumulate(nob.begin()+10,nob.end(),0.)/50.;
  nob.clear();
  if (std::abs(target-n)<std::max(1,target)/200.) adjustment_done=true;
  else set_adjusted_parameter(double(parms[adjust_parameter])+(n<target ? correction : -correction));
}

void WRun::make_meas()
{
  std::vector<std::vector<double> > const& matrix_element_n = QMCRun<>::diagonal_matrix_element["n"];

  // measure energy and particle numbers
  double energy=0.;
  std::vector<state_type> local(num_sites());
  std::valarray<double> chain_n(num_chains);

  for (int i=0;i<num_sites();i++) {
    // determine energy
    state_type s=initial_state(i);
    if (nonlocal)
      energy += H0(site(i)) - kinks[i].size()/(2.*beta);
    else
      energy += onsite_energy(s,i) -kinks[i].size()/(2.*beta);

    local[i]=s;
    if (chain_kappa)
      chain_n[chain_number[i]]+=matrix_element_n[site_type(i)][s];
  }

  energy-=adjustment_energy_shift;

  if (!do_common_measurements(Sign,local)) {
    return;
  }

  // measure winding numbers
  bond_iterator bi, bi_end;

  std::valarray<double> winding_number(0., winding_dimension());
  for(boost::tie(bi, bi_end) = bonds(); bi != bi_end; ++bi) {
    int s1 = source(*bi);
    int s2 = target(*bi);
    vector_type v=bond_vector_relative(*bi);


    cyclic_iterator i1 = first_kink(s1);
    cyclic_iterator i2 = first_kink(s2);
    if (i1.valid() && i2.valid()) {
      while(true) {
        if (i1->time() == i2->time()) {
          int hop_dir = (i1->state() > (i1-1)->state() ? -1 : 1);
          for(int d=0; d<winding_dimension(); ++d) {
            winding_number[d]+=hop_dir*v[d];
          }
        }
        if (i1->time()<i2->time()) {
          ++i1;
          if (i1==first_kink(s1))
            break;
        }
        else {
          ++i2;
          if (i2==first_kink(s2))
            break;
        }
      }
    }
  }
  std::valarray<double> winding_number2(0., winding_dimension());
  winding_number2 = winding_number*winding_number;

  if(use_1D_stiffness){
     // @#$br : record the winding numbers = -1, 0, 1
     //         the vector entries are 0,1,2, respectively
     std::valarray<double> winding_numberz(3);
     for(int i=0; i<3;++i){  winding_numberz[i]=0. ;}
     if( fabs( winding_number[0] ) <=1 ) {
          winding_numberz[(int)winding_number[0]+1]+=Sign ;
          record("Winding number histogram",winding_numberz ,Sign); //*Sign;
     }
  }
  else if (winding_dimension()){
     // determine stiffness
     double stiffness=0.;
     for(int d=0; d<winding_dimension(); ++d) {
        stiffness += winding_number[d]*winding_number[d];
     }
     record("Stiffness",stiffness*Sign/(winding_dimension()*beta),Sign);
  }

  double vol = num_sites();

  record("Statistics",stat,Sign);
  winding_number2 *= Sign;
  if (winding_dimension()) record("Winding number^2",winding_number2,Sign);
  record("Energy",energy*Sign,Sign);
  record("Energy Density",energy/vol*Sign,Sign);
  //record("Stiffness",stiffness*Sign/(winding_dimension()*beta),Sign);

  if(chain_kappa) {
    for (int i=0;i<num_chains;++i) chain_n[i]/=chain_size[i];

    std::valarray<double> chain_n2=chain_n*chain_n;

    chain_n  *= Sign;
    chain_n2 *= Sign;

    record("Chain density",chain_n,Sign);
    record("Chain density^2",chain_n2,Sign);
  }

}   // WRun::make_meas

