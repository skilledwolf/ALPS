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
#include <map>
#include <cmath>
#include <stdexcept>

void WRun::create_observables()
{
  create_common_observables();
  
  measurements << alps::RealVectorObservable("Statistics");
  measurements << alps::make_observable(alps::RealVectorObservable("Winding number^2"), is_signed_);

  if (use_1D_stiffness && dimension() != 1)
    throw std::invalid_argument("USE_1D_STIFFNESS requires one spatial dimension");
  if(use_1D_stiffness){   //@#$br: take care of D=1 ---> need a histogram of the winding numbers
     measurements << alps::make_observable(alps::RealVectorObservable("Winding number histogram"), is_signed_); 
  }
  
  // warn a user if using the D>1 stiffness estimator in 1D  
  //if(1==dimension() && !use_1D_stiffness){
  //  std::cout<<" *** Warning: using the <W^2> estimator for the stiffness \n *** Consider switching USE_1D_STIFFNESS=true \n" ;
  //}


 
  if (nonlocal) {
    measurements << alps::RealObservable("Ratio time intervals");
    measurements << alps::RealObservable("Ratio time intervals 2");
    measurements << alps::RealObservable("Time interval");
    measurements << alps::RealVectorObservable("Statistics time intervals");
  }
  
  site_iterator it,end;
  num_chains=0;
  if (chain_kappa) {
    measurements << alps::make_observable(alps::RealVectorObservable("Chain density"), is_signed_);
    measurements << alps::make_observable(alps::RealVectorObservable("Chain density^2"), is_signed_);
    if (!where.empty()) {
      // Coordinate labels identify chains; each chain uses its actual site
      // count, including unequal chains. Never narrow a coordinate to int.
      std::map<double, unsigned> labels;
      for (boost::tie(it,end)=sites();it!=end;++it) {
        const auto& position = coordinate(*it);
        if (position.size() < 2 || !std::isfinite(position[1]))
          throw std::invalid_argument("CHAIN_KAPPA requires finite second site coordinates");
        labels.emplace(position[1], 0);
      }
      for (auto& label : labels) {
        label.second = chain_labels.size(); chain_labels.push_back(label.first);
      }
      num_chains = chain_labels.size(); chain_sizes.assign(num_chains, 0);
      for (boost::tie(it,end)=sites();it!=end;++it) {
        const auto chain = labels.at(coordinate(*it)[1]);
        chain_number[*it] = chain; ++chain_sizes[chain];
      }
    }
  }
}

int WRun::get_particle_number() {
  double n=0.;
  for (int i=0;i<num_sites();i++)
    n += QMCRun<>::diagonal_matrix_element["n"][site_type(i)][initial_state(i)];
  return static_cast<int>(n);
}

void WRun::adjustment() {

  std::cerr << "Starting adjustment\n";
  int number_of_bosons=static_cast<int>(parms["NUMBER_OF_PARTICLES"]);
  double correction=static_cast<double>(parms["CORRECTION"]);

  if ((corrections_upwards==0||corrections_downwards==0)&&steps%25==0) {
    // coarse adjustment
    if (number_of_bosons>get_particle_number()) {
      parms[adjust_parameter] = static_cast<double>(parms[adjust_parameter]) + 10.*correction;
      corrections_upwards++;
    }
    else {
      parms[adjust_parameter] = static_cast<double>(parms[adjust_parameter]) - 10.*correction;
      corrections_downwards++;
    }
    std::cerr << "Init\n";
    initialize_hamiltonian();
  }
  else
    if (!preadjustment_done&&steps%25==0)
      preadjustment_done = true;
  
  if (!adjustment_done&&preadjustment_done) {
    // fine adjustment
    nob.push_back(get_particle_number());
    if (steps%60==0) {
      double n=0.;
      // averaging of measured number of bosons over last 50 steps
      for (int i=10;i<nob.size();i++)
        n += nob[i];
      n /= (nob.size()-10);
      nob.clear();
      // check if average is close enough to number_of_bosons
      if (fabs(number_of_bosons-n)<number_of_bosons/200.) {
        adjustment_done = true;
        std::cout << "Adjusted parameter is " << static_cast<double>(parms[adjust_parameter]) << ".\n";
        std::cout << "Adjustment took " << steps << " steps.\n";
        if (steps>thermal_sweeps/4.) {
          std::cout << "Adjustment took longer than 25 % of thermalization phase.\n";
          exit(1);
        }
      }
      else {
        if (number_of_bosons>n)
          parms[adjust_parameter] = static_cast<double>(parms[adjust_parameter]) + correction;
        else
          parms[adjust_parameter] = static_cast<double>(parms[adjust_parameter]) - correction;
    std::cerr << "Init\n";
        initialize_hamiltonian();
      }
    }
  }
  std::cerr << "Done adjustment\n";
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

  if (!do_common_measurements(Sign,local)) {
    green = 0.;
    return;
  }

  // measure winding numbers
  bond_iterator bi, bi_end;

  std::vector<int> bond_dir(dimension());
  std::valarray<double> winding_number(0., dimension());
  std::vector<double> boundary_winding(dimension(), 0.);
  for(boost::tie(bi, bi_end) = bonds(); bi != bi_end; ++bi) {
    int s1 = source(*bi);
    int s2 = target(*bi);
    vector_type v=bond_vector_relative(*bi);

    for(int d=0; d<dimension(); ++d)
      bond_dir[d] = boundary_crossing[*bi].crosses(d);
    
    // Distinct worldline events can share a floating time. Only the linked
    // pair with the same nonzero identifier represents a hop on this bond.
    std::map<unsigned, const kink_type*> target_kinks;
    for (const auto& kink : kinks[s2]) if (kink.id()) target_kinks.emplace(kink.id(), &kink);
    auto previous = initial_state(s1);
    for (const auto& kink : kinks[s1]) {
      const auto found = target_kinks.find(kink.id());
      if (kink.id() && found != target_kinks.end()) {
        if (kink.time() != found->second->time())
          throw std::logic_error("linked winding endpoints have different physical times");
        const int hop_dir = int(previous) - int(kink.state());
        for (int d=0; d<dimension(); ++d) {
          winding_number[d] += hop_dir*v[d];
          boundary_winding[d] += hop_dir*bond_dir[d];
        }
      }
      previous = kink.state();
    }
  }

  std::valarray<double> winding_number2(0., dimension());
  winding_number2 = winding_number*winding_number;

  if(use_1D_stiffness){
     // @#$br : record the winding numbers = -1, 0, 1
     //         the vector entries are 0,1,2, respectively            
     std::valarray<double> winding_numberz(3);
     for(int i=0; i<3;++i){  winding_numberz[i]=0. ;}
     if (std::fabs(boundary_winding[0]) <= 1)
       winding_numberz[int(boundary_winding[0])+1] = Sign;
     // Higher winding sectors contribute zero to these three probabilities,
     // but remain in the same sign denominator and W^2 observable.
     measurements["Winding number histogram"] << winding_numberz;
  }
  else{
     // determine stiffness
     double stiffness=0.;
     for(int d=0; d<dimension(); ++d) {
        stiffness += winding_number[d]*winding_number[d];
     }
     measurements["Stiffness"] << stiffness*Sign/(dimension()*beta);
  }

  double vol = num_sites();

  measurements["Statistics"]  << stat; 
  winding_number2 *= Sign;
  measurements["Winding number^2"] << winding_number2;
  measurements["Energy"] << energy*Sign;
  measurements["Energy Density"] << energy/vol*Sign;
  //measurements["Stiffness"] << stiffness*Sign/(dimension()*beta);

  if (measure_green_function_) {
    green /= 2.*double(skip_measurements);
    green *= Sign;
    measurements["Green's Function"] << green;
    green=0.;
  }
  if(chain_kappa) {
    for (int chain=0; chain<num_chains; ++chain) chain_n[chain] /= chain_sizes[chain];
    std::valarray<double> chain_n2 = chain_n*chain_n;

    chain_n  *= Sign;
    chain_n2 *= Sign;

    measurements["Chain density"] << chain_n;
    measurements["Chain density^2"] << chain_n2;
  }

}   // WRun::make_meas


void WRun::measure_green()
{
/*
  int delta=worm_head[1].site()-worm_head[0].site();
  if (delta<0)
    delta +=num_sites();
  for (int h=0;h<(delta==0 ? 2 :1);++h) {
    cyclic_iterator head=worm_head[1-h].kink();
    cyclic_iterator start=head-1;
    cyclic_iterator end=head+1;
    site_descriptor s=worm_head[1-h].site();
    bool do_create = (head->state() < start->state());

    if ((end->time()-start->time())>=(end->time()-worm_head[h]->time())) {
      // needs to be fixed ...
      double statistical_weight= (onsite_energy(start->state(),s)-onsite_energy(head->state(),s));
      double time_interval=end->time() - start->time();
      double time_delta = worm_head[h]->time() - start->time();
      if(time_delta==1.) 
        time_delta=0.;
      double P=exp(-beta*time_delta*statistical_weight)/integrated_weight(statistical_weight, beta*time_interval);
      green[(delta!=0 || do_create) ? delta+1 : delta]+=P;    
    }
  }
  */
}   // WRun::measure_green

