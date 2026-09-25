/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#pragma once

#include "types.hpp"
#include <boost/numeric/bindings/blas/level2/ger.hpp>

namespace alps { namespace ctint {

// The numerical workflow is shared by both scheduler adapters. All dispatch
// is resolved at compile time; RNG calls retain each runtime's ordering.
template<class X, class Y> inline Y linear_interpolate(const X x0, const X x1, const Y y0, const Y y1, const X x)
{
  return y0 + (x-x0)/(x1-x0)*(y1-y0);
}

struct update_kernel {
  template<class Run>
  static double fastupdate_up(Run& run, const int flavor, bool compute_only_weight) {
  assert(num_rows(run.M[flavor].matrix()) == num_cols(run.M[flavor].matrix()));
  unsigned int noperators = num_rows(run.M[flavor].matrix());
  //current size of M: number of vertices -1. We need to add the last vertex.
  //A pointer to the creator and annihilator is already stored in creators_ and annihilators_ at position M[flavor].size();
  double Green0_n_n=run.green0_spline(run.M[flavor].creators()[noperators],run.M[flavor].annihilators()[noperators]);
  alps::numeric::vector<double> Green0_n_j(noperators);
  alps::numeric::vector<double> Green0_j_n(noperators);
  alps::numeric::vector<double> lastcolumn(noperators);
  alps::numeric::vector<double> lastrow(noperators);
  //compute the Green's functions with interpolation
  for(unsigned int i=0;i<noperators;++i){
    Green0_n_j[i]=run.green0_spline(run.M[flavor].creators()[noperators],run.M[flavor].annihilators()[i]);
    Green0_j_n[i]=run.green0_spline(run.M[flavor].creators()[i],run.M[flavor].annihilators()[noperators]);
  }
  //compute the last row
  if(noperators!=0)
    gemv(run.M[flavor].matrix(),Green0_j_n,lastcolumn);
  //compute lambda
  double ip = (noperators==0?0:scalar_product(Green0_n_j, lastcolumn));
  double lambda = Green0_n_n - ip + run.M[flavor].alpha()[noperators];
  //return weight if we have nothing else to do
  if(compute_only_weight){
    return lambda;
  }
  //compute last column
  if(noperators!=0)
    gemv(transpose(run.M[flavor].matrix()), Green0_n_j, lastrow);
  //compute norm of vector and mean for roundoff error check
  if(noperators > 0){
    lastrow    *= 1./lambda;
    boost::numeric::bindings::blas::ger(1.0,lastcolumn,lastrow,run.M[flavor].matrix());
    lastcolumn *= 1./lambda;
    //std::cout<<lambda<<" "<<M[flavor]<<std::endl;
  }
  //add row and column to M
  resize(run.M[flavor].matrix(), noperators+1, noperators+1);
  for(unsigned int i=0;i<noperators;++i){
    run.M[flavor].matrix()(i, noperators)=-lastcolumn[i];
  }
  for(unsigned int i=0;i<noperators;++i){
    run.M[flavor].matrix()(noperators, i)=-lastrow[i];
  }
  run.M[flavor].matrix()(noperators, noperators)=1./lambda;
  return lambda;
}

  template<class Run>
  static double fastupdate_down(Run& run, const int operator_nr, const int flavor, bool compute_only_weight) {
  using std::swap;
  using alps::numeric::column_view;
  using alps::numeric::vector;
  using alps::numeric::matrix;
  assert(num_rows(run.M[flavor].matrix()) == num_cols(run.M[flavor].matrix()));
  //perform updates according to formula 21.1, 21.2
  unsigned int noperators=num_rows(run.M[flavor].matrix());  //how many operators do we have in total?
  if(compute_only_weight){
    return run.M[flavor].matrix()(operator_nr,operator_nr);
  }
  //swap rows and colums of M <-> move selected vertex to the end.
  for(unsigned int i=0;i<noperators;++i){
    swap(run.M[flavor].matrix()(i,noperators-1), run.M[flavor].matrix()(i,operator_nr));
  }
  for(unsigned int i=0;i<noperators;++i){
    swap(run.M[flavor].matrix()(noperators-1,i),run.M[flavor].matrix()(operator_nr,i));
  }
  //swap creator and annihilator
  swap(run.M[flavor].creators()[operator_nr],     run.M[flavor].creators()[noperators-1]);
  swap(run.M[flavor].annihilators()[operator_nr], run.M[flavor].annihilators()[noperators-1]);
  swap(run.M[flavor].alpha()[operator_nr],        run.M[flavor].alpha()[noperators-1]);
  double Mnn=run.M[flavor].matrix()(noperators-1,noperators-1);
  //now perform fastupdate of M
  vector<double> lastrow(noperators-1);
  vector<double> lastcolumn(column_view<matrix<double> >(run.M[flavor].matrix(),noperators-1));
  for(unsigned int j=0;j<noperators-1;++j){
    lastrow[j]=run.M[flavor].matrix()(noperators-1, j);
  }
  if(noperators>1)
    lastcolumn *= -1./Mnn;
  resize(run.M[flavor].matrix(),noperators-1,noperators-1);  //lose the last row and last column, reduce size by one, but keep contents.
  if(noperators>1)
    boost::numeric::bindings::blas::ger(1.0,lastcolumn,lastrow,run.M[flavor].matrix());
  return Mnn;  //the determinant ratio det D_{k-1}/det D_{k}
}

  template<class Run>
  static double green0_spline(const Run& run, const creator &cdagger, const annihilator &c) {
  //for the M matrix we need the bare temporal Green's function.
  //For the dressed frequency GF we need the bare frequency GF.
  //here we get the bare temporal GF, interpolated, as needed for the M matrix.
  //we will receive this as an input into our solver later.
  spin_t flavor;
  if((flavor=cdagger.flavor()) != c.flavor()){
    return 0; //the delta in spin space.
  }
  itime_t delta_t=cdagger.t()-c.t();
  site_t site1 = cdagger.s();
  site_t site2 = c.s();
  return run.green0_spline(delta_t, flavor, site1, site2);
}

  template<class Run>
  static double green0_spline(const Run& run, const itime_t delta_t, const spin_t flavor, const site_t site1, const site_t site2) {
  if(delta_t*delta_t < run.almost_zero){
    return run.bare_green_itime(0,site1, site2, flavor);
  }
  else if(delta_t>0){
    int time_index_1 = (int)(delta_t*run.n_tau*run.temperature);
    int time_index_2 = time_index_1+1;
    return linear_interpolate((double)time_index_1*run.beta*run.n_tau_inv, (double)time_index_2*run.beta*run.n_tau_inv,
                              run.bare_green_itime(time_index_1,site1, site2, flavor),
                              run.bare_green_itime(time_index_2,site1, site2, flavor),delta_t);
  }
  else{
    int time_index_1 = (int)(delta_t*run.n_tau*run.temperature+run.n_tau);
    int time_index_2 = time_index_1+1;
    return -linear_interpolate((double)time_index_1*run.beta*run.n_tau_inv, (double)time_index_2*run.beta*run.n_tau_inv,
                               run.bare_green_itime(time_index_1,site1,site2,flavor),
                               run.bare_green_itime(time_index_2,site1,site2,flavor),delta_t+run.beta);
  }
}

  template<class Run>
  static double green0_spline(const Run& run, const itime_t delta_t, const spin_t flavor) {
  if(delta_t*delta_t < run.almost_zero){
    return run.bare_green_itime(0, flavor);
  }
  else if(delta_t>0){
    int time_index_1 = (int)(delta_t*run.n_tau*run.temperature);
    int time_index_2 = time_index_1+1;
    return linear_interpolate((double)time_index_1*run.beta*run.n_tau_inv, (double)time_index_2*run.beta*run.n_tau_inv,
                              run.bare_green_itime(time_index_1,flavor),
                              run.bare_green_itime(time_index_2,flavor),delta_t);
  }
  else{
    int time_index_1 = (int)(delta_t*run.n_tau*run.temperature+run.n_tau);
    int time_index_2 = time_index_1+1;
    return -linear_interpolate((double)time_index_1*run.beta*run.n_tau_inv, (double)time_index_2*run.beta*run.n_tau_inv, //time ordering
                               run.bare_green_itime(time_index_1,flavor),
                               run.bare_green_itime(time_index_2,flavor),delta_t+run.beta);
  }
}

  template<class Run>
  static void interaction_expansion_step(Run& run) {
  int pert_order=run.vertices.size();   //current order of perturbation series
  double metropolis_weight=0.;
  //double oldsign=sign;
  if(run.ctint_uniform()<0.5){  //trying to ADD vertex
    if(run.vertices.size()>=run.max_order)
      return; //we have already reached the highest perturbation order
    metropolis_weight=run.try_add();
    if(std::fabs(metropolis_weight)> run.ctint_uniform()){
      run.measurements["VertexInsertion"]<<1.;
      run.perform_add();
      run.sign*=metropolis_weight<0?-1:1;
    }else{
      run.measurements["VertexInsertion"]<<0.;
      run.reject_add();
    }
  }else{ // try to REMOVE a vertex
    pert_order=run.vertices.size(); //choose a vertex
    if(pert_order < 1)
      return;     //we have an empty list or one with just one vertex
    //this might be the ideal place to do some cleanup, e.g. get rid of the roundoff errors and such.
    int vertex_nr=(int)(run.ctint_uniform() * pert_order);
    metropolis_weight=run.try_remove(vertex_nr); //get the determinant ratio. don't perform fastupdate yet
    if(std::fabs(metropolis_weight)> run.ctint_uniform()){ //do the actual update
      run.measurements["VertexRemoval"]<<1.;
      run.perform_remove(vertex_nr);
      run.sign*=metropolis_weight<0?-1:1;
    }else{
      run.measurements["VertexRemoval"]<<0.;
      run.reject_remove();
    }
  }//end REMOVE
  run.weight=metropolis_weight;
}

  template<class Run>
  static void reset_perturbation_series(Run& run) {
  std::vector<inverse_m_matrix> M2(run.M); //make a copy of M
  vertex_array vertices_backup;
  for(unsigned int i=0;i<run.vertices.size();++i){
    vertices_backup.push_back(run.vertices[i]);
  }
  run.vertices.clear();
  run.sign=1;
  for(spin_t z=0;z<run.n_flavors;++z){
    resize(run.M[z].matrix(),0,0);
  }
  run.green_matsubara = run.bare_green_matsubara;
  run.green_itime     = run.bare_green_itime;
  //recompute M from scratch
  for(unsigned int i=0;i<vertices_backup.size();++i){
    run.vertices.push_back(vertices_backup[i]);
    run.perform_add();
  }
  for(unsigned int z=0;z<M2.size();++z){
    double max_diff=0;
    for(unsigned int j=0;j<num_cols(M2[z].matrix());++j){
      for(unsigned int i=0;i<num_rows(M2[z].matrix());++i){
        double diff=run.M[z].matrix()(i,j)-M2[z].matrix()(i,j);
        if(std::abs(diff)>max_diff) max_diff=std::abs(diff);
      }
    }
    if(max_diff > 1.e-8)
      std::cout<<"WARNING: roundoff errors in flavor: "<<z<<" max diff "<<max_diff<<std::endl;
  }
}

  template<class Run>
  static double half_filling_try_add(Run& run) {
  assert(run.n_flavors==1);
  spin_t flavor=0;
  //construct a creator and an annihilator in flavor 0:
  double t = run.beta*run.ctint_uniform();
  double abs_w = run.beta*run.onsite_U*run.n_site; // onsite_U is const

  unsigned int site = run.ctint_site(run.n_site);
  run.M[flavor].creators().push_back(creator(up, site, t, run.n_matsubara));
  run.M[flavor].annihilators().push_back(annihilator(up, site, t, run.n_matsubara));
  run.M[flavor].alpha().push_back(run.ctint_uniform()<0.5?run.alpha:1-run.alpha); //symmetrized version

  //keep track of vertex list
  run.vertices.push_back(vertex(0, site, run.M[0].creators().size()-1, run.M[0].annihilators().size()-1, 0, 0, 0, 0, abs_w));
  //second part of list is ignored (symmetry)
  //perform fastupdate up for weight
  double lambda=run.fastupdate_up(0, true); // true means compute_only_weight
  double metropolis_weight=abs_w/(run.vertices.size())*lambda*lambda;
  //return weight
  return metropolis_weight;
}

  template<class Run>
  static void half_filling_perform_add(Run& run) {
  //perform the fastupdate up move
  run.fastupdate_up(0,false);
}

  template<class Run>
  static void half_filling_reject_add(Run& run) {
  //get rid of the operators
  run.M[0].creators().pop_back();
  run.M[0].annihilators().pop_back();
  run.M[0].alpha().pop_back();
  //get rid of the vertex from vertex list
  run.vertices.pop_back();
}

  template<class Run>
  static double half_filling_try_remove(Run& run, unsigned int vertex_nr) {
  assert(num_rows(run.M[0].matrix()) == num_cols(run.M[0].matrix()));
  //get weight
  double lambda_1 = run.fastupdate_down(vertex_nr, 0, true);  // true means compute_only_weight
  double pert_order=num_rows(run.M[0].matrix());
  //return weight
  return  pert_order/(run.beta*run.onsite_U*run.n_site)*lambda_1*lambda_1;
}

  template<class Run>
  static void half_filling_perform_remove(Run& run, unsigned int vertex_nr) {
  //perform fastupdate down
  run.fastupdate_down(vertex_nr, 0, false);  // false means really perform, not only compute weight
  //get rid of operators
  run.M[0].creators().pop_back();
  run.M[0].annihilators().pop_back();
  run.M[0].alpha().pop_back();
  //get rid of vertex list entries
  run.vertices.pop_back();
}

  template<class Run>
  static void half_filling_reject_remove(Run& run) {
  //do nothing
  return;
}

  template<class Run>
  static double hubbard_try_add(Run& run) {
  assert(run.n_flavors==2);
  double t = run.beta*run.ctint_uniform();
  double abs_w = run.beta*run.onsite_U*run.n_site; // onsite_U is const
  unsigned int site = run.ctint_site(run.n_site);
  double alpha0 = run.ctint_uniform()<0.5?run.alpha:1-run.alpha;
  double alpha1 = 1 - alpha0;
  spin_t flavor0=0;
  spin_t flavor1=1;
  run.M[0].creators().push_back(creator(flavor0, site, t, run.n_matsubara));
  run.M[0].annihilators().push_back(annihilator(flavor0, site, t, run.n_matsubara));
  run.M[0].alpha().push_back(alpha0); //symmetrized version
  run.M[1].creators().push_back(creator(flavor1, site, t, run.n_matsubara));
  run.M[1].annihilators().push_back(annihilator(flavor1, site, t,run.n_matsubara));
  run.M[1].alpha().push_back(alpha1); //symmetrized version
  //keep track of vertex list
  run.vertices.push_back(vertex(flavor0, site, run.M[0].creators().size()-1, run.M[0].annihilators().size()-1,
                        flavor1, site, run.M[1].creators().size()-1, run.M[1].annihilators().size()-1, abs_w));
  //perform fastupdate up for weight
  double lambda0=run.fastupdate_up(flavor0, true); // true means compute_only_weight
  double lambda1=run.fastupdate_up(flavor1, true);
  //std::cout<<"lambda: "<<lambda<<std::endl;
  double metropolis_weight=-abs_w/(run.vertices.size())*lambda0*lambda1;
  //return weight
  return metropolis_weight;
}

  template<class Run>
  static void hubbard_perform_add(Run& run) {
  //perform the fastupdate up move
  run.fastupdate_up(0,false);
  run.fastupdate_up(1,false);
}

  template<class Run>
  static void hubbard_reject_add(Run& run) {
  //get rid of the operators
  run.M[0].creators().pop_back();
  run.M[0].annihilators().pop_back();
  run.M[0].alpha().pop_back();
  run.M[1].creators().pop_back();
  run.M[1].annihilators().pop_back();
  run.M[1].alpha().pop_back();
  //get rid of the vertex from vertex list
  run.vertices.pop_back();
}

  template<class Run>
  static double hubbard_try_remove(Run& run, unsigned int vertex_nr) {
  assert(num_rows(run.M[0].matrix()) == num_cols(run.M[0].matrix()));
  //get weight
  double lambda0 = run.fastupdate_down(vertex_nr, 0, true);  // true means compute_only_weight
  double lambda1 = run.fastupdate_down(vertex_nr, 1, true);
  double pert_order=num_rows(run.M[0].matrix());
  //return weight
  return  -pert_order/(run.beta*run.onsite_U*run.n_site)*lambda0*lambda1;
}

  template<class Run>
  static void hubbard_perform_remove(Run& run, unsigned int vertex_nr) {
  //perform fastupdate down
  run.fastupdate_down(vertex_nr, 0, false);  // false means really perform, not only compute weight
  run.fastupdate_down(vertex_nr, 1, false);  // false means really perform, not only compute weight
  //get rid of operators
  run.M[0].creators().pop_back();
  run.M[0].annihilators().pop_back();
  run.M[0].alpha().pop_back();
  run.M[1].creators().pop_back();
  run.M[1].annihilators().pop_back();
  run.M[1].alpha().pop_back();
  //get rid of vertex list entries
  run.vertices.pop_back();
}

  template<class Run>
  static void hubbard_reject_remove(Run& run) {
  //do nothing
  return;
}

  template<class Run>
  static double multiband_try_add(Run& run) {
  assert(run.n_site==1 && run.n_flavors >1);
  spin_t flavor1=(spin_t)(run.ctint_uniform()*run.n_flavors);
  spin_t flavor2;

  do{ flavor2=(spin_t)(run.ctint_uniform()*run.n_flavors);
  } while (run.U(flavor1, flavor2)==0);

  double t=run.beta*run.ctint_uniform();

  double abs_w=run.beta*run.U(flavor1, flavor2)/2;
  double alpha1 = run.ctint_uniform()<0.5 ? run.alpha : 1-run.alpha;
  double alpha2 = 1-alpha1;
  site_t site1=0;
  site_t site2=0;
  run.M[flavor1].creators().    push_back(creator(flavor1,site1,t, run.n_matsubara));
  run.M[flavor1].annihilators().push_back(annihilator(flavor1,site1,t,run.n_matsubara));
  run.M[flavor1].alpha().       push_back(alpha1); //symmetrized version
  run.M[flavor2].creators().    push_back(creator(flavor2,site2,t, run.n_matsubara));
  run.M[flavor2].annihilators().push_back(annihilator(flavor2,site2,t,run.n_matsubara));
  run.M[flavor2].alpha().       push_back(alpha2); //symmetrized version

  run.vertices.push_back(vertex(flavor1, site1, run.M[flavor1].creators().size()-1, run.M[flavor1].annihilators().size()-1,
                flavor2, site2, run.M[flavor2].creators().size()-1, run.M[flavor2].annihilators().size()-1, abs_w));

  double lambda1=run.fastupdate_up(flavor1, true);
  double lambda2=run.fastupdate_up(flavor2, true);
  double sym_factor = run.U.n_nonzero();
  //minus sign because we're not working with down holes but down electrons -> everything picks up a minus sign.
  double metropolis_weight=-abs_w*sym_factor/(run.vertices.size())*lambda1*lambda2;
  return metropolis_weight;
}

  template<class Run>
  static void multiband_perform_add(Run& run) {
  //find flavors
  spin_t flavor1=run.vertices.back().flavor1();
  spin_t flavor2=run.vertices.back().flavor2();
  //perform the fastupdate up move
  run.fastupdate_up(flavor1,false);
  run.fastupdate_up(flavor2,false);
}

  template<class Run>
  static void multiband_reject_add(Run& run) {
  //find flavors
  spin_t flavor1=run.vertices.back().flavor1();
  spin_t flavor2=run.vertices.back().flavor2();
  //get rid of the operators
  run.M[flavor1].creators().pop_back();
  run.M[flavor1].annihilators().pop_back();
  run.M[flavor1].alpha().pop_back();
  run.M[flavor2].creators().pop_back();
  run.M[flavor2].annihilators().pop_back();
  run.M[flavor2].alpha().pop_back();
  //get rid of the vertex from vertex list
  run.vertices.pop_back();
}

  template<class Run>
  static double multiband_try_remove(Run& run, unsigned int vertex_nr) {
  spin_t flavor1=run.vertices[vertex_nr].flavor1();
  spin_t flavor2=run.vertices[vertex_nr].flavor2();
  if(flavor1==flavor2) {
    run.ctint_invalid_flavors("bug: flavor1 and flavor2 are equal, we'd require a two vertex removal move for that!");
  }
  //find operator positions in that flavor
  unsigned int operator_nr_1=run.vertices[vertex_nr].c_dagger_1();
  unsigned int operator_nr_2=run.vertices[vertex_nr].c_dagger_2();
  //get weight
  double abs_w=run.vertices[vertex_nr].abs_w();
  double lambda_1 = run.fastupdate_down(operator_nr_1, flavor1, true);
  double lambda_2 = run.fastupdate_down(operator_nr_2, flavor2, true);
  double pert_order=run.vertices.size();

  double sym_factor = run.U.n_nonzero();
  double metropolis_weight = -pert_order/abs_w/sym_factor*lambda_1*lambda_2;

  return  metropolis_weight;
}

  template<class Run>
  static void multiband_perform_remove(Run& run, unsigned int vertex_nr) {
  //find flavors
  spin_t flavor1=run.vertices[vertex_nr].flavor1();
  spin_t flavor2=run.vertices[vertex_nr].flavor2();
  //find operator positions in that flavor
  unsigned int operator_nr_1=run.vertices[vertex_nr].c_dagger_1();
  unsigned int operator_nr_2=run.vertices[vertex_nr].c_dagger_2();
  //perform fastupdate down
  run.fastupdate_down(operator_nr_1, flavor1, false);
  run.fastupdate_down(operator_nr_2, flavor2, false);
  assert(num_rows(run.M[flavor1].matrix()) == num_cols(run.M[flavor1].matrix()));
  //take care of vertex list
  for(int i=run.vertices.size()-1;i>=0;--i){
    //this operator pointed to the last row/column of M[flavor1], which has just been moved to vertex_nr.
    if(run.vertices[i].flavor1()==flavor1 && run.vertices[i].c_dagger_1()==num_rows(run.M[flavor1].matrix())){
      run.vertices[i].c_dagger_1()=operator_nr_1;
      run.vertices[i].c_1()=operator_nr_1;
      break;
    }
    if(run.vertices[i].flavor2()==flavor1 && run.vertices[i].c_dagger_2()==num_rows(run.M[flavor1].matrix())){
      run.vertices[i].c_dagger_2()=operator_nr_1;
      run.vertices[i].c_2()=operator_nr_1;
      break;
    }
  }
  for(int i=run.vertices.size()-1;i>=0;--i){
    if(run.vertices[i].flavor1()==flavor2 && run.vertices[i].c_dagger_1()==num_rows(run.M[flavor2].matrix())){
      run.vertices[i].c_dagger_1()=operator_nr_2;
      run.vertices[i].c_1()=operator_nr_2;
      break;
    }
    if(run.vertices[i].flavor2()==flavor2 && run.vertices[i].c_dagger_2()==num_rows(run.M[flavor2].matrix())){
      run.vertices[i].c_dagger_2()=operator_nr_2;
      run.vertices[i].c_2()=operator_nr_2;
      break;
    }
  }
  run.vertices[vertex_nr]=run.vertices.back();

  //get rid of operators
  run.M[flavor1].creators().pop_back();
  run.M[flavor1].annihilators().pop_back();
  run.M[flavor1].alpha().pop_back();
  run.M[flavor2].creators().pop_back();
  run.M[flavor2].annihilators().pop_back();
  run.M[flavor2].alpha().pop_back();
  //get rid of vertex list entries
  run.vertices.pop_back();

}

  template<class Run>
  static void multiband_reject_remove(Run& run) {
  //do nothing
  return;
}

};

} } // namespace alps::ctint
