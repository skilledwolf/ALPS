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

#include "interaction_expansion.hpp"


//the following functions work for a single site hubbard model, one site, one flavor.
double HalfFillingHubbardInteractionExpansion::try_add()
{
  assert(n_flavors==1);
  spin_t flavor=0;
  //construct a creator and an annihilator in flavor 0:
  double t = beta*random(); 
  double abs_w = beta*onsite_U*n_site; // onsite_U is const

  unsigned int site = (int)(random()*n_site);
  M[flavor].creators().push_back(creator(up, site, t, n_matsubara));
  M[flavor].annihilators().push_back(annihilator(up, site, t, n_matsubara));
  M[flavor].alpha().push_back(random()<0.5?alpha:1-alpha); //symmetrized version

  //keep track of vertex list
  vertices.push_back(vertex(0, site, M[0].creators().size()-1, M[0].annihilators().size()-1, 0, 0, 0, 0, abs_w)); 
  //second part of list is ignored (symmetry)
  //perform fastupdate up for weight
  double lambda=fastupdate_up(0, true); // true means compute_only_weight
  double metropolis_weight=abs_w/(vertices.size())*lambda*lambda;
  //return weight
  return metropolis_weight;
}


void HalfFillingHubbardInteractionExpansion::perform_add()
{
  //perform the fastupdate up move
  fastupdate_up(0,false);
}


void HalfFillingHubbardInteractionExpansion::reject_add()
{
  //get rid of the operators
  M[0].creators().pop_back();
  M[0].annihilators().pop_back();
  M[0].alpha().pop_back();
  //get rid of the vertex from vertex list
  vertices.pop_back();
}


double HalfFillingHubbardInteractionExpansion::try_remove(unsigned int vertex_nr)
{
  //get weight
  double lambda_1 = fastupdate_down(vertex_nr, 0, true);  // true means compute_only_weight
  double pert_order=num_rows(M[0].matrix());
  //return weight
  return  pert_order/(beta*onsite_U*n_site)*lambda_1*lambda_1;
}


void HalfFillingHubbardInteractionExpansion::perform_remove(unsigned int vertex_nr)
{
  //perform fastupdate down
  fastupdate_down(vertex_nr, 0, false);  // false means really perform, not only compute weight
  //get rid of operators
  M[0].creators().pop_back();
  M[0].annihilators().pop_back();
  M[0].alpha().pop_back();
  //get rid of vertex list entries
  vertices.pop_back();
}


void HalfFillingHubbardInteractionExpansion::reject_remove()
{
  //do nothing
  return;
}


double HubbardInteractionExpansion::try_add()
{
  assert(n_flavors==2);
  double t = beta*random(); 
  double abs_w = beta*onsite_U*n_site; // onsite_U is const
  unsigned int site = (int)(random()*n_site);
  double alpha0 = random()<0.5?alpha:1-alpha;
  double alpha1 = 1 - alpha0;
  spin_t flavor0=0;
  spin_t flavor1=1;
  M[0].creators().push_back(creator(flavor0, site, t, n_matsubara));
  M[0].annihilators().push_back(annihilator(flavor0, site, t, n_matsubara));
  M[0].alpha().push_back(alpha0); //symmetrized version
  M[1].creators().push_back(creator(flavor1, site, t, n_matsubara));
  M[1].annihilators().push_back(annihilator(flavor1, site, t,n_matsubara));
  M[1].alpha().push_back(alpha1); //symmetrized version
  //keep track of vertex list
  vertices.push_back(vertex(flavor0, site, M[0].creators().size()-1, M[0].annihilators().size()-1, 
                        flavor1, site, M[1].creators().size()-1, M[1].annihilators().size()-1, abs_w)); 
  //perform fastupdate up for weight
  double lambda0=fastupdate_up(flavor0, true); // true means compute_only_weight
  double lambda1=fastupdate_up(flavor1, true);
  //std::cout<<"lambda: "<<lambda<<std::endl;
  double metropolis_weight=-abs_w/(vertices.size())*lambda0*lambda1;
  //return weight
  return metropolis_weight;
}


void HubbardInteractionExpansion::perform_add()
{
  //perform the fastupdate up move
  fastupdate_up(0,false);
  fastupdate_up(1,false);
}


void HubbardInteractionExpansion::reject_add()
{
  //get rid of the operators
  M[0].creators().pop_back();
  M[0].annihilators().pop_back();
  M[0].alpha().pop_back();
  M[1].creators().pop_back();
  M[1].annihilators().pop_back();
  M[1].alpha().pop_back();
  //get rid of the vertex from vertex list
  vertices.pop_back();
}


double HubbardInteractionExpansion::try_remove(unsigned int vertex_nr)
{
  //get weight
  double lambda0 = fastupdate_down(vertex_nr, 0, true);  // true means compute_only_weight
  double lambda1 = fastupdate_down(vertex_nr, 1, true);  
  double pert_order=num_rows(M[0].matrix());
  //return weight
  return  -pert_order/(beta*onsite_U*n_site)*lambda0*lambda1;
}


void HubbardInteractionExpansion::perform_remove(unsigned int vertex_nr)
{
  //perform fastupdate down
  fastupdate_down(vertex_nr, 0, false);  // false means really perform, not only compute weight
  fastupdate_down(vertex_nr, 1, false);  // false means really perform, not only compute weight
  //get rid of operators
  M[0].creators().pop_back();
  M[0].annihilators().pop_back();
  M[0].alpha().pop_back();
  M[1].creators().pop_back();
  M[1].annihilators().pop_back();
  M[1].alpha().pop_back();
  //get rid of vertex list entries
  vertices.pop_back();
}


void HubbardInteractionExpansion::reject_remove()
{
  //do nothing
  return;
}


// Adapted from ALPS reconciliation 90a36c0d4c64ee7308dcb0651264197ffdf36449,
// interaction_expansion2/{interaction_expansion.cpp,model.cpp}: propose each
// nonzero unordered density-interaction pair uniformly. The old row-conditioned
// rejection loop hung on isolated flavors and used an incorrect pair probability
// for sparse matrices. Keep the retained model hierarchy and inverse kernels.
MultiBandDensityHubbardInteractionExpansion::MultiBandDensityHubbardInteractionExpansion(
    const alps::params& p,int rank) : InteractionExpansion(p,rank) {
  if (n_site!=1 || n_flavors<2) throw std::invalid_argument("Multiband density interactions require one site and at least two flavors");
  for (spin_t first=0; first<n_flavors; ++first) {
    if (U(first,first)!=0.) throw std::invalid_argument("Multiband density interaction diagonal must vanish");
    for (spin_t second=first+1; second<n_flavors; ++second) {
      const double interaction=U(first,second);
      if (!std::isfinite(interaction) || interaction!=U(second,first))
        throw std::invalid_argument("Multiband density interaction matrix must be finite and symmetric");
      if (interaction!=0.) interaction_pairs_.emplace_back(first,second);
    }
  }
  for (const auto& [first,second] : interaction_pairs_)
    if (!std::isfinite(beta*U(first,second)*interaction_pairs_.size()))
      throw std::invalid_argument("Multiband pair proposal factor is not representable");
}
std::pair<spin_t,spin_t> MultiBandDensityHubbardInteractionExpansion::pair_at(double draw) const {
  if (!std::isfinite(draw) || draw<0. || draw>=1.)
    throw std::domain_error("Multiband pair draw must be in [0,1)");
  return interaction_pairs_.at(static_cast<std::size_t>(draw*interaction_pairs_.size()));
}
//A term U n_i n_j of the series expansion. i != j. E.g. for
//orbitals with G0(c_i,c_j)=0.
double MultiBandDensityHubbardInteractionExpansion::try_add()
{
  assert(n_site==1 && n_flavors >1);
  if (interaction_pairs_.empty()) return 0.; // Exact free model; no vertex appended.
  const auto [flavor1,flavor2]=pair_at(random());
  
  double t=beta*random();

  double abs_w=beta*U(flavor1, flavor2)*interaction_pairs_.size();
  double alpha1 = random()<0.5 ? alpha : 1-alpha;
  double alpha2 = 1-alpha1;
  site_t site1=0;
  site_t site2=0;
  M[flavor1].creators().    push_back(creator(flavor1,site1,t, n_matsubara));
  M[flavor1].annihilators().push_back(annihilator(flavor1,site1,t,n_matsubara));
  M[flavor1].alpha().       push_back(alpha1); //symmetrized version
  M[flavor2].creators().    push_back(creator(flavor2,site2,t, n_matsubara));
  M[flavor2].annihilators().push_back(annihilator(flavor2,site2,t,n_matsubara));
  M[flavor2].alpha().       push_back(alpha2); //symmetrized version
  
  vertices.push_back(vertex(flavor1, site1, M[flavor1].creators().size()-1, M[flavor1].annihilators().size()-1, 
                flavor2, site2, M[flavor2].creators().size()-1, M[flavor2].annihilators().size()-1, abs_w));
  
  double lambda1=fastupdate_up(flavor1, true);
  double lambda2=fastupdate_up(flavor2, true);
  //minus sign because we're not working with down holes but down electrons -> everything picks up a minus sign.
  double metropolis_weight=-abs_w/(vertices.size())*lambda1*lambda2;
  return metropolis_weight;
}


void MultiBandDensityHubbardInteractionExpansion::perform_add()
{    
  //find flavors
  spin_t flavor1=vertices.back().flavor1();
  spin_t flavor2=vertices.back().flavor2();
  //perform the fastupdate up move
  fastupdate_up(flavor1,false);
  fastupdate_up(flavor2,false);
}


void MultiBandDensityHubbardInteractionExpansion::reject_add()
{
  if (interaction_pairs_.empty()) return;
  //find flavors
  spin_t flavor1=vertices.back().flavor1();
  spin_t flavor2=vertices.back().flavor2();
  //get rid of the operators
  M[flavor1].creators().pop_back();
  M[flavor1].annihilators().pop_back();
  M[flavor1].alpha().pop_back();
  M[flavor2].creators().pop_back();
  M[flavor2].annihilators().pop_back();
  M[flavor2].alpha().pop_back();
  //get rid of the vertex from vertex list
  vertices.pop_back();
}


double MultiBandDensityHubbardInteractionExpansion::try_remove(unsigned int vertex_nr)
{
  spin_t flavor1=vertices[vertex_nr].flavor1();
  spin_t flavor2=vertices[vertex_nr].flavor2();
  if(flavor1==flavor2) {throw std::logic_error("bug: flavor1 and flavor2 are equal, we'd require a two vertex removal move for that!");}
  //find operator positions in that flavor
  unsigned int operator_nr_1=vertices[vertex_nr].c_dagger_1();
  unsigned int operator_nr_2=vertices[vertex_nr].c_dagger_2();
  //get weight
  double abs_w=vertices[vertex_nr].abs_w();
  double lambda_1 = fastupdate_down(operator_nr_1, flavor1, true);
  double lambda_2 = fastupdate_down(operator_nr_2, flavor2, true);
  double pert_order=vertices.size();
  
  double metropolis_weight = -pert_order/abs_w*lambda_1*lambda_2;

  return  metropolis_weight; 
}


void MultiBandDensityHubbardInteractionExpansion::perform_remove(unsigned int vertex_nr)
{
  //find flavors
  spin_t flavor1=vertices[vertex_nr].flavor1();
  spin_t flavor2=vertices[vertex_nr].flavor2();
  //find operator positions in that flavor
  unsigned int operator_nr_1=vertices[vertex_nr].c_dagger_1();
  unsigned int operator_nr_2=vertices[vertex_nr].c_dagger_2();
  //perform fastupdate down
  fastupdate_down(operator_nr_1, flavor1, false);
  fastupdate_down(operator_nr_2, flavor2, false);
  //take care of vertex list
  for(int i=vertices.size()-1;i>=0;--i){
    //this operator pointed to the last row/column of M[flavor1], which has just been moved to vertex_nr.
    if(vertices[i].flavor1()==flavor1 && vertices[i].c_dagger_1()==num_rows(M[flavor1].matrix())){ 
      vertices[i].c_dagger_1()=operator_nr_1;
      vertices[i].c_1()=operator_nr_1;
      break;
    }
    if(vertices[i].flavor2()==flavor1 && vertices[i].c_dagger_2()==num_rows(M[flavor1].matrix())){
      vertices[i].c_dagger_2()=operator_nr_1;
      vertices[i].c_2()=operator_nr_1;
      break;
    }
  }
  for(int i=vertices.size()-1;i>=0;--i){
    if(vertices[i].flavor1()==flavor2 && vertices[i].c_dagger_1()==num_rows(M[flavor2].matrix())){
      vertices[i].c_dagger_1()=operator_nr_2;
      vertices[i].c_1()=operator_nr_2;
      break;
    }
    if(vertices[i].flavor2()==flavor2 && vertices[i].c_dagger_2()==num_rows(M[flavor2].matrix())){
      vertices[i].c_dagger_2()=operator_nr_2;
      vertices[i].c_2()=operator_nr_2;
      break;
    }
  }
  vertices[vertex_nr]=vertices.back();

  //get rid of operators
  M[flavor1].creators().pop_back();
  M[flavor1].annihilators().pop_back();
  M[flavor1].alpha().pop_back();
  M[flavor2].creators().pop_back();
  M[flavor2].annihilators().pop_back();
  M[flavor2].alpha().pop_back();
  //get rid of vertex list entries
  vertices.pop_back();

}


void MultiBandDensityHubbardInteractionExpansion::reject_remove()
{
  //do nothing
  return;
}


