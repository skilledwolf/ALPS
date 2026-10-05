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


//A term U n_i n_j of the series expansion. i != j.  E.g. for
//orbitals with G0(c_i, c_j)=0.
double InteractionExpansion::try_add()
{
  assert(n_site==1 && n_flavors >1);
  // Every nonzero unordered pair has proposal probability 1 / pair_count.
  // Row-conditioned proposals would overweight pairs touching sparse rows.
  const auto pair = interaction_pairs.at(interaction_pairs.size() == 1 ? 0 :
      static_cast<std::size_t>(random() * interaction_pairs.size()));
  const auto flavor1 = pair.first, flavor2 = pair.second;

  double t=beta*random();

  double abs_w=beta*U(flavor1, flavor2)*interaction_pairs.size();
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


void InteractionExpansion::perform_add()
{
  //find flavors
  spin_t flavor1=vertices.back().flavor1();
  spin_t flavor2=vertices.back().flavor2();
  //perform the fastupdate up move
  fastupdate_up(flavor1,false);
  fastupdate_up(flavor2,false);
}


void InteractionExpansion::reject_add()
{
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


double InteractionExpansion::try_remove(unsigned int vertex_nr)
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


void InteractionExpansion::perform_remove(unsigned int vertex_nr)
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

