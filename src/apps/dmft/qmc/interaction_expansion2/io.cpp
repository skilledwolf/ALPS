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

void InteractionExpansion::print(std::ostream &os){
  os<<"***********************************************************************************************************"<<std::endl;
  os<<"*** ALPS InteractionExpansion solver                                                                    ***"<<std::endl;
  os<<"*** Emanuel Gull, Philipp Werner, Sebastian Fuchs, Brigitte Surer, Thomas Pruschke, and Matthias Troyer ***"<<std::endl;
  os<<"*** Please cite the interaction-expansion method described in:                                         ***"<<std::endl;
  os<<"*** ***************** Computer Physics Communications 182, 1078 (2011). ******************************* ***"<<std::endl;
  os<<"***********************************************************************************************************"<<std::endl;
  os<<"***                                                                                                     ***"<<std::endl;
  os<<"*** implementing the interaction expansion algorithm by Rubtsov et al., JETP Letters 80, 61.            ***"<<std::endl;
  os<<"***                                                                                                     ***"<<std::endl;
  os<<"***********************************************************************************************************"<<std::endl;
  os<<"max order\t"<<max_order<<"\tn_flavors: "
    <<n_flavors<<"\tn_site: "<<n_site
    <<"\tn_matsubara: "<<n_matsubara<<std::endl;
  os<<"n_tau: "<<n_tau<<"\tmc steps: "<<mc_steps
    <<"\ttherm steps: "<<therm_steps<<std::endl;
  
  os<<"beta: "<<beta<<"\talpha: "<<alpha<<"\tU: "<<onsite_U<<std::endl;
  
  os<<"recalc period: "<<recalc_period<<"\tmeasurement period: "<< measurement_period<<std::endl;
  os<<"almost zero: "<<almost_zero<<std::endl;
}
