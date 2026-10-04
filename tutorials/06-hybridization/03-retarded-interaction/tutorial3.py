 #############################################################################/
 #
 # ALPS Project: Algorithms and Libraries for Physics Simulations
 #
 # ALPS Libraries
 #
 # Copyright (C) 2012 by Hartmut Hafermann <hafermann@cpht.polytechnique.fr>
 #
 #
 # ALPS Project: https://alps.comp-phys.org/
 # SPDX-License-Identifier: MIT
 #
 #############################################################################/

 # This tutorial implements the DMFT self-consistency on the Bethe-lattice
 # for the Hubbard-Holstein model with a single phonon mode at w0 and coupling
 # lambda. The latter gives rise to a retarded interaction in the impurity model.
 # The script can be used to reproduce the results of Fig.3 of the ALPS CT-HYB paper.
 #
 # The impurity model is solved using the cthyb module (python interface to the
 # hybridization solver). For details on how to use the solver see its documentation.
 #
 # The tutorial shows how to flexibly implement a selfconsistency loop without using
 # the ALPS DMFT framework. The selfconsistency is based on G(tau) and additional
 # quantities (such as self-energy etc.) are measured in the last iteration only.
 # These are written to stdout in human-readable (plottable) format.
 #
 # Run this script as:
 # alpspython tutorial2.py
 #
 # This python script is MPI aware and can hence be called using mpirun:
 #
 # mpirun -np 2 alpspython tutorial2.py
 #
 # In case this does not work, try:
 #
 # mpirun -np 2 sh alpspython tutorial2.py

import shutil
import pyalps.mpi as mpi                # mpi library
from pyalps.hdf5 import archive            # hdf5 interface
from pyalps.ngs import params
import pyalps.cthyb as cthyb            # the solver module
from numpy import sqrt,cosh,sinh,exp,pi #some math
from numpy import array,zeros,append
##################################################################################################################
#                                                                                                                #
#                                               P A R A M E T E R S                                              #
#                                                                                                                #
##################################################################################################################

# This script takes some time to run to get converged results. Make sure the runtime is large enough to get
# sensible results (depends on the number of processes).
runtime_dmft=60        # runtime for each DMFT iteration
runtime_dmft_final=600 # increase runtime on final iteration for additional measurements
dmft_iterations=10     # number of DMFT iterations

# perform calculations for fixed U and Uscr and given w0
U      = 8.0
Uscr   = 3.0
w0     = 12.0
Lambda = sqrt((U-Uscr)*w0/2.0) # choose lambda as to get the given Uscr = U - 2*lambda^2/w0

# for simplicity, this script is for a single set of parameters only
# we list all solver parameters here for completeness
parms = {
# general
'SWEEPS'                     : 1000000000,                         #sweeps to be done
'THERMALIZATION'             : 1000,                               #thermalization sweeps to be done
'N_MEAS'                     : 100,                                #number of sweeps after which a measurement is done
'N_ORBITALS'                 : 2,                                  #number of 'orbitals', i.e. number of spin-orbital degrees of freedom or segments
'VERBOSE'                    : True,                               #whether to output extra information
'COMPUTE_VERTEX'             : False,                              #whether to compute the vertex function
# physics parameters
'U'                          : U,                                  #Hubbard repulsion
'MU'                         : U/2.-2*Lambda**2/w0,                #chemical potential (MU=U/2-2K'(0) corresponds to half-filling; here K(0)=Lambda^2/w0)
'BETA'                       : 50.0,                               #inverse temperature
# measurements
'MEASURE_freq'               : False,                              #whether to measure single-particle Green's function on Matsubara frequencies
'MEASURE_legendre'           : False,                              #whether to measure single-particle Green's function in Legendre polynomial basis
'MEASURE_g2w'                : False,                              #whether to measure two-particle Green's function on Matsubara frequencies
'MEASURE_h2w'                : False,                              #whether to measure the higher-order correlation function for the vertex on Matsubara frequencies
'MEASURE_nn'                 : False,                              #whether to measure equal-time density-density correlations
'MEASURE_nnt'                : False,                              #whether to measure the density-density correlation function (local susceptibility) in imaginary time
'MEASURE_nnw'                : False,                              #whether to measure the density-density correlation function (local susceptibility) on Matsubara frequencies
'MEASURE_sector_statistics'  : False,                              #whether to measure sector statistics
# measurement parameters
'N_HISTOGRAM_ORDERS'         : 50,                                 #maximum order for the perturbation order histogram
'N_TAU'                      : 5000,                               #number of imaginary time points (tau_0=0, tau_N_TAU=BETA)
'N_MATSUBARA'                : 512,                                #number of Matsubara frequencies
'N_nn'                       : 5000,                               #number of imaginary time points for the density-density correlation function
'N_W'                        : 20,                                 #number of bosonic Matsubara frequencies for the two-particle Green's function or local susceptibility
'N_w2'                       : 20,                                 #number of fermionic Matsubara frequencies for the two-particle Green's function
'N_LEGENDRE'                 : 80,                                 #number of Legendre coefficients
}# parms

# solver input files: hybridization function and retarded interaction K(tau), K'(tau)
inputs = {'delta': "Delta.h5", 'delta_format': "hdf5",
          'retarded_interaction': "K_tau.h5", 'retarded_interaction_format': "hdf5"}
results = "hyb.param.out.h5" # name of the h5 output file
runtime = runtime_dmft       # runtime of the solver per iteration
text_output = False          # whether to write results in human readable (text) format

# additional parameters (used outside the solver only)
hopping = 1.   # hopping
mix     = 0.5  # mixing parameter for hybridization update

if mpi.rank==0:
  print("generating initial hybridization...")
  g=[]
  I=complex(0.,1.)
  mu=0.0
  for n in range(parms['N_MATSUBARA']):
    w=(2*n+1)*pi/parms['BETA']
    g.append(2.0/(I*w+mu+I*sqrt(4*hopping**2-(I*w+mu)**2))) # noninteracting Green's function on Bethe lattice
  delta=[]
  for i in range(parms['N_TAU']+1):
    tau=i*parms['BETA']/parms['N_TAU']
    g0tau=0.0;
    for n in range(parms['N_MATSUBARA']):
      iw=complex(0.0,(2*n+1)*pi/parms['BETA'])
      g0tau+=((g[n]-1.0/iw)*exp(-iw*tau)).real # Fourier transform with tail subtracted
    g0tau *= 2.0/parms['BETA']
    g0tau += -1.0/2.0 # add back contribution of the tail
    delta.append(hopping**2*g0tau) # delta=t**2 g

  # write hybridization function to hdf5 archive (solver input)
  ar=archive(inputs['delta'],'w')
  for m in range(parms['N_ORBITALS']):
    ar['/Delta_%i'%m]=delta
  del ar

  print("generating retarded interaction...")
  l   =Lambda
  beta=parms['BETA']

  K  = lambda tau: - (l**2)*(cosh(w0*(beta/2.0-tau))/sinh(w0*beta/2.0) - cosh(w0*beta/2.0)/sinh(w0*beta/2.0) )/(w0*w0)
  Kp = lambda tau: + (l**2)*(sinh(w0*(beta/2.0-tau))/sinh(w0*beta/2.0))/w0

  print("U_screened =", parms['U'] - 2*l**2/w0)

  k_tau=[]
  kp_tau=[]
  for i in range(parms['N_TAU']+1):
    tau=i*parms['BETA']/parms['N_TAU']
    k_tau.append(K(tau))
    kp_tau.append(Kp(tau))

  # write retarded interaction function K(tau) and its derivative to file (solver input)
  ar=archive(inputs['retarded_interaction'],'w')
  ar['/Ret_int_K']=k_tau
  ar['/Ret_int_Kp']=kp_tau
  del ar

  if text_output:
    f=open('Ktau.dat','w')
    for i in range(len(k_tau)):
      tau=i*parms['BETA']/parms['N_TAU']
      f.write("%f %f %f\n"%(tau,k_tau[i],kp_tau[i]))

mpi.world.barrier() # wait until solver input is written to file

###################################################################################################################
#                                                                                                                 #
#                                D M F T   S E L F C O N S I S T E N C Y    L O O P                               #
#                                                                                                                 #
###################################################################################################################
for it in range(dmft_iterations):

  if mpi.rank==0:
    print("****************************************************************************")
    print("*                           DMFT iteration %3i                             *"%(it))
    print("****************************************************************************")

  # !always make sure that parameters are changed on all threads equally!
  # (i.e. don't wrap this into an 'if mpi.rank==0' statement)
  if it==dmft_iterations-1:
    runtime = runtime_dmft_final
    # turn on additional measurements for the final dmft interation
    parms['MEASURE_freq']=True # turn of Matsubara measurement
    parms['MEASURE_legendre']=True # turn on Legendre measurement
    text_output=True  # this will write results of the final iteration in text format

  # write parameters for reference (on master only)
  if mpi.rank==0:
    ar=archive('hyb.param.h5','a')
    ar['/parameters']=params(parms)
    ar['/parameters%i'%it]=params(parms) # this is a backup for each iteration
    del ar

  # solve the impurity model in parallel
  cthyb.solve(cthyb.prepare(parms, input=inputs, output={'results': results, 'text': text_output},
                            execution={'time_limit': runtime, 'seed': 42}))

  # self-consistency on the master
  if mpi.rank==0:
    if text_output:
      shutil.copy("Gt.dat", "Gt%i.dat"%it) # keep Green's function for each iteration for monitoring
      shutil.copy("simulation.dat", "simulation%i.dat"%it) # keep some basic information for each iteration

    # read Green's function from file
    ar=archive(results,'r')
    # symmetrize G(tau)
    # in this case all orbitals and spins are degenerate
    gt=array(zeros(parms['N_TAU']+1))
    for m in range(parms['N_ORBITALS']):
      gt+=ar['G_tau/%i/mean/value'%m]
    gt/=parms['N_ORBITALS']
    del ar

    # Bethe lattice self-consistency: delta(tau)=t**2 g(tau)
    # read delta_old
    ar=archive(inputs['delta'],'rw')
    for m in range(parms['N_ORBITALS']):
      delta_old=ar['/Delta_%i'%m]
      delta_new=array(zeros(parms['N_TAU']+1))
      delta_new=(1.-mix)*(hopping**2 * gt) + mix*delta_old # mix old and new delta

      # write hybridization to the h5 archive (this is solver input)
      ar['/Delta_%i'%m]=delta_new
    del ar

    # write input hybridization function for reference
    if text_output:
      f=open('Delta%i.dat'%it,'w')
      for i in range(parms['N_TAU']+1):
        f.write("%f"%(i*parms['BETA']/parms['N_TAU']))
        for m in range(parms['N_ORBITALS']):
          f.write(" %f"%delta_old[i])
        f.write("\n")
      f.close()

  mpi.world.barrier() # wait until solver input is written

# go back and loop
