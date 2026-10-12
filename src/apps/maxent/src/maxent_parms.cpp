/*****************************************************************************
 *
 * ALPS Project Applications
 *
 * Copyright (C) 2010 by Sebastian Fuchs <fuchs@comp-phys.org>
 *                       Thomas Pruschke <pruschke@comp-phys.org>
 *                       Matthias Troyer <troyer@comp-phys.org>
 *
 * ALPS Project: https://alps.comp-phys.org/
 * SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

#include <utility> // Boost.Math 1.76 includes this inside a namespace.
#include <boost/math/constants/constants.hpp>
#include "maxent.hpp"
#include <alps/config.h> // needed to set up correct bindings
#include <alps/hdf5/vector.hpp>
#include <boost/lexical_cast.hpp>
#include <limits>
#include <sstream>
#include <boost/numeric/bindings/ublas.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/bindings/lapack/driver/syev.hpp>
#include <boost/numeric/bindings/lapack/driver/gesvd.hpp>
#include <boost/numeric/bindings/upper.hpp>

#define MAXIMUM(a,b) ((a>b) ? a : a)

namespace {
int frequency_count(const alps::params& p) {
  const int count = p["NFREQ"];
  if (count < 1 || count == std::numeric_limits<int>::max())
    throw std::invalid_argument("NFREQ must be positive and leave room for an endpoint");
  return count;
}
}

ContiParameters::ContiParameters(const alps::params& p)
: ContiParameters(p, maxent_detail::legacy_data(p)) {}

ContiParameters::ContiParameters(const alps::params& p, const maxent_detail::data& data) :
Default_(make_default_model(p, "DEFAULT_MODEL", data.prior_omega, data.prior_density)),
T_(p.defined("T") ? static_cast<double>(p["T"]) : 1./static_cast<double>(p["BETA"])),
ndat_(maxent_detail::data_count(data)), nfreq_(frequency_count(p)),
y_(ndat_),sigma_(ndat_), x_(ndat_),K_(),t_array_(nfreq_+1),tau_(data.tau)
{
  if (!std::isfinite(T_) || T_ <= 0 || !std::isfinite(1./T_))
    throw std::invalid_argument("temperature and inverse temperature must be finite and positive");
  const double norm = p["NORM"] | 1.;
  if (!std::isfinite(norm) || norm <= 0)
    throw std::invalid_argument("NORM must be finite and positive");
  for (double tau : tau_)
    if (tau < 0 || tau > 1./T_) throw std::invalid_argument("tau must lie in [0,BETA]");
  std::string p_f_grid = p["FREQUENCY_GRID"]|"";
  if (p_f_grid=="Lorentzian" || p_f_grid=="half Lorentzian") {
    const double cut = p["CUT"]|0.01;
    if (!(cut > 0. && cut < 0.5))
      throw std::invalid_argument("CUT must lie strictly between 0 and 0.5");
    const bool half = p_f_grid=="half Lorentzian";
    std::vector<double> temp(nfreq_+1);
    for (int i=0; i<=nfreq_; ++i) {
      const double position = double(i)/nfreq_;
      const double angle = half ? position*(0.5-cut) : position*(1.-2*cut)+cut-0.5;
      temp[i] = std::tan(boost::math::constants::pi<double>() * angle);
    }
    for (int i=0; i<=nfreq_; ++i)
      t_array_[i] = (temp[i] - temp.front())/(temp.back() - temp.front());
  }
  else if (p_f_grid=="quadratic") {
    const double spread = p["SPREAD"]|4;
    if (!std::isfinite(spread) || spread < 1.)
      throw std::invalid_argument("SPREAD must be finite and at least 1");
    // Integrate all NFREQ positive interval weights, then normalize the knots.
    t_array_[0] = 0.;
    for (int i=0; i<nfreq_; ++i) {
      const double a = nfreq_ == 1 ? 0.5 : double(i)/(nfreq_-1);
      const double weight = 1. + (spread-1.) * std::pow(2.*a-1., 2);
      t_array_[i+1] = t_array_[i] + weight;
    }
    const double total = t_array_[nfreq_];
    for (int i=1; i<=nfreq_; ++i) t_array_[i] /= total;
  }
  else if (p_f_grid=="log") {
    const double t_min = p["LOG_MIN"]|1.0e-4;
    if (!(t_min > 0. && t_min < 0.5))
      throw std::invalid_argument("LOG_MIN must lie strictly between 0 and 0.5");
    // Even NFREQ has a central knot; odd NFREQ has a central interval.
    const int pairs = (nfreq_+1)/2;
    if (nfreq_%2 == 0) t_array_[nfreq_/2] = 0.5;
    for (int i=0; i<pairs; ++i) {
      const double distance = pairs == 1 ? 0.5 :
          std::exp(std::log(t_min) + double(i)/(pairs-1)*(std::log(0.5)-std::log(t_min)));
      t_array_[pairs-1-i] = 0.5-distance;
      t_array_[nfreq_-(pairs-1-i)] = 0.5+distance;
    }
  }
  else if (p_f_grid=="linear") {
    // Every grid supplies NFREQ+1 knots for NFREQ integration intervals.
    for (int i=0; i<=nfreq_; ++i) t_array_[i] = double(i)/nfreq_;
  }
  else {
    throw std::invalid_argument("No valid frequency grid specified");
  }
  t_array_[0] = 0.;
  t_array_[nfreq_] = 1.;
  for (int i=1; i<=nfreq_; ++i) {
    if (!std::isfinite(t_array_[i]) || !(t_array_[i] > t_array_[i-1]))
      throw std::invalid_argument("Frequency grid knots must be finite and strictly increasing");
  }
  for (int i = 0; i < ndat_; ++i) {
    y_(i) = data.values[i] / norm;
    if (data.covariance.empty()) sigma_(i) = data.errors[i] / norm;
    if (!std::isfinite(y_(i)) || (data.covariance.empty() &&
        (!std::isfinite(sigma_(i)) || sigma_(i) <= 0)))
      throw std::invalid_argument("data normalization is not representable");
  }
  if (!data.covariance.empty()) {
    cov_.resize(ndat_, ndat_);
    for (int i = 0; i < ndat_; ++i) for (int j = 0; j < ndat_; ++j)
      cov_(i,j) = data.covariance[std::size_t(i)*ndat_+j];
  }
}


void ContiParameters::setup_kernel(const alps::params& p, const int ntab, const vector_type& freq)
{
  using namespace boost::numeric;
  K_.resize(ndat_, ntab);
  std::string p_data = p["DATASPACE"]|"time";
  std::string p_kernel = p["KERNEL"]|"fermionic";
//  for (int i=0; i<ndat(); ++i)
//    y_(i) = static_cast<double>(p["X_"+boost::lexical_cast<std::string>(i)])/static_cast<double>(p["NORM"] | 1.);
  if(p_data=="time") {
    if (alps::is_master())
      std::cerr << "assume time space data" << std::endl;
    if (p_kernel == "fermionic") {
      if (alps::is_master())
        std::cerr << "Using fermionic kernel" << std::endl;
      for (int i=0; i<ndat(); ++i) {
        double tau;
        if (!tau_.empty())
          tau = tau_[i];
        else
          tau = i / ((ndat()-1)* T_);
        for (int j=0; j<ntab; ++j) {
          double omega = freq[j]; //Default().omega_of_t(double(j)/(ntab-1));
          K_(i,j) =  -1. / (std::exp(omega*tau) + std::exp(-omega*(1./T_-tau)));
        }
      }
    }
    else if (p_kernel == "bosonic") {
      if (alps::is_master())
        std::cerr << "Using bosonic kernel" << std::endl;
      for (int i=0; i<ndat(); ++i) {
        double tau;
        if (!tau_.empty())
          tau = tau_[i];
        else
          tau = i / ((ndat()-1) * T_);
        for (int j=0; j<ntab; ++j) {
          // This kernel is even in omega. Use nonpositive exponent arguments
          // and expm1 for a stable zero-frequency limit in every column.
          const double omega = std::abs(freq[j]);
          K_(i,j) = omega == 0 ? T_ : 0.5*omega *
              (std::exp(-omega*tau) + std::exp(-omega*(1./T_-tau))) /
              -std::expm1(-omega/T_);
        }
      }
    }
    //for zero temperature, only positive frequency matters
    else if (p_kernel == "Boris") {
      if (alps::is_master())
        std::cerr << "Using Boris' kernel" << std::endl;
      for (int i=0; i<ndat(); ++i) {
        if (tau_.empty()) throw std::invalid_argument("Boris kernel requires explicit tau");
        double tau = tau_[i];
        for (int j=0; j<ntab; ++j) {
          double omega = freq[j];
          K_(i,j) = -std::exp(-omega*tau);
        }
      }
    }
    else
      boost::throw_exception(std::invalid_argument("unknown integration kernel"));
  }
  else if (p_data == "frequency" && p_kernel == "fermionic" &&
           (p["PARTICLE_HOLE_SYMMETRY"]|false)) {
    std::cerr << "using particle hole symmetric kernel for fermionic data" << std::endl;
    for (int i=0; i<ndat(); ++i) {
      double omegan = (2*i+1)*boost::math::constants::pi<double>()*T_;
      for (int j=0; j<ntab; ++j) {
        double omega = freq[j];
        K_(i,j) =  -omegan / (omegan*omegan + omega*omega);
      }
    }
  }
  else if (p_data == "frequency" && p_kernel == "bosonic" &&
           (p["PARTICLE_HOLE_SYMMETRY"]|false)) {
    //std::cerr << "using particle hole symmetric kernel for bosonic data" << std::endl;
    //std::cerr<<"ndat is: "<<ndat()<<" ntab: "<<ntab<<std::endl;
    //std::cerr<<"freqs: "<<freq[0]<<" "<<freq[ntab-1]<<std::endl;

    for (int i=0; i<ndat(); ++i) {
      double Omegan = (2*i)*boost::math::constants::pi<double>()*T_;
      for (int j=0; j<ntab; ++j) {
        double Omega = freq[j];
        K_(i,j) = i == 0 ? -1. : -Omega*Omega / (Omegan*Omegan + Omega*Omega);
      }
    }
    //double z=0;
    //for (int i=0; i<ndat(); ++i) {
    //  for (int j=0; j<ntab; ++j) {
    //    z+=K_(i,j);
    //  }
    //}
    //std::cout<<"debug kernel checksum is: "<<z<<std::endl;
  }
  else if (p_data == "frequency" && p_kernel == "anomalous" &&
           (p["PARTICLE_HOLE_SYMMETRY"]|false)) {
    std::cerr << "using particle hole symmetric kernel for anomalous fermionic data" << std::endl;
    for(int i=0;i<ndat();++i){
      double omegan = (2*i+1)*boost::math::constants::pi<double>()*T_;
      for (int j=0; j<ntab; ++j) {
        double omega = freq[j];
        K_(i,j) =  omega*omega / (omegan*omegan + omega*omega);
      }
    }
  }
  else if (p_data == "frequency") {
    if (alps::is_master())
      std::cerr << "assume frequency space data" << std::endl;
    if (ndat_ % 2)
      throw std::invalid_argument("complex frequency data require real/imaginary pairs");
    // First-column corrections also identified in ALPS reconciliation
    // 90a36c0d4, src/apps/maxent/src/maxent_parms.cpp. Keep the legacy API.
    ublas::matrix<std::complex<double>, ublas::column_major> Kc(ndat_/2, ntab);
    if (p_kernel == "fermionic") {
      if (alps::is_master())
        std::cerr << "Using fermionic kernel" << std::endl;
      for (int i=0; i<ndat()/2; ++i) {
        std::complex<double> iomegan(0, (2*i+1)*boost::math::constants::pi<double>()*T_);
        for (int j=0; j<ntab; ++j) {
          double omega = freq[j];
          Kc(i,j) =  1. / (iomegan - omega);
        }
      }
    }
    else if (p_kernel == "bosonic") {
      if (alps::is_master())
        std::cerr << "Using bosonic kernel" << std::endl;
      for (int i=0; i<ndat()/2; ++i) {
        std::complex<double> iomegan(0, 2*i*boost::math::constants::pi<double>()*T_);
        for (int j=0; j<ntab; ++j) {
          double omega = freq[j];
          //Kc(i,j) =  -1. / (iomegan - omega);
          Kc(i,j) = i == 0 ? -1. : omega / (iomegan - omega);
        }
      }
    }
    else if (p_kernel == "anomalous"){
      std::cerr<<"Using general anomalous kernel omega / (iomega_n - omega) for, e.g., omega*Delta"<<std::endl;
      for (int i=0; i<ndat()/2; ++i) {
        std::complex<double> iomegan(0, (2*i+1)*boost::math::constants::pi<double>()*T_);
        for (int j=0; j<ntab; ++j) {
          double omega = freq[j];
          Kc(i,j) =  -omega / (iomegan - omega);
        }
      }
    }
    else
      boost::throw_exception(std::invalid_argument("unknown integration kernel"));
    for (int i=0; i<ndat(); i+=2) {
      for (int j=0; j<ntab; ++j) {
        K_(i,j) = Kc(i/2,j).real();
        K_(i+1,j) = Kc(i/2,j).imag();
      }
    }
  }
  else
    boost::throw_exception(std::invalid_argument("unknown value for parameter DATASPACE"));
  for (int i = 0; i < ndat(); ++i)
    for (int j = 0; j < ntab; ++j)
      if (!std::isfinite(K_(i,j)))
        throw std::invalid_argument("nonfinite integration kernel");
//  vector_type sigma(ndat());
  if (cov_.size1() != 0) {
    vector_type var(ndat());
    if (bindings::lapack::syev('V', bindings::upper(cov_), var, bindings::lapack::optimal_workspace()) != 0)
      throw std::runtime_error("covariance eigendecomposition failed");
    const double tolerance = 64 * std::numeric_limits<double>::epsilon() * ndat() * std::abs(var[var.size()-1]);
    if (!std::isfinite(var[0]) || var[0] < -tolerance)
      throw std::invalid_argument("covariance must be positive semidefinite");
    matrix_type cov_trans = ublas::trans(cov_);
    matrix_type K_loc = ublas::prec_prod(cov_trans, K_);
    vector_type y_loc = ublas::prec_prod(cov_trans, y_);
    if (alps::is_master() && p["VERBOSE"]|false)
      std::cout << "# Eigenvalues of the covariance matrix:\n";
    // We drop eigenvalues of the covariance matrix which are smaller than 1e-10
    // as they represent bad data directions (usually there is a steep drop
    // below that value)
    int new_ndat_,old_ndat_=ndat();
    for (new_ndat_ =0;new_ndat_<ndat();new_ndat_++)
      if (var[new_ndat_]>1e-10) break;
    // This is the number of good data
    ndat_ = old_ndat_ - new_ndat_;
    if (ndat_ == 0) throw std::invalid_argument("covariance has no usable directions above 1e-10");
    if (alps::is_master())
      std::cout << "# Ignoring singular eigenvalues (0-" << new_ndat_-1 << " out of " << old_ndat_ << ")\n";
    // Now resize kernel and data matrix and fill it with the values for the
    // good data directions
    K_.resize(ndat_,ntab);
    y_.resize(ndat_);
    sigma_.resize(ndat_);
    for (int i=0; i<ndat(); i++) {
      y_(i) = y_loc(new_ndat_+i);
      for (int j=0; j<ntab; j++) {
        K_(i,j)=K_loc(new_ndat_+i,j);
      }
    }
    for (int i=0; i<ndat(); ++i) {
        sigma_[i] = std::sqrt(std::abs(var(new_ndat_+i)))/static_cast<double>(p["NORM"] | 1.);
      if (alps::is_master() && p["VERBOSE"]|false)
        std::cout << "# " << var(new_ndat_+i) << "\n";
    }
  }
  //else {
  //  for (int i=0; i<ndat(); ++i)
  //    sigma[i] = static_cast<double>(p["SIGMA_"+boost::lexical_cast<std::string>(i)])/static_cast<double>(p["NORM"] | 1.);
  //}
  //Look around Eq. D.5 in Sebastian's thesis. We have sigma = sqrt(eigenvalues of covariance matrix) or, in case of a diagonal covariance matrix, we have sigma=SIGMA_X. The then define y := \bar{G}/sigma and K := (1/sigma)\tilde{K}
  if (p_data == "hillibilli" && !(p["PARTICLE_HOLE_SYMMETRY"]|true))  {
//      if (p["DATASPACE"]=="frequency" && !p.value_or_default("PARTICLE_HOLE_SYMMETRY",true))  {
    if (alps::is_master())
          std::cerr << "Kernel for complex data\n";
    for (int i=0; i<ndat()/2; i+=2) {
      std::complex<double> y(y_[i],y_[i+1]),s(sigma_[i],sigma_[i+1]);
      y /= s;
      y_[i] = y.real();
      y_[i+1] = y.imag();
      for (int j=0; j<ntab; ++j) {
        std::complex<double> K(K_(i,j),K_(i+1,j));
        K /= s;
        K_(i,j) = K.real();
        K_(i+1,j) = K.imag();
      }
    }
  } else {
    for (int i=0; i<ndat(); i++) {
      y_[i] /= sigma_[i];
      for (int j=0; j<ntab; ++j) {
        K_(i,j) /= sigma_[i];
      }
    }
  }

  //this enforces a strict normalization if needed.
  //not sure that this is done properly. recheck!
  if(p["ENFORCE_NORMALIZATION"]|false) {
    std::cout<<"enforcing strict normalization."<<std::endl;
    double artificial_norm_enforcement_sigma=static_cast<double>(p["SIGMA_NORMALIZATION"])/static_cast<double>(p["NORM"] | 1.);
    for(int j=0;j<ntab;++j){
      K_(ndat()-1,j) = 1./artificial_norm_enforcement_sigma;
    }
    y_[ndat()-1]=1./artificial_norm_enforcement_sigma;
  }
    for (int i=0; i<ndat(); ++i) {
      if (!std::isfinite(y_[i])) throw std::invalid_argument("nonfinite whitened data");
      for (int j=0; j<ntab; ++j)
        if (!std::isfinite(K_(i,j))) throw std::invalid_argument("nonfinite whitened kernel");
    }
    std::cerr << "Kernel set up\n";
}


MaxEntParameters::MaxEntParameters(const alps::params& p) :
ContiParameters(p),
U_(ndat(), ndat()), Vt_(ndat(), nfreq()), Sigma_(ndat(), ndat()),
omega_coord_(nfreq()), delta_omega_(nfreq()), ns_(0)
{
  using namespace boost::numeric;
//  if (ndat() > nfreq())
//    boost::throw_exception(std::invalid_argument("NDAT should be smaller than NFREQ"));
/*  for (int i=0; i<nfreq(); ++i) {
    omega_coord_[i] = Default().omega_of_t(t_array_[i]); //(Default().omega_of_t(t_array_[i]) + Default().omega_of_t(t_array_[i+1]))/2.;
    if (i>0 && i<nfreq()-1)
      delta_omega_[i] = (Default().omega_of_t(t_array_[i+1]) - Default().omega_of_t(t_array_[i-1]))/2.0; //    delta_omega_[i] = (Default().omega_of_t(t_array_[MAXIMUM(i+1,nfreq()-1)]) - Default().omega_of_t(t_array_[MAXIMUM(0,i-1)]))/2.0; //Default().omega_of_t(t_array_[i+1]) - Default().omega_of_t(t_array_[i]);
    else if (i==0)
      delta_omega_[i] = (Default().omega_of_t(t_array_[i+1]) - Default().omega_of_t(t_array_[i]))/2.0; //    delta_omega_[i] = (Default().omega_of_t(t_array_[MAXIMUM(i+1,nfreq()-1)]) - Default().omega_of_t(t_array_[MAXIMUM(0,i-1)]))/2.0; //Default().omega_of_t(t_array_[i+1]) - Default().omega_of_t(t_array_[i]);
    else
      delta_omega_[i] = (Default().omega_of_t(t_array_[i]) - Default().omega_of_t(t_array_[i-1]))/2.0; //    delta_omega_[i] = (Default().omega_of_t(t_array_[MAXIMUM(i+1,nfreq()-1)]) - Default().omega_of_t(t_array_[MAXIMUM(0,i-1)]))/2.0; //Default().omega_of_t(t_array_[i+1]) - Default().omega_of_t(t_array_[i]);
//
//      delta_omega_[i] = (Default().omega_of_t(t_array_[MAXIMUM(i+1,nfreq()-1)]) - Default().omega_of_t(t_array_[MAXIMUM(0,i-1)]))/2.0; //    delta_omega_[i] = (Default().omega_of_t(t_array_[MAXIMUM(i+1,nfreq()-1)]) - Default().omega_of_t(t_array_[MAXIMUM(0,i-1)]))/2.0; //Default().omega_of_t(t_array_[i+1]) - Default().omega_of_t(t_array_[i]);
*/
//      std::cout << i << " " << omega_coord_[i] << " " << delta_omega_[i] << std::endl;
//  }
  for (int i=0; i<nfreq(); ++i) {
    omega_coord_[i] = (Default().omega_of_t(t_array_[i]) + Default().omega_of_t(t_array_[i+1]))/2.;
    delta_omega_[i] = Default().omega_of_t(t_array_[i+1]) - Default().omega_of_t(t_array_[i]);
  }
  //std::cerr<<"debug: omega and delta: "<<std::endl;
  //for(int i=0;i<nfreq();++i){
  //  std::cout<<omega_coord_[i]<<" "<<delta_omega_[i]<<std::endl;
  //}
  setup_kernel(p, nfreq(), omega_coord_);

  //perform the SVD decomposition K = U Sigma V^T
  const int rank = std::min(ndat(), nfreq());
  vector_type S(rank);
  U_.resize(ndat(), rank); Vt_.resize(rank, nfreq());
  matrix_type Kt = K_; // gesvd destroys K!
  if (bindings::lapack::gesvd('S','S',Kt, S, U_, Vt_) != 0)
    throw std::runtime_error("kernel singular value decomposition failed");
  if (p["VERBOSE"]|false) std::cout << "# Singular values of the Kernel:\n";
    const double prec = std::sqrt(std::numeric_limits<double>::epsilon())*nfreq()*S[0];
  if (p["VERBOSE"]|false) std::cout << "# eps = " << sqrt(std::numeric_limits<double>::epsilon()) << std::endl << "# prec = " << prec << std::endl;
  for (unsigned int s=0; s<S.size(); ++s) {
    if (p["VERBOSE"]|false)std::cout << "# " << s << "\t" << S[s] <<"\n";
    ns_ = (std::isfinite(S[s]) && S[s] > 0 && S[s] >= prec) ? s+1 : ns_;
  }
  if (ns() == 0)
    boost::throw_exception(std::logic_error("all singular values smaller than the precision"));

  //(truncated) U has dimension ndat() * ns_; ndat() is # of input (matsubara frequency/imag time) points
  //(truncated) Sigma has dimension ns_*ns_ (number of singular eigenvalues)
  //(truncated) V^T has dimensions ns_* nfreq(); nfreq() is number of output (real) frequencies
  //ns_ is the dimension of the singular space.

  U_.resize(ndat(), ns_, true);
  Vt_.resize(ns_, nfreq(), true);
  Sigma_.resize(ns_, ns_);
  Sigma_.clear();
  for (int s=0; s<ns_; ++s) {
    std::cout << "# " << s << "\t" << S[s] <<"\n";
    Sigma_(s,s) = S[s];
  }
  //compute Ut and
  matrix_type Ut = ublas::trans(U_);             //U^T has dimension ns_*ndat()
  vector_type t = ublas::prec_prod(Ut, y_);      //t has dimension ns_
  vector_type y2 = ublas::prec_prod(U_, t);      //y2 has dimension ndat(), which is dimension of y
  double chi = ublas::norm_2(y_-y2);             //this measures the loss of precision when transforming to singular space and back.
  std::cout << "minimal chi2: " << chi*chi/y_.size() << std::endl;
}



