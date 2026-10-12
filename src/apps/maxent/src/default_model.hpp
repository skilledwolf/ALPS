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


#ifndef ALPS_TOOL_DEFAULT_MODEL_HPP
#define ALPS_TOOL_DEFAULT_MODEL_HPP

#include <utility> // Boost.Math 1.76 includes this inside a namespace.
#include <boost/math/constants/constants.hpp>
#include <math.h>
#include <cmath>
#include <alps/ngs/params.hpp>
#include <alps/osiris/comm.h>
#include <boost/shared_ptr.hpp>
#include <boost/throw_exception.hpp>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>


namespace maxent_detail {
inline double finite_parameter(const alps::params& p, const char* key) {
  const double x = p[key];
  if (!std::isfinite(x)) throw std::invalid_argument(std::string(key)+" must be finite");
  return x;
}
inline double positive_parameter(const alps::params& p, const char* key) {
  const double x = finite_parameter(p,key);
  if (x <= 0) throw std::invalid_argument(std::string(key)+" must be positive");
  return x;
}
}

//Note the slightly crooked structure here:

//class DefaultModel
// -> FlatDefaultModel : public DefaultModel
// -> GeneralDefaultModel : public DefaultModel
//  and GeneralDefaultModel contains a 'Model' object.

//class Model
// -> Gaussian : public Model
//    -> ShiftedGaussian : public Gaussian
//      -> DoubleGaussian : public ShiftedGaussian
//        -> GeneralDoubleGaussian : public ShiftedGaussian
// ->  TabFunction : public Model
// ->  LinearRiseExpDecay : public Model
// ->  QuadraticRiseExpDecay : public Model

class DefaultModel
{
public:
//  DefaultModel(const alps::Parameters& p) :
  DefaultModel(const alps::params& p) :
    omega_max(p["OMEGA_MAX"]),
    omega_min(static_cast<double>(p["OMEGA_MIN"]|(-omega_max))), //we had a 0 here in the bosonic case. That's not a good idea if you're continuing symmetric functions like chi(omega)/omega. Change omega_min to zero manually if you need it.
    blow_up_(p["BLOW_UP"]|1.)
  {
    if (!std::isfinite(omega_min) || !std::isfinite(omega_max) ||
        !std::isfinite(omega_max-omega_min) || omega_min >= omega_max ||
        !std::isfinite(blow_up_) || blow_up_ <= 0)
      throw std::invalid_argument("invalid default-model frequency range or BLOW_UP");
    //std::cout<<"found omega_min: "<<omega_min<<std::endl;
    //std::cout<<"found omega_max: "<<omega_max<<std::endl;
    //std::cout<<"found blowup:    "<<blow_up_<<std::endl;
  }

  virtual ~DefaultModel(){}

  //omega returns a frequency point, given x between 0 and 1.
  virtual double omega(const double x) const = 0;

  //D returns the derivative of the integrated default model
  virtual double D(const double omega) const = 0;

  //returns the integrated default model
  virtual double x(const double t=0) const = 0;

  //equidistant mapping from [0,1] to [omega_min, omega_max]
  double omega_of_t(const double t) const { return omega_min + (omega_max-omega_min)*t; }
  //equidistant mapping from [omega_min, omega_max] to [0,1]
  double t_of_omega(const double omega) const { return (omega-omega_min)/(omega_max-omega_min); }
  double blow_up() const { return blow_up_; }

protected:
  const double omega_max;
  const double omega_min;
  const double blow_up_;
};



class FlatDefaultModel : public DefaultModel
{
public:

//    FlatDefaultModel(const alps::Parameters& p) : DefaultModel(p) {}
  FlatDefaultModel(const alps::params& p) : DefaultModel(p) {}

  double omega(const double x) const {
    return x/blow_up()*(omega_max-omega_min) + omega_min;
  }

  double D(const double) const {
    return 1./(omega_max-omega_min);
  }

  double x(const double t) const {
    return t*blow_up();
  }

};



class Model
{
public:
  virtual double operator()(const double omega)=0;
  virtual ~Model(){}
};



class Gaussian : public Model
{
public:
//  Gaussian(const alps::Parameters& p) : sigma(maxent_detail::positive_parameter(p,"SIGMA")) {}
  Gaussian(const alps::params& p) : sigma(maxent_detail::positive_parameter(p,"SIGMA")) {}

  virtual double operator()(const double omega) {
    return std::exp(-omega*omega/2./sigma/sigma)/sqrt(2*boost::math::constants::pi<double>())/sigma;
  }

private:
  const double sigma;
};

class TwoGaussians : public Model
{
public:
//    TwoGaussians(const alps::Parameters& p) : sigma1(maxent_detail::positive_parameter(p,"SIGMA1")),
    TwoGaussians(const alps::params& p) : sigma1(maxent_detail::positive_parameter(p,"SIGMA1")),
    sigma2(maxent_detail::positive_parameter(p,"SIGMA2")),
    shift1(static_cast<double>(p["SHIFT1"]|0.0)),
    shift2(maxent_detail::finite_parameter(p,"SHIFT2")),
    norm1(static_cast<double>(p["NORM1"]|0.5)) {
      if (!std::isfinite(shift1) || !std::isfinite(norm1) || norm1 < 0 || norm1 > 1)
        throw std::invalid_argument("invalid two-Gaussian shift or weight");
    }

    virtual double operator()(const double omega) {
        return norm1*std::exp(-(omega-shift1)*(omega-shift1)/2./sigma1/sigma1)/sqrt(2*boost::math::constants::pi<double>())/sigma1+(1.0-norm1)*std::exp(-(omega-shift2)*(omega-shift2)/2./sigma2/sigma2)/sqrt(2*boost::math::constants::pi<double>())/sigma2;
    }

private:
    const double sigma1,sigma2,shift1,shift2,norm1;
};


class ShiftedGaussian : public Gaussian
{
public:
//  ShiftedGaussian(const alps::Parameters& p) :
  ShiftedGaussian(const alps::params& p) :
    Gaussian(p), shift(maxent_detail::finite_parameter(p,"SHIFT")){}

  double operator()(const double omega) {
    return Gaussian::operator()(omega-shift);
  }

protected:
  const double shift;
};


class DoubleGaussian : public ShiftedGaussian
{
public:
//  DoubleGaussian(const alps::Parameters& p) :
  DoubleGaussian(const alps::params& p) :
    ShiftedGaussian(p){}

  double operator()(const double omega) {
    return 0.5*(Gaussian::operator()(omega-shift) + Gaussian::operator()(omega+shift));
  }
};

class LinearRiseExpDecay : public Model{
public:
//  LinearRiseExpDecay(const alps::Parameters &p): lambda_(maxent_detail::positive_parameter(p,"LAMBDA")){}
  LinearRiseExpDecay(const alps::params &p): lambda_(maxent_detail::positive_parameter(p,"LAMBDA")){}
  double operator()(const double omega) {
    return lambda_*lambda_*omega*std::exp(-lambda_*omega);
  }

private:
  const double lambda_;
};

class QuadraticRiseExpDecay : public Model{
public:
//  QuadraticRiseExpDecay(const alps::Parameters &p): lambda_(maxent_detail::positive_parameter(p,"LAMBDA")){}
  QuadraticRiseExpDecay(const alps::params &p): lambda_(maxent_detail::positive_parameter(p,"LAMBDA")){}
  double operator()(const double omega) {
    return (lambda_*lambda_*lambda_)/2.*(omega*omega)*std::exp(-lambda_*omega);
  }

private:
  const double lambda_;
};

class GeneralDoubleGaussian : public ShiftedGaussian
{
public:
//  GeneralDoubleGaussian(const alps::Parameters& p) :
  GeneralDoubleGaussian(const alps::params& p) :
    ShiftedGaussian(p), bnorm(maxent_detail::positive_parameter(p,"BOSE_NORM")) {}

  double operator()(const double omega) {
    if (omega > 0)
      return Gaussian::operator()(omega);
    else
      return bnorm*Gaussian::operator()(omega+shift);
  }

private:
  const double bnorm;
};


class TabFunction : public Model
{
public:
//  TabFunction(const alps::Parameters& p, std::string const& name) //: index(0)
  TabFunction(const alps::params& p, std::string const& name) //: index(0)
  {
    std::string p_name = p[name].cast<std::string>();
    std::ifstream defstream(p_name.c_str());
    if (!defstream)
      boost::throw_exception(std::invalid_argument("could not open default model file: "+p[name]));
    while (defstream >> std::ws && !defstream.eof()) {
      double om, density;
      if (!(defstream >> om >> density)) throw std::invalid_argument("malformed tabulated prior");
      Omega.push_back(om); Def.push_back(density);
      defstream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    validate(p);
  }

  TabFunction(const alps::params& p, std::vector<double> omega, std::vector<double> density)
    : Omega(std::move(omega)), Def(std::move(density)) { validate(p); }


  double operator()(const double omega) {
    if (!std::isfinite(omega) || omega < Omega.front() || omega > Omega.back())
      throw std::invalid_argument("frequency outside tabulated prior");
    if (omega == Omega.front()) return Def.front();
    std::vector<double>::const_iterator ub = std::upper_bound(Omega.begin(), Omega.end(), omega);
    int index = ub - Omega.begin();
    if (ub==Omega.end())
      index = Omega.end()-Omega.begin()-1;
    double om1 = Omega[index-1];
    double om2 = Omega[index];
    double D1 = Def[index-1];
    double D2 = Def[index];
    return -(D2-D1)/(om2-om1)*(om2-omega)+D2;
  }

private:
  void validate(const alps::params& p) const {
    if (Omega.size() < 2 || Omega.size() != Def.size())
      throw std::invalid_argument("tabulated prior requires matching arrays with at least two points");
    bool positive = false;
    for (std::size_t i = 0; i < Omega.size(); ++i) {
      if (!std::isfinite(Omega[i]) || !std::isfinite(Def[i]) || Def[i] < 0 ||
          (i && Omega[i] <= Omega[i-1]))
        throw std::invalid_argument("prior coordinates must increase and densities be finite and nonnegative");
      positive |= Def[i] > 0;
    }
    const double hi = p["OMEGA_MAX"], lo = p["OMEGA_MIN"] | -hi;
    if (!positive || Omega.front() > lo || Omega.back() < hi)
      throw std::invalid_argument("tabulated prior must have positive mass and cover the frequency interval");
  }
  std::vector<double> Omega;
  std::vector<double> Def;
};



class GeneralDefaultModel : public DefaultModel
{
public:

//  GeneralDefaultModel(const alps::Parameters& p, boost::shared_ptr<Model> mod)
  GeneralDefaultModel(const alps::params& p, boost::shared_ptr<Model> mod)
   : DefaultModel(p)
   , Mod(mod)
   , ntab(5001)
   , xtab(ntab)
  {
    double sum = 0;
    xtab[0] = 0.;
    //this is an evaluation on an equidistant grid; sum integrated by trapezoidal rule
    double delta_omega = (omega_max-omega_min)/(ntab-1);
    for (int o=1; o<ntab; ++o) {
      double omega1 = omega_min + (o-1)*delta_omega;
      double omega2 = omega_min + o*delta_omega;
      const double left = (*Mod)(omega1), right = (*Mod)(omega2);
      if (!std::isfinite(left) || !std::isfinite(right) || left < 0 || right < 0)
        throw std::invalid_argument("default-model density must be finite and nonnegative");
      sum += (left+right)/2.*delta_omega;
      xtab[o] = sum;
    }
    if (!std::isfinite(sum) || sum <= 0 || !std::isfinite(blow_up()/sum))
      throw std::invalid_argument("default-model integral must be finite and positive");
    for (int o=0; o<ntab; ++o) {
      xtab[o] *= blow_up()/sum;
    }
    /*std::cout<<std::setprecision(14)<<"tabulated a model, these are the values: "<<std::endl;
    for(int i=0;i<ntab;++i){
      std::cout<<i<<" "<<xtab[i]<<std::endl;
    }
    std::cout<<"total sum is: "<<sum<<std::endl;*/

    /*for(int N=0;N<=1000;++N){
      double d=1./1000.*N;
      std::cout<<omega(d)<<" "<<D(omega(d))<<" "<<x(d)<<" "<<d<<std::endl;
    }*/
  }

  double omega(const double x) const {
    if(!(x<=blow_up() && x>=0.)) throw std::logic_error("parameter x is out of bounds!"); //DNDEBUG switches off debug assertions
    if (x == 0) return omega_min;
    if (x == blow_up()) return omega_max;
    std::vector<double>::const_iterator ub = std::upper_bound(xtab.begin(), xtab.end(), x);
    int omega_index = ub - xtab.begin();
    if (ub==xtab.end())
      omega_index = xtab.end() - xtab.begin() - 1;
    double om1 = omega_min + (omega_index-1)*(omega_max-omega_min)/(ntab-1);
    double om2 = omega_min + omega_index*(omega_max-omega_min)/(ntab-1);
    double x1 = xtab[omega_index-1];
    double x2 = xtab[omega_index];
    return -(om2-om1)/(x2-x1)*(x2-x)+om2;
  }

  //this returns the value of the model function at frequency omega
  double D(const double omega) const {
    return (*Mod)(omega);
  }

  //I have no idea what this does.
  double x(const double t) const {
    if(!std::isfinite(t) || t>1. || t<0.) throw std::logic_error("parameter t is out of bounds!");
    int od = (int)(t*(ntab-1));
    if (od==(ntab-1))
      return blow_up();
    double x1 = xtab[od];
    double x2 = xtab[od+1];
    return x1 + (x2-x1)*(t*(ntab-1)-od);
  }

private:
  boost::shared_ptr<Model> Mod;
  const int ntab;
  std::vector<double> xtab; //xtab has an equidistantly tabulated discretized model function
};



//inline boost::shared_ptr<DefaultModel> make_default_model(const alps::Parameters& parms, std::string const& name)
inline boost::shared_ptr<DefaultModel> make_default_model(const alps::params& parms, std::string const& name,
    const std::vector<double>& prior_omega = {}, const std::vector<double>& prior_density = {})
{
    std::string p_name = parms[name]|"flat";
  if (p_name == "tabulated") {
    boost::shared_ptr<Model> mod(new TabFunction(parms, prior_omega, prior_density));
    return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, mod));
  }
  if (!prior_omega.empty() || !prior_density.empty())
    throw std::invalid_argument("prior arrays require DEFAULT_MODEL=tabulated");
  if (p_name == "flat") {
    if (alps::is_master())
      std::cerr << "Using flat default model" << std::endl;
    return boost::shared_ptr<DefaultModel>(new FlatDefaultModel(parms));
  }
  else if (p_name == "gaussian") {
    if (alps::is_master())
      std::cerr << "Using Gaussian default model" << std::endl;
    boost::shared_ptr<Model> Mod(new Gaussian(parms));
    return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "twogaussians") {
      if (alps::is_master())
          std::cerr << "Using sum of two Gaussians default model" << std::endl;
      boost::shared_ptr<Model> Mod(new TwoGaussians(parms));
      return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "shifted gaussian") {
    if (alps::is_master())
      std::cerr << "Using shifted Gaussian default model" << std::endl;
    boost::shared_ptr<Model> Mod(new ShiftedGaussian(parms));
    return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "double gaussian") {
    if (alps::is_master())
      std::cerr << "Using double Gaussian default model" << std::endl;
    boost::shared_ptr<Model> Mod(new DoubleGaussian(parms));
    return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "general double gaussian") {
    if (alps::is_master())
      std::cerr << "Using general double Gaussian default model" << std::endl;
    boost::shared_ptr<Model> Mod(new GeneralDoubleGaussian(parms));
    return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "linear rise exp decay") {
    if (alps::is_master())
      std::cerr << "Using linear rise exponential decay default model" << std::endl;
    boost::shared_ptr<Model> Mod(new LinearRiseExpDecay(parms));
    return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else if (p_name == "quadratic rise exp decay") {
    if (alps::is_master())
      std::cerr << "Using quadratic rise exponential decay default model" << std::endl;
    boost::shared_ptr<Model> Mod(new QuadraticRiseExpDecay(parms));
    return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
  else {
    if (alps::is_master())
      std::cerr << "Using tabulated default model" << std::endl;
    boost::shared_ptr<Model> Mod(new TabFunction(parms, name));
    return boost::shared_ptr<DefaultModel>(new GeneralDefaultModel(parms, Mod));
  }
}



#endif
