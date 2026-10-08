/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2009 by Matthias Troyer <troyer@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_MODEL_MEASUREMENT_OPERATORS_H
#define ALPS_MODEL_MEASUREMENT_OPERATORS_H

#include <alps/hdf5.hpp>
#include <alps/hdf5/vector.hpp>

#include <alps/parameter.h>
#include <alps/xml.h>
#include <alps/config.h>
#include <alps/utility/encode.hpp>
#include <alps/utility/numeric_cast.hpp>
#include <alps/numeric/real.hpp>
#include <boost/foreach.hpp>
#include <vector>
#include <string>
#include <map> 

namespace alps{

class ALPS_DECL MeasurementOperators
{
public:
  MeasurementOperators (Parameters const& p);
  
  bool calc_averages() const 
  { 
    return !(average_expressions.empty() && local_expressions.empty()
              && correlation_expressions.empty() && structurefactor_expressions.empty() );
  }

protected:
  bool calc_labels() const 
  { 
    return !(local_expressions.empty()
              && correlation_expressions.empty() && structurefactor_expressions.empty() );
  }
  
  std::map<std::string,std::string> average_expressions;
  std::map<std::string,std::string> local_expressions;
  std::map<std::string,std::pair<std::string,std::string> > correlation_expressions;
  std::map<std::string,std::pair<std::string,std::string> > structurefactor_expressions;
};

class ALPS_DECL MeasurementLabels : public MeasurementOperators
{
public:

  template <class LatticeModel>
  MeasurementLabels(LatticeModel const&, int);

protected:
  std::vector<std::string> distlabel_;
  std::vector<std::string> momentumlabel_;
  std::vector<std::string> bondlabel_;
  std::vector<std::string> sitelabel_;
  mutable std::map<std::string,bool> bond_operator_; // mutable to allow operator[]
};


template <class ValueType>
class EigenvectorMeasurements 
 : public MeasurementLabels
{
public:
  typedef ValueType value_type;
  
  template <class LatticeModel>
  EigenvectorMeasurements(LatticeModel const&);
  virtual ~EigenvectorMeasurements() {}

  
  virtual void save(hdf5::archive &) const;
  virtual void load(hdf5::archive &);

  std::map<std::string,std::vector<value_type> > average_values;
  std::map<std::string,std::vector<std::vector<value_type> > > local_values;
  std::map<std::string,std::vector<std::vector<value_type> > > correlation_values;
  std::map<std::string,std::vector<std::vector<value_type> > > structurefactor_values;
  
  bool empty() const { return average_values.empty() && local_values.empty() && correlation_values.empty() && structurefactor_values.empty();}
};


template<class LatticeModel>
MeasurementLabels::MeasurementLabels(LatticeModel const& lattice_model, int /* unused*/)
 : MeasurementOperators(lattice_model.get_parameters())
{
  if (calc_labels()) 
  {
    distlabel_ = lattice_model.distance_labels();
    momentumlabel_ = lattice_model.momenta_labels();
    bondlabel_ = lattice_model.bond_labels();
    sitelabel_ = lattice_model.site_labels();
    typedef std::pair<std::string,std::string> string_pair;
    BOOST_FOREACH(string_pair const& x, local_expressions)
      bond_operator_[x.first] = lattice_model.has_bond_operator(x.second);
  }
}


template <class ValueType>
template<class LatticeModel>
EigenvectorMeasurements<ValueType>::EigenvectorMeasurements(LatticeModel const& lattice_model)
 : MeasurementLabels(lattice_model,0)
{
}

template <class ValueType>
void EigenvectorMeasurements<ValueType>::save(alps::hdf5::archive& ar) const
{
  using alps::numeric::real;

  for (typename std::map<std::string,std::vector<value_type> >::const_iterator
    it=average_values.begin();it!=average_values.end();++it) 
  {
    std::string path = "results/"+ hdf5_name_encode(it->first);
    ar << make_pvp(path+"/mean/value", real(it->second));
  }    

  for (typename std::map<std::string,std::vector<std::vector<value_type> > >::const_iterator 
          it=local_values.begin();it!=local_values.end();++it) {
    std::string path = "results/"+ hdf5_name_encode(it->first);
    ar << make_pvp(path+"/mean/value", real(it->second));
    if (bond_operator_[it->first])
      ar << make_pvp(path+"/labels", bondlabel_);
    else
      ar << make_pvp(path+"/labels", sitelabel_);
  }

  for (typename std::map<std::string,std::vector<std::vector<value_type> > >::const_iterator 
        it=correlation_values.begin();it!=correlation_values.end();++it) {
    std::string path = "results/"+ hdf5_name_encode(it->first);
    ar << make_pvp(path+"/mean/value", real(it->second));
    ar << make_pvp(path+"/labels", distlabel_);
  }

  for (typename std::map<std::string,std::vector<std::vector<value_type> > >::const_iterator 
        it=structurefactor_values.begin();it!=structurefactor_values.end();++it) {
    std::string path = "results/"+ hdf5_name_encode(it->first);
    ar << make_pvp(path+"/mean/value", real(it->second));
    ar << make_pvp(path+"/labels", momentumlabel_);
  }

}

template <class ValueType>
void EigenvectorMeasurements<ValueType>::load(alps::hdf5::archive & ar)
{
  std::vector<std::string> list = ar.list_children(ar.get_context()+"/results");
  for (std::vector<std::string>::const_iterator it = list.begin(); it != list.end(); ++it) {
    std::string name = hdf5_name_decode(*it);
    std::string path = "results/"+*it;
    if (average_expressions.find(name) != average_expressions.end() || name == "Energy") {
      std::vector<double> vals;
      ar >> make_pvp(path+"/mean/value",vals);
      average_values[name] = numeric_cast<std::vector<ValueType> > (vals);
    }
    else {
      std::vector<std::vector<double> > vals;
      std::vector<std::vector<ValueType> >  converted = numeric_cast<std::vector<std::vector<ValueType> > > (vals);
      ar >> make_pvp(path+"/mean/value", vals);
      std::vector<std::string> labels;
      if (local_expressions.find(name) != local_expressions.end()) {
        if (bond_operator_[local_expressions.find(name)->first]) {
          if (bondlabel_.empty() && ar.is_data(path+"/labels"))
            ar >> make_pvp(path+"/labels", bondlabel_);
        }
        else {
          if (sitelabel_.empty() && ar.is_data(path+"/labels"))
            ar >> make_pvp(path+"/labels", sitelabel_);
        }
        local_values[name]=converted;
      }
      else if (correlation_expressions.find(name) != correlation_expressions.end()) {
        if (distlabel_.empty() && ar.is_data(path+"/labels"))
          ar >> make_pvp(path+"/labels", distlabel_);
        correlation_values[name]=converted;
      }
      else if (structurefactor_expressions.find(name) != structurefactor_expressions.end()) {
        if (momentumlabel_.empty() && ar.is_data(path+"/labels"))
          ar >> make_pvp(path+"/labels", momentumlabel_);
        structurefactor_values[name]=converted;
      }
      else
        boost::throw_exception(std::runtime_error("cannot decide whether " + name + " is local or correlation measurement "));
    }
  }
}


}

#endif
