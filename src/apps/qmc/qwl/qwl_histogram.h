/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2004 by Stefan Wessel <wessel@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef ALPS_QWL_HISTOGRAM_H
#define ALPS_QWL_HISTOGRAM_H

#include <alps/hdf5/archive.hpp>
#include <cmath>
#include <vector>
#include <algorithm>
#include <valarray>

template<typename T=double> 
class histogram {
 public:
  histogram() : data(), left_(0), right_(0) {}
  void resize(unsigned int newsize,unsigned int newleft=0) {
    data.resize(newsize);
    left_=newleft;
    right_=left_+size()-1;
    fill(0);
  }
  void fill(T x) {
    std::fill(data.begin(), data.end(), x);
  }
  void subtract() {
    T x=data[0];
    for (unsigned int i=0;i<data.size();++i)
      data[i]-=x;
  }
  std::valarray<T> getvalarray(unsigned int min,unsigned int max) const {
    std::valarray<T> dval;
    dval.resize(max-min+1);
    int count=0;
    for (int i=min-left();i<=max-left();++i) {
      dval[count]=data[i];
      ++count;
    }
    return dval;
  }
  T& operator[](unsigned int i) {
    return data[i-left()];
  }
  unsigned int left() const {
    return left_;
  }  
  unsigned int right() const {
    return right_;
  }  
  unsigned int size() const {
    return data.size();
  }
  double flatness() const {  
    double av=data[0];
    for (int i=1;i<size();++i)
      av+=data[i];
    av/=size();
    double diff=fabs(data[0]-av);
    for (int i=1;i<size();++i)
      if (diff<fabs(data[i]-av))  
        diff=fabs(data[i]-av);
    return ( av>0 ? diff/av : 0);
  }  
  T min() const {  
   T min=data[0];
    for (int i=1;i<size();++i)
      if (data[i]<min)
         min=data[i];
    return min;
  }  
  void save(alps::hdf5::archive& ar) const { ar["values"] << data; }
  void load(alps::hdf5::archive& ar) {
    std::vector<T> values;
    ar["values"] >> values;
    if (values.size()!=data.size() || !std::all_of(values.begin(),values.end(),
        [](T value) { return std::isfinite(value); }))
      throw std::invalid_argument("Invalid QWL histogram checkpoint");
    data=std::move(values);
  }
 private:
  std::vector<T> data;
  unsigned int left_;
  unsigned int right_;
};

#endif


