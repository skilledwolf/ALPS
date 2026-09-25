/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef __DMTK_GSLICE_ITER_H__
#define __DMTK_GSLICE_ITER_H__

// gslice iterator class. From Stroustrup's book.

#include <vector>
#include "slice_iter.h" 

namespace dmtk
{

#include "meta.h"

namespace detail {
template<class T, class View>
class gslice_cursor {
public:
  using value_type = std::remove_const_t<T>;
  using _V = std::vector<value_type>;
  using storage = std::conditional_t<std::is_const_v<T>, const _V, _V>;
  using reference = std::conditional_t<std::is_const_v<T>, value_type, T&>;
  using row_view = std::conditional_t<std::is_const_v<T>, cslice_iter<value_type>, slice_iter<value_type>>;

  gslice_cursor() = default;
  gslice_cursor(storage* vv, std::slice first, std::slice second): v(vv), s1(first), s2(second) {}
  View begin() const { View t = view(); t.curr1 = t.curr2 = 0; return t; }
  View end() const { View t = view(); t.curr1 = s1.size(); t.curr2 = s2.size(); return t; }
  View& operator++() { ++curr1; return view(); }
  View operator++(int) { View t = view(); ++curr1; return t; }
  row_view operator[](size_t i) { return column(curr1 = i); }
  reference operator()(size_t i, size_t j) { curr1 = i; curr2 = j; return (*v)[index(i,j)]; }
  row_view operator*() { return column(curr1); }
  row_view column(size_t i) const {
    return row_view(v, std::slice(index(i,0), s2.size(), s2.stride()));
  }
  row_view row(size_t j) const {
    return row_view(v, std::slice(index(0,j), s1.size(), s1.stride()));
  }
  const value_type* get_pointer(size_t i, size_t j) const { return &(*v)[index(i,j)]; }
  size_t size() const { return s1.size(); }
  size_t size1() const { return s1.size(); }
  size_t size2() const { return s2.size(); }
  size_t stride1() const { return s1.stride(); }
  size_t stride2() const { return s2.stride(); }
  size_t start1() const { return s1.start(); }
  size_t start2() const { return s2.start(); }
  View transpose() const { return View(v, s2, s1); }

protected:
  storage* v = nullptr;
  std::slice s1{0,0,0}, s2{0,0,0};
  size_t curr1 = 0, curr2 = 0;

private:
  View& view() { return static_cast<View&>(*this); }
  const View& view() const { return static_cast<const View&>(*this); }
  size_t index(size_t i, size_t j) const {
    return s1.start()+i*s1.stride()+s2.start()+j*s2.stride();
  }
};
} // namespace detail

template<class T>
class cgslice_iter : public detail::gslice_cursor<const T, cgslice_iter<T>> {
  using base = detail::gslice_cursor<const T, cgslice_iter<T>>;
public:
  using base::base;
};

template<class T>
class gslice_iter : public detail::gslice_cursor<T, gslice_iter<T>> {
  using base = detail::gslice_cursor<T, gslice_iter<T>>;
  using base::s1;
  using base::s2;
public:
  using typename base::_V;
  using base::base;
  using base::operator();
    gslice_iter& operator=(gslice_iter<T> ss)
      { 
        size_t n1 = std::min(ss.s1.size(), s1.size());
        size_t n2 = std::min(ss.s2.size(), s2.size());
        for(size_t i = 0; i < n1; i++)
          for(size_t j = 0; j < n2; j++)
            this->operator()(i,j) = ss(i,j);
        return *this; 
      }

    gslice_iter& operator=(cgslice_iter<T> ss)
      { 
        size_t n1 = std::min(ss.size1(), s1.size());
        size_t n2 = std::min(ss.size2(), s2.size());
        for(size_t i = 0; i < n1; i++)
          for(size_t j = 0; j < n2; j++)
            this->operator()(i,j) = ss(i,j);
        return *this; 
      }

//  Operators

    template<class Expr>
    gslice_iter& operator=(const IterExpr<T,Expr>&);

    gslice_iter& operator+=(gslice_iter<T>);
    gslice_iter& operator-=(gslice_iter<T>);
    gslice_iter& operator*=(gslice_iter<T>);
    gslice_iter& operator/=(gslice_iter<T>);

    gslice_iter& operator+=(cgslice_iter<T>);
    gslice_iter& operator-=(cgslice_iter<T>);
    gslice_iter& operator*=(cgslice_iter<T>);
    gslice_iter& operator/=(cgslice_iter<T>);

    gslice_iter& operator=(_V);
    gslice_iter& operator+=(_V);
    gslice_iter& operator-=(_V);
    gslice_iter& operator*=(_V);
    gslice_iter& operator/=(_V);

    gslice_iter& operator=(const T*);
    gslice_iter& operator+=(const T*);
    gslice_iter& operator-=(const T*);
    gslice_iter& operator*=(const T*);
    gslice_iter& operator/=(const T*);

    template<class Expr>
    gslice_iter& operator+=(const IterExpr<T,Expr>&);
    template<class Expr>
    gslice_iter& operator-=(const IterExpr<T,Expr>&);
    template<class Expr>
    gslice_iter& operator*=(const IterExpr<T,Expr>&);
    template<class Expr>
    gslice_iter& operator/=(const IterExpr<T,Expr>&);

};


#include "gslice_implement.h"

template<class T>
template<class Expr>
inline gslice_iter<T>&
gslice_iter<T>::operator=(const IterExpr<T,Expr>& vv)
{
  for(size_t i = 0; i < s1.size(); i++)
    for(size_t j = 0; j < s2.size(); j++) 
      this->operator()(i,j) = vv(i,j); 

  return *this;
}

//////////////////////////////////////////////////////////

#define BINARY_OP(op,ap) \
template<class T> \
template<class Expr> \
inline gslice_iter<T>& \
gslice_iter<T>::op(const IterExpr<T,Expr>& mm) \
{ \
  for(size_t i = 0; i < s1.size(); i++) \
    for(size_t j = 0; j < s2.size(); j++){ \
      operator()(i,j) ap mm(i,j); \
  } \
 \
  return *this; \
}

BINARY_OP(operator+=,+=);
BINARY_OP(operator-=,-=);
BINARY_OP(operator/=,/=);
BINARY_OP(operator*=,*=);
#undef BINARY_OP

#define BINARY_OP(op,ap) \
template<class T> \
inline gslice_iter<T>& \
gslice_iter<T>::op(gslice_iter<T> mm) \
{ \
  for(size_t i = 0; i < s1.size(); i++) \
    for(size_t j = 0; j < s2.size(); j++){ \
      operator()(i,j) ap mm(i,j); \
  } \
 \
  return *this; \
}

BINARY_OP(operator+=,+=);
BINARY_OP(operator-=,-=);
BINARY_OP(operator/=,/=);
BINARY_OP(operator*=,*=);
#undef BINARY_OP

#define BINARY_OP(op,ap) \
template<class T> \
inline gslice_iter<T>& \
gslice_iter<T>::op(cgslice_iter<T> mm) \
{ \
  for(size_t i = 0; i < s1.size(); i++) \
    for(size_t j = 0; j < s2.size(); j++){ \
      operator()(i,j) ap mm(i,j); \
  } \
 \
  return *this; \
}

BINARY_OP(operator+=,+=);
BINARY_OP(operator-=,-=);
BINARY_OP(operator/=,/=);
BINARY_OP(operator*=,*=);
#undef BINARY_OP

#define BINARY_OP(op,ap) \
template<class T> \
inline gslice_iter<T>& \
gslice_iter<T>::op(_V v) \
{ \
  T* _v = &v[0]; \
  for(size_t i = 0; i < s1.size(); i++) \
    for(size_t j = 0; j < s2.size(); j++){ \
      operator()(i,j) ap *_v++; \
  } \
 \
  return *this; \
}

BINARY_OP(operator=,=);
BINARY_OP(operator+=,+=);
BINARY_OP(operator-=,-=);
BINARY_OP(operator/=,/=);
BINARY_OP(operator*=,*=);
#undef BINARY_OP

#define BINARY_OP(op,ap) \
template<class T> \
inline gslice_iter<T>& \
gslice_iter<T>::op(const T* _v) \
{ \
  for(size_t i = 0; i < s1.size(); i++) \
    for(size_t j = 0; j < s2.size(); j++){ \
      operator()(i,j) ap *_v++; \
  } \
 \
  return *this; \
}

BINARY_OP(operator=,=);
BINARY_OP(operator+=,+=);
BINARY_OP(operator-=,-=);
BINARY_OP(operator/=,/=);
BINARY_OP(operator*=,*=);
#undef BINARY_OP
} // namespace dmtk 

#endif // __DMTK_GSLICE_ITER_H__


