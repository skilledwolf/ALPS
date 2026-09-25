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

#ifndef __DMTK_SLICE_ITER_H__
#define __DMTK_SLICE_ITER_H__

// slice iterator class. From Stroustrup's book.

#include <iosfwd>
#include <vector>
#include <type_traits>
#include <valarray>
#include "conj.h"
#include "meta.h"
#include "array_util.h"

namespace dmtk
{

namespace detail {
// Share cursor and addressing logic; writable views keep their copy-through
// assignment while readonly views retain ordinary cursor assignment.
template<class T, class View>
class slice_cursor {
public:
  using value_type = std::remove_const_t<T>;
  using _V = std::vector<value_type>;
  using storage = std::conditional_t<std::is_const_v<T>, const _V, _V>;
  using reference = std::conditional_t<std::is_const_v<T>, value_type, T&>;
  using const_reference = std::conditional_t<std::is_const_v<T>, value_type, const T&>;

  slice_cursor() = default;
  slice_cursor(storage* vv, std::slice ss): v(vv), s(ss) {}
  size_t size() const { return s.size(); }
  size_t size1() const { return 1; }
  size_t size2() const { return size(); }
  size_t current() const { return curr; }
  View begin() const { View t = view(); t.curr = 0; return t; }
  View end() const { View t = view(); t.curr = size(); return t; }
  View& operator++() { ++curr; return view(); }
  View operator++(int) { View t = view(); ++curr; return t; }
  View& operator--() { --curr; return view(); }
  View operator--(int) { View t = view(); --curr; return t; }
  reference operator[](size_t i) { return ref(curr = i); }
  const_reference operator[](size_t i) const { return ref(i); }
  reference operator()(size_t i) { return ref(curr = i); }
  const_reference operator()(size_t i) const { return ref(i); }
  reference operator*() const { return ref(curr); }
  const value_type* get_pointer(size_t i) const { return &(*v)[s.start()+i*s.stride()]; }
  bool operator==(const View& q) const { return same_slice(q) && curr == q.curr; }
  bool operator!=(const View& q) const { return !(*this == q); }
  bool operator<(const View& q) const { return same_slice(q) && curr < q.curr; }
  bool operator>(const View& q) const { return same_slice(q) && curr > q.curr; }

protected:
  storage* v = nullptr;
  std::slice s{0, 0, 0};
  size_t curr = 0;
  reference ref(size_t i) const { return (*v)[s.start()+i*s.stride()]; }

private:
  View& view() { return static_cast<View&>(*this); }
  const View& view() const { return static_cast<const View&>(*this); }
  bool same_slice(const View& q) const {
    return s.stride() == q.s.stride() && s.start() == q.s.start();
  }
};
} // namespace detail

template<class T>
class cslice_iter : public detail::slice_cursor<const T, cslice_iter<T>> {
  using base = detail::slice_cursor<const T, cslice_iter<T>>;
public:
  using base::base;
};

template<class T>
class slice_iter : public detail::slice_cursor<T, slice_iter<T>> {
  using base = detail::slice_cursor<T, slice_iter<T>>;
  friend base;
  using base::curr;
  using base::ref;
public:
  using typename base::_V;
  using base::base;
  using base::size;
    slice_iter& operator=(slice_iter<T> ss)
      { 
        int n = std:: min(ss.size(), size());
        array_copy(n, ss, *this);
//        for(int i = 0; i < n; i++) ref(i) = ss[i];
        return *this;
      }
    slice_iter& operator=(cslice_iter<T> ss)
      { 
        int n = std:: min(ss.size(), size());
        array_copy2(n, ss, *this);
//        for(int i = 0; i < n; i++) ref(i) = ss[i];
        return *this;
      }

    slice_iter& operator=(const T& v){ref(curr) = v; return *this;}

    template<class Expr>
    slice_iter& operator=(const IterExpr<T,Expr>&);

    slice_iter& operator+=(slice_iter<T>);
    slice_iter& operator-=(slice_iter<T>);
    slice_iter& operator*=(slice_iter<T>);
    slice_iter& operator/=(slice_iter<T>);

    slice_iter& operator=(_V);
    slice_iter& operator+=(_V);
    slice_iter& operator-=(_V);
    slice_iter& operator*=(_V);
    slice_iter& operator/=(_V);

    slice_iter& operator=(const T*);
    slice_iter& operator+=(const T*);
    slice_iter& operator-=(const T*);
    slice_iter& operator*=(const T*);
    slice_iter& operator/=(const T*);

    template<class Expr>
    slice_iter& operator+=(const IterExpr<T,Expr>&);
    template<class Expr>
    slice_iter& operator-=(const IterExpr<T,Expr>&);
    template<class Expr>
    slice_iter& operator*=(const IterExpr<T,Expr>&);
    template<class Expr>
    slice_iter& operator/=(const IterExpr<T,Expr>&);

};


#include "slice_implement.h"

template<class T>
template<class Expr>
inline slice_iter<T>&
slice_iter<T>::operator=(const IterExpr<T,Expr>& vv)
{
  curr = 0;
  array_copy2(vv.size2(), vv, *this);
//  for(uint i = 0; i < size(); i++) { *this->operator++(0) = vv[i]; }
  curr = 0;
  return *this;
}

#define BINARY_OP(op,ap) \
template<class T> \
template<class Expr> \
inline slice_iter<T>& \
slice_iter<T>::op(const IterExpr<T,Expr>& vv) \
{ \
  curr = 0; \
  for(int i = 0; i < size(); i++) { *this->operator++(0) ap vv[i]; } \
  curr = 0; \
  return *this; \
} 

BINARY_OP(operator+=,+=);
BINARY_OP(operator-=,-=);
BINARY_OP(operator/=,/=);
BINARY_OP(operator*=,*=);
#undef BINARY_OP

#define BINARY_OP(op,ap) \
template<class T> \
inline slice_iter<T>& \
slice_iter<T>::op(slice_iter<T> vv) \
{ \
  curr = 0; \
  for(int i = 0; i < size(); i++) { *this->operator++(0) ap vv[i]; } \
  curr = 0; \
  return *this; \
} 

BINARY_OP(operator+=,+=);
BINARY_OP(operator-=,-=);
BINARY_OP(operator/=,/=);
BINARY_OP(operator*=,*=);
#undef BINARY_OP

#define BINARY_OP(op,ap) \
template<class T> \
inline slice_iter<T>& \
slice_iter<T>::op(_V vv) \
{ \
  curr = 0; \
  int n = size(); \
  T* _v = &vv[0]; \
  while(n--) { *this->operator++(0) ap *_v++; } \
  curr = 0; \
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
inline slice_iter<T>& \
slice_iter<T>::op(const T* _v) \
{ \
  curr = 0; \
  int n = size(); \
  while(n--) { *this->operator++(0) ap *_v++; } \
  curr = 0; \
  return *this; \
} 

BINARY_OP(operator=,=);
BINARY_OP(operator+=,+=);
BINARY_OP(operator-=,-=);
BINARY_OP(operator/=,/=);
BINARY_OP(operator*=,*=);
#undef BINARY_OP

template<class T>
inline T
product(cslice_iter<T> a, cslice_iter<T> b)
{
  if(a.size() != b.size())
     cerr << "** Warning: sizes do not comform\n";

  T r = T(0);

  for(size_t i = 0; i < std::min(a.size(), b.size()); i++) 
    r += std::conj(a[i])*b[i]; 

  return r;
}

template<class T>
inline T
product(slice_iter<T> a, slice_iter<T> b)
{
  if(a.size() != b.size())
     cerr << "** Warning: sizes do not comform\n";

  T r = T(0);

  for(size_t i = 0; i < std::min(a.size(), b.size()); i++) 
    r += std::conj(a[i])*b[i]; 

  return r;
}

} // namespace dmtk 

#endif // __DMTK_SLICE_ITER_H__

