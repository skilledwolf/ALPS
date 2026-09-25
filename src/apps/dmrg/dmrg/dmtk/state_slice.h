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

#ifndef __DMTK_STATE_SLICE_H__
#define __DMTK_STATE_SLICE_H__

#include <vector>
#include <array>
#include "slice_iter.h" 
#include "gslice_iter.h" 

namespace dmtk
{

namespace detail {
template<class T>
class state_slice_view {
public:
  using value_type = std::remove_const_t<T>;
  using storage = std::conditional_t<std::is_const_v<T>, const std::vector<value_type>, std::vector<value_type>>;
  using reference = std::conditional_t<std::is_const_v<T>, value_type, T&>;
  using single_view = std::conditional_t<std::is_const_v<T>, cslice_iter<value_type>, slice_iter<value_type>>;
  using pair_view = std::conditional_t<std::is_const_v<T>, cgslice_iter<value_type>, gslice_iter<value_type>>;

  state_slice_view() = default;
  state_slice_view(storage* vv, std::slice s1, std::slice s2, std::slice s3, std::slice s4)
    : v(vv), slices{s1,s2,s3,s4} {}

  single_view operator()(std::slice r1, size_t i2, size_t i3, size_t i4) const {
    return single({0, i2, i3, i4}, 0, r1);
  }
  single_view operator()(size_t i1, std::slice r2, size_t i3, size_t i4) const {
    return single({i1, 0, i3, i4}, 1, r2);
  }
  single_view operator()(size_t i1, size_t i2, std::slice r3, size_t i4) const {
    return single({i1, i2, 0, i4}, 2, r3);
  }
  single_view operator()(size_t i1, size_t i2, size_t i3, std::slice r4) const {
    return single({i1, i2, i3, 0}, 3, r4);
  }
  pair_view operator()(std::slice r1, std::slice r2, size_t i3, size_t i4) const {
    return pair({0, 0, i3, i4}, 0, 1, r1, r2);
  }
  pair_view operator()(std::slice r1, size_t i2, std::slice r3, size_t i4) const {
    return pair({0, i2, 0, i4}, 0, 2, r1, r3);
  }
  pair_view operator()(std::slice r1, size_t i2, size_t i3, std::slice r4) const {
    return pair({0, i2, i3, 0}, 0, 3, r1, r4);
  }
  pair_view operator()(size_t i1, std::slice r2, std::slice r3, size_t i4) const {
    return pair({i1, 0, 0, i4}, 1, 2, r2, r3);
  }
  pair_view operator()(size_t i1, std::slice r2, size_t i3, std::slice r4) const {
    return pair({i1, 0, i3, 0}, 1, 3, r2, r4);
  }
  pair_view operator()(size_t i1, size_t i2, std::slice r3, std::slice r4) const {
    return pair({i1, i2, 0, 0}, 2, 3, r3, r4);
  }

  size_t start1() const { return slices[0].start(); }
  size_t start2() const { return slices[1].start(); }
  size_t start3() const { return slices[2].start(); }
  size_t start4() const { return slices[3].start(); }
  size_t stride1() const { return slices[0].stride(); }
  size_t stride2() const { return slices[1].stride(); }
  size_t stride3() const { return slices[2].stride(); }
  size_t stride4() const { return slices[3].stride(); }
  size_t size1() const { return slices[0].size(); }
  size_t size2() const { return slices[1].size(); }
  size_t size3() const { return slices[2].size(); }
  size_t size4() const { return slices[3].size(); }
  size_t size() const { return size1()*size2()*size3()*size4(); }
  reference operator()(size_t i, size_t j, size_t k, size_t l) const { return (*v)[index(i,j,k,l)]; }
  reference operator()(size_t i) const { return (*v)[i]; }
  size_t index(size_t i, size_t j, size_t k, size_t l) const { return offset({i,j,k,l}); }

private:
  storage* v = nullptr;
  std::array<std::slice,4> slices{};

  size_t offset(std::array<size_t,4> const& indices) const {
    size_t result = 0;
    for (size_t axis = 0; axis < 4; ++axis)
      result += slices[axis].start() + indices[axis]*slices[axis].stride();
    return result;
  }
  single_view single(std::array<size_t,4> indices, size_t axis, std::slice range) const {
    indices[axis] = range.start();
    return single_view(v, std::slice(offset(indices), range.size(), range.stride()*slices[axis].stride()));
  }
  pair_view pair(std::array<size_t,4> indices, size_t first, size_t second,
                 std::slice a, std::slice b) const {
    indices[first] = a.start();
    indices[second] = b.start();
    size_t second_start = slices[second].start() + b.start()*slices[second].stride();
    // Keep the second axis's offset separate, as in the public layout accessors.
    return pair_view(v, std::slice(offset(indices)-second_start, a.size(), a.stride()*slices[first].stride()),
                        std::slice(second_start, b.size(), b.stride()*slices[second].stride()));
  }
};
} // namespace detail

template<class T>
using cstate_slice = detail::state_slice_view<const T>;

template<class T>
class state_slice : public detail::state_slice_view<T> {
  using base = detail::state_slice_view<T>;
public:
  using base::base;
  state_slice& operator=(state_slice other) {
    for (size_t i = 0; i < std::min(this->size1(), other.size1()); ++i)
      for (size_t j = 0; j < std::min(this->size2(), other.size2()); ++j)
        for (size_t k = 0; k < std::min(this->size3(), other.size3()); ++k)
          for (size_t l = 0; l < std::min(this->size4(), other.size4()); ++l)
            (*this)(i,j,k,l) = other(i,j,k,l);
    return *this;
  }
};

} // namespace dmtk
#endif // __DMTK_STATE_SLICE_H__
