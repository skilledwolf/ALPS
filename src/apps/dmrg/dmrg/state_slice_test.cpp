// SPDX-License-Identifier: MIT
#define WITH_LAPACK
#include <dmtk/dmtk.h>
#include <array>
#include <stdexcept>

namespace {
using namespace dmtk;
using coordinates = std::array<size_t, 4>;

double value(coordinates const& i) {
  return 1000*i[0] + 100*i[1] + 10*i[2] + i[3];
}
void require(bool ok) {
  if (!ok) throw std::runtime_error("state slice differs from coordinate reference");
}

template<int A, class State, class Selector>
auto one(State& v, Selector q, coordinates i) {
  if constexpr (A == 0) return v(q, i[1], i[2], i[3]);
  if constexpr (A == 1) return v(i[0], q, i[2], i[3]);
  if constexpr (A == 2) return v(i[0], i[1], q, i[3]);
  if constexpr (A == 3) return v(i[0], i[1], i[2], q);
}

template<int A, int B, class State, class Selector>
auto two(State& v, Selector q, Selector r, coordinates i) {
  if constexpr (A == 0 && B == 1) return v(q, r, i[2], i[3]);
  if constexpr (A == 0 && B == 2) return v(q, i[1], r, i[3]);
  if constexpr (A == 0 && B == 3) return v(q, i[1], i[2], r);
  if constexpr (A == 1 && B == 2) return v(i[0], q, r, i[3]);
  if constexpr (A == 1 && B == 3) return v(i[0], q, i[2], r);
  if constexpr (A == 2 && B == 3) return v(i[0], i[1], q, r);
}

template<int A, class Selector>
void check_one(VectorState<double>& v, Selector q, Range range, coordinates i) {
  auto writable = one<A>(v, q, i);
  auto readable = one<A>(static_cast<VectorState<double> const&>(v), q, i);
  require(writable.size() == range.size() && readable.size() == range.size());
  for (size_t a = 0; a < range.size(); ++a) {
    i[A] = range.start() + a*range.stride();
    require(writable[a] == value(i) && readable[a] == value(i));
    writable[a] += 0.5;
    require(v(i[0], i[1], i[2], i[3]) == value(i)+0.5);
    writable[a] -= 0.5;
  }
}

template<int A, int B, class Selector>
void check_two(VectorState<double>& v, Selector q, Selector r,
               Range first, Range second, coordinates i) {
  auto writable = two<A,B>(v, q, r, i);
  auto readable = two<A,B>(static_cast<VectorState<double> const&>(v), q, r, i);
  require(writable.size1() == first.size() && writable.size2() == second.size());
  require(readable.size1() == first.size() && readable.size2() == second.size());
  for (size_t a = 0; a < first.size(); ++a)
    for (size_t b = 0; b < second.size(); ++b) {
      i[A] = first.start() + a*first.stride();
      i[B] = second.start() + b*second.stride();
      require(writable(a,b) == value(i) && readable(a,b) == value(i));
      writable(a,b) += 0.5;
      require(v(i[0], i[1], i[2], i[3]) == value(i)+0.5);
      writable(a,b) -= 0.5;
    }
}

void check_ranges() {
  VectorState<double> v(6,7,8,9);
  for (size_t a=0; a<6; ++a) for (size_t b=0; b<7; ++b)
    for (size_t c=0; c<8; ++c) for (size_t d=0; d<9; ++d)
      v(a,b,c,d) = value({a,b,c,d});
  coordinates i{2,3,4,5};
  Range r(1,5,2), s(2,5);
  check_one<0>(v,r,r,i); check_one<1>(v,r,r,i);
  check_one<2>(v,r,r,i); check_one<3>(v,r,r,i);
  check_two<0,1>(v,r,s,r,s,i); check_two<0,2>(v,r,s,r,s,i);
  check_two<0,3>(v,r,s,r,s,i); check_two<1,2>(v,r,s,r,s,i);
  check_two<1,3>(v,r,s,r,s,i); check_two<2,3>(v,r,s,r,s,i);
  auto view = v(i[0], r, s, i[3]);
  require(view.start1() == i[0]*7*8*9 + r.start()*8*9);
  require(view.start2() == s.start()*9 + i[3]);
  require(view.stride1() == r.stride()*8*9 && view.stride2() == 9);
}

void check_quantum(bool constrained) {
  QN::init(); QN::add_qn_index("N", true); QN::set_qn_mask(1);
  PackedBasis basis;
  basis.push_back(SubSpace(QN(0),0,1));
  basis.push_back(SubSpace(QN(1),2,4));
  VectorState<double> v = constrained
    ? VectorState<double>(basis,basis,basis,basis,QN(2),1)
    : VectorState<double>(basis,basis,basis,basis);
  for (size_t a=0; a<5; ++a) for (size_t b=0; b<5; ++b)
    for (size_t c=0; c<5; ++c) for (size_t d=0; d<5; ++d)
      if (!constrained || (a>=2)+(b>=2)+(c>=2)+(d>=2) == 2)
        v(a,b,c,d) = value({a,b,c,d});
  for (int bits : {3,5,6,9,10,12}) {
    std::array<QN,4> q;
    std::array<Range,4> r{Range(0,1),Range(0,1),Range(0,1),Range(0,1)};
    coordinates i;
    for (int axis=0; axis<4; ++axis) {
      q[axis] = QN((bits>>axis)&1);
      r[axis] = basis(q[axis]);
      i[axis] = r[axis].start()+1;
    }
    check_one<0>(v,q[0],r[0],i); check_one<1>(v,q[1],r[1],i);
    check_one<2>(v,q[2],r[2],i); check_one<3>(v,q[3],r[3],i);
    check_two<0,1>(v,q[0],q[1],r[0],r[1],i);
    check_two<1,2>(v,q[1],q[2],r[1],r[2],i);
    check_two<2,3>(v,q[2],q[3],r[2],r[3],i);
  }
}
}

int main() {
  check_ranges(); check_quantum(false); check_quantum(true);
}
