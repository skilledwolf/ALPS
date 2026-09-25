// SPDX-License-Identifier: MIT
#define WITH_LAPACK
#include <dmtk/dmtk.h>
#include <array>
#include <numeric>
#include <random>
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

template<class Value, class Reference>
void check_order(Value const& a, Value const& b, Reference const& left, Reference const& right) {
  require((a == b) == (left == right)); require((a != b) == (left != right));
  require((a < b) == (left < right)); require((a <= b) == (left <= right));
  require((a > b) == (left > right)); require((a >= b) == (left >= right));
}

void check_ordering() {
  QN::init();
  for (int i = 0; i < 6; ++i) QN::add_qn_index(std::to_string(i));
  std::mt19937 random(1729);
  for (int mask = 0; mask < 64; ++mask) {
    QN::set_qn_mask(mask);
    for (int sample = 0; sample < 128; ++sample) {
      StateSpace a, b;
      std::vector<int> state_left, state_right;
      for (int site = 1; site <= 4; ++site) {
        QN qa, qb;
        std::vector<int> left, right;
        for (int i = 0; i < 6; ++i) {
          int x = int(random() % 9) - 4, y = int(random() % 9) - 4;
          if (sample % 4 == 0) y = x;
          qa[i].set_half(x); qb[i].set_half(y);
          if (mask & (1 << i)) { left.push_back(x); right.push_back(y); }
        }
        check_order(qa, qb, left, right);
        require(qa.equal(qb, mask) == (left == right));
        // Ranges do not participate in StateSpace ordering.
        a[site] = SubSpace(qa, 0, 2); b[site] = SubSpace(qb, 3, 8);
        state_left.insert(state_left.end(), left.begin(), left.end());
        state_right.insert(state_right.end(), right.begin(), right.end());
      }
      check_order(a, b, state_left, state_right);
    }
  }
}

void check_copy_ranges() {
  std::array<int, 24> source, target;
  std::iota(source.begin(), source.end(), -7);
  for (int count = 0; count <= 24; ++count) {
    target.fill(99);
    array_copy(count, static_cast<int const*>(source.data()), target.data());
    for (int i = 0; i < 24; ++i) require(target[i] == (i < count ? source[i] : 99));
    target.fill(99);
    array_copy(count, source, target);
    for (int i = 0; i < 24; ++i) require(target[i] == (i < count ? source[i] : 99));
    std::array<double, 24> converted;
    converted.fill(99);
    array_copy2(count, static_cast<std::array<int, 24> const&>(source), converted);
    for (int i = 0; i < 24; ++i) require(converted[i] == (i < count ? source[i] : 99));
    converted.fill(99);
    array_copy2(count, source, converted);
    for (int i = 0; i < 24; ++i) require(converted[i] == (i < count ? source[i] : 99));
  }
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

template<int A, class View>
void check_subslice(View& view) {
  coordinates i{1,1,1,1};
  auto sub = one<A>(view, std::slice(1,2,2), i);
  for (size_t n = 0; n < 2; ++n) {
    i[A] = 1+2*n;
    require(sub[n] == view(i[0],i[1],i[2],i[3]));
  }
}

template<int A, int B, class View>
void check_subslice_pair(View& view) {
  coordinates i{1,1,1,1};
  auto sub = two<A,B>(view, std::slice(1,2,2), std::slice(0,2,3), i);
  for (size_t a = 0; a < 2; ++a) for (size_t b = 0; b < 2; ++b) {
    i[A] = 1+2*a; i[B] = 3*b;
    require(sub(a,b) == view(i[0],i[1],i[2],i[3]));
    require(sub.column(a)[b] == sub(a,b));
    require(sub.row(b)[a] == sub(a,b));
    require(*sub.get_pointer(a,b) == sub(a,b));
    require(sub.transpose()(b,a) == sub(a,b));
  }
}

template<class View>
void check_subslices(View& view) {
  require(view.size() == 4*5*6*7);
  check_subslice<0>(view); check_subslice<1>(view);
  check_subslice<2>(view); check_subslice<3>(view);
  check_subslice_pair<0,1>(view); check_subslice_pair<0,2>(view);
  check_subslice_pair<0,3>(view); check_subslice_pair<1,2>(view);
  check_subslice_pair<1,3>(view); check_subslice_pair<2,3>(view);
}

void check_view_storage() {
  std::vector<double> data(4000), target(4000, -1);
  std::iota(data.begin(), data.end(), 0.0);
  state_slice<double> writable(&data, {1,4,1000}, {2,5,100}, {3,6,10}, {4,7,1});
  cstate_slice<double> readable(&data, {1,4,1000}, {2,5,100}, {3,6,10}, {4,7,1});
  check_subslices(writable); check_subslices(readable);
  state_slice<double> destination(&target, {1,4,1000}, {2,5,100}, {3,6,10}, {4,7,1});
  destination = writable;
  destination(1,2,3,4) = -7;
  require(target[1244] == -7 && data[1244] == 1244 && target[0] == -1);

  gslice_iter<double> matrix(&data, {3,3,20}, {5,4,2});
  gslice_iter<double> copy(&target, {1,2,20}, {2,3,2});
  cgslice_iter<double> source(&data, {3,3,20}, {5,4,2});
  copy = source;
  require(target[3] == data[8] && target[27] == data[32]);
  copy(0,0) = -9;
  require(data[8] == 8 && target[3] == -9);
  auto cursor = matrix.begin();
  require((*cursor)[0] == data[8]);
  ++cursor;
  require((*cursor)[1] == data[30]);

  slice_iter<double> line(&data, {3,5,2}), line_copy(&target, {4,3,3});
  line_copy = line;
  require(target[4] == 3 && target[10] == 7);
  auto it = line.begin();
  require(*it++ == 3 && *it == 5 && it.current() == 1);
  require(line.end().current() == 5 && line.begin() < line.end());
  cslice_iter<double> read(&data, {3,5,2}), rebound(&target, {4,3,3});
  rebound = read;
  require(*rebound.get_pointer(2) == 7 && rebound.size() == 5);
  static_assert(!std::is_assignable_v<decltype(read[0]), double>);
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
  for (size_t a=0; a<5; ++a) for (size_t b=0; b<5; ++b)
    for (size_t c=0; c<5; ++c) for (size_t d=0; d<5; ++d)
      if (!constrained || (a>=2)+(b>=2)+(c>=2)+(d>=2) == 2)
        require(static_cast<VectorState<double> const&>(v)(a,b,c,d) == value({a,b,c,d}));
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
  check_view_storage();
  check_ordering();
  check_copy_ranges();
  check_ranges(); check_quantum(false); check_quantum(true);
}
