// SPDX-License-Identifier: MIT
#include <alps/alea.h>
#include <alps/alea/mcdata.hpp>
#include <alps/ngs/mcresult.hpp>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <type_traits>

void close(double actual, double expected) {
  if (std::abs(actual - expected) > 1e-12 || !std::isfinite(actual))
    throw std::runtime_error("Result arithmetic changed a value or uncertainty");
}
void close(std::vector<double> const& actual, std::vector<double> const& expected) {
  if (actual.size() != expected.size()) throw std::runtime_error("Wrong result size");
  for (std::size_t i = 0; i < actual.size(); ++i) close(actual[i], expected[i]);
}
template<class T>
void check(alps::mcresult const& actual, alps::alea::mcdata<T> const& expected) {
  if (!actual.is_type<T>()) throw std::runtime_error("Wrong result type");
  close(actual.mean<T>(), expected.mean());
  close(actual.error<T>(), expected.error());
}
template<class Exception, class F>
void rejects(F operation) {
  try { operation(); } catch (Exception const&) { return; }
  throw std::runtime_error("Invalid result arithmetic was accepted");
}

int main() {
  alps::RealObservable scalar("scalar"), other("other");
  alps::RealVectorObservable vector("vector"), other_vector("other vector");
  for (int i = 0; i < 32; ++i) {
    scalar << 2. + .1*(i % 3);
    other << 4. + .2*(i % 5);
    vector << std::valarray<double>{3. + .1*(i % 3), 7. + .3*(i % 4)};
    other_vector << std::valarray<double>{5. + .2*(i % 5), 8. + .1*(i % 7)};
  }
  alps::mcresult s(&scalar), t(&other), v(&vector), w(&other_vector);
  alps::alea::mcdata<double> sd(scalar), td(other);
  alps::alea::mcdata<std::vector<double>> vd(vector), wd(other_vector);
  std::vector<double> raw{2., 4.};
  auto exercise = [&](auto op, auto assign) {
    check(op(s, t), op(sd, td));
    check(op(v, w), op(vd, wd));
    if constexpr (std::is_same_v<decltype(op), std::plus<>> || std::is_same_v<decltype(op), std::minus<>>) {
      // The numerical backend does not implement mixed-result addition/subtraction.
      rejects<std::logic_error>([&] { op(s, v); });
      rejects<std::logic_error>([&] { op(v, s); });
    } else {
      check(op(s, v), op(sd, vd));
      check(op(v, s), op(vd, sd));
    }
    check(op(s, 2.), op(sd, 2.));
    check(op(2., s), op(2., sd));
    check(op(v, 2.), op(vd, 2.));
    check(op(2., v), op(2., vd));
    check(op(v, raw), op(vd, raw));
    check(op(raw, v), op(raw, vd));
    auto inplace = [&](auto const& obs, auto const& rhs, auto expected, auto const& rhs_data) {
      alps::mcresult actual(&obs);
      assign(actual, rhs);
      assign(expected, rhs_data);
      check(actual, expected);
    };
    inplace(scalar, 2., sd, 2.);
    inplace(scalar, t, sd, td);
    inplace(vector, 2., vd, 2.);
    inplace(vector, raw, vd, raw);
    inplace(vector, w, vd, wd);
    rejects<std::bad_cast>([&] { op(s, raw); });
    rejects<std::bad_cast>([&] { op(raw, s); });
    rejects<std::bad_cast>([&] { assign(s, raw); });
    rejects<std::bad_cast>([&] { assign(s, v); });
    rejects<std::bad_cast>([&] { assign(v, s); });
  };
  exercise(std::plus<>{}, [](auto& a, auto const& b) { a += b; });
  exercise(std::minus<>{}, [](auto& a, auto const& b) { a -= b; });
  exercise(std::multiplies<>{}, [](auto& a, auto const& b) { a *= b; });
  exercise(std::divides<>{}, [](auto& a, auto const& b) { a /= b; });
}
