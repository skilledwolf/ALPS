// SPDX-License-Identifier: MIT
#include "selfenergy.hpp"
#include <map>
#include <stdexcept>

namespace {
void require(bool ok) { if (!ok) throw std::runtime_error("CT-INT measurement contract failed"); }
void near(double actual, double expected) { if (std::abs(actual-expected)>=1e-12) { std::cerr << "actual " << actual << " expected " << expected << "\n"; require(false); } }
struct vertex {
  double time;
  std::vector<std::complex<double>> exponents;
  double t() const { return time; }
  unsigned s() const { return 0; }
  void compute_exp(unsigned n, int direction) {
    exponents.resize(n);
    for (unsigned w=0;w<n;++w)
      exponents[w]=std::exp(std::complex<double>(0,double(direction)*(2*w+1)*std::acos(-1.)*time));
  }
  const std::complex<double>* exp_iomegat() const { return exponents.data(); }
};
struct matrix_state {
  alps::numeric::matrix<double> values{1,1};
  std::vector<vertex> creation{{.25,{}}}, annihilation{{.75,{}}};
  auto& matrix() { return values; }
  auto& creators() { return creation; }
  auto& annihilators() { return annihilation; }
};
struct scalar_sample { double value = 0; void operator<<(double x) { value=x; } };
template<bool SingleSiteName> struct run {
  static constexpr bool ctint_single_site_density = SingleSiteName;
  unsigned n_flavors=2,n_site=1,n_matsubara=4,n_matsubara_measurements=4,n_self=4;
  double beta=1,sign=-1;
  unsigned draws=0;
  bool thermalized=true;
  std::vector<matrix_state> M{2};
  std::map<std::string,scalar_sample> measurements;
  std::map<std::string,std::valarray<double>> vectors;
  run() { for (auto& m:M) m.values(0,0)=1.2; }
  double ctint_uniform() { return .1 + .08*(draws++ % 10); }
  bool is_thermalized() { return thermalized; }
  double green0_spline(double t, unsigned, unsigned, unsigned) const { return t>=0 ? -.5 : .5; }
  void ctint_measure(std::string const& name, std::valarray<double> const& values) { vectors[name]=values; }
};

template<bool Name> void measurements() {
  using kernel=alps::ctint::measurement_kernel;
  run<Name> r;
  kernel::compute_W_matsubara(r);
  require(r.draws==1);
  for (unsigned w=0;w<4;++w) {
    auto expected=-1.2*std::exp(std::complex<double>(0,(2*w+1)*std::acos(-1.)*(-.5)));
    near(r.vectors["Wk_real_0_0_0"][w],expected.real());
    near(r.vectors["Wk_imag_0_0_0"][w],expected.imag());
  }
  // At tau=.1, the Green functions are -.5 and +.5, with inverse matrix 1.2.
  near(r.vectors["densities"][0],-.8);
  near(r.measurements["density_correlation"].value,-.64);
  r.draws=0;
  kernel::compute_W_itime(r);
  require(r.draws==10);
  std::valarray<double> expected(0.,5);
  for (unsigned i=0;i<10;++i) {
    double delta=.25-(.1+.08*i);
    double value=1.2*(delta>=0 ? -.5 : .5);
    if (delta<0) { delta+=1; value=-value; }
    expected[unsigned(delta*4+.5)]+=value*(-.1);
  }
  for(unsigned i=0;i<5;++i) near(r.vectors["W_0_0_0"][i],expected[i]);
  require(r.measurements.count(Name ? "density_0" : "density_0_0")==1);
  r.thermalized=false;r.vectors.clear();r.measurements.clear();
  kernel::compute_W_itime(r);
  require(r.draws==20 && r.vectors.empty() && r.measurements.empty());
}
}
int main() {
  measurements<false>();measurements<true>();
  itime_green_function_t bare(5,1,1), result(5,1,1);
  for (int k=0;k<=4;++k) bare(k,0,0,0)=-.2-.6*k/4;
  near(alps::ctint::green0_spline(bare,.25,0,0,0,4,1),-.35);
  near(alps::ctint::green0_spline(bare,-.25,0,0,0,4,1),.65);
  near(alps::ctint::green0_spline(bare,0,0,0,0,4,1),-.2);
  alps::ctint::evaluate_itime_green(result,bare,1,1,1,4,4,
    [](std::string const&) { return std::make_pair(std::vector<double>(5,0),std::vector<double>(5,0)); });
  for (int k=0;k<=4;++k) near(result(k,0,0,0),bare(k,0,0,0));
  alps::ctint::evaluate_itime_green(result,bare,1,1,1,4,4,
    [](std::string const&) { return std::make_pair(std::vector<double>{0,0,2,0,0},std::vector<double>(5,0)); });
  // A single nonzero convolution bin at tau=.5.
  near(result(0,0,0,0),-.2-2*.5);
  near(result(3,0,0,0),-.65+2*.35);
  near(result(4,0,0,0),-1-result(0,0,0,0));
}
