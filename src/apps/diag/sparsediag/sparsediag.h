/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1994-2005 by Matthias Troyer <troyer@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include "../diag.h"
#include <alps/numeric/real.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <ietl/interface/ublas.h>
#include <ietl/vectorspace.h>
#include <ietl/lanczos.h>
#include <boost/random.hpp>
#include <complex>

template <class T>
class SparseDiagMatrix : public DiagMatrix<T, boost::numeric::ublas::mapped_vector_of_mapped_vector<T, boost::numeric::ublas::row_major> >
{
public:
  typedef DiagMatrix<T, boost::numeric::ublas::mapped_vector_of_mapped_vector<T, boost::numeric::ublas::row_major> > super_type;
  typedef T value_type;
  typedef typename super_type::magnitude_type magnitude_type;
  typedef typename super_type::matrix_type matrix_type;
  typedef typename super_type::site_iterator site_iterator;
  typedef typename super_type::site_descriptor site_descriptor;
  typedef typename boost::numeric::ublas::vector<value_type> vector_type;
  typedef typename super_type::size_type size_type;
  typedef typename super_type::mag_vector_type mag_vector_type;
  typedef typename super_type::half_integer_type half_integer_type;
  typedef typename super_type::operator_matrix_type operator_matrix_type;
  
  explicit SparseDiagMatrix (alps::Parameters const& p) : super_type(p) {}
  void do_subspace();
  void print_eigenvectors(std::ostream& os) const;
private:
  
  std::vector<value_type> calculate(operator_matrix_type const& m) const;
  std::vector<vector_type> eigenvectors;
};

template <class T>
void SparseDiagMatrix<T>::do_subspace()
{
  typedef ietl::vectorspace<vector_type> vectorspace_type;
  using alps::numeric::real;
  boost::lagged_fibonacci607 generator;
  mag_vector_type ev;
  this->build();
  if (this->dimension()==0)
    return;

  vectorspace_type vec(this->dimension());
  ietl::lanczos<matrix_type,vectorspace_type> lanczos(this->matrix(),vec);
  int n;
  
  if (this->dimension()>1) {
    int max_iter = this->get_parameters().value_or_default("MAX_ITERATIONS",std::min(int(10*this->dimension()),1000));  
    int num_eigenvalues = this->get_parameters().value_or_default("NUMBER_EIGENVALUES",1);
    ietl::lanczos_iteration_nlowest<double> iter(max_iter,num_eigenvalues);
    std::cerr << "Starting Lanczos \n";
    lanczos.calculate_eigenvalues(iter,generator);
    std::cerr << "Finished Lanczos\n";
    n=std::min(num_eigenvalues,int(lanczos.eigenvalues().size()));
    ev.resize(n);
    for (int i=0;i<n;++i) 
      ev[i]=lanczos.eigenvalues()[i];
  }
  else
  {
    ev.resize(1);
    ev[0]=real(value_type(this->matrix()(0,0)));
  }
  if (this->calc_vectors()) {
    if  (this->dimension()>1) {
      // calculate eigen vectors
      ietl::Info<magnitude_type> info; // (m1, m2, ma, eigenvalue, residualm, status).
  
      try {
        eigenvectors.clear();
        lanczos.eigenvectors(lanczos.eigenvalues().begin(),lanczos.eigenvalues().begin()+n,
                             std::back_inserter(eigenvectors),info,generator); 
      }
      catch (std::runtime_error& e) {
        std::cout <<"Exception during eigenvector calculation: " <<  e.what() << "\n";
      }  
    }
    else {
      vector_type v(1);
      v[0]=1.;
      eigenvectors.push_back(v);
    }
  }
  this->perform_measurements();
  eigenvectors.clear();
  std::copy(ev.begin(),ev.end(),std::back_inserter(this->measurements_.rbegin()->average_values["Energy"]));
  this->eigenvalues_.push_back(ev);
}

template <class T>
std::vector<T> SparseDiagMatrix<T>::calculate(operator_matrix_type const& m) const
{
  using namespace boost::numeric::ublas;
  std::vector<value_type> av;
  for(typename std::vector<vector_type>::const_iterator it = eigenvectors.begin();it!=eigenvectors.end();it++)
    av.push_back(inner_prod(boost::numeric::ublas::conj(*it),prod(m,*it)));
  return av;
}

template <class T>
void SparseDiagMatrix<T>::print_eigenvectors(std::ostream& os) const
{
  unsigned int n=0;
  for(typename std::vector<vector_type>::const_iterator it = eigenvectors.begin();it!=eigenvectors.end();it++)
    os << "Eigenvector# " << n++ << ":\n" << *it << "\n";
}
