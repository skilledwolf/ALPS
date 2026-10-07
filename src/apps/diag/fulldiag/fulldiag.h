/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1994-2009 by Matthias Troyer <troyer@comp-phys.org>,
*                            Andreas Honecker <ahoneck@uni-goettingen.de>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <cassert>

#include <alps/config.h> // needed to set up correct bindings
#include <boost/numeric/bindings/lapack/driver/syev.hpp>
#include <boost/numeric/bindings/lapack/driver/heev.hpp>

#include "../diag.h"

#include <alps/numeric/matrix.hpp>
#include <alps/numeric/matrix/ublas_sparse_functions.hpp>

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/bindings/ublas.hpp>
#include <boost/numeric/bindings/upper.hpp>
#include <boost/numeric/bindings/lower.hpp>

template <class T> 
struct pick_real_complex 
{
  template <class Matrix, class Vector>
  static int heev (const char jobz, Matrix const& m, Vector& v)
   { return boost::numeric::bindings::lapack::syev(jobz,m,v);}
};

template <class T>
struct pick_real_complex<std::complex<T> > 
{
  template <class Matrix, class Vector>
  static int heev (const char jobz, Matrix const& m, Vector& v)
   { return boost::numeric::bindings::lapack::heev(jobz,m,v);}
};

template <class T>
class FullDiagMatrix : public DiagMatrix<T,alps::numeric::matrix<T> >
{
public:
  typedef DiagMatrix<T,alps::numeric::matrix<T> > super_type;
  typedef T value_type;
  typedef typename super_type::magnitude_type magnitude_type;
  typedef typename super_type::matrix_type matrix_type;
  typedef typename super_type::vector_type vector_type;
  typedef typename super_type::size_type size_type;
  typedef typename super_type::mag_vector_type mag_vector_type;
  typedef typename super_type::half_integer_type half_integer_type;
  typedef typename super_type::site_iterator site_iterator;
  typedef typename super_type::site_descriptor site_descriptor;
  typedef typename super_type::operator_matrix_type operator_matrix_type;
  
  explicit FullDiagMatrix (alps::Parameters const& p) : super_type(p) {}
  void print_eigenvectors(std::ostream& os) const;
private:
  void do_subspace();
  std::vector<value_type> calculate(operator_matrix_type const& m) const;
};

template <class T>
std::vector<T> FullDiagMatrix<T>::calculate(operator_matrix_type const& m) const
{
  std::vector<value_type> av;
  for (unsigned i=0; i < num_cols(this->matrix()); ++i) {
    alps::numeric::column_view<matrix_type const> const v(this->matrix(),i);
    //av.push_back(scalar_product(conj(v),m*v)); // scalar_product uses conj_mult [M. Pikulski]
    av.push_back(scalar_product(v,m*v));
  }
  return av;
}


template <class T>
void FullDiagMatrix<T>::print_eigenvectors(std::ostream& os) const
{
  for (unsigned i=0;i<num_cols(this->matrix());++i)
    os << "Eigenvector# " << i << ":\n" << alps::numeric::column_view<matrix_type const>(this->matrix(),i) << "\n";
}

template <class T>
void FullDiagMatrix<T>::do_subspace()
{
  using std::copy;
  this->build();
  if (this->dimension()) {
    mag_vector_type eigenvalues(this->dimension());
    pick_real_complex<T>::heev(this->calc_vectors() ? 'V' : 'N',boost::numeric::bindings::upper(this->matrix()),eigenvalues);
    this->perform_measurements();
    copy(eigenvalues.begin(),eigenvalues.end(),std::back_inserter(this->measurements_.rbegin()->average_values["Energy"]));
    this->eigenvalues_.push_back(eigenvalues);
  }
}
