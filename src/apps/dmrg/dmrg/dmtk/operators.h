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

#ifndef __DMTK_OPERATORS_H__
#define __DMTK_OPERATORS_H__

#include "product.h"
#include <memory>

namespace dmtk {

template<class T>
void
product(const BasicOp<T> &the_op,
        const VectorState<T> &v, VectorState<T> &res,
        size_t m, T coef = T(1), bool hc = false, DMTKglobals<T> *globals = NULL)
{
  if (the_op.name() == "I") {
    VectorState<T> aux(v);
    aux *= coef;
    res += aux;
  }
  std::unique_ptr<DMTKglobals<T> > local;
  if (!globals) {
    local.reset(new DMTKglobals<T>);
    globals = local.get();
  }
  for_each_product_term(the_op, v, res, m, coef, hc,
    [&](const ProductTerm<T>& term) { product_term1(term, v, res, MASK_PRODUCT_DEFAULT, globals); });
}

template<class T>
void
product(const BasicOp<T> &op1, const BasicOp<T>& op2,
        const VectorState<T> &v, VectorState<T> &res,
        size_t m1, size_t m2, T coef = T(1), bool hc = false, DMTKglobals<T> *globals = NULL, bool use_condensed = false)
{
  std::unique_ptr<DMTKglobals<T> > local;
  if (!globals) {
    local.reset(new DMTKglobals<T>);
    globals = local.get();
  }
  for_each_product_term(op1, op2, v, res, m1, m2, coef, hc,
    [&](const ProductTerm<T>& term) { product_term2(term, v, res, MASK_PRODUCT_DEFAULT | MASK_PRODUCT_HC, globals, use_condensed); });
}

template<class T>
void
product(const BasicOp<T> &op1, const BasicOp<T>& op2,
        const BasicOp<T> &op3, const BasicOp<T>& op4,
        const VectorState<T> &v, VectorState<T> &res,
        size_t m1, size_t m2, size_t m3, size_t m4,
        T coef = T(1), bool hc = false, DMTKglobals<T> *globals = NULL)
{
  for_each_product_term(op1, op2, op3, op4, v, res, m1, m2, m3, m4, coef, hc,
    [&](const ProductTerm<T>& term) { product_term4(term, v, res, MASK_PRODUCT_DEFAULT | MASK_PRODUCT_HC, globals); });
}

template<class T>
void
product(const BasicOp<T> &op1, const BasicOp<T>& op2,
        const BasicOp<T> &op3, 
        const VectorState<T> &v, VectorState<T> &res,
        size_t m1, size_t m2, size_t m3,
        T coef = T(1), bool hc = false, DMTKglobals<T> *globals = NULL)
{
  for_each_product_term(op1, op2, op3, v, res, m1, m2, m3, coef, hc,
    [&](const ProductTerm<T>& term) { product_term3(term, v, res, MASK_PRODUCT_DEFAULT | MASK_PRODUCT_HC, globals); });
}

template<class T>
BasicOp<T>
product(const BasicOp<T> &op2, const BasicOp<T> &op1)
{
  BasicOp<T> res;
  QN dqn1 = op1.dqn;
  QN dqn2 = op2.dqn;
  res.dqn = dqn1+dqn2;
  res.repack(op2.subspaces());

  typename BMatrix<T>::const_iterator iter;
  for(iter = op1.begin(); iter != op1.end(); iter++){
    const SubMatrix<T> &sm1 = (*iter);
    QN qn1 = sm1.qn();
    const SubMatrix<T> *_sm2 = op2.block(qn1+dqn1);
    if(_sm2){
      const SubMatrix<T> &sm2 = *_sm2;
      if(sm1.rows() == sm2.cols()){
        SubMatrix<T> sm(qn1,sm1.col_range(),sm2.row_range());
        sm=(product(sm2,sm1));
        res.push_back(sm); 
      } else {
        cout << "ERROR : op1 and op2 are inconsistent\n";
      }
    }  
  }
  return res;
}
//////////////////////////////////////////////////////////////////////////
template<class T>
BasicOp<T> 
BasicOp<T>::reshape(int qn_mask) const
{
  BasicOp<T> new_op;
  new_op = this->internals(); 
  //for(int i = 0; i < QN::QN_LAST; i++)
    //if(IBITS(qn_mask,i) == 0) new_op.dqn[i] = 0;

  PackedBasis::const_iterator biter;
  PackedBasis new_subspaces;

  new_subspaces = this->subspaces().reshape(qn_mask);
  new_op.resize(new_subspaces);

/*
  cout << "=========================\n";
  for(biter = new_subspaces.begin(); biter != new_subspaces.end(); biter++)
   cout << "HOLA BEGIN " << (*biter).start() << " END " << (*biter).end() << " QN " << (*biter).qn().n() << " " << (*biter).dim() << endl;
  cout << "=========================\n";
*/

  typename BMatrix<T>::iterator iter;
  typename BMatrix<T>::const_iterator citer;
  for(iter = new_op.begin(); iter != new_op.end(); iter++){
    SubMatrix<T> &sub = (*iter);
    QN new_qn = sub.qn();

    for(citer = this->begin(); citer != this->end(); citer++){
      const SubMatrix<T> &m = (*citer);
      if(new_qn.equal(m.qn(), qn_mask)){
        for(int col = 0; col < m.cols(); col++){
          int new_col = col + m.col_range().begin()-sub.col_range().begin();
          for(int row = 0; row < m.rows(); row++){
            int new_row = row + m.row_range().begin()-sub.row_range().begin();
            sub(new_col,new_row) = m(col,row);
          }
        }
      } 
    }
  } 

  return new_op;
}
//////////////////////////////////////////////////////////////////////////

} // dmtk

#endif // __DMTK_OPERATORS_H__
