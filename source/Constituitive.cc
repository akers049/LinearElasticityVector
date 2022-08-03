#ifndef CONSTITUITIVE_CC_
#define CONSTITUITIVE_CC_

#include "Constituitive.h"


using namespace dealii;

  inline
  double LinearElastic::get_energy(Tensor<2, DIM> &grad_u)
  {
    double W = 0.0;
    Tensor<2, DIM> eps = symmetrize(grad_u);

    W = 0.5*double_contract<0,0,1,1>(double_contract<2, 0, 3, 1>(C, eps), eps);

    return W;
  }

  inline
  void LinearElastic::get_sigma(Tensor<2, DIM> &grad_u, Tensor<2, DIM> &sigma)
  {

    Tensor<2, DIM> eps = symmetrize(grad_u);
    sigma = double_contract<2, 0, 3, 1>(C, eps);
  }


#endif
