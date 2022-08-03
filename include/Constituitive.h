#ifndef CONSTITUITIVE_H_
#define CONSTITUITIVE_H_

#include <iostream>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/symmetric_tensor.h>

#define DIM 2

using namespace dealii;


  class LinearElastic
  {
    public:

    virtual ~LinearElastic(){};
    LinearElastic(double E, double nu)
    {
      mu = E/(2.0*(1 + nu));
      lambda = E*nu/((1 + nu)*(1 - 2.0*nu));

      for (unsigned int i = 0; i < DIM; i ++)
        for (unsigned int j = 0; j < DIM; j ++)
          for (unsigned int k = 0; k < DIM; k ++)
            for (unsigned int l = 0; l < DIM; l ++)
              C[i][j][k][l] = ((i == j) && (k ==l) ? lambda : 0.0) +
                ((i == k) && (j ==l) ? mu : 0.0) + ((i == l) && (j ==k) ? mu : 0.0);
    };


    double get_energy(Tensor<2, DIM> &grad_u);
    void get_sigma(Tensor<2, DIM> &grad_u, Tensor<2, DIM> &sigma);

    Tensor<4, DIM> C;

    private:
      double lambda;
      double mu;
  };

  #endif
