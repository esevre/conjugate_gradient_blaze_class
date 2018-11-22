//
// Created by Erik Sevre on 11/22/18.
//
// This code is takes notes from blaze-iterative
//   https://github.com/tjolsen/BlazeIterative

#pragma once

#include <iostream>

#include <blaze/Math.h>
#include "matrix_setup.hpp"

/// conjugate_gradient will solve the system Ax=b for
/// x where A is a symetric positive definite matrix A
///
/// \tparam DataType
/// \param A input matrix
/// \param b input vector
/// \param x output vector
template <class DataType>
void conjugate_gradient(
        const blaze::DynamicMatrix<DataType> &A,
        const blaze::DynamicVector<DataType> &b,
              blaze::DynamicVector<DataType> &x
        )
{
    using namespace blaze;

    DynamicVector<DataType> r = b - A*x;
    DynamicVector<DataType> p = r;
    DynamicVector<DataType> Ap = A*p;

    DataType residual_norm_0 = trans(r)*r;
    DataType residual_norm_prev = residual_norm_0;
    DataType residual_norm = residual_norm_0;

    size_t iteration{0};
    const size_t iteration_max{100};

    while (true) {

        auto alpha = residual_norm / (trans(p)*Ap);
        x += alpha * p;
        r -= alpha * Ap;

        if (iteration >= iteration_max) break;
        if (isConverged(residual_norm, residual_norm/residual_norm_0)) break;

        residual_norm_prev = residual_norm;
        residual_norm = trans(r)*r;

        auto beta = residual_norm / residual_norm_prev;

        p = r + beta * p;
        Ap = A*p;

        ++iteration;
    }


    std::cout << "CG Converged in " << iteration << " iterations\n";

}


