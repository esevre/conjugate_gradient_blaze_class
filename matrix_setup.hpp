//
// Created by Erik Sevre on 11/22/18.
//

#pragma once

#include <blaze/Math.h>


bool isConverged(double absolute_residual, double relative_residual)
{
    double relative_residual_tolerance = 1e-7;
    double absolute_residual_tolerance = 1e-7;

    if (std::abs(relative_residual) < relative_residual_tolerance) {
        return true;
    } else if (std::abs(absolute_residual) < absolute_residual_tolerance) {
        return true;
    } else {
        return false;
    }
}

template <class DynamicMatrixType>
void setup_matrix(DynamicMatrixType &A, size_t cols, size_t rows)
{
    const size_t COL_SZ = cols;
    const size_t ROW_SZ = rows;
    const size_t NDIM = COL_SZ * ROW_SZ;

    A.resize(NDIM, NDIM);

    for (size_t pos = 0; pos<NDIM; ++pos) {
        A(pos, pos) = -4;
        size_t sub_col_index = pos % COL_SZ;
        size_t sub_row_index = pos / COL_SZ;
        //
        // set up tri-diagonal indeces
        //
        if (sub_col_index > 0) {
            A(pos, pos-1) = 1;
        }
        if (sub_col_index < COL_SZ - 1) {
            A(pos, pos+1) = 1;
        }
        //
        // set up Identity blocks
        //
        if (sub_row_index > 0) {
            A(pos, pos-COL_SZ) = 1;
        }
        if (sub_row_index < COL_SZ - 1) {
            A(pos, pos+COL_SZ) = 1;
        }
    }
}



