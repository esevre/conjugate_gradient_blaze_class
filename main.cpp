#include <iostream>
#include <cmath>

#include <blaze/Math.h>

#include "conjugate_gradient.hpp"
#include "matrix_setup.hpp"


int main()
{
    const size_t rows = 50;
    const size_t cols = 50;
    const size_t dims = rows * cols;

    blaze::DynamicMatrix<double> A(dims, dims);
    blaze::DynamicVector<double> x(dims);
    blaze::DynamicVector<double> b(dims);

    setup_matrix(A, cols, rows);
    for (auto &elt : b) {
        elt = 1;
    }

    conjugate_gradient(A, b, x);

    blaze::DynamicVector<double> err = b - A*x;
    double err_sum = 0;
    for (auto &e : err) {
        err_sum += e*e;
    }
    std::cout << "error sum:\n" << err_sum << "\n";


}