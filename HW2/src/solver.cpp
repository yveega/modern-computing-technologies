#include <vector>
#include <functional>
#include <iostream>
#include <cmath>
#include "inmost.h"

using namespace INMOST;

Sparse::Vector
solve_DE(std::function<double(double, double)> f,
            std::function<double(double)> g_bottom,
            std::function<double(double)> g_top,
            std::function<double(double)> g_left,
            std::function<double(double)> g_right, size_t n)
{
    double h = 1.0 / n;
    Sparse::Matrix A;
    Sparse::Vector b;
    Sparse::Vector x;
    // Set their size
    A.SetInterval(0, (n - 1) * (n - 1));
    b.SetInterval(0, (n - 1) * (n - 1));
    x.SetInterval(0, (n - 1) * (n - 1));

    std::function<void(size_t, size_t, size_t)> set_coeff =
        [&] (size_t row, size_t i_set, size_t j_set) {
            if (i_set == 0) {
                b[row] += g_bottom(j_set * h) / h / h;
            } else if (i_set == n) {
                b[row] += g_top(j_set * h) / h / h;
            } else if (j_set == 0) {
                b[row] += g_left(i_set * h) / h / h;
            } else if (j_set == n) {
                b[row] += g_right(i_set * h) / h / h;
            } else {
                A[row][(i_set - 1) * (n - 1) + j_set - 1] = -1.0;
            }
        };
    for (size_t i = 1; i < n; i++) {
        for (size_t j = 1; j < n; j++) {
            size_t idx = (i - 1) * (n - 1) + j - 1;
            A[idx][idx] = 4.0;
            b[idx] = f(i * h, j * h);
            set_coeff(idx, i - 1, j);
            set_coeff(idx, i + 1, j);
            set_coeff(idx, i, j - 1);
            set_coeff(idx, i, j + 1);
            
            // for (size_t k = 0; k < (n - 1) * (n - 1); k++) {
            //     std::cout << A[idx][k] << ' ';
            // }
            // std::cout << "|  " << b[idx] << std::endl;
        }
    }
    Solver S(Solver::INNER_ILU2);
    S.SetParameter("absolute_tolerance", "1e-10");
    S.SetParameter("relative_tolerance", "1e-6");
    S.SetMatrix(A);
    bool solve = S.Solve(b, x);
    for (size_t i = 0; i < (n - 1) * (n - 1); i++) {
        x[i] *= h * h;
    }
    return x;
}

double f(double x, double y)
{
    return 50.0 * sin(5*x) * sin(5*y);
}

double u(double x, double y)
{
    return sin(5*x) * sin(5*y);
}

double g(double x) {
    return 1.0;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    size_t n = 100;

    Sparse::Vector x = solve_DE(f,
                [] (double x) { return u(x, 0.0); },
                [] (double x) { return u(x, 1.0); },
                [] (double y) { return u(0.0, y); },
                [] (double y) { return u(1.0, y); }, n);

    double h = 1.0 / n;
    double norm = 0.0;
    for (size_t i = 1; i < n; i++) {
        for (size_t j = 1; j < n; j++) {
            size_t idx = (i - 1) * (n - 1) + j - 1;
            double diff = u(i * h, j * h) - x[idx];
            norm += diff * diff;
        }
    }
    std::cout << "ERROR NORM: " << norm << std::endl;
    return 0;
}