#include <vector>
#include <functional>
#include <iostream>
#include <cmath>
#include "inmost.h"

using namespace INMOST;

/* The function solves differential equation for u(x, y)
div(-D∇u) = f, (x, y) in (0; 1)×(0; 1)
u = g on the border of square (0; 1)×(0; 1)
D = diag{dx, dy} - diagonal diffusion tensor

f takes 2 doubles (x and y)
dx and dy are two doubles - coefficients of diffusion tensor
g is passed as 4 functions:
g_bottom(x) = g(x, 0)
g_top(x) = g(x, 1)
g_left(y) = g(0, y)
g_right(y) = g(1, y)

n is the size of grid

the result is vector of values in the inner points in grid, ordered by rows
*/
Sparse::Vector
solve_DE(std::function<double(double, double)> f, double dx, double dy,
            std::function<double(double)> g_bottom,
            std::function<double(double)> g_top,
            std::function<double(double)> g_left,
            std::function<double(double)> g_right, size_t n)
{
    double h = 1.0 / n;
    // Create matrix and vectors for system Ax = b
    Sparse::Matrix A;
    Sparse::Vector b;
    Sparse::Vector x;
    // Set their size
    A.SetInterval(0, (n - 1) * (n - 1));
    b.SetInterval(0, (n - 1) * (n - 1));
    x.SetInterval(0, (n - 1) * (n - 1));

    // function which modifies b if (i_set, j_set) is a border element
    // and puts -1 in matrix A for non-border elements
    // row is index of current row of system Ax = b
    std::function<void(size_t, size_t, size_t, double)> set_coeff =
        [&] (size_t row, size_t i_set, size_t j_set, double d) {
            if (i_set == 0) {
                b[row] += d * g_bottom(j_set * h) / h / h;
            } else if (i_set == n) {
                b[row] += d * g_top(j_set * h) / h / h;
            } else if (j_set == 0) {
                b[row] += d * g_left(i_set * h) / h / h;
            } else if (j_set == n) {
                b[row] += d * g_right(i_set * h) / h / h;
            } else {
                A[row][(i_set - 1) * (n - 1) + j_set - 1] = -d;
            }
        };

    // Building matrix A and vector b
    for (size_t i = 1; i < n; i++) {
        for (size_t j = 1; j < n; j++) {
            size_t idx = (i - 1) * (n - 1) + j - 1;
            A[idx][idx] = 2.0 * (dx + dy); // diagonal element
            b[idx] = f(i * h, j * h);
            set_coeff(idx, i - 1, j, dy);
            set_coeff(idx, i + 1, j, dy);
            set_coeff(idx, i, j - 1, dx);
            set_coeff(idx, i, j + 1, dx);
        }
    }

    // Solving system Ax = b
    Solver S(Solver::INNER_ILU2);
    S.SetParameter("absolute_tolerance", "1e-12");
    S.SetParameter("relative_tolerance", "1e-12");
    S.SetParameter("drop_tolerance", "0.005");
    S.SetMatrix(A);
    S.Solve(b, x);
    std::cout << S.Iterations() << ' ' << S.IterationsTime() << ' ';

    // Multiplying resulted vector x by h*h
    for (size_t i = 0; i < (n - 1) * (n - 1); i++) {
        x[i] *= h * h;
    }
    return x;
}