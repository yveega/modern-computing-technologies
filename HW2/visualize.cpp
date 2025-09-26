#include "inmost.h"
#include "solver.cpp"
#include <functional>
#include <matplot/matplot.h>
#include <vector>
#include <cmath>

using namespace INMOST;

double f(double x, double y)
{
    return 50.0 * sin(5 * x) * sin(5 * y);
}

double u(double x, double y)
{
    return sin(5 * x) * sin(5 * y);
}

int main() {
    size_t n = 40;
    std::vector<std::vector<double>> data(n + 1, std::vector<double>(n + 1, 0));
    Sparse::Vector res = solve_DE(f, 1.0, 1.0,
                [] (double x) { return u(x, 0.0); },
                [] (double x) { return u(x, 1.0); },
                [] (double y) { return u(0.0, y); },
                [] (double y) { return u(1.0, y); }, n);
    std::vector<std::vector<double> > Z(n + 1, std::vector<double>(n + 1, 0));
    for (size_t i = 1; i < n; i++) {
        for (size_t j = 1; j < n; j++) {
            size_t idx_res = (i - 1) * (n - 1) + j - 1;
            Z[i][j] = res[idx_res];
        }
    }
    double h = 1.0 / n;
    for (size_t i = 0; i < n + 1; i++) {
        Z[0][i] = u(0.0, i * h);
        Z[n][i] = u(1.0, i * h);
        Z[i][0] = u(i * h, 0.0);
        Z[i][n] = u(i * h, 1.0);
    }
    auto [X, Y] = matplot::meshgrid(matplot::linspace(0.0, 1.0, n + 1));
    auto plot = matplot::mesh(X, Y, Z)->palette_map_at_surface(true).face_alpha(0.5);
    matplot::rotate(23, 33);

    matplot::show();
}