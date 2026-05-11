#include "GaussianSolver.h"
#include <cmath>
#include <stdexcept>
#include <string>

namespace truss_fem {

Vector GaussianSolver::solve(const Matrix& K, const Vector& F) {
    int n = K.size();

    // Copy K into working matrix A
    Matrix A(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A(i, j) = K(i, j);
        }
    }

    Vector b = F;

    // Forward elimination with partial pivoting
    for (int k = 0; k < n; ++k) {
        if (std::abs(A(k, k)) < 1e-15) {
            bool found = false;
            for (int i = k + 1; i < n; ++i) {
                if (std::abs(A(i, k)) > 1e-15) {
                    for (int j = k; j < n; ++j) {
                        std::swap(A(k, j), A(i, j));
                    }
                    std::swap(b[k], b[i]);
                    found = true;
                    break;
                }
            }
            if (!found) {
                std::cerr << "Warning: near-zero pivot at DOF " << k
                          << " (possible mechanism)" << std::endl;
                throw std::runtime_error("Singular matrix: zero pivot at row " + std::to_string(k));
            }
        }

        for (int i = k + 1; i < n; ++i) {
            double factor = A(i, k) / A(k, k);
            for (int j = k; j < n; ++j) {
                A(i, j) -= factor * A(k, j);
            }
            b[i] -= factor * b[k];
        }
    }

    // Back substitution
    Vector u(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        double sum = b[i];
        for (int j = i + 1; j < n; ++j) {
            sum -= A(i, j) * u[j];
        }
        u[i] = sum / A(i, i);
    }

    return u;
}

} // namespace truss_fem
