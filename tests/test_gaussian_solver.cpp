#include <gtest/gtest.h>
#include "GaussianSolver.h"
#include <stdexcept>

using namespace truss_fem;

constexpr double EPS = 1e-9;

TEST(GaussianSolverTest, Solve2x2System) {
    Matrix K(2);
    K(0, 0) = 2.0; K(0, 1) = 1.0;
    K(1, 0) = 1.0; K(1, 1) = 3.0;

    Vector F = {5.0, 6.0};

    GaussianSolver solver;
    Vector u = solver.solve(K, F);

    // 2x + y = 5, x + 3y = 6 → x=1.8, y=1.4
    EXPECT_NEAR(u[0], 1.8, EPS);
    EXPECT_NEAR(u[1], 1.4, EPS);
}

TEST(GaussianSolverTest, Solve4x4Identity) {
    Matrix K(4);
    for (int i = 0; i < 4; ++i) K(i, i) = 1.0;

    Vector F = {1.0, 2.0, 3.0, 4.0};

    GaussianSolver solver;
    Vector u = solver.solve(K, F);

    for (int i = 0; i < 4; ++i) {
        EXPECT_NEAR(u[i], F[i], EPS);
    }
}

TEST(GaussianSolverTest, SolveDiagonalSystem) {
    Matrix K(3);
    K(0, 0) = 3.0; K(1, 1) = 5.0; K(2, 2) = 7.0;
    Vector F = {9.0, 10.0, 14.0};

    GaussianSolver solver;
    Vector u = solver.solve(K, F);

    EXPECT_NEAR(u[0], 3.0, EPS);
    EXPECT_NEAR(u[1], 2.0, EPS);
    EXPECT_NEAR(u[2], 2.0, EPS);
}

TEST(GaussianSolverTest, ThrowsOnSingularMatrix) {
    // Linearly dependent rows: row1 = 2 * row0
    Matrix K(2);
    K(0, 0) = 1.0; K(0, 1) = 2.0;
    K(1, 0) = 2.0; K(1, 1) = 4.0;
    Vector F = {3.0, 6.0};

    GaussianSolver solver;
    EXPECT_THROW(solver.solve(K, F), std::runtime_error);
}
