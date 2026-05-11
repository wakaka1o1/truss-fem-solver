#include <gtest/gtest.h>
#include "BoundaryConditionHandler.h"
#include "TrussModel.h"

using namespace truss_fem;

constexpr double EPS = 1e-9;

TEST(BCHandlerTest, ApplySingleBC) {
    TrussModel model;
    model.nodes = {{0, 0.0, 0.0}, {1, 1.0, 0.0}};
    model.bcs = {{0, DoF::UX, 0.0}};

    Matrix K(4);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            K(i, j) = 1.0;

    Vector F = {5.0, 3.0, 7.0, 2.0};

    BoundaryConditionHandler::apply(K, F, model);

    int d = 0; // DOF 0 (Node 0, UX)
    for (int j = 0; j < 4; ++j) {
        if (j == d) {
            EXPECT_NEAR(K(d, j), 1.0, EPS);
        } else {
            EXPECT_NEAR(K(d, j), 0.0, EPS);
        }
    }
    for (int i = 0; i < 4; ++i) {
        if (i != d) {
            EXPECT_NEAR(K(i, d), 0.0, EPS);
        }
    }
    EXPECT_NEAR(F[d], 0.0, EPS);
}

TEST(BCHandlerTest, ApplyBCOffsetsForceVector) {
    TrussModel model;
    model.nodes = {{0, 0.0, 0.0}, {1, 1.0, 0.0}};
    model.bcs = {{1, DoF::UY, 0.005}};

    Matrix K(4);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            K(i, j) = 100.0;

    Vector F = {0.0, 0.0, 0.0, 0.0};

    BoundaryConditionHandler::apply(K, F, model);

    int d = 3;
    double v = 0.005;
    for (int j = 0; j < 4; ++j) {
        if (j != d) {
            EXPECT_NEAR(F[j], -100.0 * v, EPS);
        }
    }
    EXPECT_NEAR(F[d], v, EPS);
}

TEST(BCHandlerTest, NoBCsDoesNothing) {
    TrussModel model;
    model.nodes = {{0, 0.0, 0.0}, {1, 1.0, 0.0}};

    Matrix K(4);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            K(i, j) = 42.0;

    Vector F = {1.0, 2.0, 3.0, 4.0};
    auto F_copy = F;

    BoundaryConditionHandler::apply(K, F, model);

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            EXPECT_NEAR(K(i, j), 42.0, EPS);
        }
        EXPECT_NEAR(F[i], F_copy[i], EPS);
    }
}
