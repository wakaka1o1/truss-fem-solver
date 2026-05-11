#include <gtest/gtest.h>
#include "Assembler.h"
#include "TrussModel.h"
#include "Node.h"
#include "Element.h"
#include "BoundaryCondition.h"

using namespace truss_fem;

constexpr double EPS = 1e-9;

TrussModel makeTwoBarModel() {
    TrussModel m;
    m.nodes = {{0, 0.0, 0.0}, {1, 1.0, 0.0}, {2, 0.5, 1.0}};
    m.elements = {
        {0, 0, 2, 210e9, 0.001},
        {1, 1, 2, 210e9, 0.001}
    };
    return m;
}

TEST(AssemblerTest, GlobalStiffnessDimensions) {
    auto model = makeTwoBarModel();
    Matrix K = Assembler::assembleStiffness(model);

    int n = model.totalDofs(); // 3 nodes -> 6 DOFs
    EXPECT_EQ(K.size(), n);
}

TEST(AssemblerTest, StiffnessMatrixSymmetry) {
    auto model = makeTwoBarModel();
    Matrix K = Assembler::assembleStiffness(model);

    for (int i = 0; i < K.size(); ++i) {
        for (int j = 0; j < K.size(); ++j) {
            EXPECT_NEAR(K(i, j), K(j, i), EPS);
        }
    }
}

TEST(AssemblerTest, ForceVectorSize) {
    auto model = makeTwoBarModel();
    model.loads = {{2, DoF::UY, -10000.0}};

    Vector F = Assembler::assembleForceVector(model);
    EXPECT_EQ(static_cast<int>(F.size()), model.totalDofs());
}

TEST(AssemblerTest, ForceVectorAssembly) {
    auto model = makeTwoBarModel();
    model.loads = {{2, DoF::UY, -10000.0}};

    Vector F = Assembler::assembleForceVector(model);
    // Node 2 UY = DOF index 5
    EXPECT_NEAR(F[5], -10000.0, EPS);
    // All other entries should be 0
    EXPECT_NEAR(F[0], 0.0, EPS);
    EXPECT_NEAR(F[1], 0.0, EPS);
    EXPECT_NEAR(F[2], 0.0, EPS);
    EXPECT_NEAR(F[3], 0.0, EPS);
    EXPECT_NEAR(F[4], 0.0, EPS);
}
