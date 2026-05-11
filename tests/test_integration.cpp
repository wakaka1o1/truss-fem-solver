#include <gtest/gtest.h>
#include "JsonParser.h"
#include "Assembler.h"
#include "BoundaryConditionHandler.h"
#include "GaussianSolver.h"
#include "PostProcessor.h"
#include "Matrix.h"
#include "TrussModel.h"
#include <cmath>
#include <string>

using namespace truss_fem;

constexpr double EPS = 1e-6;

std::string fixturePath(const std::string& name) {
    return std::string(TEST_FIXTURES_DIR) + "/" + name;
}

struct SolverResult {
    Vector displacements;
    Vector reactions;
};

SolverResult runSolver(const std::string& inputFile) {
    TrussModel model = JsonParser::parse(inputFile);
    Matrix K_original = Assembler::assembleStiffness(model);
    Vector F = Assembler::assembleForceVector(model);
    Matrix K = K_original;
    BoundaryConditionHandler::apply(K, F, model);
    GaussianSolver solver;
    Vector u = solver.solve(K, F);
    Vector reactions = PostProcessor::computeReactions(model, u, K_original);
    return {u, reactions};
}

TEST(IntegrationTest, TwoBarSymmetricTruss) {
    auto [u, reactions] = runSolver(fixturePath("two_bar_truss.json"));

    // 3 nodes -> 6 DOFs
    EXPECT_EQ(u.size(), 6u);

    // Node 0 fully fixed -> zero displacement
    EXPECT_NEAR(u[0], 0.0, EPS); // ux
    EXPECT_NEAR(u[1], 0.0, EPS); // uy

    // Node 1 fully fixed -> zero displacement
    EXPECT_NEAR(u[2], 0.0, EPS); // ux
    EXPECT_NEAR(u[3], 0.0, EPS); // uy

    // Node 2 loaded downward -> should deflect downward
    EXPECT_LT(u[5], 0.0); // uy < 0 (downward)

    // Symmetric truss: each support takes half the vertical load
    double expected_uy_reaction = 5000.0; // 10000 / 2
    EXPECT_NEAR(reactions[1], expected_uy_reaction, EPS * expected_uy_reaction);
    EXPECT_NEAR(reactions[3], expected_uy_reaction, EPS * expected_uy_reaction);

    // Horizontal reactions should be equal and opposite (symmetric)
    EXPECT_NEAR(reactions[0] + reactions[2], 0.0, EPS);
}

TEST(IntegrationTest, ThreeBarIndeterminateTruss) {
    auto [u, reactions] = runSolver(fixturePath("three_bar_truss.json"));

    // 4 nodes -> 8 DOFs
    EXPECT_EQ(u.size(), 8u);

    // Node 3 loaded downward
    EXPECT_LT(u[7], 0.0); // Node 3 uy < 0

    // Vertical reactions should sum to applied load (15,000 N)
    double sum_uy = 0.0;
    for (const auto& bc : JsonParser::parse(fixturePath("three_bar_truss.json")).bcs) {
        if (bc.dof == DoF::UY) {
            int idx = 2 * bc.nodeId + 1;
            sum_uy += reactions[idx];
        }
    }
    EXPECT_NEAR(sum_uy, 15000.0, EPS * 15000.0);

    // Horizontal reactions should sum to zero
    double sum_ux = 0.0;
    for (const auto& bc : JsonParser::parse(fixturePath("three_bar_truss.json")).bcs) {
        if (bc.dof == DoF::UX) {
            int idx = 2 * bc.nodeId;
            sum_ux += reactions[idx];
        }
    }
    EXPECT_NEAR(sum_ux, 0.0, EPS);
}

TEST(IntegrationTest, StressNonZeroForLoadedTruss) {
    auto [u, reactions] = runSolver(fixturePath("two_bar_truss.json"));

    TrussModel model = JsonParser::parse(fixturePath("two_bar_truss.json"));
    std::vector<Node> nodesWithDisp = model.nodes;
    for (size_t i = 0; i < nodesWithDisp.size(); ++i) {
        nodesWithDisp[i].ux = u[2 * i];
        nodesWithDisp[i].uy = u[2 * i + 1];
    }

    for (const auto& e : model.elements) {
        const auto& n1 = nodesWithDisp[static_cast<size_t>(e.n1)];
        const auto& n2 = nodesWithDisp[static_cast<size_t>(e.n2)];
        double stress = e.computeStress(n1, n2);
        // Stress should be non-zero (the truss is loaded)
        EXPECT_NE(std::abs(stress), 0.0);
    }
}
