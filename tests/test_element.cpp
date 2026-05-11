#include <gtest/gtest.h>
#include <cmath>
#include "Element.h"
#include "Node.h"

using namespace truss_fem;

constexpr double EPS = 1e-9;

TEST(ElementTest, LocalStiffnessHorizontalBar) {
    Node n0{0, 0.0, 0.0};
    Node n1{1, 1.0, 0.0};
    Element e{0, 0, 1, 210e9, 0.001};

    auto k_local = e.computeLocalStiffness(n0, n1);
    double expected = 210e9 * 0.001 / 1.0; // E*A/L, L=1.0
    EXPECT_NEAR(k_local(0, 0),  expected, EPS);
    EXPECT_NEAR(k_local(0, 1), -expected, EPS);
    EXPECT_NEAR(k_local(1, 0), -expected, EPS);
    EXPECT_NEAR(k_local(1, 1),  expected, EPS);
}

TEST(ElementTest, GlobalStiffnessHorizontalBar) {
    Node n0{0, 0.0, 0.0};
    Node n1{1, 1.0, 0.0};
    Element e{0, 0, 1, 210e9, 0.001};

    auto k_global = e.computeGlobalStiffness(n0, n1);
    double EA_L = 210e9 * 0.001 / 1.0;
    // Horizontal: c=1, s=0
    EXPECT_NEAR(k_global(0, 0),  EA_L, EPS);  // c²
    EXPECT_NEAR(k_global(0, 2), -EA_L, EPS);  // -c²
    EXPECT_NEAR(k_global(2, 0), -EA_L, EPS);  // -c²
    EXPECT_NEAR(k_global(2, 2),  EA_L, EPS);  // c²
    // s terms should be 0
    EXPECT_NEAR(k_global(0, 1), 0.0, EPS);
    EXPECT_NEAR(k_global(1, 1), 0.0, EPS);
}

TEST(ElementTest, GlobalStiffnessDiagonalBar) {
    Node n0{0, 0.0, 0.0};
    Node n1{1, 1.0, 1.0};
    Element e{0, 0, 1, 210e9, 0.001};

    auto k_global = e.computeGlobalStiffness(n0, n1);
    double L = std::sqrt(2.0);
    double EA_L = 210e9 * 0.001 / L;
    double c = 1.0 / std::sqrt(2.0);
    double s = 1.0 / std::sqrt(2.0);

    EXPECT_NEAR(k_global(0, 0), EA_L * c * c, EPS);
    EXPECT_NEAR(k_global(0, 1), EA_L * c * s, EPS);
    EXPECT_NEAR(k_global(1, 0), EA_L * c * s, EPS);
    EXPECT_NEAR(k_global(1, 1), EA_L * s * s, EPS);
    EXPECT_NEAR(k_global(2, 2), EA_L * c * c, EPS);
}

TEST(ElementTest, LengthComputation) {
    Node n0{0, 0.0, 0.0};
    Node n1{1, 3.0, 4.0};
    Element e{0, 0, 1, 210e9, 0.001};

    double L = e.computeLength(n0, n1);
    EXPECT_NEAR(L, 5.0, EPS);
}

TEST(ElementTest, StressComputation) {
    Node n0{0, 0.0, 0.0, 0.0, 0.0};
    Node n1{1, 1.0, 0.0, 0.001, 0.0};
    Element e{0, 0, 1, 210e9, 0.001};

    double strain = e.computeStrain(n0, n1);
    EXPECT_NEAR(strain, 0.001, EPS);

    double stress = e.computeStress(n0, n1);
    // stress = E * strain; compare against computed strain to avoid floating-point
    // amplification (~1e-16 strain error * 2.1e11 = ~2e-5 stress error)
    EXPECT_NEAR(stress, e.E * strain, EPS);
}

TEST(ElementTest, ThrowsOnZeroLength) {
    Node n0{0, 0.0, 0.0};
    Node n1{1, 0.0, 0.0}; // coincident -> L=0
    Element e{0, 0, 1, 210e9, 0.001};

    EXPECT_THROW(e.computeLocalStiffness(n0, n1), std::runtime_error);
    EXPECT_THROW(e.computeGlobalStiffness(n0, n1), std::runtime_error);
    EXPECT_THROW(e.computeStrain(n0, n1), std::runtime_error);
}

TEST(ElementTest, ZeroStrainForNoDisplacement) {
    Node n0{0, 0.0, 0.0, 0.0, 0.0};
    Node n1{1, 1.0, 0.0, 0.0, 0.0};
    Element e{0, 0, 1, 210e9, 0.001};

    double strain = e.computeStrain(n0, n1);
    EXPECT_NEAR(strain, 0.0, EPS);
}
