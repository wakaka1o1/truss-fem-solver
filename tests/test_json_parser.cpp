#include <gtest/gtest.h>
#include "JsonParser.h"
#include <string>

using namespace truss_fem;

constexpr double EPS = 1e-9;

std::string fixturePath(const std::string& name) {
    return std::string(TEST_FIXTURES_DIR) + "/" + name;
}

TEST(JsonParserTest, ParseValidTwoBarTruss) {
    auto model = JsonParser::parse(fixturePath("two_bar_truss.json"));

    EXPECT_EQ(model.nodes.size(), static_cast<size_t>(3));
    EXPECT_EQ(model.elements.size(), static_cast<size_t>(2));
    EXPECT_EQ(model.bcs.size(), static_cast<size_t>(3));
    EXPECT_EQ(model.loads.size(), static_cast<size_t>(1));

    EXPECT_EQ(model.nodes[0].id, 0);
    EXPECT_NEAR(model.nodes[0].x, 0.0, EPS);
    EXPECT_NEAR(model.nodes[0].y, 0.0, EPS);

    EXPECT_EQ(model.elements[0].n1, 0);
    EXPECT_EQ(model.elements[0].n2, 2);
    EXPECT_NEAR(model.elements[0].E, 210e9, EPS);
    EXPECT_NEAR(model.elements[0].A, 0.001, EPS);

    EXPECT_EQ(model.bcs[1].nodeId, 0);
    EXPECT_EQ(model.bcs[1].dof, DoF::UY);
    EXPECT_NEAR(model.bcs[1].value, 0.0, EPS);

    EXPECT_EQ(model.loads[0].nodeId, 2);
    EXPECT_EQ(model.loads[0].dof, DoF::UY);
    EXPECT_NEAR(model.loads[0].value, -10000.0, EPS);
}

TEST(JsonParserTest, ThrowsOnMissingFile) {
    EXPECT_THROW(JsonParser::parse("nonexistent.json"), std::runtime_error);
}

TEST(JsonParserTest, ParseThreeBarTruss) {
    auto model = JsonParser::parse(fixturePath("three_bar_truss.json"));

    EXPECT_EQ(model.nodes.size(), static_cast<size_t>(4));
    EXPECT_EQ(model.elements.size(), static_cast<size_t>(3));
    EXPECT_EQ(model.bcs.size(), static_cast<size_t>(4));
    EXPECT_EQ(model.loads.size(), static_cast<size_t>(1));
}
