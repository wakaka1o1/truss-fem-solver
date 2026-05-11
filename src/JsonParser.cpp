#include "JsonParser.h"
#include <fstream>
#include <stdexcept>
#include <nlohmann/json.hpp>

namespace truss_fem {

namespace {

using json = nlohmann::json;

DoF dofFromString(const std::string& s) {
    if (s == "ux") return DoF::UX;
    if (s == "uy") return DoF::UY;
    throw std::runtime_error("Invalid DOF: " + s + " (expected 'ux' or 'uy')");
}

void validateElement(const Element& e, int numNodes) {
    if (e.n1 < 0 || e.n1 >= numNodes || e.n2 < 0 || e.n2 >= numNodes) {
        throw std::runtime_error("Element " + std::to_string(e.id)
            + " references out-of-range node");
    }
    if (e.n1 == e.n2) {
        throw std::runtime_error("Element " + std::to_string(e.id)
            + " has identical nodes (zero length)");
    }
    if (e.E <= 0.0 || e.A <= 0.0) {
        throw std::runtime_error("Element " + std::to_string(e.id)
            + " has non-positive E or A");
    }
}

void validateBC(const BoundaryCondition& bc, int numNodes,
                const std::vector<bool>& constrainedDofs) {
    if (bc.nodeId < 0 || bc.nodeId >= numNodes) {
        throw std::runtime_error("BC references out-of-range node " + std::to_string(bc.nodeId));
    }
    int idx = 2 * bc.nodeId + (bc.dof == DoF::UY ? 1 : 0);
    if (constrainedDofs[idx]) {
        throw std::runtime_error("Duplicate BC on node " + std::to_string(bc.nodeId)
            + " DOF " + (bc.dof == DoF::UX ? "ux" : "uy"));
    }
}

} // anonymous namespace

TrussModel JsonParser::parse(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filepath);
    }

    json j;
    try {
        file >> j;
    } catch (const json::parse_error& e) {
        throw std::runtime_error("JSON parse error in " + filepath + ": " + e.what());
    }

    TrussModel model;

    // Parse nodes
    if (!j.contains("nodes") || !j["nodes"].is_array()) {
        throw std::runtime_error("Missing or invalid 'nodes' array");
    }
    for (const auto& nj : j["nodes"]) {
        model.nodes.push_back({
            nj.at("id").get<int>(),
            nj.at("x").get<double>(),
            nj.at("y").get<double>()
        });
    }

    int numNodes = static_cast<int>(model.nodes.size());

    // Parse elements
    if (!j.contains("elements") || !j["elements"].is_array()) {
        throw std::runtime_error("Missing or invalid 'elements' array");
    }
    for (const auto& ej : j["elements"]) {
        Element e{
            ej.at("id").get<int>(),
            ej.at("n1").get<int>(),
            ej.at("n2").get<int>(),
            ej.at("E").get<double>(),
            ej.at("A").get<double>()
        };
        validateElement(e, numNodes);
        model.elements.push_back(e);
    }

    // Parse boundary conditions
    std::vector<bool> constrainedDofs(2 * numNodes, false);
    if (j.contains("boundary_conditions")) {
        for (const auto& bj : j["boundary_conditions"]) {
            BoundaryCondition bc{
                bj.at("node").get<int>(),
                dofFromString(bj.at("dof").get<std::string>()),
                bj.at("value").get<double>()
            };
            validateBC(bc, numNodes, constrainedDofs);
            int idx = 2 * bc.nodeId + (bc.dof == DoF::UY ? 1 : 0);
            constrainedDofs[idx] = true;
            model.bcs.push_back(bc);
        }
    }

    // Parse loads
    if (j.contains("loads")) {
        for (const auto& lj : j["loads"]) {
            model.loads.push_back({
                lj.at("node").get<int>(),
                dofFromString(lj.at("dof").get<std::string>()),
                lj.at("value").get<double>()
            });
        }
    }

    return model;
}

} // namespace truss_fem
