#include <iostream>
#include "JsonParser.h"
#include "Assembler.h"
#include "BoundaryConditionHandler.h"
#include "GaussianSolver.h"
#include "PostProcessor.h"

using namespace truss_fem;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: truss_solver <input.json>\n";
        return 1;
    }

    try {
        TrussModel model = JsonParser::parse(argv[1]);
        Matrix K_original = Assembler::assembleStiffness(model);
        Vector F = Assembler::assembleForceVector(model);

        // Copy K for BC modification (K_original preserved for reactions)
        Matrix K = K_original;
        BoundaryConditionHandler::apply(K, F, model);

        GaussianSolver solver;
        Vector u = solver.solve(K, F);

        PostProcessor::printResults(model, u);
        PostProcessor::computeReactions(model, u, K_original);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
