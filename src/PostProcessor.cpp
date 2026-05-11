#include "PostProcessor.h"
#include "Assembler.h"
#include <iostream>
#include <iomanip>
#include <cmath>

namespace truss_fem {

void PostProcessor::printResults(const TrussModel& model, const Vector& displacements) {
    std::cout << std::scientific << std::setprecision(6);

    // Apply displacements to model nodes for stress computation
    std::vector<Node> nodesWithDisp = model.nodes;
    for (size_t i = 0; i < nodesWithDisp.size(); ++i) {
        nodesWithDisp[i].ux = displacements[2 * i];
        nodesWithDisp[i].uy = displacements[2 * i + 1];
    }

    std::cout << "=== Nodal Displacements ===\n";
    for (const auto& n : nodesWithDisp) {
        std::cout << "Node " << n.id << ": ux = " << std::setw(12) << n.ux
                  << ", uy = " << std::setw(12) << n.uy << "\n";
    }

    std::cout << "\n=== Element Stresses ===\n";
    for (const auto& e : model.elements) {
        const auto& n1 = nodesWithDisp[static_cast<size_t>(e.n1)];
        const auto& n2 = nodesWithDisp[static_cast<size_t>(e.n2)];
        double stress = e.computeStress(n1, n2);
        std::cout << "Element " << e.id << ": sigma = " << std::setw(12) << stress << " Pa\n";
    }
}

Vector PostProcessor::computeReactions(const TrussModel& model, const Vector& displacements,
                                        const Matrix& K) {
    Vector R = K * displacements;

    // Subtract applied loads to isolate reactions
    Vector appliedLoads(model.totalDofs(), 0.0);
    for (const auto& load : model.loads) {
        int idx = model.dofIndex(load.nodeId, load.dof);
        appliedLoads[idx] = load.value;
    }

    int n = model.totalDofs();
    Vector reactions(n, 0.0);
    for (int i = 0; i < n; ++i) {
        reactions[i] = R[i] - appliedLoads[i];
    }

    std::cout << "\n=== Support Reactions ===\n";
    for (const auto& bc : model.bcs) {
        int idx = model.dofIndex(bc.nodeId, bc.dof);
        std::string dofName = (bc.dof == DoF::UX) ? "Fx" : "Fy";
        std::cout << "Node " << bc.nodeId << ": " << dofName
                  << " = " << std::setw(12) << reactions[idx] << "\n";
    }

    return reactions;
}

} // namespace truss_fem
