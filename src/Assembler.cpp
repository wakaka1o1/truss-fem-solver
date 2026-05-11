#include "Assembler.h"

namespace truss_fem {

Matrix Assembler::assembleStiffness(const TrussModel& model) {
    int n = model.totalDofs();
    Matrix K(n);

    for (const auto& elem : model.elements) {
        const Node& node1 = model.getNode(elem.n1);
        const Node& node2 = model.getNode(elem.n2);
        Matrix ke = elem.computeGlobalStiffness(node1, node2);

        int dofs[4] = {
            model.dofIndex(elem.n1, DoF::UX),
            model.dofIndex(elem.n1, DoF::UY),
            model.dofIndex(elem.n2, DoF::UX),
            model.dofIndex(elem.n2, DoF::UY)
        };

        for (int a = 0; a < 4; ++a) {
            for (int b = 0; b < 4; ++b) {
                K(dofs[a], dofs[b]) += ke(a, b);
            }
        }
    }

    return K;
}

Vector Assembler::assembleForceVector(const TrussModel& model) {
    int n = model.totalDofs();
    Vector F(n, 0.0);

    for (const auto& load : model.loads) {
        int idx = model.dofIndex(load.nodeId, load.dof);
        F[idx] += load.value;
    }

    return F;
}

} // namespace truss_fem
