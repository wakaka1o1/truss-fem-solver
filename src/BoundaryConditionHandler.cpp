#include "BoundaryConditionHandler.h"

namespace truss_fem {

void BoundaryConditionHandler::apply(Matrix& K, Vector& F, const TrussModel& model) {
    int n = K.size();

    for (const auto& bc : model.bcs) {
        int d = model.dofIndex(bc.nodeId, bc.dof);
        double v = bc.value;

        // Adjust force vector: F[j] -= K[j][d] * v for all j != d
        for (int j = 0; j < n; ++j) {
            if (j != d) {
                F[j] -= K(j, d) * v;
            }
        }

        // Zero row d and column d
        for (int j = 0; j < n; ++j) {
            K(d, j) = 0.0;
            K(j, d) = 0.0;
        }

        // Set diagonal and force
        K(d, d) = 1.0;
        F[d] = v;
    }
}

} // namespace truss_fem
