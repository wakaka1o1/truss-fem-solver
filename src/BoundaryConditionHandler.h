#ifndef TRUSS_FEM_BC_HANDLER_H
#define TRUSS_FEM_BC_HANDLER_H

#include "Matrix.h"
#include "TrussModel.h"

namespace truss_fem {

class BoundaryConditionHandler {
public:
    static void apply(Matrix& K, Vector& F, const TrussModel& model);
};

} // namespace truss_fem

#endif
