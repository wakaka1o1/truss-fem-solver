#ifndef TRUSS_FEM_BOUNDARY_CONDITION_H
#define TRUSS_FEM_BOUNDARY_CONDITION_H

#include "Types.h"

namespace truss_fem {

struct BoundaryCondition {
    int nodeId;
    DoF dof;
    double value;
};

} // namespace truss_fem

#endif
