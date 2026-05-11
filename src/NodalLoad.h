#ifndef TRUSS_FEM_NODAL_LOAD_H
#define TRUSS_FEM_NODAL_LOAD_H

#include "Types.h"

namespace truss_fem {

struct NodalLoad {
    int nodeId;
    DoF dof;
    double value;
};

} // namespace truss_fem

#endif
