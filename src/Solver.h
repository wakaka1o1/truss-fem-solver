#ifndef TRUSS_FEM_SOLVER_H
#define TRUSS_FEM_SOLVER_H

#include "Matrix.h"

namespace truss_fem {

class Solver {
public:
    virtual Vector solve(const Matrix& K, const Vector& F) = 0;
    virtual ~Solver() = default;
};

} // namespace truss_fem

#endif
