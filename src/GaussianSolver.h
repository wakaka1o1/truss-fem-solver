#ifndef TRUSS_FEM_GAUSSIAN_SOLVER_H
#define TRUSS_FEM_GAUSSIAN_SOLVER_H

#include "Solver.h"
#include <iostream>

namespace truss_fem {

class GaussianSolver : public Solver {
public:
    Vector solve(const Matrix& K, const Vector& F) override;
};

} // namespace truss_fem

#endif
