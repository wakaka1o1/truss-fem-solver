#ifndef TRUSS_FEM_POST_PROCESSOR_H
#define TRUSS_FEM_POST_PROCESSOR_H

#include "TrussModel.h"
#include "Matrix.h"

namespace truss_fem {

class PostProcessor {
public:
    static void printResults(const TrussModel& model, const Vector& displacements);
    static Vector computeReactions(const TrussModel& model, const Vector& displacements,
                                   const Matrix& K);
};

} // namespace truss_fem

#endif
