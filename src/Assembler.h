#ifndef TRUSS_FEM_ASSEMBLER_H
#define TRUSS_FEM_ASSEMBLER_H

#include "Matrix.h"
#include "TrussModel.h"

namespace truss_fem {

class Assembler {
public:
    static Matrix assembleStiffness(const TrussModel& model);
    static Vector assembleForceVector(const TrussModel& model);
};

} // namespace truss_fem

#endif
