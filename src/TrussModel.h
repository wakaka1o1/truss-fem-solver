#ifndef TRUSS_FEM_TRUSS_MODEL_H
#define TRUSS_FEM_TRUSS_MODEL_H

#include <vector>
#include <stdexcept>
#include "Node.h"
#include "Element.h"
#include "BoundaryCondition.h"
#include "NodalLoad.h"
#include "Types.h"

namespace truss_fem {

class TrussModel {
public:
    std::vector<Node> nodes;
    std::vector<Element> elements;
    std::vector<BoundaryCondition> bcs;
    std::vector<NodalLoad> loads;

    int totalDofs() const;
    int dofIndex(int nodeId, DoF dof) const;
    const Node& getNode(int id) const;
};

} // namespace truss_fem

#endif
