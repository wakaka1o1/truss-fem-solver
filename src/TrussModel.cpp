#include "TrussModel.h"
#include <string>

namespace truss_fem {

int TrussModel::totalDofs() const {
    return 2 * static_cast<int>(nodes.size());
}

int TrussModel::dofIndex(int nodeId, DoF dof) const {
    return 2 * nodeId + (dof == DoF::UY ? 1 : 0);
}

const Node& TrussModel::getNode(int id) const {
    for (const auto& n : nodes) {
        if (n.id == id) return n;
    }
    throw std::runtime_error("Node " + std::to_string(id) + " not found");
}

} // namespace truss_fem
