#ifndef TRUSS_FEM_ELEMENT_H
#define TRUSS_FEM_ELEMENT_H

#include "Matrix.h"
#include "Node.h"

namespace truss_fem {

struct Element {
    int id;
    int n1, n2;
    double E, A;

    double computeLength(const Node& node1, const Node& node2) const;
    Matrix computeLocalStiffness(const Node& node1, const Node& node2) const;
    Matrix computeGlobalStiffness(const Node& node1, const Node& node2) const;
    double computeStrain(const Node& node1, const Node& node2) const;
    double computeStress(const Node& node1, const Node& node2) const;
};

} // namespace truss_fem

#endif
