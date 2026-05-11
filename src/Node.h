#ifndef TRUSS_FEM_NODE_H
#define TRUSS_FEM_NODE_H

namespace truss_fem {

struct Node {
    int id;
    double x{0.0}, y{0.0};
    double ux{0.0}, uy{0.0};
};

} // namespace truss_fem

#endif
