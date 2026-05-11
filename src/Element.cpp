#include "Element.h"
#include <cmath>
#include <stdexcept>
#include <string>

namespace truss_fem {

double Element::computeLength(const Node& node1, const Node& node2) const {
    double dx = node2.x - node1.x;
    double dy = node2.y - node1.y;
    return std::sqrt(dx * dx + dy * dy);
}

Matrix Element::computeLocalStiffness(const Node& node1, const Node& node2) const {
    double L = computeLength(node1, node2);
    if (L < 1e-15) {
        throw std::runtime_error("Element " + std::to_string(id) + " has zero length");
    }
    double k = E * A / L;
    Matrix ke(2);
    ke(0, 0) =  k;   ke(0, 1) = -k;
    ke(1, 0) = -k;   ke(1, 1) =  k;
    return ke;
}

Matrix Element::computeGlobalStiffness(const Node& node1, const Node& node2) const {
    double dx = node2.x - node1.x;
    double dy = node2.y - node1.y;
    double L = std::sqrt(dx * dx + dy * dy);

    if (L < 1e-15) {
        throw std::runtime_error("Element " + std::to_string(id) + " has zero length");
    }

    double c = dx / L;
    double s = dy / L;
    double k = E * A / L;

    Matrix ke(4);
    ke(0, 0) = k * c * c;   ke(0, 1) = k * c * s;   ke(0, 2) = -k * c * c;  ke(0, 3) = -k * c * s;
    ke(1, 0) = k * c * s;   ke(1, 1) = k * s * s;   ke(1, 2) = -k * c * s;  ke(1, 3) = -k * s * s;
    ke(2, 0) = -k * c * c;  ke(2, 1) = -k * c * s;  ke(2, 2) = k * c * c;   ke(2, 3) = k * c * s;
    ke(3, 0) = -k * c * s;  ke(3, 1) = -k * s * s;  ke(3, 2) = k * c * s;   ke(3, 3) = k * s * s;

    return ke;
}

double Element::computeStrain(const Node& node1, const Node& node2) const {
    double dx = node2.x - node1.x;
    double dy = node2.y - node1.y;
    double L0 = std::sqrt(dx * dx + dy * dy);

    double dx_new = (node2.x + node2.ux) - (node1.x + node1.ux);
    double dy_new = (node2.y + node2.uy) - (node1.y + node1.uy);
    double L_new = std::sqrt(dx_new * dx_new + dy_new * dy_new);

    return (L_new - L0) / L0;
}

double Element::computeStress(const Node& node1, const Node& node2) const {
    return E * computeStrain(node1, node2);
}

} // namespace truss_fem
