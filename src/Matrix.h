#ifndef TRUSS_FEM_MATRIX_H
#define TRUSS_FEM_MATRIX_H

#include <vector>
#include <cassert>

namespace truss_fem {

using Vector = std::vector<double>;

class Matrix {
public:
    explicit Matrix(int n) : data_(n * n, 0.0), n_(n) {
        assert(n > 0);
    }

    double& operator()(int i, int j) {
        assert(i >= 0 && i < n_ && j >= 0 && j < n_);
        return data_[static_cast<size_t>(i) * n_ + j];
    }

    double operator()(int i, int j) const {
        assert(i >= 0 && i < n_ && j >= 0 && j < n_);
        return data_[static_cast<size_t>(i) * n_ + j];
    }

    int size() const { return n_; }

    Vector operator*(const Vector& v) const {
        assert(static_cast<int>(v.size()) == n_);
        Vector result(n_, 0.0);
        for (int i = 0; i < n_; ++i) {
            for (int j = 0; j < n_; ++j) {
                result[i] += (*this)(i, j) * v[j];
            }
        }
        return result;
    }

private:
    std::vector<double> data_;
    int n_;
};

} // namespace truss_fem

#endif
