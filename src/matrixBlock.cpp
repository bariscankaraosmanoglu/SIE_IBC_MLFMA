#include "matrixBlock.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace {
bool finite(const std::complex<double>& value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}
}

void matrixBlock::initmatrixBlock(int rows, int columns) {
    if (rows <= 0 || rows != columns)
        throw std::runtime_error("expected a nonempty square self block");
    const std::size_t n = rows;
    if (n > std::numeric_limits<std::size_t>::max()/n ||
        n*n > factors.max_size())
        throw std::runtime_error("block dimensions exceed storage capacity");
    matrixBlock next;
    next.sizerow = rows; next.sizecol = columns;
    next.indrows.resize(n); next.indcols.resize(n);
    next.factors.assign(n*n, {}); next.pivots.resize(n);
    *this = std::move(next);
}

void matrixBlock::assignvalue(int row, int column, std::complex<double> value) {
    factors[std::size_t(row)*sizecol+column] = value;
    factorized_ = false;
}

void matrixBlock::factorize() {
    factorized_ = false;
    const int n = sizerow;
    if (n <= 0 || sizecol != n || factors.size() != std::size_t(n)*n)
        throw std::runtime_error("uninitialized or invalid block storage");
    for (const auto& value : factors)
        if (!finite(value)) throw std::runtime_error("nonfinite block coefficient");
    for (int k = 0; k < n; ++k) {
        int pivot = k;
        double largest = std::abs(factors[std::size_t(k)*n+k]);
        for (int row = k+1; row < n; ++row) {
            const double candidate = std::abs(factors[std::size_t(row)*n+k]);
            if (candidate > largest) { largest = candidate; pivot = row; }
        }
        // No absolute epsilon cutoff: a uniformly small, nonsingular block
        // should receive the same pivot decisions as its rescaled counterpart.
        if (largest == 0.0)
            throw std::runtime_error("singular block at pivot " + std::to_string(k));
        if (!std::isfinite(largest))
            throw std::runtime_error("unrepresentable pivot magnitude at column " + std::to_string(k));
        pivots[k] = pivot;
        if (pivot != k)
            for (int column = 0; column < n; ++column)
                std::swap(factors[std::size_t(k)*n+column], factors[std::size_t(pivot)*n+column]);
        for (int row = k+1; row < n; ++row) {
            auto& multiplier = factors[std::size_t(row)*n+k];
            multiplier /= factors[std::size_t(k)*n+k];
            if (!finite(multiplier))
                throw std::runtime_error("nonfinite multiplier at pivot " + std::to_string(k));
            for (int column = k+1; column < n; ++column) {
                auto& value = factors[std::size_t(row)*n+column];
                value -= multiplier*factors[std::size_t(k)*n+column];
                if (!finite(value))
                    throw std::runtime_error("nonfinite elimination at pivot " + std::to_string(k));
            }
        }
    }
    factorized_ = true;
}

void matrixBlock::solve_inplace(std::complex<double>* rhs) const {
    if (!factorized_ || !rhs) throw std::runtime_error("block factors are not ready for application");
    const int n = sizerow;
    // P*B = L*U: apply the same sequence of row exchanges to the RHS.
    for (int k = 0; k < n; ++k)
        if (pivots[k] != k) std::swap(rhs[k], rhs[pivots[k]]);
    for (int row = 0; row < n; ++row)
        for (int column = 0; column < row; ++column)
            rhs[row] -= factors[std::size_t(row)*n+column]*rhs[column];
    for (int row = n-1; row >= 0; --row) {
        for (int column = row+1; column < n; ++column)
            rhs[row] -= factors[std::size_t(row)*n+column]*rhs[column];
        rhs[row] /= factors[std::size_t(row)*n+row];
    }
    // Arithmetic overflow in an application propagates to the existing Krylov
    // finite-value checks; it must not bypass their termination/status reporting.
}
