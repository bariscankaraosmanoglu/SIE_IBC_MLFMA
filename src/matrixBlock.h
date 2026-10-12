#ifndef MATRIXBLOCK_H
#define MATRIXBLOCK_H
#define MLFMA_PIVOTED_BLOCK_LU 1
#include <complex>
#include <vector>

// Owned self-block factors and RWG maps. No symmetry/Hermitian assumption.
class matrixBlock {
public:
    int sizerow = 0;
    int sizecol = 0;
    std::vector<int> indrows, indcols;
    std::vector<std::complex<double>> factors;
    std::vector<int> pivots;

    matrixBlock() = default;
    matrixBlock(const matrixBlock&) = delete;
    matrixBlock& operator=(const matrixBlock&) = delete;
    matrixBlock(matrixBlock&&) = default;
    matrixBlock& operator=(matrixBlock&&) = default;
    void initmatrixBlock(int rows, int columns);
    void assignvalue(int row, int column, std::complex<double> value);
    void factorize();
    void solve_inplace(std::complex<double>* rhs) const;
    bool ready() const { return factorized_; }
private:
    bool factorized_ = false;
};
#endif
