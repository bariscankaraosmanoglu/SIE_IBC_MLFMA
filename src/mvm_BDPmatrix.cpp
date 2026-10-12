#include "matrixBlock.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <stdexcept>
#include <vector>

extern matrixBlock *blockdiagonalmatrix11, *blockdiagonalmatrix22;
extern int val_pc;
extern int *num_filledboxes;
extern int num_mlfmalevels;

void mvm_BDPmatrix(std::complex<double> x[], std::complex<double> y[]) {
    try {
        if (val_pc != 1 || !blockdiagonalmatrix11 || !blockdiagonalmatrix22 || !x || !y)
            throw std::runtime_error("mode-1 factors or application vectors are unavailable");
        const int num_boxes = num_filledboxes[num_mlfmalevels-1];
        int largest = 0;
        for (int box = 0; box < num_boxes; ++box)
            largest = std::max({largest,blockdiagonalmatrix11[box].sizerow,blockdiagonalmatrix22[box].sizerow});
        // Per-call workspace keeps applications independent. Disjoint RWG
        // partitions permit x == y: gather a complete block before writing it.
        std::vector<std::complex<double>> rhs(largest);
        for (int current = 0; current < 2; ++current) {
        for (int box = 0; box < num_boxes; ++box) {
            const matrixBlock& block = (current == 0 ? blockdiagonalmatrix11 : blockdiagonalmatrix22)[box];
            for (int row = 0; row < block.sizerow; ++row)
                rhs[row] = x[block.indrows[row]];
            block.solve_inplace(rhs.data());
            for (int column = 0; column < block.sizecol; ++column)
                y[block.indcols[column]] = rhs[column];
        }
        }
    } catch (const std::exception& error) {
        std::fprintf(stderr,"Block preconditioner application error: %s\n",error.what());
        std::exit(EXIT_FAILURE);
    }
}
