#include "matrixBlock.h"
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <memory>
#include <stdexcept>
#include <vector>

extern int val_pc;
extern int *children, *boxes_near, *num_filledboxes, *num_children, *rwg_ind;
extern int num_mlfmalevels, num_rwgs;
extern std::complex<double> *matrix_Znear11, *matrix_Znear22;
extern long int *nearindex, size_Znear;
extern matrixBlock *blockdiagonalmatrix11, *blockdiagonalmatrix22;

void calculate_blockdiagonalpreconditioner() {
    delete[] blockdiagonalmatrix11; blockdiagonalmatrix11=nullptr;
    delete[] blockdiagonalmatrix22; blockdiagonalmatrix22=nullptr;
    // Mode 2 uses the approximate-operator inner GMRES without these factors.
    // Keep the guard here as well as in main for direct callers/test drivers.
    if (val_pc != 1) return;
    int box_test = -1;
    try {
        if (num_mlfmalevels < 1 || num_rwgs <= 0)
            throw std::runtime_error("invalid tree or RWG dimensions");
        const int num_boxes = num_filledboxes[num_mlfmalevels-1];
        const int box_start = num_filledboxes[num_mlfmalevels*2-1];
        const int box_end = num_filledboxes[num_mlfmalevels*3-1]+1;
        if (num_boxes <= 0 || box_start < 0 || box_end-box_start != num_boxes)
            throw std::runtime_error("inconsistent leaf-box range");
        std::unique_ptr<matrixBlock[]> blocks11(new matrixBlock[num_boxes]);
        std::unique_ptr<matrixBlock[]> blocks22(new matrixBlock[num_boxes]);
        std::vector<unsigned char> covered(num_rwgs, 0);
        long long coefficients = 0;
        for (box_test = box_start; box_test < box_end; ++box_test) {
            const int n = num_children[box_test];
            const int first = children[box_test*2];
            const int last = children[box_test*2+1];
            if (n <= 0 || first < 0 || last < first || last >= num_rwgs || last-first+1 != n)
                throw std::runtime_error("invalid leaf RWG range");
            int self = -1;
            for (int neighbor = 0; neighbor < 27; ++neighbor)
                if (boxes_near[box_test*27+neighbor] == box_test) {
                    if (self != -1) throw std::runtime_error("duplicate self interaction");
                    self = neighbor;
                }
            if (self == -1) throw std::runtime_error("missing self interaction");
            const long int offset = nearindex[(box_test-box_start)*27+self];
            const long long count = 1LL*n*n;
            if (offset < 0 || offset > size_Znear || count > size_Znear-offset)
                throw std::runtime_error("self-block offset exceeds near matrix");
            for (int row = 0; row < n; ++row) {
                const int rwg = rwg_ind[first+row];
                if (rwg < 0 || rwg >= num_rwgs || covered[rwg])
                    throw std::runtime_error("invalid or repeated RWG index");
                covered[rwg] = 1;
            }
            // Keep the legacy block-diagonal method on each current separately.
            // Factoring a combined 2n-by-2n block would change that preconditioner.
            for (int current = 0; current < 2; ++current) {
                matrixBlock& block=(current == 0 ? blocks11 : blocks22)[box_test-box_start];
                const auto* matrix=current == 0 ? matrix_Znear11 : matrix_Znear22;
                block.initmatrixBlock(n,n);
                for (int row = 0; row < n; ++row) {
                    block.indrows[row]=block.indcols[row]=rwg_ind[first+row]+current*num_rwgs;
                    for (int column = 0; column < n; ++column)
                        block.assignvalue(row,column,matrix[offset+1LL*row*n+column]);
                }
                block.factorize();
            }
            coefficients += count;
        }
        for (const auto value : covered)
            if (!value) throw std::runtime_error("leaf blocks do not cover every RWG");
        blockdiagonalmatrix11=blocks11.release();
        blockdiagonalmatrix22=blocks22.release();
        std::printf("Block preconditioner: %d pivoted LU blocks, %lld complex coefficients\n",2*num_boxes,2*coefficients);
    } catch (const std::exception& error) {
        std::fprintf(stderr,"Block preconditioner error at leaf box %d: %s\n",box_test,error.what());
        std::exit(EXIT_FAILURE);
    }
}
