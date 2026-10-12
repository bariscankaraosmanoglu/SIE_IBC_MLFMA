#include "solver_globals.h"
#include "matrixBlock.h"
#include "mlfma_workspace.h"
#include <complex>
using namespace std;

// Simulation-owned arrays. Local owners reset borrowed globals when releasing
// their arrays (direct far fields, topology temporaries, per-illumination RHS).
extern double *node_pos, *rwg_lengths, *tri_areas, *vec_trinormals;
extern double *ffdir_phi, *ffdir_theta, *dir_PW, *size_boxes, *box_pos;
extern int *tri_nodes, *tri_media, *nodetrisindex, *rwg_tris, *rwg_nodes;
extern int *tri2rwgindex, *tri2rwgnode, *tri2rwgsign, *num_truncation, *num_truncation_pc;
extern int *num_filledboxes, *num_children, *children, *children2parent, *rwg_ind;
extern int *boxes_far, *boxes_near, *rwg2box, *rwg_reverseind, *smallestbox2tri;
extern int *trans_index, *trans_index_pc;
extern long int *nearindex, *agg_index, *agg_index_pc;
extern complex<double> *trans_op, *trans_op_pc, *matrix_Znear11, *matrix_Znear12, *matrix_Znear21, *matrix_Znear22, *aggregation, *disaggregation, *vec_farfield, *vec_sol, *vec_RHS;
extern matrixBlock *blockdiagonalmatrix11, *blockdiagonalmatrix22;

namespace {
template<class T> void release(T*& pointer) {
    delete[] pointer;
    pointer = nullptr;
}
}
void clear_memory(){
    // Nulling every owner makes repeated and partially initialized cleanup safe.
    // num_source borrows fixed static storage and must not be deleted.
    mlfma_workspace().clear();
    release(node_pos);
    release(rwg_lengths);
    release(tri_areas);
    release(vec_trinormals);
    release(ffdir_phi);
    release(ffdir_theta);
    release(dir_PW);
    release(size_boxes);
    release(box_pos);
    release(tri_nodes);
    release(tri_media);
    release(nodetrisindex);
    release(rwg_tris);
    release(rwg_nodes);
    release(tri2rwgindex);
    release(tri2rwgnode);
    release(tri2rwgsign);
    release(num_truncation);
    release(num_truncation_pc);
    release(num_filledboxes);
    release(num_children);
    release(children);
    release(children2parent);
    release(rwg_ind);
    release(boxes_far);
    release(boxes_near);
    release(rwg2box);
    release(rwg_reverseind);
    release(smallestbox2tri);
    release(trans_index);
    release(trans_index_pc);
    release(nearindex);
    release(agg_index);
    release(agg_index_pc);
    release(trans_op);
    release(trans_op_pc);
    release(matrix_Znear11);
    release(matrix_Znear12);
    release(matrix_Znear21);
    release(matrix_Znear22);
    release(aggregation);
    release(disaggregation);
    release(vec_farfield);
    release(vec_sol);
    release(vec_RHS);
    release(blockdiagonalmatrix11);
    release(blockdiagonalmatrix22);
}
