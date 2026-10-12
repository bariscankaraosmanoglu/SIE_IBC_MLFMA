#include "mlfma_workspace.h"
#include "solver_globals.h"
#include <stdio.h>
#include <complex>
#include <chrono>
#include <omp.h>
#include <cmath>

using namespace std;
using namespace std::chrono;

extern double pi;
extern long int *agg_index;
extern long int size_agg;
extern double *box_pos;
extern double *size_boxes;
extern int *children;
extern double val_kO;
extern int *boxes_far;
extern int *num_truncation;
extern int num_interpol;
extern int num_mlfmalevels;
extern int *num_filledboxes;
extern int num_rwgs;
extern int *rwg_ind;
extern double *rwg_lengths;
extern double *node_pos;
extern int *rwg_nodes;
extern int *rwg_tris;
extern double *vec_trinormals;
extern int *trans_index;
extern complex<double> *trans_op;
extern complex<double> *aggregation;
extern complex<double> *disaggregation;
extern int parallel_omp_on;
extern int num_threads;
extern FILE *file_timer;

void phi_locations_weights(int,double*,double*);
void theta_locations_weights(int,double*,double*);
void calculate_sin_cos(int, int, double*, double*, double*, double*, double*, double*);
void interpolationweights(int,int,int,int,double*,double*,double*,int*);
void interpolation_c2f(int,int,int,int,int,int*,int*,double*,double*,complex<double>*,complex<double>*);
void interpolation_f2c(int,int,int,int,int,int*,int*,double*,double*,complex<double>*,complex<double>*);
void crossProduct3d(double vec1[], double vec2[], double resvec[]);
double dotProduct3d(double vec1[], double vec2[]);

void mvm_Zfar1121dotx(complex<double> x[], complex<double> y[]){
    auto& workspace = mlfma_workspace();
    //
    //
    // mvm_Zfardotx
    // Calculate the far-zone interactions via MLFMA.
    // Zfar dot x operation is performed.
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    int box_start, box_end;
    complex<double> complex0=complex<double>(0.0, 0.0);
    complex<double> i=complex<double>(0.0, 1.0);
    complex<double> y_constant=(i*val_kO/(4.0*pi))*(i*val_kO/(4.0*pi));

    for(int ind0=0; ind0<(size_agg*2);ind0++){
        aggregation[ind0]=complex0;
    }
    for(int ind0=0; ind0<(size_agg*2);ind0++){
        disaggregation[ind0]=complex0;
    }


    //
    //
    // Aggregation operation at the last level
    // Aggregation is performed on every basis function at each box.
    // Radiation patterns of each basis function is shifted from sample position to the box center.
    // Shifted patterns are multiplied with the given complex coefficients.
    box_start=num_filledboxes[num_mlfmalevels*2-1];
    box_end=num_filledboxes[num_mlfmalevels*3-1]+1;
    const auto& leaf_grid = workspace.grid(num_truncation[num_mlfmalevels-1]+1,2*(num_truncation[num_mlfmalevels-1]+1));
    omp_set_num_threads(num_threads);

    #pragma omp parallel if(parallel_omp_on) shared(aggregation,x)
    {
        double shiftvector[3];
        double khat[3]{0.0, 0.0, 0.0};
        double dotval;
        complex<double> expfactor;
        int numberofterms=num_truncation[num_mlfmalevels-1];
        int numberofthetas=numberofterms+1;
        int numberofphis=numberofthetas*2;
        const double* sintheta = leaf_grid.sin_theta.data();
        const double* costheta = leaf_grid.cos_theta.data();
        const double* sinphi = leaf_grid.sin_phi.data();
        const double* cosphi = leaf_grid.cos_phi.data();
        int basisnode1,basisnode2,basisnode3;
        double basisnode1loc[3],basisnode2loc[3],basisnode3loc[3],basislocation[3];
        double basislength,basissign,basisvector[3];
        double basisradiationtheta,basisradiationphi;
        double rwgsign[2]{1.0,-1.0};
        int rwglastnodeindex[2]{2,3};
        int indtheta,indphi;
        complex<double> val_x;


        int tid=omp_get_thread_num();
        double timer_start = omp_get_wtime();
        double timer_total = 0.0;

        #pragma omp for schedule(dynamic,5)
        for(int box_basis=box_start;box_basis<box_end;box_basis++){
            double timer_in = omp_get_wtime();
            int start_basis=children[box_basis*2+0];
            int end_basis=children[box_basis*2+1]+1;
            for(int ind_basis=start_basis;ind_basis<end_basis;ind_basis++){
                int basisfunction=rwg_ind[ind_basis];
                val_x=x[basisfunction];
                double basisboxcenter[3]={box_pos[box_basis*3+0],box_pos[box_basis*3+1],box_pos[box_basis*3+2]};
                for(int rwgtri=0; rwgtri<2; rwgtri++){
                    basisnode1=rwg_nodes[basisfunction*4+0];
                    basisnode2=rwg_nodes[basisfunction*4+1];
                    basisnode3=rwg_nodes[basisfunction*4+rwglastnodeindex[rwgtri]];
                    basisnode1loc[0]=node_pos[basisnode1*3+0];
                    basisnode1loc[1]=node_pos[basisnode1*3+1];
                    basisnode1loc[2]=node_pos[basisnode1*3+2];
                    basisnode2loc[0]=node_pos[basisnode2*3+0];
                    basisnode2loc[1]=node_pos[basisnode2*3+1];
                    basisnode2loc[2]=node_pos[basisnode2*3+2];
                    basisnode3loc[0]=node_pos[basisnode3*3+0];
                    basisnode3loc[1]=node_pos[basisnode3*3+1];
                    basisnode3loc[2]=node_pos[basisnode3*3+2];
                    basislocation[0]=(basisnode1loc[0]+basisnode2loc[0]+basisnode3loc[0])/3;
                    basislocation[1]=(basisnode1loc[1]+basisnode2loc[1]+basisnode3loc[1])/3;
                    basislocation[2]=(basisnode1loc[2]+basisnode2loc[2]+basisnode3loc[2])/3;
                    basislength=rwg_lengths[basisfunction];
                    basissign=rwgsign[rwgtri];
                    shiftvector[0]=basislocation[0]-basisboxcenter[0];
                    shiftvector[1]=basislocation[1]-basisboxcenter[1];
                    shiftvector[2]=basislocation[2]-basisboxcenter[2];
                    basisvector[0]=basislocation[0]-basisnode3loc[0];
                    basisvector[1]=basislocation[1]-basisnode3loc[1];
                    basisvector[2]=basislocation[2]-basisnode3loc[2];
                    long int indagg=agg_index[box_basis]*2;
                    for(indtheta=0;indtheta<numberofthetas;indtheta++){
                        for(indphi=0;indphi<numberofphis;indphi++){
                            khat[0]=sintheta[indtheta]*cosphi[indphi];
                            khat[1]=sintheta[indtheta]*sinphi[indphi];
                            khat[2]=costheta[indtheta];
                            dotval=shiftvector[0]*khat[0]+shiftvector[1]*khat[1]+shiftvector[2]*khat[2];
                            expfactor=exp(-i*val_kO*dotval);
                            basisradiationtheta=costheta[indtheta]*cosphi[indphi]*basisvector[0]+costheta[indtheta]*sinphi[indphi]*basisvector[1]-sintheta[indtheta]*basisvector[2];
                            basisradiationphi=-sinphi[indphi]*basisvector[0]+cosphi[indphi]*basisvector[1];
                            aggregation[indagg]+=val_x*basissign*basislength*expfactor*basisradiationtheta/2.0;
                            indagg++;
                            aggregation[indagg]+=val_x*basissign*basislength*expfactor*basisradiationphi/2.0;
                            indagg++;
                        }
                    }
                }
            }
            timer_total+=omp_get_wtime()-timer_in;
        }
        double timer_end = omp_get_wtime();
//        fprintf(file_timer,"MVMFarAgg%d %d %1.5f %1.5f\n", num_mlfmalevels-1, omp_get_thread_num(), timer_end-timer_start, timer_total );
    }
    
    


    //
    //
    // Aggregation operation at the higher levels
    // Aggregation is performed on every child of each parent box.
    // Fields of the boxes are interpolated (upsampled) and shifted to the parent box center
    for(int indlev=(num_mlfmalevels-2); indlev>1; indlev--){
        int numberofterms=num_truncation[indlev];
        int num_finesample_theta=numberofterms+1;
        int num_finesample_phi=num_finesample_theta*2;
        numberofterms=num_truncation[indlev+1];
        int num_coarsesample_theta=numberofterms+1;
        int num_coarsesample_phi=num_coarsesample_theta*2;
        auto& fine_grid = workspace.grid(num_finesample_theta,num_finesample_phi);
        const auto& phases = workspace.child_shifts(num_finesample_theta,num_finesample_phi,
                                                    size_boxes[indlev],val_kO);
        auto& weights = workspace.table(num_coarsesample_theta,num_coarsesample_phi,
                                        num_finesample_theta,num_finesample_phi,num_interpol);
        double* sintheta = fine_grid.sin_theta.data();
        double* costheta = fine_grid.cos_theta.data();
        double* sinphi = fine_grid.sin_phi.data();
        double* cosphi = fine_grid.cos_phi.data();
        double* inter_weights_theta = weights.theta_weights.data();
        double* inter_weights_phi = weights.phi_weights.data();
        int* indexloc_theta = weights.theta_indices.data();
        int* indexloc_phi = weights.phi_indices.data();
        complex<double> *temp_fine_theta, *temp_fine_phi;
        complex<double> *temp_coarse_theta, *temp_coarse_phi;

        box_start=num_filledboxes[num_mlfmalevels*1+indlev];
        box_end=num_filledboxes[num_mlfmalevels*2+indlev]+1;

        omp_set_num_threads(num_threads);
        #pragma omp parallel if(parallel_omp_on) shared(aggregation) private(temp_fine_theta,temp_fine_phi,\
                                                                             temp_coarse_theta,temp_coarse_phi)
        {
            workspace.prepare_team();
            double shiftvector[3];
            double khat[3]{0.0, 0.0, 0.0};
            double dotval;
            complex<double> expfactor;
            double timer_start = omp_get_wtime();
            double timer_total = 0.0;
            auto& scratch = workspace.worker();
            scratch.prepare(num_coarsesample_theta,num_coarsesample_phi,num_finesample_theta,num_finesample_phi);
            temp_fine_theta = scratch.fine_theta.data();
            temp_fine_phi = scratch.fine_phi.data();
            temp_coarse_theta = scratch.coarse_theta.data();
            temp_coarse_phi = scratch.coarse_phi.data();
            #pragma omp for schedule(dynamic)
            for(int box_parent=box_start;box_parent<box_end;box_parent++){
                double timer_in = omp_get_wtime();
                double parentbasisboxcenter[3]={box_pos[box_parent*3+0],box_pos[box_parent*3+1],box_pos[box_parent*3+2]};
                int box_child_start=children[box_parent*2+0];
                int box_child_end=children[box_parent*2+1]+1;
                for(int box_child=box_child_start;box_child<box_child_end;box_child++){
                    double childbasisboxcenter[3]={box_pos[box_child*3+0],box_pos[box_child*3+1],box_pos[box_child*3+2]};
                    shiftvector[0]=childbasisboxcenter[0]-parentbasisboxcenter[0];
                    shiftvector[1]=childbasisboxcenter[1]-parentbasisboxcenter[1];
                    shiftvector[2]=childbasisboxcenter[2]-parentbasisboxcenter[2];
                    const complex<double>* shift_phase = phases.for_shift(shiftvector);
                    long int ind_child=agg_index[box_child];
                    int ind_temp=0;
                    for(int ind=0;ind<(num_finesample_theta*num_finesample_phi);ind++){
                        temp_fine_theta[ind]=complex0;
                        temp_fine_phi[ind]=complex0;
                    }
                    for(int ind=0;ind<(num_coarsesample_theta*num_coarsesample_phi);ind++){
                        temp_coarse_theta[ind]=complex0;
                        temp_coarse_phi[ind]=complex0;
                    }
                    for(int indtheta=0;indtheta<num_coarsesample_theta;indtheta++){
                        for(int indphi=0;indphi<num_coarsesample_phi;indphi++){
                            temp_coarse_theta[ind_temp]=aggregation[ind_child*2+0];
                            temp_coarse_phi[ind_temp]=aggregation[ind_child*2+1];
                            ind_temp++;
                            ind_child++;
                        }
                    }
                    interpolation_c2f(num_coarsesample_theta,num_coarsesample_phi,num_finesample_theta,num_finesample_phi,num_interpol,indexloc_theta,indexloc_phi,inter_weights_theta,inter_weights_phi,temp_coarse_theta,temp_fine_theta,scratch.interpolation);
                    interpolation_c2f(num_coarsesample_theta,num_coarsesample_phi,num_finesample_theta,num_finesample_phi,num_interpol,indexloc_theta,indexloc_phi,inter_weights_theta,inter_weights_phi,temp_coarse_phi,temp_fine_phi,scratch.interpolation);
                    long int ind_parent=agg_index[box_parent];

                    ind_temp=0;
                    for(int indtheta=0;indtheta<num_finesample_theta;indtheta++){
                        for(int indphi=0;indphi<num_finesample_phi;indphi++){
                            if (shift_phase) expfactor=shift_phase[ind_temp];
                            else {
                                khat[0]=sintheta[indtheta]*cosphi[indphi];
                                khat[1]=sintheta[indtheta]*sinphi[indphi];
                                khat[2]=costheta[indtheta];
                                dotval=dotProduct3d(shiftvector,khat);
                                expfactor=exp(-i*val_kO*dotval);
                            }
                            aggregation[ind_parent*2+0]+=expfactor*temp_fine_theta[ind_temp];
                            aggregation[ind_parent*2+1]+=expfactor*temp_fine_phi[ind_temp];
                            ind_temp++;
                            ind_parent++;
                        }
                    }
                }

                timer_total+=omp_get_wtime()-timer_in;
            }

            double timer_end = omp_get_wtime();
//            fprintf(file_timer,"MVMFarAgg%d %d %1.5f %1.5f \n", indlev, omp_get_thread_num(), timer_end-timer_start, timer_total );
        }
    }



    //
    //
    // Translation operation at every level
    // Far boxes of each basis box is considered.
    // Basis to test box translation is performed.
    for(int indlev=2; indlev<num_mlfmalevels; indlev++){
        int numberofterms=num_truncation[indlev];
        int numberofthetas=numberofterms+1;
        int numberofphis=numberofthetas*2;
        box_start=num_filledboxes[num_mlfmalevels*1+indlev];
        box_end=num_filledboxes[num_mlfmalevels*2+indlev]+1;

        omp_set_num_threads(num_threads);
        #pragma omp parallel if(parallel_omp_on)
        {
            double timer_start = omp_get_wtime();
            double timer_total = 0.0;
            #pragma omp for schedule(dynamic)
            for(int box_test=box_start;box_test<box_end;box_test++){
                double timer_in = omp_get_wtime();
                double testboxcenter[3]={box_pos[box_test*3+0],box_pos[box_test*3+1],box_pos[box_test*3+2]};
                for(int ind_farbox=0;ind_farbox<189;ind_farbox++){
                    int box_basis=boxes_far[box_test*189+ind_farbox];
                    if(box_basis!=-1){
                        double basisboxcenter[3]={box_pos[box_basis*3+0],box_pos[box_basis*3+1],box_pos[box_basis*3+2]};
                        int indx=lround(3.0+(testboxcenter[0]-basisboxcenter[0])/size_boxes[indlev]);
                        int indy=lround(3.0+(testboxcenter[1]-basisboxcenter[1])/size_boxes[indlev]);
                        int indz=lround(3.0+(testboxcenter[2]-basisboxcenter[2])/size_boxes[indlev]);
                        int indexval=indlev*343+indx*49+indy*7+indz;
                        int ind_trans_op=trans_index[indexval];
                        long int indagg=agg_index[box_basis];
                        long int inddisagg=agg_index[box_test];
                        for(int indtheta=0;indtheta<numberofthetas;indtheta++){
                            for(int indphi=0;indphi<numberofphis;indphi++){
                                disaggregation[inddisagg*2+0]+=trans_op[ind_trans_op]*aggregation[indagg*2+0];
                                disaggregation[inddisagg*2+1]+=trans_op[ind_trans_op]*aggregation[indagg*2+1];
                                ind_trans_op++;
                                indagg++;
                                inddisagg++;
                            }
                        }
                    }
                }
                timer_total+=omp_get_wtime()-timer_in;
            }
            double timer_end = omp_get_wtime();
//            fprintf(file_timer,"MVMFarTra%d %d %1.5f %1.5f \n", indlev, omp_get_thread_num(), timer_end-timer_start, timer_total );
        }
    }



    //
    //
    // Disaggregation operation at the higher levels
    // Disaggregation is performed in the direction of parent box to child boxes.
    // Fields of the parent boxes are shifted and interpolated (downsampled) to the child box center.
    for(int indlev=2; indlev<(num_mlfmalevels-1); indlev++){
        int numberofterms=num_truncation[indlev];
        int num_finesample_theta=numberofterms+1;
        int num_finesample_phi=num_finesample_theta*2;
        numberofterms=num_truncation[indlev+1];
        int num_coarsesample_theta=numberofterms+1;
        int num_coarsesample_phi=num_coarsesample_theta*2;
        auto& fine_grid = workspace.grid(num_finesample_theta,num_finesample_phi);
        const auto& phases = workspace.child_shifts(num_finesample_theta,num_finesample_phi,
                                                    size_boxes[indlev],val_kO);
        auto& weights = workspace.table(num_coarsesample_theta,num_coarsesample_phi,
                                        num_finesample_theta,num_finesample_phi,num_interpol);
        double* sintheta = fine_grid.sin_theta.data();
        double* costheta = fine_grid.cos_theta.data();
        double* sinphi = fine_grid.sin_phi.data();
        double* cosphi = fine_grid.cos_phi.data();
        double* inter_weights_theta = weights.theta_weights.data();
        double* inter_weights_phi = weights.phi_weights.data();
        int* indexloc_theta = weights.theta_indices.data();
        int* indexloc_phi = weights.phi_indices.data();
        complex<double> *temp_fine_theta, *temp_fine_phi;
        complex<double> *temp_coarse_theta, *temp_coarse_phi;

        box_start=num_filledboxes[num_mlfmalevels*1+indlev];
        box_end=num_filledboxes[num_mlfmalevels*2+indlev]+1;

        omp_set_num_threads(num_threads);
        #pragma omp parallel if(parallel_omp_on) shared(disaggregation) private(temp_fine_theta,temp_fine_phi,\
                                                                                temp_coarse_theta,temp_coarse_phi)
        {
            workspace.prepare_team();
            double shiftvector[3];
            double khat[3]{0.0, 0.0, 0.0};
            double dotval;
            complex<double> expfactor;

            double timer_start = omp_get_wtime();
            double timer_total = 0.0;
            auto& scratch = workspace.worker();
            scratch.prepare(num_coarsesample_theta,num_coarsesample_phi,num_finesample_theta,num_finesample_phi);
            temp_fine_theta = scratch.fine_theta.data();
            temp_fine_phi = scratch.fine_phi.data();
            temp_coarse_theta = scratch.coarse_theta.data();
            temp_coarse_phi = scratch.coarse_phi.data();
            #pragma omp for  schedule(dynamic)
            for(int box_parent=box_start;box_parent<box_end;box_parent++){
                double timer_in = omp_get_wtime();
                double parenttestboxcenter[3]={box_pos[box_parent*3+0],box_pos[box_parent*3+1],box_pos[box_parent*3+2]};
                int box_child_start=children[box_parent*2+0];
                int box_child_end=children[box_parent*2+1]+1;
                for(int box_child=box_child_start;box_child<box_child_end;box_child++){
                    double childtestboxcenter[3]={box_pos[box_child*3+0],box_pos[box_child*3+1],box_pos[box_child*3+2]};
                    shiftvector[0]=childtestboxcenter[0]-parenttestboxcenter[0];
                    shiftvector[1]=childtestboxcenter[1]-parenttestboxcenter[1];
                    shiftvector[2]=childtestboxcenter[2]-parenttestboxcenter[2];
                    const complex<double>* shift_phase = phases.for_shift(shiftvector);
                    long int ind_parent=agg_index[box_parent];

                    int ind_temp=0;
                    for(int ind=0;ind<(num_finesample_theta*num_finesample_phi);ind++){
                        temp_fine_theta[ind]=complex0;
                        temp_fine_phi[ind]=complex0;
                    }
                    for(int ind=0;ind<(num_coarsesample_theta*num_coarsesample_phi);ind++){
                        temp_coarse_theta[ind]=complex0;
                        temp_coarse_phi[ind]=complex0;
                    }
                    for(int indtheta=0;indtheta<num_finesample_theta;indtheta++){
                        for(int indphi=0;indphi<num_finesample_phi;indphi++){
                            if (shift_phase) expfactor=conj(shift_phase[ind_temp]);
                            else {
                                khat[0]=sintheta[indtheta]*cosphi[indphi];
                                khat[1]=sintheta[indtheta]*sinphi[indphi];
                                khat[2]=costheta[indtheta];
                                dotval=dotProduct3d(shiftvector,khat);
                                expfactor=exp(i*val_kO*dotval);
                            }
                            temp_fine_theta[ind_temp]=expfactor*disaggregation[ind_parent*2+0];
                            temp_fine_phi[ind_temp]=expfactor*disaggregation[ind_parent*2+1];
                            ind_temp++;
                            ind_parent++;
                        }
                    }
                    interpolation_f2c(num_coarsesample_theta,num_coarsesample_phi,num_finesample_theta,num_finesample_phi,num_interpol,indexloc_theta,indexloc_phi,inter_weights_theta,inter_weights_phi,temp_fine_theta,temp_coarse_theta,scratch.interpolation);
                    interpolation_f2c(num_coarsesample_theta,num_coarsesample_phi,num_finesample_theta,num_finesample_phi,num_interpol,indexloc_theta,indexloc_phi,inter_weights_theta,inter_weights_phi,temp_fine_phi,temp_coarse_phi,scratch.interpolation);
                    long int ind_child=agg_index[box_child];
                    ind_temp=0;
                    for(int indtheta=0;indtheta<num_coarsesample_theta;indtheta++){
                        for(int indphi=0;indphi<num_coarsesample_phi;indphi++){
                            disaggregation[ind_child*2+0]+=temp_coarse_theta[ind_temp];
                            disaggregation[ind_child*2+1]+=temp_coarse_phi[ind_temp];
                            ind_temp++;
                            ind_child++;
                        }
                    }
                }

                timer_total+=omp_get_wtime()-timer_in;
            }

            double timer_end = omp_get_wtime();
//            fprintf(file_timer,"MVMFarDis%d %d %1.5f %1.5f \n", indlev, omp_get_thread_num(), timer_end-timer_start, timer_total );
        }
    }



    //
    //
    // Disaggregation operation at the last level
    // Disaggregation is performed on every testing function at each box.
    // Receiving patterns of each testing function is shifted from box center to the sampling position.
    // Integration over the fields are finalized and assigned on the testing function index of y.
    box_start=num_filledboxes[num_mlfmalevels*2-1];
    box_end=num_filledboxes[num_mlfmalevels*3-1]+1;
    omp_set_num_threads(num_threads);

    #pragma omp parallel if(parallel_omp_on) shared(disaggregation,y)
    {
        double shiftvector[3];
        double khat[3]{0.0, 0.0, 0.0};
        double dotval;
        complex<double> expfactor;
        int numberofterms=num_truncation[num_mlfmalevels-1];
        int numberofthetas=numberofterms+1;
        int numberofphis=numberofthetas*2;
        const double* sintheta = leaf_grid.sin_theta.data();
        const double* costheta = leaf_grid.cos_theta.data();
        const double* sinphi = leaf_grid.sin_phi.data();
        const double* cosphi = leaf_grid.cos_phi.data();
        int testnode1,testnode2,testnode3;
        double testnode1loc[3],testnode2loc[3],testnode3loc[3],testlocation[3];
        double testlength,testsign,testvector[3],ncrosstestvector[3],testglobalnormal[3];
        double testreceivingtheta,testreceivingphi;
        double rwgsign[2]{1.0,-1.0};
        int rwglastnodeindex[2]{2,3};
        int testfunction, testtri;

        int tid=omp_get_thread_num();
        double timer_start = omp_get_wtime();
        double timer_total = 0.0;

        #pragma omp for schedule(dynamic, 5)
        for(int box_test=box_start;box_test<box_end;box_test++){
            double timer_in = omp_get_wtime();
            double testboxcenter[3]={box_pos[box_test*3+0],box_pos[box_test*3+1],box_pos[box_test*3+2]};
            int start_test=children[box_test*2+0];
            int end_test=children[box_test*2+1]+1;
            for(int ind_test=start_test;ind_test<end_test;ind_test++){
                complex<double> val_integral11=complex0;
                complex<double> val_integral21=complex0;
                for(int rwgtri=0; rwgtri<2; rwgtri++){
                    testfunction=rwg_ind[ind_test];
                    testnode1=rwg_nodes[testfunction*4+0];
                    testnode2=rwg_nodes[testfunction*4+1];
                    testnode3=rwg_nodes[testfunction*4+rwglastnodeindex[rwgtri]];
                    testnode1loc[0]=node_pos[testnode1*3+0];
                    testnode1loc[1]=node_pos[testnode1*3+1];
                    testnode1loc[2]=node_pos[testnode1*3+2];
                    testnode2loc[0]=node_pos[testnode2*3+0];
                    testnode2loc[1]=node_pos[testnode2*3+1];
                    testnode2loc[2]=node_pos[testnode2*3+2];
                    testnode3loc[0]=node_pos[testnode3*3+0];
                    testnode3loc[1]=node_pos[testnode3*3+1];
                    testnode3loc[2]=node_pos[testnode3*3+2];
                    testlocation[0]=(testnode1loc[0]+testnode2loc[0]+testnode3loc[0])/3;
                    testlocation[1]=(testnode1loc[1]+testnode2loc[1]+testnode3loc[1])/3;
                    testlocation[2]=(testnode1loc[2]+testnode2loc[2]+testnode3loc[2])/3;
                    testlength=rwg_lengths[testfunction];
                    testsign=rwgsign[rwgtri];
                    testtri=rwg_tris[testfunction*2+rwgtri];
                    testglobalnormal[0]=vec_trinormals[testtri*3+0];
                    testglobalnormal[1]=vec_trinormals[testtri*3+1];
                    testglobalnormal[2]=vec_trinormals[testtri*3+2];
                    shiftvector[0]=testlocation[0]-testboxcenter[0];
                    shiftvector[1]=testlocation[1]-testboxcenter[1];
                    shiftvector[2]=testlocation[2]-testboxcenter[2];
                    testvector[0]=testlocation[0]-testnode3loc[0];
                    testvector[1]=testlocation[1]-testnode3loc[1];
                    testvector[2]=testlocation[2]-testnode3loc[2];
                    crossProduct3d(testglobalnormal,testvector,ncrosstestvector);
                    long int inddisagg=agg_index[box_test]*2;
                    for(int indtheta=0;indtheta<numberofthetas;indtheta++){
                        for(int indphi=0;indphi<numberofphis;indphi++){
                            khat[0]=sintheta[indtheta]*cosphi[indphi];
                            khat[1]=sintheta[indtheta]*sinphi[indphi];
                            khat[2]=costheta[indtheta];
                            dotval=dotProduct3d(shiftvector,khat);
                            expfactor=exp(i*val_kO*dotval);
                            testreceivingtheta=costheta[indtheta]*cosphi[indphi]*testvector[0]+costheta[indtheta]*sinphi[indphi]*testvector[1]-sintheta[indtheta]*testvector[2];
                            testreceivingphi=-sinphi[indphi]*testvector[0]+cosphi[indphi]*testvector[1];
                            val_integral11+=disaggregation[inddisagg]*testsign*testlength*expfactor*testreceivingtheta/2.0;
                            val_integral11+=disaggregation[inddisagg+1]*testsign*testlength*expfactor*testreceivingphi/2.0;

                            val_integral21+=disaggregation[inddisagg]*testsign*testlength*expfactor*testreceivingphi/2.0;
                            val_integral21-=disaggregation[inddisagg+1]*testsign*testlength*expfactor*testreceivingtheta/2.0;
                            inddisagg++;
                            inddisagg++;
                        }
                    }
                }
                y[testfunction]+=val_integral11*y_constant;
                y[testfunction+num_rwgs]+=val_integral21*y_constant;
            }
            timer_total+=omp_get_wtime()-timer_in;
        }
        double timer_end = omp_get_wtime();
//        fprintf(file_timer,"MVMFarDis%d %d %1.5f %1.5f \n", num_mlfmalevels-1, omp_get_thread_num(), timer_end-timer_start, timer_total );
    }



    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
//    printf("Zfar MVM Calculation Takes %f s \n", time_span.count());

}

