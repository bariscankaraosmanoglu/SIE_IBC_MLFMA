#include "mlfma_workspace.h"
#include <stdio.h>
#include <complex>
using namespace std;

void interpolation_c2f(int num_coarsesample_theta,int num_coarsesample_phi,int num_finesample_theta,int num_finesample_phi,int num_interpolpoint,int *indexloc_theta,int *indexloc_phi,double *inter_weights_theta,double *inter_weights_phi,complex<double> *field_coarse,complex<double> *field_fine,InterpolationScratch& scratch){
    //
    //
    // interpolation_c2f
    // Interpolation for MLFMA aggregation part.
    // This function interpolates the coarse input into the fine output
    // Two-step interpolation is applied
    complex<double> *sample_inflated_theta;
    complex<double> *field_fine_half;
    complex<double> *sample_inflated_phi;
    scratch.prepare(num_coarsesample_theta,num_coarsesample_phi,num_finesample_theta);
    sample_inflated_theta = scratch.inflated_theta.data();
    field_fine_half = scratch.half.data();
    sample_inflated_phi = scratch.inflated_phi.data();


    // Preprocess part for interpolation on theta direction
    // Inflation on theta direction is performed
    for(int ind0=0;ind0<num_coarsesample_theta;ind0++){
        for(int ind1=0;ind1<(num_coarsesample_phi/2);ind1++){
            int index0=(num_coarsesample_theta-ind0-1)*num_coarsesample_phi+(ind1+num_coarsesample_phi/2);
            int index1=ind0*num_coarsesample_phi+ind1;
            sample_inflated_theta[index1]=-field_coarse[index0];
        }

        for(int ind1=(num_coarsesample_phi/2);ind1<num_coarsesample_phi;ind1++){
            int index0=(num_coarsesample_theta-ind0-1)*num_coarsesample_phi+(ind1-num_coarsesample_phi/2);
            int index1=ind0*num_coarsesample_phi+ind1;
            sample_inflated_theta[index1]=-field_coarse[index0];
        }
    }
    for(int ind0=0;ind0<num_coarsesample_theta;ind0++){
        for(int ind1=0;ind1<num_coarsesample_phi;ind1++){
            int index0=ind0*num_coarsesample_phi+ind1;
            int index1=ind0*num_coarsesample_phi+num_coarsesample_theta*num_coarsesample_phi+ind1;
            sample_inflated_theta[index1]=field_coarse[index0];
        }
    }
    for(int ind0=0;ind0<num_coarsesample_theta;ind0++){
        for(int ind1=0;ind1<(num_coarsesample_phi/2);ind1++){
            int index0=(num_coarsesample_theta-ind0-1)*num_coarsesample_phi+(ind1+num_coarsesample_phi/2);
            int index1=ind0*num_coarsesample_phi+num_coarsesample_theta*num_coarsesample_phi*2+ind1;
            sample_inflated_theta[index1]=-field_coarse[index0];
        }
        for(int ind1=(num_coarsesample_phi/2);ind1<num_coarsesample_phi;ind1++){
            int index0=(num_coarsesample_theta-ind0-1)*num_coarsesample_phi+(ind1-num_coarsesample_phi/2);
            int index1=ind0*num_coarsesample_phi+num_coarsesample_theta*num_coarsesample_phi*2+ind1;
            sample_inflated_theta[index1]=-field_coarse[index0];
        }
    }

    // Interpolation on theta direction
    for(int ind_coarse_phi=0;ind_coarse_phi<num_coarsesample_phi;ind_coarse_phi++){
        for(int ind_fine_theta=0;ind_fine_theta<num_finesample_theta;ind_fine_theta++){
            complex<double> sample_sum=complex<double>(0.0, 0.0);
            for(int ind_sample=0;ind_sample<(2*num_interpolpoint);ind_sample++){
                int index0=ind_fine_theta*num_interpolpoint*2+ind_sample;
                int index1=indexloc_theta[index0]*num_coarsesample_phi+ind_coarse_phi;
                complex<double> sample_val=sample_inflated_theta[index1];
                double sample_wei=inter_weights_theta[index0];
                sample_sum+=sample_wei*sample_val;
            }
            field_fine_half[ind_fine_theta*num_coarsesample_phi+ind_coarse_phi]=sample_sum;
        }
    }

    // Preprocess part for interpolation on phi direction
    // Inflation on phi direction is performed
    for(int ind0=0;ind0<num_finesample_theta;ind0++){
        for(int ind1=0;ind1<num_coarsesample_phi;ind1++){
            int index0=ind0*num_coarsesample_phi+ind1;
            int index1=ind0*num_coarsesample_phi*3+ind1;
            sample_inflated_phi[index1]=field_fine_half[index0];
            sample_inflated_phi[index1+num_coarsesample_phi]=field_fine_half[index0];
            sample_inflated_phi[index1+num_coarsesample_phi*2]=field_fine_half[index0];
        }
    }

    // Interpolation on phi direction
    for(int ind_fine_theta=0;ind_fine_theta<num_finesample_theta;ind_fine_theta++){
        for(int ind_fine_phi=0;ind_fine_phi<num_finesample_phi;ind_fine_phi++){
            complex<double> sample_sum=complex<double>(0.0, 0.0);
            for(int ind_sample=0;ind_sample<(2*num_interpolpoint);ind_sample++){
                int index0=ind_fine_phi*num_interpolpoint*2+ind_sample;
                int index1=ind_fine_theta*num_coarsesample_phi*3+indexloc_phi[index0];
                complex<double> sample_val=sample_inflated_phi[index1];
                double sample_wei=inter_weights_phi[index0];
                sample_sum+=sample_wei*sample_val;
            }
            field_fine[ind_fine_theta*num_finesample_phi+ind_fine_phi]=sample_sum;

        }
    }



}







void interpolation_c2f(int num_coarsesample_theta,int num_coarsesample_phi,int num_finesample_theta,int num_finesample_phi,int num_interpolpoint,int *indexloc_theta,int *indexloc_phi,double *inter_weights_theta,double *inter_weights_phi,complex<double> *field_coarse,complex<double> *field_fine){
    InterpolationScratch scratch;
    interpolation_c2f(num_coarsesample_theta,num_coarsesample_phi,num_finesample_theta,num_finesample_phi,num_interpolpoint,indexloc_theta,indexloc_phi,inter_weights_theta,inter_weights_phi,field_coarse,field_fine,scratch);
}
