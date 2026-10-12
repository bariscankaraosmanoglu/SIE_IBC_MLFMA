#include <stdio.h>
using namespace std;

extern double pi;

void interpolationweights(int num_coarsesample,int num_finesample,int num_interpolpoint,int flagphitheta,double *loc_coarse,double *loc_fine,double *inter_weights,int *indexloc){
    //
    //
    // interpolationweights
    // Lagrange interpolation weights and index locations

    double *loc_coarse_ex;
    loc_coarse_ex=new double[3*num_coarsesample];

    if (flagphitheta==1){ // Phi direction
        for(int ind1=0;ind1<num_coarsesample;ind1++){
            loc_coarse_ex[ind1+num_coarsesample*0]=loc_coarse[ind1]-2*pi;
            loc_coarse_ex[ind1+num_coarsesample*1]=loc_coarse[ind1];
            loc_coarse_ex[ind1+num_coarsesample*2]=loc_coarse[ind1]+2*pi;
        }
    }
    else{ // Theta direction
        for(int ind1=0;ind1<num_coarsesample;ind1++){
            loc_coarse_ex[ind1+num_coarsesample*0]=-loc_coarse[num_coarsesample-1-ind1];
            loc_coarse_ex[ind1+num_coarsesample*1]=loc_coarse[ind1];
            loc_coarse_ex[ind1+num_coarsesample*2]=2*pi-loc_coarse[num_coarsesample-1-ind1];
        }
    }

    for(int ind_fine=0;ind_fine<num_finesample;ind_fine++){
        int ind_coarse=0;
        while((loc_coarse_ex[ind_coarse+num_coarsesample]<loc_fine[ind_fine])&&(ind_coarse<=num_coarsesample)){
            ind_coarse++;
        }
        int ind_coarse_start=ind_coarse+num_coarsesample-num_interpolpoint;
        int ind_coarse_end=ind_coarse+num_coarsesample+num_interpolpoint;
        int ind_sample=0;
        for(int ind_coarse_ex=ind_coarse_start;ind_coarse_ex<ind_coarse_end;ind_coarse_ex++){
            double weight_val=1.0;
            for(int ind_coarse_ex_s=ind_coarse_start;ind_coarse_ex_s<ind_coarse_end;ind_coarse_ex_s++){
                if(ind_coarse_ex_s!=ind_coarse_ex){
                    weight_val=weight_val*(loc_fine[ind_fine]-loc_coarse_ex[ind_coarse_ex_s])/(loc_coarse_ex[ind_coarse_ex]-loc_coarse_ex[ind_coarse_ex_s]);
                }
            }
            inter_weights[ind_fine*num_interpolpoint*2+ind_sample]=weight_val;
            indexloc[ind_fine*num_interpolpoint*2+ind_sample]=ind_coarse_ex;
            ind_sample++;
        }
    }

    delete[] loc_coarse_ex;



}



