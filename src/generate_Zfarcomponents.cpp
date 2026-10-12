#include "generate_Zfarcomponents.h"

void generate_Zfarcomponents(){
    //
    //
    // generate_Zfarcomponents
    // Determine the size of the aggregation and disaggregation array size.
    // The number of samples of each level is used for the calculation.
    int startbox=0;
    int endbox=0;
    int num_theta=0;
    int num_phi=0;

    size_agg=0;
    agg_index=new long int[num_filledboxesall];
    for(int ind_box=0; ind_box<num_filledboxesall; ind_box++)
        agg_index[ind_box]=0;

    for(int indlev=2; indlev<num_mlfmalevels; indlev++){
        startbox=num_filledboxes[num_mlfmalevels*1+indlev];
        endbox=num_filledboxes[num_mlfmalevels*2+indlev]+1;
        num_theta=num_truncation[indlev]+1;
        num_phi=num_theta*2;
        for(int ind_box=startbox; ind_box<endbox; ind_box++){
            agg_index[ind_box]=size_agg;
            size_agg+=long(num_theta)*long(num_phi);
        }
    }
    printf("Aggregation And Disaggregation Operator Size: %ld  [%f GB]\n", size_agg*4, size_agg*4.0*16.0/1024.0/1024/1024);

    if(val_pc==2){
        size_agg_pc=0;
        agg_index_pc=new long int[num_filledboxesall];
        for(int ind_box=0; ind_box<num_filledboxesall; ind_box++)
            agg_index_pc[ind_box]=0;

        for(int indlev=2; indlev<num_mlfmalevels; indlev++){
            startbox=num_filledboxes[num_mlfmalevels*1+indlev];
            endbox=num_filledboxes[num_mlfmalevels*2+indlev]+1;
            num_theta=num_truncation_pc[indlev]+1;
            num_phi=num_theta*2;
            for(int ind_box=startbox; ind_box<endbox; ind_box++){
                agg_index_pc[ind_box]=size_agg_pc;
                size_agg_pc+=long(num_theta)*long(num_phi);
            }
        }
        printf("PC Aggregation And Disaggregation Operator Size: %ld  [%f GB]\n", size_agg_pc*4, size_agg_pc*4.0*16.0/1024.0/1024/1024);

    }

    std::cout.flush();

}

