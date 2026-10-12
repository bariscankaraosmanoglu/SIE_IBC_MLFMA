#include "generate_Zfarcomponents.h"
#include "input_validation.h"
#include <climits>
#include <vector>
#include <algorithm>

void generate_Zfarcomponents(){
    try {
    //
    //
    // generate_Zfarcomponents
    // Determine the size of the aggregation and disaggregation array size.
    // The number of samples of each level is used for the calculation.
    int startbox=0;
    int endbox=0;
    int num_theta=0;
    int num_phi=0;

    // IBC far interactions start at level 2. Keep leaf storage for a two-level
    // near-only tree, without allocating unused root/coarsest grids.
    const int first_active_level = std::min(2,num_mlfmalevels-1);
    size_agg=0;
    agg_index=new long int[num_filledboxesall];
    for(int ind_box=0; ind_box<num_filledboxesall; ind_box++)
        agg_index[ind_box]=0;

    for(int indlev=first_active_level; indlev<num_mlfmalevels; indlev++){
        startbox=num_filledboxes[num_mlfmalevels*1+indlev];
        endbox=num_filledboxes[num_mlfmalevels*2+indlev]+1;
        num_theta=num_truncation[indlev]+1;
        num_phi=num_theta*2;
        for(int ind_box=startbox; ind_box<endbox; ind_box++){
            agg_index[ind_box]=size_agg;
            size_agg=checked_index_sum(size_agg,long(num_theta)*num_phi,INT_MAX/2,"aggregation");
        }
    }
    printf("Aggregation And Disaggregation Operator Size: %ld  [%f GB]\n", size_agg*4, size_agg*4.0*16.0/1024.0/1024/1024);

    if(val_pc==2){
        size_agg_pc=0;
        agg_index_pc=new long int[num_filledboxesall];
        for(int ind_box=0; ind_box<num_filledboxesall; ind_box++)
            agg_index_pc[ind_box]=0;

        for(int indlev=first_active_level; indlev<num_mlfmalevels; indlev++){
            startbox=num_filledboxes[num_mlfmalevels*1+indlev];
            endbox=num_filledboxes[num_mlfmalevels*2+indlev]+1;
            num_theta=num_truncation_pc[indlev]+1;
            num_phi=num_theta*2;
            for(int ind_box=startbox; ind_box<endbox; ind_box++){
                agg_index_pc[ind_box]=size_agg_pc;
                size_agg_pc=checked_index_sum(size_agg_pc,long(num_theta)*num_phi,INT_MAX/2,"approximate aggregation");
            }
        }
        printf("PC Aggregation And Disaggregation Operator Size: %ld  [%f GB]\n", size_agg_pc*4, size_agg_pc*4.0*16.0/1024.0/1024/1024);

    }

    std::cout.flush();
    } catch (const InputError& error) { report_input_error(error); }
      catch (const std::bad_alloc&) { report_input_error(InputError("aggregation","unable to allocate index storage")); }
}

