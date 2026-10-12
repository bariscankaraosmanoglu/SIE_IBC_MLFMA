#include "generate_Znearcomponent.h"
#include "input_validation.h"
#include <climits>
#include <vector>

void generate_Znearcomponent(){
    try {
    //
    //
    // generate_Znearcomponent
    // Determine the number of near interactions, construct the reverse index the rwg functions for Znear calculation,
    // and construct triangle-based indexing for each least level box.
    int ind_rwg=0;
    int maxrwgsperbox=0;
    size_Znear=0;
    int startbox=num_filledboxes[num_mlfmalevels*2-1];
    int endbox=num_filledboxes[num_mlfmalevels*3-1]+1;
    int num_leastlevelbox=num_filledboxes[num_mlfmalevels-1];

    // Determine the number of interactions in Znear matrix
    nearindex=new long int[num_leastlevelbox*27];
    for(int box_test=startbox; box_test<endbox; box_test++){
        for(int ind_nearbox=0; ind_nearbox<27; ind_nearbox++){
            if(boxes_near[box_test*27+ind_nearbox]!=-1){
                int box_basis=boxes_near[box_test*27+ind_nearbox];
                nearindex[(box_test-startbox)*27+ind_nearbox]=size_Znear;
                const long block=checked_index_product(num_children[box_test],num_children[box_basis],INT_MAX,"near matrix");
                size_Znear=checked_index_sum(size_Znear,block,INT_MAX,"near matrix");
            }
        }
    }

    // Generate the inverse index of the position-based sorted rwg_tris
    ind_rwg=0;
    rwg_reverseind=new int[num_rwgs];
    rwg2box=new int[num_rwgs];
    for(int ind_box=startbox; ind_box<endbox; ind_box++){
        for(int ind_rwginbox=0; ind_rwginbox<num_children[ind_box]; ind_rwginbox++){
            int ind_rev=rwg_ind[ind_rwg];
            rwg_reverseind[ind_rev]=ind_rwg;
            rwg2box[ind_rev]=ind_box;
            ind_rwg++;
        }
    }

    maxrwgsperbox=0;
    for(int ind_box=startbox; ind_box<endbox; ind_box++){
        if(maxrwgsperbox<num_children[ind_box]){
            maxrwgsperbox=num_children[ind_box];
        }
    }
    ind_rwg=0;
    val_maxtri=0;
    int ind_newtriinbox=0;
    std::vector<int> boxtri_storage(size_t(maxrwgsperbox)*2+1);
    int* boxtri=boxtri_storage.data();
    // Calculating the maximum number of triangles in boxes for allocation of "smallestbox2tri"
    for(int ind_box=startbox; ind_box<endbox; ind_box++){
        for(int ind_tri=0; ind_tri<(maxrwgsperbox*2+1); ind_tri++){
            boxtri[ind_tri]=-1;
        }
        for(int indrwginbox=0; indrwginbox<num_children[ind_box]; indrwginbox++){
            int ind_rev=rwg_ind[ind_rwg];
            boxtri[2*indrwginbox+0]=rwg_tris[ind_rev*2+0];
            boxtri[2*indrwginbox+1]=rwg_tris[ind_rev*2+1];
            ind_rwg++;
        }
        qsort(boxtri, (maxrwgsperbox*2+1), sizeof(int), compareint);
        ind_newtriinbox=0;
        for(int ind_triinbox=0; ind_triinbox<(maxrwgsperbox*2+1); ind_triinbox++){
            if(boxtri[ind_triinbox]!=-1){
                if(boxtri[ind_triinbox]!=boxtri[ind_triinbox-1]){
                    ind_newtriinbox++;
                }
            }
        }
        if(val_maxtri<ind_newtriinbox){
            val_maxtri=ind_newtriinbox;
        }
    }
//    printf("max number of triangles in boxes: %d \n", val_maxtri);

    ind_rwg=0;
    ind_newtriinbox=0;
    checked_index_product(num_leastlevelbox,val_maxtri,INT_MAX,"triangle-to-box indices");
    smallestbox2tri=new int[num_leastlevelbox*val_maxtri];
    for(int indtri=0; indtri<(num_leastlevelbox*val_maxtri); indtri++)
        smallestbox2tri[indtri]=-1;
    // Generate the triangle to box for Znear calculation
    for(int ind_box=startbox; ind_box<endbox; ind_box++){
        for(int ind_tri=0; ind_tri<(maxrwgsperbox*2+1); ind_tri++){
            boxtri[ind_tri]=-1;
        }
        for(int indrwginbox=0; indrwginbox<num_children[ind_box]; indrwginbox++){
            int ind_rev=rwg_ind[ind_rwg];
            boxtri[2*indrwginbox+0]=rwg_tris[ind_rev*2+0];
            boxtri[2*indrwginbox+1]=rwg_tris[ind_rev*2+1];
            ind_rwg++;
        }
        qsort(boxtri, (maxrwgsperbox*2+1), sizeof(int), compareint);
        ind_newtriinbox=0;
        for(int ind_triinbox=0; ind_triinbox<(maxrwgsperbox*2+1); ind_triinbox++){
            if(boxtri[ind_triinbox]!=-1){
                if(boxtri[ind_triinbox]!=boxtri[ind_triinbox-1]){
                    smallestbox2tri[(ind_box-startbox)*val_maxtri+ind_newtriinbox]=boxtri[ind_triinbox];
                    ind_newtriinbox++;
                }
            }
        }
    }
    } catch (const InputError& error) { report_input_error(error); }
      catch (const std::bad_alloc&) { report_input_error(InputError("near matrix","unable to allocate index storage")); }
}

int compareint (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}
