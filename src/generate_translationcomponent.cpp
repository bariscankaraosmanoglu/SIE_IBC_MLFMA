#include "generate_translationcomponent.h"
#include "input_validation.h"
#include <climits>
#include <vector>

void generate_translationcomponent(){
    try {
    //
    //
    // generate_translationcomponent
    // Calculate the size of the translation operator array and determine its index for each translation direction.
    // The farboxes are checked and the necessary translation directions are indexed.
    int startbox;
    int endbox;
    size_trans=0;
    trans_index=new int[num_mlfmalevels*7*7*7]; //num_mlfmalevels*343
    for(int ind1=0;ind1<(num_mlfmalevels*7*7*7);ind1++){
        trans_index[ind1]=-1;
    }
    for(int indlev=2; indlev<num_mlfmalevels; indlev++){
        startbox=num_filledboxes[num_mlfmalevels*1+indlev];
        endbox=num_filledboxes[num_mlfmalevels*2+indlev]+1;
        for(int testbox=startbox; testbox<endbox; testbox++){
            double testboxcenter[3]={box_pos[testbox*3+0], box_pos[testbox*3+1], box_pos[testbox*3+2]};
            for(int indfarbox=0; indfarbox<189; indfarbox++){
                if(boxes_far[testbox*189+indfarbox]!=-1){
                    int basisbox=boxes_far[testbox*189+indfarbox];
                    double basisboxcenter[3]={box_pos[basisbox*3+0], box_pos[basisbox*3+1], box_pos[basisbox*3+2]};
                    int delta[3];
                    delta[0]=round((testboxcenter[0]-basisboxcenter[0])/size_boxes[indlev])+3;
                    delta[1]=round((testboxcenter[1]-basisboxcenter[1])/size_boxes[indlev])+3;
                    delta[2]=round((testboxcenter[2]-basisboxcenter[2])/size_boxes[indlev])+3;
                    int indexval=indlev*343+delta[0]*49+delta[1]*7+delta[2];
                    if(trans_index[indexval]==-1){
                        trans_index[indexval]=size_trans;
                        size_trans=checked_index_sum(size_trans,2L*(num_truncation[indlev]+1)*(num_truncation[indlev]+1),INT_MAX,"translation");
                    }
                }
            }
        }
    }



    if(val_pc==2){
        size_trans_pc=0;
        trans_index_pc=new int[num_mlfmalevels*7*7*7]; //num_mlfmalevels*343
        for(int ind1=0;ind1<(num_mlfmalevels*7*7*7);ind1++){
            trans_index_pc[ind1]=-1;
        }
        for(int indlev=2; indlev<num_mlfmalevels; indlev++){
            startbox=num_filledboxes[num_mlfmalevels*1+indlev];
            endbox=num_filledboxes[num_mlfmalevels*2+indlev]+1;
            for(int testbox=startbox; testbox<endbox; testbox++){
                double testboxcenter[3]={box_pos[testbox*3+0], box_pos[testbox*3+1], box_pos[testbox*3+2]};
                for(int indfarbox=0; indfarbox<189; indfarbox++){
                    if(boxes_far[testbox*189+indfarbox]!=-1){
                        int basisbox=boxes_far[testbox*189+indfarbox];
                        double basisboxcenter[3]={box_pos[basisbox*3+0], box_pos[basisbox*3+1], box_pos[basisbox*3+2]};
                        int delta[3];
                        delta[0]=round((testboxcenter[0]-basisboxcenter[0])/size_boxes[indlev])+3;
                        delta[1]=round((testboxcenter[1]-basisboxcenter[1])/size_boxes[indlev])+3;
                        delta[2]=round((testboxcenter[2]-basisboxcenter[2])/size_boxes[indlev])+3;
                        int indexval=indlev*343+delta[0]*49+delta[1]*7+delta[2];
                        if(trans_index_pc[indexval]==-1){
                            trans_index_pc[indexval]=size_trans_pc;
                            size_trans_pc=checked_index_sum(size_trans_pc,2L*(num_truncation_pc[indlev]+1)*(num_truncation_pc[indlev]+1),INT_MAX,"translation");
                        }
                    }
                }
            }
        }
    }


    } catch (const InputError& error) { report_input_error(error); }
      catch (const std::bad_alloc&) { report_input_error(InputError("translation","unable to allocate index storage")); }
}






