#include "generate_nearfarboxes.h"

void generate_nearfarboxes(){
    boxes_near=new int[num_filledboxesall*27];
    boxes_far=new int[num_filledboxesall*189];
    for(int indbox=0; indbox<num_filledboxesall*27; indbox++){
        boxes_near[indbox]=-1;
    }
    for(int indbox=0; indbox<num_filledboxesall*189; indbox++){
        boxes_far[indbox]=-1;
    }
    boxes_near[0]=0;
    for(int indlev=1; indlev<num_mlfmalevels; indlev++){
        int startbox=num_filledboxes[num_mlfmalevels*1+indlev];
        int endbox=num_filledboxes[num_mlfmalevels*2+indlev]+1;
        for(int indbox=startbox; indbox<endbox; indbox++){
            int nearfieldcounter = 0;
            int farfieldcounter = 0;
            int parentbox = children2parent[indbox];
            double boxcenter[3]={box_pos[indbox*3+0], box_pos[indbox*3+1], box_pos[indbox*3+2]};
            for(int indparentboxnear=0; indparentboxnear<27; indparentboxnear++){
                if (boxes_near[parentbox*27+indparentboxnear]!= -1){
                    int unclebox=boxes_near[parentbox*27+indparentboxnear];
                    for(int indcousinbox=children[unclebox*2+0]; indcousinbox<(children[unclebox*2+1]+1); indcousinbox++){
                        double cousincenter[3]={box_pos[indcousinbox*3+0], box_pos[indcousinbox*3+1], box_pos[indcousinbox*3+2]};
                        double distance=sqrt((cousincenter[0]-boxcenter[0])*(cousincenter[0]-boxcenter[0])\
                                                  +(cousincenter[1]-boxcenter[1])*(cousincenter[1]-boxcenter[1])\
                                                  +(cousincenter[2]-boxcenter[2])*(cousincenter[2]-boxcenter[2]));
                        if(distance<1.89*size_boxes[indlev]){
                            boxes_near[indbox*27+nearfieldcounter]=indcousinbox;
                            nearfieldcounter++;
                        }
                        else{
                            boxes_far[indbox*189+farfieldcounter]=indcousinbox;
                            farfieldcounter++;
                        }
                    }
                }
                else{
                    break;
                }
            }
        }
    }
}
