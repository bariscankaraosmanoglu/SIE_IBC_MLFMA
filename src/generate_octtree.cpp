#include "mlfma_workspace.h"
#include "generate_octtree.h"
#include "input_validation.h"
#include <climits>

void generate_octtree(){
    // A new tree starts a new sampling/workspace lifetime.
    mlfma_workspace().clear();
    try {
    require_input(num_rwgs > 0, "RWG", "mesh has no full RWG unknowns; isolated triangles are unsupported by this discretization");
    //
    //
    // generate_octtree
    // Construct the oct-tree structure of MLFMA. Determine the clusters/boxes of RWG functions.
    double objectmin[3];
    double objectmax[3];
    double objectcenter[3];
    double objectsize;
    double boxmin[3];
    int rwgnode1;
    int rwgnode2;
    double rwgnode1pos[3];
    double rwgnode2pos[3];
    double rwgcenter[3];
    int rwgpositionTemp[3];
    int poslimit;
    int *rwgposition;
    long int *rwglocationsTemp;
    sortedElements *rwglocations;


    // Determine the min-max location, the center point, and the object size of the mesh-geometry.
    objectmin[0]=node_pos[0];
    objectmin[1]=node_pos[1];
    objectmin[2]=node_pos[2];
    objectmax[0]=node_pos[0];
    objectmax[1]=node_pos[1];
    objectmax[2]=node_pos[2];
    for(int ind1=1; ind1<num_nodes; ind1++){
        if(objectmin[0]>node_pos[ind1*3+0])
            objectmin[0]=node_pos[ind1*3+0];
        if(objectmin[1]>node_pos[ind1*3+1])
            objectmin[1]=node_pos[ind1*3+1];
        if(objectmin[2]>node_pos[ind1*3+2])
            objectmin[2]=node_pos[ind1*3+2];
        if(objectmax[0]<node_pos[ind1*3+0])
            objectmax[0]=node_pos[ind1*3+0];
        if(objectmax[1]<node_pos[ind1*3+1])
            objectmax[1]=node_pos[ind1*3+1];
        if(objectmax[2]<node_pos[ind1*3+2])
            objectmax[2]=node_pos[ind1*3+2];
    }
    objectcenter[0]=(objectmax[0]+objectmin[0])/2;
    objectcenter[1]=(objectmax[1]+objectmin[1])/2;
    objectcenter[2]=(objectmax[2]+objectmin[2])/2;
    for (int axis=0; axis<3; ++axis)
        require_input(std::isfinite(objectcenter[axis]) && std::isfinite(objectmax[axis]-objectmin[axis]),
                      "POINTS", "bounding box is not representable");
    objectsize=0.0;
    for(int ind1=0; ind1<3; ind1++){
        if(objectsize<(objectmax[ind1]-objectmin[ind1]))
            objectsize=objectmax[ind1]-objectmin[ind1];
    }
    printf("Min: %f %f %f \n", objectmin[0],objectmin[1],objectmin[2]);
    printf("Max: %f %f %f \n", objectmax[0],objectmax[1],objectmax[2]);
    printf("Center: %f %f %f \n", objectcenter[0],objectcenter[1],objectcenter[2]);
    printf("Model Size: %f \n", objectsize);
    printf("\n");

    // Resolve and validate every dimension before allocating tree/operator arrays.
    extern int num_interpol;
    extern double dir_ff;
    const MlfmaDimensions dimensions = validate_mlfma_dimensions(num_mlfmalevels,
        objectsize,val_lambdaO,val_kO,num_digits,val_pc,val_approx,num_interpol,dir_ff);
    num_mlfmalevels=dimensions.levels;
    printf("MLFMA Levels: %d \n", num_mlfmalevels);
    size_boxes=new double[num_mlfmalevels];
    num_truncation=new int[num_mlfmalevels];
    for(int level=0;level<num_mlfmalevels;++level) {
        size_boxes[level]=dimensions.boxes[level];
        num_truncation[level]=dimensions.terms[level];
        printf("Level: %d, Num. Harmonics: %d \n",level,num_truncation[level]);
    }
    const double smallestboxsize=size_boxes[num_mlfmalevels-1];
    for(int axis=0;axis<3;++axis) boxmin[axis]=objectcenter[axis]-objectsize/2;
    if(val_pc==2) {
        num_truncation_pc=new int[num_mlfmalevels];
        for(int level=0;level<num_mlfmalevels;++level) {
            num_truncation_pc[level]=dimensions.approximate_terms[level];
            printf("Level: %d, Num. Harmonics PC: %d \n",level,num_truncation_pc[level]);
        }
    }

    // Calculate RWG location indices at the least level boxes. Storing the variables in XYZ coordinates and in combined-form.
    // The RWG box location indices are sorted for the calculation of box locations at the higher levels.
    rwgposition=new int[num_rwgs*3];
    rwglocations=new sortedElements[num_rwgs];


    for(int ind1=0; ind1<num_rwgs; ind1++){
        rwgposition[ind1*3+0]=0;
        rwgposition[ind1*3+1]=0;
        rwgposition[ind1*3+2]=0;
        rwglocations[ind1].elemID=ind1;
        rwglocations[ind1].elemVal=0;
    }
    poslimit=pow(2,num_mlfmalevels-1)-1;
    for(int ind1=0; ind1<num_rwgs; ind1++){
        rwgnode1=rwg_nodes[ind1*4+0];
        rwgnode2=rwg_nodes[ind1*4+1];
        rwgnode1pos[0]=node_pos[rwgnode1*3+0];
        rwgnode1pos[1]=node_pos[rwgnode1*3+1];
        rwgnode1pos[2]=node_pos[rwgnode1*3+2];
        rwgnode2pos[0]=node_pos[rwgnode2*3+0];
        rwgnode2pos[1]=node_pos[rwgnode2*3+1];
        rwgnode2pos[2]=node_pos[rwgnode2*3+2];
        rwgcenter[0]=(rwgnode1pos[0]+rwgnode2pos[0])/2;
        rwgcenter[1]=(rwgnode1pos[1]+rwgnode2pos[1])/2;
        rwgcenter[2]=(rwgnode1pos[2]+rwgnode2pos[2])/2;
        rwgpositionTemp[0]=floor((rwgcenter[0]-boxmin[0])/smallestboxsize);
        rwgpositionTemp[1]=floor((rwgcenter[1]-boxmin[1])/smallestboxsize);
        rwgpositionTemp[2]=floor((rwgcenter[2]-boxmin[2])/smallestboxsize);
        if(rwgpositionTemp[0]>poslimit)
            rwgpositionTemp[0]=poslimit;
        if(rwgpositionTemp[0]<0)
            rwgpositionTemp[0]=0;
        if(rwgpositionTemp[1]>poslimit)
            rwgpositionTemp[1]=poslimit;
        if(rwgpositionTemp[1]<0)
            rwgpositionTemp[1]=0;
        if(rwgpositionTemp[2]>poslimit)
            rwgpositionTemp[2]=poslimit;
        if(rwgpositionTemp[2]<0)
            rwgpositionTemp[2]=0;
        rwgposition[ind1*3+0]=rwgpositionTemp[0];
        rwgposition[ind1*3+1]=rwgpositionTemp[1];
        rwgposition[ind1*3+2]=rwgpositionTemp[2];
        rwglocations[ind1].elemVal=dec2bin2dec(rwgpositionTemp[0],rwgpositionTemp[1],rwgpositionTemp[2],num_mlfmalevels);
    }
    qsort(rwglocations, num_rwgs, sizeof(struct sortedElements), compare);


    // Resorting according to the index values
    int ind_rwg=0;
    while(ind_rwg<num_rwgs){
        int ind_start=ind_rwg;
        int num_size1=0;
        while(ind_rwg<num_rwgs&&rwglocations[ind_rwg].elemVal==rwglocations[ind_start].elemVal){
            ind_rwg++;
            num_size1++;
        }
        int ind_end=ind_rwg;
        if(num_size1>0){
            struct sortedElements rwglocationsbuffer[num_size1];
            int ind_b=0;
            for(int ind_s=ind_start;ind_s<ind_end;ind_s++){
                rwglocationsbuffer[ind_b].elemVal=rwglocations[ind_s].elemVal;
                rwglocationsbuffer[ind_b].elemID=rwglocations[ind_s].elemID;
                ind_b++;
            }
            qsort(rwglocationsbuffer, num_size1, sizeof(struct sortedElements), compareind);
            ind_b=0;
            for(int ind_s=ind_start;ind_s<ind_end;ind_s++){
                rwglocations[ind_s].elemVal=rwglocationsbuffer[ind_b].elemVal;
                rwglocations[ind_s].elemID=rwglocationsbuffer[ind_b].elemID;
                ind_b++;
            }
        }
        else{
            num_size1++;
        }
    }


    // Determine the total number of filled boxes "ind_filledboxes"
    // for the allocation of arrays "num_children" and "box_pos"
    rwglocationsTemp=new long int[num_rwgs];
    rwg_ind=new int[num_rwgs];
    for(int ind1=0;ind1<num_rwgs;ind1++){
        rwglocationsTemp[ind1]=rwglocations[ind1].elemVal;
        rwg_ind[ind1]=rwglocations[ind1].elemID;

    }
    num_filledboxes=new int[num_mlfmalevels*3];
	num_filledboxesall=0;
    for(int indlev=num_mlfmalevels; indlev>0; indlev--){
        num_filledboxes[indlev]=1;
        for(int indrwg=1; indrwg<num_rwgs; indrwg++){
            if (rwglocationsTemp[indrwg]!=rwglocationsTemp[indrwg-1]){
                num_filledboxes[indlev]++;
            }
        }
        num_filledboxesall=checked_index_sum(num_filledboxesall,num_filledboxes[indlev],INT_MAX/189,"num_mlfmalevels box indices");
        for(int indrwg=0; indrwg<num_rwgs; indrwg++){
            rwglocationsTemp[indrwg]/=8;
        }
    }

    // Calculate the box positions in cartesian coordinates "box_pos" and determine the number of children at each box "num_children"
    for(int ind1=0;ind1<num_rwgs;ind1++){
        rwglocationsTemp[ind1]=rwglocations[ind1].elemVal;
    }
    num_children=new int[num_filledboxesall];
    box_pos=new double[num_filledboxesall*3];
    double boxsizesvec=0.0;
    int ind_filledboxes=-1;
    for(int indlev=num_mlfmalevels-1; indlev>=0; indlev--){
        boxsizesvec=size_boxes[indlev]/2;
        num_filledboxes[indlev]=1;
        num_children[ind_filledboxes+num_filledboxes[indlev]]=1;
        box_pos[(ind_filledboxes+num_filledboxes[indlev])*3+0] = boxmin[0]+rwgposition[rwglocations[0].elemID*3+0]*size_boxes[indlev]+boxsizesvec;
        box_pos[(ind_filledboxes+num_filledboxes[indlev])*3+1] = boxmin[1]+rwgposition[rwglocations[0].elemID*3+1]*size_boxes[indlev]+boxsizesvec;
        box_pos[(ind_filledboxes+num_filledboxes[indlev])*3+2] = boxmin[2]+rwgposition[rwglocations[0].elemID*3+2]*size_boxes[indlev]+boxsizesvec;
        for(int indrwg=1; indrwg<num_rwgs; indrwg++){
            if (rwglocationsTemp[indrwg]==rwglocationsTemp[indrwg-1]){
                if((rwglocations[indrwg].elemVal!=rwglocations[indrwg-1].elemVal) || (indlev==num_mlfmalevels-1))
                    num_children[ind_filledboxes+num_filledboxes[indlev]]++;
            } else{
                num_filledboxes[indlev]++;
                num_children[ind_filledboxes+num_filledboxes[indlev]]=1;
                box_pos[(ind_filledboxes+num_filledboxes[indlev])*3+0] = boxmin[0]+rwgposition[rwglocations[indrwg].elemID*3+0]*size_boxes[indlev]+boxsizesvec;
                box_pos[(ind_filledboxes+num_filledboxes[indlev])*3+1] = boxmin[1]+rwgposition[rwglocations[indrwg].elemID*3+1]*size_boxes[indlev]+boxsizesvec;
                box_pos[(ind_filledboxes+num_filledboxes[indlev])*3+2] = boxmin[2]+rwgposition[rwglocations[indrwg].elemID*3+2]*size_boxes[indlev]+boxsizesvec;
            }
        }
        ind_filledboxes+=num_filledboxes[indlev];
        for(int indrwg=0; indrwg<num_rwgs; indrwg++)
            rwglocationsTemp[indrwg]/=8;
        for(int indrwg=0; indrwg<num_rwgs; indrwg++){
            rwgposition[indrwg*3+0]/=2;
            rwgposition[indrwg*3+1]/=2;
            rwgposition[indrwg*3+2]/=2;
        }
        if(indlev!=num_mlfmalevels-1){
            for(int indrwg=0; indrwg<num_rwgs; indrwg++)
                rwglocations[indrwg].elemVal/=8;
        }
    }
    delete[] rwglocations;
    delete[] rwgposition;
    delete[] rwglocationsTemp;


    // Prepare the start and end index of the boxes at each level
    num_filledboxes[num_mlfmalevels*1+0]=0;
    num_filledboxes[num_mlfmalevels*2+0]=0;
    for(int indlev=1; indlev<num_mlfmalevels; indlev++){
        num_filledboxes[num_mlfmalevels*2+indlev]=num_filledboxes[num_mlfmalevels*0+indlev]+num_filledboxes[num_mlfmalevels*2+indlev-1];
        num_filledboxes[num_mlfmalevels*1+indlev]=num_filledboxes[num_mlfmalevels*2+indlev]-num_filledboxes[num_mlfmalevels*0+indlev]+1;
    }

    // Flip the arrays "num_children" and "box_pos" for MLFMA operations
    int *num_childrenTemp=new int[num_filledboxesall];
    for(int ind1=0; ind1<num_filledboxesall; ind1++){
        num_childrenTemp[ind1]=num_children[num_filledboxesall-ind1-1];
    }
    for(int indlev=0; indlev<num_mlfmalevels; indlev++){
        int ind1_start=num_filledboxes[num_mlfmalevels*1+indlev];
        int ind1_end=num_filledboxes[num_mlfmalevels*2+indlev];
        for(int ind1=0; ind1<num_filledboxes[num_mlfmalevels*0+indlev]; ind1++){
            num_children[ind1_start]=num_childrenTemp[ind1_end];
            ind1_start++;
            ind1_end--;
        }
    }
    delete[] num_childrenTemp;

    double *box_posTemp=new double[num_filledboxesall*3];
    for(int ind1=0; ind1<num_filledboxesall; ind1++){
        box_posTemp[ind1*3+0]=box_pos[(num_filledboxesall-ind1-1)*3+0];
        box_posTemp[ind1*3+1]=box_pos[(num_filledboxesall-ind1-1)*3+1];
        box_posTemp[ind1*3+2]=box_pos[(num_filledboxesall-ind1-1)*3+2];
    }
    for(int indlev=0; indlev<num_mlfmalevels; indlev++){
        int ind1_start=num_filledboxes[num_mlfmalevels*1+indlev];
        int ind1_end=num_filledboxes[num_mlfmalevels*2+indlev];
        for(int ind1=0; ind1<num_filledboxes[num_mlfmalevels*0+indlev]; ind1++){
            box_pos[ind1_start*3+0]=box_posTemp[ind1_end*3+0];
            box_pos[ind1_start*3+1]=box_posTemp[ind1_end*3+1];
            box_pos[ind1_start*3+2]=box_posTemp[ind1_end*3+2];
            ind1_start++;
            ind1_end--;
        }
    }
    delete[] box_posTemp;

    // Relate the children of the parent boxes "children" and parent of children "children2parent"
    children=new int[num_filledboxesall*2];
    children2parent=new int[num_filledboxesall];
    for(int indbox=0; indbox<num_filledboxesall; indbox++){
        children[indbox*2+0]=0;
        children[indbox*2+1]=0;
        children2parent[indbox]=-1;
    }
    for(int indlev=0; indlev<num_mlfmalevels-1; indlev++){
        int firstchild=num_filledboxes[num_mlfmalevels*1+indlev+1];
        int startbox=num_filledboxes[num_mlfmalevels*1+indlev];
        int endbox=num_filledboxes[num_mlfmalevels*2+indlev]+1;
        for(int indbox=startbox; indbox<endbox; indbox++){
            children[indbox*2+0]=firstchild;
            children[indbox*2+1]=firstchild+num_children[indbox]-1;
            int startparent=firstchild;
            int endparent=firstchild+num_children[indbox];
            for(int indparent=startparent; indparent<endparent; indparent++){
                children2parent[indparent]=indbox;
            }
            firstchild+=num_children[indbox];
        }
    }
    int firstchild=0;
    int startbox=num_filledboxes[num_mlfmalevels*2-1];
    int endbox=num_filledboxes[num_mlfmalevels*3-1]+1;
    for(int indbox=startbox; indbox<endbox; indbox++){
        children[indbox*2+0]=firstchild;
        children[indbox*2+1]=firstchild+num_children[indbox]-1;
        firstchild+=num_children[indbox];
    }


    printf("\n");
    for(int indlev=0;indlev<num_mlfmalevels;indlev++){
        int num_box=num_filledboxes[indlev];
        printf("Level: %d, Number of Boxes: %d \n", indlev, num_box);
    }
    printf("\n");

    std::cout.flush();
    } catch (const InputError& error) { report_input_error(error); }
      catch (const std::bad_alloc&) { report_input_error(InputError("num_mlfmalevels", "unable to allocate tree storage")); }

}

long int dec2bin2dec(int posx, int posy, int posz, int num_lev){
    long int power2=1;
    long int rwg_pos=0;
    for(int ind1=0;ind1<num_lev;ind1++){
        rwg_pos+=long(posz%2)*power2;
        posz/=2;
        power2*=2;
        rwg_pos+=long(posy%2)*power2;
        posy/=2;
        power2*=2;
        rwg_pos+=long(posx%2)*power2;
        posx/=2;
        power2*=2;

    }
    return rwg_pos;
}

int compare (const void *a, const void *b)
{
//    return ( (*(struct sortedElements*)a).elemVal - (*(struct sortedElements*)b).elemVal );
    if(( (*(struct sortedElements*)a).elemVal - (*(struct sortedElements*)b).elemVal )<0)
        return -1;
    if(( (*(struct sortedElements*)a).elemVal - (*(struct sortedElements*)b).elemVal )>0)
        return 1;
    return 0;

}

int compareind (const void *a, const void *b)
{
    return ( (*(struct sortedElements*)a).elemID - (*(struct sortedElements*)b).elemID );
}

