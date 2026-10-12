#include <stdio.h>

using namespace std;

extern int parallel_omp_on;
extern int num_threads;
extern int *num_filledboxes;
extern int num_mlfmalevels;
extern int *thread_boxes;
extern int *thread_boxes_Rad;
extern long int size_Znear;
extern long int *nearindex;
extern int num_rwgs;
extern int *children;

void calculate_loadbalance(){

    // Load balance operation is performed by calculating average number of interactions per thread
    // Each load is parsed in test boxes in Znear MVM

    long int thread_load=size_Znear/long(num_threads);
    int thread_id=1;

    int box_start=num_filledboxes[num_mlfmalevels*2-1];
    int box_end=num_filledboxes[num_mlfmalevels*3-1]+1;

    thread_boxes=new int[num_threads*2];
    thread_boxes[0]=box_start;
    thread_boxes[num_threads*2-1]=box_end;

    for(int box_test=box_start;box_test<box_end;box_test++){
        long int ind_Znear=nearindex[(box_test-box_start)*27];
        if (ind_Znear>=(thread_id*thread_load)){
            thread_boxes[thread_id*2+0]=box_test;
            thread_boxes[(thread_id-1)*2+1]=box_test;
            thread_id++;
        }
    }



    int thread_load_rad=num_rwgs/num_threads;
    thread_id=1;

    thread_boxes_Rad=new int[num_threads*2];
    thread_boxes_Rad[0]=box_start;
    thread_boxes_Rad[num_threads*2-1]=box_end;
    int ind_box_rwgs=0;
    for(int box_basis=box_start;box_basis<box_end;box_basis++){
        int start_basis=children[box_basis*2+0];
        int end_basis=children[box_basis*2+1]+1;
        ind_box_rwgs+=(end_basis-start_basis);
        if (ind_box_rwgs>=(thread_id*thread_load_rad)){
            thread_boxes_Rad[thread_id*2+0]=box_basis;
            thread_boxes_Rad[(thread_id-1)*2+1]=box_basis;
            thread_id++;
        }

    }



//    int box_start=num_filledboxes[num_mlfmalevels*2-1];
//    int box_end=num_filledboxes[num_mlfmalevels*3-1]+1;
//    long int thread_load=(box_end-box_start+1)/long(num_threads);
//    int thread_id=1;
//
//    thread_boxes=new int[num_threads*2];
//    thread_boxes[0]=box_start;
//    thread_boxes[num_threads*2-1]=box_end;
//
//    for(int box_test=box_start;box_test<box_end;box_test++){
//        if ((box_test-box_start+1)>=(thread_load*thread_id)){
//            thread_boxes[thread_id*2+0]=box_test;
//            thread_boxes[(thread_id-1)*2+1]=box_test;
//            thread_id++;
//        }
//    }




//    for(int tid=0; tid<num_threads; tid++){
//        printf("Thread ID: %d, Start Box: %d, End Box: %d, Znear Index: %ld \n", tid, thread_boxes[tid*2+0], thread_boxes[tid*2+1], nearindex[(thread_boxes[tid*2+0]-box_start)*27]);
//        printf("Thread ID: %d, Start Box: %d, End Box: %d, Number of Box: %d \n \n", tid, thread_boxes_Rad[tid*2+0], thread_boxes_Rad[tid*2+1], thread_boxes_Rad[tid*2+1]-thread_boxes_Rad[tid*2+0] );
//    }
//    printf("\n");

}
