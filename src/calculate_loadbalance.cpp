#include "calculate_loadbalance.h"
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

extern int num_threads;
extern int *num_filledboxes;
extern int num_mlfmalevels;
extern int *num_children;
extern int *boxes_near;
extern long int size_Znear;

namespace {
// Geometry-dependent weights are immutable throughout assembly and solution.
std::vector<std::uint64_t> near_work_prefix;
int first_leaf_box;

int partition_boundary(int worker, int team_size){
    const int box_count=static_cast<int>(near_work_prefix.size())-1;
    if(worker==0) return first_leaf_box;
    if(worker==team_size) return first_leaf_box+box_count;

    const std::uint64_t total=near_work_prefix.back();
    if(total==0){
        return first_leaf_box+static_cast<int>(std::uint64_t(box_count)*worker/team_size);
    }
    // floor(total*worker/team_size), avoiding overflow of total*worker.
    const std::uint64_t target=(total/team_size)*worker
        +((total%team_size)*worker)/team_size;
    const auto boundary=std::lower_bound(near_work_prefix.begin(),near_work_prefix.end(),target);
    return first_leaf_box+static_cast<int>(boundary-near_work_prefix.begin());
}
}

void calculate_loadbalance(){
    if(num_threads<=0){
        std::fprintf(stderr,"Invalid thread count: %d; expected a positive integer.\n",num_threads);
        std::exit(EXIT_FAILURE);
    }
    first_leaf_box=num_filledboxes[num_mlfmalevels*2-1];
    const int end_box=num_filledboxes[num_mlfmalevels*3-1]+1;
    near_work_prefix.assign(end_box-first_leaf_box+1,0);
    for(int box=first_leaf_box;box<end_box;++box){
        std::uint64_t work=0;
        for(int neighbor=0;neighbor<27;++neighbor){
            const int basis_box=boxes_near[box*27+neighbor];
            if(basis_box!=-1){
                work+=std::uint64_t(num_children[box])*std::uint64_t(num_children[basis_box]);
            }
        }
        const int index=box-first_leaf_box;
        near_work_prefix[index+1]=near_work_prefix[index]+work;
    }
    assert(near_work_prefix.back()==std::uint64_t(size_Znear));
}

std::pair<int,int> near_box_range(int thread_id,int team_size){
    assert(!near_work_prefix.empty());
    assert(team_size>0 && thread_id>=0 && thread_id<team_size);
    const int begin=partition_boundary(thread_id,team_size);
    const int end=partition_boundary(thread_id+1,team_size);
    assert(first_leaf_box<=begin && begin<=end);
    assert(end<=first_leaf_box+static_cast<int>(near_work_prefix.size())-1);
    return {begin,end};
}
