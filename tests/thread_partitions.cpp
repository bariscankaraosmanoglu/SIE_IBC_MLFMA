#include "calculate_loadbalance.h"
#include <omp.h>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

int num_threads=8;
int num_mlfmalevels=2;
int *num_filledboxes;
int *num_children;
int *boxes_near;
long int size_Znear;

static void check(const std::vector<int>& counts){
    const int box_count=static_cast<int>(counts.size());
    std::vector<int> levels{1,box_count,0,1,0,box_count};
    std::vector<int> children(box_count+1);
    std::vector<int> neighbors((box_count+1)*27,-1);
    num_filledboxes=levels.data();
    num_children=children.data();
    boxes_near=neighbors.data();
    size_Znear=0;
    for(int b=1;b<=box_count;++b){
        children[b]=counts[b-1];
        // Self and previous box exercise unequal interaction weights.
        neighbors[b*27]=b;
        if(b>1) neighbors[b*27+1]=b-1;
        size_Znear+=std::int64_t(children[b])*children[b];
        if(b>1) size_Znear+=std::int64_t(children[b])*children[b-1];
    }
    calculate_loadbalance();
    for(int team : {1,2,3,8,16,32,64,512}){
        std::vector<int> visits(box_count);
        int previous_end=1;
        for(int worker=0;worker<team;++worker){
            const auto range=near_box_range(worker,team);
            assert(range.first==previous_end);
            assert(1<=range.first && range.first<=range.second && range.second<=box_count+1);
            for(int b=range.first;b<range.second;++b) ++visits[b-1];
            previous_end=range.second;
        }
        assert(previous_end==box_count+1);
        for(int visits_to_box : visits) assert(visits_to_box==1);
    }
    // Large team sizes also exercise target arithmetic without a huge allocation.
    for(int worker : {0,1,1073741823,2147483646}){
        const auto range=near_box_range(worker,2147483647);
        assert(1<=range.first && range.first<=range.second && range.second<=box_count+1);
    }
    for(int parallel : {0,1}){
        std::vector<int> visits(box_count);
        int actual_team=0;
        #pragma omp parallel num_threads(32) if(parallel)
        {
            #pragma omp single
            actual_team=omp_get_num_threads();
            const auto range=near_box_range(omp_get_thread_num(),omp_get_num_threads());
            for(int b=range.first;b<range.second;++b){
                #pragma omp atomic update
                visits[b-1]++;
            }
        }
        for(int visits_to_box : visits) assert(visits_to_box==1);
        if(!parallel) assert(actual_team==1);
    }
}

int main(int argc,char** argv){
    if(argc==2){
        num_threads=std::atoi(argv[1]);
        calculate_loadbalance(); // Invalid counts must fail before accessing geometry.
        return 2;
    }
    check({});
    check({0});
    check({1});
    check({0,0,0,0});
    check({1,0,1,0,1});
    check({1,1,10000,1,1});
    check(std::vector<int>(32,1)); // Original radiation partition overflow case.
    check({100000,100000}); // Interaction counts exceed 32 bits.
    check({2000000000,1}); // total*team_index would overflow 64 bits.
    std::cout<<"Partition coverage and actual-team checks passed.\n";
}
