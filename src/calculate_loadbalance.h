#ifndef CALCULATE_LOADBALANCE_H
#define CALCULATE_LOADBALANCE_H

#include <utility>

void calculate_loadbalance();

// Call inside the parallel region with its actual thread ID and team size.
// The returned range is [first box, exclusive end box).
std::pair<int, int> near_box_range(int thread_id, int team_size);

#endif
