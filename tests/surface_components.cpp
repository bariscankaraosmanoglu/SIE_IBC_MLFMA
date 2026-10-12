#include "surface_components.h"
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <stdexcept>
#include <vector>
int num_nodes = 0;
int num_tris = 0;
int* tri_nodes = nullptr;
static bool rejects(int nodes, const std::vector<int>& triangles) {
    try { classify_open_surfaces(nodes, int(triangles.size()/3), triangles.data()); }
    catch (const std::invalid_argument&) { return true; }
    return false;
}
int main() {
    assert(classify_open_surfaces(0, 0, nullptr).empty());
    std::vector<int> tetrahedron{0,2,1, 0,1,3, 1,2,3, 2,0,3};
    assert(classify_open_surfaces(4, 4, tetrahedron.data()) == std::vector<unsigned char>(4, 0));
    std::vector<int> open = tetrahedron; open.resize(9);
    assert(classify_open_surfaces(4, 3, open.data()) == std::vector<unsigned char>(3, 1));
    std::vector<int> mixed = tetrahedron;
    for (int node : open) mixed.push_back(node + 4);
    std::vector<unsigned char> expected{0,0,0,0,1,1,1};
    assert(classify_open_surfaces(8, 7, mixed.data()) == expected);
    // Triangle order and winding do not affect boundary propagation.
    std::reverse(mixed.begin(), mixed.end()); std::reverse(expected.begin(), expected.end());
    assert(classify_open_surfaces(8, 7, mixed.data()) == expected);
    assert(rejects(5, {0,1,2, 1,0,3, 0,1,4}));
    assert(rejects(5, {0,1,2, 0,3,4}));
    std::vector<int> pinched = tetrahedron;
    for (int node : tetrahedron) pinched.push_back(node == 0 ? 0 : node+3);
    assert(rejects(7, pinched));
    assert(rejects(3, {0,0,1}));
    assert(rejects(3, {0,1,3}));
    assert(rejects(3, {0,1,-1}));
    printf("Surface component checks passed\n");
}
