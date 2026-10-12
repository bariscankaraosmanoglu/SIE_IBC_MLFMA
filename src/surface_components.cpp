#include "surface_components.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <stdexcept>
#include <string>

extern int num_nodes;
extern int num_tris;
extern int* tri_nodes;

namespace {
struct Edge { int first, second, triangle; };
std::vector<unsigned char> open_triangles;
}

std::vector<unsigned char> classify_open_surfaces(int node_count, int triangle_count,
                                                const int* triangle_nodes) {
    if (node_count < 0 || triangle_count < 0 || (triangle_count && !triangle_nodes))
        throw std::invalid_argument("Invalid surface dimensions");
    std::vector<Edge> edges;
    edges.reserve(size_t(triangle_count) * 3);
    std::vector<int> parent(triangle_count), rank(triangle_count, 0);
    std::iota(parent.begin(), parent.end(), 0);
    // Separate union-find for triangle corners checks each vertex's incident
    // triangle fan; edge incidence alone misses surfaces pinched at a vertex.
    std::vector<size_t> corner_parent(size_t(triangle_count) * 3);
    std::iota(corner_parent.begin(), corner_parent.end(), size_t(0));
    auto corner_root = [&](size_t corner) {
        while (corner_parent[corner] != corner) {
            corner_parent[corner] = corner_parent[corner_parent[corner]];
            corner = corner_parent[corner];
        }
        return corner;
    };
    auto root = [&](int triangle) {
        while (parent[triangle] != triangle) {
            parent[triangle] = parent[parent[triangle]];
            triangle = parent[triangle];
        }
        return triangle;
    };
    for (int triangle = 0; triangle < triangle_count; ++triangle) {
        const int* nodes = triangle_nodes + size_t(triangle) * 3;
        for (int j = 0; j < 3; ++j) {
            if (nodes[j] < 0 || nodes[j] >= node_count)
                throw std::invalid_argument("Invalid node index on triangle " + std::to_string(triangle));
            int a = nodes[j], b = nodes[(j + 1) % 3];
            if (a == b)
                throw std::invalid_argument("Repeated node on triangle " + std::to_string(triangle));
            edges.push_back({std::min(a, b), std::max(a, b), triangle});
        }
    }
    std::sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) {
        return a.first != b.first ? a.first < b.first : a.second < b.second;
    });
    std::vector<unsigned char> boundary(triangle_count, 0);
    for (size_t begin = 0; begin < edges.size();) {
        size_t end = begin + 1;
        while (end < edges.size() && edges[end].first == edges[begin].first &&
               edges[end].second == edges[begin].second) ++end;
        if (end - begin > 2)
            throw std::invalid_argument("Nonmanifold edge " + std::to_string(edges[begin].first) +
                                        " " + std::to_string(edges[begin].second));
        if (end - begin == 1) boundary[edges[begin].triangle] = 1;
        else {
            for (int node : {edges[begin].first, edges[begin].second}) {
                size_t a = size_t(edges[begin].triangle) * 3;
                size_t b = size_t(edges[begin + 1].triangle) * 3;
                while (triangle_nodes[a] != node) ++a;
                while (triangle_nodes[b] != node) ++b;
                a = corner_root(a); b = corner_root(b);
                corner_parent[a] = b;
            }
            int a = root(edges[begin].triangle), b = root(edges[begin + 1].triangle);
            if (a != b) {
                if (rank[a] < rank[b]) std::swap(a, b);
                parent[b] = a;
                if (rank[a] == rank[b]) ++rank[a];
            }
        }
        begin = end;
    }
    std::vector<size_t> vertex_fan(node_count, corner_parent.size());
    for (size_t corner = 0; corner < corner_parent.size(); ++corner) {
        int node = triangle_nodes[corner];
        size_t fan = corner_root(corner);
        if (vertex_fan[node] == corner_parent.size()) vertex_fan[node] = fan;
        else if (vertex_fan[node] != fan)
            throw std::invalid_argument("Nonmanifold vertex " + std::to_string(node));
    }
    std::vector<unsigned char> open(triangle_count, 0);
    for (int triangle = 0; triangle < triangle_count; ++triangle)
        if (boundary[triangle]) open[root(triangle)] = 1;
    for (int triangle = 0; triangle < triangle_count; ++triangle)
        boundary[triangle] = open[root(triangle)];
    return boundary;
}

void generate_surface_components() {
    try {
        open_triangles = classify_open_surfaces(num_nodes, num_tris, tri_nodes);
    } catch (const std::invalid_argument& error) {
        std::fprintf(stderr, "Invalid surface topology: %s\n", error.what());
        std::exit(EXIT_FAILURE);
    }
}

bool triangle_on_open_surface(int triangle) {
    if (triangle < 0 || size_t(triangle) >= open_triangles.size()) {
        std::fprintf(stderr, "Surface components must be initialized before CFIE coefficients\n");
        std::exit(EXIT_FAILURE);
    }
    return open_triangles[triangle] != 0;
}
