#ifndef SURFACE_COMPONENTS_H
#define SURFACE_COMPONENTS_H

#include <vector>

// Validate edge incidence before RWG construction. The returned flags are
// indexed by triangle; all triangles in a component with a boundary are open.
std::vector<unsigned char> classify_open_surfaces(int node_count, int triangle_count,
                                                const int* triangle_nodes);
void generate_surface_components();
bool triangle_on_open_surface(int triangle);

#endif
