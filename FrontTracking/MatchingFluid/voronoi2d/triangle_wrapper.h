#ifndef TRIANGLE_WRAPPER_H
#define TRIANGLE_WRAPPER_H

#include "vec.h"
#include <vector>

void compute_delaunay(const std::vector<Vec2f>& points, std::vector<Vec3ui>& tris);

#endif