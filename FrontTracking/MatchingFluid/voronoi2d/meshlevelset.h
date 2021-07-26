/*
 *  meshlevelset.h
 *  eltopo2d_project
 *
 *  Created by tyson on 11/11/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <vec.h>

namespace eltopo2d
{
class DynamicSurface;
}

void make_mesh_level_set( const eltopo2d::DynamicSurface& surface,
                          const std::vector<Vec3ui> &mesh_triangles, 
                          const std::vector<Vec2d> &mesh_vertices,
                          std::vector<double>& phi );