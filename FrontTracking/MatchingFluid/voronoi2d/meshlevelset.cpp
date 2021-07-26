/*
 *  meshlevelset.cpp
 *  eltopo2d_project
 *
 *  Created by tyson on 11/11/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "meshlevelset.h"
#include "dynamicsurface.h"

// simplest/stupidest version:
// get the distance for each point, then use ray casting to get the sign

void make_mesh_level_set( const eltopo2d::DynamicSurface& surface,
                          const std::vector<Vec3ui> &mesh_triangles, const std::vector<Vec2d> &mesh_vertices,
                          std::vector<double>& phi )
{
   phi.clear();
   phi.resize( mesh_vertices.size(), 1e+30 );
   
   for ( unsigned int i = 0; i < mesh_vertices.size(); ++i )
   {
      // get absolute distance
      phi[i] = surface.get_distance_to_surface( mesh_vertices[i] );
      
      // get sign
      if ( surface.point_is_inside_volume( mesh_vertices[i] ) )
      {
         phi[i] = -phi[i];
      }
   }

}


