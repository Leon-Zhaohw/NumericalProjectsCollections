//#####################################################################
// Copyright 2006-2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header COLLISIONS_FORWARD
//#####################################################################
#ifndef __COLLISIONS_FORWARD__
#define __COLLISIONS_FORWARD__

#include <Particles/PARTICLES_FORWARD.h>
namespace PhysBAM{

template<class TV> class BOX_HIERARCHY;
template<class TV> class SEGMENT_HIERARCHY;
template<class T> class TRIANGLE_HIERARCHY;
template<class T> class TRIANGLE_HIERARCHY_2D;
template<class T> class TETRAHEDRON_HIERARCHY;
template<class T> class PARTICLE_PARTITION;
template<class TV,class T_PARTICLE=SOLIDS_PARTICLE<TV> > class PARTICLE_HIERARCHY;

template<class TV> class COLLISION_BODY;
template<class TV> class COLLISION_BODY_LIST;
template<class T> class TETRAHEDRON_COLLISION_BODY;
template<class TV> class COLLISION_PENALTY_FORCES;
template<class TV> class DEFORMABLE_OBJECT_FLUID_COLLISIONS;
template<class TV> class COLLISION_BODY_IMPULSE_ACCUMULATOR;
template<class TV> class RIGID_DEFORMABLE_COLLISIONS;

template<class T_COLLISION_BODY> class COLLISION_BODY_SPATIAL_PARTITION;
template<class T_COLLISION_BODY> class COLLISION_BODY_SPATIAL_PARTITION_UNION;
enum SPATIAL_PARTITION_VOXEL_SIZE_HEURISTIC {SPATIAL_PARTITION_SCENE_SIZE,SPATIAL_PARTITION_MAX_BOX_SIZE,SPATIAL_PARTITION_AVERAGE_BOX_SIZE};

template<class TV> struct POINT_FACE_REPULSION_PAIR;
template<class TV> struct EDGE_EDGE_REPULSION_PAIR;
template<class TV> class SEGMENT_REPULSIONS_AND_COLLISIONS_GEOMETRY;
template<class TV> class TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY;
template<class TV> class SEGMENT_REPULSIONS;
template<class TV> class TRIANGLE_REPULSIONS;
template<class TV> class SEGMENT_COLLISIONS;
template<class TV> class TRIANGLE_COLLISIONS;

template<class TV> struct PRECOMPUTE_PROJECT_POINT_FACE;
template<class TV> struct PRECOMPUTE_PROJECT_EDGE_EDGE;

}
#endif
