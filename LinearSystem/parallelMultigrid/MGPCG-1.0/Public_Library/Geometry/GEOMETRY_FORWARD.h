//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GEOMETRY_FORWARD 
//#####################################################################
#ifndef __GEOMETRY_FORWARD__
#define __GEOMETRY_FORWARD__

namespace PhysBAM{

template<class TV> class BOX;
template<class TV> class RAY;
template<class T,int d> class VECTOR;
template<class TV> class ORIENTED_BOX;
template<class T> class POINT_SIMPLEX_1D;
template<class T> class LINE_2D;
template<class T> class PLANE;
template<class T> class SEGMENT_1D;
template<class T> class SEGMENT_2D;
template<class T> class SEGMENT_3D;
template<class T> class TRIANGLE_2D;
template<class T> class TRIANGLE_3D;
template<class T> class TETRAHEDRON;
template<class T> class POINT_SIMPLICES_1D;
template<class T> class SEGMENTED_CURVE_2D;
template<class TV> class SEGMENTED_CURVE;
template<class T> class TRIANGULATED_AREA;
template<class T> class TRIANGULATED_SURFACE;
template<class T> class TETRAHEDRALIZED_VOLUME;
template<class T> class HEXAHEDRALIZED_VOLUME;
template<class T> class CYLINDER;
template<class T> class RING;
template<class TV> class SPHERE;
template<class T> class TORUS;
template<class TV> class IMPLICIT_OBJECT;
template<class TV,class TRANSFORM> class IMPLICIT_OBJECT_TRANSFORMED;
template<class T> class POINT_2D;
template<class TV> class BOUNDED_HORIZONTAL_PLANE;

enum POINT_SIMPLEX_COLLISION_TYPE {POINT_SIMPLEX_NO_COLLISION,POINT_SIMPLEX_COLLISION_ENDS_OUTSIDE,POINT_SIMPLEX_COLLISION_ENDS_INSIDE,POINT_SIMPLEX_UNKNOWN_COLLISION};

}
#endif
