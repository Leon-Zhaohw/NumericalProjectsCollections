//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GEOMETRY_POLICY 
//#####################################################################
#ifndef __GEOMETRY_POLICY__
#define __GEOMETRY_POLICY__

#include <Geometry/GEOMETRY_FORWARD.h>
#include <Collisions_And_Interactions/COLLISIONS_FORWARD.h>
namespace PhysBAM{

template<class TV> struct GEOMETRY_POLICY;

//#####################################################################
// 0D
//#####################################################################
template<class T>
struct GEOMETRY_POLICY<VECTOR<T,0> >
{
    struct UNUSABLE{};
};
//#####################################################################
// 1D
//#####################################################################
template<class T>
struct GEOMETRY_POLICY<VECTOR<T,1> >
{
    struct UNUSABLE{};
    typedef VECTOR<T,1> POINT;
    typedef BOX<VECTOR<T,1> > ORIENTED_BOX;
    typedef PhysBAM::RAY<VECTOR<T,1> > RAY;
    typedef POINT_SIMPLEX_1D<T> HYPERPLANE;
    typedef POINT_SIMPLEX_1D<T> POINT_SIMPLEX;
    typedef POINT_SIMPLICES_1D<T> POINT_SIMPLICES;
    typedef SEGMENT_1D<T> SEGMENT;

    typedef PhysBAM::SEGMENTED_CURVE<VECTOR<T,1> > SEGMENTED_CURVE;
    typedef UNUSABLE TRIANGULATED_OBJECT;
};
//#####################################################################
// 2D
//#####################################################################
template<class T>
struct GEOMETRY_POLICY<VECTOR<T,2> >
{
    typedef POINT_2D<T> POINT;
    typedef PhysBAM::RAY<VECTOR<T,2> > RAY;
    typedef PhysBAM::ORIENTED_BOX<VECTOR<T,2> > ORIENTED_BOX;
    typedef LINE_2D<T> HYPERPLANE;
    typedef SEGMENT_2D<T> SEGMENT;
    typedef TRIANGLE_2D<T> TRIANGLE;

    typedef SEGMENTED_CURVE_2D<T> SEGMENTED_CURVE;
    typedef TRIANGULATED_AREA<T> TRIANGULATED_OBJECT;
    typedef TRIANGLE_HIERARCHY_2D<T> TRIANGLE_HIERARCHY;
};
//#####################################################################
// 3D
//#####################################################################
template<class T>
struct GEOMETRY_POLICY<VECTOR<T,3> >
{
public:
    typedef VECTOR<T,3> POINT;
    typedef PhysBAM::RAY<VECTOR<T,3> > RAY;
    typedef PhysBAM::ORIENTED_BOX<VECTOR<T,3> > ORIENTED_BOX;
    typedef PLANE<T> HYPERPLANE;
    typedef SEGMENT_3D<T> SEGMENT;
    typedef TRIANGLE_3D<T> TRIANGLE;

    typedef PhysBAM::SEGMENTED_CURVE<VECTOR<T,3> > SEGMENTED_CURVE;
    typedef TRIANGULATED_SURFACE<T> TRIANGULATED_OBJECT;
    typedef PhysBAM::TRIANGLE_HIERARCHY<T> TRIANGLE_HIERARCHY;
};
//#####################################################################
}
#endif
