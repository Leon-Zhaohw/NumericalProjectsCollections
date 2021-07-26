//#####################################################################
// Copyright 2006, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_POLICY 
//#####################################################################
#ifndef __GRID_POLICY__
#define __GRID_POLICY__

#include <Matrices_And_Vectors/VECTOR_POLICY.h>
namespace PhysBAM{

template<class T> class GRID_0D;
template<class T> class GRID_1D;
template<class T> class GRID_2D;
template<class T> class GRID_3D;
template<class T> class BINTREE_GRID;
template<class T> class QUADTREE_GRID;
template<class T> class OCTREE_GRID;
template<class T> class RLE_GRID_1D;
template<class T> class RLE_GRID_2D;
template<class T> class RLE_GRID_3D;
template<class T> class UNIFORM_GRID_ITERATOR_1D;
template<class T> class UNIFORM_GRID_ITERATOR_2D;
template<class T> class UNIFORM_GRID_ITERATOR_3D;
template<class T> class UNIFORM_GRID_ITERATOR_NODE_1D;
template<class T> class UNIFORM_GRID_ITERATOR_NODE_2D;
template<class T> class UNIFORM_GRID_ITERATOR_NODE_3D;
template<class T> class UNIFORM_GRID_ITERATOR_CELL_1D;
template<class T> class UNIFORM_GRID_ITERATOR_CELL_2D;
template<class T> class UNIFORM_GRID_ITERATOR_CELL_3D;
template<class T> class UNIFORM_GRID_ITERATOR_FACE_1D;
template<class T> class UNIFORM_GRID_ITERATOR_FACE_2D;
template<class T> class UNIFORM_GRID_ITERATOR_FACE_3D;
template<class TV> struct DYADIC_TAG;
template<class TV> struct UNIFORM_TAG;
template<class TV> struct RLE_TAG;

template<class TV> struct GRID_POLICY;

//#####################################################################
// 0D
//#####################################################################
template<class T>
struct GRID_POLICY<VECTOR<T,0> >
{
    typedef GRID_0D<T> UNIFORM_GRID;
};
//#####################################################################
// 1D
//#####################################################################
template<class T>
struct GRID_POLICY<VECTOR<T,1> >
{
    typedef GRID_1D<T> UNIFORM_GRID;
    typedef BINTREE_GRID<T> DYADIC_GRID;
    typedef RLE_GRID_1D<T> RLE_GRID;
};
//#####################################################################
// 2D
//#####################################################################
template<class T>
struct GRID_POLICY<VECTOR<T,2> >
{
    typedef GRID_2D<T> UNIFORM_GRID;
    typedef QUADTREE_GRID<T> DYADIC_GRID;
    typedef RLE_GRID_2D<T> RLE_GRID;
};
//#####################################################################
// 3D
//#####################################################################
template<class T>
struct GRID_POLICY<VECTOR<T,3> >
{
    typedef GRID_3D<T> UNIFORM_GRID;
    typedef OCTREE_GRID<T> DYADIC_GRID;
    typedef RLE_GRID_3D<T> RLE_GRID;
};
//#####################################################################
}
#endif
