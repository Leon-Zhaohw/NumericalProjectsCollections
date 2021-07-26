//#####################################################################
// Copyright 2007, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_0D
//#####################################################################
#ifndef __GRID_0D__
#define __GRID_0D__

#include <Grids/POLICY_UNIFORM.h>
#include <Arrays/ARRAYS_FORWARD.h>
#include <Geometry/BOX.h>
namespace PhysBAM{

template<class T> class UNIFORM_GRID_ITERATOR_NODE_0D;
template<class T> class UNIFORM_GRID_ITERATOR_CELL_0D;

template<class T>
class GRID_0D
{
    STATIC_ASSERT((IS_SCALAR<T>::value));
    typedef VECTOR<T,0> TV;typedef VECTOR<int,0> TV_INT;
public:
    enum REGION {WHOLE_REGION,GHOST_REGION,BOUNDARY_REGION,INTERIOR_REGION,BOUNDARY_INTERIOR_REGION}; // for iterators
    static const int dimension=0;
    typedef UNIFORM_GRID_ITERATOR_NODE_0D<T> NODE_ITERATOR;
    typedef UNIFORM_GRID_ITERATOR_CELL_0D<T> CELL_ITERATOR;
    typedef TV_INT VECTOR_INT;

    GRID_0D()
    {}

    GRID_0D(const TV_INT& counts,const BOX<VECTOR<T,0> >& box,const bool MAC_grid=false)
    {}

    GRID_0D<T>& operator=(const GRID_0D<T>& grid_input)
    {return *this;}
//#####################################################################
};
}
#endif
