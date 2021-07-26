//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Craig Schroeder, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VECTOR_POLICY
//#####################################################################
#ifndef __VECTOR_POLICY__
#define __VECTOR_POLICY__

#include <Matrices_And_Vectors/MATRIX_FORWARD.h>
#include <Matrices_And_Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class TV> struct VECTOR_POLICY;

//#####################################################################
// 0D
//#####################################################################
template<class T>
struct VECTOR_POLICY<VECTOR<T,0> >
{
    typedef PhysBAM::MATRIX<T,0> MATRIX;
    typedef MATRIX SYMMETRIC_MATRIX;
    typedef MATRIX DIAGONAL_MATRIX;
};
//#####################################################################
// 1D
//#####################################################################
template<class T>
struct VECTOR_POLICY<VECTOR<T,1> >
{
    typedef PhysBAM::MATRIX<T,1> MATRIX;
    typedef MATRIX SYMMETRIC_MATRIX;
    typedef MATRIX DIAGONAL_MATRIX;
    typedef PhysBAM::MATRIX<T,2> TRANSFORMATION_MATRIX;
    typedef VECTOR<T,0> SPIN;
};
//#####################################################################
// 2D
//#####################################################################
template<class T>
struct VECTOR_POLICY<VECTOR<T,2> >
{
    typedef PhysBAM::MATRIX<T,2> MATRIX;
    typedef PhysBAM::SYMMETRIC_MATRIX<T,2> SYMMETRIC_MATRIX;
    typedef PhysBAM::DIAGONAL_MATRIX<T,2> DIAGONAL_MATRIX;
    typedef PhysBAM::MATRIX<T,3> TRANSFORMATION_MATRIX;
    typedef VECTOR<T,1> SPIN;
};
//#####################################################################
// 3D
//#####################################################################
template<class T>
struct VECTOR_POLICY<VECTOR<T,3> >
{
    typedef PhysBAM::MATRIX<T,3> MATRIX;
    typedef PhysBAM::SYMMETRIC_MATRIX<T,3> SYMMETRIC_MATRIX;
    typedef PhysBAM::DIAGONAL_MATRIX<T,3> DIAGONAL_MATRIX;
    typedef PhysBAM::MATRIX<T,4> TRANSFORMATION_MATRIX;
    typedef VECTOR<T,3> SPIN;
};
//#####################################################################
}
#endif
