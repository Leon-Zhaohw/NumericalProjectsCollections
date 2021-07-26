//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header VECTOR_FORWARD
//#####################################################################
#ifndef __VECTOR_FORWARD__
#define __VECTOR_FORWARD__

#include <Read_Write/READ_WRITE_FORWARD.h>
#include <Utilities/STATIC_ASSERT.h>
#include <Utilities/TYPE_UTILITIES.h>
namespace PhysBAM{

template<class T,class T_VECTOR> class VECTOR_BASE;
template<class T,int d> class VECTOR;

template<class T> class VECTOR_ND;
template<class T> class SPARSE_VECTOR_ND;

struct ZERO;
template<class T> class INTERVAL;
template<class T> class COMPLEX;
template<class T> class QUATERNION;
template<class TV> class ROTATION;
template<class T> class TETRAHEDRAL_GROUP;
template<class TV> class FRAME;
template<class TV> class TWIST;

template<class T_VECTOR,class T_NEW> struct REBIND;
template<class T,int d,class T_NEW> struct REBIND<VECTOR<T,d>,T_NEW>{typedef VECTOR<T_NEW,d> TYPE;};

template<class T,class T_EXPRESSION> class VECTOR_EXPRESSION;
template<class T_VECTOR1,class T_VECTOR2> class VECTOR_SUM;
template<class T_VECTOR1,class T_VECTOR2> class VECTOR_DIFFERENCE;
template<class T_VECTOR,class T2> class VECTOR_SCALE;
template<class T_VECTOR> class VECTOR_NEGATION;

template<class T> struct IS_VECTOR:public mpl::false_{};
template<class T,int d> struct IS_VECTOR<VECTOR<T,d> >:public mpl::true_{};

template<class T> struct IS_SCALAR_BLOCK;
template<class T> struct IS_SCALAR_VECTOR_SPACE;
template<class T,int d> struct IS_SCALAR_BLOCK<VECTOR<T,d> >:public AND<(d>0),IS_SCALAR_BLOCK<T>::value>{};
template<class T,int d> struct IS_SCALAR_VECTOR_SPACE<VECTOR<T,d> >:public AND<(d>0),IS_SCALAR_VECTOR_SPACE<T>::value>{};
template<class T,int d,class RW> struct IS_BINARY_IO_SAFE<VECTOR<T,d>,RW>:public AND<(d>0),IS_BINARY_IO_SAFE<T,RW>::value>{};

}
#endif
