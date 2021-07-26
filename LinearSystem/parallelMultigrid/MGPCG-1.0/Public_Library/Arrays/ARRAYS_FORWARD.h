//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Nipun Kwatra, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header ARRAYS_FORWARD
//#####################################################################
#ifndef __ARRAYS_FORWARD__
#define __ARRAYS_FORWARD__

namespace PhysBAM{

template<class T,class T_ARRAY,class ID=int> class ARRAY_BASE;
template<class T,class T_ARRAY,class ID=int> class ARRAY_EXPRESSION;
template<class T,class T_ARRAYS> class ARRAYS_ND_BASE;

template<class T,int d> class VECTOR;
template<class T,class ID=int> class ARRAY;
template<class T,class ID=int> class LIST_ARRAY;
template<class T_ARRAY,class T_INDICES=LIST_ARRAY<int>&> class INDIRECT_ARRAY;
template<class T,int length=1> class ARRAYS_1D;
template<class T,int length=1> class ARRAYS_2D;
template<class T,int length=1> class ARRAYS_3D;
template<class ID=int> class IDENTITY_MAP;
template<class T,class ID=int> class CONSTANT_MAP;

template<class T,int length=1> class FACE_ARRAYS_1D;
template<class T,int length=1> class FACE_ARRAYS_2D;
template<class T,int length=1> class FACE_ARRAYS_3D;

template<class T> class BLOCKED_ARRAY_2D;
template<class T> class BLOCKED_ARRAY_3D;

template<class T,class ID=int> class RAW_ARRAY;
template<class T_ARRAY,class T_PROJECTOR> class PROJECTED_ARRAY;
template<class T_STRUCT,class T_FIELD,T_FIELD T_STRUCT::* field> struct FIELD_PROJECTOR;
struct INDEX_PROJECTOR;
struct POINTER_PROJECTOR;

class FLOOD_FILL_1D;
class FLOOD_FILL_2D;
class FLOOD_FILL_3D;
class FLOOD_FILL_GRAPH;
template<class T> class FLOOD_FILL_OCTREE;
template<class T> class FLOOD_FILL_QUADTREE;

template<class T_ARRAY,class T_NEW> struct REBIND;
template<class T,class T_NEW> struct REBIND<ARRAY<T>,T_NEW>{typedef ARRAY<T_NEW> TYPE;};
template<class T,class T_NEW> struct REBIND<LIST_ARRAY<T>,T_NEW>{typedef LIST_ARRAY<T_NEW> TYPE;};
template<class T,class T_NEW> struct REBIND<ARRAYS_1D<T>,T_NEW>{typedef ARRAYS_1D<T_NEW> TYPE;};
template<class T,class T_NEW> struct REBIND<ARRAYS_2D<T>,T_NEW>{typedef ARRAYS_2D<T_NEW> TYPE;};
template<class T,class T_NEW> struct REBIND<ARRAYS_3D<T>,T_NEW>{typedef ARRAYS_3D<T_NEW> TYPE;};
template<class T,class T_NEW> struct REBIND<FACE_ARRAYS_1D<T>,T_NEW>{typedef FACE_ARRAYS_1D<T_NEW> TYPE;};
template<class T,class T_NEW> struct REBIND<FACE_ARRAYS_2D<T>,T_NEW>{typedef FACE_ARRAYS_2D<T_NEW> TYPE;};
template<class T,class T_NEW> struct REBIND<FACE_ARRAYS_3D<T>,T_NEW>{typedef FACE_ARRAYS_3D<T_NEW> TYPE;};

template<class T_ARRAY,int length_new> struct REBIND_LENGTH;
template<class T,int length_new> struct REBIND_LENGTH<ARRAY<T>,length_new>{typedef ARRAY<VECTOR<T,length_new> > TYPE;};
template<class T> struct REBIND_LENGTH<ARRAY<T>,1>{typedef ARRAY<T> TYPE;};
template<class T,int length_new> struct REBIND_LENGTH<ARRAYS_1D<T>,length_new>{typedef ARRAYS_1D<T,length_new> TYPE;};
template<class T,int length_new> struct REBIND_LENGTH<ARRAYS_2D<T>,length_new>{typedef ARRAYS_2D<T,length_new> TYPE;};
template<class T,int length_new> struct REBIND_LENGTH<ARRAYS_3D<T>,length_new>{typedef ARRAYS_3D<T,length_new> TYPE;};
template<class T,int length_new> struct REBIND_LENGTH<FACE_ARRAYS_1D<T>,length_new>{typedef FACE_ARRAYS_1D<T,length_new> TYPE;};
template<class T,int length_new> struct REBIND_LENGTH<FACE_ARRAYS_2D<T>,length_new>{typedef FACE_ARRAYS_2D<T,length_new> TYPE;};
template<class T,int length_new> struct REBIND_LENGTH<FACE_ARRAYS_3D<T>,length_new>{typedef FACE_ARRAYS_3D<T,length_new> TYPE;};

template<class T_ARRAY> struct ARRAY_RESULT_TYPE{typedef typename T_ARRAY::RESULT_TYPE TYPE;};
template<class T_ARRAY> struct ARRAY_RESULT_TYPE<const T_ARRAY>{typedef typename T_ARRAY::CONST_RESULT_TYPE TYPE;};

}
#endif
