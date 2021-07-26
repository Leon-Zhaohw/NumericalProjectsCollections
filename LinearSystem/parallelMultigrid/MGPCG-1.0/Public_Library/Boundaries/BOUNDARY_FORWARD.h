//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header BOUNDARY_FORWARD
//#####################################################################
#ifndef __BOUNDARY_FORWARD__
#define __BOUNDARY_FORWARD__

namespace PhysBAM{

template<class T,class T2> class BOUNDARY;
template<class T_GRID,class T2,int d=1> class BOUNDARY_UNIFORM;
template<class T_GRID,class T2> class BOUNDARY_DYADIC;
template<class T_GRID,class T2,int d=1> class BOUNDARY_RLE;

template<class T_GRID,class T2> class BOUNDARY_LINEAR_EXTRAPOLATION;
template<class T_GRID,class T2> class BOUNDARY_MAC_GRID_PERIODIC;
template<class T_GRID,class T2> class BOUNDARY_PERIODIC;
template<class T_GRID,class T2> class BOUNDARY_REFLECTION;
template<class T_GRID,class T2> class BOUNDARY_REFLECTION_WATER;
template<class T_GRID,class T2> class BOUNDARY_CONSTANT_EXTRAPOLATION_COLLIDABLE;
template<class T_GRID,int d=T_GRID::dimension+2> class BOUNDARY_EULER_EQUATIONS_SOLID_WALL_PERIODIC;
template<class T_GRID> class BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP;
template<class T_GRID> class BOUNDARY_MAC_GRID_SOLID_WALL_SLIP;
template<class T_GRID> class BOUNDARY_PHI_WATER;
template<class T_GRID> class BOUNDARY_SOLID_WALL_SLIP_OUTFLOW;
template<class T> class BOUNDARY_EULER_EQUATIONS_CYLINDRICAL;
template<class T> class BOUNDARY_EULER_EQUATIONS_SPHERICAL;
template<class TV> class BOUNDARY_SOLID_WALL_SLIP;

template<class T_BOUNDARY,class T_NEW> struct REBIND;
template<class T_BOUNDARY,int d_new> struct REBIND_LENGTH;
template<class T_GRID,class T2,int d,class T_NEW> struct REBIND<BOUNDARY_UNIFORM<T_GRID,T2,d>,T_NEW>{typedef BOUNDARY_UNIFORM<T_GRID,T_NEW,d> TYPE;};
template<class T_GRID,class T2,class T_NEW> struct REBIND<BOUNDARY_DYADIC<T_GRID,T2>,T_NEW>{typedef BOUNDARY_DYADIC<T_GRID,T_NEW> TYPE;};
template<class T_GRID,class T2,class T_NEW> struct REBIND<BOUNDARY_RLE<T_GRID,T2>,T_NEW>{typedef BOUNDARY_RLE<T_GRID,T_NEW> TYPE;};
template<class T_GRID,class T2,int d,int d_new> struct REBIND_LENGTH<BOUNDARY_UNIFORM<T_GRID,T2,d>,d_new>{typedef BOUNDARY_UNIFORM<T_GRID,T2,d_new> TYPE;};
template<class T_GRID,class T2,int d_new> struct REBIND_LENGTH<BOUNDARY_DYADIC<T_GRID,T2>,d_new>{typedef BOUNDARY_DYADIC<T_GRID,T2> TYPE;};
template<class T_GRID,class T2,int d_new> struct REBIND_LENGTH<BOUNDARY_RLE<T_GRID,T2>,d_new>{typedef BOUNDARY_RLE<T_GRID,T2> TYPE;};
}
#endif
