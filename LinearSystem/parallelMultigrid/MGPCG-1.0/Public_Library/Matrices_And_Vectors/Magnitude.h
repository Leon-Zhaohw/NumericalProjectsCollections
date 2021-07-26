//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function Magnitude 
//#####################################################################
#ifndef __Magnitude__
#define __Magnitude__

#include <cmath>
#include <cstdlib>
namespace PhysBAM{

template<class T,class T_ARRAY,class ID> class ARRAY_BASE;

using ::std::abs;

template<class TV>
inline typename TV::SCALAR Magnitude(const TV& v)
{return v.Magnitude();}

template<class TV>
inline typename TV::SCALAR Magnitude_Squared(const TV& v)
{return v.Magnitude_Squared();}

inline int Magnitude(const int a)
{return abs(a);}

inline float Magnitude(const float a)
{return abs(a);}

inline double Magnitude(const double a)
{return abs(a);}

inline float Magnitude_Squared(const float a)
{return a*a;}

inline double Magnitude_Squared(const double a)
{return a*a;}

}
#endif
