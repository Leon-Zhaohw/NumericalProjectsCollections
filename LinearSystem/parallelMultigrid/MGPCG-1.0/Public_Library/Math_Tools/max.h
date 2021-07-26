//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Sergey Koltakov, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function max  
//#####################################################################
#ifndef __max__
#define __max__

#include <algorithm>
namespace PhysBAM{

template<class T>
inline T max(const T a,const T b)
{return ::std::max(a,b);}

template<class T>
inline T max(const T a,const T b,const T c)
{return ::std::max(a,max(b,c));}

template<class T>
inline T max(const T a,const T b,const T c,const T d)
{return ::std::max(a,max(b,c,d));}

template<class T>
inline T max(const T a,const T b,const T c,const T d,const T e)
{return ::std::max(a,max(b,c,d,e));}

template<class T>
inline T max(const T a,const T b,const T c,const T d,const T e,const T f)
{return ::std::max(a,max(b,c,d,e,f));}

template<class T>
inline T max(const T a,const T b,const T c,const T d,const T e,const T f,const T g)
{return ::std::max(a,max(b,c,d,e,f,g));}

template<class T>
inline T max(const T a,const T b,const T c,const T d,const T e,const T f,const T g,const T h)
{return ::std::max(a,max(b,c,d,e,f,g,h));}

template<class T>
inline T max(const T a,const T b,const T c,const T d,const T e,const T f,const T g,const T h,const T i)
{return ::std::max(a,max(b,c,d,e,f,g,h,i));}

}
#endif
