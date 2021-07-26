//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Sergey Koltakov, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function min 
//#####################################################################
#ifndef __min__
#define __min__

#include <algorithm>
namespace PhysBAM{

template<class T>
inline T min(const T a,const T b)
{return ::std::min(a,b);}

template<class T>
inline T min(const T a,const T b,const T c)
{return ::std::min(a,min(b,c));}

template<class T>
inline T min(const T a,const T b,const T c,const T d)
{return ::std::min(a,min(b,c,d));}

template<class T>
inline T min(const T a,const T b,const T c,const T d,const T e)
{return ::std::min(a,min(b,c,d,e));}

template<class T>
inline T min(const T a,const T b,const T c,const T d,const T e,const T f)
{return ::std::min(a,min(b,c,d,e,f));}

template<class T>
inline T min(const T a,const T b,const T c,const T d,const T e,const T f,const T g)
{return ::std::min(a,min(b,c,d,e,f,g));}

template<class T>
inline T min(const T a,const T b,const T c,const T d,const T e,const T f,const T g,const T h)
{return ::std::min(a,min(b,c,d,e,f,g,h));}

}
#endif
