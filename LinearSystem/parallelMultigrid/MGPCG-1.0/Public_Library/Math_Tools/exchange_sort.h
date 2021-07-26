//#####################################################################
// Copyright 2002, Robert Bridson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function exchange_sort
//#####################################################################
//
// exchanges the values passed to be in sorted (ascending) order 
//               
//#####################################################################
#ifndef __exchange_sort__
#define __exchange_sort__

#include <Math_Tools/exchange.h>
namespace PhysBAM{

template<class T>
inline void exchange_sort(T& a,T& b)
{if(b<a) exchange(a,b);}

template<class T>
inline void exchange_sort(T& a,T& b,T& c)
{if(b<a) exchange(a,b);if(c<b) exchange(b,c);if(b<a) exchange(a,b);}

template<class T>
inline void exchange_sort(T& a,T& b,T& c,T& d)
{if(b<a) exchange(a,b);if(d<c) exchange(c,d);if(c<a) exchange(a,c);if(b>d) exchange(b,d);if(c<b) exchange(b,c);}

}
#endif

