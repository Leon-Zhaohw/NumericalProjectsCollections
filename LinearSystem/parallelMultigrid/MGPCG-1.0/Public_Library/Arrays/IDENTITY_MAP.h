//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Craig Schroeder, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IDENTITY_MAP
//#####################################################################
#ifndef __IDENTITY_MAP__
#define __IDENTITY_MAP__

#include <Arrays/ARRAY_EXPRESSION.h>
#include <cassert>
namespace PhysBAM{

template<class ID> struct IS_ARRAY<IDENTITY_MAP<ID> >:public mpl::true_{};
template<class ID> struct IS_ARRAY_VIEW<IDENTITY_MAP<ID> >:public mpl::true_{};

template<class ID> // ID=int
class IDENTITY_MAP:public ARRAY_EXPRESSION<ID,IDENTITY_MAP<ID>,ID>
{
public:
    typedef ID ELEMENT;
    typedef ID INDEX;
private:
    ID m;
public:

    explicit IDENTITY_MAP(const ID m)
        :m(m)
    {}

    ID Size() const
    {return m;}

    ID operator()(const ID i) const
    {assert(ID(1)<=i && i<=m);return i;}

//#####################################################################
};
}
#endif
