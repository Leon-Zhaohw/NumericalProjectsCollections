//#####################################################################
// Copyright 2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_EXPRESSION
//#####################################################################
#ifndef __ARRAY_EXPRESSION__
#define __ARRAY_EXPRESSION__

#include <Arrays/ARRAYS_FORWARD.h>
#include <Utilities/TYPE_UTILITIES.h>
#include <boost/type_traits/remove_reference.hpp>
#include <iterator>
namespace PhysBAM{

template<class T_ARRAY> class ARRAY_STL_ITERATOR
{
    typedef typename T_ARRAY::INDEX ID;

    T_ARRAY* array;
    ID index;
public:
    typedef std::random_access_iterator_tag iterator_category;
    typedef typename ARRAY_RESULT_TYPE<T_ARRAY>::TYPE reference;
    typedef typename boost::remove_reference<reference>::type value_type;
    typedef int difference_type;
    typedef value_type* pointer;

    ARRAY_STL_ITERATOR()
        :array(0),index(0)
    {}

    ARRAY_STL_ITERATOR(T_ARRAY& array,const ID index)    
        :array(&array),index(index)
    {}

    void operator++()
    {index++;}

    bool operator==(const ARRAY_STL_ITERATOR& other)
    {return index==other.index;} // assume array==other.array

    bool operator!=(const ARRAY_STL_ITERATOR& other)
    {return index!=other.index;} // assume array==other.array

    typename ARRAY_RESULT_TYPE<T_ARRAY>::TYPE operator*() const
    {return (*array)(index);}
};

template<class T,class T_ARRAY,class ID> // ID = int
class ARRAY_EXPRESSION:public ARRAY_BASE<T,T_ARRAY,ID>
{
    using ARRAY_BASE<T,T_ARRAY,ID>::Derived;
public:
    typedef const T RESULT_TYPE;
    typedef const T CONST_RESULT_TYPE;
    typedef ARRAY_STL_ITERATOR<T_ARRAY> iterator;
    typedef ARRAY_STL_ITERATOR<const T_ARRAY> const_iterator;

    ARRAY_STL_ITERATOR<T_ARRAY> begin() // for stl
    {return ARRAY_STL_ITERATOR<T_ARRAY>(Derived(),ID(1));}

    ARRAY_STL_ITERATOR<const T_ARRAY> begin() const // for stl
    {return ARRAY_STL_ITERATOR<const T_ARRAY>(Derived(),ID(1));}

    ARRAY_STL_ITERATOR<T_ARRAY> end() // for stl
    {return ARRAY_STL_ITERATOR<T_ARRAY>(Derived(),ID(Derived().Size()+1));}

    ARRAY_STL_ITERATOR<const T_ARRAY> end() const // for stl
    {return ARRAY_STL_ITERATOR<const T_ARRAY>(Derived(),ID(Derived().Size()+1));}

//#####################################################################
};
}
#endif
