//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PROJECTED_ARRAY
//#####################################################################
#ifndef __PROJECTED_ARRAY__
#define __PROJECTED_ARRAY__

#include <Arrays/ARRAY_BASE.h>
#include <Matrices_And_Vectors/VECTOR_FORWARD.h>
#include <boost/pointee.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/type_traits/remove_reference.hpp>
namespace PhysBAM{

template<class T_ARRAY,class T_PROJECTOR> struct IS_ARRAY<PROJECTED_ARRAY<T_ARRAY,T_PROJECTOR> >:public mpl::true_{};
template<class T_ARRAY,class T_PROJECTOR> struct IS_ARRAY_VIEW<PROJECTED_ARRAY<T_ARRAY,T_PROJECTOR> >:public mpl::true_{};

template<class T_ARRAY> struct PROJECTED_ARRAY_BASE{};
template<int d> struct PROJECTED_ARRAY_BASE<VECTOR<int,d>&>{enum {m=VECTOR<int,d>::m};};
template<class T_ARRAY,class T_PROJECTOR> struct PROJECTED_ARRAY_BASE<PROJECTED_ARRAY<T_ARRAY,T_PROJECTOR>&>:public PROJECTED_ARRAY_BASE<T_ARRAY>{};

template<class T_ARRAY,class T_PROJECTOR> struct PROJECTED_ARRAY_ELEMENT
{
    typedef typename boost::remove_const<typename boost::remove_reference<typename T_PROJECTOR::template RESULT<typename T_ARRAY::ELEMENT>::TYPE>::type>::type TYPE;
};

//#####################################################################
// Class PROJECTED_ARRAY
//#####################################################################
template<class T_ARRAY,class T_PROJECTOR>
class PROJECTED_ARRAY:public PROJECTED_ARRAY_BASE<T_ARRAY>,public ARRAY_BASE<typename PROJECTED_ARRAY_ELEMENT<T_ARRAY,T_PROJECTOR>::TYPE,PROJECTED_ARRAY<T_ARRAY,T_PROJECTOR> >,
    private T_PROJECTOR
{
    typedef typename PROJECTED_ARRAY_ELEMENT<T_ARRAY,T_PROJECTOR>::TYPE T;
    typedef typename T_PROJECTOR::template RESULT<typename T_ARRAY::ELEMENT>::TYPE RESULT_NONCONST;
    typedef typename T_PROJECTOR::template RESULT<const typename T_ARRAY::ELEMENT>::TYPE RESULT_CONST;
public:
    typedef T ELEMENT;typedef typename T_ARRAY::INDEX INDEX;

    T_ARRAY& array;

    PROJECTED_ARRAY(T_ARRAY& array)
        :array(array)
    {
        STATIC_ASSERT((IS_EMPTY<T_PROJECTOR>::value));
    }

    PROJECTED_ARRAY(T_ARRAY& array,const T_PROJECTOR& projector)
        :T_PROJECTOR(projector),array(array)
    {}

    PROJECTED_ARRAY(const PROJECTED_ARRAY<typename REMOVE_CONST<T_ARRAY>::TYPE,T_PROJECTOR>& projected_array)
        :T_PROJECTOR(projected_array.Projector()),array(projected_array.array)
    {}

    const T_PROJECTOR& Projector() const
    {return *this;}

    INDEX Size() const
    {return array.Size();}

    typename IF<IS_CONST<T_ARRAY>::value,RESULT_CONST,RESULT_NONCONST>::TYPE operator()(const INDEX i)
    {return T_PROJECTOR::operator()(array(i));}

    RESULT_CONST operator()(const INDEX i) const
    {return T_PROJECTOR::operator()(array(i));}

    PROJECTED_ARRAY operator=(const PROJECTED_ARRAY& source)
    {return ARRAY_BASE<T,PROJECTED_ARRAY>::operator=(source);}

    template<class T_ARRAY2>
    PROJECTED_ARRAY operator=(const T_ARRAY2& source)
    {return ARRAY_BASE<T,PROJECTED_ARRAY>::operator=(source);}
};
//#####################################################################
// Class FIELD_PROJECTOR
//#####################################################################
template<class T_STRUCT,class T_FIELD,T_FIELD T_STRUCT::* field>
struct FIELD_PROJECTOR
{
    template<class U> struct RESULT:public IF<IS_CONST<U>::value,const T_FIELD&,T_FIELD&>{};

    const T_FIELD& operator()(const T_STRUCT& element) const
    {return element.*field;}

    T_FIELD& operator()(T_STRUCT& element) const
    {return element.*field;}
};
//#####################################################################
// Class INDEX_PROJECTOR
//#####################################################################
struct INDEX_PROJECTOR
{
    template<class T_ARRAY> struct RESULT:public IF<IS_CONST<T_ARRAY>::value,const typename T_ARRAY::ELEMENT&,typename T_ARRAY::ELEMENT&>{};
    int index;

    INDEX_PROJECTOR(const int index)
        :index(index)
    {}

    template<class T_ARRAY> const typename T_ARRAY::ELEMENT& operator()(const T_ARRAY& array) const
    {return array(index);}

    template<class T_ARRAY> typename T_ARRAY::ELEMENT& operator()(T_ARRAY& array) const
    {return array(index);}
};
//#####################################################################
// Class POINTER_PROJECTOR
//#####################################################################
struct POINTER_PROJECTOR // turns shared_ptr<T>, etc., into T*
{
    template<class T_POINTER> struct RESULT{typedef typename boost::pointee<typename boost::remove_const<T_POINTER>::type>::type* TYPE;};

    template<class T_POINTER> typename RESULT<T_POINTER>::TYPE operator()(const T_POINTER& pointer) const
    {return &*pointer;}
};
//#####################################################################
}
#endif
