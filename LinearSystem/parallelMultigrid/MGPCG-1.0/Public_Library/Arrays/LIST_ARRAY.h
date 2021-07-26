//#####################################################################
// Copyright 2004-2007, Ronald Fedkiw, Geoffrey Irving, Tamar Shinar, Eftychios Sifakis, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LIST_ARRAY
//#####################################################################
#ifndef __LIST_ARRAY__
#define __LIST_ARRAY__

#include <Arrays/ARRAY.h>
#include <Utilities/override.h>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/has_trivial_destructor.hpp>
namespace PhysBAM{

template<class T,class ID> struct IS_ARRAY<LIST_ARRAY<T,ID> >:public mpl::true_{};
template<class T,class ID> struct IS_ARRAY_VIEW<LIST_ARRAY<T,ID> >:public mpl::false_{};
template<class T,class ID> struct CANONICALIZE_CONST_ARRAY<LIST_ARRAY<T,ID> >:public FIRST<RAW_ARRAY<const T,ID> >{};

template<class T,class ID>
class LIST_ARRAY:public ARRAY_BASE<T,LIST_ARRAY<T,ID>,ID>
{
public:
    template<class T2> struct REBIND{typedef LIST_ARRAY<T2,ID> TYPE;};
    typedef T ELEMENT;typedef ID INDEX;
private:
    struct UNUSABLE{};
public:

    ARRAY<T,ID> array;
    ID m; // the current size of the array (array.m may be larger for elbow room)

    LIST_ARRAY()
        :m(0)
    {}

    explicit LIST_ARRAY(const ID m_input,const bool initialize_using_default_constructor=true)
        :array(m_input,initialize_using_default_constructor),m(m_input)
    {}

    LIST_ARRAY(const LIST_ARRAY& array_input)
        :array(array_input),m(array.m)
    {}

    template<class T_ARRAY>
    //explicit LIST_ARRAY(const T_ARRAY& array_input,typename boost::enable_if<IS_SAME<T,typename T_ARRAY::ELEMENT>,UNUSABLE>::type unused=UNUSABLE())
    explicit LIST_ARRAY(const T_ARRAY& array_input,typename boost::enable_if<boost::is_same<T,typename T_ARRAY::ELEMENT> >::type* =0)
        :array(array_input),m(array.m)
    {}

    LIST_ARRAY& operator=(const LIST_ARRAY& source)
    {ID source_m=source.m;
    if(array.m<source_m) array.Resize(source_m,false,false);
    else if(Same_Array(*this,source)) return *this;
    for(ID i(1);i<=source_m;i++) array(i)=source(i);
    m=source_m;return *this;}

    template<class T_ARRAY>
    LIST_ARRAY& operator=(const T_ARRAY& source)
    {STATIC_ASSERT_SAME(T,typename T_ARRAY::ELEMENT);
    ID source_m=source.Size();
    if(array.m<source_m) array.Resize(source_m,false,false);
    else if(Same_Array(*this,source)) return *this;
    for(ID i(1);i<=source_m;i++) array(i)=source(i);
    m=source_m;return *this;}

    ID Size() const
    {return m;}

    T& operator()(const ID i)
    {assert(i<=m);return array(i);}

    const T& operator()(const ID i) const
    {assert(i<=m);return array(i);}

    T* Get_Array_Pointer()
    {return array.Get_Array_Pointer();}

    const T* Get_Array_Pointer() const
    {return array.Get_Array_Pointer();}

    ID Max_Size() const
    {return array.m;}

    void Compact()
    {if(m<array.m) Exact_Resize(m);}

private:
    void Ensure_Enough_Space(const ID m_new,const bool copy_existing_elements=true) PHYSBAM_ALWAYS_INLINE
    {if(array.m<m_new) array.Resize_Grow(ID(4*Value(m_new)/3+2),false,copy_existing_elements);}

public:
    void Preallocate(const ID max_size)
    {if(array.m<max_size){
        ARRAY<T,ID> new_array(max_size,false);for(ID i(1);i<=m;i++) new_array(i)=array(i);
        ARRAY<T,ID>::Exchange_Arrays(array,new_array);}}

    void Resize(const ID m_new,const bool initialize_new_elements=true,const bool copy_existing_elements=true)
    {Ensure_Enough_Space(m_new,copy_existing_elements);if(initialize_new_elements && m_new>m) for(ID i=m+1;i<=m_new;i++) array(i)=T();
    m=m_new;}

    void Exact_Resize(const ID m_new,const bool initialize_new_elements=true) // zero elbow room
    {array.Resize(m_new,initialize_new_elements);m=array.m;}

    ID Append(const T& element) PHYSBAM_ALWAYS_INLINE
    {m++;Ensure_Enough_Space(m);array(m)=element;return m;}

    template<class T_ARRAY>
    void Append_Elements(const T_ARRAY& append_array)
    {STATIC_ASSERT_SAME(ELEMENT,typename T_ARRAY::ELEMENT);m+=Value(append_array.Size());Ensure_Enough_Space(m);
    for(typename T_ARRAY::INDEX i(1);i<=append_array.Size();i++) array(m-Value(append_array.Size())+Value(i))=append_array(i);}

    void Append_Unique(const T& element)
    {for(ID i(1);i<=m;i++) if(array(i)==element) return;Append(element);}

    template<class T_ARRAY>
    void Append_Unique_Elements(const T_ARRAY& append_array)
    {STATIC_ASSERT_SAME(T,typename T_ARRAY::ELEMENT);
    typename T_ARRAY::INDEX append_m=append_array.Size();for(typename T_ARRAY::INDEX i(1);i<=append_m;i++) Append_Unique(append_array(i));}

    void Remove_End()
    {assert(m>ID());m--;}

    void Remove_Index(const ID index) // preserves ordering of remaining elements
    {assert(ID(1)<=index && index<=m);for(ID i=index;i<m;i++) array(i)=array(i+1);Remove_End();}

    void Remove_Index_Lazy(const ID index)
    {assert(ID(1)<=index && index<=m);
    if(index<m) array(index)=array(m);
    Remove_End();}

    void Remove_All() // if elements are non-primitive this may waste memory
    {m=ID();}

    void Clean_Memory()
    {Exact_Resize(ID());}

    void Delete_Pointers_And_Clean_Memory() // only valid if T is a pointer type
    {for(ID i(1);i<=m;i++) delete array(i);Clean_Memory();}

    void Insert(const T& element,const ID index)
    {m++;Ensure_Enough_Space(m);for(ID i=m;i>index;i--) array(i)=array(i-1);array(index)=element;}

    T Pop()
    {Remove_End();return array(m+1);}

    RAW_ARRAY<const T> Pop_Elements(const int count) // return value should be copied immediately, not kept around
    {BOOST_MPL_ASSERT((boost::has_trivial_destructor<T>));
    assert(m-count>=ID());m-=count;
    return RAW_ARRAY<const T>(count,array.Get_Array_Pointer()+m);}

    static void Exchange_Arrays(LIST_ARRAY<T,ID>& a,LIST_ARRAY<T,ID>& b)
    {ARRAY<T,ID>::Exchange_Arrays(a.array,b.array);exchange(a.m,b.m);}

    template<class RW>
    void Read(std::istream& input_stream)
    {array.template Read<RW>(input_stream);m=array.m;}

    template<class RW>
    void Write(std::ostream& output_stream) const
    {array.template Write_Prefix<RW>(output_stream,m);}

//#####################################################################
};
}
#endif
