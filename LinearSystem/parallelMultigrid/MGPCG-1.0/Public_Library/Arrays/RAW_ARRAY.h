//#####################################################################
// Copyright 2007-2008, Geoffrey Irving, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RAW_ARRAY
//#####################################################################
#ifndef __RAW_ARRAY__
#define __RAW_ARRAY__

#include <Arrays/ARRAY_BASE.h>
#include <Math_Tools/exchange.h>
#include <Read_Write/READ_WRITE_FUNCTIONS.h>
#include <Utilities/EXCEPTIONS.h>
#include <Utilities/TYPE_UTILITIES.h>
#include <boost/format.hpp>
#include <boost/type_traits/add_const.hpp>
namespace PhysBAM{

template<class T,class ID> struct IS_ARRAY<RAW_ARRAY<T,ID> >:public mpl::true_{};
template<class T,class ID> struct IS_ARRAY_VIEW<RAW_ARRAY<T,ID> >:public mpl::true_{};
template<class T,class ID,class T_NEW> struct REBIND<RAW_ARRAY<T,ID>,T_NEW>{typedef RAW_ARRAY<typename IF<IS_CONST<T>::value,const T_NEW,T_NEW>::TYPE,ID> TYPE;};

template<class T,class ID> struct IS_CONST<RAW_ARRAY<const T,ID> >:public mpl::true_{}; // RAW_ARRAY<const T,ID> is equivalent to const RAW_ARRAY<const T,ID>

template<class T_ARRAY> struct CANONICALIZE_CONST_ARRAY<T_ARRAY,typename boost::enable_if<IS_BASE_OF<RAW_ARRAY<typename T_ARRAY::ELEMENT,typename T_ARRAY::INDEX>,T_ARRAY> >::type>
    :public CANONICALIZE_CONST_ARRAY<RAW_ARRAY<typename boost::add_const<typename T_ARRAY::ELEMENT>::type,typename T_ARRAY::INDEX> >{};

template<class T,class ID>
class RAW_ARRAY:public ARRAY_BASE<typename REMOVE_CONST<T>::TYPE,RAW_ARRAY<T,ID>,ID>
{
    struct UNUSABLE{};
    template<class S> struct COPY_CONST:public IF<IS_CONST<T>::value,typename boost::add_const<S>::type,S>{};
    typedef ARRAY_BASE<typename REMOVE_CONST<T>::TYPE,RAW_ARRAY<T,ID>,ID> BASE;
public:
    typedef typename REMOVE_CONST<T>::TYPE ELEMENT;typedef ID INDEX;
    typedef T& RESULT_TYPE;

protected:
    // m and base_pointer inherit constness of T
    typename COPY_CONST<ID>::TYPE m;
private:
    friend class RAW_ARRAY<typename IF<IS_CONST<T>::value,ELEMENT,const ELEMENT>::TYPE,ID>;
    typename COPY_CONST<T*>::TYPE base_pointer;

public:
    using BASE::Same_Array;

    RAW_ARRAY(const ID m,T* raw_data)
        :m(m),base_pointer(raw_data)
    {}

    RAW_ARRAY(const RAW_ARRAY<typename REMOVE_CONST<T>::TYPE,ID>& array)
        :m(array.m),base_pointer(array.base_pointer)
    {}

    template<class T_ARRAY>
    RAW_ARRAY(T_ARRAY& array,typename boost::enable_if<AND<IS_SAME<ELEMENT,typename T_ARRAY::ELEMENT>::value,NOT<IS_ARRAY_VIEW<T_ARRAY>::value>::value>,UNUSABLE>::type unusable=UNUSABLE())
        :m(array.Size()),base_pointer(array.Get_Array_Pointer())
    {}

    template<class T_ARRAY>
    RAW_ARRAY(T_ARRAY array,typename boost::enable_if<AND<IS_SAME<ELEMENT,typename T_ARRAY::ELEMENT>::value,IS_ARRAY_VIEW<T_ARRAY>::value>,UNUSABLE>::type unusable=UNUSABLE())
        :m(array.Size()),base_pointer(array.Get_Array_Pointer())
    {}

    ID Size() const
    {return m;}

    T& operator()(const ID i)
    {assert(ID(1)<=i && i<=m);return base_pointer[Value(i)-1];}

    const T& operator()(const ID i) const
    {assert(ID(1)<=i && i<=m);return base_pointer[Value(i)-1];}

    bool Valid_Index(const ID i) const
    {return ID(1)<=i && i<=m;}

    RAW_ARRAY& operator=(const RAW_ARRAY& source)
    {return BASE::operator=(source);}

    template<class T_ARRAY2>
    RAW_ARRAY& operator=(const T_ARRAY2& source)
    {return BASE::operator=(source);}

    T* Get_Array_Pointer()
    {return base_pointer;}

    const T* Get_Array_Pointer() const
    {return base_pointer;}

    void Exchange(RAW_ARRAY& other)
    {STATIC_ASSERT(!IS_CONST<T>::value); // make RAW_ARRAY<const T> equivalent to const RAW_ARRAY<const T>
    exchange(m,other.m);exchange(base_pointer,other.base_pointer);}

    static void Exchange_Arrays(RAW_ARRAY& array1,RAW_ARRAY& array2)
    {STATIC_ASSERT(!IS_CONST<T>::value); // make RAW_ARRAY<const T> equivalent to const RAW_ARRAY<const T>
    exchange(array1.m,array2.m);exchange(array1.base_pointer,array2.base_pointer);}

    RAW_ARRAY<typename REMOVE_CONST<T>::TYPE>& Const_Cast() const // return reference to allow Exchange
    {return reinterpret_cast<RAW_ARRAY<typename REMOVE_CONST<T>::TYPE>&>(const_cast<RAW_ARRAY&>(*this));}

    static bool Same_Array(const RAW_ARRAY& array1,const RAW_ARRAY& array2)
    {return array1.Get_Array_Pointer()==array2.Get_Array_Pointer();}

    template<class RW>
    void Read(std::istream& input)
    {ID read_size;Read_Binary<RW>(input,read_size);
    if(read_size!=Size()) throw READ_ERROR(str(boost::format("Expected size %d, read %d")%Size()%read_size));
    Read_Binary_Array<RW>(input,Get_Array_Pointer(),Value(Size()));}

    template<class RW>
    void Write(std::ostream& output) const
    {ID m=Size();Write_Binary<RW>(output,m);Write_Binary_Array<RW>(output,Get_Array_Pointer(),Value(m));}

//#####################################################################
};
template<class T> inline void exchange(RAW_ARRAY<T>& a,RAW_ARRAY<T>& b) // TODO: replace Exchange_Arrays with specialization of exchange
{STATIC_ASSERT((T)false);} // use Exchange_Arrays instead
}
#endif
