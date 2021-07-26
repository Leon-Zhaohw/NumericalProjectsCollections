//#####################################################################
// Copyright 2007-2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAYS_ND_BASE
//#####################################################################
#ifndef __ARRAYS_ND_BASE__
#define __ARRAYS_ND_BASE__

#include <Arrays/RAW_ARRAY.h>
namespace PhysBAM{

template<class T,class T_ARRAYS>
class ARRAYS_ND_BASE
{
public:
    typedef T ELEMENT;
    typedef typename SCALAR_POLICY<T>::TYPE SCALAR;

    RAW_ARRAY<T> array; // one-dimensional data storage

protected:
    ARRAYS_ND_BASE()
        :array(0,0)
    {}

    explicit ARRAYS_ND_BASE(const int size)
        :array(size,new T[size])
    {}

    ~ARRAYS_ND_BASE()
    {delete[] array.Get_Array_Pointer();}
public:

    T_ARRAYS& Derived()
    {return static_cast<T_ARRAYS&>(*this);}

    const T_ARRAYS& Derived() const
    {return static_cast<const T_ARRAYS&>(*this);}

    void Delete_Pointers_And_Clean_Memory() // only valid if T is a pointer type
    {for(int i=1;i<=array.Size();i++) delete array(i);Derived().Clean_Memory();}

    bool operator==(const T_ARRAYS& v) const
    {return T_ARRAYS::Equal_Dimensions(Derived(),v) && array==v.array;}

    bool operator!=(const T_ARRAYS& v) const
    {return !(*this==v);}

    T_ARRAYS& operator+=(const T_ARRAYS& v)
    {assert(T_ARRAYS::Equal_Dimensions(Derived(),v));array+=v.array;return Derived();}

    T_ARRAYS& operator+=(const T& a)
    {array+=a;return Derived();}

    T_ARRAYS& operator-=(const T_ARRAYS& v)
    {assert(T_ARRAYS::Equal_Dimensions(Derived(),v));array-=v.array;return Derived();}

    T_ARRAYS& operator-=(const T& a)
    {array-=a;return Derived();}

    template<class T2,class T_ARRAYS_T2>
    T_ARRAYS& operator*=(const ARRAYS_ND_BASE<T2,T_ARRAYS_T2>& v)
    {assert(T_ARRAYS::Equal_Dimensions(Derived(),v.Derived()));array*=v.array;return Derived();}

    template<class T2> typename boost::enable_if<OR<IS_SCALAR<T2>::value,IS_SAME<T,T2>::value>,T_ARRAYS&>::type
    operator*=(const T2 a)
    {array*=a;return Derived();}

    template<class T2>
    T_ARRAYS& operator/=(const T2 a)
    {return *this*=Inverse(a);}

    int Number_True() const
    {return array.Number_True();}

    void Fill(const T& constant)
    {array.Fill(constant);}

    static void Copy(const T_ARRAYS& old_copy,T_ARRAYS& new_copy)
    {assert(T_ARRAYS::Equal_Dimensions(old_copy,new_copy));
    RAW_ARRAY<T>::Copy(old_copy.array,new_copy.array);}

    template<class T2>
    static void Copy(const T2 constant,const T_ARRAYS& old_copy,T_ARRAYS& new_copy)
    {assert(T_ARRAYS::Equal_Dimensions(old_copy,new_copy));
    new_copy.array=constant*old_copy.array;}

    template<class T2>
    static void Copy(const T2 c1,const T_ARRAYS& v1,const T_ARRAYS& v2,T_ARRAYS& result)
    {assert(T_ARRAYS::Equal_Dimensions(v1,v2)&&T_ARRAYS::Equal_Dimensions(v2,result));
    result.array=c1*v1.array+v2.array;}

    template<class T2>
    static void Copy(const T2 c1,const T_ARRAYS& v1,const T2 c2,const T_ARRAYS& v2,T_ARRAYS& result)
    {assert(T_ARRAYS::Equal_Dimensions(v1,v2)&&T_ARRAYS::Equal_Dimensions(v2,result));
    result.array=c1*v1.array+c2*v2.array;}

    template<class T2>
    static void Copy(const T2 c1,const T_ARRAYS& v1,const T2 c2,const T_ARRAYS& v2,const T2 c3,const T_ARRAYS& v3,T_ARRAYS& result)
    {assert(T_ARRAYS::Equal_Dimensions(v1,v2)&&T_ARRAYS::Equal_Dimensions(v2,v3)&&T_ARRAYS::Equal_Dimensions(v3,result));
    result.array=c1*v1.array+c2*v2.array+c3*v3.array;}

    void Clamp_Below(const T& value)
    {array.Clamp_Below(value);}

    T Average() const
    {return array.Average();}

    T Max() const
    {return array.Max();}

    T Maxabs() const
    {return array.Maxabs();}

    T Maxmag() const
    {return array.Maxmag();}

    T Min() const
    {return array.Min();}

    T Minmag() const
    {return array.Minmag();}

    T Sum() const
    {return array.Sum();}

    T Sumabs() const
    {return array.Sumabs();}

    T Componentwise_Maxabs() const
    {return array.Componentwise_Maxabs();}

    static T Dot_Product(const T_ARRAYS& a1,const T_ARRAYS& a2)
    {STATIC_ASSERT(T_ARRAYS::length==1);assert(T_ARRAYS::Equal_Dimensions(a1,a2));
    return RAW_ARRAY<T>::Dot_Product(a1.array,a2.array);}

    static double Dot_Product_Double_Precision(const T_ARRAYS& a1,const T_ARRAYS& a2)
    {STATIC_ASSERT(T_ARRAYS::length==1);assert(T_ARRAYS::Equal_Dimensions(a1,a2));
	return RAW_ARRAY<T>::Dot_Product_Double_Precision(a1.array,a2.array);}

    template<class TV>
    static typename SCALAR_POLICY<TV>::TYPE Maximum_Magnitude(const ARRAYS_ND_BASE<TV,T_ARRAYS>& a)
    {STATIC_ASSERT((IS_SAME<T,TV>::value && T_ARRAYS::length==1));
    return a.array.Maximum_Magnitude();}

//#####################################################################
};
}
#endif

