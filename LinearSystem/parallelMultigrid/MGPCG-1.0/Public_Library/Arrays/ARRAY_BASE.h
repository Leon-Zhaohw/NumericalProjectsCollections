//#####################################################################
// Copyright 2004-2007, Kevin Der, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_BASE
//#####################################################################
#ifndef __ARRAY_BASE__
#define __ARRAY_BASE__

#include <Arrays/ARRAYS_FORWARD.h>
#include <Arrays/ARRAY_DIFFERENCE.h>
#include <Arrays/ARRAY_LEFT_MULTIPLE.h>
#include <Arrays/ARRAY_NEGATION.h>
#include <Arrays/ARRAY_PLUS_SCALAR.h>
#include <Arrays/ARRAY_PRODUCT.h>
#include <Arrays/ARRAY_SUM.h>
#include <Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <Data_Structures/ELEMENT_ID.h>
#include <Math_Tools/max.h>
#include <Math_Tools/maxabs.h>
#include <Math_Tools/maxmag.h>
#include <Math_Tools/min.h>
#include <Math_Tools/minmag.h>
#include <Math_Tools/Inverse.h>
#include <Matrices_And_Vectors/Magnitude.h>
#include <Matrices_And_Vectors/Dot_Product.h>
#include <Matrices_And_Vectors/SCALAR_POLICY.h>
#include <Utilities/TYPE_UTILITIES.h>
#include <Utilities/STATIC_ASSERT.h>
#include <boost/utility/enable_if.hpp>
#include <iostream>
namespace PhysBAM{

template<class T,class T_ARRAY,class ID> class ARRAY_BASE;

template<class T_ARRAY,class ENABLER=void> struct CANONICALIZE_CONST_ARRAY:public FIRST<T_ARRAY>{};

template<class T_ARRAY1,class T_ARRAY2> struct SAME_ARRAY_CANONICAL{static bool Same_Array(const T_ARRAY1& array1,const T_ARRAY2& array2)
{STATIC_ASSERT(!IS_SAME<T_ARRAY1,T_ARRAY2>::value);return false;}};

template<class T_ARRAY> struct SAME_ARRAY_CANONICAL<T_ARRAY,T_ARRAY>{static bool Same_Array(const T_ARRAY& array1,const T_ARRAY& array2)
{return T_ARRAY::Same_Array(array1,array2);}};

template<class TA1,class TA2> struct SAME_ARRAY:public SAME_ARRAY_CANONICAL<typename CANONICALIZE_CONST_ARRAY<TA1>::TYPE,typename CANONICALIZE_CONST_ARRAY<TA2>::TYPE>{};

template<class T,class T_ARRAY,class ID>
class ARRAY_BASE
{
    struct UNUSABLE{};
public:
    typedef T ELEMENT;
    typedef ID INDEX;
    typedef T value_type; // for stl
    typedef T* iterator; // for stl
    typedef const T* const_iterator; // for stl
    typedef int difference_type; // for stl

    typedef T& RESULT_TYPE;
    typedef const T& CONST_RESULT_TYPE;
    typedef typename SCALAR_POLICY<T>::TYPE SCALAR;

protected:
    ARRAY_BASE(){}
    ARRAY_BASE(const ARRAY_BASE&){}
    ~ARRAY_BASE(){}
public:

    T_ARRAY& Derived()
    {return static_cast<T_ARRAY&>(*this);}

    const T_ARRAY& Derived() const
    {return static_cast<const T_ARRAY&>(*this);}

protected:
    template<class T_ARRAY2> struct IS_ARRAY_BASE:public mpl::false_{};
    template<class T2,class T_ARRAY2> struct IS_ARRAY_BASE<ARRAY_BASE<T2,T_ARRAY2,ID> >:public mpl::true_{};
public:

    template<class T_ARRAY1> static bool Same_Array(const T_ARRAY1& array1,const T_ARRAY1& array2)
    {return &array1==&array2;}

    template<class T_ARRAY1,class T_ARRAY2> static bool
    Same_Array(const T_ARRAY1& array1,const T_ARRAY2& array2)
    {return SAME_ARRAY<T_ARRAY1,T_ARRAY2>::Same_Array(array1,array2);}

protected:
    T_ARRAY& operator=(const ARRAY_BASE& source)
    {T_ARRAY& self=Derived();ID m=self.Size();const T_ARRAY& source_=source.Derived();assert(m==source_.Size());
    if(!T_ARRAY::Same_Array(self,source_)) for(ID i(1);i<=m;i++) self(i)=source_(i);
    return self;}

    template<class T_ARRAY1>
    T_ARRAY& operator=(const T_ARRAY1& source)
    {STATIC_ASSERT_SAME(T,typename T_ARRAY1::ELEMENT);
    T_ARRAY& self=Derived();ID m=self.Size();assert(m==source.Size());
    if(!T_ARRAY::Same_Array(self,source)) for(ID i(1);i<=m;i++) self(i)=source(i);
    return self;}
public:

    template<class T_INDICES>
    INDIRECT_ARRAY<T_ARRAY,T_INDICES&> Subset(const T_INDICES& indices)
    {return INDIRECT_ARRAY<T_ARRAY,T_INDICES&>(Derived(),indices);}

    template<class T_INDICES>
    INDIRECT_ARRAY<const T_ARRAY,T_INDICES&> Subset(const T_INDICES& indices) const
    {return INDIRECT_ARRAY<const T_ARRAY,T_INDICES&>(Derived(),indices);}

    INDIRECT_ARRAY<T_ARRAY,IDENTITY_MAP<> > Prefix(const ID prefix_size)
    {assert(prefix_size<=Derived().Size());return INDIRECT_ARRAY<T_ARRAY,IDENTITY_MAP<> >(Derived(),IDENTITY_MAP<>(prefix_size));}

    INDIRECT_ARRAY<const T_ARRAY,IDENTITY_MAP<> > Prefix(const ID prefix_size) const
    {assert(prefix_size<=Derived().Size());return INDIRECT_ARRAY<const T_ARRAY,IDENTITY_MAP<> >(Derived(),IDENTITY_MAP<>(prefix_size));}

private:
    typedef typename IF<IS_CLASS<T>::value,T,UNUSABLE>::TYPE T_IF_CLASS;
public:

    template<class T_FIELD,T_FIELD T_IF_CLASS::* field>
    PROJECTED_ARRAY<T_ARRAY,FIELD_PROJECTOR<T_IF_CLASS,T_FIELD,field> > Project()
    {return PROJECTED_ARRAY<T_ARRAY,FIELD_PROJECTOR<ELEMENT,T_FIELD,field> >(Derived());}

    template<class T_FIELD,T_FIELD T_IF_CLASS::* field>
    PROJECTED_ARRAY<const T_ARRAY,FIELD_PROJECTOR<T_IF_CLASS,T_FIELD,field> > Project() const
    {return PROJECTED_ARRAY<const T_ARRAY,FIELD_PROJECTOR<ELEMENT,T_FIELD,field> >(Derived());}

    PROJECTED_ARRAY<T_ARRAY,INDEX_PROJECTOR> Project(const ID index)
    {return PROJECTED_ARRAY<T_ARRAY,INDEX_PROJECTOR>(Derived(),INDEX_PROJECTOR(index));}

    PROJECTED_ARRAY<const T_ARRAY,INDEX_PROJECTOR> Project(const ID index) const
    {return PROJECTED_ARRAY<const T_ARRAY,INDEX_PROJECTOR>(Derived(),INDEX_PROJECTOR(index));}

private:
    template<class S> struct ELEMENT_OF{typedef typename S::ELEMENT TYPE;};
    typedef typename IF<IS_VECTOR<T>::value,ELEMENT_OF<T>,FIRST<UNUSABLE> >::TYPE::TYPE ELEMENT_OF_T;
public:

    RAW_ARRAY<ELEMENT_OF_T> Flattened() // valid only for contiguous arrays of VECTOR<T,d>
    {T_ARRAY& self=Derived();return RAW_ARRAY<typename T::ELEMENT>(T::m*self.Size(),self.Get_Array_Pointer()->begin());}

    RAW_ARRAY<const ELEMENT_OF_T> Flattened() const // valid only for contiguous arrays of VECTOR<T,d>
    {const T_ARRAY& self=Derived();return RAW_ARRAY<const typename T::ELEMENT>(T::m*self.Size(),self.Get_Array_Pointer()->begin());}

    PROJECTED_ARRAY<const T_ARRAY,POINTER_PROJECTOR> Raw_Pointers() const
    {return PROJECTED_ARRAY<const T_ARRAY,POINTER_PROJECTOR>(Derived());}

    template<class T_ARRAY1>
    bool operator==(const T_ARRAY1& v) const
    {STATIC_ASSERT_SAME(T,typename T_ARRAY1::ELEMENT);
    const T_ARRAY& self=Derived();ID m=self.Size();
    if(m!=v.Size()) return false;for(ID i(1);i<=m;i++) if(self(i)!=v(i)) return false;return true;}

    template<class T_ARRAY1>
    bool operator!=(const T_ARRAY1& v) const
    {return !(*this==v);}

    template<class T_ARRAY1>
    T_ARRAY& operator+=(const ARRAY_BASE<T,T_ARRAY1,ID>& v)
    {T_ARRAY& self=Derived();ID m=self.Size();const T_ARRAY1& v_=v.Derived();assert(m==v_.Size());
    for(ID i(1);i<=m;i++) self(i)+=v_(i);return self;}

    T_ARRAY& operator+=(const T& a) // this could be merged with the version below if windows wasn't broken
    {T_ARRAY& self=Derived();ID m=self.Size();
    for(ID i(1);i<=m;i++) self(i)+=a;return self;}

    template<class T2> typename boost::enable_if<IS_SCALAR<T2>,T_ARRAY&>::type
    operator+=(const T2& a)
    {T_ARRAY& self=Derived();ID m=self.Size();
    for(ID i(1);i<=m;i++) self(i)+=a;return self;}

    template<class T_ARRAY1>
    T_ARRAY& operator-=(const ARRAY_BASE<T,T_ARRAY1,ID>& v)
    {T_ARRAY& self=Derived();ID m=self.Size();const T_ARRAY1& v_=v.Derived();assert(m==v_.Size());
    for(ID i(1);i<=m;i++) self(i)-=v_(i);return self;}

    template<class T2> typename boost::enable_if<OR<IS_SCALAR<T2>::value,IS_SAME<T,T2>::value>,T_ARRAY&>::type
    operator-=(const T2& a)
    {T_ARRAY& self=Derived();ID m=self.Size();
    for(ID i(1);i<=m;i++) self(i)-=a;return self;}

    template<class T2,class T_ARRAY_T2>
    T_ARRAY& operator*=(const ARRAY_BASE<T2,T_ARRAY_T2,ID>& v)
    {T_ARRAY& self=Derived();ID m=self.Size();const T_ARRAY_T2& v_=v.Derived();assert(m==v_.Size());
    for(ID i(1);i<=m;i++) self(i)*=v_(i);return self;}

    template<class T2> typename boost::enable_if<OR<IS_SCALAR<T2>::value,IS_SAME<T,T2>::value>,T_ARRAY&>::type
    operator*=(const T2& a)
    {T_ARRAY& self=Derived();ID m=self.Size();
    for(ID i(1);i<=m;i++) self(i)*=a;return self;}

    template<class T2,class T_ARRAY_T2>
    T_ARRAY& operator/=(const ARRAY_BASE<T2,T_ARRAY_T2,ID>& a)
    {T_ARRAY& self=Derived();ID m=self.Size();const T_ARRAY_T2& a_=a.Derived();assert(m==a_.Size());
    for(ID i(1);i<=m;i++){assert(a_(i));self(i)/=a_(i);}return self;}

    template<class T2> typename boost::enable_if<OR<IS_SCALAR<T2>::value,IS_SAME<T,T2>::value>,T_ARRAY&>::type
    operator/=(const T2& a)
    {return *this*=Inverse(a);}

    T& Last()
    {T_ARRAY& self=Derived();return self(self.Size());}

    const T& Last() const
    {const T_ARRAY& self=Derived();return self(self.Size());}

    ID Find(const T& element) const
    {const T_ARRAY& self=Derived();ID m=self.Size();
    for(ID i(1);i<=m;i++) if(self(i)==element) return i;return 0;}

    bool Find(const T& element,ID& index) const // returns the first occurence of an element in an array
    {return Find(element,1,index);}

    bool Find(const T& element,const ID start_index,ID& index) const // returns the first occurence after start_index of an element in an array
    {const T_ARRAY& self=Derived();ID m=self.Size();
    for(ID i=start_index;i<=m;i++) if(self(i)==element){index=i;return true;}return false;}

    bool Contains(const T& element) const
    {const T_ARRAY& self=Derived();ID m=self.Size();
    for(ID i(1);i<=m;i++) if(self(i)==element) return true;return false;}

    bool Contains_Only(const T& element) const
    {const T_ARRAY& self=Derived();ID m=self.Size();
    for(ID i(1);i<=m;i++) if(self(i)!=element) return false;return true;}

    int Count_Matches(const T& value) const
    {const T_ARRAY& self=Derived();ID m=self.Size();
    int count=0;for(ID i(1);i<=m;i++) if(self(i)==value) count++;return count;}

    int Number_True() const
    {STATIC_ASSERT_SAME(T,bool);return Count_Matches(true);}

    int Number_False() const
    {STATIC_ASSERT_SAME(T,bool);return Count_Matches(false);}

    void Get_Unique(LIST_ARRAY<T>& array) const
    {const T_ARRAY& self=Derived();HASHTABLE<T> hash(Value(self.Size())*3/2);hash.Set_All(self);hash.Get_Keys(array);}

    void Prune_Duplicates()
    {T_ARRAY& self=Derived();HASHTABLE<T> hash(Value(self.Size())*3/2);int j=0;for(int i=1;i<=self.Size();i++) if(hash.Set(self(i))) self(++j)=self(i);self.Resize(j);}

    void Fill(const T& constant)
    {T_ARRAY& self=Derived();ID m=self.Size();for(ID i(1);i<=m;i++) self(i)=constant;}

    template<class T_ARRAY1,class T_ARRAY2>
    static void Copy(const T_ARRAY1& old_copy,T_ARRAY2& new_copy,typename boost::disable_if<IS_ARRAY_VIEW<T_ARRAY2>,UNUSABLE>::type unusable=UNUSABLE())
    {STATIC_ASSERT_SAME(T,typename T_ARRAY1::ELEMENT);
    STATIC_ASSERT_SAME(T,typename T_ARRAY2::ELEMENT);
    if(T_ARRAY::Same_Array(old_copy,new_copy)) return;
    ID m=old_copy.Size();assert(m==new_copy.Size());
    for(ID i(1);i<=m;i++) new_copy(i)=old_copy(i);}

    template<class T_ARRAY1,class T_ARRAY2>
    static typename boost::enable_if<IS_ARRAY_VIEW<T_ARRAY2> >::type
    Copy(const T_ARRAY1& old_copy,T_ARRAY2 new_copy)
    {STATIC_ASSERT_SAME(T,typename T_ARRAY1::ELEMENT);
    STATIC_ASSERT_SAME(T,typename T_ARRAY2::ELEMENT);
    ID m=old_copy.Size();assert(m==new_copy.Size());
    for(ID i(1);i<=m;i++) new_copy(i)=old_copy(i);}

    template<class T2,class T_ARRAY1,class T_ARRAY2>
    static void Copy(const T2 constant,const T_ARRAY1& array,T_ARRAY2& result)
    {Copy(constant*array,result);}

    template<class T2,class T_ARRAY1,class T_ARRAY2,class T_ARRAY3>
    static void Copy(const T2 c1,const T_ARRAY1& v1,const T_ARRAY2& v2,T_ARRAY3& result)
    {Copy(c1*v1+v2,result);}

    template<class T2,class T_ARRAY1,class T_ARRAY2,class T_ARRAY3>
    static void Copy(const T2 c1,const T_ARRAY1& v1,const T2 c2,const T_ARRAY2& v2,T_ARRAY3& result)
    {Copy(c1*v1+c2*v2,result);}

    template<class T2,class T_ARRAY1,class T_ARRAY2,class T_ARRAY3,class T_ARRAY4>
    static void Copy(const T2 c1,const T_ARRAY1& v1,const T2 c2,const T_ARRAY2& v2,const T2 c3,const T_ARRAY3& v3,T_ARRAY4& result)
    {Copy(c1*v1+c2*v2+c3*v3,result);}

    void Copy_Or_Fill(const T& constant) // for occasional templatization purposes
    {Fill(constant);}

    template<class T_ARRAY1>
    void Copy_Or_Fill(const T_ARRAY1& old_copy) // for occasional templatization purposes
    {Copy(old_copy,Derived());}

    static void Get(T_ARRAY& new_copy,const T_ARRAY& old_copy)
    {if(&old_copy!=&new_copy) Copy(old_copy.Prefix(new_copy.Size()),new_copy);}

    static void Put(const T_ARRAY& old_copy,T_ARRAY& new_copy)
    {if(&old_copy!=&new_copy) Copy(old_copy,new_copy.Prefix(old_copy.Size()));}

    template<class T2>
    static void Put(const T2 constant,const T_ARRAY& old_copy,T_ARRAY& new_copy)
    {Copy(constant*old_copy,new_copy.Prefix(old_copy.Size()));}

    void Clamp_Below(const T& value)
    {T_ARRAY& self=Derived();ID m=self.Size();for(ID i(1);i<=m;i++) self(i)=clamp_min(self(i),value);}

    T Max() const
    {const T_ARRAY& self=Derived();T result=self(ID(1));ID m=self.Size();for(ID i(2);i<=m;i++) result=PhysBAM::max(result,self(i));return result;}

    T Maxabs() const
    {const T_ARRAY& self=Derived();T result=abs(self(ID(1)));ID m=self.Size();for(ID i(2);i<=m;i++) result=PhysBAM::max(result,abs(self(i)));return result;}

    T Maxmag() const
    {const T_ARRAY& self=Derived();T result=self(ID(1));ID m=self.Size();for(ID i(2);i<=m;i++) result=PhysBAM::maxmag(result,self(i));return result;}

    ID Argmax() const
    {const T_ARRAY& self=Derived();ID result(1),m=self.Size();for(ID i(2);i<=m;i++) if(self(i)>self(result)) result=i;return result;}

    T Min() const
    {const T_ARRAY& self=Derived();T result=self(ID(1));ID m=self.Size();for(ID i(2);i<=m;i++) result=PhysBAM::min(result,self(i));return result;}

    T Minmag() const
    {const T_ARRAY& self=Derived();T result=self(ID(1));ID m=self.Size();for(ID i(2);i<=m;i++) result=PhysBAM::minmag(result,self(i));return result;}

    ID Argmin() const
    {const T_ARRAY& self=Derived();ID result(1),m=self.Size();for(ID i(2);i<=m;i++) if(self(i)<self(result)) result=i;return result;}

    T Sum() const
    {const T_ARRAY& self=Derived();T result=T();ID m=self.Size();for(ID i(1);i<=m;i++) result+=self(i);return result;}

    double Sum_Double_Precision() const
    {const T_ARRAY& self=Derived();double result=0;ID m=self.Size();for(ID i(1);i<=m;i++) result+=self(i);return result;}

    T Sumabs() const
    {const T_ARRAY& self=Derived();T result=T();ID m=self.Size();for(ID i(1);i<=m;i++) result+=abs(self(i));return result;}

    T Componentwise_Maxabs() const
    {const T_ARRAY& self=Derived();T result=abs(self(ID(1)));ID m=self.Size();for(ID i(2);i<=m;i++) result=T::Componentwise_Max(result,abs(self(i)));return result;}

    T Average() const
    {const T_ARRAY& self=Derived();return self.Size()?Sum()/SCALAR(self.Size()):T();}

    template<class T_ARRAY1,class T_ARRAY2>
    static T Linear_Combination(const T_ARRAY1& w,const T_ARRAY2& a)
    {STATIC_ASSERT_SAME(T,typename T_ARRAY2::ELEMENT);assert(w.Size()==a.Size());
    T result=T();ID m=a.Size();for(ID i(1);i<=m;i++) result+=w(i)*a(i);return result;}

    template<class T_ARRAY1,class T_ARRAY2> static typename SCALAR_POLICY<typename T_ARRAY1::ELEMENT>::TYPE
    Dot_Product(const T_ARRAY1& a1,const T_ARRAY2& a2)
    {STATIC_ASSERT_SAME(T,typename T_ARRAY1::ELEMENT);
    STATIC_ASSERT_SAME(T,typename T_ARRAY2::ELEMENT);
    assert(a1.Size()==a2.Size());
    SCALAR result=(SCALAR)0;ID m=a1.Size();for(ID i(1);i<=m;i++) result+=PhysBAM::Dot_Product(a1(i),a2(i));return result;}

    template<class T_ARRAY1,class T_ARRAY2> static double
    Dot_Product_Double_Precision(const T_ARRAY1& a1,const T_ARRAY2& a2)
    {STATIC_ASSERT_SAME(T,typename T_ARRAY1::ELEMENT);
    STATIC_ASSERT_SAME(T,typename T_ARRAY2::ELEMENT);
    assert(a1.Size()==a2.Size());
    double result=0;ID m=a1.Size();for(ID i(1);i<=m;i++) result+=PhysBAM::Dot_Product(a1(i),a2(i));return result;}

    template<class T_ARRAY1,class T_ARRAY2> static typename T_ARRAY1::ELEMENT
    Inner_Product(const T_ARRAY1& m,const T_ARRAY2& a1,const T_ARRAY2& a2,typename boost::enable_if<IS_SCALAR<typename T_ARRAY1::ELEMENT>,UNUSABLE>::type unusable=UNUSABLE())
    {assert(a1.Size()==a2.Size());
    STATIC_ASSERT_SAME(SCALAR,typename T_ARRAY1::ELEMENT);
    STATIC_ASSERT_SAME(T,typename T_ARRAY2::ELEMENT);
    SCALAR result=(SCALAR)0;ID size=a1.Size();for(ID i(1);i<=size;i++) result+=m(i)*PhysBAM::Dot_Product(a1(i),a2(i));return result;}

    template<class T_ARRAY1,class T_ARRAY2> static typename T_ARRAY2::ELEMENT::SCALAR
    Inner_Product(const T_ARRAY1& m,const T_ARRAY2& a1,const T_ARRAY2& a2,typename boost::disable_if<IS_SCALAR<typename T_ARRAY1::ELEMENT>,UNUSABLE>::type unusable=UNUSABLE())
    {STATIC_ASSERT_SAME(T,typename T_ARRAY2::ELEMENT);assert(a1.Size()==a2.Size());
    SCALAR result=(SCALAR)0;ID size=a1.Size();for(ID i(1);i<=size;i++) result+=m(i).Inner_Product(a1(i),a2(i));return result;}

    template<class T_ARRAY1,class T_ARRAY2> static double
    Inner_Product_Double_Precision(const T_ARRAY1& m,const T_ARRAY2& a1,const T_ARRAY2& a2,typename boost::enable_if<IS_SCALAR<typename T_ARRAY1::ELEMENT>,UNUSABLE>::type unusable=UNUSABLE())
    {assert(a1.Size()==a2.Size());
    STATIC_ASSERT_SAME(SCALAR,typename T_ARRAY1::ELEMENT);
    STATIC_ASSERT_SAME(T,typename T_ARRAY2::ELEMENT);
    double result=0;ID size=a1.Size();for(ID i(1);i<=size;i++) result+=m(i)*PhysBAM::Dot_Product(a1(i),a2(i));return result;}

    template<class T_ARRAY1,class T_ARRAY2> static double
    Inner_Product_Double_Precision(const T_ARRAY1& m,const T_ARRAY2& a1,const T_ARRAY2& a2,typename boost::disable_if<IS_SCALAR<typename T_ARRAY1::ELEMENT>,UNUSABLE>::type unusable=UNUSABLE())
    {STATIC_ASSERT_SAME(T,typename T_ARRAY2::ELEMENT);assert(a1.Size()==a2.Size());
    double result=0;ID size=a1.Size();for(ID i(1);i<=size;i++) result+=m(i).Inner_Product(a1(i),a2(i));return result;}

    SCALAR Magnitude_Squared() const
    {const T_ARRAY& self=Derived();
    SCALAR result=(SCALAR)0;ID m=self.Size();for(ID i(1);i<=m;i++) result+=PhysBAM::Magnitude_Squared(self(i));return result;}

    SCALAR Magnitude() const
    {return sqrt(Magnitude_Squared());}

    template<class T_ARRAY1> static typename SCALAR_POLICY<typename T_ARRAY1::ELEMENT>::TYPE
    Distance(const T_ARRAY1& a,const T_ARRAY1& b)
    {STATIC_ASSERT_SAME(T,typename T_ARRAY1::ELEMENT);assert(a.Size()==b.Size());
    SCALAR distance=(SCALAR)0;ID m=a.Size();for(ID i(1);i<=m;i++) distance+=PhysBAM::Magnitude_Squared(a(i)-b(i));return sqrt(distance);}

private:
    SCALAR Maximum_Magnitude_Helper(mpl::true_ is_scalar) const
    {const T_ARRAY& self=Derived();
    T result=(T)0;ID m=self.Size();for(ID i(1);i<=m;i++) result=PhysBAM::max(result,abs(self(i)));return result;}

    SCALAR Maximum_Magnitude_Helper(mpl::false_ is_scalar) const
    {return sqrt(Maximum_Magnitude_Squared());}

    SCALAR Maximum_Magnitude_Squared_Helper(mpl::true_ is_scalar) const
    {return sqr(Maximum_Magnitude());}

    SCALAR Maximum_Magnitude_Squared_Helper(mpl::false_ is_scalar) const
    {const T_ARRAY& self=Derived();
    SCALAR result=(SCALAR)0;ID m=self.Size();for(ID i(1);i<=m;i++) result=PhysBAM::max(result,self(i).Magnitude_Squared());return result;}

    ID Arg_Maximum_Magnitude_Helper(mpl::true_ is_scalar) const
    {const T_ARRAY& self=Derived();ID m=self.Size();
    SCALAR maximum=-1;ID argmax=ID();
    for(ID i(1);i<=m;i++){
        SCALAR current=abs(self(i));
        if(maximum<current){maximum=current;argmax=i;}}
    return argmax;}

    ID Arg_Maximum_Magnitude_Helper(mpl::false_ is_scalar) const
    {const T_ARRAY& self=Derived();ID m=self.Size();
    SCALAR maximum=-1;ID argmax=ID();
    for(ID i(1);i<=m;i++){
        SCALAR current=self(i).Magnitude_Squared();
        if(maximum<current){maximum=current;argmax=i;}}
    return argmax;}
public:

    SCALAR Maximum_Magnitude() const
    {return Maximum_Magnitude_Helper(IS_SCALAR<T>());}

    SCALAR Maximum_Magnitude_Squared() const
    {return Maximum_Magnitude_Squared_Helper(IS_SCALAR<T>());}

    ID Arg_Maximum_Magnitude() const
    {return Arg_Maximum_Magnitude_Helper(IS_SCALAR<T>());}

    template<class T_ARRAY1,class T_ARRAY2> static T
    Maximum_Pairwise_Distance(const T_ARRAY1& a,const T_ARRAY2& b,typename boost::enable_if<IS_SCALAR<typename T_ARRAY1::ELEMENT>,UNUSABLE>::TYPE unusable=UNUSABLE())
    {STATIC_ASSERT_SAME(T,typename T_ARRAY1::ELEMENT);
    STATIC_ASSERT_SAME(T,typename T_ARRAY2::ELEMENT);
    assert(a.Size()==b.Size());
    T result=(T)0;ID m=a.Size();for(ID i(1);i<=m;i++) result=PhysBAM::max(result,abs(a(i)-b(i)));return result;}

    template<class T_ARRAY_TV1,class T_ARRAY_TV2> static typename T_ARRAY_TV1::ELEMENT::SCALAR
    Maximum_Pairwise_Distance(const T_ARRAY_TV1& a,const T_ARRAY_TV2& b)
    {STATIC_ASSERT_SAME(T,typename T_ARRAY_TV1::ELEMENT);
    STATIC_ASSERT_SAME(T,typename T_ARRAY_TV2::ELEMENT);
    assert(a.Size()==b.Size());
    SCALAR result=(SCALAR)0;ID m=a.Size();for(ID i(1);i<=m;i++) result=PhysBAM::max(result,(a(i)-b(i)).Magnitude_Squared());return sqrt(result);}

    static void Find_Common_Elements(const T_ARRAY& a,const T_ARRAY& b,T_ARRAY& result)
    {assert(&a!=&result);assert(&b!=&result);result.Remove_All();
    ID m=a.Size();for(ID i(1);i<=m;i++) if(b.Contains(a(i))) result.Append(a(i));}

    template<class T_ARRAY1,class T_ARRAY2>
    static bool Equal_Dimensions(const T_ARRAY1& a,const T_ARRAY2& b)
    {return a.Size()==b.Size();}

    template<class T_ARRAY1,class T_ARRAY_INT>
    static void Permute(const T_ARRAY1& source,T_ARRAY1& destination,const T_ARRAY_INT& permutation)
    {STATIC_ASSERT_SAME(T,typename T_ARRAY1::ELEMENT);
    STATIC_ASSERT_SAME(ID,typename T_ARRAY_INT::ELEMENT);
    ID m=permutation.Size();for(ID i(1);i<=m;i++) destination(i)=source(permutation(i));}

    template<class T_ARRAY1,class T_ARRAY_INT>
    static void Unpermute(const T_ARRAY1& source,T_ARRAY1& destination,const T_ARRAY_INT& permutation)
    {STATIC_ASSERT_SAME(T,typename T_ARRAY1::ELEMENT);
    STATIC_ASSERT_SAME(ID,typename T_ARRAY_INT::ELEMENT);
    ID m=permutation.Size();for(ID i(1);i<=m;i++) destination(permutation(i))=source(i);}

    template<class T_ARRAY1>
    void Remove_Sorted_Indices(const T_ARRAY1& index)
    {STATIC_ASSERT_SAME(ID,typename T_ARRAY1::ELEMENT);
    T_ARRAY& self=Derived();ID m=self.Size(),index_m=index.Size();
    if(index_m==0) return;
    for(ID kk(1);kk<=index_m-1;kk++){
        assert(1<=index(kk) && index(kk)<=m);
        for(ID i=index(kk)+1-kk;i<=index(kk+1)-1-kk;i++) self(i)=self(i+kk);}
    for(ID i=index(index_m)+1-index_m;i<=m-index_m;i++) self(i)=self(i+index_m);
    self.Resize(m-index_m);}

    template<class T_ARRAY1>
    void Remove_Sorted_Indices_Lazy(const T_ARRAY1& index)
    {STATIC_ASSERT_SAME(ID,typename T_ARRAY1::ELEMENT);
    T_ARRAY& self=Derived();ID index_m=index.Size();
    if(index_m==0) return;
    ID curr=0;
    for(ID k=index_m;k>=ID(1);k--)if(index(k)!=curr){curr=index(k);self.Remove_Index_Lazy(curr);}
    self.Compact();}

    template<class T_ARRAY1>
    static void Reverse_In_Place(T_ARRAY1& input)
    {for(ID i(1);i<=ID(Value(input.m)/2);i++) exchange(input(i),input(input.m+1-i));}

    T* begin() // for stl
    {return Derived().Get_Array_Pointer();}

    const T* begin() const // for stl
    {return Derived().Get_Array_Pointer();}

    T* end() // for stl
    {return Derived().Get_Array_Pointer()+Derived().Size();}

    const T* end() const // for stl
    {return Derived().Get_Array_Pointer()+Derived().Size();}

//#####################################################################
};
template<class T,class T_ARRAY,class ID> inline std::ostream& operator<<(std::ostream& output,const ARRAY_BASE<T,T_ARRAY,ID>& a)
{
    const T_ARRAY& a_=a.Derived();ID m=a_.Size();
    if(m){
        for(ID i(1);i<m;i++) output<<a_(i)<<" ";
        output<<a_(m);}
    output<<std::endl;return output;
}
}
#include <Arrays/IDENTITY_MAP.h>
#include <Arrays/RAW_ARRAY.h>
#endif
