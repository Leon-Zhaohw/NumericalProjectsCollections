//#####################################################################
// Copyright 2002-2007, Geoffrey Irving, Frank Losasso, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VECTOR_1D
//#####################################################################
#ifndef __VECTOR_1D__
#define __VECTOR_1D__

#include <Matrices_And_Vectors/SCALAR_POLICY.h>
#include <Matrices_And_Vectors/VECTOR_0D.h>
#include <Matrices_And_Vectors/VECTOR_BASE.h>
#include <Math_Tools/Inverse.h>
#include <Math_Tools/clamp.h>
#include <Math_Tools/argmax.h>
#include <Math_Tools/argmin.h>
#include <Math_Tools/max.h>
#include <Math_Tools/min.h>
#include <Math_Tools/sign.h>
#include <Math_Tools/sqr.h>
#include <Read_Write/READ_WRITE_FUNCTIONS.h>
#include <Utilities/STATIC_ASSERT.h>
#include <Utilities/DEBUG_UTILITIES.h>
#include <cmath>
namespace PhysBAM{

using ::std::floor;
using ::std::sin;
using ::std::cos;

template<class T>
class VECTOR<T,1>:public VECTOR_BASE<T,VECTOR<T,1> >
{
    struct UNUSABLE{};
public:
    template<class T2> struct REBIND{typedef VECTOR<T2,1> TYPE;};
    typedef typename IF<IS_SCALAR<T>::value,T,UNUSABLE>::TYPE SCALAR;
    typedef T ELEMENT;
    typedef int INDEX;
    typedef T value_type; // for stl
    typedef T* iterator; // for stl
    enum WORKAROUND1 {dimension=1};
    enum WORKAROUND2 {m=1};

    T x;

    explicit VECTOR(INITIAL_SIZE n=INITIAL_SIZE(1))
        :x()
    {
        STATIC_ASSERT(sizeof(VECTOR)==sizeof(T));assert(n==INITIAL_SIZE(1));
    }

    explicit VECTOR(const T& x_input)
        :x(x_input)
    {}

    template<class T2> explicit VECTOR(const VECTOR<T2,1>& vector_input)
        :x((T)vector_input.x)
    {}

    VECTOR(const VECTOR& vector_input)
        :x(vector_input.x)
    {}

    template<class T_VECTOR>
    explicit VECTOR(const T_VECTOR& v,typename boost::enable_if<AND<IS_SAME<T,typename T_VECTOR::ELEMENT>::value,INTS_EQUAL<T_VECTOR::m,1>::value>,UNUSABLE>::type unusable=UNUSABLE())
        :x(v(1))
    {}

    VECTOR(const VECTOR_ND<T>& vector_input)
        :x(vector_input(1))
    {
        assert(vector_input.n==1);
    }

    template<int n>
    VECTOR(const VECTOR<T,n>& v1,const VECTOR<T,1-n>& v2)
    {
        for(int i=1;i<=n;i++) (*this)(i)=v1(i);for(int i=n+1;i<=1;i++) (*this)(i)=v2(i-n);
    }

    template<class T_VECTOR> typename boost::enable_if<boost::mpl::and_<IS_SAME<T,typename T_VECTOR::ELEMENT>,INTS_EQUAL<T_VECTOR::m,1> >,VECTOR&>::type
    operator=(const T_VECTOR& v)
    {
        x=v(1);return *this;
    }

    VECTOR& operator=(const VECTOR& v)
    {
        x=v(1);return *this;
    }

    int Size() const
    {return 1;}

    const T& operator[](const int i) const
    {assert(i==1);return x;}

    T& operator[](const int i)
    {assert(i==1);return x;}

    const T& operator()(const int i) const
    {assert(i==1);return x;}

    T& operator()(const int i)
    {assert(i==1);return x;}

    bool operator==(const VECTOR& v) const
    {return x==v.x;}

    bool operator!=(const VECTOR& v) const
    {return x!=v.x;}

    VECTOR operator-() const
    {return VECTOR(-x);}

    VECTOR& operator+=(const VECTOR& v)
    {x+=v.x;return *this;}

    VECTOR& operator-=(const VECTOR& v)
    {x-=v.x;return *this;}

    VECTOR& operator*=(const VECTOR& v)
    {x*=v.x;return *this;}

    VECTOR& operator+=(const T& a)
    {x+=a;return *this;}

    VECTOR& operator-=(const T& a)
    {x-=a;return *this;}

    VECTOR& operator*=(const T& a)
    {x*=a;return *this;}

    VECTOR& operator*=(const INT_INVERSE a)
    {x*=a;return *this;}

    VECTOR& operator/=(const T& a)
    {return *this*=Inverse(a);}

    VECTOR& operator/=(const VECTOR& v)
    {x/=v.x;return *this;}

    VECTOR operator+(const VECTOR& v) const
    {return VECTOR(x+v.x);}

    VECTOR operator-(const VECTOR& v) const
    {return VECTOR(x-v.x);}

    VECTOR operator*(const VECTOR& v) const
    {return VECTOR(x*v.x);}

    VECTOR operator/(const VECTOR& v) const
    {return VECTOR(x/v.x);}

    VECTOR operator+(const T& a) const
    {return VECTOR(x+a);}

    VECTOR operator-(const T& a) const
    {return VECTOR(x-a);}

    VECTOR operator*(const T& a) const
    {return VECTOR(x*a);}

    VECTOR operator*(const INT_INVERSE a) const
    {return VECTOR(x*a);}

    VECTOR operator/(const T& a) const
    {return *this*Inverse(a);}

    T Magnitude_Squared() const
    {return sqr(x);}

    T Magnitude() const
    {return abs(x);}

    T L1_Norm() const
    {return abs(x);}

    T Normalize()
    {T magnitude=abs(x);x=(T)(x>=0?1:-1);return magnitude;}

    VECTOR Normalized() const
    {return VECTOR((T)(x>=0?1:-1));}

    T Min() const
    {return x;}

    T Max() const
    {return x;}

    T Max_Abs() const
    {return abs(x);}

    int Arg_Min() const
    {return 1;}

    int Arg_Max() const
    {return 1;}

    bool Elements_Equal() const
    {return true;}

    int Dominant_Axis() const
    {return 1;}

    static T Dot_Product(const VECTOR& v1,const VECTOR& v2)
    {return v1.x*v2.x;}

    static VECTOR Componentwise_Min(const VECTOR& v1,const VECTOR& v2)
    {return VECTOR(min(v1.x,v2.x));}

    static VECTOR Componentwise_Max(const VECTOR& v1,const VECTOR& v2)
    {return VECTOR(max(v1.x,v2.x));}

    bool All_Greater(const VECTOR& v) const
    {return x>v.x;}

    bool All_Less(const VECTOR& v) const
    {return x<v.x;}

    bool All_Greater_Equal(const VECTOR& v) const
    {return x>=v.x;}

    bool All_Less_Equal(const VECTOR& v) const
    {return x<=v.x;}

    VECTOR Projected_On_Unit_Direction(const VECTOR& direction) const
    {return Dot_Product(*this,direction)*direction;}

    VECTOR Projected(const VECTOR& direction) const // un-normalized direction
    {return Dot_Product(*this,direction)/direction.Magnitude_Squared()*direction;}

    void Project_On_Unit_Direction(const VECTOR& direction)
    {*this=Dot_Product(*this,direction)*direction;}

    void Project(const VECTOR& direction) // un-normalized direction
    {*this=Dot_Product(*this,direction)/direction.Magnitude_Squared()*direction;}

    VECTOR Projected_Orthogonal_To_Unit_Direction(const VECTOR& direction) const
    {return *this-Dot_Product(*this,direction)*direction;}

    static VECTOR<T,1> Cross_Product(const VECTOR<T,0>,const VECTOR&)
    {return VECTOR<T,1>();}

    static VECTOR<T,1> Cross_Product(const VECTOR&,const VECTOR<T,0>)
    {return VECTOR<T,1>();}

    static VECTOR<T,0> Cross_Product(const VECTOR&,const VECTOR&)
    {return VECTOR<T,0>();}

    T Sum() const
    {return x;}

    T Average() const
    {return x;}

    T Product() const
    {return x;}

    const VECTOR& Column_Sum() const
    {return *this;}

    int Number_True() const
    {STATIC_ASSERT((IS_SAME<T,bool>::value));return x;}

    static VECTOR Axis_Vector(const int axis)
    {assert(axis==1);return VECTOR((T)1);}

    static VECTOR All_Ones_Vector()
    {return VECTOR((T)1);}

    void Fill(const T& constant)
    {x=constant;}

    void Get(T& element1) const
    {element1=x;}

    void Set(const T& element1)
    {x=element1;}

    template<class T_FUNCTION>
    static VECTOR Map(const T_FUNCTION& f,const VECTOR& v)
    {return VECTOR(f(v.x));}

    int Find(const T& element) const
    {return x==element?1:0;}

    bool Contains(const T& element) const
    {return x==element;}

    template<class T_ARRAY>
    bool Contains_All(const T_ARRAY& elements) const
    {STATIC_ASSERT((IS_SAME<typename T_ARRAY::ELEMENT,T>::value));
    for(int i=1;i<=elements.Size();i++) if(!Contains(elements(i))) return false;
    return true;}

    template<class T_ARRAY>
    bool Contains_Any(const T_ARRAY& elements) const
    {STATIC_ASSERT((IS_SAME<typename T_ARRAY::ELEMENT,T>::value));
    for(int i=1;i<=elements.Size();i++) if(Contains(elements(i))) return true;
    return false;}

    VECTOR<T,0> Remove_Index(const int index) const
    {assert(1==index);return VECTOR<T,0>();}

    VECTOR<T,2> Insert(const T& element,const int index) const
    {VECTOR<T,2> r;r[index]=element;r[3-index]=x;return r;}

    VECTOR<T,2> Append(const T& element) const
    {return VECTOR<T,2>(x,element);}

    template<int d2> VECTOR<T,1+d2> Append_Elements(const VECTOR<T,d2>& elements) const
    {VECTOR<T,1+d2> r;r[1]=x;for(int i=1;i<=d2;i++) r[i+1]=elements[i];return r;}

    VECTOR Reversed() const
    {return *this;}

    template<int d1,int d2> VECTOR<int,d2-d1+1> Slice() const
    {STATIC_ASSERT((mpl::and_<mpl::less_equal<mpl::int_<1>,mpl::int_<d1> >,mpl::less_equal<mpl::int_<d2>,mpl::int_<1> > >::value));
    VECTOR<T,d2-d1+1> r;for(int i=d1;i<=d2;i++) r[i-d1+1]=(*this)[i];return r;}

    template<int n> void Split(VECTOR<T,n>& v1,VECTOR<T,1-n>& v2) const
    {for(int i=1;i<=n;i++) v1(i)=(*this)(i);
    for(int i=n+1;i<=1;i++) v2(i-n)=(*this)(i);}

    template<class RW>
    void Read(std::istream& input_stream)
    {Read_Binary<RW>(input_stream,x);}

    template<class RW>
    void Write(std::ostream& output_stream) const
    {Write_Binary<RW>(output_stream,x);}

    T* begin() // for stl
    {return &x;}

    const T* begin() const // for stl
    {return &x;}

    T* end() // for stl
    {return &x+1;}

    const T* end() const // for stl
    {return &x+1;}

//#####################################################################
};

//#####################################################################
// Miscellaneous free operators and functions
//#####################################################################
template<class T> inline VECTOR<T,1>
operator+(const T& a,const VECTOR<T,1>& v)
{return VECTOR<T,1>(a+v.x);}

template<class T> inline VECTOR<T,1>
operator-(const T& a,const VECTOR<T,1>& v)
{return VECTOR<T,1>(a-v.x);}

template<class T> inline VECTOR<T,1>
operator*(const T& a,const VECTOR<T,1>& v)
{return VECTOR<T,1>(a*v.x);}

template<class T> inline VECTOR<T,1>
operator/(const T& a,const VECTOR<T,1>& v)
{return VECTOR<T,1>(a/v.x);}

template<class T> inline VECTOR<T,1>
abs(const VECTOR<T,1>& v)
{return VECTOR<T,1>(abs(v.x));}

template<class T> inline VECTOR<T,1>
floor(const VECTOR<T,1>& v)
{return VECTOR<T,1>(floor(v.x));}

template<class T> inline VECTOR<T,1>
exp(const VECTOR<T,1>& v)
{return VECTOR<T,1>(exp(v.x));}

template<class T> inline VECTOR<T,1>
sin(const VECTOR<T,1>& v)
{return VECTOR<T,1>(sin(v.x));}

template<class T> inline VECTOR<T,1>
cos(const VECTOR<T,1>& v)
{return VECTOR<T,1>(cos(v.x));}

template<class T> inline VECTOR<T,1>
sqrt(const VECTOR<T,1>& v)
{return VECTOR<T,1>(sqrt(v.x));}

template<class T> inline VECTOR<T,1>
Inverse(const VECTOR<T,1>& v)
{return VECTOR<T,1>(1/v.x);}

//#####################################################################
// Functions clamp, clamp_min, clamp_max, in_bounds
//#####################################################################
template<class T> inline VECTOR<T,1>
clamp(const VECTOR<T,1>& v,const VECTOR<T,1>& vmin,const VECTOR<T,1>& vmax)
{return VECTOR<T,1>(clamp(v.x,vmin.x,vmax.x));}

template<class T> inline VECTOR<T,1>
clamp(const VECTOR<T,1>& v,T min,T max)
{return VECTOR<T,1>(clamp(v.x,min,max));}

template<class T> inline VECTOR<T,1>
clamp_min(const VECTOR<T,1>& v,const VECTOR<T,1>& vmin)
{return VECTOR<T,1>(clamp_min(v.x,vmin.x));}

template<class T> inline VECTOR<T,1>
clamp_min(const VECTOR<T,1>& v,const T& min)
{return VECTOR<T,1>(clamp_min(v.x,min));}

template<class T> inline VECTOR<T,1>
clamp_max(const VECTOR<T,1>& v,const VECTOR<T,1>& vmax)
{return VECTOR<T,1>(clamp_max(v.x,vmax.x));}

template<class T> inline VECTOR<T,1>
clamp_max(const VECTOR<T,1>& v,const T& max)
{return VECTOR<T,1>(clamp_max(v.x,max));}

template<class T> inline bool
in_bounds(const VECTOR<T,1>& v,const VECTOR<T,1>& vmin,const VECTOR<T,1>& vmax)
{return in_bounds(v.x,vmin.x,vmax.x);}

//#####################################################################
// Stream input and output
//#####################################################################
template<class T> inline std::istream&
operator>>(std::istream& input,VECTOR<T,1>& v)
{return input>>v.x;}

template<class T> inline std::ostream&
operator<<(std::ostream& output,const VECTOR<T,1>& v)
{return output<<v.x;}

//#####################################################################
}
#endif
