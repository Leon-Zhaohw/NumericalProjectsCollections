//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Neil Molino, Igor Neverov, Duc Nguyen, Craig Schroeder, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RANDOM_NUMBERS
//#####################################################################
#ifndef __RANDOM_NUMBERS__
#define __RANDOM_NUMBERS__

#include <Geometry/BOX.h>
#include <Geometry/GEOMETRY_POLICY.h>
#include <Utilities/NONCOPYABLE.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>
#include <ctime>
#include <boost/cstdint.hpp>
namespace PhysBAM{

using ::std::log;

template<class T,class GENERATOR=boost::mt19937>
class RANDOM_NUMBERS:public NONCOPYABLE
{
public:
    typedef int HAS_TYPED_READ_WRITE;

private:
    int gaussian_iset; // Used to force Get_Gaussian to reset
    T gset; // used internally by Get_Gaussian
    GENERATOR random_number_generator;
    boost::uniform_real<T> random_real;
    boost::variate_generator<GENERATOR&,boost::uniform_real<T> > variate;

public:
    explicit RANDOM_NUMBERS(const unsigned int seed=time(0))
        :gaussian_iset(0),random_real((T)0,(T)1),variate(random_number_generator,random_real)
    {
        Set_Seed(seed);
    }

    virtual ~RANDOM_NUMBERS()
    {}

    void Set_Seed(const unsigned int seed_input=time(0))
    {random_number_generator.seed((boost::uint32_t)seed_input);}

    int Get_Uniform_Integer(const int a,const int b)
    {return min(b,(int)(a+(b+1-a)*Get_Number()));} // in [a,b]

    T Get_Uniform_Number(const T a,const T b)
    {STATIC_ASSERT((NOT<std::numeric_limits<T>::is_integer>::value));
    return a+(b-a)*Get_Number();} // in [a,b)

    template<int d>
    VECTOR<T,d> Get_Uniform_Vector(const VECTOR<T,d>& v0,const VECTOR<T,d>& v1)
    {VECTOR<T,d> r;for(int i=1;i<=d;i++) r(i)=Get_Uniform_Number(v0(i),v1(i));return r;}

    template<int d>
    VECTOR<T,d> Get_Uniform_Vector(const T a,const T b)
    {VECTOR<T,d> r;Fill_Uniform_Vector(r,a,b);return r;}

    template<class T_VECTOR>
    void Fill_Uniform_Vector(VECTOR_BASE<T,T_VECTOR>& v,const T a,const T b)
    {for(int i=1;i<=v.Size();i++) v(i)=Get_Uniform_Number(a,b);}

    template<class T_MATRIX>
    void Fill_Uniform_Matrix(MATRIX_BASE<T,T_MATRIX>& m,const T a,const T b)
    {for(int i=1;i<=m.Rows();i++) for(int j=1;j<=m.Columns();j++) m(i,j)=Get_Uniform_Number(a,b);}

    template<int d>
    void Fill_Uniform_Matrix(DIAGONAL_MATRIX<T,d>& m,const T a,const T b)
    {for(int i=1;i<=m.Rows();i++) m(i,i)=Get_Uniform_Number(a,b);}

    template<int d>
    void Fill_Uniform_Matrix(SYMMETRIC_MATRIX<T,d>& m,const T a,const T b)
    {for(int i=1;i<=m.Rows();i++) for(int j=1;j<=i;j++) m(i,j)=Get_Uniform_Number(a,b);}

    template<int d>
    void Fill_Uniform_Matrix(UPPER_TRIANGULAR_MATRIX<T,d>& m,const T a,const T b)
    {for(int j=1;j<=d;j++) for(int i=1;i<=j;i++) m(i,j)=Get_Uniform_Number(a,b);}

    template<class TV>
    TV Get_Uniform_Vector(const BOX<TV>& box)
    {return Get_Uniform_Vector(box.min_corner,box.max_corner);}

    T Get_Gaussian()
    {T fac,rsq,v1,v2;
    if(gaussian_iset==0){
        do{v1=2*Get_Uniform_Number((T)0,(T)1)-1;v2=2*Get_Uniform_Number((T)0,(T)1)-1;rsq=sqr(v1)+sqr(v2);}while(rsq>=1 || rsq==0);
        fac=sqrt(-2*log(rsq)/rsq);gset=v1*fac;gaussian_iset=1;return v2*fac;}
    else{gaussian_iset=0;return gset;}}

    template<class TV>
    TV Get_Vector_In_Unit_Sphere()
    {for(;;){
        TV v=Get_Uniform_Vector(BOX<TV>());
        if(v.Magnitude_Squared()<=1) return v;}}

    template<class TV>
    TV Get_Direction()
    {if(!TV::m) return TV();
    for(;;){
        TV v=Get_Uniform_Vector(BOX<TV>());
        typename TV::SCALAR magnitude_squared=v.Magnitude_Squared();
        if(magnitude_squared>0 && magnitude_squared<=1) return v/sqrt(magnitude_squared);}}

    template<class TV>
    ROTATION<TV> Get_Rotation()
    {return Get_Rotation_Helper(Get_Direction<VECTOR<T,2*TV::m-2> >());}

    template<class TV>
    FRAME<TV> Get_Frame(const TV& v0,const TV& v1)
    {TV v=Get_Uniform_Vector(v0,v1);return FRAME<TV>(v,Get_Rotation<TV>());}

    template<class TV>
    TWIST<TV> Get_Twist(const T& a)
    {TWIST<TV> tw;tw.Set_Vector(Get_Uniform_Vector<T,TWIST<TV>::m>(-a,a));return tw;}

private:
    static ROTATION<VECTOR<T,1> > Get_Rotation_Helper(const VECTOR<T,0>&)
    {return ROTATION<VECTOR<T,1> >();}

    static ROTATION<VECTOR<T,2> > Get_Rotation_Helper(const VECTOR<T,2>& v)
    {return ROTATION<VECTOR<T,2> >::From_Complex(COMPLEX<T>(v));}

    static ROTATION<VECTOR<T,3> > Get_Rotation_Helper(const VECTOR<T,4>& v)
    {return ROTATION<VECTOR<T,3> >::From_Quaternion(QUATERNION<T>(v));}
public:

    virtual void Read(TYPED_ISTREAM& input)
    {Read_Binary(input,gaussian_iset,gset);
#if !defined(WIN32) || !defined(_MSC_VER) || _MSC_VER > 1400
    input.stream>>random_number_generator>>random_real; // does not compile on VS2005
#else
    PHYSBAM_NOT_IMPLEMENTED();
#endif
    }

    virtual void Write(TYPED_OSTREAM& output) const
    {Write_Binary(output,gaussian_iset,gset);
#if !defined(WIN32) || !defined(_MSC_VER) || _MSC_VER > 1400
    output.stream<<random_number_generator<<random_real; // does not compile on VS2005
#else
    PHYSBAM_NOT_IMPLEMENTED();
#endif
    }

    T Get_Number()
    {return variate();} // in [0,1)
//#####################################################################
};
}
#endif
