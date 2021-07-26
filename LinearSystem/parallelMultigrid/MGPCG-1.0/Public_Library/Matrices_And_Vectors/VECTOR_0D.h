//#####################################################################
// Copyright 2007, Nipun Kwatra, Craig Schroeder, Tamar Shinar, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VECTOR_0D
//#####################################################################
#ifndef __VECTOR_0D__
#define __VECTOR_0D__

#include <Matrices_And_Vectors/SCALAR_POLICY.h>
#include <Matrices_And_Vectors/VECTOR_BASE.h>
#include <Utilities/DEBUG_UTILITIES.h>
#include <cmath>
namespace PhysBAM{

template<class T>
class VECTOR<T,0>:public VECTOR_BASE<T,VECTOR<T,0> >
{
    struct UNUSABLE{};
public:
    enum WORKAROUND1 {dimension=0};
    enum WORKAROUND2 {m=0};
    typedef typename IF<IS_SCALAR<T>::value,T,UNUSABLE>::TYPE SCALAR;
    template<class T2> struct REBIND{typedef VECTOR<T2,0> TYPE;};
    typedef T ELEMENT;
    typedef int INDEX;
    typedef VECTOR_BASE<T,VECTOR<T,0> > BASE;

    explicit VECTOR(INITIAL_SIZE n=INITIAL_SIZE())
    {
        assert(n==INITIAL_SIZE());
    }

    VECTOR(const VECTOR& vector_input)
        :BASE()
    {}

    template<class T2> explicit VECTOR(const VECTOR<T2,0>& vector_input)
    {}

    VECTOR& operator=(const VECTOR& v)
    {return *this;}

    const T& operator[](const int) const
    {PHYSBAM_FATAL_ERROR();}

    T& operator[](const int)
    {PHYSBAM_FATAL_ERROR();}

    const T& operator()(const int i) const
    {PHYSBAM_FATAL_ERROR();}

    T& operator()(const int i)
    {PHYSBAM_FATAL_ERROR();}

    bool operator==(const VECTOR& v) const
    {return true;}

    bool operator!=(const VECTOR& v) const
    {return false;}

    VECTOR operator-() const
    {return *this;}

    VECTOR& operator+=(const VECTOR&)
    {return *this;}

    VECTOR& operator-=(const VECTOR&)
    {return *this;}

    VECTOR& operator*=(const VECTOR&)
    {return *this;}

    VECTOR& operator*=(const T&)
    {return *this;}

    VECTOR& operator/=(const T&)
    {return *this;}

    VECTOR& operator/=(const VECTOR&)
    {return *this;}

    VECTOR operator+(const VECTOR&) const
    {return *this;}

    VECTOR operator+(const T&) const
    {return *this;}

    VECTOR operator-(const VECTOR&) const
    {return *this;}

    VECTOR operator*(const VECTOR&) const
    {return *this;}

    VECTOR operator/(const VECTOR&) const
    {return *this;}

    VECTOR operator*(const T&) const
    {return *this;}

    VECTOR operator/(const T&) const
    {return *this;}

    bool Contains(const T&) const
    {return false;}

    T Magnitude_Squared() const
    {return 0;}

    T Magnitude() const
    {return 0;}

    T Normalize()
    {return T();}

    VECTOR Normalized() const
    {return *this;}

    static T Dot_Product(const VECTOR&,const VECTOR&)
    {return T();}

    VECTOR<T,1> Insert(const T& element,const int index) const
    {VECTOR<T,1> r;r[index]=element;return r;}

    static VECTOR All_Ones_Vector()
    {return VECTOR();}

    template<class RW>
    void Read(std::istream&)
    {}

    template<class RW>
    void Write(std::ostream&) const
    {}
};

template<class T> inline VECTOR<T,0>
operator*(const T&,const VECTOR<T,0>& v)
{return v;}

template<class T> inline VECTOR<T,0>
operator/(const T&,const VECTOR<T,0>& v)
{return v;}

template<class T> inline std::istream&
operator>>(std::istream& input,VECTOR<T,0>&)
{return input;}

template<class T> inline std::ostream&
operator<<(std::ostream& output,const VECTOR<T,0>&)
{return output;}

}
#endif
