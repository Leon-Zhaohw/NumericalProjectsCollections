//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis,Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
#ifndef __Combined_Saxpy_Helper__
#define __Combined_Saxpy_Helper__
namespace PhysBAM{

template<class T>
class Combined_Saxpy_Size_Specific_Helper_Base
{
protected:
    // x := x + alpha*p
    // p := z + beta*p
    T* const x;
    T* const p;
    const T* const z;
    const T alpha;
    const T beta;

public:
    explicit Combined_Saxpy_Size_Specific_Helper_Base(T* const x_input,T* const p_input,const T* const z_input,const T alpha_input,const T beta_input)
        :x(x_input),p(p_input),z(z_input),alpha(alpha_input),beta(beta_input)
    {}

    virtual ~Combined_Saxpy_Size_Specific_Helper_Base() {};

//#####################################################################
    virtual void Run()=0;
    virtual void Run_Parallel(const int number_of_partitions)=0;
//#####################################################################
};

template<class T>
class Combined_Saxpy_Helper{
    const int x_size,y_size,z_size;
    Combined_Saxpy_Size_Specific_Helper_Base<T>* derived;

public:
    ~Combined_Saxpy_Helper() {delete derived;}
    void Run() {derived->Run();}
    void Run_Parallel(const int number_of_partitions) {derived->Run_Parallel(number_of_partitions);}

//#####################################################################
    Combined_Saxpy_Helper(const int x_size_input,const int y_size_input,const int z_size_input,T* const x_input,T* const p_input,const T* const z_input,const T alpha_input,const T beta_input);
//#####################################################################
};

template<class T,int y_size,int z_size=y_size>
class Combined_Saxpy_Size_Specific_Helper:public Combined_Saxpy_Size_Specific_Helper_Base<T>
{
    typedef Combined_Saxpy_Size_Specific_Helper_Base<T> Base;
    using Base::x;using Base::p;using Base::z;using Base::alpha;using Base::beta;

    const int x_size,padded_x_size;

    enum WORKAROUND{
        x_block_size=4,
        y_block_size=4,
        z_block_size=4,
        padded_y_size=y_size+2,
        padded_z_size=z_size+2,
        x_shift=padded_y_size*padded_z_size,
        y_shift=padded_z_size,
        z_shift=1
    };

public:
    explicit Combined_Saxpy_Size_Specific_Helper(const int x_size_input,T* const x_input,T* const p_input,const T* const z_input,const T alpha_input,const T beta_input)
        :Base(x_input,p_input,z_input,alpha_input,beta_input)
	,x_size(x_size_input),padded_x_size(x_size_input+2)
    {}

    void Run()
    {Run_X_Range(1,x_size);}
  
    // For debugging purposes only

    static void Allocate_Data(T*& x,T*& y,T*& diagonal_part)
    {
        int padded_length=padded_y_size*padded_y_size*padded_z_size;
	x=new T[padded_length];y=new T[padded_length];diagonal_part=new T[padded_length];}

    static void Initialize_Data(T* const x,T* const y,T* const diagonal_part)
    {
        int padded_length=padded_y_size*padded_y_size*padded_z_size;
	for(int i=0;i<padded_length;i++) x[i]=y[i]=diagonal_part[i]=(T)i;}

//#####################################################################
    void Run_Parallel(const int number_of_partitions);
    void Run_X_Range(const int xmin, const int xmax);
//#####################################################################
};
}
#endif
