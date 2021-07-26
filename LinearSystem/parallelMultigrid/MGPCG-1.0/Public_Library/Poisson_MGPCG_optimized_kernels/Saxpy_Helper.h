//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis,Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
#ifndef __Saxpy_Helper__
#define __Saxpy_Helper__
namespace PhysBAM{

template<class T>
class Saxpy_Size_Specific_Helper_Base
{
protected:
    // result(i) += c1*u1(i)
    T* const result;
    const T* const u1;
    const T& c1;

public:
    explicit Saxpy_Size_Specific_Helper_Base(const T& c1_input,const T* const u1_input,T* const result_input)
	:result(result_input),u1(u1_input),c1(c1_input)
    {}

    virtual ~Saxpy_Size_Specific_Helper_Base() {};

//#####################################################################
    virtual void Run()=0;
    virtual void Run_Parallel(const int number_of_partitions)=0;
//#####################################################################
};

template<class T>
class Saxpy_Helper{
    const int x_size,y_size,z_size;
    Saxpy_Size_Specific_Helper_Base<T>* derived;

public:
    ~Saxpy_Helper() {delete derived;}
    void Run() {derived->Run();}
    void Run_Parallel(const int number_of_partitions) {derived->Run_Parallel(number_of_partitions);}

//#####################################################################
    Saxpy_Helper(const int x_size_input,const int y_size_input,const int z_size_input,const T& c1,const T* const u1,T* const result);
//#####################################################################
};

template<class T,int y_size,int z_size=y_size>
class Saxpy_Size_Specific_Helper:public Saxpy_Size_Specific_Helper_Base<T>
{
    typedef Saxpy_Size_Specific_Helper_Base<T> Base;
    using Base::result;using Base::u1;using Base::c1;

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
    explicit Saxpy_Size_Specific_Helper(const int x_size_input,const T& c1_input,const T* const u1_input,T* const result_input)
	:Base(c1_input,u1_input,result_input)
	,x_size(x_size_input),padded_x_size(x_size_input+2)
    {}

    void Run()
    {Run_X_Range(1,x_size);}
  
    // For debugging purposes only

    static void Allocate_Data(T*& result,T*& u1)
    {
	int padded_length=padded_y_size*padded_y_size*padded_z_size;
        result=new T[padded_length];
        u1=new T[padded_length];
    }

    static void Initialize_Data(T* const result,T* const u1)
    {
	int padded_length=padded_y_size*padded_y_size*padded_z_size;
        for(int i=0;i<padded_length;i++) result[i]=u1[i]=(T)i;
    }

//#####################################################################
    void Run_Parallel(const int number_of_partitions);
    void Run_X_Range(const int xmin,const int xmax);
//#####################################################################
};
}
#endif
