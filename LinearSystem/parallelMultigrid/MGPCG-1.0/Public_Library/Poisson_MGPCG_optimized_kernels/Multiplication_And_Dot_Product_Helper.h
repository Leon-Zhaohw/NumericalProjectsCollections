//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis, Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
#ifndef __Multiplication_And_Dot_Product_Helper__
#define __Multiplication_And_Dot_Product_Helper__
namespace PhysBAM{

template<class T>
class Multiplication_And_Dot_Product_Size_Specific_Helper_Base
{
protected:
    //Compute v=Lu
    const T* const u;
    T* const v;
    const T* const diagonal_part;
    double* u_dot_v_partial_results;
    const T scale_factor;

public:
    explicit Multiplication_And_Dot_Product_Size_Specific_Helper_Base(const T* const u_input,T* const v_input,const T* const diagonal_part_input,const T scale_factor_input)
	:u(u_input),v(v_input),diagonal_part(diagonal_part_input),u_dot_v_partial_results(0),scale_factor(scale_factor_input)
    {}
    
    virtual ~Multiplication_And_Dot_Product_Size_Specific_Helper_Base() {}

//#####################################################################
    virtual double Run()=0;
    virtual double Run_Parallel(const int number_of_partitions)=0;
//#####################################################################
};

template<class T>
class Multiplication_And_Dot_Product_Helper{
    const int x_size,y_size,z_size;
    Multiplication_And_Dot_Product_Size_Specific_Helper_Base<T>* derived;

public:
    ~Multiplication_And_Dot_Product_Helper() {delete derived;}
    double Run() {return derived->Run();}
    double Run_Parallel(const int number_of_partitions) {return derived->Run_Parallel(number_of_partitions);}

//#####################################################################
    Multiplication_And_Dot_Product_Helper(const int x_size_input,const int y_size_input,const int z_size_input,const T* const u,T* const v,const T* const diagonal_part,const T scale_factor=(T)1);
//#####################################################################
};

template<class T,int y_size,int z_size=y_size>
class Multiplication_And_Dot_Product_Size_Specific_Helper:public Multiplication_And_Dot_Product_Size_Specific_Helper_Base<T>
{
    typedef Multiplication_And_Dot_Product_Size_Specific_Helper_Base<T> Base;
    using Base::u;using Base::v;using Base::diagonal_part;using Base::u_dot_v_partial_results;using Base::scale_factor;

    int x_size,padded_x_size;

    enum WORKAROUND{
        x_block_size=4,
        y_block_size=4,
        z_block_size=4,
        padded_y_size=y_size+2,
        padded_z_size=z_size+2,
        x_shift=padded_y_size*padded_z_size,
        y_shift=padded_z_size,
        z_shift=1,
        x_plus_one_shift=x_shift,
        x_minus_one_shift=-x_shift,
        y_plus_one_shift=y_shift,
        y_minus_one_shift=-y_shift,
        z_plus_one_shift=z_shift,
        z_minus_one_shift=-z_shift
    };

public:
    explicit Multiplication_And_Dot_Product_Size_Specific_Helper(const int x_size_input,const T* const u_input,T* const v_input,const T* const diagonal_part_input,const T scale_factor_input=(T)1)
        :Base(u_input,v_input,diagonal_part_input,scale_factor_input)
	,x_size(x_size_input),padded_x_size(x_size_input+2)
    {}
    
    double Run()
    {u_dot_v_partial_results=new double[1];
    Run_X_Range(1,x_size,0);
    double u_dot_v=u_dot_v_partial_results[0];
    delete[] u_dot_v_partial_results;u_dot_v_partial_results=0;
    return u_dot_v;}
    
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
    double Run_Parallel(const int number_of_partitions);
    void Run_X_Range(const int xmin,const int xmax,const int partition_number);
//#####################################################################
};
}
#endif
