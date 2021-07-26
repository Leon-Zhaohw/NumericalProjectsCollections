//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis, Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
#ifndef __Multiply_And_Compute_Sum_And_Extrema_Helper__
#define __Multiply_And_Compute_Sum_And_Extrema_Helper__
namespace PhysBAM{

template<class T>
class Multiply_And_Compute_Sum_And_Extrema_Size_Specific_Helper_Base
{
protected:
    const T* const u;
    T* const v;
    const T* const diagonal_part;
    const T scale_factor;
    double& sum;
    T& minimum;
    T& maximum;
    double* sum_partial_results;
    T* minimum_partial_results;
    T* maximum_partial_results;

public:
    explicit Multiply_And_Compute_Sum_And_Extrema_Size_Specific_Helper_Base(const T* const u_input,T* const v_input,const T* const diagonal_part_input,
        double &sum_input,T& minimum_input,T& maximum_input,const T scale_factor_input)
	:u(u_input),v(v_input),diagonal_part(diagonal_part_input),scale_factor(scale_factor_input),sum(sum_input),minimum(minimum_input),maximum(maximum_input),
         sum_partial_results(0),minimum_partial_results(0),maximum_partial_results(0)
    {}
    
    virtual ~Multiply_And_Compute_Sum_And_Extrema_Size_Specific_Helper_Base() {}

//#####################################################################
    virtual void Run()=0;
    virtual void Run_Parallel(const int number_of_partitions)=0;
//#####################################################################
};

template<class T>
class Multiply_And_Compute_Sum_And_Extrema_Helper{
    const int x_size,y_size,z_size;
    Multiply_And_Compute_Sum_And_Extrema_Size_Specific_Helper_Base<T>* derived;

public:
    ~Multiply_And_Compute_Sum_And_Extrema_Helper() {delete derived;}
    void Run() {derived->Run();}
    void Run_Parallel(const int number_of_partitions) {derived->Run_Parallel(number_of_partitions);}

//#####################################################################
    Multiply_And_Compute_Sum_And_Extrema_Helper(const int x_size_input,const int y_size_input,const int z_size_input,const T* const u,T* const v,
        const T* const diagonal_part,double &sum_input,T& minimum_input,T& maximum_input,const T scale_factor=(T)1);
//#####################################################################
};

template<class T,int y_size,int z_size=y_size>
class Multiply_And_Compute_Sum_And_Extrema_Size_Specific_Helper:public Multiply_And_Compute_Sum_And_Extrema_Size_Specific_Helper_Base<T>
{
    typedef Multiply_And_Compute_Sum_And_Extrema_Size_Specific_Helper_Base<T> Base;    
    using Base::v;using Base::u;using Base::diagonal_part;using Base::scale_factor;using Base::sum;using Base::minimum;using Base::maximum;
    using Base::sum_partial_results;using Base::minimum_partial_results;using Base::maximum_partial_results;

    const int x_size,padded_x_size;

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
    explicit Multiply_And_Compute_Sum_And_Extrema_Size_Specific_Helper(const int x_size_input,const T* const u_input,T* const v_input,const T* const diagonal_part_input,
        double &sum_input,T& minimum_input,T& maximum_input,const T scale_factor_input)
	:Base(u_input,v_input,diagonal_part_input,sum_input,minimum_input,maximum_input,scale_factor_input)
	,x_size(x_size_input),padded_x_size(x_size_input+2)
    {}
    
    void Run()
    {sum_partial_results=new double[1];minimum_partial_results=new T[1];maximum_partial_results=new T[1];
    Run_X_Range(1,x_size,0);
    sum=sum_partial_results[0];minimum=minimum_partial_results[0];maximum=maximum_partial_results[0];
    delete[] sum_partial_results;sum_partial_results=0;
    delete[] minimum_partial_results;minimum_partial_results=0;
    delete[] maximum_partial_results;maximum_partial_results=0;}
    
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
    void Run_X_Range(const int xmin,const int xmax,const int partition_number);
//#####################################################################
};
}
#endif
