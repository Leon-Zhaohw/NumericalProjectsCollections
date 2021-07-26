//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis,Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
#ifndef __Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema_Helper__
#define __Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema_Helper__
namespace PhysBAM{

template<class T>
class Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema_Size_Specific_Helper_Base
{
protected:
    T* const result;
    const T* const u1;
    const T& c1;
    double &sum;
    T& minimum;
    T& maximum;
    double* sum_partial_results;
    T* minimum_partial_results;
    T* maximum_partial_results;

public:
    explicit Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema_Size_Specific_Helper_Base(const T& c1_input,const T* const u1_input,
        T* const result_input,double &sum_input,T& minimum_input,T& maximum_input)
	:result(result_input),u1(u1_input),c1(c1_input),sum(sum_input),minimum(minimum_input),maximum(maximum_input),sum_partial_results(0),
        minimum_partial_results(0),maximum_partial_results(0)
    {}

    virtual ~Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema_Size_Specific_Helper_Base() {}

//#####################################################################
    virtual void Run()=0;
    virtual void Run_Parallel(const int number_of_partitions)=0;
//#####################################################################
};

template<class T>
class Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema_Helper{
    const int x_size,y_size,z_size;
    Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema_Size_Specific_Helper_Base<T>* derived;

public:
    ~Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema_Helper() {delete derived;}
    void Run() {derived->Run();}
    void Run_Parallel(const int number_of_partitions) {derived->Run_Parallel(number_of_partitions);}

//#####################################################################
    Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema_Helper(const int x_size_input,const int y_size_input,const int z_size_input,const T& c1,
        const T* const u1,T* const result,double &sum_input,T& minimum_input,T& maximum_input);
//#####################################################################
};

template<class T,int y_size,int z_size=y_size>
class Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema_Size_Specific_Helper
    :public Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema_Size_Specific_Helper_Base<T>
{
    typedef Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema_Size_Specific_Helper_Base<T> Base;
    using Base::result;using Base::u1;using Base::c1;using Base::sum;using Base::minimum;using Base::maximum;using Base::sum_partial_results;
    using Base::minimum_partial_results;using Base::maximum_partial_results;

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
    explicit Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema_Size_Specific_Helper(const int x_size_input,const T& c1_input,const T* const u1_input,T* const result_input,
        double &sum_input,T& minimum_input,T& maximum_input)
        :Base(c1_input,u1_input,result_input,sum_input,minimum_input,maximum_input)
	,x_size(x_size_input),padded_x_size(x_size_input+2)
    {}

    void Run()
    {sum_partial_results=new double[1];minimum_partial_results=new T[1];maximum_partial_results=new T[1];
    Run_X_Range(1,x_size,0);
    sum=sum_partial_results[0];minimum=minimum_partial_results[0];maximum=maximum_partial_results[0];
    delete[] sum_partial_results;sum_partial_results=0;
    delete[] minimum_partial_results;minimum_partial_results=0;
    delete[] maximum_partial_results;maximum_partial_results=0;}
  

//#####################################################################
    void Run_Parallel(const int number_of_partitions);
    void Run_X_Range(const int xmin,const int xmax,const int partition_number);
//#####################################################################
};
}
#endif
