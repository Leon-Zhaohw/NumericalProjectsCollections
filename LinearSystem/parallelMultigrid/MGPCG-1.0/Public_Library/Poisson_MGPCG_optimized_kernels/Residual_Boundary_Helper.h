//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis, Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
#ifndef __Residual_Boundary_Helper__
#define __Residual_Boundary_Helper__
namespace PhysBAM{

template<class T>
class Residual_Boundary_Size_Specific_Helper_Base
{
protected:
    const T* const u;
    const T* const b;
    T* const r;
    const T* const diagonal_part;
    const int* const boundary_index;
    const int number_of_indices;

public:
    explicit Residual_Boundary_Size_Specific_Helper_Base(const T* const u_input,const T* const b_input,T* const r_input,const T* const diagonal_part_input,
	const int* const boundary_index_input,const int number_of_indices_input)
        :u(u_input),b(b_input),r(r_input),diagonal_part(diagonal_part_input), boundary_index(boundary_index_input),number_of_indices(number_of_indices_input)
    {}

    virtual ~Residual_Boundary_Size_Specific_Helper_Base() {}

//#####################################################################
    virtual void Run()=0;
    virtual void Run_Parallel(const int number_of_partitions)=0;
//#####################################################################
};

template<class T>
class Residual_Boundary_Helper{
    const int x_size,y_size,z_size;
    Residual_Boundary_Size_Specific_Helper_Base<T>* derived;

public:
    ~Residual_Boundary_Helper() {delete derived;}
    void Run() {derived->Run();}
    void Run_Parallel(const int number_of_partitions) {derived->Run_Parallel(number_of_partitions);}

//#####################################################################
    Residual_Boundary_Helper(const int x_size_input,const int y_size_input,const int z_size_input,const T* const u,const T* const b,T* const r,
        const T* const diagonal_part,const int* const boundary_index, const int number_of_indices);
//#####################################################################
};

template<class T,int y_size,int z_size=y_size>
class Residual_Boundary_Size_Specific_Helper:public Residual_Boundary_Size_Specific_Helper_Base<T>
{
    typedef Residual_Boundary_Size_Specific_Helper_Base<T> Base;
    using Base::u;using Base::b;using Base::r;using Base::diagonal_part;using Base::boundary_index;using Base::number_of_indices;

    const int x_size,padded_x_size;

    enum WORKAROUND{
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
    explicit Residual_Boundary_Size_Specific_Helper(const int x_size_input,const T* const u_input,const T* const b_input,T* const r_input,const T* const diagonal_part_input,
	const int* const boundary_index_input,const int number_of_indices_input)
        :Base(u_input,b_input,r_input,diagonal_part_input,boundary_index_input,number_of_indices_input)
	,x_size(x_size_input),padded_x_size(x_size_input+2)
    {}

    void Run()
    {Run_Index_Range(0,number_of_indices);}

//#####################################################################
    void Run_Parallel(const int number_of_partitions);
    void Run_Index_Range(const int index_start, const int index_end);
//#####################################################################
};
}
#endif
