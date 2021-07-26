//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis, Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
#ifndef __Relaxation_Boundary_Helper__
#define __Relaxation_Boundary_Helper__
namespace PhysBAM{

template<class T>
class Relaxation_Boundary_Size_Specific_Helper_Base
{
protected:
    T* const u;
    const T* const b;
    const T* const one_over_diagonal_part;
    const int* const boundary_index;
    const int* const block_start;
    const int* const block_end;
    const int number_of_red_blocks;
    const int number_of_black_blocks;
    const int loops;
    const bool reverse_order;

public:
    explicit Relaxation_Boundary_Size_Specific_Helper_Base(T* const u_input,const T* const b_input,const T* const one_over_diagonal_part_input,
        const int* const boundary_index_input,const int* const block_start_input,const int* const block_end_input,const int number_of_red_blocks_input,
        const int number_of_black_blocks_input,const int loops_input,const bool reverse_order_input)
        :u(u_input),b(b_input),one_over_diagonal_part(one_over_diagonal_part_input),boundary_index(boundary_index_input),block_start(block_start_input),
         block_end(block_end_input),number_of_red_blocks(number_of_red_blocks_input),number_of_black_blocks(number_of_black_blocks_input),loops(loops_input),
         reverse_order(reverse_order_input)
    {}

    virtual ~Relaxation_Boundary_Size_Specific_Helper_Base() {}

//#####################################################################
    virtual void Run()=0;
    virtual void Run_Parallel(const int number_of_partitions)=0;
//#####################################################################
};

template<class T>
class Relaxation_Boundary_Helper{
    const int x_size,y_size,z_size;
    Relaxation_Boundary_Size_Specific_Helper_Base<T>* derived;

public:
    ~Relaxation_Boundary_Helper() {delete derived;}
    void Run() {derived->Run();}
    void Run_Parallel(const int number_of_partitions) {derived->Run_Parallel(number_of_partitions);}

//#####################################################################
    Relaxation_Boundary_Helper(const int x_size_input,const int y_size_input,const int z_size_input,T* const u,const T* const b,
        const T* const one_over_diagonal_part,const int* const boundary_index,const int* const block_start,const int* const block_end,
        const int number_of_red_blocks,const int number_of_black_blocks,const int loops,const bool reverse_order);
//#####################################################################
};

template<class T,int y_size,int z_size=y_size>
class Relaxation_Boundary_Size_Specific_Helper:public Relaxation_Boundary_Size_Specific_Helper_Base<T>
{
    typedef Relaxation_Boundary_Size_Specific_Helper_Base<T> Base;
    using Base::u;using Base::b;using Base::one_over_diagonal_part;using Base::boundary_index;using Base::block_start;using Base::block_end;
    using Base::number_of_red_blocks;using Base::number_of_black_blocks;using Base::loops;using Base::reverse_order;

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
    explicit Relaxation_Boundary_Size_Specific_Helper(const int x_size_input,T* const u_input,const T* const b_input,const T* const one_over_diagonal_part_input,
	const int* const boundary_index_input,const int* const block_start_input,const int* const block_end_input,const int number_of_red_blocks_input,
        const int number_of_black_blocks_input,const int loops_input,const bool reverse_order_input)
        :Base(u_input,b_input,one_over_diagonal_part_input,boundary_index_input,block_start_input,block_end_input,number_of_red_blocks_input,
            number_of_black_blocks_input,loops_input,reverse_order_input)
	,x_size(x_size_input),padded_x_size(x_size_input+2)
    {}

    void Run()
    {
	Run_Block_Range(0,number_of_red_blocks+number_of_black_blocks-1);
    }

    // For debugging purposes only

//     static void Allocate_Data(T*& u,T*& b,T*& r,T*& one_over_diagonal_part)
//     {u=new T[padded_length];b=new T[padded_length];r=new T[padded_length];one_over_diagonal_part=new T[padded_length];}

//     static void Initialize_Data(T* const u,T* const b,T* const r,T* const one_over_diagonal_part)
//     {for(int i=0;i<padded_length;i++) u[i]=b[i]=r[i]=one_over_diagonal_part[i]=(T)i;}

//#####################################################################
    void Run_Parallel(const int number_of_partitions);
    void Run_Block_Range(const int first_block,const int last_block);
//#####################################################################
};
}
#endif
