//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis, Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
#ifndef __Block_Enumeration_Helper__
#define __Block_Enumeration_Helper__

namespace PhysBAM{

template<class T>
class Block_Enumeration_Size_Specific_Helper_Base
{
protected:
    const unsigned char* const is_boundary_bytemask;
    const unsigned char* const is_extended_boundary_bytemask;
    const int* const boundary_nodes_in_block;
    const int* const extended_boundary_nodes_in_block;
    int* const red_boundary_block_start;
    int* const red_boundary_block_end;
    int* const black_boundary_block_start;
    int* const black_boundary_block_end;
    int* const red_boundary_indices;
    int* const black_boundary_indices;
    const int red_indices_count;
    int* const extended_boundary_indices;
public:
    explicit Block_Enumeration_Size_Specific_Helper_Base(const unsigned char* const is_boundary_bytemask_input,
        const unsigned char* const is_extended_boundary_bytemask_input,const int* const boundary_nodes_in_block_input,const int* const extended_boundary_nodes_in_block_input,
	int* const red_boundary_block_start_input,int* const red_boundary_block_end_input,int* const black_boundary_block_start_input,int* const black_boundary_block_end_input,
	int* const red_boundary_indices_input,int* const black_boundary_indices_input,const int red_indices_count_input,int* const extended_boundary_indices_input)
        :is_boundary_bytemask(is_boundary_bytemask_input),is_extended_boundary_bytemask(is_extended_boundary_bytemask_input),
        boundary_nodes_in_block(boundary_nodes_in_block_input),extended_boundary_nodes_in_block(extended_boundary_nodes_in_block_input),
	red_boundary_block_start(red_boundary_block_start_input),red_boundary_block_end(red_boundary_block_end_input),
	black_boundary_block_start(black_boundary_block_start_input),black_boundary_block_end(black_boundary_block_end_input),
	red_boundary_indices(red_boundary_indices_input),black_boundary_indices(black_boundary_indices_input),
	red_indices_count(red_indices_count_input),extended_boundary_indices(extended_boundary_indices_input)
    {}

    virtual ~Block_Enumeration_Size_Specific_Helper_Base() {}

//#####################################################################
    virtual void Run()=0;
    virtual void Run_Parallel(const int number_of_partitions)=0;
//#####################################################################
};

template<class T>
class Block_Enumeration_Helper{
    const int x_size,y_size,z_size;
    Block_Enumeration_Size_Specific_Helper_Base<T>* derived;

public:
    ~Block_Enumeration_Helper() {delete derived;}
    void Run() {derived->Run();}
    void Run_Parallel(const int number_of_partitions) {derived->Run_Parallel(number_of_partitions);}

//#####################################################################
    Block_Enumeration_Helper(const int x_size_input,const int y_size_input,const int z_size_input,const unsigned char* const is_boundary_bytemask_input,
        const unsigned char* const is_extended_boundary_bytemask_input,const int* const boundary_nodes_in_block_input,const int* const extended_boundary_nodes_in_block_input,
	int* const red_boundary_block_start_input,int* const red_boundary_block_end_input,int* const black_boundary_block_start_input,int* const black_boundary_block_end_input,
	int* const red_boundary_indices_input,int* const black_boundary_indices_input,const int red_indices_count_input,int* const extended_boundary_indices_input);
//#####################################################################
};

template<class T,int y_size,int z_size=y_size>
class Block_Enumeration_Size_Specific_Helper:public Block_Enumeration_Size_Specific_Helper_Base<T>
{
    typedef Block_Enumeration_Size_Specific_Helper_Base<T> Base;
    using Base::is_boundary_bytemask;using Base::is_extended_boundary_bytemask;using Base::boundary_nodes_in_block;using Base::extended_boundary_nodes_in_block;
    using Base::red_boundary_block_start;
    using Base::red_boundary_block_end;
    using Base::black_boundary_block_start;
    using Base::black_boundary_block_end;
    using Base::red_boundary_indices;
    using Base::black_boundary_indices;
    using Base::red_indices_count;
    using Base::extended_boundary_indices;


    const int x_size;
    const int padded_x_size;
    const int number_of_x_blocks;
    enum WORKAROUND{
        INTERIOR_CELL_TYPE=1,
        DIRICHLET_CELL_TYPE=2,
        NEUMANN_CELL_TYPE=3,
        x_block_size=4,
        y_block_size=4,
        z_block_size=4,
        number_of_y_blocks=y_size/y_block_size,
        number_of_z_blocks=z_size/z_block_size,

        padded_y_size=y_size+2,
        padded_z_size=z_size+2,
	
	x_shift=padded_y_size*padded_z_size,
        y_shift=padded_z_size,
        z_shift=1,
    };


public:
    explicit Block_Enumeration_Size_Specific_Helper(const int x_size_input,const unsigned char* const is_boundary_bytemask_input,
        const unsigned char* const is_extended_boundary_bytemask_input,const int* const boundary_nodes_in_block_input,const int* const extended_boundary_nodes_in_block_input,
	int* const red_boundary_block_start_input,int* const red_boundary_block_end_input,int* const black_boundary_block_start_input,int* const black_boundary_block_end_input,
	int* const red_boundary_indices_input,int* const black_boundary_indices_input,const int red_indices_count_input,int* const extended_boundary_indices_input)
        :Base(is_boundary_bytemask_input,is_extended_boundary_bytemask_input,boundary_nodes_in_block_input,extended_boundary_nodes_in_block_input,
	     red_boundary_block_start_input, red_boundary_block_end_input, black_boundary_block_start_input, black_boundary_block_end_input,
	    red_boundary_indices_input, black_boundary_indices_input, red_indices_count_input,extended_boundary_indices_input)
	,x_size(x_size_input),padded_x_size(x_size_input+2),number_of_x_blocks(x_size_input/x_block_size)
    {}

    void Run()
    {
	Run_X_Range(1,x_size,0,0,0,0,0);
    }
    
    // For debugging purposes only

//     static void Allocate_Data(T*& u,T*& u_coarse,unsigned char*& bit_writemask)
//     {
// 	int padded_length = padded_y_size*padded_y_size*padded_z_size;
// 	int coarse_padded_length=coarse_padded_y_size*coarse_padded_y_size*coarse_padded_z_size;u=new T[padded_length];
// 	u_coarse=new T[coarse_padded_length];bit_writemask=new unsigned char[coarse_padded_length];}

//     static void Initialize_Data(T* const u,T* const u_coarse,unsigned char* const bit_writemask)
//     {
// 	int padded_length = padded_y_size*padded_y_size*padded_z_size;
// 	int coarse_padded_length=coarse_padded_y_size*coarse_padded_y_size*coarse_padded_z_size;
// 	for(int i=0;i<padded_length;i++) u[i]=(T)i;
// 	for(int i=0;i<coarse_padded_length;i++) u_coarse[i]=(T)i;
// 	for(int i=0;i<coarse_padded_length;i++) bit_writemask[i]=i%256;}

//#####################################################################
    void Run_Parallel(const int number_of_partitions);
    void Run_X_Range(const int xmin,const int xmax,const int red_block_offset,const int black_block_offset,const int red_index_offset, const int black_index_offset, const int extended_index_offset);

    void Compute_Offsets(int& xmin,int& xmax,int& red_block_offset,int& black_block_offset,int& red_index_offset,int& black_index_offset,int& extended_index_offset);

//#####################################################################
};
}
#endif
