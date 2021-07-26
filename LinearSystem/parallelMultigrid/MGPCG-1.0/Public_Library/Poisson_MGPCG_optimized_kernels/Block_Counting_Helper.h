//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis,Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
#ifndef __Block_Counting_Helper__
#define __Block_Counting_Helper__

namespace PhysBAM{

template<class T>
class Block_Counting_Size_Specific_Helper_Base
{
protected:
    const unsigned char* const is_boundary_bytemask;
    const unsigned char* const is_extended_boundary_bytemask;
    int* const boundary_nodes_in_block;
    int* const extended_boundary_nodes_in_block;
    int& total_red_boundary_blocks;
    int& total_black_boundary_blocks;
    int& total_extended_boundary_blocks;
    int& total_red_boundary_indices;
    int& total_black_boundary_indices;
    int& total_extended_boundary_indices;
    int* red_boundary_blocks_partial_results;
    int* black_boundary_blocks_partial_results;
    int* extended_boundary_blocks_partial_results;
    int* red_boundary_indices_partial_results;
    int* black_boundary_indices_partial_results;
    int* extended_boundary_indices_partial_results;

public:
    explicit Block_Counting_Size_Specific_Helper_Base(const unsigned char* const is_boundary_bytemask_input,
        const unsigned char* const is_extended_boundary_bytemask_input,int* const boundary_nodes_in_block_input,int* const extended_boundary_nodes_in_block_input,
	int& total_red_boundary_blocks_input,int& total_black_boundary_blocks_input,int& total_extended_boundary_blocks_input,
	int& total_red_boundary_indices_input,int& total_black_boundary_indices_input,int& total_extended_boundary_indices_input)
        :is_boundary_bytemask(is_boundary_bytemask_input),is_extended_boundary_bytemask(is_extended_boundary_bytemask_input),
        boundary_nodes_in_block(boundary_nodes_in_block_input),extended_boundary_nodes_in_block(extended_boundary_nodes_in_block_input),
	total_red_boundary_blocks(total_red_boundary_blocks_input),total_black_boundary_blocks(total_black_boundary_blocks_input),
	total_extended_boundary_blocks(total_extended_boundary_blocks_input),
	total_red_boundary_indices(total_red_boundary_indices_input),total_black_boundary_indices(total_black_boundary_indices_input),
	total_extended_boundary_indices(total_extended_boundary_indices_input)
    {}

    virtual ~Block_Counting_Size_Specific_Helper_Base() {}

//#####################################################################
    virtual void Run()=0;
    virtual void Run_Parallel(const int number_of_partitions)=0;
//#####################################################################
};

template<class T>
class Block_Counting_Helper{
    const int x_size,y_size,z_size;
    Block_Counting_Size_Specific_Helper_Base<T>* derived;

public:
    ~Block_Counting_Helper() {delete derived;}
    void Run() {derived->Run();}
    void Run_Parallel(const int number_of_partitions) {derived->Run_Parallel(number_of_partitions);}

//#####################################################################
    Block_Counting_Helper(const int x_size_input,const int y_size_input,const int z_size_input,const unsigned char* const is_boundary_bytemask_input,
        const unsigned char* const is_extended_boundary_bytemask_input,int* const boundary_nodes_in_block_input,int* const extended_boundary_nodes_in_block_input,
	int& total_red_boundary_blocks_input,int& total_black_boundary_blocks_input,int& total_extended_boundary_blocks_input,
	int& total_red_boundary_indices_input,int& total_black_boundary_indices_input,int& total_extended_boundary_indices_input);
//#####################################################################
};

template<class T,int y_size,int z_size=y_size>
class Block_Counting_Size_Specific_Helper:public Block_Counting_Size_Specific_Helper_Base<T>
{
    typedef Block_Counting_Size_Specific_Helper_Base<T> Base;
    using Base::is_boundary_bytemask;using Base::is_extended_boundary_bytemask;using Base::boundary_nodes_in_block;using Base::extended_boundary_nodes_in_block;

    using Base::total_red_boundary_blocks;
    using Base::total_black_boundary_blocks;
    using Base::total_extended_boundary_blocks;
    using Base::total_red_boundary_indices;
    using Base::total_black_boundary_indices;
    using Base::total_extended_boundary_indices;
    using Base::red_boundary_blocks_partial_results;
    using Base::black_boundary_blocks_partial_results;
    using Base::extended_boundary_blocks_partial_results;
    using Base::red_boundary_indices_partial_results;
    using Base::black_boundary_indices_partial_results;
    using Base::extended_boundary_indices_partial_results;

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
    explicit Block_Counting_Size_Specific_Helper(const int x_size_input,const unsigned char* const is_boundary_bytemask_input,
        const unsigned char* const is_extended_boundary_bytemask_input,int* const boundary_nodes_in_block_input,int* const extended_boundary_nodes_in_block_input,
	int& total_red_boundary_blocks_input,int& total_black_boundary_blocks_input,int& total_extended_boundary_blocks_input,
	int& total_red_boundary_indices_input,int& total_black_boundary_indices_input,int& total_extended_boundary_indices_input)
        :Base(is_boundary_bytemask_input,is_extended_boundary_bytemask_input,boundary_nodes_in_block_input,extended_boundary_nodes_in_block_input,
	    total_red_boundary_blocks_input,total_black_boundary_blocks_input,total_extended_boundary_blocks_input,
	    total_red_boundary_indices_input,total_black_boundary_indices_input,total_extended_boundary_indices_input)
	,x_size(x_size_input),padded_x_size(x_size_input+2),number_of_x_blocks(x_size_input/x_block_size)
    {}

    void Run()
    {
	
	red_boundary_blocks_partial_results = new int[1];
	black_boundary_blocks_partial_results = new int[1];
	extended_boundary_blocks_partial_results = new int[1];
	red_boundary_indices_partial_results = new int[1];
	black_boundary_indices_partial_results = new int[1];
	extended_boundary_indices_partial_results = new int[1];
	Run_X_Range(1,x_size,0);

	total_red_boundary_blocks = red_boundary_blocks_partial_results[0];
	total_black_boundary_blocks = black_boundary_blocks_partial_results[0];
	total_extended_boundary_blocks=extended_boundary_blocks_partial_results[0];
	total_red_boundary_indices=red_boundary_indices_partial_results[0];
	total_black_boundary_indices=black_boundary_indices_partial_results[0];
	total_extended_boundary_indices=extended_boundary_indices_partial_results[0];
	
	delete[] red_boundary_blocks_partial_results; red_boundary_blocks_partial_results=0;
	delete[] black_boundary_blocks_partial_results; black_boundary_blocks_partial_results=0;
	delete[] extended_boundary_blocks_partial_results; extended_boundary_blocks_partial_results=0;
	delete[] red_boundary_indices_partial_results; red_boundary_indices_partial_results=0;
	delete[] black_boundary_indices_partial_results; black_boundary_indices_partial_results=0;
	delete[] extended_boundary_indices_partial_results; extended_boundary_indices_partial_results=0;
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
    void Run_X_Range(const int xmin,const int xmax,const int partition);
//#####################################################################
};
}
#endif
