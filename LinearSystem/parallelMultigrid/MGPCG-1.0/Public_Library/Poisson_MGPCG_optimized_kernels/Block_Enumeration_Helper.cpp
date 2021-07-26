//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis, Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
#include "Block_Enumeration_Helper.h"
#include "PTHREAD_QUEUE.h"
#include "Instantiation_Helpers.h"
using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;
//#####################################################################
// Constructor/Desctructor of non-size-specific class
//#####################################################################
CONSTRUCTOR_INSTANTIATION_HELPERX(Block_Enumeration,(const int x_size_input,const int y_size_input,const int z_size_input,const unsigned char* const is_boundary_bytemask_input,
        const unsigned char* const is_extended_boundary_bytemask_input,const int* const boundary_nodes_in_block_input,const int* const extended_boundary_nodes_in_block_input,
	int* const red_boundary_block_start_input,int* const red_boundary_block_end_input,int* const black_boundary_block_start_input,int* const black_boundary_block_end_input,
	int* const red_boundary_indices_input,int* const black_boundary_indices_input,const int red_indices_count_input,int* const extended_boundary_indices_input),(x_size_input,is_boundary_bytemask_input,is_extended_boundary_bytemask_input,boundary_nodes_in_block_input,extended_boundary_nodes_in_block_input,
	     red_boundary_block_start_input, red_boundary_block_end_input, black_boundary_block_start_input, black_boundary_block_end_input,
	    red_boundary_indices_input, black_boundary_indices_input, red_indices_count_input,extended_boundary_indices_input))
//#####################################################################
// Function Run_Parallel
//#####################################################################
namespace{
template<class T,int y_size,int z_size>
struct Block_Enumeration_Size_Specific_Thread_Helper:public PhysBAM::PTHREAD_QUEUE::TASK
{
    Block_Enumeration_Size_Specific_Helper<T,y_size,z_size>* const obj;
    const int xmin,xmax,red_block_offset,black_block_offset,red_index_offset, black_index_offset, extended_index_offset;
    Block_Enumeration_Size_Specific_Thread_Helper(Block_Enumeration_Size_Specific_Helper<T,y_size,z_size>* const obj_input,const int xmin_input,const int xmax_input,const int red_block_offset_input,const int black_block_offset_input,const int red_index_offset_input, const int black_index_offset_input, const int extended_index_offset_input)
        :obj(obj_input),xmin(xmin_input),xmax(xmax_input),red_block_offset(red_block_offset_input),black_block_offset(black_block_offset_input)
	,red_index_offset(red_index_offset_input),black_index_offset(black_index_offset_input),extended_index_offset(extended_index_offset_input){}
    void Run(){obj->Run_X_Range(xmin,xmax,red_block_offset,black_block_offset,red_index_offset,black_index_offset,extended_index_offset);}
};
}
template<class T,int y_size,int z_size> void Block_Enumeration_Size_Specific_Helper<T,y_size,z_size>::
Run_Parallel(const int number_of_partitions)
{
    if(x_size*y_size*z_size<=32*32*32 || number_of_partitions ==1) {Run();return;} // Run serially for small problems
    
    int red_block_offset=0,black_block_offset=0,red_index_offset=0, black_index_offset=0, extended_index_offset=0;
    for(int partition=0;partition<number_of_partitions;partition++){
        int first_block_of_partition=(number_of_x_blocks/number_of_partitions)*(partition)+std::min(number_of_x_blocks%number_of_partitions,partition)+1;
        int last_block_of_partition=(number_of_x_blocks/number_of_partitions)*(partition+1)+std::min(number_of_x_blocks%number_of_partitions,partition+1);
        int xmin=(first_block_of_partition-1)*x_block_size+1;
        int xmax=last_block_of_partition*x_block_size;
        Block_Enumeration_Size_Specific_Thread_Helper<T,y_size,z_size>* task=new Block_Enumeration_Size_Specific_Thread_Helper<T,y_size,z_size>(this,xmin,xmax,red_block_offset,black_block_offset,red_index_offset,black_index_offset,extended_index_offset);
        pthread_queue->Queue(task);
	Compute_Offsets(xmin,xmax,red_block_offset,black_block_offset,red_index_offset, black_index_offset,extended_index_offset);

    }
    pthread_queue->Wait();
}
//#####################################################################
// Function Run_X_Range
//#####################################################################
template<class T,int y_size,int z_size> void Block_Enumeration_Size_Specific_Helper<T,y_size,z_size>::
Run_X_Range(const int xmin,const int xmax,const int red_block_offset,const int black_block_offset,const int red_index_offset, const int black_index_offset, const int extended_index_offset)
{    

    int* local_red_block_start = &red_boundary_block_start[red_block_offset];
    int* local_red_block_end=&red_boundary_block_end[red_block_offset];
    int* local_black_block_start = &black_boundary_block_start[black_block_offset];
    int* local_black_block_end=&black_boundary_block_end[black_block_offset];
    int* local_extended_indices = &extended_boundary_indices[extended_index_offset];
    
    int red_index=red_index_offset;
    int black_index=black_index_offset;

    for(int block_i=xmin;block_i<=xmax;block_i+=x_block_size)
    for(int block_j=1;block_j<=y_size;block_j+=y_block_size)
    for(int block_k=1;block_k<=z_size;block_k+=z_block_size){
 
        int bi=(block_i-1)/x_block_size;
        int bj=(block_j-1)/y_block_size;
        int bk=(block_k-1)/z_block_size;
        
        int block_index=(bi*number_of_y_blocks+bj)*number_of_z_blocks+bk;
	int block_color=(bi+bj+bk)%2;
	int* block_indices = 0;
	if(boundary_nodes_in_block[block_index]){
	    if(block_color){
		*local_black_block_start=black_index+red_indices_count;
		local_black_block_start++;
		*local_black_block_end=black_index+red_indices_count+boundary_nodes_in_block[block_index];
		local_black_block_end++;
		block_indices=&black_boundary_indices[black_index];
		black_index+=boundary_nodes_in_block[block_index];
	    }else{
		*local_red_block_start=red_index;
		local_red_block_start++;
		*local_red_block_end=red_index+boundary_nodes_in_block[block_index];
		local_red_block_end++;
		block_indices=&red_boundary_indices[red_index];	
		red_index+=boundary_nodes_in_block[block_index];
	    }
	}
	    
	int c=0,e=0;
        for(int i=block_i;i<block_i+x_block_size;i++)
        for(int j=block_j;j<block_j+y_block_size;j++)
        for(int k=block_k;k<block_k+z_block_size;k++){

            int index=i*x_shift+j*y_shift+k*z_shift;	    

            if(is_boundary_bytemask[index]){ *block_indices = index; block_indices++;c++;}
            if(is_extended_boundary_bytemask[index]){ *local_extended_indices=index; local_extended_indices++;e++;}
        }
	if(c!=boundary_nodes_in_block[block_index]){ std::cout<<"c("<<bi<<","<<bj<<","<<bk<<")="<<c<<":"<<boundary_nodes_in_block[block_index]<<std::endl;}
	if(e!=extended_boundary_nodes_in_block[block_index]){std::cout<<"e("<<bi<<","<<bj<<","<<bk<<")="<<e<<":"<<extended_boundary_nodes_in_block[block_index]<<std::endl;}
    }
}
//#####################################################################
// Function Compute_Offsets
//#####################################################################
template<class T,int y_size,int z_size> void Block_Enumeration_Size_Specific_Helper<T,y_size,z_size>::
Compute_Offsets(int& xmin,int& xmax,int& red_block_offset,int& black_block_offset,int& red_index_offset,int& black_index_offset,int& extended_index_offset)
{
    for(int block_i=xmin;block_i<=xmax;block_i+=x_block_size)
    for(int block_j=1;block_j<=y_size;block_j+=y_block_size)
    for(int block_k=1;block_k<=z_size;block_k+=z_block_size){
 
        int bi=(block_i-1)/x_block_size;
        int bj=(block_j-1)/y_block_size;
        int bk=(block_k-1)/z_block_size;
        
        int block_index=(bi*number_of_y_blocks+bj)*number_of_z_blocks+bk;
	int block_color=(bi+bj+bk)%2;
	if(block_color){
	    black_block_offset+=boundary_nodes_in_block[block_index]?1:0;
	    black_index_offset+=boundary_nodes_in_block[block_index];
	}else{
	    red_block_offset+=boundary_nodes_in_block[block_index]?1:0;
	    red_index_offset+=boundary_nodes_in_block[block_index];
	}
	extended_index_offset+=extended_boundary_nodes_in_block[block_index];
    }
		
}

//#####################################################################
template class Block_Enumeration_Helper<void>;

