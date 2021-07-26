//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis,Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
#include "Restriction_Helper.h"
#include "PTHREAD_QUEUE.h"
#include "Instantiation_Helpers.h"
using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;
//#####################################################################
// Constructor of non-size-specific class
//#####################################################################
CONSTRUCTOR_INSTANTIATION_HELPERX(Restriction,(const int x_size_input,const int y_size_input,const int z_size_input,const T* const r_input,T* const b_coarse_input,const unsigned char* const cell_type_coarse_input),(x_size_input,r_input,b_coarse_input,cell_type_coarse_input))
//#####################################################################
// Function Run_Parallel
//#####################################################################
namespace{
template<class T,int y_size,int z_size>
struct Restriction_Size_Specific_Thread_Helper:public PhysBAM::PTHREAD_QUEUE::TASK
{
    Restriction_Size_Specific_Helper<T,y_size,z_size>* const obj;
    const int xmin,xmax;
    Restriction_Size_Specific_Thread_Helper(Restriction_Size_Specific_Helper<T,y_size,z_size>* const obj_input,const int xmin_input,const int xmax_input)
        :obj(obj_input),xmin(xmin_input),xmax(xmax_input) {}
    void Run(){obj->Run_X_Range(xmin,xmax);}
};
}
template<class T,int y_size,int z_size> void Restriction_Size_Specific_Helper<T,y_size,z_size>::
Run_Parallel(const int number_of_partitions)
{
    if(x_size*y_size*z_size<=32*32*32) {Run();return;} // Run serially for small problems

    int number_of_x_blocks=coarse_x_size/x_block_size;
    for(int partition=1;partition<=number_of_partitions;partition++){
        int first_block_of_partition=(number_of_x_blocks/number_of_partitions)*(partition-1)+std::min(number_of_x_blocks%number_of_partitions,partition-1)+1;
        int last_block_of_partition=(number_of_x_blocks/number_of_partitions)*partition+std::min(number_of_x_blocks%number_of_partitions,partition);
        int xmin=(first_block_of_partition-1)*x_block_size+1;
        int xmax=last_block_of_partition*x_block_size;
        Restriction_Size_Specific_Thread_Helper<T,y_size,z_size>* task=new Restriction_Size_Specific_Thread_Helper<T,y_size,z_size>(this,xmin,xmax);
        pthread_queue->Queue(task);}
    pthread_queue->Wait();
}
//#####################################################################
// Function Run_X_Range
//#####################################################################
template<class T,int y_size,int z_size> void Restriction_Size_Specific_Helper<T,y_size,z_size>::
Run_X_Range(const int xmin,const int xmax)
{    
    for(int block_i=xmin;block_i<=xmax;block_i+=x_block_size)
    for(int block_j=1;block_j<=coarse_y_size;block_j+=y_block_size)
    for(int block_k=1;block_k<=coarse_z_size;block_k+=z_block_size)
 
        for(int i=block_i;i<block_i+x_block_size;i++)
        for(int j=block_j;j<block_j+y_block_size;j++)
        for(int k=block_k;k<block_k+z_block_size;k++)
        {    
            int coarse_index=i*coarse_x_shift+j*coarse_y_shift+k*coarse_z_shift;
            int index=2*(i*x_shift+j*y_shift+k*z_shift)-x_shift-y_shift-z_shift;            

            b_coarse[coarse_index]=
                r[index+x_minus_one_y_minus_one_z_minus_one_shift]
                +r[index+x_minus_one_y_minus_one_z_plus_two_shift]
                +r[index+x_minus_one_y_plus_two_z_minus_one_shift]
                +r[index+x_minus_one_y_plus_two_z_plus_two_shift]
                +r[index+x_plus_two_y_minus_one_z_minus_one_shift]
                +r[index+x_plus_two_y_minus_one_z_plus_two_shift]
                +r[index+x_plus_two_y_plus_two_z_minus_one_shift]
                +r[index+x_plus_two_y_plus_two_z_plus_two_shift];
            b_coarse[coarse_index]+=(T)3*(
                r[index+x_minus_one_y_minus_one_shift]
                +r[index+x_minus_one_y_minus_one_z_plus_one_shift]
                +r[index+x_minus_one_z_minus_one_shift]
                +r[index+x_minus_one_z_plus_two_shift]
                +r[index+x_minus_one_y_plus_one_z_minus_one_shift]
                +r[index+x_minus_one_y_plus_one_z_plus_two_shift]
                +r[index+x_minus_one_y_plus_two_shift]
                +r[index+x_minus_one_y_plus_two_z_plus_one_shift]
                +r[index+y_minus_one_z_minus_one_shift]
                +r[index+y_minus_one_z_plus_two_shift]
                +r[index+y_plus_two_z_minus_one_shift]
                +r[index+y_plus_two_z_plus_two_shift]
                +r[index+x_plus_one_y_minus_one_z_minus_one_shift]
                +r[index+x_plus_one_y_minus_one_z_plus_two_shift]
                +r[index+x_plus_one_y_plus_two_z_minus_one_shift]
                +r[index+x_plus_one_y_plus_two_z_plus_two_shift]
                +r[index+x_plus_two_y_minus_one_shift]
                +r[index+x_plus_two_y_minus_one_z_plus_one_shift]
                +r[index+x_plus_two_z_minus_one_shift]
                +r[index+x_plus_two_z_plus_two_shift]
                +r[index+x_plus_two_y_plus_one_z_minus_one_shift]
                +r[index+x_plus_two_y_plus_one_z_plus_two_shift]
                +r[index+x_plus_two_y_plus_two_shift]
                +r[index+x_plus_two_y_plus_two_z_plus_one_shift]);
            b_coarse[coarse_index]+=(T)9*(
                r[index+x_minus_one_shift]
                +r[index+x_minus_one_z_plus_one_shift]
                +r[index+x_minus_one_y_plus_one_shift]
                +r[index+x_minus_one_y_plus_one_z_plus_one_shift]
                +r[index+y_minus_one_shift]
                +r[index+y_minus_one_z_plus_one_shift]
                +r[index+z_minus_one_shift]
                +r[index+z_plus_two_shift]
                +r[index+y_plus_one_z_minus_one_shift]
                +r[index+y_plus_one_z_plus_two_shift]
                +r[index+y_plus_two_shift]
                +r[index+y_plus_two_z_plus_one_shift]
                +r[index+x_plus_one_y_minus_one_shift]
                +r[index+x_plus_one_y_minus_one_z_plus_one_shift]
                +r[index+x_plus_one_z_minus_one_shift]
                +r[index+x_plus_one_z_plus_two_shift]
                +r[index+x_plus_one_y_plus_one_z_minus_one_shift]
                +r[index+x_plus_one_y_plus_one_z_plus_two_shift]
                +r[index+x_plus_one_y_plus_two_shift]
                +r[index+x_plus_one_y_plus_two_z_plus_one_shift]
                +r[index+x_plus_two_shift]
                +r[index+x_plus_two_z_plus_one_shift]
                +r[index+x_plus_two_y_plus_one_shift]
                +r[index+x_plus_two_y_plus_one_z_plus_one_shift]);
            b_coarse[coarse_index]+=(T)27*(
                r[index]
                +r[index+z_plus_one_shift]
                +r[index+y_plus_one_shift]
                +r[index+y_plus_one_z_plus_one_shift]
                +r[index+x_plus_one_shift]
                +r[index+x_plus_one_z_plus_one_shift]
                +r[index+x_plus_one_y_plus_one_shift]
                +r[index+x_plus_one_y_plus_one_z_plus_one_shift]);

	    b_coarse[coarse_index]*=(T)0.0078125;

            if(cell_type_coarse[coarse_index]!=INTERIOR_CELL_TYPE) b_coarse[coarse_index]=T();
        }
}
//#####################################################################
template class Restriction_Helper<float>;
template class Restriction_Helper<double>;
