//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis, Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
#include "Relaxation_Interior_Helper.h"
#include "PTHREAD_QUEUE.h"
#include "Instantiation_Helpers.h"
using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;
//#####################################################################
// Constructor of non-size-specific class
//#####################################################################
CONSTRUCTOR_INSTANTIATION_HELPERX(Relaxation_Interior,(const int x_size_input,const int y_size_input,const int z_size_input,T* const u,const T* const b,T* const delta,const unsigned char* const bit_writemask),(x_size_input,u,b,delta,bit_writemask))
//#####################################################################
// Function Run_Parallel
//#####################################################################
namespace{
template<class T,int y_size,int z_size>
struct Relaxation_Interior_Run_Major_Size_Specific_Thread_Helper:public PhysBAM::PTHREAD_QUEUE::TASK
{
    Relaxation_Interior_Size_Specific_Helper<T,y_size,z_size>* const obj;
    const int xmin,xmax;
    Relaxation_Interior_Run_Major_Size_Specific_Thread_Helper(
        Relaxation_Interior_Size_Specific_Helper<T,y_size,z_size>* const obj_input,const int xmin_input,const int xmax_input)
        :obj(obj_input),xmin(xmin_input),xmax(xmax_input) {}
    void Run(){obj->Run_Major_X_Range(xmin,xmax);}
};
template<class T,int y_size,int z_size>
struct Relaxation_Interior_Run_Minor_Size_Specific_Thread_Helper:public PhysBAM::PTHREAD_QUEUE::TASK
{
    Relaxation_Interior_Size_Specific_Helper<T,y_size,z_size>* const obj;
    const int xmin,xmax;
    Relaxation_Interior_Run_Minor_Size_Specific_Thread_Helper(
        Relaxation_Interior_Size_Specific_Helper<T,y_size,z_size>* const obj_input,const int xmin_input,const int xmax_input)
        :obj(obj_input),xmin(xmin_input),xmax(xmax_input) {}
    void Run(){obj->Run_Minor_X_Range(xmin,xmax);}
};
}
template<class T,int y_size,int z_size> void Relaxation_Interior_Size_Specific_Helper<T,y_size,z_size>::
Run_Parallel(const int number_of_partitions)
{
    if(x_size*y_size*z_size<=32*32*32) {Run();return;} // Run serially for small problems

    int number_of_x_blocks=coarse_x_size/x_block_size;
    for(int partition=1;partition<=number_of_partitions;partition++){
        int first_block_of_partition=(number_of_x_blocks/number_of_partitions)*(partition-1)+std::min(number_of_x_blocks%number_of_partitions,partition-1)+1;
        int last_block_of_partition=(number_of_x_blocks/number_of_partitions)*partition+std::min(number_of_x_blocks%number_of_partitions,partition);
        int xmin=(first_block_of_partition-1)*x_block_size+1;
        int xmax=last_block_of_partition*x_block_size;
        Relaxation_Interior_Run_Major_Size_Specific_Thread_Helper<T,y_size,z_size>* task=new Relaxation_Interior_Run_Major_Size_Specific_Thread_Helper<T,y_size,z_size>(this,xmin,xmax);
        pthread_queue->Queue(task);}
    pthread_queue->Wait();

    for(int partition=1;partition<=number_of_partitions;partition++){
        int first_block_of_partition=(number_of_x_blocks/number_of_partitions)*(partition-1)+std::min(number_of_x_blocks%number_of_partitions,partition-1)+1;
        int last_block_of_partition=(number_of_x_blocks/number_of_partitions)*partition+std::min(number_of_x_blocks%number_of_partitions,partition);
        int xmin=(first_block_of_partition-1)*x_block_size+1;
        int xmax=last_block_of_partition*x_block_size;
        Relaxation_Interior_Run_Minor_Size_Specific_Thread_Helper<T,y_size,z_size>* task=new Relaxation_Interior_Run_Minor_Size_Specific_Thread_Helper<T,y_size,z_size>(this,xmin,xmax);
        pthread_queue->Queue(task);}
    pthread_queue->Wait();
}
//#####################################################################
// Function Run_Major_X_Range
//#####################################################################
template<class T,int y_size,int z_size> void Relaxation_Interior_Size_Specific_Helper<T,y_size,z_size>::
Run_Major_X_Range(const int xmin,const int xmax)
{ 
    const T one_sixth=1./6.;
    const T two_thirds=2./3.;

    int i=xmin;
    for(int j=1;j<=coarse_y_size;j++)
    for(int k=1;k<=coarse_z_size;k++)
    {
        int index=2*(i*x_shift+j*y_shift+k*z_shift)-x_plus_one_y_plus_one_z_plus_one_shift;

        delta[index]=one_sixth*(
            u[index+x_plus_one_shift]
            +u[index+x_minus_one_shift]
            +u[index+y_plus_one_shift]
            +u[index+y_minus_one_shift]
            +u[index+z_plus_one_shift]
            +u[index+z_minus_one_shift]
            -b[index]);
        delta[index]-=u[index];

        delta[index+z_plus_one_shift]=one_sixth*(
            u[index+x_plus_one_z_plus_one_shift]
            +u[index+x_minus_one_z_plus_one_shift]
            +u[index+y_plus_one_z_plus_one_shift]
            +u[index+y_minus_one_z_plus_one_shift]
            +u[index+z_plus_two_shift]
            +u[index]
            -b[index+z_plus_one_shift]);
        delta[index+z_plus_one_shift]-=u[index+z_plus_one_shift];

        delta[index+y_plus_one_shift]=one_sixth*(
            u[index+x_plus_one_y_plus_one_shift]
            +u[index+x_minus_one_y_plus_one_shift]
            +u[index+y_plus_two_shift]
            +u[index]
            +u[index+y_plus_one_z_plus_one_shift]
            +u[index+y_plus_one_z_minus_one_shift]
            -b[index+y_plus_one_shift]);
        delta[index+y_plus_one_shift]-=u[index+y_plus_one_shift];
    
        delta[index+y_plus_one_z_plus_one_shift]=one_sixth*(
            u[index+x_plus_one_y_plus_one_z_plus_one_shift]
            +u[index+x_minus_one_y_plus_one_z_plus_one_shift]
            +u[index+y_plus_two_z_plus_one_shift]
            +u[index+z_plus_one_shift]
            +u[index+y_plus_one_z_plus_two_shift]
            +u[index+y_plus_one_shift]
            -b[index+y_plus_one_z_plus_one_shift]);
        delta[index+y_plus_one_z_plus_one_shift]-=u[index+y_plus_one_z_plus_one_shift];

        delta[index+x_plus_one_shift]=one_sixth*(
            u[index+x_plus_two_shift]
            +u[index]
            +u[index+x_plus_one_y_plus_one_shift]
            +u[index+x_plus_one_y_minus_one_shift]
            +u[index+x_plus_one_z_plus_one_shift]
            +u[index+x_plus_one_z_minus_one_shift]
            -b[index+x_plus_one_shift]);
        delta[index+x_plus_one_shift]-=u[index+x_plus_one_shift];
    
        delta[index+x_plus_one_z_plus_one_shift]=one_sixth*(
            u[index+x_plus_two_z_plus_one_shift]
            +u[index+z_plus_one_shift]
            +u[index+x_plus_one_y_plus_one_z_plus_one_shift]
            +u[index+x_plus_one_y_minus_one_z_plus_one_shift]
            +u[index+x_plus_one_z_plus_two_shift]
            +u[index+x_plus_one_shift]
            -b[index+x_plus_one_z_plus_one_shift]);
        delta[index+x_plus_one_z_plus_one_shift]-=u[index+x_plus_one_z_plus_one_shift];

        delta[index+x_plus_one_y_plus_one_shift]=one_sixth*(
            u[index+x_plus_two_y_plus_one_shift]
            +u[index+y_plus_one_shift]
            +u[index+x_plus_one_y_plus_two_shift]
            +u[index+x_plus_one_shift]
            +u[index+x_plus_one_y_plus_one_z_plus_one_shift]
            +u[index+x_plus_one_y_plus_one_z_minus_one_shift]
            -b[index+x_plus_one_y_plus_one_shift]);
        delta[index+x_plus_one_y_plus_one_shift]-=u[index+x_plus_one_y_plus_one_shift];

        delta[index+x_plus_one_y_plus_one_z_plus_one_shift]=one_sixth*(
            u[index+x_plus_two_y_plus_one_z_plus_one_shift]
            +u[index+y_plus_one_z_plus_one_shift]
            +u[index+x_plus_one_y_plus_two_z_plus_one_shift]
            +u[index+x_plus_one_z_plus_one_shift]
            +u[index+x_plus_one_y_plus_one_z_plus_two_shift]
            +u[index+x_plus_one_y_plus_one_shift]
            -b[index+x_plus_one_y_plus_one_z_plus_one_shift]);
        delta[index+x_plus_one_y_plus_one_z_plus_one_shift]-=u[index+x_plus_one_y_plus_one_z_plus_one_shift];
    }

    for(int i=xmin+1;i<=xmax;i++)
    for(int j=1;j<=coarse_y_size;j++)
    for(int k=1;k<=coarse_z_size;k++)
    {
        int index=2*(i*x_shift+j*y_shift+k*z_shift)-x_plus_one_y_plus_one_z_plus_one_shift;

        delta[index]=one_sixth*(
            u[index+x_plus_one_shift]
            +u[index+x_minus_one_shift]
            +u[index+y_plus_one_shift]
            +u[index+y_minus_one_shift]
            +u[index+z_plus_one_shift]
            +u[index+z_minus_one_shift]
            -b[index]);
        delta[index]-=u[index];

        delta[index+z_plus_one_shift]=one_sixth*(
            u[index+x_plus_one_z_plus_one_shift]
            +u[index+x_minus_one_z_plus_one_shift]
            +u[index+y_plus_one_z_plus_one_shift]
            +u[index+y_minus_one_z_plus_one_shift]
            +u[index+z_plus_two_shift]
            +u[index]
            -b[index+z_plus_one_shift]);
        delta[index+z_plus_one_shift]-=u[index+z_plus_one_shift];

        delta[index+y_plus_one_shift]=one_sixth*(
            u[index+x_plus_one_y_plus_one_shift]
            +u[index+x_minus_one_y_plus_one_shift]
            +u[index+y_plus_two_shift]
            +u[index]
            +u[index+y_plus_one_z_plus_one_shift]
            +u[index+y_plus_one_z_minus_one_shift]
            -b[index+y_plus_one_shift]);
        delta[index+y_plus_one_shift]-=u[index+y_plus_one_shift];
    
        delta[index+y_plus_one_z_plus_one_shift]=one_sixth*(
            u[index+x_plus_one_y_plus_one_z_plus_one_shift]
            +u[index+x_minus_one_y_plus_one_z_plus_one_shift]
            +u[index+y_plus_two_z_plus_one_shift]
            +u[index+z_plus_one_shift]
            +u[index+y_plus_one_z_plus_two_shift]
            +u[index+y_plus_one_shift]
            -b[index+y_plus_one_z_plus_one_shift]);
        delta[index+y_plus_one_z_plus_one_shift]-=u[index+y_plus_one_z_plus_one_shift];

        delta[index+x_plus_one_shift]=one_sixth*(
            u[index+x_plus_two_shift]
            +u[index]
            +u[index+x_plus_one_y_plus_one_shift]
            +u[index+x_plus_one_y_minus_one_shift]
            +u[index+x_plus_one_z_plus_one_shift]
            +u[index+x_plus_one_z_minus_one_shift]
            -b[index+x_plus_one_shift]);
        delta[index+x_plus_one_shift]-=u[index+x_plus_one_shift];
    
        delta[index+x_plus_one_z_plus_one_shift]=one_sixth*(
            u[index+x_plus_two_z_plus_one_shift]
            +u[index+z_plus_one_shift]
            +u[index+x_plus_one_y_plus_one_z_plus_one_shift]
            +u[index+x_plus_one_y_minus_one_z_plus_one_shift]
            +u[index+x_plus_one_z_plus_two_shift]
            +u[index+x_plus_one_shift]
            -b[index+x_plus_one_z_plus_one_shift]);
        delta[index+x_plus_one_z_plus_one_shift]-=u[index+x_plus_one_z_plus_one_shift];

        delta[index+x_plus_one_y_plus_one_shift]=one_sixth*(
            u[index+x_plus_two_y_plus_one_shift]
            +u[index+y_plus_one_shift]
            +u[index+x_plus_one_y_plus_two_shift]
            +u[index+x_plus_one_shift]
            +u[index+x_plus_one_y_plus_one_z_plus_one_shift]
            +u[index+x_plus_one_y_plus_one_z_minus_one_shift]
            -b[index+x_plus_one_y_plus_one_shift]);
        delta[index+x_plus_one_y_plus_one_shift]-=u[index+x_plus_one_y_plus_one_shift];

        delta[index+x_plus_one_y_plus_one_z_plus_one_shift]=one_sixth*(
            u[index+x_plus_two_y_plus_one_z_plus_one_shift]
            +u[index+y_plus_one_z_plus_one_shift]
            +u[index+x_plus_one_y_plus_two_z_plus_one_shift]
            +u[index+x_plus_one_z_plus_one_shift]
            +u[index+x_plus_one_y_plus_one_z_plus_two_shift]
            +u[index+x_plus_one_y_plus_one_shift]
            -b[index+x_plus_one_y_plus_one_z_plus_one_shift]);
        delta[index+x_plus_one_y_plus_one_z_plus_one_shift]-=u[index+x_plus_one_y_plus_one_z_plus_one_shift];

        int shifted_index=index-2*x_shift;
        int coarse_shifted_index=(i-1)*coarse_x_shift+j*coarse_y_shift+k*coarse_z_shift;
	
        if(bit_writemask[coarse_shifted_index] & 0x01)
            u[shifted_index]+=two_thirds*delta[shifted_index];
        if(bit_writemask[coarse_shifted_index] & 0x02)            
            u[shifted_index+z_plus_one_shift]+=two_thirds*delta[shifted_index+z_plus_one_shift];
        if(bit_writemask[coarse_shifted_index] & 0x04)
            u[shifted_index+y_plus_one_shift]+=two_thirds*delta[shifted_index+y_plus_one_shift];
        if(bit_writemask[coarse_shifted_index] & 0x08)
            u[shifted_index+y_plus_one_z_plus_one_shift]+=two_thirds*delta[shifted_index+y_plus_one_z_plus_one_shift];
        if(bit_writemask[coarse_shifted_index] & 0x10)
            u[shifted_index+x_plus_one_shift]+=two_thirds*delta[shifted_index+x_plus_one_shift];
        if(bit_writemask[coarse_shifted_index] & 0x20)
            u[shifted_index+x_plus_one_z_plus_one_shift]+=two_thirds*delta[shifted_index+x_plus_one_z_plus_one_shift];
        if(bit_writemask[coarse_shifted_index] & 0x40)
            u[shifted_index+x_plus_one_y_plus_one_shift]+=two_thirds*delta[shifted_index+x_plus_one_y_plus_one_shift];
        if(bit_writemask[coarse_shifted_index] & 0x80)
            u[shifted_index+x_plus_one_y_plus_one_z_plus_one_shift]+=two_thirds*delta[shifted_index+x_plus_one_y_plus_one_z_plus_one_shift];
    }
}
//#####################################################################
// Function Run_Minor_X_Range
//#####################################################################
template<class T,int y_size,int z_size> void Relaxation_Interior_Size_Specific_Helper<T,y_size,z_size>::
Run_Minor_X_Range(const int xmin,const int xmax)
{    
    const T two_thirds=2./3.;

    int i=xmax;
    for(int block_j=1;block_j<=coarse_y_size;block_j+=y_block_size)
    for(int block_k=1;block_k<=coarse_z_size;block_k+=z_block_size)

        for(int j=block_j;j<block_j+y_block_size;j++)
        for(int k=block_k;k<block_k+z_block_size;k++)
        {
	    int index=2*(i*x_shift+j*y_shift+k*z_shift)-x_shift-y_shift-z_shift;
            int coarse_index=i*coarse_x_shift+j*coarse_y_shift+k*coarse_z_shift;

            if(bit_writemask[coarse_index] & 0x01)
                u[index]+=two_thirds*delta[index];
            if(bit_writemask[coarse_index] & 0x02)            
                u[index+z_plus_one_shift]+=two_thirds*delta[index+z_plus_one_shift];
            if(bit_writemask[coarse_index] & 0x04)
                u[index+y_plus_one_shift]+=two_thirds*delta[index+y_plus_one_shift];
            if(bit_writemask[coarse_index] & 0x08)
                u[index+y_plus_one_z_plus_one_shift]+=two_thirds*delta[index+y_plus_one_z_plus_one_shift];
            if(bit_writemask[coarse_index] & 0x10)
                u[index+x_plus_one_shift]+=two_thirds*delta[index+x_plus_one_shift];
            if(bit_writemask[coarse_index] & 0x20)
                u[index+x_plus_one_z_plus_one_shift]+=two_thirds*delta[index+x_plus_one_z_plus_one_shift];
            if(bit_writemask[coarse_index] & 0x40)
                u[index+x_plus_one_y_plus_one_shift]+=two_thirds*delta[index+x_plus_one_y_plus_one_shift];
            if(bit_writemask[coarse_index] & 0x80)
                u[index+x_plus_one_y_plus_one_z_plus_one_shift]+=two_thirds*delta[index+x_plus_one_y_plus_one_z_plus_one_shift];
	}
}
//#####################################################################
template class Relaxation_Interior_Helper<float>;
template class Relaxation_Interior_Helper<double>;

