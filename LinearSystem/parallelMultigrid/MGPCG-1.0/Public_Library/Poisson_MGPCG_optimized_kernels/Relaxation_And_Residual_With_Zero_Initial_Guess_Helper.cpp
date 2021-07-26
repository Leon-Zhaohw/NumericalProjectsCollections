//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis, Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
#include "Relaxation_And_Residual_With_Zero_Initial_Guess_Helper.h"
#include "PTHREAD_QUEUE.h"
#include "Instantiation_Helpers.h"
using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;
//#####################################################################
// Constructor of non-size-specific class
//#####################################################################
CONSTRUCTOR_INSTANTIATION_HELPERX(Relaxation_And_Residual_With_Zero_Initial_Guess,(const int x_size_input,const int y_size_input,const int z_size_input,T* const u,T* const b,T* const r,const unsigned char* const bit_writemask,const unsigned char* const bit_interiormask,const T nullspace_component),(x_size_input,u,b,r,bit_writemask,bit_interiormask,nullspace_component))
//#####################################################################
// Function Run_Parallel
//#####################################################################
namespace{
    
template<class T,int y_size,int z_size>
struct Relaxation_And_Residual_With_Zero_Initial_Guess_Size_Specific_Thread_Helper:public PhysBAM::PTHREAD_QUEUE::TASK
{
    Relaxation_And_Residual_With_Zero_Initial_Guess_Size_Specific_Helper<T,y_size,z_size>* const obj;
    const int xmin,xmax;
    Relaxation_And_Residual_With_Zero_Initial_Guess_Size_Specific_Thread_Helper(
        Relaxation_And_Residual_With_Zero_Initial_Guess_Size_Specific_Helper<T,y_size,z_size>* const obj_input,const int xmin_input,const int xmax_input)
        :obj(obj_input),xmin(xmin_input),xmax(xmax_input) {}
    void Run(){obj->Run_X_Range(xmin,xmax);}
};

template<class T,int y_size,int z_size>
struct Relaxation_With_Zero_Initial_Guess_Size_Specific_Thread_Helper:public PhysBAM::PTHREAD_QUEUE::TASK
{
    Relaxation_And_Residual_With_Zero_Initial_Guess_Size_Specific_Helper<T,y_size,z_size>* const obj;
    const int xmin,xmax;
    Relaxation_With_Zero_Initial_Guess_Size_Specific_Thread_Helper(
        Relaxation_And_Residual_With_Zero_Initial_Guess_Size_Specific_Helper<T,y_size,z_size>* const obj_input,const int xmin_input,const int xmax_input)
        :obj(obj_input),xmin(xmin_input),xmax(xmax_input) {}
    void Run(){obj->Run_X_Range_Relax(xmin,xmax);}
};

template<class T,int y_size,int z_size>
struct Residual_With_Zero_Initial_Guess_Size_Specific_Thread_Helper:public PhysBAM::PTHREAD_QUEUE::TASK
{
    Relaxation_And_Residual_With_Zero_Initial_Guess_Size_Specific_Helper<T,y_size,z_size>* const obj;
    const int xmin,xmax;
    Residual_With_Zero_Initial_Guess_Size_Specific_Thread_Helper(
        Relaxation_And_Residual_With_Zero_Initial_Guess_Size_Specific_Helper<T,y_size,z_size>* const obj_input,const int xmin_input,const int xmax_input)
        :obj(obj_input),xmin(xmin_input),xmax(xmax_input) {}
    void Run(){obj->Run_X_Range_Residual(xmin,xmax);}
};
}
template<class T,int y_size,int z_size> void Relaxation_And_Residual_With_Zero_Initial_Guess_Size_Specific_Helper<T,y_size,z_size>::
Run_Parallel(const int number_of_partitions)
{
    if(x_size*y_size*z_size<=32*32*32) {Run();return;} // Run serially for small problems

    if(nullspace_component){ // need to do relaxation (and nullspace projection) separately from residual computation
	int number_of_x_blocks=coarse_x_size/x_block_size;
	for(int partition=1;partition<=number_of_partitions;partition++){
	    int first_block_of_partition=(number_of_x_blocks/number_of_partitions)*(partition-1)+std::min(number_of_x_blocks%number_of_partitions,partition-1)+1;
	    int last_block_of_partition=(number_of_x_blocks/number_of_partitions)*partition+std::min(number_of_x_blocks%number_of_partitions,partition);
	    int xmin=(first_block_of_partition-1)*x_block_size+1;
	    int xmax=last_block_of_partition*x_block_size;
	    Relaxation_With_Zero_Initial_Guess_Size_Specific_Thread_Helper<T,y_size,z_size>* task=new Relaxation_With_Zero_Initial_Guess_Size_Specific_Thread_Helper<T,y_size,z_size>(this,xmin,xmax);
	    pthread_queue->Queue(task);}
	pthread_queue->Wait();
	
	for(int partition=1;partition<=number_of_partitions;partition++){
	    int first_block_of_partition=(number_of_x_blocks/number_of_partitions)*(partition-1)+std::min(number_of_x_blocks%number_of_partitions,partition-1)+1;
	    int last_block_of_partition=(number_of_x_blocks/number_of_partitions)*partition+std::min(number_of_x_blocks%number_of_partitions,partition);
	    int xmin=(first_block_of_partition-1)*x_block_size+1;
	    int xmax=last_block_of_partition*x_block_size;
	    Residual_With_Zero_Initial_Guess_Size_Specific_Thread_Helper<T,y_size,z_size>* task=new Residual_With_Zero_Initial_Guess_Size_Specific_Thread_Helper<T,y_size,z_size>(this,xmin,xmax);
	    pthread_queue->Queue(task);}
	pthread_queue->Wait();
    }else{	
	int number_of_x_blocks=coarse_x_size/x_block_size;
	for(int partition=1;partition<=number_of_partitions;partition++){
	    int first_block_of_partition=(number_of_x_blocks/number_of_partitions)*(partition-1)+std::min(number_of_x_blocks%number_of_partitions,partition-1)+1;
	    int last_block_of_partition=(number_of_x_blocks/number_of_partitions)*partition+std::min(number_of_x_blocks%number_of_partitions,partition);
	    int xmin=(first_block_of_partition-1)*x_block_size+1;
	    int xmax=last_block_of_partition*x_block_size;
	    Relaxation_And_Residual_With_Zero_Initial_Guess_Size_Specific_Thread_Helper<T,y_size,z_size>* task=new Relaxation_And_Residual_With_Zero_Initial_Guess_Size_Specific_Thread_Helper<T,y_size,z_size>(this,xmin,xmax);
	    pthread_queue->Queue(task);}
	pthread_queue->Wait();
    }
}
//#####################################################################
// Function Run_X_Range
//#####################################################################
template<class T,int y_size,int z_size> void Relaxation_And_Residual_With_Zero_Initial_Guess_Size_Specific_Helper<T,y_size,z_size>::
Run_X_Range(const int xmin,const int xmax)
{    
    const T one_third=1./3.;
    const T one_ninth=1./9.;
    const T minus_one_ninth=-1./9.;

    for(int block_i=xmin;block_i<=xmax;block_i+=x_block_size)
    for(int block_j=1;block_j<=coarse_y_size;block_j+=y_block_size)
    for(int block_k=1;block_k<=coarse_z_size;block_k+=z_block_size)
 
        for(int i=block_i;i<block_i+x_block_size;i++)
        for(int j=block_j;j<block_j+y_block_size;j++)
        for(int k=block_k;k<block_k+z_block_size;k++)
        {            
            int index=2*(i*x_shift+j*y_shift+k*z_shift)-x_shift-y_shift-z_shift;
            int coarse_index=i*coarse_x_shift+j*coarse_y_shift+k*coarse_z_shift;
            
            // Replace solution with relaxed value, or zero if not deeply interior
            if(bit_writemask[coarse_index] & 0x01)
                u[index]=minus_one_ninth*b[index];
            else
                u[index]=T();

            if(bit_writemask[coarse_index] & 0x02)            
                u[index+z_plus_one_shift]=minus_one_ninth*b[index+z_plus_one_shift];
            else
                u[index+z_plus_one_shift]=T();

            if(bit_writemask[coarse_index] & 0x04)
                u[index+y_plus_one_shift]=minus_one_ninth*b[index+y_plus_one_shift];
            else
                u[index+y_plus_one_shift]=T();

            if(bit_writemask[coarse_index] & 0x08)
                u[index+y_plus_one_z_plus_one_shift]=minus_one_ninth*b[index+y_plus_one_z_plus_one_shift];
            else
                u[index+y_plus_one_z_plus_one_shift]=T();

            if(bit_writemask[coarse_index] & 0x10)
                u[index+x_plus_one_shift]=minus_one_ninth*b[index+x_plus_one_shift];
            else
                u[index+x_plus_one_shift]=T();

            if(bit_writemask[coarse_index] & 0x20)
                u[index+x_plus_one_z_plus_one_shift]=minus_one_ninth*b[index+x_plus_one_z_plus_one_shift];
            else
                u[index+x_plus_one_z_plus_one_shift]=T();

            if(bit_writemask[coarse_index] & 0x40)
                u[index+x_plus_one_y_plus_one_shift]=minus_one_ninth*b[index+x_plus_one_y_plus_one_shift];
            else
                u[index+x_plus_one_y_plus_one_shift]=T();

            if(bit_writemask[coarse_index] & 0x80)
                u[index+x_plus_one_y_plus_one_z_plus_one_shift]=minus_one_ninth*b[index+x_plus_one_y_plus_one_z_plus_one_shift];
            else
                u[index+x_plus_one_y_plus_one_z_plus_one_shift]=T();
            
	    // update residual for interior cells
#if 1
	    r[index]=one_ninth*(
                b[index+x_plus_one_shift]
                +b[index+x_minus_one_shift]
                +b[index+y_plus_one_shift]
                +b[index+y_minus_one_shift]
                +b[index+z_plus_one_shift]
                +b[index+z_minus_one_shift]);
            r[index]+=one_third*b[index];

            r[index+z_plus_one_shift]=one_ninth*(
                b[index+x_plus_one_z_plus_one_shift]
                +b[index+x_minus_one_z_plus_one_shift]
                +b[index+y_plus_one_z_plus_one_shift]
                +b[index+y_minus_one_z_plus_one_shift]
                +b[index+z_plus_two_shift]
                +b[index]);
            r[index+z_plus_one_shift]+=one_third*b[index+z_plus_one_shift];

            r[index+y_plus_one_shift]=one_ninth*(
                b[index+x_plus_one_y_plus_one_shift]
                +b[index+x_minus_one_y_plus_one_shift]
                +b[index+y_plus_two_shift]
                +b[index]
                +b[index+y_plus_one_z_plus_one_shift]
                +b[index+y_plus_one_z_minus_one_shift]);
            r[index+y_plus_one_shift]+=one_third*b[index+y_plus_one_shift];

            r[index+y_plus_one_z_plus_one_shift]=one_ninth*(
                b[index+x_plus_one_y_plus_one_z_plus_one_shift]
                +b[index+x_minus_one_y_plus_one_z_plus_one_shift]
                +b[index+y_plus_two_z_plus_one_shift]
                +b[index+z_plus_one_shift]
                +b[index+y_plus_one_z_plus_two_shift]
                +b[index+y_plus_one_shift]);
            r[index+y_plus_one_z_plus_one_shift]+=one_third*b[index+y_plus_one_z_plus_one_shift];

            r[index+x_plus_one_shift]=one_ninth*(
                b[index+x_plus_two_shift]
                +b[index]
                +b[index+x_plus_one_y_plus_one_shift]
                +b[index+x_plus_one_y_minus_one_shift]
                +b[index+x_plus_one_z_plus_one_shift]
                +b[index+x_plus_one_z_minus_one_shift]);
            r[index+x_plus_one_shift]+=one_third*b[index+x_plus_one_shift];

            r[index+x_plus_one_z_plus_one_shift]=one_ninth*(
                b[index+x_plus_two_z_plus_one_shift]
                +b[index+z_plus_one_shift]
                +b[index+x_plus_one_y_plus_one_z_plus_one_shift]
                +b[index+x_plus_one_y_minus_one_z_plus_one_shift]
                +b[index+x_plus_one_z_plus_two_shift]
                +b[index+x_plus_one_shift]);
            r[index+x_plus_one_z_plus_one_shift]+=one_third*b[index+x_plus_one_z_plus_one_shift];

            r[index+x_plus_one_y_plus_one_shift]=one_ninth*(
                b[index+x_plus_two_y_plus_one_shift]
                +b[index+y_plus_one_shift]
                +b[index+x_plus_one_y_plus_two_shift]
                +b[index+x_plus_one_shift]
                +b[index+x_plus_one_y_plus_one_z_plus_one_shift]
                +b[index+x_plus_one_y_plus_one_z_minus_one_shift]);
            r[index+x_plus_one_y_plus_one_shift]+=one_third*b[index+x_plus_one_y_plus_one_shift];

            r[index+x_plus_one_y_plus_one_z_plus_one_shift]=one_ninth*(
                b[index+x_plus_two_y_plus_one_z_plus_one_shift]
                +b[index+y_plus_one_z_plus_one_shift]
                +b[index+x_plus_one_y_plus_two_z_plus_one_shift]
                +b[index+x_plus_one_z_plus_one_shift]
                +b[index+x_plus_one_y_plus_one_z_plus_two_shift]
                +b[index+x_plus_one_y_plus_one_shift]);
            r[index+x_plus_one_y_plus_one_z_plus_one_shift]+=one_third*b[index+x_plus_one_y_plus_one_z_plus_one_shift];

            // Replace residuals with zero outside the domain
            if(!(bit_interiormask[coarse_index] & 0x01))
                r[index]=0;

            if(!(bit_interiormask[coarse_index] & 0x02))
                r[index+z_plus_one_shift]=0;

            if(!(bit_interiormask[coarse_index] & 0x04))
                r[index+y_plus_one_shift]=0;

            if(!(bit_interiormask[coarse_index] & 0x08))
                r[index+y_plus_one_z_plus_one_shift]=0;

            if(!(bit_interiormask[coarse_index] & 0x10))
                r[index+x_plus_one_shift]=0;

            if(!(bit_interiormask[coarse_index] & 0x20))
                r[index+x_plus_one_z_plus_one_shift]=0;

            if(!(bit_interiormask[coarse_index] & 0x40))
                r[index+x_plus_one_y_plus_one_shift]=0;

            if(!(bit_interiormask[coarse_index] & 0x80))
                r[index+x_plus_one_y_plus_one_z_plus_one_shift]=0;
#else
	    r[index]=0;
            if(bit_interiormask[coarse_index] & 0x01)
		r[index]=one_ninth*(
		    b[index+x_plus_one_shift]
		    +b[index+x_minus_one_shift]
		    +b[index+y_plus_one_shift]
		    +b[index+y_minus_one_shift]
		    +b[index+z_plus_one_shift]
		    +b[index+z_minus_one_shift])
		    +one_third*b[index];

	    r[index+z_plus_one_shift]=0;
            if(bit_interiormask[coarse_index] & 0x02)
		r[index+z_plus_one_shift]=one_ninth*(
		    b[index+x_plus_one_z_plus_one_shift]
		    +b[index+x_minus_one_z_plus_one_shift]
		    +b[index+y_plus_one_z_plus_one_shift]
		    +b[index+y_minus_one_z_plus_one_shift]
		    +b[index+z_plus_two_shift]
		    +b[index])
		    +one_third*b[index+z_plus_one_shift];


	    r[index+y_plus_one_shift]=0;
            if(bit_interiormask[coarse_index] & 0x04)
		r[index+y_plus_one_shift]=one_ninth*(
		    b[index+x_plus_one_y_plus_one_shift]
		    +b[index+x_minus_one_y_plus_one_shift]
		    +b[index+y_plus_two_shift]
		    +b[index]
		    +b[index+y_plus_one_z_plus_one_shift]
		    +b[index+y_plus_one_z_minus_one_shift])
		    +one_third*b[index+y_plus_one_shift];
	    

	    r[index+y_plus_one_shift]=0;
            if(bit_interiormask[coarse_index] & 0x08)
		r[index+y_plus_one_z_plus_one_shift]=one_ninth*(
		    b[index+x_plus_one_y_plus_one_z_plus_one_shift]
		    +b[index+x_minus_one_y_plus_one_z_plus_one_shift]
		    +b[index+y_plus_two_z_plus_one_shift]
		    +b[index+z_plus_one_shift]
		    +b[index+y_plus_one_z_plus_two_shift]
		    +b[index+y_plus_one_shift])
		    +one_third*b[index+y_plus_one_z_plus_one_shift];


	    r[index+x_plus_one_shift]=0;
            if(bit_interiormask[coarse_index] & 0x10)
		r[index+x_plus_one_shift]=one_ninth*(
		    b[index+x_plus_two_shift]
		    +b[index]
		    +b[index+x_plus_one_y_plus_one_shift]
		    +b[index+x_plus_one_y_minus_one_shift]
		    +b[index+x_plus_one_z_plus_one_shift]
		    +b[index+x_plus_one_z_minus_one_shift])
		    +one_third*b[index+x_plus_one_shift];


	    r[index+x_plus_one_z_plus_one_shift]=0;
            if(bit_interiormask[coarse_index] & 0x20)
		r[index+x_plus_one_z_plus_one_shift]=one_ninth*(
		    b[index+x_plus_two_z_plus_one_shift]
		    +b[index+z_plus_one_shift]
		    +b[index+x_plus_one_y_plus_one_z_plus_one_shift]
		    +b[index+x_plus_one_y_minus_one_z_plus_one_shift]
		    +b[index+x_plus_one_z_plus_two_shift]
		    +b[index+x_plus_one_shift])
		    +one_third*b[index+x_plus_one_z_plus_one_shift];

	    
	    r[index+x_plus_one_y_plus_one_shift]=0;
            if(bit_interiormask[coarse_index] & 0x40)
		r[index+x_plus_one_y_plus_one_shift]=one_ninth*(
		    b[index+x_plus_two_y_plus_one_shift]
		    +b[index+y_plus_one_shift]
		    +b[index+x_plus_one_y_plus_two_shift]
		    +b[index+x_plus_one_shift]
		    +b[index+x_plus_one_y_plus_one_z_plus_one_shift]
		    +b[index+x_plus_one_y_plus_one_z_minus_one_shift])
		    +one_third*b[index+x_plus_one_y_plus_one_shift];
	    

	    r[index+x_plus_one_y_plus_one_z_plus_one_shift]=0;
            if(bit_interiormask[coarse_index] & 0x80)
		r[index+x_plus_one_y_plus_one_z_plus_one_shift]=one_ninth*(
		    b[index+x_plus_two_y_plus_one_z_plus_one_shift]
		    +b[index+y_plus_one_z_plus_one_shift]
		    +b[index+x_plus_one_y_plus_two_z_plus_one_shift]
		    +b[index+x_plus_one_z_plus_one_shift]
		    +b[index+x_plus_one_y_plus_one_z_plus_two_shift]
		    +b[index+x_plus_one_y_plus_one_shift])
		    +one_third*b[index+x_plus_one_y_plus_one_z_plus_one_shift];
#endif
        }
}
//#####################################################################
// Function Run_X_Range_Relax
//#####################################################################
template<class T,int y_size,int z_size> void Relaxation_And_Residual_With_Zero_Initial_Guess_Size_Specific_Helper<T,y_size,z_size>::
Run_X_Range_Relax(const int xmin,const int xmax)
{    
    const T minus_one_ninth=-1./9.;

    for(int block_i=xmin;block_i<=xmax;block_i+=x_block_size)
    for(int block_j=1;block_j<=coarse_y_size;block_j+=y_block_size)
    for(int block_k=1;block_k<=coarse_z_size;block_k+=z_block_size)
 
        for(int i=block_i;i<block_i+x_block_size;i++)
        for(int j=block_j;j<block_j+y_block_size;j++)
        for(int k=block_k;k<block_k+z_block_size;k++)
        {            
            int index=2*(i*x_shift+j*y_shift+k*z_shift)-x_shift-y_shift-z_shift;
            int coarse_index=i*coarse_x_shift+j*coarse_y_shift+k*coarse_z_shift;
            
            // Subtract nullspace component from the right hand side, for interior indices
            if(bit_interiormask[coarse_index] & 0x01)
                b[index]-=nullspace_component;

            if(bit_interiormask[coarse_index] & 0x02)            
                b[index+z_plus_one_shift]-=nullspace_component;

            if(bit_interiormask[coarse_index] & 0x04)
                b[index+y_plus_one_shift]-=nullspace_component;

            if(bit_interiormask[coarse_index] & 0x08)
                b[index+y_plus_one_z_plus_one_shift]-=nullspace_component;

            if(bit_interiormask[coarse_index] & 0x10)
                b[index+x_plus_one_shift]-=nullspace_component;

            if(bit_interiormask[coarse_index] & 0x20)
                b[index+x_plus_one_z_plus_one_shift]-=nullspace_component;

            if(bit_interiormask[coarse_index] & 0x40)
                b[index+x_plus_one_y_plus_one_shift]-=nullspace_component;

            if(bit_interiormask[coarse_index] & 0x80)
                b[index+x_plus_one_y_plus_one_z_plus_one_shift]-=nullspace_component;
            
            // Replace solution with relaxed value, or zero if not deeply interior
            if(bit_writemask[coarse_index] & 0x01)
                u[index]=minus_one_ninth*b[index];
            else
                u[index]=T();

            if(bit_writemask[coarse_index] & 0x02)            
                u[index+z_plus_one_shift]=minus_one_ninth*b[index+z_plus_one_shift];
            else
                u[index+z_plus_one_shift]=T();

            if(bit_writemask[coarse_index] & 0x04)
                u[index+y_plus_one_shift]=minus_one_ninth*b[index+y_plus_one_shift];
            else
                u[index+y_plus_one_shift]=T();

            if(bit_writemask[coarse_index] & 0x08)
                u[index+y_plus_one_z_plus_one_shift]=minus_one_ninth*b[index+y_plus_one_z_plus_one_shift];
            else
                u[index+y_plus_one_z_plus_one_shift]=T();

            if(bit_writemask[coarse_index] & 0x10)
                u[index+x_plus_one_shift]=minus_one_ninth*b[index+x_plus_one_shift];
            else
                u[index+x_plus_one_shift]=T();

            if(bit_writemask[coarse_index] & 0x20)
                u[index+x_plus_one_z_plus_one_shift]=minus_one_ninth*b[index+x_plus_one_z_plus_one_shift];
            else
                u[index+x_plus_one_z_plus_one_shift]=T();

            if(bit_writemask[coarse_index] & 0x40)
                u[index+x_plus_one_y_plus_one_shift]=minus_one_ninth*b[index+x_plus_one_y_plus_one_shift];
            else
                u[index+x_plus_one_y_plus_one_shift]=T();

            if(bit_writemask[coarse_index] & 0x80)
                u[index+x_plus_one_y_plus_one_z_plus_one_shift]=minus_one_ninth*b[index+x_plus_one_y_plus_one_z_plus_one_shift];
            else
                u[index+x_plus_one_y_plus_one_z_plus_one_shift]=T();
        }
}
//#####################################################################
// Function Run_X_Range_Residual
//#####################################################################
template<class T,int y_size,int z_size> void Relaxation_And_Residual_With_Zero_Initial_Guess_Size_Specific_Helper<T,y_size,z_size>::
Run_X_Range_Residual(const int xmin,const int xmax)
{    
    const T one_third=1./3.;
    const T one_ninth=1./9.;

    for(int block_i=xmin;block_i<=xmax;block_i+=x_block_size)
    for(int block_j=1;block_j<=coarse_y_size;block_j+=y_block_size)
    for(int block_k=1;block_k<=coarse_z_size;block_k+=z_block_size)
 
        for(int i=block_i;i<block_i+x_block_size;i++)
        for(int j=block_j;j<block_j+y_block_size;j++)
        for(int k=block_k;k<block_k+z_block_size;k++)
        {            
            int index=2*(i*x_shift+j*y_shift+k*z_shift)-x_shift-y_shift-z_shift;
            int coarse_index=i*coarse_x_shift+j*coarse_y_shift+k*coarse_z_shift;
            
            
#if 1
	    r[index]=one_ninth*(
                b[index+x_plus_one_shift]
                +b[index+x_minus_one_shift]
                +b[index+y_plus_one_shift]
                +b[index+y_minus_one_shift]
                +b[index+z_plus_one_shift]
                +b[index+z_minus_one_shift]);
            r[index]+=one_third*b[index];

            r[index+z_plus_one_shift]=one_ninth*(
                b[index+x_plus_one_z_plus_one_shift]
                +b[index+x_minus_one_z_plus_one_shift]
                +b[index+y_plus_one_z_plus_one_shift]
                +b[index+y_minus_one_z_plus_one_shift]
                +b[index+z_plus_two_shift]
                +b[index]);
            r[index+z_plus_one_shift]+=one_third*b[index+z_plus_one_shift];

            r[index+y_plus_one_shift]=one_ninth*(
                b[index+x_plus_one_y_plus_one_shift]
                +b[index+x_minus_one_y_plus_one_shift]
                +b[index+y_plus_two_shift]
                +b[index]
                +b[index+y_plus_one_z_plus_one_shift]
                +b[index+y_plus_one_z_minus_one_shift]);
            r[index+y_plus_one_shift]+=one_third*b[index+y_plus_one_shift];

            r[index+y_plus_one_z_plus_one_shift]=one_ninth*(
                b[index+x_plus_one_y_plus_one_z_plus_one_shift]
                +b[index+x_minus_one_y_plus_one_z_plus_one_shift]
                +b[index+y_plus_two_z_plus_one_shift]
                +b[index+z_plus_one_shift]
                +b[index+y_plus_one_z_plus_two_shift]
                +b[index+y_plus_one_shift]);
            r[index+y_plus_one_z_plus_one_shift]+=one_third*b[index+y_plus_one_z_plus_one_shift];

            r[index+x_plus_one_shift]=one_ninth*(
                b[index+x_plus_two_shift]
                +b[index]
                +b[index+x_plus_one_y_plus_one_shift]
                +b[index+x_plus_one_y_minus_one_shift]
                +b[index+x_plus_one_z_plus_one_shift]
                +b[index+x_plus_one_z_minus_one_shift]);
            r[index+x_plus_one_shift]+=one_third*b[index+x_plus_one_shift];

            r[index+x_plus_one_z_plus_one_shift]=one_ninth*(
                b[index+x_plus_two_z_plus_one_shift]
                +b[index+z_plus_one_shift]
                +b[index+x_plus_one_y_plus_one_z_plus_one_shift]
                +b[index+x_plus_one_y_minus_one_z_plus_one_shift]
                +b[index+x_plus_one_z_plus_two_shift]
                +b[index+x_plus_one_shift]);
            r[index+x_plus_one_z_plus_one_shift]+=one_third*b[index+x_plus_one_z_plus_one_shift];

            r[index+x_plus_one_y_plus_one_shift]=one_ninth*(
                b[index+x_plus_two_y_plus_one_shift]
                +b[index+y_plus_one_shift]
                +b[index+x_plus_one_y_plus_two_shift]
                +b[index+x_plus_one_shift]
                +b[index+x_plus_one_y_plus_one_z_plus_one_shift]
                +b[index+x_plus_one_y_plus_one_z_minus_one_shift]);
            r[index+x_plus_one_y_plus_one_shift]+=one_third*b[index+x_plus_one_y_plus_one_shift];

            r[index+x_plus_one_y_plus_one_z_plus_one_shift]=one_ninth*(
                b[index+x_plus_two_y_plus_one_z_plus_one_shift]
                +b[index+y_plus_one_z_plus_one_shift]
                +b[index+x_plus_one_y_plus_two_z_plus_one_shift]
                +b[index+x_plus_one_z_plus_one_shift]
                +b[index+x_plus_one_y_plus_one_z_plus_two_shift]
                +b[index+x_plus_one_y_plus_one_shift]);
            r[index+x_plus_one_y_plus_one_z_plus_one_shift]+=one_third*b[index+x_plus_one_y_plus_one_z_plus_one_shift];

            // Replace residuals with zero outside the domain
            if(!(bit_interiormask[coarse_index] & 0x01))
                r[index]=0;

            if(!(bit_interiormask[coarse_index] & 0x02))
                r[index+z_plus_one_shift]=0;

            if(!(bit_interiormask[coarse_index] & 0x04))
                r[index+y_plus_one_shift]=0;

            if(!(bit_interiormask[coarse_index] & 0x08))
                r[index+y_plus_one_z_plus_one_shift]=0;

            if(!(bit_interiormask[coarse_index] & 0x10))
                r[index+x_plus_one_shift]=0;

            if(!(bit_interiormask[coarse_index] & 0x20))
                r[index+x_plus_one_z_plus_one_shift]=0;

            if(!(bit_interiormask[coarse_index] & 0x40))
                r[index+x_plus_one_y_plus_one_shift]=0;

            if(!(bit_interiormask[coarse_index] & 0x80))
                r[index+x_plus_one_y_plus_one_z_plus_one_shift]=0;
#else
	    r[index]=0;
            if(bit_interiormask[coarse_index] & 0x01)
		r[index]=one_ninth*(
		    b[index+x_plus_one_shift]
		    +b[index+x_minus_one_shift]
		    +b[index+y_plus_one_shift]
		    +b[index+y_minus_one_shift]
		    +b[index+z_plus_one_shift]
		    +b[index+z_minus_one_shift])
		    +one_third*b[index];

	    r[index+z_plus_one_shift]=0;
            if(bit_interiormask[coarse_index] & 0x02)
		r[index+z_plus_one_shift]=one_ninth*(
		    b[index+x_plus_one_z_plus_one_shift]
		    +b[index+x_minus_one_z_plus_one_shift]
		    +b[index+y_plus_one_z_plus_one_shift]
		    +b[index+y_minus_one_z_plus_one_shift]
		    +b[index+z_plus_two_shift]
		    +b[index])
		    +one_third*b[index+z_plus_one_shift];


	    r[index+y_plus_one_shift]=0;
            if(bit_interiormask[coarse_index] & 0x04)
		r[index+y_plus_one_shift]=one_ninth*(
		    b[index+x_plus_one_y_plus_one_shift]
		    +b[index+x_minus_one_y_plus_one_shift]
		    +b[index+y_plus_two_shift]
		    +b[index]
		    +b[index+y_plus_one_z_plus_one_shift]
		    +b[index+y_plus_one_z_minus_one_shift])
		    +one_third*b[index+y_plus_one_shift];
	    

	    r[index+y_plus_one_shift]=0;
            if(bit_interiormask[coarse_index] & 0x08)
		r[index+y_plus_one_z_plus_one_shift]=one_ninth*(
		    b[index+x_plus_one_y_plus_one_z_plus_one_shift]
		    +b[index+x_minus_one_y_plus_one_z_plus_one_shift]
		    +b[index+y_plus_two_z_plus_one_shift]
		    +b[index+z_plus_one_shift]
		    +b[index+y_plus_one_z_plus_two_shift]
		    +b[index+y_plus_one_shift])
		    +one_third*b[index+y_plus_one_z_plus_one_shift];


	    r[index+x_plus_one_shift]=0;
            if(bit_interiormask[coarse_index] & 0x10)
		r[index+x_plus_one_shift]=one_ninth*(
		    b[index+x_plus_two_shift]
		    +b[index]
		    +b[index+x_plus_one_y_plus_one_shift]
		    +b[index+x_plus_one_y_minus_one_shift]
		    +b[index+x_plus_one_z_plus_one_shift]
		    +b[index+x_plus_one_z_minus_one_shift])
		    +one_third*b[index+x_plus_one_shift];


	    r[index+x_plus_one_z_plus_one_shift]=0;
            if(bit_interiormask[coarse_index] & 0x20)
		r[index+x_plus_one_z_plus_one_shift]=one_ninth*(
		    b[index+x_plus_two_z_plus_one_shift]
		    +b[index+z_plus_one_shift]
		    +b[index+x_plus_one_y_plus_one_z_plus_one_shift]
		    +b[index+x_plus_one_y_minus_one_z_plus_one_shift]
		    +b[index+x_plus_one_z_plus_two_shift]
		    +b[index+x_plus_one_shift])
		    +one_third*b[index+x_plus_one_z_plus_one_shift];

	    
	    r[index+x_plus_one_y_plus_one_shift]=0;
            if(bit_interiormask[coarse_index] & 0x40)
		r[index+x_plus_one_y_plus_one_shift]=one_ninth*(
		    b[index+x_plus_two_y_plus_one_shift]
		    +b[index+y_plus_one_shift]
		    +b[index+x_plus_one_y_plus_two_shift]
		    +b[index+x_plus_one_shift]
		    +b[index+x_plus_one_y_plus_one_z_plus_one_shift]
		    +b[index+x_plus_one_y_plus_one_z_minus_one_shift])
		    +one_third*b[index+x_plus_one_y_plus_one_shift];
	    

	    r[index+x_plus_one_y_plus_one_z_plus_one_shift]=0;
            if(bit_interiormask[coarse_index] & 0x80)
		r[index+x_plus_one_y_plus_one_z_plus_one_shift]=one_ninth*(
		    b[index+x_plus_two_y_plus_one_z_plus_one_shift]
		    +b[index+y_plus_one_z_plus_one_shift]
		    +b[index+x_plus_one_y_plus_two_z_plus_one_shift]
		    +b[index+x_plus_one_z_plus_one_shift]
		    +b[index+x_plus_one_y_plus_one_z_plus_two_shift]
		    +b[index+x_plus_one_y_plus_one_shift])
		    +one_third*b[index+x_plus_one_y_plus_one_z_plus_one_shift];
#endif
        }
}
//#####################################################################
//#####################################################################
template class Relaxation_And_Residual_With_Zero_Initial_Guess_Helper<float>;
template class Relaxation_And_Residual_With_Zero_Initial_Guess_Helper<double>;

