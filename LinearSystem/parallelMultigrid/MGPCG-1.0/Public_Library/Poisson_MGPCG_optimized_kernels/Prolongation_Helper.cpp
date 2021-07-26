//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis, Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
#include "Prolongation_Helper.h"
#include "PTHREAD_QUEUE.h"
#include "Instantiation_Helpers.h"
using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;
//#####################################################################
// Constructor/Desctructor of non-size-specific class
//#####################################################################
CONSTRUCTOR_INSTANTIATION_HELPERX(Prolongation,(const int x_size_input,const int y_size_input,const int z_size_input,T* const u,const T* const u_coarse,const unsigned char* const bit_writemask),(x_size_input,u,u_coarse,bit_writemask))
//#####################################################################
// Function Run_Parallel
//#####################################################################
namespace{
template<class T,int y_size,int z_size>
struct Prolongation_Size_Specific_Thread_Helper:public PhysBAM::PTHREAD_QUEUE::TASK
{
    Prolongation_Size_Specific_Helper<T,y_size,z_size>* const obj;
    const int xmin,xmax;
    Prolongation_Size_Specific_Thread_Helper(Prolongation_Size_Specific_Helper<T,y_size,z_size>* const obj_input,const int xmin_input,const int xmax_input)
        :obj(obj_input),xmin(xmin_input),xmax(xmax_input) {}
    void Run(){obj->Run_X_Range(xmin,xmax);}
};
}
template<class T,int y_size,int z_size> void Prolongation_Size_Specific_Helper<T,y_size,z_size>::
Run_Parallel(const int number_of_partitions)
{
    if(x_size*y_size*z_size<=32*32*32) {Run();return;} // Run serially for small problems

    int number_of_x_blocks=coarse_x_size/x_block_size;
    for(int partition=1;partition<=number_of_partitions;partition++){
        int first_block_of_partition=(number_of_x_blocks/number_of_partitions)*(partition-1)+std::min(number_of_x_blocks%number_of_partitions,partition-1)+1;
        int last_block_of_partition=(number_of_x_blocks/number_of_partitions)*partition+std::min(number_of_x_blocks%number_of_partitions,partition);
        int xmin=(first_block_of_partition-1)*x_block_size+1;
        int xmax=last_block_of_partition*x_block_size;
        Prolongation_Size_Specific_Thread_Helper<T,y_size,z_size>* task=new Prolongation_Size_Specific_Thread_Helper<T,y_size,z_size>(this,xmin,xmax);
        pthread_queue->Queue(task);}
    pthread_queue->Wait();
}
//#####################################################################
// Function Run_X_Range
//#####################################################################
template<class T,int y_size,int z_size> void Prolongation_Size_Specific_Helper<T,y_size,z_size>::
Run_X_Range(const int xmin,const int xmax)
{    
    const T one_over_sixty_four=1./64.;
    const T three_over_sixty_four=3./64.;
    const T nine_over_sixty_four=9./64.;
    const T twenty_seven_over_sixty_four=27./64.;

    for(int block_i=xmin;block_i<=xmax;block_i+=x_block_size)
    for(int block_j=1;block_j<=coarse_y_size;block_j+=y_block_size)
    for(int block_k=1;block_k<=coarse_z_size;block_k+=z_block_size)
 
        for(int i=block_i;i<block_i+x_block_size;i++)
        for(int j=block_j;j<block_j+y_block_size;j++)
        for(int k=block_k;k<block_k+z_block_size;k++)
        {    
            int coarse_index=i*coarse_x_shift+j*coarse_y_shift+k*coarse_z_shift;
            int index=2*(i*x_shift+j*y_shift+k*z_shift)-x_shift-y_shift-z_shift;            

            u[index]+=
                one_over_sixty_four*(
                    u_coarse[coarse_index+coarse_x_minus_one_coarse_y_minus_one_coarse_z_minus_one_shift])
                +three_over_sixty_four*(
                    u_coarse[coarse_index+coarse_y_minus_one_coarse_z_minus_one_shift]
                    +u_coarse[coarse_index+coarse_x_minus_one_coarse_z_minus_one_shift]
                    +u_coarse[coarse_index+coarse_x_minus_one_coarse_y_minus_one_shift])
                +nine_over_sixty_four*(
                    u_coarse[coarse_index+coarse_z_minus_one_shift]
                    +u_coarse[coarse_index+coarse_y_minus_one_shift]
                    +u_coarse[coarse_index+coarse_x_minus_one_shift])
                +twenty_seven_over_sixty_four*(
                    u_coarse[coarse_index]);

	    u[index+z_plus_one_shift]+=
                one_over_sixty_four*(
                    u_coarse[coarse_index+coarse_x_minus_one_coarse_y_minus_one_coarse_z_plus_one_shift])
                +three_over_sixty_four*(
                    u_coarse[coarse_index+coarse_y_minus_one_coarse_z_plus_one_shift]
                    +u_coarse[coarse_index+coarse_x_minus_one_coarse_z_plus_one_shift]
                    +u_coarse[coarse_index+coarse_x_minus_one_coarse_y_minus_one_shift])
                +nine_over_sixty_four*(
                    u_coarse[coarse_index+coarse_z_plus_one_shift]
                    +u_coarse[coarse_index+coarse_y_minus_one_shift]
                    +u_coarse[coarse_index+coarse_x_minus_one_shift])
                +twenty_seven_over_sixty_four*(
                    u_coarse[coarse_index]);

            u[index+y_plus_one_shift]+=
                one_over_sixty_four*(
                    u_coarse[coarse_index+coarse_x_minus_one_coarse_y_plus_one_coarse_z_minus_one_shift])
                +three_over_sixty_four*(
                    u_coarse[coarse_index+coarse_y_plus_one_coarse_z_minus_one_shift]
                    +u_coarse[coarse_index+coarse_x_minus_one_coarse_y_plus_one_shift]
                    +u_coarse[coarse_index+coarse_x_minus_one_coarse_z_minus_one_shift])
                +nine_over_sixty_four*(
                    u_coarse[coarse_index+coarse_y_plus_one_shift]
                    +u_coarse[coarse_index+coarse_z_minus_one_shift]
                    +u_coarse[coarse_index+coarse_x_minus_one_shift])
                +twenty_seven_over_sixty_four*(
                    u_coarse[coarse_index]);

            u[index+y_plus_one_z_plus_one_shift]+=
                one_over_sixty_four*(
                    u_coarse[coarse_index+coarse_x_minus_one_coarse_y_plus_one_coarse_z_plus_one_shift])
                +three_over_sixty_four*(
                    u_coarse[coarse_index+coarse_y_plus_one_coarse_z_plus_one_shift]
                    +u_coarse[coarse_index+coarse_x_minus_one_coarse_y_plus_one_shift]
                    +u_coarse[coarse_index+coarse_x_minus_one_coarse_z_plus_one_shift])
                +nine_over_sixty_four*(
                    u_coarse[coarse_index+coarse_y_plus_one_shift]
                    +u_coarse[coarse_index+coarse_z_plus_one_shift]
                    +u_coarse[coarse_index+coarse_x_minus_one_shift])
                +twenty_seven_over_sixty_four*(
                    u_coarse[coarse_index]);

            u[index+x_plus_one_shift]+=
                one_over_sixty_four*(
                    u_coarse[coarse_index+coarse_x_plus_one_coarse_y_minus_one_coarse_z_minus_one_shift])
                +three_over_sixty_four*(
                    u_coarse[coarse_index+coarse_x_plus_one_coarse_z_minus_one_shift]
                    +u_coarse[coarse_index+coarse_x_plus_one_coarse_y_minus_one_shift]
                    +u_coarse[coarse_index+coarse_y_minus_one_coarse_z_minus_one_shift])
                +nine_over_sixty_four*(
                    u_coarse[coarse_index+coarse_x_plus_one_shift]
                    +u_coarse[coarse_index+coarse_z_minus_one_shift]
                    +u_coarse[coarse_index+coarse_y_minus_one_shift])
                +twenty_seven_over_sixty_four*(
                    u_coarse[coarse_index]);

            u[index+x_plus_one_z_plus_one_shift]+=
                one_over_sixty_four*(
                    u_coarse[coarse_index+coarse_x_plus_one_coarse_y_minus_one_coarse_z_plus_one_shift])
                +three_over_sixty_four*(
                    u_coarse[coarse_index+coarse_x_plus_one_coarse_z_plus_one_shift]
                    +u_coarse[coarse_index+coarse_x_plus_one_coarse_y_minus_one_shift]
                    +u_coarse[coarse_index+coarse_y_minus_one_coarse_z_plus_one_shift])
                +nine_over_sixty_four*(
                    u_coarse[coarse_index+coarse_x_plus_one_shift]
                    +u_coarse[coarse_index+coarse_z_plus_one_shift]
                    +u_coarse[coarse_index+coarse_y_minus_one_shift])
                +twenty_seven_over_sixty_four*(
                    u_coarse[coarse_index]);

            u[index+x_plus_one_y_plus_one_shift]+=
                one_over_sixty_four*(
                    u_coarse[coarse_index+coarse_x_plus_one_coarse_y_plus_one_coarse_z_minus_one_shift])
                +three_over_sixty_four*(
                    u_coarse[coarse_index+coarse_x_plus_one_coarse_y_plus_one_shift]
                    +u_coarse[coarse_index+coarse_x_plus_one_coarse_z_minus_one_shift]
                    +u_coarse[coarse_index+coarse_y_plus_one_coarse_z_minus_one_shift])
                +nine_over_sixty_four*(
                    u_coarse[coarse_index+coarse_x_plus_one_shift]
                    +u_coarse[coarse_index+coarse_y_plus_one_shift]
                    +u_coarse[coarse_index+coarse_z_minus_one_shift])
                +twenty_seven_over_sixty_four*(
                    u_coarse[coarse_index]);

            u[index+x_plus_one_y_plus_one_z_plus_one_shift]+=
                one_over_sixty_four*(
                    u_coarse[coarse_index+coarse_x_plus_one_coarse_y_plus_one_coarse_z_plus_one_shift])
                +three_over_sixty_four*(
                    u_coarse[coarse_index+coarse_x_plus_one_coarse_y_plus_one_shift]
                    +u_coarse[coarse_index+coarse_x_plus_one_coarse_z_plus_one_shift]
                    +u_coarse[coarse_index+coarse_y_plus_one_coarse_z_plus_one_shift])
                +nine_over_sixty_four*(
                    u_coarse[coarse_index+coarse_x_plus_one_shift]
                    +u_coarse[coarse_index+coarse_y_plus_one_shift]
                    +u_coarse[coarse_index+coarse_z_plus_one_shift])
                +twenty_seven_over_sixty_four*(
                    u_coarse[coarse_index]);

	    if(!(bit_writemask[coarse_index]&0x01))
		u[index]=T();
	    if(!(bit_writemask[coarse_index] & 0x02))            
                u[index+z_shift]=T();
            if(!(bit_writemask[coarse_index] & 0x04))
                u[index+y_shift]=T();
            if(!(bit_writemask[coarse_index] & 0x08))
                u[index+y_plus_one_z_plus_one_shift]=T();
            if(!(bit_writemask[coarse_index] & 0x10))
                u[index+x_shift]=T();
            if(!(bit_writemask[coarse_index] & 0x20))
                u[index+x_plus_one_z_plus_one_shift]=T();
            if(!(bit_writemask[coarse_index] & 0x40))
                u[index+x_plus_one_y_plus_one_shift]=T();
            if(!(bit_writemask[coarse_index] & 0x80))
                u[index+x_plus_one_y_plus_one_z_plus_one_shift]=T();
        }
}
//#####################################################################
template class Prolongation_Helper<float>;
template class Prolongation_Helper<double>;

