//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis,Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################

#include <Matrices_And_Vectors/VECTOR.h>
#include <Read_Write/FILE_UTILITIES.h>
#include <Utilities/LOG.h>
#include "MULTIGRID_POISSON_TESTS.h"
#include <Poisson_MGPCG_optimized_kernels/PTHREAD_QUEUE.h>


using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;
int main(int argc,char* argv[])
{ 
    typedef float T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
    static const int d=3;

#ifdef MGPCG_UNOPTIMIZED
    STATIC_ASSERT((d==2)||(d==3));
#else
    STATIC_ASSERT(d==3);
#endif

    if(argc<3 || argc >4){
	std::cout<<"Usage: "<<argv[0]<<" <test_number> <number_of_threads> <resolution>(optional)"<<std::endl; return 1;
    }

    const int test_number=atoi(argv[1]);
    const int number_of_threads=atoi(argv[2]);
    int resolution =256;
    if(argc==4)
	resolution=atoi(argv[3]);


    pthread_queue=new PhysBAM::PTHREAD_QUEUE(number_of_threads);
    MULTIGRID_POISSON_TESTS<T,d> tests(test_number,number_of_threads,resolution);
    tests.Run();
    return 0;
}
