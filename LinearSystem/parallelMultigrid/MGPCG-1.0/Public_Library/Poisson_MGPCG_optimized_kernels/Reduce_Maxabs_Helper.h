//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis, Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
#ifndef __Reduce_Maxabs_Helper__
#define __Reduce_Maxabs_Helper__
namespace PhysBAM{

template<class T>
class Reduce_Maxabs_Size_Specific_Helper_Base
{
protected:
    const T* const u;
    T* maxabs_partial_results;

public:
    explicit Reduce_Maxabs_Size_Specific_Helper_Base(const T* const u_input)
	:u(u_input),maxabs_partial_results(0)
    {}
    
    virtual ~Reduce_Maxabs_Size_Specific_Helper_Base() {}

//#####################################################################
    virtual T Run()=0;
    virtual T Run_Parallel(const int number_of_partitions)=0;
//#####################################################################
};

template<class T>
class Reduce_Maxabs_Helper{
    const int x_size,y_size,z_size;
    Reduce_Maxabs_Size_Specific_Helper_Base<T>* derived;

public:
    ~Reduce_Maxabs_Helper() {delete derived;}
    T Run() {return derived->Run();}
    T Run_Parallel(const int number_of_partitions) {return derived->Run_Parallel(number_of_partitions);}

//#####################################################################
    Reduce_Maxabs_Helper(const int x_size_input,const int y_size_input,const int z_size_input,const T* const u);
//#####################################################################
};

template<class T,int y_size,int z_size=y_size>
class Reduce_Maxabs_Size_Specific_Helper:public Reduce_Maxabs_Size_Specific_Helper_Base<T>
{
    typedef Reduce_Maxabs_Size_Specific_Helper_Base<T> Base;
    using Base::u;using Base::maxabs_partial_results;
    
    const int x_size,padded_x_size;
    enum WORKAROUND{
        x_block_size=4,
        y_block_size=4,
        z_block_size=4,
        padded_y_size=y_size+2,
        padded_z_size=z_size+2,
        x_shift=padded_y_size*padded_z_size,
        y_shift=padded_z_size,
        z_shift=1
    };

public:
    explicit Reduce_Maxabs_Size_Specific_Helper(const int x_size_input,const T* const u_input)
	:Base(u_input),x_size(x_size_input),padded_x_size(x_size_input+2)
    {}
    
    T Run()
    {maxabs_partial_results=new T[1];
    Run_X_Range(1,x_size,0);
    T maxabs=maxabs_partial_results[0];
    delete[] maxabs_partial_results;maxabs_partial_results=0;
    return maxabs;}
    
    // For debugging purposes only

    static void Allocate_Data(T*& x)
    {int padded_length=padded_y_size*padded_y_size*padded_z_size;
    x=new T[padded_length];}

    static void Initialize_Data(T* const x)
    {int padded_length=padded_y_size*padded_y_size*padded_z_size;
    for(int i=0;i<padded_length;i++) x[i]=(T)i;}

//#####################################################################
    T Run_Parallel(const int number_of_partitions);
    void Run_X_Range(const int xmin,const int xmax,const int partition_number);
//#####################################################################
};
}
#endif
