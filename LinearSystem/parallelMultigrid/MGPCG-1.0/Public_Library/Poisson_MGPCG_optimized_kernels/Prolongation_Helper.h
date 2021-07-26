//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis, Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
#ifndef __Prolongation_Helper__
#define __Prolongation_Helper__

namespace PhysBAM{

template<class T>
class Prolongation_Size_Specific_Helper_Base
{
protected:
    T* const u;
    const T* const u_coarse;
    const unsigned char* const bit_writemask;

public:
    explicit Prolongation_Size_Specific_Helper_Base(T* const u_input,const T* const u_coarse_input,const unsigned char* const bit_writemask_input)
        :u(u_input),u_coarse(u_coarse_input),bit_writemask(bit_writemask_input)
    {}

    virtual ~Prolongation_Size_Specific_Helper_Base() {}

//#####################################################################
    virtual void Run()=0;
    virtual void Run_Parallel(const int number_of_partitions)=0;
//#####################################################################
};

template<class T>
class Prolongation_Helper{
    const int x_size,y_size,z_size;
    Prolongation_Size_Specific_Helper_Base<T>* derived;

public:
    ~Prolongation_Helper() {delete derived;}
    void Run() {derived->Run();}
    void Run_Parallel(const int number_of_partitions) {derived->Run_Parallel(number_of_partitions);}

//#####################################################################
    Prolongation_Helper(const int x_size_input,const int y_size_input,const int z_size_input,T* const u_input,const T* const u_coarse_input,
        const unsigned char* const bit_writemask_input);
//#####################################################################
};

template<class T,int y_size,int z_size=y_size>
class Prolongation_Size_Specific_Helper:public Prolongation_Size_Specific_Helper_Base<T>
{
    typedef Prolongation_Size_Specific_Helper_Base<T> Base;
    using Base::u;using Base::u_coarse;using Base::bit_writemask;

    const int x_size,padded_x_size;
    const int coarse_x_size,coarse_padded_x_size;

    enum WORKAROUND{
        x_block_size=4,
        y_block_size=4,
        z_block_size=4,
        padded_y_size=y_size+2,
        padded_z_size=z_size+2,
        coarse_y_size=y_size/2,
        coarse_z_size=z_size/2,
        coarse_padded_y_size=coarse_y_size+2,
        coarse_padded_z_size=coarse_z_size+2,
	x_shift=padded_y_size*padded_z_size,
        y_shift=padded_z_size,
        z_shift=1,
        coarse_x_shift=coarse_padded_y_size*coarse_padded_z_size,
        coarse_y_shift=coarse_padded_z_size,
        coarse_z_shift=1,
        z_plus_one_shift=z_shift,
        y_plus_one_shift=y_shift,
        y_plus_one_z_plus_one_shift=y_shift+z_shift,
        x_plus_one_shift=x_shift,
        x_plus_one_z_plus_one_shift=x_shift+z_shift,
        x_plus_one_y_plus_one_shift=x_shift+y_shift,
        x_plus_one_y_plus_one_z_plus_one_shift=x_shift+y_shift+z_shift,
        coarse_x_minus_one_shift=-coarse_x_shift,
        coarse_x_minus_one_coarse_y_minus_one_shift=-coarse_x_shift-coarse_y_shift,
        coarse_x_minus_one_coarse_y_minus_one_coarse_z_minus_one_shift=-coarse_x_shift-coarse_y_shift-coarse_z_shift,
        coarse_x_minus_one_coarse_y_minus_one_coarse_z_plus_one_shift=-coarse_x_shift-coarse_y_shift+coarse_z_shift,
        coarse_x_minus_one_coarse_y_plus_one_shift=-coarse_x_shift+coarse_y_shift,
        coarse_x_minus_one_coarse_y_plus_one_coarse_z_minus_one_shift=-coarse_x_shift+coarse_y_shift-coarse_z_shift,
        coarse_x_minus_one_coarse_y_plus_one_coarse_z_plus_one_shift=-coarse_x_shift+coarse_y_shift+coarse_z_shift,
        coarse_x_minus_one_coarse_z_minus_one_shift=-coarse_x_shift-coarse_z_shift,
        coarse_x_minus_one_coarse_z_plus_one_shift=-coarse_x_shift+coarse_z_shift,
        coarse_x_plus_one_shift=coarse_x_shift,
        coarse_x_plus_one_coarse_y_minus_one_shift=coarse_x_shift-coarse_y_shift,
        coarse_x_plus_one_coarse_y_minus_one_coarse_z_minus_one_shift=coarse_x_shift-coarse_y_shift-coarse_z_shift,
        coarse_x_plus_one_coarse_y_minus_one_coarse_z_plus_one_shift=coarse_x_shift-coarse_y_shift+coarse_z_shift,
        coarse_x_plus_one_coarse_y_plus_one_shift=coarse_x_shift+coarse_y_shift,
        coarse_x_plus_one_coarse_y_plus_one_coarse_z_minus_one_shift=coarse_x_shift+coarse_y_shift-coarse_z_shift,
        coarse_x_plus_one_coarse_y_plus_one_coarse_z_plus_one_shift=coarse_x_shift+coarse_y_shift+coarse_z_shift,
        coarse_x_plus_one_coarse_z_minus_one_shift=coarse_x_shift-coarse_z_shift,
        coarse_x_plus_one_coarse_z_plus_one_shift=coarse_x_shift+coarse_z_shift,
        coarse_y_minus_one_shift=-coarse_y_shift,
        coarse_y_minus_one_coarse_z_minus_one_shift=-coarse_y_shift-coarse_z_shift,
        coarse_y_minus_one_coarse_z_plus_one_shift=-coarse_y_shift+coarse_z_shift,
        coarse_y_plus_one_shift=coarse_y_shift,
        coarse_y_plus_one_coarse_z_minus_one_shift=coarse_y_shift-coarse_z_shift,
        coarse_y_plus_one_coarse_z_plus_one_shift=coarse_y_shift+coarse_z_shift,
        coarse_z_minus_one_shift=-coarse_z_shift,
        coarse_z_plus_one_shift=coarse_z_shift
    };

public:
    explicit Prolongation_Size_Specific_Helper(const int x_size_input,T* const u_input,const T* const u_coarse_input,const unsigned char* const bit_writemask_input)
        :Base(u_input,u_coarse_input,bit_writemask_input)
	,x_size(x_size_input),padded_x_size(x_size_input+2)
	,coarse_x_size(x_size/2),coarse_padded_x_size(coarse_x_size+2)
    {}

    void Run()
    {Run_X_Range(1,coarse_x_size);}
    
    // For debugging purposes only

    static void Allocate_Data(T*& u,T*& u_coarse,unsigned char*& bit_writemask)
    {
	int padded_length = padded_y_size*padded_y_size*padded_z_size;
	int coarse_padded_length=coarse_padded_y_size*coarse_padded_y_size*coarse_padded_z_size;
	u=new T[padded_length];u_coarse=new T[coarse_padded_length];bit_writemask=new unsigned char[coarse_padded_length];}

    static void Initialize_Data(T* const u,T* const u_coarse,unsigned char* const bit_writemask)
    {
	int padded_length = padded_y_size*padded_y_size*padded_z_size;
	int coarse_padded_length=coarse_padded_y_size*coarse_padded_y_size*coarse_padded_z_size;
	for(int i=0;i<padded_length;i++) u[i]=(T)i;
	for(int i=0;i<coarse_padded_length;i++) u_coarse[i]=(T)i;
	for(int i=0;i<coarse_padded_length;i++) bit_writemask[i]=i%256;}

//#####################################################################
    void Run_Parallel(const int number_of_partitions);
    void Run_X_Range(const int xmin,const int xmax);
//#####################################################################
};
}
#endif
