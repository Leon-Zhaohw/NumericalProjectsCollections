//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis,Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
#ifndef __Coarsened_Discretization_Helper__
#define __Coarsened_Discretization_Helper__

namespace PhysBAM{

template<class T>
class Coarsened_Discretization_Size_Specific_Helper_Base
{
protected:
    const unsigned char* const fine_cell_type;
    unsigned char* const coarse_cell_type;

public:
    explicit Coarsened_Discretization_Size_Specific_Helper_Base(const unsigned char* const fine_cell_type_input, unsigned char* const coarse_cell_type_input)
	: fine_cell_type(fine_cell_type_input),coarse_cell_type(coarse_cell_type_input)
    {}

    virtual ~Coarsened_Discretization_Size_Specific_Helper_Base() {}

//#####################################################################
    virtual void Run()=0;
    virtual void Run_Parallel(const int number_of_partitions)=0;
//#####################################################################
};

template<class T>
class Coarsened_Discretization_Helper{
    const int x_size,y_size,z_size;
    Coarsened_Discretization_Size_Specific_Helper_Base<T>* derived;

public:
    ~Coarsened_Discretization_Helper() {delete derived;}
    void Run() {derived->Run();}
    void Run_Parallel(const int number_of_partitions) {derived->Run_Parallel(number_of_partitions);}

//#####################################################################
    Coarsened_Discretization_Helper(const int x_size_input,const int y_size_input,const int z_size_input,const unsigned char* const fine_cell_type_input, unsigned char* const coarse_cell_type_input);
//#####################################################################
};

template<class T,int y_size,int z_size=y_size>
class Coarsened_Discretization_Size_Specific_Helper:public Coarsened_Discretization_Size_Specific_Helper_Base<T>
{
    typedef Coarsened_Discretization_Size_Specific_Helper_Base<T> Base;
    using Base::fine_cell_type; using Base::coarse_cell_type;

    const int x_size,padded_x_size;
    const int coarse_x_size,coarse_padded_x_size;
    
    enum WORKAROUND{
	INTERIOR_CELL_TYPE=1,
	DIRICHLET_CELL_TYPE=2,
	NEUMANN_CELL_TYPE=3,
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
        x_plus_one_y_plus_one_z_plus_one_shift=x_shift+y_shift+z_shift,
        x_plus_one_y_plus_one_shift=x_shift+y_shift,
        x_plus_one_z_plus_one_shift=x_shift+z_shift,
        x_plus_one_shift=x_shift,
        y_plus_one_shift=y_shift,
        z_plus_one_shift=z_shift,
        y_plus_one_z_plus_one_shift=y_shift+z_shift
    };

public:
    explicit Coarsened_Discretization_Size_Specific_Helper(const int x_size_input,const unsigned char* const fine_cell_type_input, unsigned char* const coarse_cell_type_input)
        :Base(fine_cell_type_input,coarse_cell_type_input)
	,x_size(x_size_input),padded_x_size(x_size_input+2)
	,coarse_x_size(x_size/2),coarse_padded_x_size(coarse_x_size+2)
    {}

    void Run()
    {Run_X_Range(1,coarse_x_size);}
    
    // For debugging purposes only

//     static void Allocate_Data(T*& u,T*& b,T*& r,unsigned char*& bit_writemask)
//     {u=new T[padded_length];b=new T[padded_length];r=new T[padded_length];bit_writemask=new unsigned char[coarse_padded_length];}

//     static void Initialize_Data(T* const u,T* const b,T* const r,unsigned char* const bit_writemask)
//     {for(int i=0;i<padded_length;i++) u[i]=b[i]=r[i]=(T)i;for(int i=0;i<coarse_padded_length;i++) bit_writemask[i]=0xff;}

//#####################################################################
    void Run_Parallel(const int number_of_partitions);
    void Run_X_Range(const int xmin,const int xmax);
//#####################################################################
};
}
#endif
