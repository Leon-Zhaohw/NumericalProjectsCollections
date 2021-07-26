//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis, Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
#ifndef __Initialize_Interior_Bitmaps_And_Diagonal_Entries_Helper__
#define __Initialize_Interior_Bitmaps_And_Diagonal_Entries_Helper__

namespace PhysBAM{

template<class T>
class Initialize_Interior_Bitmaps_And_Diagonal_Entries_Size_Specific_Helper_Base
{
protected:
    const unsigned char* const cell_type;
    unsigned char* const is_interior_bitmask;
    T* const diagonal;
    unsigned char* const has_full_diagonal_bitmask;

public:
    explicit Initialize_Interior_Bitmaps_And_Diagonal_Entries_Size_Specific_Helper_Base(const unsigned char* const cell_type_input,
        unsigned char* const is_interior_bitmask_input,T* const diagonal_input,unsigned char* const has_full_diagonal_bitmask_input)
        :cell_type(cell_type_input),is_interior_bitmask(is_interior_bitmask_input),diagonal(diagonal_input),
         has_full_diagonal_bitmask(has_full_diagonal_bitmask_input)
    {}

    virtual ~Initialize_Interior_Bitmaps_And_Diagonal_Entries_Size_Specific_Helper_Base() {}

//#####################################################################
    virtual void Run()=0;
    virtual void Run_Parallel(const int number_of_partitions)=0;
//#####################################################################
};

template<class T>
class Initialize_Interior_Bitmaps_And_Diagonal_Entries_Helper{
    const int x_size,y_size,z_size;
    Initialize_Interior_Bitmaps_And_Diagonal_Entries_Size_Specific_Helper_Base<T>* derived;

public:
    ~Initialize_Interior_Bitmaps_And_Diagonal_Entries_Helper() {delete derived;}
    void Run() {derived->Run();}
    void Run_Parallel(const int number_of_partitions) {derived->Run_Parallel(number_of_partitions);}

//#####################################################################
    Initialize_Interior_Bitmaps_And_Diagonal_Entries_Helper(const int x_size_input,const int y_size_input,const int z_size_input,
        const unsigned char* const cell_type_input,unsigned char* is_interior_bitmask_coarse_input,T* const diagonal_input,
        unsigned char* const has_full_diagonal_bitmask_input);
//#####################################################################
};

template<class T,int y_size,int z_size=y_size>
class Initialize_Interior_Bitmaps_And_Diagonal_Entries_Size_Specific_Helper:public Initialize_Interior_Bitmaps_And_Diagonal_Entries_Size_Specific_Helper_Base<T>
{
    typedef Initialize_Interior_Bitmaps_And_Diagonal_Entries_Size_Specific_Helper_Base<T> Base;
    using Base::cell_type;using Base::is_interior_bitmask;using Base::diagonal;using Base::has_full_diagonal_bitmask;

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
        x_plus_one_shift=x_shift,
        y_plus_one_shift=y_shift,
        z_plus_one_shift=z_shift,
        x_minus_one_shift=-x_shift,
        y_minus_one_shift=-y_shift,
        z_minus_one_shift=-z_shift
    };

public:
    explicit Initialize_Interior_Bitmaps_And_Diagonal_Entries_Size_Specific_Helper(const int x_size_input,const unsigned char* const cell_type_input,
        unsigned char* const is_interior_bitmask_input,T* const diagonal_input,unsigned char* const has_full_diagonal_bitmask_input)
        :Base(cell_type_input,is_interior_bitmask_input,diagonal_input,has_full_diagonal_bitmask_input)
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
