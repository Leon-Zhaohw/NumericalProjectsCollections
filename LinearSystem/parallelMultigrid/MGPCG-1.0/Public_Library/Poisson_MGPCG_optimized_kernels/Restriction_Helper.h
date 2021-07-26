//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis, Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
#ifndef __Restriction_Helper__
#define __Restriction_Helper__

namespace PhysBAM{

template<class T>
class Restriction_Size_Specific_Helper_Base
{
protected:
    const T* const r;
    T* const b_coarse;
    const unsigned char* const cell_type_coarse;

public:
    explicit Restriction_Size_Specific_Helper_Base(const T* const r_input,T* const b_coarse_input,const unsigned char* const cell_type_coarse_input)
        :r(r_input),b_coarse(b_coarse_input),cell_type_coarse(cell_type_coarse_input)
    {}

    virtual ~Restriction_Size_Specific_Helper_Base() {}
    
//#####################################################################
    virtual void Run()=0;
    virtual void Run_Parallel(const int number_of_partitions)=0;
//#####################################################################
};

template<class T>
class Restriction_Helper{
    const int x_size,y_size,z_size;
    Restriction_Size_Specific_Helper_Base<T>* derived;

public:
    ~Restriction_Helper() {delete derived;}
    void Run() {derived->Run();}
    void Run_Parallel(const int number_of_partitions) {derived->Run_Parallel(number_of_partitions);}

//#####################################################################
    Restriction_Helper(const int x_size_input,const int y_size_input,const int z_size_input,const T* const r_input,T* const b_coarse_input,
        const unsigned char* const cell_type_coarse_input);
//#####################################################################
};

template<class T,int y_size,int z_size=y_size>
class Restriction_Size_Specific_Helper:public Restriction_Size_Specific_Helper_Base<T>
{
    typedef Restriction_Size_Specific_Helper_Base<T> Base;
    using Base::r;using Base::b_coarse;using Base::cell_type_coarse;

    const int x_size,padded_x_size;
    const int coarse_x_size,coarse_padded_x_size;

    enum WORKAROUND{
        INTERIOR_CELL_TYPE=1,
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
        x_minus_one_y_minus_one_z_minus_one_shift=-x_shift-y_shift-z_shift,
        x_minus_one_y_minus_one_shift=-x_shift-y_shift,
        x_minus_one_y_minus_one_z_plus_one_shift=-x_shift-y_shift+z_shift,
        x_minus_one_y_minus_one_z_plus_two_shift=-x_shift-y_shift+2*z_shift,
        x_minus_one_z_minus_one_shift=-x_shift-z_shift,
        x_minus_one_shift=-x_shift,
        x_minus_one_z_plus_one_shift=-x_shift+z_shift,
        x_minus_one_z_plus_two_shift=-x_shift+2*z_shift,
        x_minus_one_y_plus_one_z_minus_one_shift=-x_shift+y_shift-z_shift,
        x_minus_one_y_plus_one_shift=-x_shift+y_shift,
        x_minus_one_y_plus_one_z_plus_one_shift=-x_shift+y_shift+z_shift,
        x_minus_one_y_plus_one_z_plus_two_shift=-x_shift+y_shift+2*z_shift,
        x_minus_one_y_plus_two_z_minus_one_shift=-x_shift+2*y_shift-z_shift,
        x_minus_one_y_plus_two_shift=-x_shift+2*y_shift,
        x_minus_one_y_plus_two_z_plus_one_shift=-x_shift+2*y_shift+z_shift,
        x_minus_one_y_plus_two_z_plus_two_shift=-x_shift+2*y_shift+2*z_shift,
        y_minus_one_z_minus_one_shift=-y_shift-z_shift,
        y_minus_one_shift=-y_shift,
        y_minus_one_z_plus_one_shift=-y_shift+z_shift,
        y_minus_one_z_plus_two_shift=-y_shift+2*z_shift,
        z_minus_one_shift=-z_shift,
        z_plus_one_shift=z_shift,
        z_plus_two_shift=2*z_shift,
        y_plus_one_z_minus_one_shift=y_shift-z_shift,
        y_plus_one_shift=y_shift,
        y_plus_one_z_plus_one_shift=y_shift+z_shift,
        y_plus_one_z_plus_two_shift=y_shift+2*z_shift,
        y_plus_two_z_minus_one_shift=2*y_shift-z_shift,
        y_plus_two_shift=2*y_shift,
        y_plus_two_z_plus_one_shift=2*y_shift+z_shift,
        y_plus_two_z_plus_two_shift=2*y_shift+2*z_shift,
        x_plus_one_y_minus_one_z_minus_one_shift=x_shift-y_shift-z_shift,
        x_plus_one_y_minus_one_shift=x_shift-y_shift,
        x_plus_one_y_minus_one_z_plus_one_shift=x_shift-y_shift+z_shift,
        x_plus_one_y_minus_one_z_plus_two_shift=x_shift-y_shift+2*z_shift,
        x_plus_one_z_minus_one_shift=x_shift-z_shift,
        x_plus_one_shift=x_shift,
        x_plus_one_z_plus_one_shift=x_shift+z_shift,
        x_plus_one_z_plus_two_shift=x_shift+2*z_shift,
        x_plus_one_y_plus_one_z_minus_one_shift=x_shift+y_shift-z_shift,
        x_plus_one_y_plus_one_shift=x_shift+y_shift,
        x_plus_one_y_plus_one_z_plus_one_shift=x_shift+y_shift+z_shift,
        x_plus_one_y_plus_one_z_plus_two_shift=x_shift+y_shift+2*z_shift,
        x_plus_one_y_plus_two_z_minus_one_shift=x_shift+2*y_shift-z_shift,
        x_plus_one_y_plus_two_shift=x_shift+2*y_shift,
        x_plus_one_y_plus_two_z_plus_one_shift=x_shift+2*y_shift+z_shift,
        x_plus_one_y_plus_two_z_plus_two_shift=x_shift+2*y_shift+2*z_shift,
        x_plus_two_y_minus_one_z_minus_one_shift=2*x_shift-y_shift-z_shift,
        x_plus_two_y_minus_one_shift=2*x_shift-y_shift,
        x_plus_two_y_minus_one_z_plus_one_shift=2*x_shift-y_shift+z_shift,
        x_plus_two_y_minus_one_z_plus_two_shift=2*x_shift-y_shift+2*z_shift,
        x_plus_two_z_minus_one_shift=2*x_shift-z_shift,
        x_plus_two_shift=2*x_shift,
        x_plus_two_z_plus_one_shift=2*x_shift+z_shift,
        x_plus_two_z_plus_two_shift=2*x_shift+2*z_shift,
        x_plus_two_y_plus_one_z_minus_one_shift=2*x_shift+y_shift-z_shift,
        x_plus_two_y_plus_one_shift=2*x_shift+y_shift,
        x_plus_two_y_plus_one_z_plus_one_shift=2*x_shift+y_shift+z_shift,
        x_plus_two_y_plus_one_z_plus_two_shift=2*x_shift+y_shift+2*z_shift,
        x_plus_two_y_plus_two_z_minus_one_shift=2*x_shift+2*y_shift-z_shift,
        x_plus_two_y_plus_two_shift=2*x_shift+2*y_shift,
        x_plus_two_y_plus_two_z_plus_one_shift=2*x_shift+2*y_shift+z_shift,
        x_plus_two_y_plus_two_z_plus_two_shift=2*x_shift+2*y_shift+2*z_shift
    };

public:
    explicit Restriction_Size_Specific_Helper(const int x_size_input,const T* const r_input,T* const b_coarse_input,const unsigned char* const cell_type_coarse_input)
        :Base(r_input,b_coarse_input,cell_type_coarse_input)
	,x_size(x_size_input),padded_x_size(x_size_input+2)
	,coarse_x_size(x_size/2),coarse_padded_x_size(coarse_x_size+2)
    {}

    void Run()
    {Run_X_Range(1,coarse_x_size);}
    
    // For debugging purposes only

    static void Allocate_Data(T*& r,T*& b_coarse,unsigned char*& cell_type_coarse)
    {
	int padded_length = padded_y_size*padded_y_size*padded_z_size;
	int coarse_padded_length=coarse_padded_y_size*coarse_padded_y_size*coarse_padded_z_size;
	r=new T[padded_length];b_coarse=new T[coarse_padded_length];cell_type_coarse=new unsigned char[coarse_padded_length];}

    static void Initialize_Data(T* const r,T* const b_coarse,unsigned char* const cell_type_coarse)
    {
	int padded_length = padded_y_size*padded_y_size*padded_z_size;
	int coarse_padded_length=coarse_padded_y_size*coarse_padded_y_size*coarse_padded_z_size;
	for(int i=0;i<padded_length;i++) r[i]=(T)i;
	for(int i=0;i<coarse_padded_length;i++) b_coarse[i]=(T)i;
	for(int i=0;i<coarse_padded_length;i++) cell_type_coarse[i]=i%256;}

//#####################################################################
    void Run_Parallel(const int number_of_partitions);
    void Run_X_Range(const int xmin,const int xmax);
//#####################################################################
};
}
#endif
