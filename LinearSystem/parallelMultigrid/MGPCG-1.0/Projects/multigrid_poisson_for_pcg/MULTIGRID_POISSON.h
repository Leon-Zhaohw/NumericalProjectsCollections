//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis,Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
// Class MULTIGRID_POISSON
//#####################################################################
#ifndef __MULTIGRID_POISSON__
#define __MULTIGRID_POISSON__

#include <Arrays/ARRAYS_2D.h>
#include <Arrays/ARRAYS_3D.h>
#include <Grids/POLICY_UNIFORM.h>
#include <Matrices_And_Vectors/VECTOR_2D.h>
#include <Matrices_And_Vectors/VECTOR_3D.h>
#include <Utilities/NONCOPYABLE.h>
#include <Utilities/LOG.h>

namespace PhysBAM{

template<class T,int d>
class MULTIGRID_POISSON:public NONCOPYABLE
{
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    typedef typename GRID_POLICY<TV>::UNIFORM_GRID T_GRID;
public:
    typedef enum {INTERIOR_CELL_TYPE=1,DIRICHLET_CELL_TYPE=2,NEUMANN_CELL_TYPE=3} CELL_TYPE;
private:
    typedef typename POLICY_UNIFORM<VECTOR<CELL_TYPE,d> >::ARRAYS_SCALAR T_CELL_TYPE_FIELD;
    typedef typename POLICY_UNIFORM<VECTOR<unsigned char,d> >::ARRAYS_SCALAR T_CELL_TYPE_AS_UNSIGNED_CHAR;
    typedef typename POLICY_UNIFORM<TV>::ARRAYS_SCALAR T_VARIABLE;
    typedef typename POLICY_UNIFORM<VECTOR<bool,d> >::ARRAYS_SCALAR T_FLAG;
    typedef typename POLICY_UNIFORM<VECTOR<unsigned char,d> >::ARRAYS_SCALAR T_BITMASK;

public:

    const T h;
    const T_INDEX n;
    const T_GRID grid;
    const T_GRID coarse_grid;

    T_CELL_TYPE_AS_UNSIGNED_CHAR cell_type;
    T_VARIABLE u;
    T_VARIABLE b;

    T_VARIABLE delta;
    LIST_ARRAY<T> one_over_diagonal_part; //defined only on boundary_indices
    T_VARIABLE diagonal_entries;

    int total_red_boundary_blocks,total_black_boundary_blocks;

    const int colors;
    const int& number_of_threads;

    // NOTE: in optimized version, boundary indices stored in flattened format
    // for unoptimized version, boundary indices stored in vector format
#ifndef MGPCG_UNOPTIMIZED

    ARRAY<int> boundary_block_start;
    ARRAY<int> boundary_block_end;
    ARRAY<int> boundary_indices;

    ARRAY<int> extended_boundary_indices;
#else
    LIST_ARRAY<int> boundary_block_start;
    LIST_ARRAY<int> boundary_block_end;
    LIST_ARRAY<T_INDEX> boundary_indices;

    LIST_ARRAY<T_INDEX> extended_boundary_indices;

#endif

    // for serial iterators
    BOX<T_INDEX> padded_domain;
    BOX<T_INDEX> unpadded_domain;
    BOX<T_INDEX> padded_coarse_domain;
    BOX<T_INDEX> unpadded_coarse_domain;

    enum WORKAROUND1 {boundary_block_size=4};
    enum WORKAROUND2 {boundary_block_padding=0};

    // bitmasks defined over 2x2x2 blocks 
    T_BITMASK index_has_full_diagonal_coarse_bitmask;
    T_BITMASK index_is_interior_coarse_bitmask;


//#####################################################################
public:
    MULTIGRID_POISSON(const T_INDEX& n,const T h,const int& number_of_threads_input);
    void Initialize();
    void Reinitialize();
    T Apply_System_Matrix(const T_INDEX& index);
    T Apply_System_Matrix(const T_INDEX& index,const T_VARIABLE& u_input) const;
    void Subtract_Multiple_Of_System_Matrix_And_Compute_Sum_And_Extrema(const T_VARIABLE& x,T_VARIABLE& y,double& sum,T& rmin,T& rmax) const;
    T Multiply_With_System_Matrix_And_Compute_Dot_Product(const T_VARIABLE&x,T_VARIABLE& y) const;

private:
    void Initialize_Boundary_Region();
    void Initialize_Boundary_Blocks();

    void Initialize_Interior_Bitmaps_And_Diagonal_Entries();

    void Build_System_Matrix();
    void Boundary_Relaxation(bool reverse_order,const int loops=1);
    T Interior_Relaxation(const int loops=1,const bool compute_dot_product=false);
    void Compute_Residuals_Boundary();

public:
    T Relaxation_Sweep(const bool reverse_order,const int boundary_pre_relaxations=2,const int interior_relaxations=1,const int boundary_post_relaxations=2,const bool compute_dot_product=false,const bool verbose=true);
    void Relax_And_Compute_Residuals(const int interior_relaxations,const int boundary_post_relaxations,const T nullspace_component);


    // Initialization helpers for debugging example
    void Initialize_Square_Domain();
    void Initialize_Square_Minus_Circle_Domain();
    void Initialize_Test_Domain();
    void Initialize_Test_Right_Hand_Side(T_VARIABLE& b_input);
    void Initialize_Test_Right_Hand_Side();
    void Initialize_Test_Initial_Guess(T_VARIABLE& u_input);
    void Initialize_Test_Initial_Guess();
//#####################################################################
};
}
#endif
