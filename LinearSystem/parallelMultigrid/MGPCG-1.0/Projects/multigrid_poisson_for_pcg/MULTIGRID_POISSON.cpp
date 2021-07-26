//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis,Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
// Class MULTIGRID_POISSON
//#####################################################################
#include <Utilities/LOG.h>
#include "MULTIGRID_POISSON.h"
#include "BOX_ITERATOR.h"

#ifndef MGPCG_UNOPTIMIZED
#include <Poisson_MGPCG_optimized_kernels/Relaxation_And_Residual_With_Zero_Initial_Guess_Helper.h>
#include <Poisson_MGPCG_optimized_kernels/Relaxation_Boundary_Helper.h>
#include <Poisson_MGPCG_optimized_kernels/Residual_Boundary_Helper.h>
#include <Poisson_MGPCG_optimized_kernels/Relaxation_Interior_Helper.h>
#include <Poisson_MGPCG_optimized_kernels/Multiplication_And_Dot_Product_Helper.h>
#include <Poisson_MGPCG_optimized_kernels/Multiply_And_Compute_Sum_And_Extrema_Helper.h>
#include <Poisson_MGPCG_optimized_kernels/Relaxation_And_Dot_Product_Interior_Helper.h>

#include <Poisson_MGPCG_optimized_kernels/Initialize_Interior_Bitmaps_And_Diagonal_Entries_Helper.h>
#include <Poisson_MGPCG_optimized_kernels/Boundary_Initialization_Helper.h>
#include <Poisson_MGPCG_optimized_kernels/Block_Enumeration_Helper.h>
#include <Poisson_MGPCG_optimized_kernels/Block_Counting_Helper.h>
#endif
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> MULTIGRID_POISSON<T,d>::
MULTIGRID_POISSON(const T_INDEX& n_input,const T h_input,const int& number_of_threads_input)
    :h(h_input),n(n_input),grid(n+2,BOX<TV>(-TV::All_Ones_Vector()*.5*h,TV(n)*h+TV::All_Ones_Vector()*.5*h)),
    coarse_grid(n/2+2,BOX<TV>(-TV::All_Ones_Vector()*h,TV(n)*h+TV::All_Ones_Vector()*h)),
    colors(2),number_of_threads(number_of_threads_input)
    ,padded_domain(T_INDEX::All_Ones_Vector(),n_input+2),unpadded_domain(2*T_INDEX::All_Ones_Vector(),n_input+1)
    ,padded_coarse_domain(T_INDEX::All_Ones_Vector(),n_input/2+2),unpadded_coarse_domain(2*T_INDEX::All_Ones_Vector(),n_input/2+1)
    
{
    cell_type.Resize(grid,0,false,false);
    u.Resize(grid,0,true,false);
    b.Resize(grid,0,true,false);
    delta.Resize(grid,0,true,false);

}
//#####################################################################
// Function Initialize 
//#####################################################################
template<class T,int d> void MULTIGRID_POISSON<T,d>::
Initialize()
{
    // Initialize state variables
    index_has_full_diagonal_coarse_bitmask.Resize(coarse_grid,0,false,false);
    index_is_interior_coarse_bitmask.Resize(coarse_grid,0,false,false);
    diagonal_entries.Resize(grid,0,false,false);

    Initialize_Interior_Bitmaps_And_Diagonal_Entries();

    // Initialize boundary-related structures
    Initialize_Boundary_Region();

    // Initialize discretization
    Build_System_Matrix();
}

//#####################################################################
// Function Reinitialize 
//#####################################################################
template<class T,int d> void MULTIGRID_POISSON<T,d>::
Reinitialize()
{
    // Reinitialize state variables
    u.Fill(0);

    one_over_diagonal_part.Remove_All();
    Initialize_Interior_Bitmaps_And_Diagonal_Entries();

    // Initialize boundary-related structures
    Initialize_Boundary_Region();

    // Initialize discretization
    Build_System_Matrix();
}


//#####################################################################
// Function Initialize_Boundary_Region
//#####################################################################
template<class T,int d> void MULTIGRID_POISSON<T,d>::
Initialize_Boundary_Region(){

    LOG::SCOPE scope("Initialize_Boundary_Region","Initializing boundary region");

#ifndef MGPCG_UNOPTIMIZED
    LOG::Time("Allocate temporaries");
    T_BITMASK index_is_boundary(grid);
    T_BITMASK index_is_extended_boundary(grid);
    for(int v=1;v<=d;v++) PHYSBAM_ASSERT(n(v)%boundary_block_size==0);
    ARRAY<int> boundary_nodes_in_block(n(1)*n(2)*n(3)/(boundary_block_size*boundary_block_size*boundary_block_size),false);
    ARRAY<int> extended_boundary_nodes_in_block(n(1)*n(2)*n(3)/(boundary_block_size*boundary_block_size*boundary_block_size),false);

    LOG::Time("Boundary_Initialization");
    Boundary_Initialization_Helper<void> boundary_initialization(n(1),n(2),n(3),&cell_type(1,1,1),&index_is_interior_coarse_bitmask(1,1,1),&index_is_boundary(1,1,1),&index_is_extended_boundary(1,1,1));
    boundary_initialization.Run_Parallel(number_of_threads);

    LOG::Time("Block_Counting"); 
    int total_extended_boundary_blocks=0,total_red_boundary_indices=0,total_black_boundary_indices=0,total_extended_boundary_indices=0;
    total_red_boundary_indices=0;total_black_boundary_indices=0;

    Block_Counting_Helper<void> block_counting(n(1),n(2),n(3),&index_is_boundary(1,1,1),&index_is_extended_boundary(1,1,1),&boundary_nodes_in_block(1),&extended_boundary_nodes_in_block(1),total_red_boundary_blocks,total_black_boundary_blocks,total_extended_boundary_blocks,total_red_boundary_indices,total_black_boundary_indices,total_extended_boundary_indices);
    block_counting.Run_Parallel(number_of_threads);

    LOG::cout<<"Found "<<total_red_boundary_indices<<" red boundary indices"<<std::endl;
    LOG::cout<<"Found "<<total_black_boundary_indices<<" black boundary indices"<<std::endl;

    LOG::Time("Allocate boundary block arrays");
    boundary_block_start.Resize(total_red_boundary_blocks+total_black_boundary_blocks,false,false);
    boundary_block_end.Resize(boundary_block_start.m,false,false);
    boundary_indices.Resize(total_red_boundary_indices+total_black_boundary_indices,false,false);
    extended_boundary_indices.Resize(total_extended_boundary_indices);

    LOG::Time("Block Enumeration");
    Block_Enumeration_Helper<void> block_enumeration(n(1),n(2),n(3),&index_is_boundary(1,1,1),&index_is_extended_boundary(1,1,1),&boundary_nodes_in_block(1),&extended_boundary_nodes_in_block(1),&boundary_block_start(1),&boundary_block_end(1),&boundary_block_start(total_red_boundary_blocks+1),&boundary_block_end(total_red_boundary_blocks+1),&boundary_indices(1),&boundary_indices(total_red_boundary_indices+1),total_red_boundary_indices,&extended_boundary_indices(1));
    block_enumeration.Run_Parallel(number_of_threads);


#else

    LOG::Time("Allocate Temporaries");    
    for(int v=1;v<=d;v++) PHYSBAM_ASSERT(n(v)%2==0);
    T_FLAG index_is_boundary(grid);
    T_FLAG index_is_extended_boundary(grid);

    index_is_boundary.Fill(0);
    index_is_extended_boundary.Fill(0);
    unsigned char full_interior_coarse_cell; // bitmask for fully interior coarse cell
    if(d==2)
	full_interior_coarse_cell=0x0f;
    else if(d==3)
	full_interior_coarse_cell=0xff;

    LOG::Time("Boundary Initialization");
    // any interior cell along the unpadded domain is boundary
    for(BOUNDARY_ITERATOR<d> iterator(unpadded_domain);iterator.Valid();iterator.Next()){
	const T_INDEX& index=iterator.Index();
	if(cell_type(index)==INTERIOR_CELL_TYPE)
	    index_is_boundary(index)=true;
    }

    for(BOX_ITERATOR<d> coarse_iterator(unpadded_coarse_domain);coarse_iterator.Valid();coarse_iterator.Next()){
        const T_INDEX& coarse_index=coarse_iterator.Index();
        const T_INDEX base_fine_index=(coarse_index-1)*2;

	// check if any 2x2(x2) children of the coarse cell are not interior cells
	// if any are not interior, mark all fine cells with coarse cell 
	// in prolongation stencil as boundary
        if(index_is_interior_coarse_bitmask(coarse_index)!=full_interior_coarse_cell)
            for(BOX_ITERATOR<d> fine_iterator(BOX<T_INDEX>(base_fine_index-1,base_fine_index+2));fine_iterator.Valid();fine_iterator.Next()){
		const T_INDEX& fine_index=fine_iterator.Index();
		if(cell_type(fine_index)==INTERIOR_CELL_TYPE)
		    index_is_boundary(fine_index)=true;}
    }

    // find extended boundary cells
    // an interior cell is an extended boundary cell if 
    // it or any of its 6 face neighbors are boundary cells
    for(BOX_ITERATOR<d> iterator(unpadded_domain);iterator.Valid();iterator.Next()){
	const T_INDEX& index=iterator.Index();
	if(cell_type(index)!=INTERIOR_CELL_TYPE) continue;
	if(index_is_boundary(index)){ 
	    index_is_extended_boundary(index)=true; 
	    continue;
	}
	for(int v=1;v<=d;v++){
	    if(index_is_boundary(index+T_INDEX::Axis_Vector(v))){
		index_is_extended_boundary(index)=true; 
		break;}
	    if(index_is_boundary(index-T_INDEX::Axis_Vector(v))){
		index_is_extended_boundary(index)=true; 
		break;}
	}
    }
    LOG::Time("Block Counting");
    for(int v=1;v<=d;v++) PHYSBAM_ASSERT(n(v)%boundary_block_size==0);
    ARRAY<LIST_ARRAY<T_INDEX> > boundary_block_base_index; 
    boundary_block_base_index.Resize(2);
    LIST_ARRAY<T_INDEX> extended_boundary_block_base_index;

    // find all blocks with boundary and extended boundary nodes.
    // store index of the first node in each block.
    for(BOX_ITERATOR<d,boundary_block_size> domain_iterator(unpadded_domain);domain_iterator.Valid();domain_iterator.Next()){
	const T_INDEX& base_index=domain_iterator.Index();
	BOX<T_INDEX> box=BOX<T_INDEX>(base_index-boundary_block_padding,base_index+boundary_block_size-1+boundary_block_padding);

	bool contains_boundary=false;
	bool contains_extended_boundary=false;
	for(BOX_ITERATOR<d> box_iterator(box);box_iterator.Valid();box_iterator.Next())
	    if(index_is_boundary(box_iterator.Index())){
		contains_boundary=true;
		contains_extended_boundary=true;
 		break;
	    }else if(index_is_extended_boundary(box_iterator.Index())){
		contains_extended_boundary=true;
	    }
	if(contains_boundary){
	    T_INDEX block_index=(base_index-2)/boundary_block_size;
	    int box_color;
	    if(d==2)
		box_color = (block_index(1)+block_index(2))%2+1;
	    else if (d==3)
		box_color = (block_index(1)+block_index(2)+block_index(3))%2+1;
	    boundary_block_base_index(box_color).Append(base_index);
	}
	if(contains_extended_boundary)
	    extended_boundary_block_base_index.Append(base_index);
    }

    total_red_boundary_blocks=boundary_block_base_index(1).m;
    total_black_boundary_blocks=boundary_block_base_index(2).m;
    
    boundary_block_start.Clean_Memory();
    boundary_block_end.Clean_Memory();
    boundary_indices.Clean_Memory();
    extended_boundary_indices.Clean_Memory();
    boundary_block_start.Preallocate(boundary_block_base_index(1).m+boundary_block_base_index(2).m);
    boundary_block_end.Preallocate(boundary_block_base_index(1).m+boundary_block_base_index(2).m);  
    
    
    // NOTE: boundary_block_start in 
    // optimized and unoptimized code will be off by 1
    LOG::Time("Block Enumeration");
    for(int c=1;c<=2;c++){
	for(int b=1;b<=boundary_block_base_index(c).m;b++){ // for each block with boundary indices
            const T_INDEX& offset=boundary_block_base_index(c)(b);
 	    boundary_block_start.Append(boundary_indices.m+1); // first boundary_index of the block (1-indexed)

	    for(BOX_ITERATOR<d> iterator(BOX<T_INDEX>(offset,offset+T_INDEX::All_Ones_Vector()*(boundary_block_size-1)));
		iterator.Valid();iterator.Next()){
		const T_INDEX& index=iterator.Index();
		if(index_is_boundary(index)){
 		    boundary_indices.Append(index);
		}
	    }
 	    boundary_block_end.Append(boundary_indices.m); // last boundary_index of the block (1-indexed)

	}
    }
    // by going through the blocks, 
    // we enumerate extended boundary indices
    // in a block-by-block order (for fast reading)
    for(int b=1;b<=extended_boundary_block_base_index.m;b++){
	const T_INDEX& offset=extended_boundary_block_base_index(b);
	
	for(BOX_ITERATOR<d> iterator(BOX<T_INDEX>(offset,offset+T_INDEX::All_Ones_Vector()*(boundary_block_size-1)));
	    iterator.Valid();iterator.Next()){
	    const T_INDEX& index=iterator.Index();
	    if(index_is_extended_boundary(index))
		extended_boundary_indices.Append(index);
	}
    }
    LOG::cout<<"Found "<<boundary_block_end(total_red_boundary_blocks)<<" red boundary indices"<<std::endl;
 #endif
    LOG::cout<<"Found "<<boundary_indices.m<<" boundary indices"<<std::endl;
    LOG::cout<<"Found "<<total_red_boundary_blocks<<" red boundary blocks"<<std::endl;
    LOG::cout<<"Found "<<total_black_boundary_blocks<<" black boundary blocks"<<std::endl;
    LOG::cout<<"Found "<<extended_boundary_indices.m<<" extended boundary indices"<<std::endl;
	 
}


//#####################################################################
// Function Build_System_Matrix
//#####################################################################
template<class T,int d> void MULTIGRID_POISSON<T,d>::
Build_System_Matrix(){
    LOG::SCOPE scope("Build_System_Matrix","Building system matrix");

    LOG::Time("Compute diagonal entries");
    // compute one_over_diagonal_part only on boundary_indices
#ifndef MGPCG_UNOPTIMIZED
    const T* const diagonal_entries_ptr=&diagonal_entries(1,1,1);
    for(int i=1;i<=boundary_indices.m;i++){
	int index=boundary_indices(i);
	if(diagonal_entries_ptr[index]!=0) 
	    one_over_diagonal_part.Append((T)1/diagonal_entries_ptr[index]);
	else one_over_diagonal_part.Append((T)0);
    }

#else
    for(int i=1;i<=boundary_indices.m;i++){
	const T_INDEX& index=boundary_indices(i);
	if(diagonal_entries(index)!=0) 
	    one_over_diagonal_part.Append((T)1/diagonal_entries(index));
	else
	    one_over_diagonal_part.Append((T)0);
    }

#endif
}

//#####################################################################
// Function Initialize_Interior_Bitmaps_And_Diagonal_Entries
//#####################################################################
template<class T,int d> void MULTIGRID_POISSON<T,d>::
Initialize_Interior_Bitmaps_And_Diagonal_Entries(){
    LOG::SCOPE scope("Initialize_Interior_Bitmaps");

    LOG::Time("Zeroing out fields");

    for(BOUNDARY_ITERATOR<d> iterator(padded_coarse_domain);iterator.Valid();iterator.Next()){
	const T_INDEX& index=iterator.Index();
	index_is_interior_coarse_bitmask(index)=0;
	index_has_full_diagonal_coarse_bitmask(index)=0;
    }

    for(BOUNDARY_ITERATOR<d> iterator(padded_domain);iterator.Valid();iterator.Next())
	diagonal_entries(iterator.Index())=0;

#ifndef MGPCG_UNOPTIMIZED
    LOG::Time("Computing bitmaps and diagonal");
    Initialize_Interior_Bitmaps_And_Diagonal_Entries_Helper<T> initialize_interior_bitmaps_and_diagonal_entries(n(1),n(2),n(3),&cell_type(1,1,1),&index_is_interior_coarse_bitmask(1,1,1),&diagonal_entries(1,1,1),&index_has_full_diagonal_coarse_bitmask(1,1,1));
    initialize_interior_bitmaps_and_diagonal_entries.Run_Parallel(number_of_threads);

#else

    // we're initializing bitmasks
    // so loop over 2x2(x2) blocks

    int full_diagonal=d*2;
    
    for(BOX_ITERATOR<d> coarse_iterator(unpadded_coarse_domain);coarse_iterator.Valid();coarse_iterator.Next()){
        const T_INDEX& coarse_index=coarse_iterator.Index();
        const T_INDEX base_fine_index=(coarse_index-1)*2;
	
	index_is_interior_coarse_bitmask(coarse_index)=0;
	index_has_full_diagonal_coarse_bitmask(coarse_index)=0;
	
	unsigned char mask=0x01;
	for(BOX_ITERATOR<d> fine_iterator(BOX<T_INDEX>(base_fine_index,base_fine_index+1));fine_iterator.Valid();fine_iterator.Next()){
	    const T_INDEX& fine_index=fine_iterator.Index();
	    int active_neighbors=0;
	    if(cell_type(fine_index)==INTERIOR_CELL_TYPE){
		index_is_interior_coarse_bitmask(coarse_index)|=mask;
		for(int v=1;v<=d;v++){
		    if(cell_type(fine_index+T_INDEX::Axis_Vector(v))!=NEUMANN_CELL_TYPE)
			active_neighbors++;
		    if(cell_type(fine_index-T_INDEX::Axis_Vector(v))!=NEUMANN_CELL_TYPE)
			active_neighbors++;
		}
	    }
	    diagonal_entries(fine_index)=-(T)active_neighbors;
	    if(active_neighbors==full_diagonal) 
		index_has_full_diagonal_coarse_bitmask(coarse_index)|=mask;

	    mask<<=1;
	}
    }


#endif
    
}
//#####################################################################
// Function Boundary_Relaxation
//#####################################################################
template<class T,int d> void MULTIGRID_POISSON<T,d>::
Boundary_Relaxation(const bool reverse_order,const int subloops)
{
    if (boundary_indices.m==0) return;
#ifndef MGPCG_UNOPTIMIZED
    Relaxation_Boundary_Helper<T> boundary_relax(n(1),n(2),n(3),&u(1,1,1),&b(1,1,1),&one_over_diagonal_part(1),&boundary_indices(1),&boundary_block_start(1),&boundary_block_end(1),total_red_boundary_blocks,total_black_boundary_blocks,10,reverse_order);
    for(int loop=0;loop<subloops;loop++)
 	boundary_relax.Run_Parallel(number_of_threads);

#else
    // NOTE: we can ignore red-black partitioning since we are running serially
    // when running in parallel, must process the first 
    // #(total_red_boundary_blocks) boundary blocks separately
    // from the last #(total_black_boundary_blocks) boundary blocks
    for(int loop=0;loop<subloops;loop++){
	if(reverse_order)
	    for(int block=boundary_block_start.m;block>=1;block--){
		int first_index=boundary_block_start(block);
		int last_index=boundary_block_end(block);	
		for(int block_loop=1;block_loop<=10;block_loop++)
		    for(int i=last_index;i>=first_index;i--){
			const T_INDEX& index=boundary_indices(i);
			u(index)=b(index);
			for(int v=1;v<=d;v++)
			    u(index)-=u(index+T_INDEX::Axis_Vector(v))
				+u(index-T_INDEX::Axis_Vector(v));
			u(index)*=one_over_diagonal_part(i);
		    }
	    }


	else
	    for(int block=1;block<=boundary_block_start.m;block++){
		int first_index=boundary_block_start(block);
		int last_index=boundary_block_end(block);
		for(int block_loop=1;block_loop<=10;block_loop++)
		    for(int i=first_index;i<=last_index;i++){
			const T_INDEX& index=boundary_indices(i);
			u(index)=b(index);
			for(int v=1;v<=d;v++)
			    u(index)-=u(index+T_INDEX::Axis_Vector(v))
				+u(index-T_INDEX::Axis_Vector(v));
			u(index)*=one_over_diagonal_part(i);
		    }
	    }
    }

#endif
}


//#####################################################################
// Function Apply_System_Matrix
//#####################################################################
template<class T,int d> T MULTIGRID_POISSON<T,d>::
Apply_System_Matrix(const T_INDEX& index)
{
    return Apply_System_Matrix(index,u);
}
template<class T,int d> T MULTIGRID_POISSON<T,d>::
Apply_System_Matrix(const T_INDEX& index,const T_VARIABLE& u_input) const
{
    if(!diagonal_entries(index))
	return (T)0;
    T lu=0;
    for(int v=1;v<=d;v++){
	T_INDEX neighbor_index = index+T_INDEX::Axis_Vector(v);
	if(diagonal_entries(neighbor_index))
	    lu+=u_input(neighbor_index);
	neighbor_index=index-T_INDEX::Axis_Vector(v);
	if(diagonal_entries(neighbor_index))
	    lu+=u_input(neighbor_index);
    }
    lu+=diagonal_entries(index)*u_input(index);
    return lu;
}
//#####################################################################
// Function Multiply_With_System_Matrix
//#####################################################################
template<class T,int d> void MULTIGRID_POISSON<T,d>::
Subtract_Multiple_Of_System_Matrix_And_Compute_Sum_And_Extrema(const T_VARIABLE& x,T_VARIABLE& y,double& sum,T& rmin,T& rmax) const
{
#ifndef MGPCG_UNOPTIMIZED
    Multiply_And_Compute_Sum_And_Extrema_Helper<T> multiplication_helper(n(1),n(2),n(3),&x(1,1,1),&y(1,1,1),&diagonal_entries(1,1,1),sum,rmin,rmax,(T)1/(h*h));
    multiplication_helper.Run_Parallel(number_of_threads);
#else
    
    sum=0;
    rmin=std::numeric_limits<T>::max();
    rmax=-std::numeric_limits<T>::max();
    T scale_factor=(T)1/(h*h);
    for(BOX_ITERATOR<d> iterator(unpadded_domain);iterator.Valid();iterator.Next()){
	const T_INDEX& index=iterator.Index();
	if(diagonal_entries(index)==0)
	    continue;
	T result=diagonal_entries(index)*x(index);
	for(int v=1;v<=d;v++){
	    T_INDEX axis_vector=T_INDEX::Axis_Vector(v);
	    if(diagonal_entries(index+axis_vector))
		result+=x(index+axis_vector);
	    if(diagonal_entries(index-axis_vector))
		result+=x(index-axis_vector);
	}
	result*=scale_factor;
	y(index)-=result;
	sum+=y(index);
	if(y(index)<rmin) 
	    rmin=y(index);
	if(y(index)>rmax)
	    rmax=y(index);
    }
#endif
}
//#####################################################################
// Function Multiply_With_System_Matrix_And_Compute_Dot_Product
//#####################################################################
template<class T,int d> T MULTIGRID_POISSON<T,d>::
Multiply_With_System_Matrix_And_Compute_Dot_Product(const T_VARIABLE& x,T_VARIABLE& y) const
{

#ifndef MGPCG_UNOPTIMIZED
    Multiplication_And_Dot_Product_Helper<T> multiplication_helper(n(1),n(2),n(3),&x(1,1,1),&y(1,1,1),&diagonal_entries(1,1,1),(T)1/(h*h));
    return multiplication_helper.Run_Parallel(number_of_threads);

#else

    T scale_factor=(T)1/(h*h);
    double dot_product=0;
    for(BOX_ITERATOR<d> iterator(unpadded_domain);iterator.Valid();iterator.Next()){
	const T_INDEX& index=iterator.Index();
	if(diagonal_entries(index)==0)
	    continue;
	T result=diagonal_entries(index)*x(index);
	for(int v=1;v<=d;v++){
	    T_INDEX axis_vector=T_INDEX::Axis_Vector(v);
	    if(diagonal_entries(index+axis_vector))
		result+=x(index+axis_vector);
	    if(diagonal_entries(index-axis_vector))
		result+=x(index-axis_vector);
	}
	result*=scale_factor;
	y(index)=result;
	dot_product+=y(index)*x(index);
    }
    return dot_product;
	

#endif
}
//#####################################################################
// Function Interior_Relaxation
//#####################################################################
template<class T,int d> T  MULTIGRID_POISSON<T,d>::
Interior_Relaxation(const int loops,const bool compute_dot_product)
{

    // compute deltas on entire unpadded domain, assuming diagonal entry is full 
    // (i.e. not on Neumann boundary)
    // apply deltas only where diagonal entry is actually full.
#ifndef MGPCG_UNOPTIMIZED
    for(int loop=1;loop<loops;loop++){
	Relaxation_Interior_Helper<T> relaxation(n(1),n(2),n(3),&u(1,1,1),&b(1,1,1),&delta(1,1,1),&index_has_full_diagonal_coarse_bitmask(1,1,1));
	relaxation.Run_Parallel(number_of_threads);}
        
    if(compute_dot_product){
	Relaxation_And_Dot_Product_Interior_Helper<T> relaxation(n(1),n(2),n(3),&u(1,1,1),&b(1,1,1),&delta(1,1,1),&index_has_full_diagonal_coarse_bitmask(1,1,1));
        return relaxation.Run_Parallel(number_of_threads);
    }
    else{
	Relaxation_Interior_Helper<T> relaxation(n(1),n(2),n(3),&u(1,1,1),&b(1,1,1),&delta(1,1,1),&index_has_full_diagonal_coarse_bitmask(1,1,1));
	relaxation.Run_Parallel(number_of_threads);}
    return T();
#else
    const T one_over_diagonal=(T)1/(2*d);
    const T two_thirds=2./3.;
    
    // NOTE: when running in parallel,
    // care must be taken to ensure 
    // any deltas dependent on u(index)
    // are computed before being applied to u(index)

    double dot_product=0;
    for(int loop=1;loop<=loops;loop++){
	

	// compute deltas on first coarse x-slice
	T_INDEX xmin_ymax_zmax=unpadded_coarse_domain.max_corner;
	xmin_ymax_zmax(1)=unpadded_coarse_domain.min_corner(1);
	for(BOX_ITERATOR<d> coarse_iterator(BOX<T_INDEX>(unpadded_coarse_domain.min_corner,xmin_ymax_zmax));
	    coarse_iterator.Valid();coarse_iterator.Next()){
	    
	    const T_INDEX& coarse_index=coarse_iterator.Index();
	    const T_INDEX base_fine_index=(coarse_index-1)*2;
	    
	    for(BOX_ITERATOR<d> fine_iterator(BOX<T_INDEX>(base_fine_index,base_fine_index+1));fine_iterator.Valid();fine_iterator.Next()){
		
		const T_INDEX& fine_index=fine_iterator.Index();
		delta(fine_index)=-b(fine_index);
		for(int v=1;v<=d;v++){
		    delta(fine_index)+=u(fine_index+T_INDEX::Axis_Vector(v));
		    delta(fine_index)+=u(fine_index-T_INDEX::Axis_Vector(v));
		}
		delta(fine_index)*=one_over_diagonal;
		delta(fine_index)-=u(fine_index);
	    }
	}
	
	// compute deltas on remainder of x-slices.
	// apply deltas, lagging by 1 x-slice
	// since write bitmask is based on coarse cells, 
	// we'll iterate over coarse domain
	for(BOX_ITERATOR<d> coarse_iterator(BOX<T_INDEX>(unpadded_coarse_domain.min_corner+T_INDEX::Axis_Vector(1),unpadded_coarse_domain.max_corner));
	    coarse_iterator.Valid();coarse_iterator.Next()){
	    const T_INDEX& coarse_index=coarse_iterator.Index();
	    T_INDEX shifted_coarse_index=coarse_index;
	    shifted_coarse_index(1)-=1;
	    const T_INDEX base_fine_index=(coarse_index-1)*2;
	    
	    unsigned char mask=0x01;
	    for(BOX_ITERATOR<d> fine_iterator(BOX<T_INDEX>(base_fine_index,base_fine_index+1));fine_iterator.Valid();fine_iterator.Next()){
		
		const T_INDEX& fine_index=fine_iterator.Index();		    
		T_INDEX shifted_fine_index=fine_index;
		shifted_fine_index(1)-=2;

		delta(fine_index)=-b(fine_index);
		for(int v=1;v<=d;v++){
		    delta(fine_index)+=u(fine_index+T_INDEX::Axis_Vector(v));
		    delta(fine_index)+=u(fine_index-T_INDEX::Axis_Vector(v));
		}
		delta(fine_index)*=one_over_diagonal;
		delta(fine_index)-=u(fine_index);
		
		if(index_has_full_diagonal_coarse_bitmask(shifted_coarse_index) & mask){
		    u(shifted_fine_index)+=two_thirds*delta(shifted_fine_index);
		}


		if(compute_dot_product && loop==loops)
		    dot_product+=u(shifted_fine_index)*b(shifted_fine_index);
		mask<<=0x01;
	    }
	}
	
	// apply deltas to last coarse x-slice
	T_INDEX xmax_ymin_zmin=unpadded_coarse_domain.min_corner;
	xmax_ymin_zmin(1)=unpadded_coarse_domain.max_corner(1);
	for(BOX_ITERATOR<d> coarse_iterator(BOX<T_INDEX>(xmax_ymin_zmin,unpadded_coarse_domain.max_corner));
	    coarse_iterator.Valid();coarse_iterator.Next()){
	    const T_INDEX& coarse_index=coarse_iterator.Index();
	    const T_INDEX base_fine_index=(coarse_index-1)*2;
	    
	    unsigned char mask=0x01;
	    for(BOX_ITERATOR<d> fine_iterator(BOX<T_INDEX>(base_fine_index,base_fine_index+1));fine_iterator.Valid();fine_iterator.Next()){

		const T_INDEX& fine_index=fine_iterator.Index();
		if(index_has_full_diagonal_coarse_bitmask(coarse_index) & mask){
		    u(fine_index)+=two_thirds*delta(fine_index);
		    
		}
		if(compute_dot_product && loop==loops)
		    dot_product+=u(fine_index)*b(fine_index);
		mask<<=0x01;
	    }
	}
    }
    return dot_product;
#endif
}
//#####################################################################
// Function Relaxation_Sweep
//#####################################################################
template<class T,int d> T MULTIGRID_POISSON<T,d>::
Relaxation_Sweep(const bool reverse_order,const int boundary_pre_relaxations, const int interior_relaxations, const int boundary_post_relaxations,const bool compute_dot_product,const bool verbose)
{
    T dot_product=T();
    // Boundary pre-relaxation
    if(boundary_pre_relaxations){
	Boundary_Relaxation(reverse_order,boundary_pre_relaxations);
    }
    
    // Interior relaxation

    dot_product=Interior_Relaxation(interior_relaxations,compute_dot_product);
    

    // Boundary post-relaxation
    if(boundary_post_relaxations){
	Boundary_Relaxation(reverse_order,boundary_post_relaxations);
    }
    return dot_product;
}
//#####################################################################
//Function Relax_And_Compute_Residuals
//#####################################################################
template<class T,int d> void MULTIGRID_POISSON<T,d>::
Relax_And_Compute_Residuals(const int interior_relaxations,const int boundary_post_relaxations,const T nullspace_component)
{

#ifndef MGPCG_UNOPTIMIZED
    Relaxation_And_Residual_With_Zero_Initial_Guess_Helper<T> relax_and_residual_helper(n(1),n(2),n(3),&u(1,1,1),&b(1,1,1),&delta(1,1,1),&index_has_full_diagonal_coarse_bitmask(1,1,1),&index_is_interior_coarse_bitmask(1,1,1),nullspace_component);
    relax_and_residual_helper.Run_Parallel(number_of_threads);
#else

    const T one_third=1./3.;
    const T two_thirds_over_diagonal_part=(T)1/(3*d);
    const T minus_two_thirds_over_diagonal_part=-(T)1/(3*d);
    
    for(BOX_ITERATOR<d> coarse_iterator(unpadded_coarse_domain);coarse_iterator.Valid();coarse_iterator.Next()){
        const T_INDEX& coarse_index=coarse_iterator.Index();
        const T_INDEX base_fine_index=(coarse_index-1)*2;
	
	unsigned char mask=0x01;
	for(BOX_ITERATOR<d> fine_iterator(BOX<T_INDEX>(base_fine_index,base_fine_index+1));fine_iterator.Valid();fine_iterator.Next()){
	    const T_INDEX& fine_index=fine_iterator.Index();

	    // subtract nullspace component from right hand side (interior indices only)
	    if(nullspace_component)
		if(index_is_interior_coarse_bitmask(coarse_index) & mask)
		    b(fine_index)-=nullspace_component;

	    // replace solution with relaxed value,
	    // or zero if it does not have full diagonal
	    if(index_has_full_diagonal_coarse_bitmask(coarse_index) & mask)
		u(fine_index)=minus_two_thirds_over_diagonal_part*b(fine_index);
	    else
		u(fine_index)=T();
	    
	    // update residual if no nullspace component.
	    // (we'll update it in a separate loop otherwise
	    if(!nullspace_component){
		delta(fine_index)=0;
		if(index_is_interior_coarse_bitmask(coarse_index) & mask){
		    delta(fine_index) = one_third*b(fine_index);
		    for(int v=1;v<=d;v++)
			delta(fine_index)+=two_thirds_over_diagonal_part*(b(fine_index+T_INDEX::Axis_Vector(v))
			    +b(fine_index-T_INDEX::Axis_Vector(v)));
		}
	    }
	    mask<<=1;
	}
    }
		
    if(nullspace_component)
	for(BOX_ITERATOR<d> coarse_iterator(unpadded_coarse_domain);coarse_iterator.Valid();coarse_iterator.Next()){
	    const T_INDEX& coarse_index=coarse_iterator.Index();
	    const T_INDEX base_fine_index=(coarse_index-1)*2;
	
	    unsigned char mask=0x01;
	    for(BOX_ITERATOR<d> fine_iterator(BOX<T_INDEX>(base_fine_index,base_fine_index+1));fine_iterator.Valid();fine_iterator.Next()){
		const T_INDEX& fine_index=fine_iterator.Index();

	    	delta(fine_index)=0;
		if(index_is_interior_coarse_bitmask(coarse_index) & mask){
		    delta(fine_index) = one_third*b(fine_index);
		    for(int v=1;v<=d;v++)
			delta(fine_index)+=two_thirds_over_diagonal_part*(b(fine_index+T_INDEX::Axis_Vector(v))
			    +b(fine_index-T_INDEX::Axis_Vector(v)));
		}
		mask<<=1;
	    }
	}
#endif

    Boundary_Relaxation(true,boundary_post_relaxations);
    Compute_Residuals_Boundary();
}
//#####################################################################
//Function Compute_Residuals_Boundary
//#####################################################################
template<class T,int d> void MULTIGRID_POISSON<T,d>::
Compute_Residuals_Boundary()
{

#ifndef MGPCG_UNOPTIMIZED
    if (extended_boundary_indices.m==0) return;
    Residual_Boundary_Helper<T> residual_helper(n(1),n(2),n(3),&u(1,1,1),&b(1,1,1),&delta(1,1,1),&diagonal_entries(1,1,1),&extended_boundary_indices(1),extended_boundary_indices.m);
    residual_helper.Run_Parallel(number_of_threads);
#else
    for(int i=1;i<=extended_boundary_indices.m;i++){
	const T_INDEX& index=extended_boundary_indices(i);
	if(diagonal_entries(index)!=0){
	    delta(index)=b(index)-diagonal_entries(index)*u(index);
	    for(int v=1;v<=d;v++)
		delta(index)-=u(index+T_INDEX::Axis_Vector(v))
		    +u(index-T_INDEX::Axis_Vector(v));
	}
	else
	    delta(index)=0;
    }

#endif
}

//#####################################################################
// INITIALIZATION HELPERS FOR DEBUGGING PURPOSES
//#####################################################################
template<class T,int d> void MULTIGRID_POISSON<T,d>::
Initialize_Square_Domain()
{
    for(BOX_ITERATOR<d> iterator(unpadded_domain);iterator.Valid();iterator.Next())
	cell_type(iterator.Index())=INTERIOR_CELL_TYPE;

    for(BOUNDARY_ITERATOR<d> iterator(padded_domain);iterator.Valid();iterator.Next())
	cell_type(iterator.Index())=NEUMANN_CELL_TYPE;
}
template<class T,int d> void MULTIGRID_POISSON<T,d>::
Initialize_Square_Minus_Circle_Domain()
{

   T radius=.125;

   TV circle_center=TV::All_Ones_Vector()*.5;
   circle_center(2)=2*radius;
    T radius_squared=radius*radius;

    for(BOUNDARY_ITERATOR<d> iterator(padded_domain);iterator.Valid();iterator.Next())
	cell_type(iterator.Index())=NEUMANN_CELL_TYPE;

    for(BOX_ITERATOR<d> iterator(unpadded_domain);iterator.Valid();iterator.Next()){
	TV X=grid.Node(iterator.Index());
	if((X-circle_center).Magnitude_Squared()<radius_squared)
            cell_type(iterator.Index())=NEUMANN_CELL_TYPE;
        else
            cell_type(iterator.Index())=INTERIOR_CELL_TYPE;
    }
}

template<class T,int d> void MULTIGRID_POISSON<T,d>::
Initialize_Test_Domain()
{

    // Sinusoidal dirichlet interface, with 2 Neumann bubbles
    for(BOX_ITERATOR<d> iterator(unpadded_domain);iterator.Valid();iterator.Next()){
	TV X=grid.Node(iterator.Index());
	if((X-TV::All_Ones_Vector()*(T).25).Magnitude()<.125){
	    cell_type(iterator.Index())=NEUMANN_CELL_TYPE;}
        else if((X-TV::All_Ones_Vector()*(T).5).Magnitude()<.125){
	    cell_type(iterator.Index())=NEUMANN_CELL_TYPE;}
	else if(d==2 && (X(2)<.5+.25*sin((X(1))*2*pi)))
            cell_type(iterator.Index())=INTERIOR_CELL_TYPE;
        else if(d==3 && (X(2)<.5+.25*sin((X(1)+X(3))*2*pi)))
            cell_type(iterator.Index())=INTERIOR_CELL_TYPE;
        else
            cell_type(iterator.Index())=DIRICHLET_CELL_TYPE;
    }


    for(BOUNDARY_ITERATOR<d> iterator(padded_domain);iterator.Valid();iterator.Next())
	cell_type(iterator.Index())=NEUMANN_CELL_TYPE;
}
template<class T,int d> void MULTIGRID_POISSON<T,d>::
Initialize_Test_Right_Hand_Side(T_VARIABLE& b_input)
{
#if 0
    // Use a discrete delta as right hand side
    T_INDEX spike_index(n);spike_index(1)/=2;for(int v=2;v<=d;v++) spike_index(v)/=8;
    b_input(spike_index)=-.5;
#else
    for(BOX_ITERATOR<d> iterator(unpadded_domain);iterator.Valid();iterator.Next()){
	TV X=grid.Node(iterator.Index());
	if(cell_type(iterator.Index())==INTERIOR_CELL_TYPE){
	    if(d==2)
		b_input(iterator.Index())=sin(X(1)*2*pi)*cos(X(2)*2*pi);
	    else if(d==3)
		b_input(iterator.Index())=sin(X(1)*2*pi)*cos(X(2)*2*pi)*sin(X(3)*2*pi);
	}
    }
#endif

}
template<class T,int d> void MULTIGRID_POISSON<T,d>::
Initialize_Test_Right_Hand_Side()
{
    Initialize_Test_Right_Hand_Side(b);
}
template<class T,int d> void MULTIGRID_POISSON<T,d>::
Initialize_Test_Initial_Guess(T_VARIABLE& u_input)
{
    #if 0
    // Randomize initial guess
    RANDOM_NUMBERS<T> random_numbers;
    random_numbers.Set_Seed(1);
    for(BOX_ITERATOR<d> iterator(unpadded_domain);iterator.Valid();iterator.Next()){
        const T_INDEX& index=iterator.Index();
        if(cell_type(index)==INTERIOR_CELL_TYPE)
            u_input(index)=random_numbers.Get_Number();}
#endif
}
template<class T,int d> void MULTIGRID_POISSON<T,d>::
Initialize_Test_Initial_Guess()
{
    Initialize_Test_Initial_Guess(u);
}
//#####################################################################
#ifdef MGPCG_UNOPTIMIZED
template class MULTIGRID_POISSON<float,2>;
#endif
template class MULTIGRID_POISSON<float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
#ifdef MGPCG_UNOPTIMIZED
template class MULTIGRID_POISSON<double,2>;
#endif
template class MULTIGRID_POISSON<double,3>;
#endif
