//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis,Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
// Class MULTIGRID_POISSON_REFINEMENT
//#####################################################################
#include "MULTIGRID_POISSON_REFINEMENT.h"

#ifndef MGPCG_UNOPTIMIZED
#include <Poisson_MGPCG_optimized_kernels/Restriction_Helper.h>
#include <Poisson_MGPCG_optimized_kernels/Prolongation_Helper.h>
#endif

using namespace PhysBAM;

//#####################################################################
// Function Transfer_Residual_To_Coarse_Grid
//#####################################################################
template<class T,int d> void MULTIGRID_POISSON_REFINEMENT<T,d>::
Transfer_Residual_To_Coarse_Grid()
{
#ifndef MGPCG_UNOPTIMIZED
    Restriction_Helper<T> restriction_helper(fine.n(1),fine.n(2),fine.n(3),&fine.delta(1,1,1),&coarse.b(1,1,1),&coarse.cell_type(1,1,1));
    restriction_helper.Run_Parallel(fine.number_of_threads);
#else
    BOX<T_INDEX> restriction_indices(-T_INDEX::All_Ones_Vector(),2*T_INDEX::All_Ones_Vector());
    T_VARIABLE restriction_stencil(restriction_indices);
    for(BOX_ITERATOR<d> iterator(restriction_indices);iterator.Valid();iterator.Next()){
	const T_INDEX& offset_index=iterator.Index();
 	T scale=(T)1/(1<<(d-2));
	for(int v=1;v<=d;v++)
	    switch(offset_index(v)){
		case -1:
		case 2:
		    scale*=.25;
		    break;
		case 0:
		case 1:
		    scale*=.75;
		    break;
	    }
	restriction_stencil(offset_index)=scale;
    }


    for(BOX_ITERATOR<d> iterator(coarse.unpadded_domain);iterator.Valid();iterator.Next()){
	const T_INDEX& coarse_index=iterator.Index();
	T_INDEX base_fine_index=(coarse_index-1)*2;
	coarse.b(coarse_index)=0;
	if(coarse.cell_type(coarse_index)!=MULTIGRID_POISSON<T,d>::INTERIOR_CELL_TYPE)
	    continue;
	for(BOX_ITERATOR<d> local_iterator(restriction_indices);local_iterator.Valid();local_iterator.Next()){
	    const T_INDEX& offset_index=local_iterator.Index();
	    coarse.b(coarse_index)+=restriction_stencil(offset_index)*fine.delta(base_fine_index+offset_index);
	}
    }
#endif    
}

//#####################################################################
// Function Transfer_Correction_To_Fine_Grid
//#####################################################################
template<class T,int d> void MULTIGRID_POISSON_REFINEMENT<T,d>::
Transfer_Correction_To_Fine_Grid()
{
#ifndef MGPCG_UNOPTIMIZED
    Prolongation_Helper<T> prolongation(fine.n(1),fine.n(2),fine.n(3),&fine.u(1,1,1),&coarse.u(1,1,1),&fine.index_is_interior_coarse_bitmask(1,1,1));
    prolongation.Run_Parallel(fine.number_of_threads);
#else
    BOX<T_INDEX> prolongation_indices(-T_INDEX::All_Ones_Vector(),T_INDEX::All_Ones_Vector());
    T_VARIABLE prolongation_stencil(prolongation_indices);
    for(BOX_ITERATOR<d> iterator(prolongation_indices);iterator.Valid();iterator.Next()){
	const T_INDEX& offset_index=iterator.Index();
	T scale=(T)1;
	for(int v=1;v<=d;v++)
	    switch(offset_index(v)){
		case -1:
		case 1:
		    scale*=.25;
		    break;
		case 0:
		    scale*=.75;
		    break;
	    }
	prolongation_stencil(offset_index)=scale;
    }

    for(BOX_ITERATOR<d> coarse_iterator(coarse.unpadded_domain);coarse_iterator.Valid();coarse_iterator.Next()){
	const T_INDEX& coarse_index=coarse_iterator.Index();
	T_INDEX base_fine_index=(coarse_index-1)*2;
	unsigned char mask=0x01;
	for(BOX_ITERATOR<d> fine_offset_iterator(BOX<T_INDEX>(T_INDEX(),T_INDEX::All_Ones_Vector()));fine_offset_iterator.Valid();fine_offset_iterator.Next()){
	    const T_INDEX& fine_offset=fine_offset_iterator.Index();
	    if(!(fine.index_is_interior_coarse_bitmask(coarse_index) & mask)){
		fine.u(base_fine_index+fine_offset)=0;
		mask<<=1;
		continue;
	    }
	    
	    for(BOX_ITERATOR<d> coarse_offset_iterator(BOX<T_INDEX>(fine_offset-1,fine_offset));coarse_offset_iterator.Valid();coarse_offset_iterator.Next()){
		const T_INDEX& coarse_offset=coarse_offset_iterator.Index();
		
		fine.u(base_fine_index+fine_offset)+=prolongation_stencil(coarse_offset)*coarse.u(coarse_index+coarse_offset);
	    }

	    mask<<=1;
	}
    }
	
#endif
}


//#####################################################################
#ifdef MGPCG_UNOPTIMIZED
template class MULTIGRID_POISSON_REFINEMENT<float,2>;
#endif
template class MULTIGRID_POISSON_REFINEMENT<float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
#ifdef MGPCG_UNOPTIMIZED
template class MULTIGRID_POISSON_REFINEMENT<double,2>;
#endif
template class MULTIGRID_POISSON_REFINEMENT<double,3>;
#endif
