//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis, Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
// Class MULTIGRID_POISSON_SOLVER
//#####################################################################
#ifndef __MULTIGRID_POISSON_SOLVER__
#define __MULTIGRID_POISSON_SOLVER__
#include "MULTIGRID_POISSON.h"
#include "MULTIGRID_POISSON_REFINEMENT.h"
namespace PhysBAM{

template<class T,int d>
class MULTIGRID_POISSON_SOLVER
{
    typedef VECTOR<int,d> T_INDEX;
    typedef typename GRID_POLICY<VECTOR<T,d> >::UNIFORM_GRID T_GRID;
    typedef typename T_GRID::NODE_ITERATOR T_NODE_ITERATOR;
    typedef typename MULTIGRID_POISSON<T,d>::CELL_TYPE T_CELL_TYPE;

public:
    const int number_of_threads;

private:
    const T_INDEX n;
    const T h;
    const int levels;

    ARRAY<MULTIGRID_POISSON<T,d>*> discretizations;
    ARRAY<MULTIGRID_POISSON_REFINEMENT<T,d>*> refinements;
public:

    MULTIGRID_POISSON_SOLVER()
	:number_of_threads(1),n(0),h(0),levels(0)
    {}


    MULTIGRID_POISSON_SOLVER(const T_INDEX n_input,const T h_input,const int levels_input, const int number_of_threads_input)
        :number_of_threads(number_of_threads_input),n(n_input),h(h_input),levels(levels_input)
    {
        discretizations.Resize(levels);
        refinements.Resize(levels-1);
        discretizations(1)=new MULTIGRID_POISSON<T,d>(n,h,number_of_threads);
    }

    ~MULTIGRID_POISSON_SOLVER()
    {
        refinements.Delete_Pointers_And_Clean_Memory();
        discretizations.Delete_Pointers_And_Clean_Memory();
    }

    void Initialize_Multigrid_Hierarchy()
    {discretizations(1)->Initialize();
    for(int level=1;level<levels;level++){
        discretizations(level+1)=MULTIGRID_POISSON_REFINEMENT<T,d>::Coarsened_Discretization(*discretizations(level));
        refinements(level)=new MULTIGRID_POISSON_REFINEMENT<T,d>(*discretizations(level),*discretizations(level+1));}}

    void Reinitialize_Multigrid_Hierarchy()
    { 
	discretizations(1)->Reinitialize();
	refinements.Delete_Pointers_And_Clean_Memory();
        refinements.Resize(levels-1);
	   for(int level=1;level<levels;level++){
	       delete discretizations(level+1);
	       discretizations(level+1)=MULTIGRID_POISSON_REFINEMENT<T,d>::Coarsened_Discretization(*discretizations(level));
	       refinements(level)=new MULTIGRID_POISSON_REFINEMENT<T,d>(*discretizations(level),*discretizations(level+1));}
    }

    T_CELL_TYPE& Cell_Type(const T_INDEX& index)
    {return discretizations(1)->cell_type(index);}

    T& U(const T_INDEX& index)
    {return discretizations(1)->u(index);}

    T& B(const T_INDEX& index)
    {return discretizations(1)->b(index);}

    const T_GRID& Grid()
    {return discretizations(1)->grid;}

    MULTIGRID_POISSON<T,d>& Discretization() const
    {
	assert(discretizations.Size());
	return *discretizations(1);}


    MULTIGRID_POISSON_REFINEMENT<T,d>& Refinement() const
    {return *refinements(1);}
    
    T V_Cycle(const T nullspace_component)
    {
	T dot_product=T();

	int min_boundary_sweeps=2;
	int max_boundary_sweeps=min(16,min_boundary_sweeps<<(levels-2));
        int boundary_sweeps_at_coarsest_level=2;

        PHYSBAM_ASSERT(levels>=1); // Needed for the removal of nullspace component

        {//LOG::SCOPE scope("V-cycle descent","V-cycle descent");
        for(int level=1;level<levels;level++){
            int boundary_sweeps=min(max_boundary_sweeps,min_boundary_sweeps<<(level-1));

            T nullspace_component_at_this_level=(level==1)?nullspace_component:0;
	    discretizations(level)->Relax_And_Compute_Residuals(1,boundary_sweeps,nullspace_component_at_this_level);
	    refinements(level)->Transfer_Residual_To_Coarse_Grid();
  
	}}
            
        {//LOG::SCOPE scope("V-cycle bottom","V-cycle bottom");
        discretizations(levels)->u.Fill(0);
        for(int coarse_loop=1;coarse_loop<=50;coarse_loop++)
            discretizations(levels)->Relaxation_Sweep(true,boundary_sweeps_at_coarsest_level,1,boundary_sweeps_at_coarsest_level,false,false);}

        {//LOG::SCOPE scope("V-cycle bottom","V-cycle bottom");
        for(int coarse_loop=1;coarse_loop<=50;coarse_loop++)
            discretizations(levels)->Relaxation_Sweep(false,boundary_sweeps_at_coarsest_level,1,boundary_sweeps_at_coarsest_level,false,false);}

        {//LOG::SCOPE scope("V-cycle ascent","V-cycle ascent");
        for(int level=levels-1;level>=1;level--){
            refinements(level)->Transfer_Correction_To_Fine_Grid();
            int boundary_sweeps=min(max_boundary_sweeps,min_boundary_sweeps<<(level-1));
	    dot_product=discretizations(level)->Relaxation_Sweep(false,boundary_sweeps,1,0,level==1);	}}

	return dot_product;
    }

//#####################################################################
};
}
#endif
