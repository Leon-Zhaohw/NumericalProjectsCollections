//#####################################################################
// Copyright 2006-2010, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis,Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
// Class MG_PRECONDITIONED_CONJUGATE_GRADIENT
//#####################################################################
#ifndef __MG_PRECONDITIONED_CONJUGATE_GRADIENT__
#define __MG_PRECONDITIONED_CONJUGATE_GRADIENT__

#include <Utilities/LOG.h>
#include <limits>
#include <cfloat>

#ifndef MGPCG_UNOPTIMIZED
#include <Poisson_MGPCG_optimized_kernels/Combined_Saxpy_Helper.h>
#include <Poisson_MGPCG_optimized_kernels/Saxpy_Helper.h>
#include <Poisson_MGPCG_optimized_kernels/Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema_Helper.h>
#include <Poisson_MGPCG_optimized_kernels/Reduce_Maxabs_Helper.h>
#endif

namespace PhysBAM{
template<class T,int d>
class MULTIGRID_SYSTEM
{
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    typedef typename GRID_POLICY<TV>::UNIFORM_GRID T_GRID;
    typedef typename T_GRID::NODE_ITERATOR T_NODE_ITERATOR;
    typedef typename POLICY_UNIFORM<TV>::ARRAYS_SCALAR T_VARIABLE;
    
    MULTIGRID_POISSON_SOLVER<T,d>& multigrid_poisson_solver;
    MULTIGRID_POISSON<T,d>& multigrid_poisson;
public:
    MULTIGRID_SYSTEM(MULTIGRID_POISSON_SOLVER<T,d>& multigrid_poisson_solver_input)
	:multigrid_poisson_solver(multigrid_poisson_solver_input),multigrid_poisson(multigrid_poisson_solver_input.Discretization())
    {}

    void Subtract_Multiple_Of_System_Matrix_And_Compute_Sum_And_Extrema(const T_VARIABLE& x,T_VARIABLE& result,double& sum,T& rmin,T& rmax) const
    {
	LOG::SCOPE scope("[MULTIGRID_SYSTEM] Subtract_Multiple_Of_System_Matrix_And_Compute_Sum_And_Extrema");
	multigrid_poisson.Subtract_Multiple_Of_System_Matrix_And_Compute_Sum_And_Extrema(x,result,sum,rmin,rmax);
    }

    void Multiply_And_Compute_Dot_Product(const T_VARIABLE& x, T_VARIABLE& result, T& dot_product)
    {
	LOG::SCOPE scope("[MULTIGRID_SYSTEM] Multiply_And_Compute_Dot_Product");
 	dot_product=multigrid_poisson.Multiply_With_System_Matrix_And_Compute_Dot_Product(x,result);
    }

    void Saxpy(const T c1, const T_VARIABLE& v1,T_VARIABLE& result)
    {
	LOG::SCOPE scope("[MULTIGRID_SYSTEM] Saxpy ");
#ifndef MGPCG_UNOPTIMIZED
	Saxpy_Helper<T> copy_helper(multigrid_poisson.n(1),multigrid_poisson.n(2),multigrid_poisson.n(3),c1,&v1(1,1,1),&result(1,1,1));
	copy_helper.Run_Parallel(multigrid_poisson_solver.number_of_threads);
#else
	for(BOX_ITERATOR<d> iterator(multigrid_poisson.unpadded_domain);iterator.Valid();iterator.Next()){
	    const T_INDEX& index=iterator.Index();
	    result(index)=c1*v1(index);
	}
#endif
    }
    
    void Combined_Saxpy(T_VARIABLE& x,T_VARIABLE& p,const T_VARIABLE& z,const T alpha,const T beta)
    {
	LOG::SCOPE scope("[MULTIGRID_SYSTEM] Combined_Saxpy ");
#ifndef MGPCG_UNOPTIMIZED
	Combined_Saxpy_Helper<T> combined_copy_helper(multigrid_poisson.n(1),multigrid_poisson.n(2),multigrid_poisson.n(3),&x(1,1,1),&p(1,1,1),&z(1,1,1),alpha,beta);
	combined_copy_helper.Run_Parallel(multigrid_poisson_solver.number_of_threads);
#else
	for(BOX_ITERATOR<d> iterator(multigrid_poisson.unpadded_domain);iterator.Valid();iterator.Next()){
	    const T_INDEX& index=iterator.Index();
	    x(index)+=alpha*p(index);
	    p(index)=z(index)+beta*p(index);
	}
#endif
    }
    
    void Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema(const T c1, const T_VARIABLE& v1, T_VARIABLE& result,double &sum,T& minimum,T& maximum)
    {
	LOG::SCOPE scope("[MULTIGRID_SYSTEM] Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema ");
#ifndef MGPCG_UNOPTIMIZED
	Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema_Helper<T> scalar_multiply_and_accumulate_and_compute_sum_and_extrema_helper(multigrid_poisson.n(1),multigrid_poisson.n(2),multigrid_poisson.n(3),c1,&v1(1,1,1),&result(1,1,1),sum,minimum,maximum);
	scalar_multiply_and_accumulate_and_compute_sum_and_extrema_helper.Run_Parallel(multigrid_poisson_solver.number_of_threads);
#else
	sum=0;
	minimum=std::numeric_limits<T>::max();
	maximum=-minimum;
	for(BOX_ITERATOR<d> iterator(multigrid_poisson.unpadded_domain);iterator.Valid();iterator.Next()){
	    const T_INDEX& index=iterator.Index();
	    result(index)+=c1*v1(index);
	    
	    sum+=result(index);
	    if(result(index)<minimum)
		minimum=result(index);
	    if(result(index)>maximum)
		maximum=result(index);
	}
#endif
    }
    
    T Convergence_Norm(const T_VARIABLE& x) const
    {
	LOG::SCOPE scope("[MULTIGRID_SYSTEM] Convergence_Norm");
#ifndef MGPCG_UNOPTIMIZED
	Reduce_Maxabs_Helper<T> maxabs_helper(multigrid_poisson.n(1),multigrid_poisson.n(2),multigrid_poisson.n(3),&x(1,1,1));
	return maxabs_helper.Run_Parallel(multigrid_poisson_solver.number_of_threads);
#else
	T maxabs=0;
	
	for(BOX_ITERATOR<d> iterator(multigrid_poisson.unpadded_domain);iterator.Valid();iterator.Next()){
	    const T_INDEX& index=iterator.Index();
	    maxabs=std::max(maxabs,u(index)<0?-u(index):u(index));
	}
	return maxabs;
#endif
    }

    // removes component of x in nullspace of A (used to project residual for stopping conditions)
    double One_Over_Interior_Count(T_VARIABLE& x) const
    {
	LOG::SCOPE scope ("[MULTIGRID_SYSTEM] Interior_Count");
	double count=0;
	bool has_dirichlet=false;
 	for(BOX_ITERATOR<d> iterator(multigrid_poisson.padded_domain);iterator.Valid();iterator.Next()){
 		const T_INDEX& index=iterator.Index();
 		if(multigrid_poisson.cell_type(index)==MULTIGRID_POISSON<T,d>::INTERIOR_CELL_TYPE) count+=(T)1.;
 		if(multigrid_poisson.cell_type(index)==MULTIGRID_POISSON<T,d>::DIRICHLET_CELL_TYPE) has_dirichlet=true;}
        if(has_dirichlet) return 0;else return ((double)1)/count;
    }

    const T_VARIABLE& Precondition_And_Compute_Dot_Product(const T_VARIABLE& r,T_VARIABLE& z,T& dot_product,const T nullspace_component) const
    {
 	LOG::SCOPE scope("[MULTIGRID_SYSTEM] Precondition_And_Compute_Dot_Product");
 	dot_product=multigrid_poisson_solver.V_Cycle(nullspace_component);
	return z;
    }
    
};



template<class T,int d>
class MG_PRECONDITIONED_CONJUGATE_GRADIENT
{
    typedef typename POLICY_UNIFORM<VECTOR<T,d> >::ARRAYS_SCALAR VECTOR_T;
public:
    bool print_diagnostics,print_residuals;
    int* iterations_used;
    int restart_iterations;


    MG_PRECONDITIONED_CONJUGATE_GRADIENT()
        :print_diagnostics(true),print_residuals(false),// nullspace_tolerance((T)1e-5),
	iterations_used(0),// residual_magnitude_squared(0),nullspace_measure(0),
	restart_iterations(0)
    {}

    ~MG_PRECONDITIONED_CONJUGATE_GRADIENT()
    {}

    bool Solve(MULTIGRID_SYSTEM<T,d>& system,VECTOR_T& x,VECTOR_T& r,VECTOR_T& z,VECTOR_T& p,const T tolerance,const int min_iterations,const int max_iterations) PHYSBAM_OVERRIDE
    {
        T rho=0,rho_new=0,p_dot_z=0,nullspace_component=0,rmax=0,rmin=0;
        double sum=0,one_over_interior_count=0;

        // This should be precomputation
        one_over_interior_count=system.One_Over_Interior_Count(r);

        system.Subtract_Multiple_Of_System_Matrix_And_Compute_Sum_And_Extrema(x,r,sum,rmin,rmax);
        nullspace_component=(T)(one_over_interior_count*sum);

        T convergence_norm=maxabs(rmin-nullspace_component,rmax-nullspace_component);

        if(print_residuals) LOG::cout<<"Norm : "<<convergence_norm<<std::endl;
        if(convergence_norm<=tolerance){
            if(print_diagnostics) LOG::Stat("cg iterations",0);return true;}

	VECTOR_T::Exchange_Arrays(z,p);
	system.Precondition_And_Compute_Dot_Product(r,p,rho,nullspace_component);
	VECTOR_T::Exchange_Arrays(z,p);

        int iterations;for(iterations=1;;iterations++){
            LOG::SCOPE scope("PCG iteration","PCG iteration %d",iterations);
	    system.Multiply_And_Compute_Dot_Product(p,z,p_dot_z);

            T alpha=rho/p_dot_z;

            system.Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema(-alpha,z,r,sum,rmin,rmax);

            nullspace_component=(T)(one_over_interior_count*sum);

            convergence_norm=maxabs(rmin-nullspace_component,rmax-nullspace_component);

            if(print_residuals) LOG::cout<<"Norm : "<<convergence_norm<<std::endl;
            if(convergence_norm<=tolerance){
                if(print_diagnostics) LOG::Stat("cg iterations",iterations);if(iterations_used) *iterations_used=iterations;
                system.Saxpy(alpha,p,x);
                return true;}
            if(iterations==max_iterations){
                system.Saxpy(alpha,p,x);
                break;
            }
	    system.Precondition_And_Compute_Dot_Product(r,z,rho_new,nullspace_component);

            T beta=rho_new/rho;
            rho=rho_new;
	    system.Combined_Saxpy(x,p,z,alpha,beta);

        }

        if(print_diagnostics) LOG::Stat("cg iterations",iterations);if(iterations_used) *iterations_used=iterations;
        if(print_diagnostics) LOG::cout<<"cg not converged after "<<max_iterations<<" iterations, error = "<<convergence_norm<<std::endl;
        return false;
    }
};
}
#endif
