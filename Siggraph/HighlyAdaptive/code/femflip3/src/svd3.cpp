/*
 *	svd3.cpp
 *	
 *	Created by Ryoichi Ando on 1/10/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "macros.h"
#include "svd3.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

bool svd3::run( FLOAT64 matrix[DIM][DIM] ) {
	// Alloc Matrix
	gsl_matrix *M = gsl_matrix_alloc(DIM,DIM);
	
	// Fill In
	for( int i=0; i<DIM; i++ ) for( int j=0; j<DIM; j++ ) {
		gsl_matrix_set(M,i,j,matrix[i][j]);
	}
	// Compute Eigen Vectors
	gsl_vector *r = gsl_vector_alloc(DIM);
	gsl_vector *w = gsl_vector_alloc(DIM);
	gsl_matrix *v = gsl_matrix_alloc(DIM,DIM);
	gsl_linalg_SV_decomp(M,v,r,w);
	gsl_eigen_symmv_sort(r,v,GSL_EIGEN_SORT_VAL_DESC );
	
	// Extract Eigen Vector Matrix
	for( int i=0; i<DIM; i++ ) for( int j=0; j<DIM; j++ ) {
		R[i][j] = gsl_matrix_get(v,i,j);
	}
	// Extract Eigen Values
	for( int i=0; i<DIM; i++ ) {
		eig[i] = gsl_vector_get(r,i);
	}
	// Deallocation
	gsl_matrix_free(M);
	gsl_matrix_free(v);
	gsl_vector_free(r);
	gsl_vector_free(w);
	return true;
}