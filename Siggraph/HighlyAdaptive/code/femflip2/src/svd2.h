/*
 *	svd2.h
 *	
 *	Created by Ryoichi Ando on 11/12/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "macros.h"
#include "array2.h"
#include <stdio.h>

#ifndef _SVD2_H
#define _SVD2_H

class svd2 {
public:
	svd2() {
		eig[0] = 1.0;
		eig[1] = 1.0;
		FOR_EACH(DIM,DIM) {
			R[i][j] = i==j;
		} END_FOR
	}
	bool run( FLOAT64 matrix[DIM][DIM] ) {
		// Make sure matrix is symmetric
		FLOAT64 symdif = fabs(matrix[0][1]-matrix[1][0]);
		if( symdif > 1e-6 ) {
			printf( "Matrix should be symmetric ! dif = %lf\n", symdif );
			return false;
		}
		// Compute M*M^T
		array2<FLOAT64> M(2,2);
		FOR_EACH(DIM,DIM) {
			for( int k=0; k<DIM; k++ ) M[i][j] += matrix[i][k]*matrix[k][j];
		} END_FOR
		// First compute eigen values, easy enough for 2x2 matrix isn't it ?
		FLOAT64 a;
		FLOAT64 b = -(M[0][0]+M[1][1]);
		FLOAT64 c = M[0][0]*M[1][1]-M[0][1]*M[1][0];
		FLOAT64 D = fmax(0,sqrt(b*b-4*c));
		FLOAT64 eig[2];
		eig[0] = (-b+D)/2;
		eig[1] = (-b-D)/2;
		// Then compute eigen vectors
		for( int c=0; c<DIM; c++ ) {
			a = M[0][0]-eig[c];
			b = M[0][1];
			FLOAT64 det = sqrt(1.0/(a*a+b*b));
			FLOAT64 sgn = b > 0 ? 1 : -1;
			R[0][c] = fabs(b)*det;
			R[1][c] = -a*sgn*det;
		}
		 
		// Apply square root and get SVD eigen components
		this->eig[0] = sqrt(fmax(0,eig[0]));
		this->eig[1] = sqrt(fmax(0,eig[1]));
		return true;
	}
	// Directly access these varables to retrive results
	FLOAT64 R[DIM][DIM];
	FLOAT64 eig[DIM];
	
	bool operator==(const svd2 svd) const {
		for( uint dim=0; dim<DIM; dim++ ) {
			if( eig[dim] != svd.eig[dim] ) return false;
		}
		for( uint i=0; i<DIM; i++ ) for( uint j=0; j<DIM; j++ ) {
			if( R[i][j] != svd.R[i][j] ) return false;
		}
		return true;
	}
};

#endif