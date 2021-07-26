/*
 *	svd3.h
 *	
 *	Created by Ryoichi Ando on 1/8/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "macros.h"
#include "array3.h"
#include <stdio.h>

#ifndef _SVD3_H
#define _SVD3_H

class svd3 {
public:
	svd3() {
		eig[0] = 1.0;
		eig[1] = 1.0;
		eig[2] = 1.0;
		for( uint i=0; i<DIM; i++ ) for( uint j=0; j<DIM; j++ ) {
			R[i][j] = i==j;
		}
	}
	bool run( FLOAT64 matrix[DIM][DIM] );
	// Directly access these varables to retrive results
	FLOAT64 R[DIM][DIM];
	FLOAT64 eig[DIM];
};

#endif