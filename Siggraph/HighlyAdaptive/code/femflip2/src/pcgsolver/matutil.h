/*
 *	matutil.h
 *	
 *	Created by Ryoichi Ando on 12/6/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "macros.h"
#include "email.h"
#include "sparse_matrix.h"
#include <stdio.h>
#include <iostream>
#include <fstream>

#ifndef MATUTIL_H
#define MATUTIL_H

static void interp( const SparseMatrix<FLOAT64> &A, const SparseMatrix<FLOAT64> &B, SparseMatrix<FLOAT64> &C, FLOAT64 a, FLOAT64 b ) {
	if( A.index.size() != B.index.size() ) {
		email::print( "Matrix size not match !\n" );
		email::send();
		exit(0);
	}
	C = SparseMatrix<FLOAT64>(A.index.size());
	// A
	for( uint i=0; i<A.index.size(); i++ ) {
		for( uint j=0; j<A.index[i].size(); j++ ) {
			C.add_to_element(i,A.index[i][j],a*A.value[i][j]);
		}
	}
	// B
	for( uint i=0; i<B.index.size(); i++ ) {
		for( uint j=0; j<B.index[i].size(); j++ ) {
			C.add_to_element(i,B.index[i][j],b*B.value[i][j]);
		}
	}
}
				   
static void normalize_row( SparseMatrix<FLOAT64> &InterpMatrix ) {
	// Normalize
	for( uint n=0; n<InterpMatrix.value.size(); n++ ) {
		FLOAT64 sum = 0.0;
		for( uint m=0; m<InterpMatrix.value[n].size(); m++ ) {
			sum += InterpMatrix.value[n][m];
		}
		if( sum ) {
			for( uint m=0; m<InterpMatrix.value[n].size(); m++ ) {
				InterpMatrix.value[n][m] /= sum;
			}
		}
	}
}

template<class T>
static void compress( SparseMatrix<T> &matrix, std::vector<T> *rhs, uint num, std::vector<int> &index_map ) {
	// Compact the matrix
	SparseMatrix<T> final_matrix;
	std::vector<T> *final_rhs = new std::vector<T>[num];
	index_map.resize(matrix.index.size());
	
	uint validsum = 0;
	uint i=0;
	for( uint n=0; n<matrix.index.size(); n++ ) {
		if( matrix(n,n) ) {
			validsum ++;
			index_map[n] = i++;
		} else {
			index_map[n] = -1;
		}
	}
	final_matrix.resize(validsum);
	for( uint k=0; k<num; k++ ) final_rhs[k].resize(validsum);
	i=0;
	for( uint n=0; n<matrix.index.size(); n++ ) {
		if( matrix(n,n) ) {
			for( uint j=0; j<matrix.index[n].size(); j++ ) {
				if( matrix.index[n][j] >= matrix.index.size() ) {
					email::print( "Compress faild !\n" );
					email::send();
				}
				int mapped_index = index_map[matrix.index[n][j]];
				if( mapped_index >= 0 )	final_matrix.set_element(i,mapped_index,matrix.value[n][j]);
			}
			for( uint k=0; k<num; k++ ) final_rhs[k][i] = rhs[k][n];
			i++;
		}
	}
	matrix = final_matrix;
	for( uint k=0; k<num; k++ ) rhs[k] = final_rhs[k];
	delete [] final_rhs;
}

template<class T>
static void compress( SparseMatrix<T> &matrix, std::vector<T> &rhs, std::vector<int> &index_map ) {
	compress<T>(matrix,&rhs,1,index_map);
}

// result = A * B
static bool multiply( SparseMatrix<FLOAT64> &result, const SparseMatrix<FLOAT64> &A, const SparseMatrix<FLOAT64> &B ) {
	result.clear();
	result.resize(A.index.size());
	for( uint i=0; i<A.index.size(); i++ ) {
		for( uint n=0; n<A.index[i].size(); n++ ) {
			uint j = A.index[i][n];
			FLOAT64 a = A.value[i][n];
			for( uint m=0; m<B.index[j].size(); m++ ) {
				uint bi = B.index[j][m];
				FLOAT64 b = B.value[j][m];
				result.add_to_element(i,bi,a*b);
			}
		}
	}
	return true;
}

// result = A^T
static void transpose( SparseMatrix<FLOAT64> &result, const SparseMatrix<FLOAT64> &A ) {
	// Proble maximum column of matrix A
	uint A_columns = 0;
	for( uint i=0; i<A.index.size(); i++ ) {
		for( uint j=0; j<A.index[i].size(); j++ ) {
			if( A_columns < A.index[i][j]+1 ) A_columns = A.index[i][j]+1;
		}
	}
	// Now begin transpose !
	result.clear();
	result.resize(A_columns);
	for( uint i=0; i<A.index.size(); i++ ) {
		for( uint j=0; j<A.index[i].size(); j++ ) {
			result.set_element(A.index[i][j],i,A.value[i][j]);
		}
	}
}

// result = diag(vec)
static void diag( SparseMatrix<FLOAT64> &result, const std::vector<FLOAT64> &vec ) {
	result.clear();
	result.resize(vec.size());
	for( uint i=0; i<vec.size(); i++ ) {
		result.set_element(i,i,vec[i]);
	}
}

static void write_matlab1d( std::ostream &output, const std::vector<FLOAT64> &v, const char *variable_name ) {
	output<<variable_name<<"=[\n";
	for( uint i=0; i<v.size(); i++ ) {
		output<<v[i]<<"\n";
	}
	output<<"];"<<std::endl;
}

template<class T>
static bool invert2x2( const T A[2][2], T result[2][2] ) {
	T det = A[0][0]*A[1][1]-A[0][1]*A[1][0];
	if( det ) {
		result[0][0] = A[1][1]/det;
		result[1][0] = -A[1][0]/det;
		result[1][1] = A[0][0]/det;
		result[0][1] = -A[0][1]/det;
	}
	return det;
}

template<class T>
static T determinant3x3( const T A[3][3] ) {
	return A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])
	-A[1][0]*(A[0][1]*A[2][2]-A[2][1]*A[0][2])
	+A[2][0]*(A[0][1]*A[1][2]-A[1][1]*A[0][2]);
}

template<class T>
static bool invert3x3( const T A[3][3], T result[3][3] ) {
	T determinant = determinant3x3(A);
	if ( ! determinant ) return false;
	T invdet = 1.0/determinant;
	result[0][0] =  (A[1][1]*A[2][2]-A[1][2]*A[2][1])*invdet;
	result[1][0] = -(A[1][0]*A[2][2]-A[2][0]*A[1][2])*invdet;
	result[2][0] =  (A[1][0]*A[2][1]-A[2][0]*A[1][1])*invdet;
	result[0][1] = -(A[0][1]*A[2][2]-A[2][1]*A[0][2])*invdet;
	result[1][1] =  (A[0][0]*A[2][2]-A[2][0]*A[0][2])*invdet;
	result[2][1] = -(A[0][0]*A[2][1]-A[0][1]*A[2][0])*invdet;
	result[0][2] =  (A[0][1]*A[1][2]-A[0][2]*A[1][1])*invdet;
	result[1][2] = -(A[0][0]*A[1][2]-A[0][2]*A[1][0])*invdet;
	result[2][2] =  (A[0][0]*A[1][1]-A[0][1]*A[1][0])*invdet;
	return true;
}

// Original: gluInvertMatrix() copied from 
// http://stackoverflow.com/questions/1148309/inverting-a-4x4-matrix
template<class T>
static bool myInvertMatrix(const T m[16], T invOut[16]) {
	T inv[16], det; int i;
	inv[0] =   m[5]*m[10]*m[15] - m[5]*m[11]*m[14] - m[9]*m[6]*m[15]
	+ m[9]*m[7]*m[14] + m[13]*m[6]*m[11] - m[13]*m[7]*m[10];
	inv[4] =  -m[4]*m[10]*m[15] + m[4]*m[11]*m[14] + m[8]*m[6]*m[15]
	- m[8]*m[7]*m[14] - m[12]*m[6]*m[11] + m[12]*m[7]*m[10];
	inv[8] =   m[4]*m[9]*m[15] - m[4]*m[11]*m[13] - m[8]*m[5]*m[15]
	+ m[8]*m[7]*m[13] + m[12]*m[5]*m[11] - m[12]*m[7]*m[9];
	inv[12] = -m[4]*m[9]*m[14] + m[4]*m[10]*m[13] + m[8]*m[5]*m[14]
	- m[8]*m[6]*m[13] - m[12]*m[5]*m[10] + m[12]*m[6]*m[9];
	inv[1] =  -m[1]*m[10]*m[15] + m[1]*m[11]*m[14] + m[9]*m[2]*m[15]
	- m[9]*m[3]*m[14] - m[13]*m[2]*m[11] + m[13]*m[3]*m[10];
	inv[5] =   m[0]*m[10]*m[15] - m[0]*m[11]*m[14] - m[8]*m[2]*m[15]
	+ m[8]*m[3]*m[14] + m[12]*m[2]*m[11] - m[12]*m[3]*m[10];
	inv[9] =  -m[0]*m[9]*m[15] + m[0]*m[11]*m[13] + m[8]*m[1]*m[15]
	- m[8]*m[3]*m[13] - m[12]*m[1]*m[11] + m[12]*m[3]*m[9];
	inv[13] =  m[0]*m[9]*m[14] - m[0]*m[10]*m[13] - m[8]*m[1]*m[14]
	+ m[8]*m[2]*m[13] + m[12]*m[1]*m[10] - m[12]*m[2]*m[9];
	inv[2] =   m[1]*m[6]*m[15] - m[1]*m[7]*m[14] - m[5]*m[2]*m[15]
	+ m[5]*m[3]*m[14] + m[13]*m[2]*m[7] - m[13]*m[3]*m[6];
	inv[6] =  -m[0]*m[6]*m[15] + m[0]*m[7]*m[14] + m[4]*m[2]*m[15]
	- m[4]*m[3]*m[14] - m[12]*m[2]*m[7] + m[12]*m[3]*m[6];
	inv[10] =  m[0]*m[5]*m[15] - m[0]*m[7]*m[13] - m[4]*m[1]*m[15]
	+ m[4]*m[3]*m[13] + m[12]*m[1]*m[7] - m[12]*m[3]*m[5];
	inv[14] = -m[0]*m[5]*m[14] + m[0]*m[6]*m[13] + m[4]*m[1]*m[14]
	- m[4]*m[2]*m[13] - m[12]*m[1]*m[6] + m[12]*m[2]*m[5];
	inv[3] =  -m[1]*m[6]*m[11] + m[1]*m[7]*m[10] + m[5]*m[2]*m[11]
	- m[5]*m[3]*m[10] - m[9]*m[2]*m[7] + m[9]*m[3]*m[6];
	inv[7] =   m[0]*m[6]*m[11] - m[0]*m[7]*m[10] - m[4]*m[2]*m[11]
	+ m[4]*m[3]*m[10] + m[8]*m[2]*m[7] - m[8]*m[3]*m[6];
	inv[11] = -m[0]*m[5]*m[11] + m[0]*m[7]*m[9] + m[4]*m[1]*m[11]
	- m[4]*m[3]*m[9] - m[8]*m[1]*m[7] + m[8]*m[3]*m[5];
	inv[15] =  m[0]*m[5]*m[10] - m[0]*m[6]*m[9] - m[4]*m[1]*m[10]
	+ m[4]*m[2]*m[9] + m[8]*m[1]*m[6] - m[8]*m[2]*m[5];	
	det = m[0]*inv[0] + m[1]*inv[4] + m[2]*inv[8] + m[3]*inv[12];
	if (det == 0) return false;
	det = 1.0 / det;	
	for (i = 0; i < 16; i++) invOut[i] = inv[i] * det;	
	return true;
}

template<class T>
static bool invert4x4( const T A[4][4], T result[4][4] ) {
	T M[16];
	T Minv[16];
	for( uint i=0; i<4; i++ ) {
		for( uint j=0; j<4; j++ ) {
			M[i+4*j] = A[i][j];
		}
	}
	if( ! myInvertMatrix(M,Minv)) return false;
	for( uint i=0; i<4; i++ ) {
		for( uint j=0; j<4; j++ ) {
			result[i][j] = Minv[i+4*j];
		}
	}
	return true;
}

#endif