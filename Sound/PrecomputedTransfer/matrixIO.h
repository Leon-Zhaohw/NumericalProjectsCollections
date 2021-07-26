#ifndef _MATRIXIO_H_
#define _MATRIXIO_H_

/*

Jernej Barbic
Carnegie Mellon University
2003-2006

A small matrix I/O library.

All matrices loaded/written to disk are stored in the following **binary** format:
<number of rows> (signed integer, 4 bytes)
<number of columns> (signed integer, 4 bytes)
data in double precision, in COLUMN-major order (first column first, then second column, and so on)

All matrices are loaded into 1-D C memory arrays, in COLUMN-major order.
                                                                                                                                                             
A m x n matrix will occupy 2*sizeof(int) + m*n*sizeof(double) bytes on disk.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// LAPACK style column-major order matrices
#define ELT(rows,i,j) (j)*(rows)+(i)
#define MIN(x,y) ((x)<(y) ? (x) : (y))
#define MAX(x,y) ((x)>(y) ? (x) : (y))

// read a matrix from 'filename' into main memory
// m, n and matrix are output parameters
// matrix need not be pre-allocated, simply pass the address of a pointer,
// and the routine will allocate space for the matrix
template <class real>
int ReadMatrixFromDisk(char* filename, int * m, int * n, real ** matrix);

// write a m x n main memory matrix to the disk file 'filename'
template <class real>
int WriteMatrixToDisk(char* filename, int m, int n, real * matrix);

// auxiliary routines:
template <class real>
int WriteMatrixToStream(FILE * file, int m, int n, real * matrix);

int WriteMatrixHeaderToStream(FILE * file, int m, int n);

template <class real>
int ReadMatrixFromStream(FILE * file, int m, int n, real * matrix);

int ReadMatrixSizeFromStream(FILE * file, int * m, int * n);

// the following functions are wrappers that will cause the program to exit in case of failure:
template <class real>
int WriteMatrixToDisk_(char* filename, int m, int n, real * matrix);

template <class real>
int ReadMatrixFromDisk_(char* filename, int * m, int * n, real ** matrix);

int OpenFile_(char * filename, FILE ** fin, char * mode);

// if (a!=b) exit program
int Assert_(int a, int b, int positionIdentifier);
int Assert_(bool cond, int positionIdentifier);

#endif

