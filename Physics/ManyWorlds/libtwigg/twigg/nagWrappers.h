#ifndef __NAGWRAPPERS_H__

#define __NAGWRAPPERS_H__

#include "twigg/linalg.h"
#include "twigg/primitive.h"

#include <iostream>
#include <nag.h>
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>
#include <vector>
#include <nag_stdlib.h>

// This is pretty ugly, I admit.  The problem is
// that the pointer itself has to be accessible to the 
// NAG library for it to work properly.
template <typename T>
struct NAG_scoped_array
{
	NAG_scoped_array()
	{
		ptr = 0;
	}

	NAG_scoped_array( unsigned int size )
	{
		ptr = NAG_ALLOC( size, T );
	}

	~NAG_scoped_array()
	{
		if( ptr != 0 )
			NAG_FREE( ptr );
	}

	void reset( unsigned int size )
	{
		// Not sure if NAG_FREE(0) is well-defined or not.
		if( ptr != 0 )
			NAG_FREE( ptr );

		ptr = NAG_ALLOC( size, T );
	}

	T* ptr;
};

class BlockMatrix;
class NAGCholeskyPreconditioner;
class NAGIncompleteLUPreconditioner;

template <typename StorageType>
void dumpSparse( const StorageType& mat, std::ostream& out )
{
	for( unsigned int i = 0; i < mat.nnz(); ++i )
	{
		out << "(" << mat.irow(i) << ", " << mat.icol(i) << ") <- " << (mat.data())[i] << std::endl;
	}
}

class SymmetricCoordinateStorage
{
public:
	typedef NAGCholeskyPreconditioner PrecondType;

private:
	typedef boost::shared_ptr<PrecondType> PrecondPtr;
	friend class NAGCholeskyPreconditioner;

	friend Integer nag_buildPreconditioner( SymmetricCoordinateStorage&, int );
	friend std::vector<double> nag_solve( SymmetricCoordinateStorage&, double*, Integer&, double& );

public:
	SymmetricCoordinateStorage( const BlockMatrix& mat, PrecondPtr precond );
	SymmetricCoordinateStorage( const BlockMatrix& mat, Integer reserve = 0 );
	SymmetricCoordinateStorage( Integer n, Integer reserve );

	const double* data() const						{ return data_.ptr; }
	double* data()									{ return data_.ptr; }
	Integer irow( unsigned int i ) const			{ return irow_.ptr[i]; }
	Integer icol( unsigned int i ) const			{ return icol_.ptr[i]; }
	unsigned int rowPtrs( unsigned int i ) const	{ return rowPtrs_[i]; }

	Integer n() const								{ return n_; }
	Integer nnz() const								{ return nnz_; }
	Integer nnzc() const							{ return nnzc_; }
	Integer la() const								{ return la_; }

	void append( Integer row, Integer col, double value );

	fortran_matrix toMatrix() const;

	bool preconditionerBuilt() const { return preconditionerBuilt_; }
	void clearPreconditioner();

private:
	// allocate all the appropriate arrays.
	void init();

	Integer n_;
	Integer nnz_;
	Integer nnzc_;
	Integer la_;
	NAG_scoped_array<double> data_;
	NAG_scoped_array<Integer> irow_;
	NAG_scoped_array<Integer> icol_;
	std::vector<unsigned int> rowPtrs_;

	std::vector<Integer> ipiv_;
	std::vector<Integer> istr_;

	Nag_Sparse_Comm comm_;

	bool preconditionerBuilt_;
};

class CoordinateStorage
{
public:
	typedef NAGIncompleteLUPreconditioner PrecondType;

private:
	typedef boost::shared_ptr<PrecondType> PrecondPtr;
	friend class NAGIncompleteLUPreconditioner;

	friend Integer nag_buildPreconditioner( CoordinateStorage&, int );
	friend std::vector<double> nag_solve( CoordinateStorage&, double*, Integer&, double& );

public:
	CoordinateStorage( const BlockMatrix& mat, Integer reserve = 0 );
	CoordinateStorage( const BlockMatrix& mat, PrecondPtr precond );
	CoordinateStorage( Integer n, Integer reserve );

	const double* data() const						{ return data_.ptr; }
	double* data()									{ return data_.ptr; }
	Integer irow( unsigned int i ) const			{ return irow_.ptr[i]; }
	Integer icol( unsigned int i ) const			{ return icol_.ptr[i]; }

	Integer n() const								{ return n_; }
	Integer nnz() const								{ return nnz_; }
	Integer nnzc() const							{ return nnzc_; }
	Integer la() const								{ return la_; }

	fortran_matrix toMatrix() const;
	bool preconditionerBuilt() const { return preconditionerBuilt_; }

	void append( Integer row, Integer col, double value );

private:
	void init();

	Integer n_;
	Integer nnz_;
	Integer nnzc_;
	Integer la_;
	NAG_scoped_array<double> data_;
	NAG_scoped_array<Integer> irow_;
	NAG_scoped_array<Integer> icol_;

	std::vector<Integer> istr_;
	std::vector<Integer> idiag_;
	std::vector<Integer> ipivp_;
	std::vector<Integer> ipivq_;

	Nag_Sparse_Comm comm_;
	bool preconditionerBuilt_;
};

class BasicPreconditioner
{
private:
	unsigned int age_;

public:
	BasicPreconditioner()
		: age_(1) {}
	virtual ~BasicPreconditioner();

	unsigned int age() const
	{
		return age_;
	}

	void use()
	{
		++age_;
	}
};

class NAGIncompleteLUPreconditioner
	: public BasicPreconditioner
{
public:
	NAGIncompleteLUPreconditioner( const CoordinateStorage& a );
	virtual ~NAGIncompleteLUPreconditioner();

	const std::vector<Integer> ipivp;
	const std::vector<Integer> ipivq;
	std::vector<Integer> istr;
	std::vector<Integer> idiag;

	const std::vector<Integer> irow;
	const std::vector<Integer> icol;
	const std::vector<double> data;
};

class NAGCholeskyPreconditioner
	: public BasicPreconditioner
{
public:
	NAGCholeskyPreconditioner( const SymmetricCoordinateStorage& a );
	virtual ~NAGCholeskyPreconditioner();

	const std::vector<Integer> ipiv;
	std::vector<Integer> istr;

	const std::vector<Integer> irow;
	const std::vector<Integer> icol;
	const std::vector<double> data;
};

// Note that this also modifies the A matrix
Integer nag_buildPreconditioner( SymmetricCoordinateStorage& A, int lfill = 0 );
Integer nag_buildPreconditioner( CoordinateStorage& A, int lfill = 0 );

std::vector<double> nag_solve( SymmetricCoordinateStorage& A, double* b, Integer& itn, double& rnorm );
std::vector<double> nag_solve( CoordinateStorage& A, double* b, Integer& itn, double& rnorm );

#endif
