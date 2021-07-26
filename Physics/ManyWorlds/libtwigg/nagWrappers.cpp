#include "stdafx.h"

#include "twigg/nagWrappers.h"
#include "twigg/blockmat.h"
#include "twigg/linalg.h"

#include <nag.h>
#include <nag_stdlib.h>
#include <nagf04.h>
#include <nagf11.h>
#include <nagd06.h>

#include <boost/bind.hpp>
#include "twigg/util.h"

using namespace std;

SymmetricCoordinateStorage::SymmetricCoordinateStorage( Integer n, Integer reserve )
{
	n_ = n;
	nnz_ = 0;
	nnzc_ = 0;
	la_ = reserve;

	init();
}

SymmetricCoordinateStorage::SymmetricCoordinateStorage( const BlockMatrix& mat, Integer in_nnzc )
{
	n_ = mat.nrows();
	nnz_ = ((mat.nnz() * 9) - n_) / 2 + n_;
	nnzc_ = 0;
	la_ = nnz_ + in_nnzc;

	init();
	mat.toSymmetricCoordinateStorage( irow_.ptr, icol_.ptr, data_.ptr, &rowPtrs_[0] );
}

SymmetricCoordinateStorage::SymmetricCoordinateStorage( const BlockMatrix& mat, 
													   SymmetricCoordinateStorage::PrecondPtr preconditioner )
{
	n_ = mat.nrows();
	nnz_ = ((mat.nnz() * 9) - n_) / 2 + n_;
	nnzc_ = 0;
	la_ = nnz_ + std::max<Integer>( nnz_ + nnz_/10, 
		(preconditioner.get() != 0 ? preconditioner->data.size() : 0) );

	init();
	mat.toSymmetricCoordinateStorage( irow_.ptr, icol_.ptr, data_.ptr, &rowPtrs_[0] );

	if( preconditioner.get() != 0 )
	{
		ipiv_ = preconditioner->ipiv;
		istr_.reserve( preconditioner->istr.size() );
		std::transform( 
			preconditioner->istr.begin(), 
			preconditioner->istr.end(), 
			back_inserter(istr_),
			boost::bind( std::plus<Integer>(), _1, nnz_ ) );

		std::copy( preconditioner->data.begin(), preconditioner->data.end(), data_.ptr + nnz_ );
		std::copy( preconditioner->irow.begin(), preconditioner->irow.end(), irow_.ptr + nnz_ );
		std::copy( preconditioner->icol.begin(), preconditioner->icol.end(), icol_.ptr + nnz_ );
		nnzc_ = preconditioner->data.size();

		preconditioner->use();
		preconditionerBuilt_ = true;
	}
}

void SymmetricCoordinateStorage::init()
{
	data_.reset( la_ );
	irow_.reset( la_ );
	icol_.reset( la_ );
	rowPtrs_.resize( n_ + 1 );

	preconditionerBuilt_ = false;
}

void SymmetricCoordinateStorage::append( Integer row, Integer col, double value )
{
	// Convert to 1-indexed.
	++row;
	++col;

	assert( nnz_ < la_ );
	assert( nnzc_ == 0 );
	assert( nnz_ == 0 || (row > irow_.ptr[ nnz_ - 1 ]) || (row == irow_.ptr[ nnz_ - 1 ] && col > icol_.ptr[ nnz_ - 1 ] ) );
	assert( row >= col );

	if( nnz_ != 0 )
	{
		for( Integer iRow = irow_.ptr[ nnz_ - 1 ] + 1; iRow <= row; ++iRow )
			rowPtrs_[ iRow - 1 ] = nnz_;
	}

	irow_.ptr[ nnz_ ] = row;
	icol_.ptr[ nnz_ ] = col;
	data_.ptr[ nnz_ ] = value;

	++nnz_;
	rowPtrs_[ row ] = nnz_;
}

fortran_matrix SymmetricCoordinateStorage::toMatrix() const
{
	// The matrix must be square since it's symmetric.
	fortran_matrix result( n_, n_ );
	for( Integer i = 0; i < nnz_; ++i )
	{
		result( irow_.ptr[i] - 1, icol_.ptr[i] - 1 ) = data_.ptr[i];
		result( icol_.ptr[i] - 1, irow_.ptr[i] - 1 ) = data_.ptr[i];
	}
	return result;
}

void SymmetricCoordinateStorage::clearPreconditioner()
{
	this->preconditionerBuilt_ = false;
}


CoordinateStorage::CoordinateStorage( const BlockMatrix& mat, Integer reserve )
{
	n_ = mat.nrows();
	nnz_ = mat.nnz() * 9;
	nnzc_ = 0;

	la_ = nnz_ + reserve;

	init();
	mat.toCoordinateStorage( irow_.ptr, icol_.ptr, data_.ptr );
}

CoordinateStorage::CoordinateStorage( const BlockMatrix& mat, PrecondPtr preconditioner )
{
	n_ = mat.nrows();
	nnz_ = mat.nnz() * 9;
	nnzc_ = 0;

	la_ = nnz_ + std::max<long>( nnz_ + nnz_/10, 
		(preconditioner.get() != 0 ? preconditioner->data.size() : 0) );

	init();
	mat.toCoordinateStorage( irow_.ptr, icol_.ptr, data_.ptr );

	if( preconditioner.get() != 0 )
	{
		ipivp_ = preconditioner->ipivp;
		ipivq_ = preconditioner->ipivq;

		istr_.reserve( preconditioner->istr.size() );
		std::transform( 
			preconditioner->istr.begin(), 
			preconditioner->istr.end(), 
			back_inserter(istr_),
			boost::bind( std::plus<Integer>(), _1, nnz_ ) );

		idiag_.reserve( preconditioner->idiag.size() );
		std::transform( 
			preconditioner->idiag.begin(), 
			preconditioner->idiag.end(), 
			back_inserter(idiag_),
			boost::bind( std::plus<Integer>(), _1, nnz_ ) );

		std::copy( preconditioner->data.begin(), preconditioner->data.end(), data_.ptr + nnz_ );
		std::copy( preconditioner->irow.begin(), preconditioner->irow.end(), irow_.ptr + nnz_ );
		std::copy( preconditioner->icol.begin(), preconditioner->icol.end(), icol_.ptr + nnz_ );
		nnzc_ = preconditioner->data.size();

		preconditioner->use();
		preconditionerBuilt_ = true;
	}
}

CoordinateStorage::CoordinateStorage( Integer n, Integer reserve )
{
	n_ = n;
	nnz_ = 0;
	nnzc_ = 0;
	la_ = reserve;

	init();
}

void CoordinateStorage::init()
{
	data_.reset( la_ );
	irow_.reset( la_ );
	icol_.reset( la_ );

	preconditionerBuilt_ = false;
}

fortran_matrix CoordinateStorage::toMatrix() const
{
	// The matrix must be square since it's symmetric.
	fortran_matrix result( n_, n_ );
	for( Integer i = 0; i < nnz_; ++i )
		result( irow_.ptr[i] - 1, icol_.ptr[i] - 1 ) = data_.ptr[i];
	return result;
}

void CoordinateStorage::append( Integer row, Integer col, double value )
{
	// Convert to 1-indexed.
	++row;
	++col;

	assert( nnz_ < la_ );
	assert( nnzc_ == 0 );
	assert( (nnz_ == 0) || (row > irow_.ptr[ nnz_ - 1 ]) || (row == irow_.ptr[ nnz_ - 1 ] && (col > icol_.ptr[ nnz_ - 1 ])) );

	irow_.ptr[ nnz_ ] = row;
	icol_.ptr[ nnz_ ] = col;
	data_.ptr[ nnz_ ] = value;

	++nnz_;
}

NAGIncompleteLUPreconditioner::NAGIncompleteLUPreconditioner( const CoordinateStorage& a )
	:	ipivp( a.ipivp_.begin(), a.ipivp_.end() ),
		ipivq( a.ipivq_.begin(), a.ipivq_.end() ),
		irow( a.irow_.ptr + a.nnz(), a.irow_.ptr + a.nnz() + a.nnzc() ),
		icol( a.icol_.ptr + a.nnz(), a.icol_.ptr + a.nnz() + a.nnzc() ),
		data( a.data_.ptr + a.nnz(), a.data_.ptr + a.nnz() + a.nnzc() )
{
	istr.reserve( a.n() + 1 );
	std::transform( a.istr_.begin(), a.istr_.end(),
		std::back_inserter( istr ), 
		boost::bind( std::minus<Integer>(), _1, a.nnz() ) );

	idiag.reserve( a.n() );
	std::transform( a.idiag_.begin(), a.idiag_.end(),
		std::back_inserter( idiag ), 
		boost::bind( std::minus<Integer>(), _1, a.nnz() ) );
}

NAGIncompleteLUPreconditioner::~NAGIncompleteLUPreconditioner()
{
}

NAGCholeskyPreconditioner::NAGCholeskyPreconditioner( const SymmetricCoordinateStorage& a )
	:	ipiv( a.ipiv_.begin(), a.ipiv_.end() ),
		irow( a.irow_.ptr + a.nnz(), a.irow_.ptr + a.nnz() + a.nnzc() ),
		icol( a.icol_.ptr + a.nnz(), a.icol_.ptr + a.nnz() + a.nnzc() ),
		data( a.data_.ptr + a.nnz(), a.data_.ptr + a.nnz() + a.nnzc() )
{
	istr.reserve( a.n() + 1 );
	std::transform( a.istr_.begin(), a.istr_.end(),
		std::back_inserter( istr ), 
		boost::bind( std::minus<Integer>(), _1, a.nnz() ) );
}

NAGCholeskyPreconditioner::~NAGCholeskyPreconditioner()
{
}

BasicPreconditioner::~BasicPreconditioner()
{
}

Integer nag_buildPreconditioner( SymmetricCoordinateStorage& A, int lfill )
{
	// Recommendations from the NAG manual
	const double dtol = 0.0;

	// Recommendations from the NAG manual
	NagError fail;
	INIT_FAIL(fail);
	Integer npivm = 100;
	const double dscale = 0.0;

	A.ipiv_.resize( A.n() );
	A.istr_.resize( A.n() + 1 );

	nag_sparse_sym_chol_fac(
		A.n(),
		A.nnz(),
		&A.data_.ptr,
		&A.la_,
		&A.irow_.ptr,
		&A.icol_.ptr,
		lfill,
		dtol,
		Nag_SparseSym_ModFact,
		dscale,
		Nag_SparseSym_MarkPiv,
		&A.ipiv_[0],
		&A.istr_[0],
		&A.nnzc_,
		&npivm,
		&A.comm_,
		&fail);

	if( fail.code != NE_NOERROR )
	{
		throw SolverException("NAG partial Cholesky decomposition failed: "
			+ std::string(fail.message) );
	}

	A.preconditionerBuilt_ = true;
	return npivm;
}

Integer nag_buildPreconditioner( CoordinateStorage& A, int lfill )
{
	// Recommendations from the NAG manual
	const double dtol = 0.001;

	A.ipivp_.resize( A.n() );
	A.ipivq_.resize( A.n() );
	A.istr_.resize( A.n() + 1 );
	A.idiag_.resize( A.n() );

	Integer npivm = 100;

	NagError fail;
	INIT_FAIL(fail);
	nag_sparse_nsym_fac(
		A.n(),
		A.nnz(),
		&A.data_.ptr,
		&A.la_,
		&A.irow_.ptr,
		&A.icol_.ptr,
		lfill,				// lfill
		dtol,				// dtol
		Nag_SparseNsym_CompletePiv,		// recommended value
		Nag_SparseNsym_UnModFact,		// don't know what this should be
		&A.ipivp_[0],
		&A.ipivq_[0],
		&A.istr_[0],
		&A.idiag_[0],
		&A.nnzc_,
		&npivm,
		&fail );

	if( fail.code != NE_NOERROR )
		throw LinearAlgebraException("NAG ILU factorization failed: "
			+ std::string(fail.message));

	return npivm;
}

std::vector<double> nag_solve( SymmetricCoordinateStorage& A, double* b, Integer& itn, double& rnorm )
{
	std::vector<double> x( A.n(), 0.0 );
	const int max_it = 200;
	const double ksp_atol = 1e-12;

	NagError fail;
	INIT_FAIL(fail);
	nag_sparse_sym_chol_sol( 
		Nag_SparseSym_CG,
		A.n(),
		A.nnz(),
		A.data_.ptr,
		A.la(),
		A.irow_.ptr,
		A.icol_.ptr,
		&(A.ipiv_[0]),
		&(A.istr_[0]),
		b,
		ksp_atol,
		max_it,
		&x[0],
		&rnorm,
		&itn,
		&A.comm_,
		&fail);

	if( fail.code != NE_NOERROR )
	{
		std::ostringstream oss;
		oss << "NAG Conjugate Gradient solver failed to converge in " << max_it << " iterations: "
			<< fail.message << std::endl
			<< "   Final residual: " << rnorm;
		throw LinearAlgebraException(oss.str());
	}

	return x;
}

std::vector<double> nag_solve( CoordinateStorage& A, double* b, Integer& itn, double& rnorm )
{
	std::vector<double> x( A.n(), 0.0 );
	const int max_it = 100;
	const double ksp_atol = 1e-12;

	NagError fail;
	INIT_FAIL(fail);
	nag_sparse_nsym_fac_sol( 
		Nag_SparseNsym_RGMRES,
		A.n(),
		A.nnz(),
		A.data_.ptr,
		A.la(),
		A.irow_.ptr,
		A.icol_.ptr,
		&(A.ipivp_[0]),
		&(A.ipivq_[0]),
		&A.istr_[0],
		&A.idiag_[0],
		b,
		10,			// restart subspace dimension
		ksp_atol,
		max_it,
		&x[0],
		&rnorm,
		&itn,
		&A.comm_,
		&fail );

	if( fail.code != NE_NOERROR )
		throw LinearAlgebraException("NAG GMRES solver failed to converge: "
			+ std::string(fail.message));

	return x;
}



/*
#define EDGE(I,J) edge.ptr[3*((J)-1)+(I)-1]
#define CONN(I,J) conn.ptr[3*((J)-1)+(I)-1]
#define COOR(I,J) coor.ptr[2*((J)-1)+(I)-1]

std::vector<Triangle> triangulate( const std::vector<vl::Vec2>& vertices )
{
	const Integer nvb = vertices.size();		// boundary verts
	const Integer nvint = 0;					// no interior vertices
	const Integer nvmax = 6000;					// no generated vertices
	const Integer nedge = nvb;					// # bdy edges

	NAG_scoped_array<Integer> edge( 3*nedge );
	for( int j = 1; j <= nvb; ++j )
	{
		EDGE(1, j) = j;
		EDGE(2, j) = j%nvb + 1;
		EDGE(3, j) = 0;
	}

	NAG_scoped_array<double> coor( 2*nvmax );
	for( unsigned int j = 1; j <= nvb; ++j )
	{
		COOR(1, j) = (vertices[j-1])[0];
		COOR(2, j) = (vertices[j-1])[1];
	}

	NAG_scoped_array<double> weight( std::max<Integer>(1, nvint) );
	weight.ptr[0] = 0.01;

	NAG_scoped_array<Integer> conn( 3*(2*nvmax + 5) );

	Integer npropa = 5; // not sure yet

	Integer nv = 0;
	Integer nelt = 0;

	Integer itrace = 0;
	char* outfile = 0;

	{
		std::ofstream ofs("test.d");
		ofs << "d06abc Example Program Data\n";
		ofs << nvb << "  " << nedge << "  :NVB NEDGE\n";
		
		ofs.precision(10);
		ofs.setf(std::ios_base::scientific);

		{
			for (unsigned int i = 1; i <= nvb; ++i)
				ofs << i << " " << COOR(1,i) << " " << COOR(2,i) << "\n";
		}

		{
			for (unsigned int i = 1; i <= nedge; ++i)
				ofs << i << " " << EDGE(1,i) << " " << EDGE(2,i) << " " << EDGE(3,i) << "\n";
		}

		ofs << "'N'               :Printing option 'Y' or 'N'" << std::endl;
	}

	NagError fail;
	INIT_FAIL(fail);

	nag_mesh2d_delaunay(
		nvb,
		nvint,
		nvmax,
		nedge,
		edge.ptr,
		&nv,
		&nelt,
		coor.ptr,
		conn.ptr,
		weight.ptr,
		npropa,
		itrace,
		outfile, 
		&fail );

	if( fail.code != NE_NOERROR )
		throw LinearAlgebraException("Delauney triangulation failed: "
			+ std::string(fail.message));

	std::vector<Triangle> result;
	for( int i = 0; i < nelt; ++i )
		result.push_back( Triangle( CONN(i, 0), CONN(i, 1), CONN(i, 2) ) );
	return result;
}
*/
