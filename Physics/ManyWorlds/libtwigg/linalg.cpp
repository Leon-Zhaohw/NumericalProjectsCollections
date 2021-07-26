#include "stdafx.h"

#include "twigg/linalg.h"
#include "twigg/random.h"

#ifdef NAG
#include "twigg/nagWrappers.h"
#endif

#include <boost/bind.hpp>
#include <boost/tuple/tuple.hpp>

#include <sstream>
#include <ctime>
#include <algorithm>

template <typename T>
int int_cast( T val )
{
    return boost::numeric_cast<int>( val );
}

#ifdef CONDOR
#include <nag.h>
#include <nag_stdlib.h>

#include <nagf06.h>
#include <nagf07.h>
#include <nagf08.h>
#include <nagf16.h>

#define cblas_dgemm f06yac
#define cblas_dgemv f06pac
#define CblasNoTrans NoTranspose
#define CblasTrans Transpose

#else
#include <mkl_lapack.h>
#include <mkl_vml.h>
#include <mkl_cblas.h>
#endif

#ifdef HDF5
#include "hdf5.h"
#endif

template <typename Matrix>
void addEpsilon( Matrix& matrix )
{
	typedef typename Matrix::value_type real_type;
	std::for_each( matrix.data(),
		matrix.data() + (matrix.nrows()*matrix.ncols()),
		boost::bind( std::plus<real_type>(), _1, std::numeric_limits<real_type>::epsilon() ) );
}

#ifndef _WIN32
NonNegativeMatrixFactorization::NonNegativeMatrixFactorization( const fortran_matrix& V, const unsigned int r, const unsigned int maxIter  )
	: W_( V.nrows(), r ), H_( r, V.ncols() )
{
	{
		// We initialize by generating a random matrix W and then fitting
		// H_ to it using NNLS.
		RandomStream rnd;
		rnd.uniform( 0.0, 1.0, W_.nrows()*W_.ncols(), W_.data() );
//		rnd.uniform( 0.0, 1.0, H_.nrows()*H_.ncols(), H_.data() );
		for( unsigned int iCol = 0; iCol < H_.ncols(); ++iCol )
		{
			fortran_matrix tempW = W_;
			std::vector<double> b( V.data() + iCol*V.nrows(), V.data() + (iCol+1)*V.nrows() );
			std::vector<double> x = leastSquares_NNLS( tempW, &b[0] );
			assert( x.size() == r );
			std::copy( x.begin(), x.end(), H_.data() + iCol*H_.nrows() );
		}
	}

	// Now, we apply our iterative procedure to find the correct factorization
	const double tol = 1e-6;
	unsigned int iter = 0;
	double prevErr = 0.0;
	while( true )
	{
#ifdef DUMPALL
		MATFile matFile( "testDump.mat", "test dump" );
		matFile.add( "H_prev", H_ );
		matFile.add( "W_prev", W_ );
		matFile.add( "V", V );
#endif

		++iter;

		{
			fortran_matrix H_numerator = prod( transposed(W_), V );
			addEpsilon( H_numerator );

			fortran_matrix WtimesW = prod( transposed(W_), W_ );
			fortran_matrix H_denominator = prod( WtimesW, H_ );
			addEpsilon( H_denominator );

			assert( H_numerator.nrows() == H_denominator.nrows() );
			assert( H_numerator.ncols() == H_denominator.ncols() );

			fortran_matrix H_update( H_.nrows(), H_.ncols() );
			assert( H_update.ncols() == H_denominator.ncols() );
			assert( H_update.nrows() == H_denominator.nrows() );

#ifdef NO_BLAS
			std::transform( H_numerator.data(), 
				H_numerator.data() + (H_.nrows()*H_.ncols()),
				H_denominator.data(),
				H_update.data(), 
				std::divides<double>() );
#else
			vdDiv( H_update.nrows() * H_update.ncols(), 
				H_numerator.data(), 
				H_denominator.data(), 
				H_update.data() );
#endif

			std::transform( H_.data(), H_.data() + (H_.nrows()*H_.ncols()),
				H_update.data(), H_.data(), std::multiplies<double>() );

#ifdef DUMPALL
			matFile.add( "H_numerator", H_numerator );
			matFile.add( "H_denominator", H_denominator );
			matFile.add( "H_update", H_update );
			matFile.add( "WtimesW", WtimesW );
#endif
		}

		{
			fortran_matrix W_numerator = prod( V, transposed(H_) );
			addEpsilon( W_numerator );

			fortran_matrix HtimesH = prod( H_, transposed(H_) );
			fortran_matrix W_denominator = prod( W_, HtimesH );
			addEpsilon( W_denominator );

			assert( W_numerator.nrows() == W_denominator.nrows() );
			assert( W_numerator.ncols() == W_denominator.ncols() );

			fortran_matrix W_update( W_.nrows(), W_.ncols() );
			assert( W_update.ncols() == W_denominator.ncols() );
			assert( W_update.nrows() == W_denominator.nrows() );

#ifdef NO_BLAS
			std::transform( 
				W_numerator.data(),
				W_numerator.data() + W_update.nrows() * W_update.ncols(),
				W_denominator.data(),
				W_update.data()
				std::divides<double>() );
#else
			vdDiv( W_update.nrows() * W_update.ncols(), 
				W_numerator.data(), 
				W_denominator.data(), 
				W_update.data() );
#endif

			std::transform( W_.data(), W_.data() + (W_.nrows()*W_.ncols()),
				W_update.data(), W_.data(), std::multiplies<double>() );

#ifdef DUMPALL
			matFile.add( "W_numerator", W_numerator );
			matFile.add( "W_denominator", W_denominator );
			matFile.add( "W_update", W_update );
			matFile.add( "HtimesH", HtimesH );
#endif
		}


#ifdef DUMPALL
		matFile.add( "H_next", H_ );
		matFile.add( "W_next", W_ );
#endif


/*
		double maxUpdate = 0.0;
		for( unsigned int iCol = 0; iCol < H_update.ncols(); ++iCol )
			for( unsigned int iRow = 0; iRow < H_update.nrows(); ++iRow )
				maxUpdate = std::max<double>( maxUpdate, 
					fabs( H_update(iRow, iCol) - 1.0 ) );

		for( unsigned int iCol = 0; iCol < W_update.ncols(); ++iCol )
			for( unsigned int iRow = 0; iRow < W_update.nrows(); ++iRow )
				maxUpdate = std::max<double>( maxUpdate, 
					fabs( W_update(iRow, iCol) - 1.0 ) );
*/
			
/*
		double maxUpdate = std::max<double>(
			std::max<double>(
				fabs(norminf( H_update ) - 1.0),
				fabs(absmin( H_update ) - 1.0) ),
			std::max<double>(
				fabs(norminf( W_update ) - 1.0),
				fabs(absmin( W_update ) - 1.0) ) );
*/

		if( iter % 50 == 0 )
		{

			// now, figure out the error
			fortran_matrix V_est = prod(W_, H_);
			assert( V_est.nrows() == V.nrows() );
			assert( V_est.ncols() == V.ncols() );

			add( scaled( V, -1.0 ), V_est );
			double err = sqrt( norm2( V_est ) );
	
			std::cout << "iter: " << iter << ", error: " << err << std::endl;
			if( fabs(prevErr - err) < tol )
				break;

			prevErr = err;
		}

		if( maxIter != 0 && iter > maxIter )
			break;
		//std::cout << "iter: " << iter << ", maxUpdate: " << maxUpdate << std::endl;
	}
}

const fortran_matrix& NonNegativeMatrixFactorization::W()
{
	return W_;
}

const fortran_matrix& NonNegativeMatrixFactorization::H()
{
	return H_;
}
#endif

#ifndef NO_BLAS
SVD::SVD( const fortran_matrix& matrix )
{
    int m = boost::numeric_cast<int>( matrix.nrows() );
    int n = boost::numeric_cast<int>( matrix.ncols() );

	fortran_matrix* A;
	char jobu, jobvt;

	if( m >= n )
	{
		// This is the standard version of the SVD
		A = &U_;
		U_.resize( matrix.nrows(), matrix.ncols() );
		U_ = matrix;
		Vt_.resize( n, n );
		jobu = 'O';
		jobvt = 'A';
		S_.resize( n );
	}
	else
	{
		A = &Vt_;
		Vt_.resize( matrix.nrows(), matrix.ncols() );
		Vt_ = matrix;
		U_.resize( m, m );
		jobu = 'A';
		jobvt = 'O';
		S_.resize( m );
	}

    int lda = boost::numeric_cast<int>( A->nrows() );
    int ldu = boost::numeric_cast<int>( U_.nrows() );
    int ldvt = boost::numeric_cast<int>( Vt_.nrows() );
	std::vector<double> work(1);
	int lwork = -1;	// first query on size of work array
	int info = 0;
	
	{
		// first query on size of work array
		dgesvd( &jobu, 
				&jobvt, 
				&m, 
				&n, 
				A->data(), 
				&lda, 
				&S_[0],
				U_.data(),
				&ldu,
				Vt_.data(),
				&ldvt,
				&work[0],
				&lwork,
				&info );
	}
	{
		// Now, know how much workspace is desired.
		lwork = static_cast<int>(work[0]);
		work.resize( lwork );
		dgesvd( &jobu, 
				&jobvt, 
				&m, 
				&n, 
				A->data(), 
				&lda, 
				&S_[0],
				U_.data(),
				&ldu,
				Vt_.data(),
				&ldvt,
				&work[0],
				&lwork,
				&info );
	}

	if( info < 0 )
	{
		std::ostringstream oss;
		oss << "SVD failed due to invalid parameter #" 
			<< -info;
		throw LinearAlgebraException( oss.str() );
	}
	else if( info > 0 )
	{
		std::ostringstream oss;
		oss << "SVD did not converge; " << info << 
			" diagonals of the intermediate bidiagonal form did not converge to 0.";
		throw LinearAlgebraException( oss.str() );
	}
}

const fortran_matrix& SVD::U() const
{
	return U_;
}

const fortran_matrix& SVD::Vt() const
{
	return Vt_;
}

const std::vector<double>& SVD::S() const
{
	return S_;
}


void test_svd(std::ostream& out)
{
	static bool seeded = false;
	if( !seeded )
		srand( time(0) );

	const unsigned int maxSide = 7;

	const unsigned int m = rand() % maxSide + 2;
	const unsigned int n = rand() % maxSide + 2;

	fortran_matrix A( m, n );
	for( unsigned int i = 0; i < m; ++i )
		for( unsigned int j = 0; j < n; ++j )
			A(i, j) = rand() % 100;

	/*
	fortran_matrix A(2, 3);
	A(0, 0) = 1;
	A(0, 1) = 0;
	A(0, 2) = 0;
	A(1, 0) = 0;
	A(1, 1) = 0;
	A(1, 2) = 1;
	*/

	SVD svd(A);
	std::vector<double> S( svd.S() );
	printDense( out, A, "A" );
	printDense( out, svd.U(), "U" );
	printDense( out, svd.Vt(), "Vt" );
	printVector( out, &S[0], S.size(), "S" );
}
#endif // NO_BLAS

fortran_matrix prod( const fortran_matrix& left, const fortran_matrix& right )
{
	assert( left.ncols() == right.nrows() );
	// please please please use return value optimization!
	fortran_matrix result( left.nrows(), right.ncols() );
	std::fill( result.data(), result.data() + result.nrows()*result.ncols(), 0.0 );

#ifdef NO_BLAS
	// want to do it a column at a time for memory locality in both left and result
	for( unsigned int iResCol = 0; iResCol < result.ncols(); ++iResCol )
	{
		for( unsigned int jLeftCol = 0; jLeftCol < left.ncols(); ++jLeftCol )
		{
			const double rij = right( jLeftCol, iResCol );
			for( unsigned int kResRow = 0; kResRow < result.nrows(); ++kResRow )
			{
				result(kResRow, iResCol) += rij * left(jLeftCol, kResRow);
			}
		}
	}
#else
	cblas_dgemm(
		CblasColMajor, 
		CblasNoTrans, 
		CblasNoTrans, 
		int_cast( left.nrows() ),       // m
		int_cast( right.ncols() ),      // n
		int_cast( left.ncols() ),       // k
		1.0,                            // alpha
		left.data(),                    // a
		int_cast( left.nrows() ),       // lda
		right.data(),                   // b
		int_cast( right.nrows() ),      // ldb
		0.0,                            // beta
		result.data(),                  // c 
		int_cast( result.nrows() )      // ldc
		);
#endif

	return result;
}

fortran_matrix prod( const Transposed<fortran_matrix>& left, const fortran_matrix& right )
{
	assert( left.mat().nrows() == right.nrows() );
	// please please please use return value optimization!
	fortran_matrix result( left.mat().ncols(), right.ncols() );

#ifdef NO_BLAS
	// TODO This is going to suck for memory locality:
	for( unsigned int iResCol = 0; iResCol < result.ncols(); ++iResCol )
	{
		for( unsigned int jLeftCol = 0; jLeftCol < left.mat().ncols(); ++jLeftCol )
		{
			const double rij = right( jLeftCol, iResCol );
			for( unsigned int kResRow = 0; kResRow < result.nrows(); ++kResRow )
			{
				result(kResRow, iResCol) += rij * left.mat()(kResRow, jLeftCol);
			}
		}
	}
#else
	cblas_dgemm(
		CblasColMajor, 
		CblasTrans, 
		CblasNoTrans, 
		int_cast( left.mat().ncols() ),     // m
		int_cast( right.ncols() ),          // n
		int_cast( left.mat().nrows() ),     // k
		1.0,                                // alpha
		left.mat().data(),                  // a
		int_cast( left.mat().nrows() ),     // lda
		right.data(),                       // b
		int_cast( right.nrows() ),          // ldb
		0.0,                                // beta
		result.data(),                      // c 
		int_cast( result.nrows() )          // ldc
		);
#endif

	return result;
}

fortran_matrix prod( const fortran_matrix& left, const Transposed<fortran_matrix>& right )
{
	assert( left.ncols() == right.mat().ncols() );
	fortran_matrix result( left.nrows(), right.mat().nrows() );

#ifdef NO_BLAS
	// want to do it a column at a time for memory locality in both left and result
	for( unsigned int iResCol = 0; iResCol < result.ncols(); ++iResCol )
	{
		for( unsigned int jLeftCol = 0; jLeftCol < left.ncols(); ++jLeftCol )
		{
			const double rij = right.mat()( iResCol, jLeftCol );
			for( unsigned int kResRow = 0; kResRow < result.nrows(); ++kResRow )
			{
				result(kResRow, iResCol) += rij * left(jLeftCol, kResRow);
			}
		}
	}
#else
	cblas_dgemm(
		CblasColMajor, 
		CblasNoTrans,
		CblasTrans,
		int_cast( left.nrows() ),
		int_cast( right.mat().nrows() ),
		int_cast( left.ncols() ),
		1.0,
		left.data(),
		int_cast( left.nrows() ),
		right.mat().data(),
		int_cast( right.mat().nrows() ),
		0.0,
		result.data(),
		int_cast( result.nrows() )
		);
#endif

	return result;
}

void prod( const fortran_matrix& left, const double* rhs, double* out )
{
	if( left.nrows() == 0 || left.ncols() == 0 )
	{
		std::fill( out, out + left.nrows(), 0.0 );
		return;
	}

#ifdef NO_BLAS
	std::fill( out, out + left.nrows(), 0.0 );
	for( unsigned int iCol = 0; iCol < left.ncols(); ++iCol )
	{
		double rhs_i = rhs[iCol];
		for( unsigned int jRow = 0; jRow < left.nrows(); ++jRow )
		{
			out[jRow] += rhs_i * left(jRow, iCol);
		}
	}
#else
	cblas_dgemv(
		CblasColMajor, 
		CblasNoTrans,
		int_cast( left.nrows() ),       // m
		int_cast( left.ncols() ),       // n
		1.0,                            // alpha
		left.data(),                    // a
		int_cast( left.nrows() ),       // lda
		rhs,                            // x
		1,                              // incx
		0,                              // beta (blas promises that if this is zero, then y need not be initialized)
		out,                            // y
		1                               // incy
		);
#endif
}

void prod( const Transposed<fortran_matrix>& left, const double* rhs, double* out )
{
	if( left.mat().nrows() == 0 || left.mat().ncols() == 0 )
	{
		std::fill( out, out + left.mat().ncols(), 0.0 );
		return;
	}

#ifdef NO_BLAS
	std::fill( out, out + left.mat().ncols(), 0.0 );
	for( unsigned int iRow = 0; iRow < left.mat().ncols(); ++iRow )
	{
		for( unsigned int jCol = 0; jCol < left.mat().nrows(); ++jCol )
		{
			out[iRow] += rhs[jCol] * left.mat()(jCol, iRow);
		}
	}
#else
	cblas_dgemv(
		CblasColMajor, 
		CblasTrans,
		int_cast( left.mat().nrows() ),     // m
		int_cast( left.mat().ncols() ),     // n
		1.0,                                // alpha
		left.mat().data(),                  // a
		int_cast( left.mat().nrows() ),     // lda
		rhs,                                // x
		1,                                  // incx
		0,                                  // beta (blas promises that if this is zero, then y need not be initialized)
		out,                                // y
		1                                   // incy
		);
#endif
}

std::vector<double> prod( const fortran_matrix& left, const double* rhs )
{
	std::vector<double> out( left.nrows() );
	prod( left, rhs, &out[0] );
	return out;
}

std::vector<double> prod( const Transposed<fortran_matrix>& left, const double* rhs )
{
	std::vector<double> out( left.mat().ncols() );
	prod( left, rhs, &out[0] );
	return out;
}

void prodAdd( const fortran_matrix& left, const Scaled< std::vector<double> >& x, std::vector<double>& y )
{
#ifdef NO_BLAS
	assert( left.ncols() == x.vec().size() );
	assert( left.nrows() == y.size() );

	// make sure this is in a register
	const double alpha = x.scale();

	// since it is column-major, it is most advantageous to walk through the columns
	for( unsigned int iCol = 0; iCol < left.ncols(); ++iCol )
	{
		const double xi = (x.vec())[iCol];
		for( unsigned int jRow = 0; jRow < left.nrows(); ++jRow )
			y[jRow] += alpha * xi * left(jRow, iCol);
	}
#else
	cblas_dgemv(  
		CblasColMajor, 
		CblasNoTrans,               // trans
		int_cast( left.nrows() ),   // m
		int_cast( left.ncols() ),   // n
		x.scale(),                  // alpha
		left.data(),                // a
		int_cast( left.nrows() ),   // lda
		&(x.vec())[0],              // x
		1,                          // incx
		1.0,                        // beta
		&y[0],                      // y
		1 );                        // incy
#endif
}

// y = a*y
// Uses the CBLAS dscal operator
void scale( double* vec, size_t n, double scale )
{
	if( n == 0 )
		return;

#ifdef NO_BLAS
	for( size_t i = 0; i < n; ++i )
		vec[i] *= scale;
#else
	cblas_dscal(
		int_cast( n ), 
		scale, 
		vec, 
		1 );
#endif
}

void scale( fortran_matrix& mat, double a )
{
	scale( mat.data(), mat.nrows() * mat.ncols(), a );
}

void scale( std::vector<double>& vec, double a )
{
	scale( &vec[0], vec.size(), a );
}

void add( const std::vector<double>& x, std::vector<double>& y )
{
	add( scaled(x, 1.0), y );
}

void add( const std::vector<float>& x, std::vector<float>& y )
{
	add( scaled(x, 1.0f), y );
}

void add( const fortran_matrix& x, fortran_matrix& y )
{
	add( scaled(x, 1.0), y );
}

// y = y + a*x
// equivalent to a CBLAS daxpy operator
void add( const Scaled<fortran_matrix>& x, fortran_matrix& y )
{
	assert( x.vec().nrows() == y.nrows() && x.vec().ncols() == y.ncols() );
	if( x.vec().nrows() == 0 || x.vec().ncols() == 0 )
		return;

#ifdef NO_BLAS
	for( size_t i = 0; i < x.vec().nrows() * x.vec().ncols(); ++i )
		(y.data())[i] += x.scale() * (x.vec().data())[i];
#else
	cblas_daxpy( int_cast( x.vec().nrows() * x.vec().ncols() ), 
		x.scale(),
		x.vec().data(), 
		1, 
		y.data(),
		1 );
#endif
}

void add( const Scaled< std::vector<double> >& x, std::vector<double>& y )
{
	assert( x.vec().size() == y.size() );
	if( x.vec().empty() )
		return;

#ifdef NO_BLAS
	for( size_t i = 0; i < y.size(); ++i )
		y[i] += x.scale() * (x.vec())[i];
#else
	cblas_daxpy( int_cast( x.vec().size() ), 
		x.scale(),
		&(x.vec())[0], 
		1, 
		&y[0],
		1 );
#endif
}

void add( const Scaled< std::vector<float> >& x, std::vector<float>& y )
{
	assert( x.vec().size() == y.size() );
	if( x.vec().empty() )
		return;

#ifdef NO_BLAS
	for( size_t i = 0; i < y.size(); ++i )
		y[i] += x.scale() * (x.vec())[i];
#else
	cblas_saxpy( int_cast( x.vec().size() ), 
		x.scale(),
		&(x.vec())[0], 
		1, 
		&y[0],
		1 );
#endif
}


#ifndef CONDOR
#ifndef _WIN32
extern "C" 
{
void nnls_( 
	double* A,
	int* MDA,
	int* M,
	int* N,
	double* B,
	double* X,
	double* RNORM,
	double* W,
	double* ZZ,
	int* INDEX,
	int* MODE);
}

std::vector<double> leastSquares_NNLS( fortran_matrix& A, const double* b_in )
{
	int m = A.nrows();
	int n = A.ncols();
	int mda = A.nrows();

	std::vector<double> b( b_in, b_in + m );
	std::vector<double> x( n );

	double rnorm;
	std::vector<double> w( n );
	std::vector<double> zz( m );
	std::vector<int> index( n );
	int mode;

	nnls_( 
		A.data(),
		&mda,
		&m,
		&n,
		&b[0],
		&x[0],
		&rnorm,
		&w[0],
		&zz[0],
		&index[0],
		&mode );

	if( mode == 1 )
	{
		// solution computed successfully
		return x;
	}
	else if( mode == 2 )
	{
		throw LinearAlgebraException( 
			"The dimensions of the problem are bad; either M < 0 or N < 0" );
	}
	else if( mode == 3 )
	{
		std::ostringstream oss;
		oss << "Iteration count exceeded; more than " << (3*n) << " iterations.";
		throw LinearAlgebraException( oss.str() );
	}
	else
	{
		std::ostringstream oss;
		oss << "Invalid MODE parameter returned from NNLS solver: " << mode;
		throw LinearAlgebraException( oss.str() );
	}
}
#endif
#endif

#ifndef NO_BLAS
std::vector<double> leastSquares_QR( fortran_matrix& A, const double* b_in )
{
	int m = int_cast( A.nrows() );
	int n = int_cast( A.ncols() );
	int nrhs = 1;
	int lda = int_cast( A.nrows() );
	std::vector<int> jpvt( n, 0 );
	double rcond = 1e-12;
	int lwork = -1;	// first query on size of work array
	std::vector<double> work(1);
	int info = 0;
	int rank = 0;

	// Again, plz use return value optimization :)
	std::vector<double> b( std::max<int>(m, n) );
	int ldb = int_cast( b.size() );
	std::copy( b_in, b_in + m, b.begin() );

	dgelsy(		&m,				// m
				&n,				// n
				&nrhs,			// nrhs
				A.data(),		// a
				&lda,			// lda
				&b[0],			// b
				&ldb,			// ldb
				&jpvt[0],		// jpvt
				&rcond,			// rcond
				&rank,			// rank
				&work[0],		// work
				&lwork,			// lwork
				&info );		// info

	// Now, know how much workspace is desired.
	lwork = static_cast<int>(work[0]);
	work.resize( lwork+1 );

	dgelsy(		&m,				// m
				&n,				// n
				&nrhs,			// nrhs
				A.data(),		// a
				&lda,			// lda
				&b[0],			// b
				&ldb,			// ldb
				&jpvt[0],		// jpvt
				&rcond,			// rcond
				&rank,			// rank
				&work[0],		// work
				&lwork,			// lwork
				&info );		// info

	if( info < 0 )
	{
		std::ostringstream oss;
		oss << "QR failed due to invalid parameter #" 
			<< -info;
		throw LinearAlgebraException( oss.str() );
	}
	else if( info > 0 )
	{
		std::ostringstream oss;
		oss << "QR did not complete; unknown error " << info;
		throw LinearAlgebraException( oss.str() );
	}

	b.resize( n );
	return b;
}


std::vector<double> leastSquares_SVD( fortran_matrix& A, const double* b_in )
{
	int m = int_cast( A.nrows() );
	int n = int_cast( A.ncols() );
	int nrhs = 1;
	double rcond = 1e-7;
	int lwork = -1;	// first query on size of work array
	int info = 0;
	int rank = 0;
	std::vector<double> s( std::max(m, n) );
	std::vector<double> work(1);
	int lda = int_cast( A.nrows() );

	std::vector<double> b( std::max(m, n) );
	int ldb = int_cast( b.size() );
	std::copy( b_in, b_in + m, b.begin() );

	dgelss(		&m,				// m
				&n,				// n
				&nrhs,			// nrhs
				A.data(),		// a
				&lda,			// lda
				&b[0],			// b
				&ldb,			// ldb
				&s[0],			// s
				&rcond,			// rcond
				&rank,			// rank
				&work[0],		// work
				&lwork,			// lwork
				&info );		// info

	// Now, know how much workspace is desired.
	lwork = static_cast<int>(work[0]);
	work.resize( lwork+1 );

	dgelss(		&m,				// m
				&n,				// n
				&nrhs,			// nrhs
				A.data(),		// a
				&lda,			// lda
				&b[0],			// b
				&ldb,			// ldb
				&s[0],			// s
				&rcond,			// rcond
				&rank,			// rank
				&work[0],		// work
				&lwork,			// lwork
				&info );		// info

	if( info < 0 )
	{
		std::ostringstream oss;
		oss << "SVD failed due to invalid parameter #" 
			<< -info;
		throw LinearAlgebraException( oss.str() );
	}
	else if( info > 0 )
	{
		std::ostringstream oss;
		oss << "SVD did not converge; " << info << 
			" diagonals of the intermediate bidiagonal form did not converge to 0.";
		throw LinearAlgebraException( oss.str() );
	}

	b.resize( n );
	return b;
}
#endif // NO_BLAS

// MATLAB .mat file format constants
const unsigned int miINT8             = 1;
const unsigned int miUINT8            = 2;
const unsigned int miINT16            = 3;
const unsigned int miUINT16           = 4;
const unsigned int miINT32            = 5;
const unsigned int miUINT32           = 6;
const unsigned int miSINGLE           = 7;
const unsigned int miDOUBLE           = 9;
const unsigned int miINT64            = 12;
const unsigned int miUINT64           = 13;
const unsigned int miMATRIX           = 14;

template <typename T>
struct MatlabTagForDataType
{
	unsigned int operator() () const
	{
		if( std::numeric_limits<T>::is_integer )
		{
			if( std::numeric_limits<T>::is_signed )
			{
				switch( sizeof(T) )
				{
				case 1:
					return miINT8;
				case 2:
					return miINT16;
				case 4:
					return miINT32;
				case 8:
					return miINT64;
				}
			}
			else
			{
				switch( sizeof(T) )
				{
				case 1:
					return miUINT8;
				case 2:
					return miUINT16;
				case 4:
					return miUINT32;
				case 8:
					return miUINT64;
				}
			}
		}
		else
		{
			switch( sizeof(T) )
			{
			case 4:
				return miSINGLE;
			case 8:
				return miDOUBLE;
			}
		}

		assert(false);
		return 0;
	}
};

const unsigned char mxCELL_CLASS      = 1;
const unsigned char mxSTRUCT_CLASS    = 2;
const unsigned char mxOBJECT_CLASS    = 3;
const unsigned char mxCHAR_CLASS      = 4;
const unsigned char mxSPARSE_CLASS    = 5;
const unsigned char mxDOUBLE_CLASS    = 6;
const unsigned char mxSINGLE_CLASS    = 7;
const unsigned char mxINT8_CLASS      = 8;
const unsigned char mxUINT8_CLASS     = 9;
const unsigned char mxINT16_CLASS     = 10;
const unsigned char mxUINT16_CLASS    = 11;
const unsigned char mxINT32_CLASS     = 12;
const unsigned char mxUINT32_CLASS    = 13;

template <typename T>
struct MatlabClassForDataType
{
	unsigned int operator() () const
	{
		if( std::numeric_limits<T>::is_integer )
		{
			if( std::numeric_limits<T>::is_signed )
			{
				switch( sizeof(T) )
				{
				case 1:
					return mxINT8_CLASS;
				case 2:
					return mxINT16_CLASS;
				case 4:
					return mxINT32_CLASS;
				case 8:
					assert(false);
				}
			}
			else
			{
				switch( sizeof(T) )
				{
				case 1:
					return mxUINT8_CLASS;
				case 2:
					return mxUINT16_CLASS;
				case 4:
					return mxUINT32_CLASS;
				case 8:
					assert(false);
				}
			}
		}
		else
		{
			switch( sizeof(T) )
			{
			case 4:
				return mxSINGLE_CLASS;
			case 8:
				return mxDOUBLE_CLASS;
			}
		}

		assert(false);
		return 0;
	}
};


class MATFileDataElement
{
public:
	virtual unsigned int byteCount() const = 0;
	virtual void dump( std::ostream& out ) const = 0;
	virtual ~MATFileDataElement() {}
};

template <typename T>
class SimpleMATFileDataElement
	: public MATFileDataElement
{
public:
	SimpleMATFileDataElement( const T* data, size_t length )
		: data_(data), length_(length)
	{
		size_t bytes = length_ * sizeof(T);
		padding_.resize( (8 - (bytes % 8)) % 8, 0 );
	}

	size_t byteCount() const
	{
		return 8 + dataByteCount() + padding_.size();
	}

	size_t dataByteCount() const
	{
		return length_ * sizeof(T);
	}

	void dump(std::ostream& out) const
	{
		// First, need to dump the tag.
		unsigned int tag = MatlabTagForDataType<T>()();
		out.write( reinterpret_cast<const char*>( &tag ), sizeof(tag) );

		size_t length = dataByteCount();
		out.write( reinterpret_cast<const char*>( &length ), sizeof(length) );

		if( data_ )
            out.write( reinterpret_cast<const char*>( data_ ), 
                boost::numeric_cast<std::streamsize>(dataByteCount()) );

		if( !padding_.empty() )
			out.write( &padding_[0], 
                boost::numeric_cast<std::streamsize>(sizeof(char) * padding_.size()) );
	}

private:
	const T* data_;
	size_t length_;
	std::vector<char> padding_;
};

class CompoundMATFileDataElement
{
public:
	CompoundMATFileDataElement( unsigned int type )
		: type_(type), byteCount_(0)
	{
	}

	void add( boost::shared_ptr<const MATFileDataElement> child )
	{
		children_.push_back( child );
		byteCount_ += child->byteCount();
	}

	unsigned int byteCount() const
	{
		return 8 + byteCount_;
	}

	void dump(std::ostream& out) const
	{
		out.write( reinterpret_cast<const char*>( &type_ ), sizeof(type_) );
		out.write( reinterpret_cast<const char*>( &byteCount_ ), sizeof(byteCount_) );

		for( ChildList::const_iterator itr = children_.begin();
			itr != children_.end();
			++itr )
		{
			(*itr)->dump( out );
		}
	}

private:
	typedef std::deque< boost::shared_ptr<const MATFileDataElement> > ChildList;
	ChildList children_;
	unsigned int type_;
	unsigned int byteCount_;
};

MATFile::MATFile( const std::string& filename, const std::string& description )
	: ofs( filename.c_str(), std::ios::out | std::ios::binary )
{
    if( !ofs )
        throw IOException( "Unable to open file '" + filename + "' for writing." );

    ofs.exceptions( std::ifstream::eofbit | std::ifstream::failbit | std::ifstream::badbit );

	const size_t descriptionLength = 116;
	struct headerFormat
	{
		char description[descriptionLength];
		char subsystemSpecific[8];
		short version;
		char endian[2];
	} header;

	std::fill( header.description, header.description+descriptionLength, 0 );
	std::copy( &description[0], &description[ std::min(description.size(), descriptionLength) ], header.description );
	std::fill( header.subsystemSpecific, header.subsystemSpecific+8, 0 );
	header.version = 0x100;
	header.endian[0] = 'I';
	header.endian[1] = 'M';
	ofs.write( (char*) &header, sizeof( headerFormat ) );
    ofs.flush();
}

template <typename T>
void addMat( std::ofstream& ofs, const char* name,
	const boost::numeric::ublas::matrix< T, boost::numeric::ublas::column_major >& matrix )
{
	CompoundMATFileDataElement matElement(miMATRIX);

	std::vector<unsigned char> arrayFlags(8, 0);
	arrayFlags[0] = MatlabClassForDataType<T>()();
	boost::shared_ptr<MATFileDataElement > arrayFlagsElt(
		new SimpleMATFileDataElement<unsigned int>( 
			reinterpret_cast<const unsigned int *>(&arrayFlags[0]), 2 ) );
	matElement.add( arrayFlagsElt );

	// Dimensions array
	std::vector<int> dimensions( 2 );
	dimensions[0] = int_cast( matrix.size1() );
	dimensions[1] = int_cast( matrix.size2() );
	boost::shared_ptr<MATFileDataElement> dimensionsElt(
		new SimpleMATFileDataElement<int>( &dimensions[0], 2 ) );
	matElement.add( dimensionsElt );

	// Array name element
	std::vector<signed char> nameSigned;
    while( *name )
        nameSigned.push_back( *name++ );
	boost::shared_ptr<MATFileDataElement> nameElt(
		new SimpleMATFileDataElement<signed char>( &nameSigned[0], nameSigned.size() ) );
	matElement.add( nameElt );

	// Real part element
	boost::shared_ptr<MATFileDataElement> matElt;
	if( matrix.size1() == 0 || matrix.size2() == 0 )
	{
		matElt.reset(
			new SimpleMATFileDataElement<T>( 0, 0 ) );
	}
	else
	{
		matElt.reset(
			new SimpleMATFileDataElement<T>( &(matrix.data())[0], matrix.size1() * matrix.size2() ) );
	}
	matElement.add( matElt );

	matElement.dump( ofs );
    ofs.flush();
}

template <typename T>
void addMat( std::ofstream& ofs, const char* name,
			const T* a, size_t cols, size_t lda, size_t rows )
{
	// would be nice to do this faster, but it looks like I'm going to have to pack it into a 
	//   dense matrix.
	boost::numeric::ublas::matrix< T, boost::numeric::ublas::column_major > mat( rows, cols );

	for( size_t iRow = 0; iRow < rows; ++iRow )
		for( size_t iCol = 0; iCol < cols; ++iCol )
			mat( iRow, iCol ) = a[ iRow*lda + iCol ];

	addMat( ofs, name, mat );
}

template <typename T>
void addMat( std::ofstream& ofs, const char* name,
			const std::vector<T>& a, size_t cols, size_t lda )
{
	if( lda == 0 )
	{
		assert( a.size() == 0 );
		addMat( ofs, name, &a[0], 0, 0, 0 );
		return;
    }

	size_t rows = a.size() / lda;
	assert( (a.size() % lda) == 0 );

	addMat( ofs, name, &a[0], cols, lda, rows );
}

void MATFile::add( const char* name, const std::vector<double>& data, size_t cols, size_t lda )
{
	addMat( this->ofs, name, data, cols, lda );
}

void MATFile::add( const char* name, const std::vector<float>& data, size_t cols, size_t lda )
{
	addMat( this->ofs, name, data, cols, lda );
}

void MATFile::add( const char* name, const double* data, size_t cols, size_t lda, size_t rows )
{
	addMat( this->ofs, name, data, cols, lda, rows );
}

void MATFile::add( const char* name, const float* data, size_t cols, size_t lda, size_t rows )
{
	addMat( this->ofs, name, data, cols, lda, rows );
}

void MATFile::add( const char* name, const fortran_matrix& matrix )
{
	addMat( this->ofs, name, matrix.array() );
}

void MATFile::add( const char* name, const int_matrix& matrix )
{
	addMat( this->ofs, name, matrix.array() );
}

void MATFile::add( const char* name, const uint_matrix& matrix )
{
	addMat( this->ofs, name, matrix.array() );
}

void MATFile::add( const char* name, const generic_matrix<float, boost::numeric::ublas::column_major>& matrix )
{
	addMat( this->ofs, name, matrix.array() );
}

#ifdef NAG
void MATFile::add( const std::string& name, const CoordinateStorage& matrix )
{
	CompoundMATFileDataElement matElement(miMATRIX);

	// Array flags
	std::vector<unsigned char> arrayFlags(8, 0);
	arrayFlags[0] = mxSPARSE_CLASS;
	boost::shared_ptr<MATFileDataElement> arrayFlagsElt(
		new SimpleMATFileDataElement<unsigned int>( 
			reinterpret_cast<const unsigned int *>(&arrayFlags[0]), 2 ) );
	unsigned int* nz = reinterpret_cast<unsigned int*>(&arrayFlags[4]);
	*nz = matrix.nnz();
	matElement.add( arrayFlagsElt );

	// Dimensions array
	std::vector<int> dimensions( 2 );
	dimensions[0] = matrix.n();
	dimensions[1] = matrix.n();
	boost::shared_ptr<MATFileDataElement> dimensionsElt(
		new SimpleMATFileDataElement<int>( &dimensions[0], 2 ) );
	matElement.add( dimensionsElt );

	// Array name element
	std::vector<signed char> nameSigned( name.size() );
	std::copy( name.begin(), name.end(), nameSigned.begin() );
	boost::shared_ptr<MATFileDataElement> nameElt(
		new SimpleMATFileDataElement<signed char>( &nameSigned[0], nameSigned.size() ) );
	matElement.add( nameElt );

	// First, we need to figure out where all the rows are
	std::vector<int> rowStarts; rowStarts.reserve( matrix.n() );
	for( int i = 0; i < matrix.nnz(); ++i )
	{
		for( int j = rowStarts.size(); j < matrix.irow(i); ++j )
			rowStarts.push_back(i);
	}

	// Now, we need to go through in column-major order.
	std::vector<double> pr;		pr.reserve( matrix.nnz() );
	std::vector<int> ir;		ir.reserve( matrix.nnz() );
	std::vector<int> jc;		jc.reserve( matrix.n() + 1 );
	for( int iCol = 0; iCol < matrix.n(); ++iCol )
	{
		jc.push_back( ir.size() );
		for( int iRow = 0; iRow < rowStarts.size(); ++iRow )
		{
			while( (matrix.icol( rowStarts[iRow] )-1) < iCol && (matrix.irow( rowStarts[iRow] )-1) == iRow )
				++rowStarts[iRow];
			if( (matrix.icol( rowStarts[iRow] )-1) == iCol &&
				(matrix.irow( rowStarts[iRow] )-1) == iRow )
			{
				ir.push_back( iRow );
				pr.push_back( *(matrix.data() + rowStarts[iRow]) );
			}
		}
	}
	jc.push_back( ir.size() );

	// Now, we have all these matrices, so in theory we just dump
	// them to disk
	boost::shared_ptr<MATFileDataElement> irElt(
		new SimpleMATFileDataElement<int>( &ir[0], ir.size() ) );
	matElement.add( irElt );

	boost::shared_ptr<MATFileDataElement> jcElt(
		new SimpleMATFileDataElement<int>( &jc[0], jc.size() ) );
	matElement.add( jcElt );

	boost::shared_ptr<MATFileDataElement> prElt(
		new SimpleMATFileDataElement<double>( &pr[0], pr.size() ) );
	matElement.add( prElt );

	matElement.dump( ofs );
}

void MATFile::add( const std::string& name, const SymmetricCoordinateStorage& matrix )
{
	assert( false );
}
#endif

void MATFile::add( const char* name, const std::vector<unsigned int>& vector )
{
	add( name, &vector[0], vector.size() );
}

void MATFile::add( const char* name, const std::vector<int>& vector )
{
	add( name, &vector[0], vector.size() );
}

void MATFile::add( const char* name, const std::vector<double>& vector )
{
	add( name, &vector[0], vector.size() );
}

void MATFile::add( const char* name, const std::vector<float>& vector )
{
	add( name, &vector[0], vector.size() );
}

void MATFile::add( const char* name, const double* vector, size_t length )
{
	fortran_matrix mat( length, 1 );
	std::copy( vector, vector+length, mat.data() );
	add( name, mat );
}

void MATFile::add( const char* name, const float* vector, size_t length )
{
	generic_matrix< float, boost::numeric::ublas::column_major > mat( length, 1 );
	std::copy( vector, vector+length, mat.data() );
	addMat( this->ofs, name, mat.array() );
}

void MATFile::add( const char* name, const unsigned int* vector, size_t length )
{
	uint_matrix mat( length, 1 );
	std::copy( vector, vector+length, mat.data() );
	add( name, mat );
}

void MATFile::add( const char* name, const int* vector, size_t length )
{
	int_matrix mat( length, 1 );
	std::copy( vector, vector+length, mat.data() );
	add( name, mat );
}

void MATFile::add( const char* name, double scalar )
{
	fortran_matrix mat( 1, 1 );
	mat(0, 0) = scalar;
	add( name, mat );
}

#ifdef HDF5
herr_t getErrorMsg( int n, H5E_error_t *err_desc, void *client_data )
{
	std::deque<std::string>* tempVec = reinterpret_cast< std::deque<std::string>* >( client_data );
	std::deque<std::string>& vec = *tempVec;
	H5E_major_t major = err_desc->maj_num;
	H5E_minor_t minor = err_desc->min_num;
	vec.push_back( std::string(H5Eget_major(major)) + ": " + std::string(H5Eget_minor(minor)) );
	if( err_desc->desc )
		vec.back().append( std::string(": ") + std::string(err_desc->desc) );

	return 0;
}

HDF5File::HDF5File( const std::string& filename )
{
	file_ = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
	if( file_ < 0 )
	{
		std::deque<std::string> messages;
		H5Ewalk( H5E_WALK_UPWARD, getErrorMsg, &messages );
		std::ostringstream oss;
		oss << "Unable to open file '" + filename + "' for writing:\n";
		for( std::deque<std::string>::const_iterator iter = messages.begin();
			iter != messages.end(); ++iter )
		{
			oss << *iter << "\n";
		}

		throw IOException( oss.str() );
	}
}

HDF5File::~HDF5File()
{
	H5Fclose(file_);
}

void HDF5File::add( const std::string& name, const fortran_matrix& matrix )
{
	// H5F insists on row-major:
	c_matrix columnMajor( matrix.nrows(), matrix.ncols() );
	for( size_t i = 0; i < matrix.nrows(); ++i )
		for( size_t j = 0; j < matrix.ncols(); ++j )
			columnMajor(i, j) = matrix(i, j);

	hsize_t dims[] = { matrix.nrows(), matrix.ncols() };
	hsize_t maxDims[] = { matrix.nrows(), matrix.ncols() };
	hid_t dataspace = H5Screate_simple( 
		2,
		dims,
		maxDims );

	hid_t datatype = H5Tcopy( H5T_NATIVE_DOUBLE );
	herr_t status = H5Tset_order(datatype, H5T_ORDER_LE);

	hid_t dataset = H5Dcreate(file_, name.c_str(), datatype, dataspace, H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, columnMajor.data());

	H5Sclose(dataspace);
	H5Tclose(datatype);
	H5Dclose(dataset);
}
#endif // HDF5

template <typename T, typename IntType>
T norm2_imp( const T* x, IntType n )
{
	T result(0);
	for( IntType i = 0; i < n; ++i )
	{
		T val = x[i];
		result += val*val;
	}

	return result;
}

template <typename T, typename IntType>
T norminf_imp( const T* x, IntType n )
{
	T result(0);
	for( IntType i = 0; i < n; ++i )
		result = std::max( result, std::abs(x[i]) );
	return result;
}

template <typename T, typename IntType>
T absmin_imp( const T* x, IntType n )
{
	assert( n > 0 );
	T result( std::abs(x[0]) );
	for( IntType i = 1; i < n; ++i )
		result = std::min( std::abs(x[i]), result );
	return result;
}

float norm2( const float* x, int n )
{
#ifdef NO_BLAS
	return norm2_imp( x, n );
#else
	return cblas_snrm2( n, x, 1 );
#endif
}

float norminf( const float* x, int n )
{
#ifdef NO_BLAS
	return norminf_imp( x, n );
#else
	size_t index = cblas_isamax( n, x, 1 );
	return std::abs(x[index]);
#endif
}

float absmin( const float* x, int n )
{
#ifdef NO_BLAS
	return absmin_imp( x, n );
#else
	size_t index = cblas_isamin( int_cast( n ), x, 1 );
	return std::abs(x[index]);
#endif
}

float norm2( const std::vector<float>& x )
{
#ifdef NO_BLAS
	return norm2_imp( &x[0], x.size() );
#else
	return cblas_snrm2( int_cast( x.size() ), &x[0], 1 );
#endif
}

float norminf( const std::vector<float>& x )
{
#ifdef NO_BLAS
	return norminf_imp( &x[0], x.size() );
#else
	size_t index = cblas_isamax( int_cast( x.size() ), &x[0], 1 );
	return std::abs(x[index]);
#endif
}

float absmin( const std::vector<float>& x )
{
#ifdef NO_BLAS
	return absmin_imp( &x[0], x.size() );
#else
	size_t index = cblas_isamin( int_cast( x.size() ), &x[0], 1 );
	return std::abs(x[index]);
#endif
}

double norm2( const fortran_matrix& x )
{
#ifdef NO_BLAS
	return norm2_imp( x.data(), x.nrows()*x.ncols() );
#else
	return cblas_dnrm2( int_cast(x.nrows()*x.ncols()), x.data(), 1 );
#endif
}

double norminf( const fortran_matrix& x )
{
#ifdef NO_BLAS
	return norminf_imp( x.data(), x.nrows()*x.ncols() );
#else
	size_t index = cblas_idamax( int_cast(x.nrows()*x.ncols()), x.data(), 1 );
	return std::abs( *(x.data() + index) );
#endif
}

double absmin( const fortran_matrix& x )
{
#ifdef NO_BLAS
	return absmin_imp( x.data(), x.nrows()*x.ncols() );
#else
	size_t index = cblas_idamin( int_cast(x.nrows()*x.ncols()), x.data(), 1 );
	return fabs( *(x.data() + index) );
#endif
}

double norm2( const std::vector<double>& x )
{
#ifdef NO_BLAS
	return norm2_imp( &x[0], x.size() );
#else
	return cblas_dnrm2( int_cast(x.size()), &x[0], 1 );
#endif
}

double norminf( const std::vector<double>& x )
{
#ifdef NO_BLAS
	return norminf_imp( &x[0], x.size() );
#else
	size_t index = cblas_idamax( int_cast(x.size()), &x[0], 1 );
	return fabs(x[index]);
#endif
}

double absmin( const std::vector<double>& x )
{
#ifdef NO_BLAS
	return absmin_imp( &x[0], x.size() );
#else
	size_t index = cblas_idamin( int_cast(x.size()), &x[0], 1 );
	return fabs(x[index]);
#endif
}

float absmax( const vl::Mat3f& mat )
{
	float result = 0.0f;
    for(vl::Int i = 0; i < 3; ++i )
		for( vl::Int j = 0; j < 3; ++j )
			result = std::max( fabs(mat[i][j]), result );

	return result;
}

float absmax( const vl::Mat3d& mat )
{
	double result = 0.0f;
	for(vl::Int i = 0; i < 3; ++i )
		for( vl::Int j = 0; j < 3; ++j )
			result = std::max( fabs(mat[i][j]), result );

	return result;
}

template <typename MatType>
std::pair<MatType, MatType> polarDecompositionDriver( const MatType& matrix )
{
	const float epsilon = 1e-4;
	MatType Q = matrix;
	while( true )
	{
		MatType Q_next = 0.5*(Q + trans(inv(Q)));
		if( absmax(Q - Q_next) < epsilon )
			break;

		Q = Q_next;
	}

	MatType S = inv(Q) * matrix;
	return std::make_pair( Q, S );
}

std::pair<vl::Mat3f, vl::Mat3f> polarDecomposition( const vl::Mat3f& matrix )
{
	return polarDecompositionDriver( matrix );
}

std::pair<vl::Mat3d, vl::Mat3d> polarDecomposition( const vl::Mat3d& matrix )
{
	return polarDecompositionDriver( matrix );
}

/*
 * A function for creating a rotation matrix that rotates a vector called
 * "from" into another vector called "to".
 * Input : from[3], to[3] which both must be *normalized* non-zero vectors
 * Output: mtx[3][3] -- a 3x3 matrix in colum-major form
 * Authors: Tomas Möller, John Hughes
 *          "Efficiently Building a Matrix to Rotate One Vector to Another"
 *          Journal of Graphics Tools, 4(4):1-4, 1999
 */
vl::Mat3f fromToRotation(const vl::Vec3f& from, const vl::Vec3f& to)
{
	vl::Mat3f mtx;

	const float EPSILON = 0.000001;

	vl::Vec3f v = vl::cross(from, to);
	float e = vl::dot(from, to);
	float f = (e < 0)? -e:e;
	if (f > 1.0 - EPSILON)       // "from" and "to"-vector almost parallel
	{
		vl::Vec3f x;             // vector most nearly orthogonal to "from"
		x[0] = (from[0] > 0.0)? from[0] : -from[0];
		x[1] = (from[1] > 0.0)? from[1] : -from[1];
		x[2] = (from[2] > 0.0)? from[2] : -from[2];

		if (x[0] < x[1])
		{
			if (x[0] < x[2])
			{
				x[0] = 1.0; x[1] = x[2] = 0.0;
			}
			else
			{
				x[2] = 1.0; x[0] = x[1] = 0.0;
			}
		}
		else
		{
			if (x[1] < x[2])
			{
				x[1] = 1.0; x[0] = x[2] = 0.0;
			}
			else
			{
				x[2] = 1.0; x[0] = x[1] = 0.0;
			}
		}

		vl::Vec3f u, v;          // temporary storage vectors
		u[0] = x[0] - from[0]; u[1] = x[1] - from[1]; u[2] = x[2] - from[2];
		v[0] = x[0] - to[0];   v[1] = x[1] - to[1];   v[2] = x[2] - to[2];

		float c1 = 2.0 / vl::dot(u, u);
		float c2 = 2.0 / vl::dot(v, v);
		float c3 = c1 * c2  * vl::dot(u, v);

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				mtx[i][j] =  - c1 * u[i] * u[j]
								- c2 * v[i] * v[j]
								+ c3 * v[i] * u[j];
			}
			mtx[i][i] += 1.0;
		}
	}
	else  /* the most common case, unless "from"="to", or "from"=-"to" */
	{
		/* unoptimized version - a good compiler will optimize this. */
		/* h = (1.0 - e)/DOT(v, v); old code */
		float h = 1.0/(1.0 + e);      /* optimization by Gottfried Chen */
		mtx[0][0] = e + h * v[0] * v[0];
		mtx[0][1] = h * v[0] * v[1] - v[2];
		mtx[0][2] = h * v[0] * v[2] + v[1];

		mtx[1][0] = h * v[0] * v[1] + v[2];
		mtx[1][1] = e + h * v[1] * v[1];
		mtx[1][2] = h * v[1] * v[2] - v[0];

		mtx[2][0] = h * v[0] * v[2] - v[1];
		mtx[2][1] = h * v[1] * v[2] + v[0];
		mtx[2][2] = e + h * v[2] * v[2];
	}

	return mtx;
}

// helper functions for below
double maximumOffDiagonal( const fortran_matrix& A )
{
	double maxValue = 0.0;
	for( size_t i = 0; i < A.nrows(); ++i )
		for( size_t j = 0; j < A.ncols(); ++j )
		{
			if( i == j )
				continue;

			maxValue = std::max( maxValue, std::abs( A(i,j) ) );
		}

	return maxValue;
}

// returns c, s
std::pair< double, double > symmetricShur2( const fortran_matrix& A, size_t p, size_t q )
{
	const double epsilon = 10.0*std::numeric_limits<double>::epsilon();

	if( std::abs( A(p, q) ) > epsilon )
	{
		double tau = (A(q, q) - A(p, p)) / (2.0*A(p,q));

		double t;
		if( tau >= 0.0 )
			t = 1.0 / (tau + sqrt(1 + tau*tau));
		else
			t = -1.0 / (-tau + sqrt(1 + tau*tau));

		double c = 1.0 / sqrt(1 + t*t);
		double s = t*c;
		return std::make_pair( c, s );
	}
	else
	{
		double c = 1.0;
		double s = 0.0;
		return std::make_pair( c, s );
	}
}

double frobeniusNorm( const fortran_matrix& A )
{
	double result = 0.0;
	for( size_t i = 0; i < A.nrows(); ++i )
	{
		for( size_t j = 0; j < A.ncols(); ++j )
		{
			const double a = A(i,j);
			result += a*a;
		}
	}

	return sqrt(result);
}

// Compute the eigendecomposition using the cylic-by-row Jacobi algorithm from p. 430 of 
//   Golub & Van Loan
std::pair< fortran_matrix, std::vector<double> > symmetricJacobi( fortran_matrix A )
{
	assert( A.nrows() == A.ncols() );
	const size_t n = A.nrows();

	fortran_matrix V( n, n, 0.0 );
	for( size_t i = 0; i < n; ++i )
		V(i, i) = 1.0;

	const double epsilon = 1e-10 * frobeniusNorm(A);
	while( maximumOffDiagonal(A) > epsilon )
	{
		for( size_t p = 0; p < n-1; ++p )
		{
			for( size_t q = p+1; q < n; ++q )
			{
				double c, s;
				boost::tie(c, s) = symmetricShur2( A, p, q );
				assert( std::abs(c*c + s*s - 1.0) < 1e-10 );

				// todo: this would be much faster if we didn't construct
				//   J explicitly, but this first version will be useful
				//   for testing later
				fortran_matrix J( n, n, 0.0 );
				for( size_t i = 0; i < n; ++i )
					J(i, i) = 1.0;

				J( p, p ) = c;
				J( q, q ) = c;
				J( p, q ) = s;
				J( q, p ) = -s;

				fortran_matrix tmp = prod(A, J);
				V = prod(V, J);

				// transpose J:
				std::swap( J(q, p), J(p, q) );
				A = prod( J, tmp );
			}
		}
	}

	std::vector<double> eigs( n );
	for( size_t i = 0; i < n; ++i )
		eigs[i] = A(i,i);

	return std::make_pair( V, eigs );
}