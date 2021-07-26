#ifndef __LINALG_H__
#define __LINALG_H__

// wrappers for useful linear algebra stuff, like SVD

#include "exception.h"
#include "util.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp> // for io

#include <iostream>
#include <fstream>

// This to make the boost matrix (or whatever matrix type we use)
// look like an MTL matrix for code backwards compatibility.
template <typename T, typename Order>
class generic_matrix
{
public:
	typedef boost::numeric::ublas::matrix< T, Order > matrix_type;
	typedef boost::numeric::ublas::matrix_range< matrix_type > matrix_range;
	typedef generic_matrix<T, Order> submatrix_type;
	typedef T value_type;

	explicit generic_matrix( const matrix_type& mat )
		: array_(mat) {}
	explicit generic_matrix(size_t rows, size_t cols)
		: array_( rows, cols ) {}

	explicit generic_matrix(size_t rows, size_t cols, T value)
		: array_( rows, cols )
	{
		std::fill( data(), data() + nrows()*ncols(), value );
	}

	explicit generic_matrix()
		: array_() {}

	const T& operator()(size_t i, size_t j) const		{ return array_(i, j); }
	T& operator()(size_t i, size_t j)					{ return array_(i, j); }
	const T* data() const
	{
		if(array_.data().size() == 0)
			return 0;
		return &(array_.data())[0];
	}

	T* data()
	{
		if(array_.data().size() == 0)
			return 0;
		return &(array_.data())[0];
	}

	size_t nrows() const												{ return array_.size1(); }
	size_t ncols() const												{ return array_.size2(); }
	void resize( size_t rows, size_t cols )								{ array_.resize( rows, cols, true ); }

	submatrix_type submatrix( unsigned int rowStart, unsigned int rowEnd, unsigned int colStart, unsigned int colEnd )
	{
		return submatrix_type( 
			boost::numeric::ublas::matrix_range<matrix_type>( 
				array_, 
				boost::numeric::ublas::range(rowStart, rowEnd),
				boost::numeric::ublas::range(colStart, colEnd) ) );
	}

	const matrix_type& array() const
	{
		return array_;
	}

private:
	matrix_type array_;
};

template <typename T, typename Order>
inline bool isBad(const generic_matrix<T, Order>& mat)
{
	for( unsigned int i = 0; i < mat.nrows(); ++i )
		for( unsigned int j = 0; j < mat.ncols(); ++j )
			if(isBad(mat(i,j)))
				return true;

	return false;
}

typedef generic_matrix< double, boost::numeric::ublas::column_major > fortran_matrix;
typedef generic_matrix< int, boost::numeric::ublas::column_major > int_matrix;
typedef generic_matrix< unsigned int, boost::numeric::ublas::column_major > uint_matrix;
typedef generic_matrix< double, boost::numeric::ublas::row_major > c_matrix;

template <typename matrix_type>
class Transposed
{
public:
	explicit Transposed( const matrix_type& matrix )
		: matrix_(matrix) {}

	const matrix_type& mat() const	{ return matrix_; }

private:
	const matrix_type& matrix_;
};

template <typename vec_type>
class Scaled
{
public:
	explicit Scaled( const vec_type& vec, double scale = 1.0 )
		: vec_( vec ), scale_(scale) {}

	double scale() const { return scale_; }
	const vec_type& vec() const { return vec_; }

private:
	const vec_type& vec_;
	const double scale_;
};

template <typename vec_type>
Scaled<vec_type> scaled( const vec_type& x, const double scale )
{
	return Scaled<vec_type>(x, scale);
}

template <typename mat_type>
Transposed<mat_type> transposed( const mat_type& mat )
{
	return Transposed<mat_type>( mat );
}


class LinearAlgebraException
	: public Exception
{
public:
	explicit LinearAlgebraException( const std::string& message )
		: Exception(message) {}

	virtual ~LinearAlgebraException() {}
};

#ifndef NO_BLAS
class SVD
{
public:
	explicit SVD( const fortran_matrix& matrix );
	const fortran_matrix& U() const;
	const fortran_matrix& Vt() const;
	const std::vector<double>& S() const;

private:
	fortran_matrix U_;
	fortran_matrix Vt_;
	std::vector<double> S_;
};

// Runs SVD on a (smallish) random matrix.
void test_svd(std::ostream& out);
#endif // NO_BLAS

#ifndef NO_BLAS
#ifndef _WIN32
class NonNegativeMatrixFactorization
{
public:
	explicit NonNegativeMatrixFactorization( const fortran_matrix& matrix, const unsigned int nBasis, const unsigned int maxIter = 0 );
	const fortran_matrix& W();
	const fortran_matrix& H();

private:
	fortran_matrix W_;
	fortran_matrix H_;
};
#endif
#endif

template <typename Mat>
void printRow(std::ostream& os, const int rowNum, const int numberLength, const Mat& mat)
{
	const static int entrySeparation = 1;
	for (unsigned int i = 0; i < mat.ncols(); ++i) {
//		std::stringstream sstream;
//		sstream << mat(rowNum, i) << std::flush;
//		std::string s(sstream.str());
//		while (s.length() < 10 + entrySeparation)
//			s += " ";
//		os << s;
		os << mat(rowNum, i) << std::string(entrySeparation, ' ');
	}
}

template <typename Mat>
void printDense(std::ostream& os, const Mat& matrix, std::string name)
{
	os.precision(20);

	unsigned int length = 10;
	os << name << " = [";
	for (unsigned int i = 0; i < matrix.nrows(); ++i) {
		if (i != 0)
			os << std::string( name.size() + 4, ' ' ); 

		printRow(os, i, length, matrix);

		if( i == matrix.nrows() - 1 )
			os << "];" << std::endl;
		else
			os << ";" << std::endl;
	}
}

template <typename T>
void printVector( std::ostream& os, const T& vec, size_t size, const std::string& name )
{
	os.precision(20);

	os << name << " = [ ";
	for( int i = 0; i < size; ++i )
	{
		os << vec[i] << " ";
	}
	os << "]; " << std::endl;
}

template <typename Matrix>
void resetMatrix( Matrix& A )
{
	typename Matrix::iterator i;
	typename Matrix::OneD::iterator j;
	for (i = A.begin(); i != A.end(); ++i)
		for (j = (*i).begin(); j != (*i).end(); ++j)
			*j = 0.0;
}

template <typename Matrix>
double symmetricRatio( Matrix& A )
{
	double numerator = 0.0;
	double denominator = 0.0;
	typename Matrix::iterator i;
	typename Matrix::OneD::iterator j;
	for (i = A.begin(); i != A.end(); ++i)
	{
		for (j = (*i).begin(); j != (*i).end(); ++j)
		{
			double symmetricVal = A( j.column(), j.row() );
			double plus = 0.5 * ((*j) + symmetricVal);
			double minus = 0.5 * ((*j) - symmetricVal);
			numerator += minus*minus;
			denominator += plus*plus;
		}
	}

	return sqrt(numerator)/sqrt(denominator);
}

// y = a*y
// Uses the CBLAS dscal operator
void scale( fortran_matrix& mat, double scale );
void scale( std::vector<double>& vec, double scale );
void scale( double* vec, size_t n, double scale );

// y = y + a*x
// equivalent to a CBLAS daxpy operator
void add( const Scaled<fortran_matrix>& x, fortran_matrix& y );
void add( const Scaled< std::vector<double> >& x, std::vector<double>& y );
void add( const Scaled< std::vector<float> >& x, std::vector<float>& y );
void add( const fortran_matrix& x, fortran_matrix& y );
void add( const std::vector<double>& x, std::vector<double>& y );
void add( const std::vector<float>& x, std::vector<float>& y );

// Use a fast matrix-matrix product here:
fortran_matrix prod( const fortran_matrix& left, const fortran_matrix& right );
fortran_matrix prod( const fortran_matrix& left, const Transposed<fortran_matrix>& right );
fortran_matrix prod( const Transposed<fortran_matrix>& left, const fortran_matrix& right );

// Simple matrix-vector product.  out need not be initialized.
void prod( const fortran_matrix& left, const double* rhs, double* out );
void prod( const Transposed<fortran_matrix>& left, const double* rhs, double* out );

std::vector<double> prod( const fortran_matrix& left, const double* rhs );
std::vector<double> prod( const Transposed<fortran_matrix>& left, const double* rhs );

void prodAdd( const fortran_matrix& left, const Scaled< std::vector<double> >& x, std::vector<double>& y );

// Solve the equation Ax = b in least squares (minimize (Ax-b)'(Ax-b))
// x is returned.
// note that the matrix A will be destroyed as a result.
#ifndef NO_BLAS
std::vector<double> leastSquares_QR( fortran_matrix& A, const double* b );
std::vector<double> leastSquares_SVD( fortran_matrix& A, const double* b );
#ifndef _WIN32
std::vector<double> leastSquares_NNLS( fortran_matrix& A, const double* b );
#endif // _WIN32
#endif // NO_BLAS

#ifndef CONDOR
float norm2( const float* x, int n );
float norminf( const float* x, int n );
float absmin( const float* x, int n );

float norm2( const std::vector<float>& x );
float norminf( const std::vector<float>& x );
float absmin( const std::vector<float>& x );
#endif

double norm2( const fortran_matrix& x );
double norminf( const fortran_matrix& x );
double absmin( const fortran_matrix& x );

double norm2( const std::vector<double>& x );
double norminf( const std::vector<double>& x );
double absmin( const std::vector<double>& x );

#ifdef NAG
class CoordinateStorage;
class SymmetricCoordinateStorage;
#endif

// Wraps a .mat file, so we can output matrices to Matlab.
// Uses this document:
// http://www.mathworks.com/access/helpdesk/help/pdf_doc/matlab/matfile_format.pdf
class MATFile
{
public:
	explicit MATFile( const std::string& filename, const std::string& description );

	void add( const char* name, const generic_matrix<float, boost::numeric::ublas::column_major>& matrix );
	void add( const char* name, const fortran_matrix& matrix );
	void add( const char* name, const int_matrix& matrix );
	void add( const char* name, const uint_matrix& matrix );

	void add( const char* name, const std::vector<double>& vector );
	void add( const char* name, const std::vector<float>& vector );
	void add( const char* name, const std::vector<int>& vector );
	void add( const char* name, const std::vector<unsigned int>& vector );

	void add( const char* name, const std::vector<double>& data, size_t cols, size_t lda );
	void add( const char* name, const std::vector<float>& data, size_t cols, size_t lda );
	void add( const char* name, const double* data, size_t cols, size_t lda, size_t rows );
	void add( const char* name, const float* data, size_t cols, size_t lda, size_t rows );

	void add( const char* name, const double* vector, size_t length );
	void add( const char* name, const float* vector, size_t length );
	void add( const char* name, const int* vector, size_t length );
	void add( const char* name, const unsigned int* vector, size_t length );


	void add( const char* name, double scalar );

#ifdef NAG
	void add( const char* name, const CoordinateStorage& matrix );
	void add( const char* name, const SymmetricCoordinateStorage& matrix );
#endif

private:
	std::ofstream ofs;
};

#ifdef HDF5
// NCSA H5F format
// Useful for getting data into the R statistical package.
class HDF5File
{
public:
	explicit HDF5File( const std::string& filename );
	~HDF5File();

	void add( const std::string& name, const fortran_matrix& matrix );

private:
	int file_;
};
#endif


template <typename Container>
void dumpArray( const Container& arr, std::ostream& os )
{
	typedef typename Container::value_type T;
	int size = arr.size();
	os.write( reinterpret_cast<const char*>(&size), sizeof(int) );

	if( size > 0 )
		os.write( reinterpret_cast<const char*>(&arr[0]), size*sizeof(T) );
}

template <typename Container>
void dumpArray( const Container& array, std::ostream& os, int size )
{
    typedef typename Container::value_type T;
	if( size > 0 )
	    os.write( reinterpret_cast<const char*>(&array[0]), size*sizeof(T) );
}

template <typename T, typename Order>
void dumpMatrix( const generic_matrix<T, Order>& mat, std::ostream& os )
{
	int nrows = mat.nrows();
	int ncols = mat.ncols();
	os.write( reinterpret_cast<const char*>(&nrows), sizeof(int) );
	os.write( reinterpret_cast<const char*>(&ncols), sizeof(int) );

	if( nrows > 0 && ncols > 0 )
		os.write( reinterpret_cast<const char*>(mat.data()), nrows*ncols*sizeof(T) );
}

template <typename T, typename Order>
void dumpMatrix( const generic_matrix<T, Order>& mat, std::ostream& os, int nrows, int ncols )
{
	assert( nrows == mat.nrows() && ncols == mat.ncols() );
	if( nrows > 0 && ncols > 0 )
		os.write( reinterpret_cast<const char*>(mat.data()), nrows*ncols*sizeof(T) );
}

template <typename Container>
void readArray( Container& array, std::istream& is )
{
	typedef typename Container::value_type T;
	int size;
	is.read( reinterpret_cast<char*>(&size), sizeof(int) );
	array.resize( size );

	if( size > 0 )
		is.read( reinterpret_cast<char*>(&array[0]), size*sizeof(T) );
}

template <typename Container>
void readArray( Container& array, std::istream& is, int size )
{
    typedef typename Container::value_type T;
    array.resize( size );

	if( size > 0 )
	    is.read( reinterpret_cast<char*>(&array[0]), size*sizeof(T) );
}

template <typename T, typename Order>
void readMatrix( generic_matrix<T, Order>& mat, std::istream& is )
{
	int nrows;
	int ncols;
	is.read( reinterpret_cast<char*>(&nrows), sizeof(int) );
	is.read( reinterpret_cast<char*>(&ncols), sizeof(int) );

	mat.resize( nrows, ncols );
	if( nrows > 0 && ncols > 0 )
		is.read( reinterpret_cast<char*>(mat.data()), nrows*ncols*sizeof(T) );
}

template <typename T, typename Order>
void readMatrix( generic_matrix<T, Order>& mat, std::istream& is, int nrows, int ncols )
{
    mat.resize( nrows, ncols );
	if( nrows > 0 && ncols > 0 )
	    is.read( reinterpret_cast<char*>(mat.data()), nrows*ncols*sizeof(T) );
}

// first is Q, second is S
std::pair<vl::Mat3f, vl::Mat3f> polarDecomposition( const vl::Mat3f& matrix );
std::pair<vl::Mat3d, vl::Mat3d> polarDecomposition( const vl::Mat3d& matrix );

vl::Mat3f fromToRotation(const vl::Vec3f& from, const vl::Vec3f& to);

std::pair< fortran_matrix, std::vector<double> > symmetricJacobi( fortran_matrix A );

#endif
