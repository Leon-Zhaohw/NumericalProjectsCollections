#include "stdafx.h"

#ifndef USE_GSL
#include <mkl_lapack.h>
#endif

#include "twigg/interpolation.h"
#include "twigg/linalg.h"

Interpolator::Interpolator( const double* x, const double* y, int N
#ifdef USE_GSL
						   , const gsl_interp_type *T
#endif
						   )
#ifdef USE_GSL
	:	type_( T )
#ifndef CONDOR
		, lock_( new boost::mutex() )
#endif // CONDOR
#endif // USE_GSL
{
	min_ = *std::min_element( x, x + N );
	max_ = *std::max_element( x, x + N );

#ifdef USE_GSL
#ifndef CONDOR
	boost::mutex::scoped_lock scoped_lock(*lock_);
#endif // CONDOR

	const unsigned int minPoints = T->min_size;
	if( N < minPoints )
	{
		std::ostringstream oss;
		oss << "Iterpolation method '" << T->name << "' requires a minimum of " << minPoints << " points." 
			<< std::endl;
		throw CreationException( oss.str() );
	}

	accel_ = gsl_interp_accel_alloc();
	spline_ = gsl_spline_alloc(T, N);
	gsl_spline_init(spline_, x, y, N);
#else
	values_.reserve(N);
	for( unsigned int i = 0; i < N; ++i )
		values_.push_back( InterpValue( x[i], y[i] ) );

	// make sure the xs are in sorted order:
	std::sort( values_.begin(), values_.end() );

	if( values_.size() > 2 )
	{
		// form the tridiagonal matrix for the solve:
		std::vector<double> b( N-2 );

		std::vector<double> du( b.size()-1 );
		std::vector<double> dl( b.size()-1 );
		std::vector<double> d( b.size() );

		for( unsigned int j = 1; j < (N-1); ++j )
		{
			if( (j-1) != 0 )
				dl[j-2] = (x[j] - x[j-1]) / 6.0;

			d[j-1] = (x[j+1] - x[j-1]) / 3.0;

			if( (j+1) != (N-1) )
				du[j-1] = (x[j+1] - x[j]) / 6.0;

			b[j-1] = ((y[j+1] - y[j]) / (x[j+1] - x[j])) - 
				((y[j] - y[j-1]) / (x[j] - x[j-1]));
		}

		{

#ifdef NO_BLAS
			// Implemented off this document:
			// http://math.berkeley.edu/~syazdani/math128a-05s/Prog3.pdf

			// First, compute the LU decomposition
			for( size_t i = 0; i < du.size(); ++i )
			{
				dl[i] /= d[i];
				d[i+1] -= dl[i]*du[i];
			}

			// Now solve Ly = b
			for( size_t i = 0; i < dl.size(); ++i )
			{
				b[i+1] -= b[i]*dl[i];
			}

			// Now solve Ux = y
			b[ b.size() - 1 ] /= d[ b.size() - 1 ];
			for( size_t i = 0; i < b.size() - 1; ++i )
			{
				size_t iVec = b.size() - i;
				assert( iVec > 0 );
				b[iVec - 1] = (b[iVec - 1] - du[iVec - 1]*b[iVec]) / d[iVec - 1];
			}
#else
			int n = N-2;
			int nrhs = 1;
			int ldb = n;
			int info = 0;

			dgtsv(	&n,
					&nrhs,
					&dl[0],
					&d[0],
					&du[0],
					&b[0],
					&ldb,
					&info );
			if( info < 0 )
			{
				std::ostringstream oss;
				oss << "Tridiagonal solve failed due to invalid parameter #" 
					<< -info;
				throw LinearAlgebraException( oss.str() );
			}
			else if( info > 0 )
			{
				std::ostringstream oss;
				oss << "Tridiagonal solve failed because A(" << info << ", " << info << ") is exactly 0.";
				throw LinearAlgebraException( oss.str() );
			}
#endif
		}

		for( unsigned int j = 1; j < N-1; ++j )
			values_[j].d2y = b[j-1];
	}
#endif
}

Interpolator::~Interpolator()
{
#ifdef USE_GSL
	boost::mutex::scoped_lock scoped_lock(*lock_);

	gsl_spline_free(spline_);
	gsl_interp_accel_free(accel_);
#endif // USE_GSL
}

double Interpolator::operator()( double x ) const
{
	// GSL only does interpolation, won't do extrapolation unfortunately.
	x = std::min<double>( std::max<double>( x, min_ ), max_ );

#ifdef USE_GSL
	boost::mutex::scoped_lock scoped_lock(*lock_);

	return gsl_spline_eval(spline_, x, accel_);
#else // USE_GSL
	typedef InterpValueArray::const_iterator Iter;

	Iter nextIter = std::lower_bound( values_.begin(), values_.end(),
		InterpValue( x, 0.0 ) );

	if( nextIter == values_.begin() )
		return values_.front().y;
	if( nextIter == values_.end() )
		return values_.back().y;

	Iter prevIter = nextIter - 1;
	const double delta_x = nextIter->x - prevIter->x;
	const double A = (nextIter->x - x) / delta_x;
	const double B = (x - prevIter->x) / delta_x;
	const double C = (1.0/6.0)*(A*A*A - A)*(delta_x*delta_x);
	const double D = (1.0/6.0)*(B*B*B - B)*(delta_x*delta_x);

	return A*prevIter->y + 
		B*nextIter->y + 
		C*prevIter->d2y + 
		D*nextIter->d2y;
#endif
}

double Interpolator::deriv( double x ) const
{
#ifdef USE_GSL
	x = std::min<double>( std::max<double>( x, min_ ), max_ );

#ifndef CONDOR
	boost::mutex::scoped_lock scoped_lock(*lock_);
#endif // CONDOR
	return gsl_spline_eval_deriv(spline_, x, accel_);
#else
	typedef InterpValueArray::const_iterator Iter;

	Iter nextIter = std::lower_bound( values_.begin(), values_.end(),
		InterpValue( x, 0.0 ) );

	if( nextIter == values_.begin() )
		return 0.0;
	if( nextIter == values_.end() )
		return 0.0;

	Iter prevIter = nextIter - 1;
	const double delta_x = prevIter->x - nextIter->x;
	const double A = (nextIter->x - x) / delta_x;
	const double dA = (nextIter->x - 1.0) / delta_x;
	const double B = (x - prevIter->x) / delta_x;
	const double dB = (1.0 - prevIter->x) / delta_x;
	const double dC = (1.0/6.0)*(3.0*A*A*dA - dA)*(delta_x*delta_x);
	const double dD = (1.0/6.0)*(3.0*B*B*dB - dB)*(delta_x*delta_x);

	return dA*prevIter->y + 
		dB*nextIter->y + 
		dC*prevIter->d2y + 
		dD*nextIter->d2y;
#endif
}

