#ifdef WIN32
#pragma once
#endif

#ifdef USE_GSL
#include <gsl/gsl_spline.h>
#endif

#ifndef __INTERPOLATION_H__
#define __INTERPOLATION_H__

class Interpolator
{
public:
	Interpolator( const double* x, const double* y, int N
#ifdef USE_GSL
		, const gsl_interp_type *T = gsl_interp_akima
#endif
		);
	~Interpolator();

	double operator()( double x ) const;
	double deriv( double x ) const;

private:
#ifdef USE_GSL
	const gsl_interp_type* type_;
	gsl_interp_accel* accel_;
	gsl_spline* spline_;
#else
	struct InterpValue
	{
		InterpValue( double x_, double y_ )
			: x(x_), y(y_), d2y(0.0) {}

		bool operator<( const InterpValue& rhs ) const
		{
			return x < rhs.x;
		}

		double x;
		double y;
		double d2y;
	};

	typedef std::vector<InterpValue> InterpValueArray;
	InterpValueArray values_;
#endif

	double min_;
	double max_;

#ifndef CONDOR
#ifdef USE_GSL
	mutable boost::scoped_ptr<boost::mutex> lock_;
#endif
#endif
};


#endif

