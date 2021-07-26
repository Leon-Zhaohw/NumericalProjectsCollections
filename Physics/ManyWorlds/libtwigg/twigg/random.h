#ifdef WIN32
#pragma once
#endif

#ifndef __RANDOM_H__
#define __RANDOM_H__

// For some reason VSL doesn't work on Linux, I 
//  think it might be just the version of the MKL
//  we're using, but I can use the GSL instead for
//  high quality random number streams
#ifdef USE_VSL
#include "mkl_vsl_functions.h"
#include <mkl_vsl_types.h>
#elif defined(USE_GSL)
#include <gsl/gsl_rng.h>
#else
#include <boost/random.hpp>
#endif

class RandomStream
{
public:
	// initialize with a prepared seed
	RandomStream(unsigned int seed);

	// generate a random seed using the system's entropy
	RandomStream();

	~RandomStream();

	// Generates numbers uniformly in [a, b)
	double uniform();
	double uniform( double a, double b );
	void uniform( double a, double b, int n, double* out );

	void uniform( int a, int b, int n, int* out );
	int uniform( int a, int b );

	unsigned int uniformBits();
	void uniformBits( int n, unsigned int* out );

	unsigned int poisson( double mu );
	void poisson( double mu, int n, int* out );

	unsigned int binomial( unsigned int n, double p );

	void normal( double mean, double stddev, int n, double* out );
	void normal( double mean, double stddev, int n, float* out );
	double normal( double mean = 0.0, double stddev = 1.0 );

	std::vector<double> uniformSpherical( unsigned int n );

private:
#ifdef USE_VSL
	typedef unsigned int SeedType;
#elif defined(USE_GSL)
	typedef unsigned long int SeedType;
#else
	typedef boost::hellekalek1995 RngType;
	typedef RngType::result_type SeedType;
#endif

	void init( SeedType seed );

#ifdef USE_VSL
	VSLStreamStatePtr state_;
#elif defined(USE_GSL)
	gsl_rng* state_;
#else
	RngType state_;
#endif
};

#endif

