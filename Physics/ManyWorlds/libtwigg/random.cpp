#include "stdafx.h"

#include "twigg/random.h"

#include <boost/bind.hpp>
#include <fstream>
#include <ctime>

#ifdef USE_VSL
#include "twigg/linalg.h"
#include <mkl_vsl.h>
#else
#include <gsl/gsl_randist.h>
#endif

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

RandomStream::RandomStream(unsigned int seed)
{
	init(seed);
}


RandomStream::RandomStream()
{
	SeedType seed = 0;

	while( seed == 0 )
	{
	#ifdef WIN32
		HCRYPTPROV hProv = 0;
		if( CryptAcquireContext(&hProv, 0, 0, PROV_RSA_FULL, CRYPT_VERIFYCONTEXT) )
		{
			CryptGenRandom(hProv, sizeof(SeedType), reinterpret_cast<BYTE*>( &seed ));
		}
		else
		{
			// this is okay because we don't really need it to be cryptographically secure
			std::cerr << "Warning: unable to obtain random seed." << std::endl;
			seed = time(NULL);
		}
	#else
		std::ifstream ifs( "/dev/random", std::ios::binary );
		ifs.read( reinterpret_cast<char*>( &seed ), sizeof( SeedType ) );
	#endif
	}

	init(seed);
}

void RandomStream::init(RandomStream::SeedType seed)
{
#ifdef USE_VSL
	const int brng = VSL_BRNG_MRG32K3A;
	vslNewStream( &state_, brng, seed );
#elif defined(USE_GSL)
	static const gsl_rng_type *T = gsl_rng_env_setup();
	state_ = gsl_rng_alloc(T);
	gsl_rng_set( state_, seed );
#else
	if( seed < 0 )
		seed = -seed;

	state_.seed( seed );
#endif
}

RandomStream::~RandomStream()
{
#ifdef USE_VSL
	vslDeleteStream( &state_ );
#elif defined(USE_GSL)
	gsl_rng_free( state_ );
#endif
}

void RandomStream::uniform( double a, double b, int n, double* out )
{
#ifdef USE_VSL
	vdRngUniform( 0, state_, n, out, a, b );
#elif defined(USE_GSL)
	for( int i = 0; i < n; ++i )
		out[i] = a + (b-a)*gsl_rng_uniform( state_ );
#else
	typedef boost::uniform_real<double> RandType;
	boost::variate_generator< RngType&, boost::uniform_real<> > rng( 
		this->state_, RandType(a, b) );
	std::generate( out, out+n, rng );
#endif
}

double RandomStream::uniform( double a, double b )
{
#ifdef USE_VSL
	double result;
	vdRngUniform( 0, state_, 1, &result, a, b );
	return result;
#elif defined(USE_GSL)
	return a + (b-a)*gsl_rng_uniform( state_ );
#else
	typedef boost::uniform_real<double> RandType;
	boost::variate_generator< RngType&, boost::uniform_real<> > rng( 
		this->state_, RandType(a, b) );
	return rng();
#endif
}

double RandomStream::uniform()
{
#ifdef USE_GSL
	return gsl_rng_uniform( state_ );
#else
	return uniform(0.0, 1.0);
#endif
}

void RandomStream::uniform( int a, int b, int n, int* out )
{
#ifdef USE_VSL
	viRngUniform( 0, state_, n, out, a, b );
#elif defined(USE_GSL)
	for( int i = 0; i < n; ++i )
		out[i] = a + gsl_rng_uniform_int( state_, b-a );
#else
	// should not be inclusive
	typedef boost::uniform_int<int> RandType;
	boost::variate_generator< RngType&, RandType > rng( 
		this->state_, RandType(a, b-1) );
	std::generate( out, out+n, rng );
#endif
}

int RandomStream::uniform( int a, int b )
{
#ifdef USE_VSL
	int result;
	viRngUniform( 0, state_, 1, &result, a, b );
	return result;
#elif defined(USE_GSL)
	return a + gsl_rng_uniform_int( state_, b-a );
#else
	typedef boost::uniform_int<int> RandType;
	boost::variate_generator< RngType&, RandType > rng( 
		this->state_, RandType(a, b-1) );
	return rng();
#endif
}

void RandomStream::uniformBits( int n, unsigned int* out )
{
#ifdef USE_VSL
	viRngUniformBits( 0, state_, n, out );
#elif defined(USE_GSL)
	std::generate( out, out+n, boost::bind( &gsl_rng_get, state_ ) );
#else
	std::generate( out, out+n, this->state_ );
#endif
}

unsigned int RandomStream::uniformBits()
{
#ifdef USE_VSL
	unsigned int result;
	viRngUniformBits( 0, state_, 1, &result );
	return result;
#elif defined(USE_GSL)
	return gsl_rng_get( state_ );
#else
	return this->state_();
#endif
}

unsigned int RandomStream::binomial( unsigned int n, double p )
{
#ifdef USE_VSL
	int result;
	viRngBinomial( 0, state_, 1, &result, n, p );
	return boost::numeric_cast<unsigned int>(result);
#elif defined(USE_GSL)
	return gsl_ran_binomial( state_, p, n );
#else
	typedef boost::binomial_distribution<unsigned int> RandType;
	boost::variate_generator< RngType&, RandType > rng( this->state_, RandType(n, p) );
	return rng();
#endif
}

unsigned int RandomStream::poisson( double mu )
{
#ifdef USE_VSL
	int result;
	viRngPoisson( 0, state_, 1, &result, mu );
	return boost::numeric_cast<unsigned int>(result);
#elif defined(USE_GSL)
	return gsl_ran_poisson( state_, mu );
#else
	typedef boost::poisson_distribution<unsigned int, double > RandType;
	boost::variate_generator< RngType&, RandType > rng( this->state_, RandType(mu) );
	return rng();
#endif
}

void RandomStream::poisson( double mu, int n, int* out )
{
#ifdef USE_VSL
	viRngPoisson( 0, state_, n, out, mu );
#elif defined(USE_GSL)
	std::generate( out, out+n, boost::bind( &gsl_ran_poisson, state_, mu ) );
#else
	typedef boost::poisson_distribution<int, double > RandType;
	boost::variate_generator< RngType&, RandType > rng( this->state_, RandType(mu) );
	std::generate( out, out+n, rng );
#endif
}

void RandomStream::normal( double mean, double stddev, int n, double* out )
{
#ifdef USE_VSL
	vdRngGaussian( 0, state_, n, out, mean, stddev );
#elif defined(USE_GSL)
	std::generate( out, out+n, boost::bind( &gsl_ran_gaussian, state_, stddev ) );
	std::transform( out, out+n, out, boost::bind( std::plus<double>(), _1, mean ) );
#else
	typedef boost::normal_distribution<double> RandType;
	boost::variate_generator< RngType&, RandType > rng( this->state_, RandType(mean, stddev) );
	std::generate( out, out+n, rng );
#endif
}

void RandomStream::normal( double mean, double stddev, int n, float* out )
{
#ifdef USE_VSL
	vsRngGaussian( 0, state_, n, out, mean, stddev );
#elif defined(USE_GSL)
	std::generate( out, out+n, boost::bind( &gsl_ran_gaussian, state_, stddev ) );
	std::transform( out, out+n, out, boost::bind( std::plus<double>(), _1, mean ) );
#else
	typedef boost::normal_distribution<float> RandType;
	boost::variate_generator< RngType&, RandType > rng( this->state_, RandType(mean, stddev) );
	std::generate( out, out+n, rng );
#endif
}

double RandomStream::normal( double mean, double stddev )
{
#ifdef USE_VSL
	double result;
	vdRngGaussian( 0, state_, 1, &result, mean, stddev );
	return result;
#elif defined(USE_GSL)
	return mean + gsl_ran_gaussian( state_, stddev );
#else
	typedef boost::normal_distribution<float> RandType;
	boost::variate_generator< RngType&, RandType > rng( this->state_, RandType(mean, stddev) );
	return rng();
#endif
}

std::vector<double> RandomStream::uniformSpherical( unsigned int n )
{
	std::vector<double> result(n);

#ifdef USE_VSL
	// we use the property that a multivariate Gaussian is symmetric
	vdRngGaussian(
		0, 
		state_, 
		result.size(), 
		&result[0],
		0.0, 
		1.0 );
	double norm = norm2( result );
	scale( result, 1.0/sqrt(norm) );
	return result;
#elif defined(USE_GSL)
	gsl_ran_dir_nd( state_, n, &result[0] );
	return result;
#else
	typedef boost::uniform_on_sphere<double> RandType;
	boost::variate_generator< RngType&, RandType > rng( this->state_, RandType(n) );
	return rng();
#endif
}


