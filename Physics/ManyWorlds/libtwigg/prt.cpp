#include "stdafx.h"

#include "twigg/gslWrappers.h"
#include "twigg/prt.h"
#include "twigg/linalg.h"

#ifdef USE_GSL
#include <gsl/gsl_sf.h>
#else

#include <cmath>
/* maps (l,m) to an index i, where i = l(l+1)+m   */
#define LMINDEX(l,m)  ((l)*((l)+1)+(m))

/* increments to next (l,m) pair */
#define NEXTLM(l,m)  {if ((l)==(m)) {(l)++;(m)=-1*(l);} else {(m)++;}}
double  SH_coeff(int l, int m);
double  SH_eval(int l, int m, double theta, double psi);
double  legendre_eval(int l, int m, double x);
#endif

#include <fstream>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif




PRTDirectionalLight::PRTDirectionalLight( const vl::Vec3f color, double theta, double phi )
	:	color_(color),
		theta_(theta),
		phi_(phi)
{
}

PRTDirectionalLight::~PRTDirectionalLight()
{
}

PRTEnvironment::SHCoeffArray PRTDirectionalLight::shCoeffs(const vl::Mat4f& transform, unsigned int nCoeffs) const
{
	// Use y-up mode
	vl::Vec3d direction( 
		cos(theta_)*sin(phi_), 
		cos(phi_),
		sin(theta_)*sin(phi_) );

	vl::Mat3d rotMat;
	for( unsigned int i = 0; i < 3; ++i )
		for( unsigned int j = 0; j < 3; ++j )
			rotMat[i][j] = transform[i][j];
	rotMat = trans(rotMat);
	direction = vl::norm(rotMat * direction);

	SHCoeffArray result;
	unsigned int shCoeff = 0;
	for( int l = 0; shCoeff < nCoeffs; ++l )
	{
		for( int m = -l; m <= l; ++m )
		{
			double shmEval = evaluateSphericalHarmonic(
				kayvon_sphericalDirectionForDirection(direction), l, m);

			for( unsigned int iColor = 0; iColor < 3; ++iColor )
				result[iColor].push_back( color_[iColor] * shmEval );

			++shCoeff;
		}
	}

	return result;
}

static double thetaInc = M_PI/30;
static double phiInc = M_PI/30;

void PRTDirectionalLight::left()
{
	theta_ += thetaInc;
}

void PRTDirectionalLight::right()
{
	theta_ -= thetaInc;
}

void PRTDirectionalLight::up()
{
	phi_ = std::max<double>( 0.0, phi_ - phiInc );
}

void PRTDirectionalLight::down()
{
	phi_ = std::min<double>( M_PI, phi_+phiInc );
}

double evaluateSphericalHarmonic( SphericalDirection dir, int l, int m )
{
#ifdef USE_GSL
	if( m > 0 )
	{
		return M_SQRT2 
			* cos(m*dir.phi()) 
			* gsl_sf_legendre_sphPlm(l, m, cos(dir.theta()));
	}
	else if( m < 0 )
	{
		return M_SQRT2 
			* sin(-m*dir.phi()) 
			* gsl_sf_legendre_sphPlm(l, -m, cos(dir.theta()));
	}
	else
	{
		// m == 0
		return gsl_sf_legendre_sphPlm(l, 0, cos(dir.theta()));
	}
#else
	return SH_eval( l, m, dir.theta(), dir.phi() );
#endif
}

SphericalDirection kayvon_sphericalDirectionForDirection( const vl::Vec3d direction )
{
	double theta = acos( direction[2] );
	double phi = atan2(direction[1], direction[0]);

	return SphericalDirection( phi, theta );
}

PRTSkyLight::PRTSkyLight( const std::string& filename )
	: frame_(0)
{
	std::ifstream ifs( filename.c_str(), std::ios::binary );

	typedef generic_matrix< float, boost::numeric::ublas::column_major > matrix_type;
	matrix_type mat;
	readMatrix( mat, ifs );

	double scale = 1.0;

	// Now, dump all the coeffs in the appropriate spots
	for( unsigned int iFrame = 0; iFrame < mat.ncols(); ++iFrame )
	{
		coeffs_.push_back( PRTEnvironment::SHCoeffArray() );
		PRTEnvironment::SHCoeffArray& frame = coeffs_.back();

		for( unsigned int i = 0; i < 3; ++i )
			frame[i].reserve( mat.nrows()/3 );

		for( unsigned int iRow = 0; iRow < mat.nrows(); ++iRow )
			frame[iRow%3].push_back( scale*mat(iRow, iFrame) );
	}
}

PRTSkyLight::~PRTSkyLight()
{
}

PRTEnvironment::SHCoeffArray PRTSkyLight::shCoeffs(const vl::Mat4f& transform, unsigned int nCoeffs) const
{
	// not sure how to transform in this case
	assert( nCoeffs <= (this->coeffs_[frame_])[0].size() );
	return this->coeffs_[frame_];
}

void PRTSkyLight::left()
{
	this->frame_ = (frame_ + coeffs_.size() - 1) % coeffs_.size();
}

void PRTSkyLight::right()
{
	this->frame_ = (frame_ + 1) % coeffs_.size();
}

static int jump = 10;
void PRTSkyLight::up()
{
	this->frame_ = (frame_ + coeffs_.size() - jump) % coeffs_.size();
}

void PRTSkyLight::down()
{
	this->frame_ = (frame_ + jump) % coeffs_.size();
}

void PRTSumEnvironment::add( PRTEnvPtr env )
{
	children_.push_back( env );
}

PRTEnvironment::SHCoeffArray PRTSumEnvironment::shCoeffs(const vl::Mat4f& transform, unsigned int nCoeffs) const
{
	SHCoeffArray result;
	for( unsigned int i = 0; i < result.size(); ++i )
		result[i].resize( nCoeffs, 0.0 );
	for( std::deque< PRTEnvPtr >::const_iterator iter = children_.begin();
		iter != children_.end(); ++iter )
	{
		SHCoeffArray child = (*iter)->shCoeffs(transform, nCoeffs);
		for( unsigned int i = 0; i < result.size(); ++i )
		{
			result[i].resize( std::max<size_t>(child[i].size(), result[i].size()), 0.0 );
			std::transform( child[i].begin(), child[i].end(), 
				result[i].begin(), result[i].begin(), std::plus<float>() );
		}
	}

	return result;
}

void PRTSumEnvironment::left()
{
	std::for_each( children_.begin(), children_.end(), 
		boost::mem_fn( &PRTEnvironment::left ) );
}

void PRTSumEnvironment::right()
{
	std::for_each( children_.begin(), children_.end(), 
		boost::mem_fn( &PRTEnvironment::right ) );
}

void PRTSumEnvironment::up()
{
	std::for_each( children_.begin(), children_.end(), 
		boost::mem_fn( &PRTEnvironment::up ) );
}

void PRTSumEnvironment::down()
{
	std::for_each( children_.begin(), children_.end(), 
		boost::mem_fn( &PRTEnvironment::down ) );
}

PRTAmbientLight::PRTAmbientLight( const vl::Vec3f& color )
	: color_(color)
{
}

PRTEnvironment::SHCoeffArray PRTAmbientLight::shCoeffs(const vl::Mat4f& transform, unsigned int nCoeffs) const
{
	SHCoeffArray result;
	for( unsigned int i = 0; i < 3; ++i )
	{
		result[i].push_back( color_[i] );
		result[i].resize( nCoeffs, 0.0 );
	}

	return result;

}


// Kayvon code for computing spherical harmonics.  Swap in for
// the more robust GSL function when necessary.

/* returns true if given l & m values define a valid SH basis function */
#define SH_VALID(l,m)  ( (((l) >= 0) && (abs(m) <= (l))) ? 1 : 0 )

static int factorial(int n) {
  
  int total;
  
    if (n == 0)
        return 1;
    
    total = n--;
    
    while (n > 1) {
      total *= n;
      n--;
    }
    
    return total;
}

/* calculates n!!, the product of all off integers <= n
   KAYVON: what is value for 0!! ?? */
static int dblfactorial(int n) {

  int total;
  
  if (n == 0)
    return 1;
  
  if ( !(n & 1) )
    n--;
  
    total = n;
    n-=2;
    
    while (n > 1) {
      total *= n;
        n-=2;
    }
    
    return total;
}

/* this should be somewhat efficient, the worst case is P_l,0 */
double legendre_eval(int l, int m, double x) {
  
  int i, m2;
  double a, b, result;

  if (l == 0 && m == 0)
        return 1;
  
  /* force m > 0 when we calculate, we can exploit property
     relating P(l,m) to P(l,-m) at the end */
  
  m2 = (m < 0) ? -1 * m : m;
  
  /* we know closed form for P_ll(x)
     This is the only legendre we actually compute, we just build
     all the others by working up from these */
  if (l == m2) {
    
    /* (1-x^2)^(l/2) */
    result = pow(1.0 - (x*x), .5 * l);        
    
    /* multiply by (2l-1)!! */
    result *= dblfactorial((2*l) - 1);       
    
    /* multiply by (-1)^l */
    if (l & 1)
      result *= -1;
    
    }
  else {
    /* make use of recurrence properties of legendre polys to
       determine P_l,m */
    a = 0.0;
    b = legendre_eval(m2, m2, x);
        
    i = m2 + 1;                
    
    /* using P_lm = (1/(l-m)) * (x(2l-1)P_l,m_1 - (l+m-1)P_l,m-2 ) */
    while (i <= l) {
      
      result = ( (x * b * ((2*i)-1) ) - (a * (i+m2-1)) ) / ((double)(i-m2));
      a = b;
      b = result;
	    
      i++;
    }
  }
  
  
  /* use identity P_l,-m = (-1)^m * (l-m)!/(l+m)! * P_lm */
  if (m < 0) {
    
    result *= ((double)factorial(l-m2)) / ((double)(factorial(l+m2)));
    
    if (m2 & 1)
      result *= -1;
    
  }
  
  return result;
  
}


/* produces spherical harmonic coefficient */
double SH_coeff(int l, int m) {
  
  
  double top, bottom;
  
  if (m < 0)
    m *= -1;

  top = (double)( ((2*l)+1) * factorial(l-m) );
  bottom = 4.0 * M_PI * factorial(l+m);
  
  return sqrt( top / bottom);
  
}


/* evaluates the spherical harmonic at (theta,psi) */
double SH_eval(int l, int m, double theta, double psi) {
  
  
  int m2;
  
  if (!SH_VALID(l,m))
        return 0.0;
  
  m2 = (m < 0) ? -1*m : m;
  
  if (m2 == 0) {
    /* KAYVON: this assumes values are in RADIANS!!!!! */
    return SH_coeff(l,0) * legendre_eval(l,0, cos(theta) );
  }
  else if (m > 0)
    return sqrt(2.0) * SH_coeff(l,m2) * cos(psi*m2) * legendre_eval(l, m2, cos(theta) );
  else
    return sqrt(2.0) * SH_coeff(l,m2) * sin(psi*m2) * legendre_eval(l, m2, cos(theta) );    
  
}
