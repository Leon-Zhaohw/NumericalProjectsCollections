#ifdef WIN32
#pragma once
#endif

#ifndef __PRT_H__
#define __PRT_H__

#include <boost/shared_ptr.hpp>

#include <deque>

class SphericalDirection
{
public:
	SphericalDirection( double phi, double theta )
		: phi_(phi), theta_(theta) {}

	double phi() const		{ return phi_; }
	double theta() const	{ return theta_; }

private:
	double phi_;
	double theta_;
};

SphericalDirection kayvon_sphericalDirectionForDirection( const vl::Vec3d direction );
double evaluateSphericalHarmonic( SphericalDirection dir, int l, int m );

class PRTEnvironment
{
public:
	typedef boost::array< std::vector<float>, 3 > SHCoeffArray;

	virtual ~PRTEnvironment() {}

	virtual SHCoeffArray shCoeffs(const vl::Mat4f& transform, unsigned int nCoeffs) const = 0;

	virtual void left()			{}
	virtual void right()		{}
	virtual void up()			{}
	virtual void down()			{}
};

class PRTSumEnvironment
	: public PRTEnvironment
{
public:
	typedef boost::shared_ptr<PRTEnvironment> PRTEnvPtr;
	virtual ~PRTSumEnvironment() {}
	void add( PRTEnvPtr env );

	SHCoeffArray shCoeffs(const vl::Mat4f& transform, unsigned int nCoeffs) const;

	void left();
	void right();
	void up();
	void down();

private:
	std::deque< PRTEnvPtr > children_;
};

class PRTDirectionalLight
	: public PRTEnvironment
{
public:
	PRTDirectionalLight( const vl::Vec3f color, double theta, double phi );
	virtual ~PRTDirectionalLight();

	SHCoeffArray shCoeffs(const vl::Mat4f& transform, unsigned int nCoeffs) const;
	void left();
	void right();
	void up();
	void down();

private:
	double theta_;
	double phi_;
	vl::Vec3f color_;
};

class PRTSkyLight
	: public PRTEnvironment
{
public:
	PRTSkyLight( const std::string& filename );
	virtual ~PRTSkyLight();

	SHCoeffArray shCoeffs(const vl::Mat4f& transform, unsigned int nCoeffs) const;
	void up();
	void down();
	void left();
	void right();

private:
	std::deque<PRTEnvironment::SHCoeffArray> coeffs_;
	int frame_;
};

class PRTAmbientLight
	: public PRTEnvironment
{
public:
	PRTAmbientLight( const vl::Vec3f& color );
	virtual ~PRTAmbientLight() {}

	SHCoeffArray shCoeffs(const vl::Mat4f& transform, unsigned int nCoeffs) const;

private:
	vl::Vec3f color_;
};

#endif

