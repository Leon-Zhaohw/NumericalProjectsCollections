/*
 *	levelset3.h
 *	
 *	Created by Ryoichi Ando on 11/26/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "macros.h"
#include "vec3.h"
#include "array3.h"
#include "octree3.h"

#ifndef _LEVELSET3_H
#define _LEVELSET3_H

class particle3;
class levelset3 {
public:
	levelset3() {
		dx = 1e-8;
		gsize = 0;
	}
	virtual void resize( uint gsize );
	virtual FLOAT64 evalLevelset(vec3d p) const = 0;
	virtual FLOAT64 evalCurvature(vec3d p) const {
		return (evalLevelset(p+vec3d(dx,0,0))+evalLevelset(p-vec3d(dx,0,0))+
				evalLevelset(p+vec3d(0,dx,0))+evalLevelset(p-vec3d(0,dx,0))+
				evalLevelset(p+vec3d(0,0,dx))+evalLevelset(p-vec3d(0,0,dx))
				-6*evalLevelset(p))/(dx*dx);
	}
	virtual vec3d evalGradient(vec3d p) const;
	virtual bool getClosestSurfacePos(vec3d &pos) const {
		pos = pos - evalLevelset(pos)*evalGradient(pos).normal();
		return true;
	}
	virtual void draw() const {};
	
	uint gsize;
	FLOAT64 dx;
	// Optional
	virtual const std::vector<octree3::sphere3 *> *getSpheres () const { return NULL; }
	
	// Utility functions...
	static FLOAT64 box( vec3d p, vec3d p0, vec3d p1 );
	static FLOAT64 sphere( vec3d p, vec3d c, FLOAT64 r );
	
	FLOAT64 getParticleConvextHullLevelset( vec3d &pos, const particle3 &p0, FLOAT64 dpx, FLOAT64 shrink, bool consistent ) const;
	FLOAT64 getParticleConvextHullLevelset( vec3d &pos, const particle3 &p0, const particle3 &p1, FLOAT64 dpx, FLOAT64 shrink, bool consistent ) const;
	FLOAT64 getParticleConvextHullLevelset( vec3d &pos, const particle3 &p0, const particle3 &p1, const particle3 &p2, FLOAT64 dpx, FLOAT64 shrink, bool consistent ) const;
	FLOAT64 getParticleLevelset( const std::vector<particle3 *> &neighbors, FLOAT64 dpx, vec3d pos, vec3d &outpos, const levelset3 *solid, FLOAT64 shrink=0.75, bool consistent=true ) const;
	FLOAT64 getAdams07ParticleLevelset( const std::vector<particle3 *> &neighbors, FLOAT64 dpx, vec3d pos, vec3d &outpos, const levelset3 *solid, FLOAT64 shrink=0.75 ) const;
};

template <class T> class gridLevelset3 : public levelset3 {
public:
	gridLevelset3() {
		centered = false;
		array = NULL;
	}
	void bind( const array3<T> *array ) {
		this->array = array;
	}
	void setCellCentered( bool centered ) {
		this->centered = centered;
	}
	virtual FLOAT64 evalLevelset(vec3d p) const {
		size3 size = array->size();
		T e = 1.0e-4;
		if( centered ) {
			p[0] = fmin(size.w-1-e,fmax(0,size.w*p[0]-0.5));
			p[1] = fmin(size.h-1-e,fmax(0,size.h*p[1]-0.5));
			p[2] = fmin(size.d-1-e,fmax(0,size.d*p[2]-0.5));
		} else {
			p[0] = (size.w-1)*fmin(1.0-e,fmax(0.0,p[0]));
			p[1] = (size.w-1)*fmin(1.0-e,fmax(0.0,p[1]));
			p[2] = (size.d-1)*fmin(1.0-e,fmax(0.0,p[2]));
		}
		int i = (int)p[0];
		T x = p[0]-i;
		int j = (int)p[1];
		T y = p[1]-j;
		int k = (int)p[2];
		T z = p[2]-k;
		const array3<T> &q = *array;
		return  z*(1.0-x)*(1.0-y)*q[i][j][k+1]+z*(1.0-x)*y*q[i][j+1][k+1]+z*x*y*q[i+1][j+1][k+1]+z*x*(1.0-y)*q[i+1][j][k+1] +
				(1.0-z)*(1.0-x)*(1.0-y)*q[i][j][k]+(1.0-z)*(1.0-x)*y*q[i][j+1][k]+(1.0-z)*x*y*q[i+1][j+1][k]+(1.0-z)*x*(1.0-y)*q[i+1][j][k];
	}
protected:
	const array3<T> *array;
	bool centered;
};

#endif