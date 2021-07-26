/*
 *	levelset2.h
 *	
 *	Created by Ryoichi Ando on 11/3/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "macros.h"
#include "vec2.h"
#include "array2.h"
#include "octree2.h"

#ifndef _LEVELSET2_H
#define _LEVELSET2_H

class particle2;
class levelset2 {
public:
	levelset2() {
		dx = 1e-4;
		gsize = 0;
	}
	virtual void resize( uint gsize );
	virtual FLOAT64 evalLevelset(vec2d p) const = 0;
	virtual FLOAT64 evalCurvature(vec2d p) const {
		return (evalLevelset(p+vec2d(dx,0))+evalLevelset(p-vec2d(dx,0))+
				evalLevelset(p+vec2d(0,dx))+evalLevelset(p-vec2d(0,dx))
				-4*evalLevelset(p))/(dx*dx);
	}
	virtual vec2d evalGradient(vec2d p) const;
	virtual bool getClosestSurfacePos(vec2d &pos) const {
		pos = pos - evalLevelset(pos)*evalGradient(pos).normal();
		return true;
	}
	// Optional
	virtual const std::vector<octree2::sphere2 *> *getSpheres () const { return NULL; }
	uint gsize;
	FLOAT64 dx;
public:
	// Utility functions...
	static FLOAT64 box( vec2d p, vec2d p0, vec2d p1 );
	static FLOAT64 circle( vec2d p, vec2d c, FLOAT64 r );
	static void marchPoints( const std::vector<vec2d> &nodes, const std::vector<FLOAT64> &levelsets, std::vector<vec2d> &points, bool fill=true );
    static void marchPoints( FLOAT64 L[2][2], vec2d p[8], int &pnum );
	static void marchPoints( int i, int j, const array2<FLOAT64> &L, vec2d p[8], int &pnum );
	FLOAT64 getParticleConvextHullLevelset( vec2d &pos, const particle2 &p0, FLOAT64 dpx, FLOAT64 shrink, bool consistent ) const;
	FLOAT64 getParticleConvextHullLevelset( vec2d &pos, const particle2 &p0, const particle2 &p1, FLOAT64 dpx, FLOAT64 shrink, bool consistent ) const;
	FLOAT64 getParticleLevelset( const std::vector<particle2 *> &neighbors, FLOAT64 dpx, vec2d pos, vec2d &outpos, const levelset2 *solid, FLOAT64 shrink=0.75, bool consistent=true ) const;
};

template <class T> class gridLevelset2 : public levelset2 {
public:
	gridLevelset2() {
		centered = false;
		array = NULL;
	}
	void bind( const array2<T> *array ) {
		this->array = array;
	}
	void setCellCentered( bool centered ) {
		this->centered = centered;
	}
	virtual FLOAT64 evalLevelset(vec2d p) const {
		size2 size = array->size();
		T e = 1.0e-4;
		if( centered ) {
			p[0] = fmin(size.w-1-e,fmax(0,size.w*p[0]-0.5));
			p[1] = fmin(size.h-1-e,fmax(0,size.h*p[1]-0.5));
		} else {
			p[0] = (size.w-1)*fmin(1.0-e,fmax(0.0,p[0]));
			p[1] = (size.w-1)*fmin(1.0-e,fmax(0.0,p[1]));
		}
		int i = (int)p[0];
		T x = p[0]-i;
		int j = (int)p[1];
		T y = p[1]-j;
		const array2<T> &q = *array;
		return (1.0-x)*(1.0-y)*q[i][j]+(1.0-x)*y*q[i][j+1]+x*y*q[i+1][j+1]+x*(1.0-y)*q[i+1][j];
	}
protected:
	const array2<T> *array;
	bool centered;
};

#endif