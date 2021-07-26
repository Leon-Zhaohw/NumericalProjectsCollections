/*
 *	flycamera3.h
 *
 *	Created by Ryoichi Ando on Dec 30 2012
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "vec3.h"
#include "camera3.h"

#ifndef _FLYCAMERA3_H
#define _FLYCAMERA3_H

class flycamera3 : public camera3 {
public:
	virtual void setTime( FLOAT64 time ) = 0;
	virtual bool isStatic() const { return false; }
};

class static_camera3 : public flycamera3 {
public:
	static_camera3() {
		lookAt(vec3d(0.5,0.3,0.5), vec3d(0.5,2.5,-1.3), vec3d(0.0,1.0,0.0), 32.0, 1280 / 720.0 );
	}
	virtual void setTime( FLOAT64 time ) {
	}
	virtual bool isStatic() const { return true; }
};

#endif