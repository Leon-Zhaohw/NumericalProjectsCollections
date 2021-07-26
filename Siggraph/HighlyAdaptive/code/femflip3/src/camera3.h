/*
 *	camera3.h
 *
 *	Created by Ryoichi Ando on Dec 25 2012
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "vec3.h"

#ifndef _CAMERA3_H
#define _CAMERA3_H

class camera3 {
public:
	camera3();
	virtual void lookAt( vec3d target, vec3d origin, vec3d up, FLOAT64 fov, FLOAT64 aspect );
	virtual void setEnabled( bool enabled );
	virtual bool isEnabled() const { return enabled; }
	virtual bool isStatic() const { return true; }
	vec3d getProjectedCoord( vec3d pos ) const;
	void draw() const;

	vec3d target;
	vec3d origin;
	vec3d up;
	FLOAT64 fov;
	FLOAT64 aspect;
	///
	vec3d forward;
	vec3d side;
	bool enabled;
};

#endif
