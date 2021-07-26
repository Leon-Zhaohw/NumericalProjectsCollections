/*
 *	camera2.h
 *
 *	Created by Ryoichi Ando on Dec 25 2012
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "vec2.h"

#ifndef _CAMERA2_H
#define _CAMERA2_H

class camera2 {
public:
	void lookAt( vec2d target, vec2d origin, vec2d up, FLOAT64 fov );
	vec2d getProjectedCoord( vec2d pos ) const;
	void draw() const;
	
	vec2d target;
	vec2d origin;
	vec2d up;
	FLOAT64 fov;
	///
	vec2d forward;
};

#endif