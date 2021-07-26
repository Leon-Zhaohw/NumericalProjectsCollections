/*
 *	camera2.cpp
 *
 *	Created by Ryoichi Ando on Dec 25 2012
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include <math.h>
#include "camera2.h"
#include "opengl.h"

void camera2::lookAt( vec2d target, vec2d origin, vec2d up, FLOAT64 fov ) {
	this->target = target;
	this->origin = origin;
	this->up = up;
	this->fov = fov / 180.0 * PI;
	
	forward = (target-origin).normal();
	vec2d up_tmp = forward.rotate();
	if( up_tmp * up < 0 ) up_tmp = -1.0 * up_tmp;
	this->up = up_tmp;
}

vec2d camera2::getProjectedCoord( vec2d pos ) const {
	pos -= origin;
	FLOAT64 z_coord = forward * pos;
	return vec2d( z_coord, up*pos / (z_coord*tan(fov)) );
}

void camera2::draw() const {
	// Draw camera center
	uint dth = 20;
	FLOAT64 r = 0.01;
	glColor4d(1.0,0.6,0.6,1.0);
	glBegin(GL_LINE_LOOP);
	for( int t=0; t<dth; t++ ) {
		FLOAT64 theta = 2.0*t*PI/dth;
		vec2d rpos( r*cos(theta),r*sin(theta) );
		glVertex2dv((rpos+origin).v);
	}
	glEnd();
	
	// Draw a fov line
	FLOAT64 len = (origin-target).len();
	glBegin(GL_LINE_LOOP);
	glColor4d(1.0,1.0,1.0,0.8);
	glVertex2dv(origin.v);
	glVertex2dv((origin+len*forward+len*tan(fov)*up).v);
	glVertex2dv((origin+len*forward-len*tan(fov)*up).v);
	glEnd();
	glBegin(GL_LINES);
	glColor4d(1.0,0.6,0.6,0.8);
	glVertex2dv(origin.v);
	glVertex2dv(target.v);
	glColor4d(0.6,0.6,1.0,0.8);
	glVertex2dv(origin.v);
	glVertex2dv((origin+0.05*up).v);
	glEnd();
}