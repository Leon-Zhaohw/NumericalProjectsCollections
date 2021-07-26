/*
 *	camera3.cpp
 *
 *	Created by Ryoichi Ando on Dec 25 2012
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include <math.h>
#include "camera3.h"
#include "opengl.h"

camera3::camera3() {
	enabled = false;
}

void camera3::lookAt( vec3d target, vec3d origin, vec3d up, FLOAT64 fov, FLOAT64 aspect ) {
	forward = (target-origin).normal();
	side = (forward ^ up).normal();
	up = side ^ forward;
	this->target = target;
	this->origin = origin;
	this->up = up;
	this->fov = fov / 180.0 * PI;
	this->aspect = aspect;
}

void camera3::setEnabled( bool enabled ) {
	this->enabled = enabled;
}

vec3d camera3::getProjectedCoord( vec3d pos ) const {
	pos -= origin;
	FLOAT64 x = side * pos;
	FLOAT64 y = up * pos;
	FLOAT64 z = forward * pos;
	FLOAT64 z_tan_fov = z*tan(fov);
	x = 2.0 * x / z_tan_fov;
	y = 2.0 * aspect * y / z_tan_fov;
	return vec3d(x,y,z);
}

void camera3::draw() const {
	// Draw center
	glColor4d(1.0,0.5,0.5,1.0);
	glPushMatrix();
	glTranslated(origin[0],origin[1],origin[2]);
	glutWireSphere(0.01,20,20);
	glPopMatrix();
	
	// Draw view frustum
	glColor4d(1.0,0.8,0.8,1.0);
	glBegin(GL_LINES);
	glVertex3dv(origin.v);
	glVertex3dv(target.v);
	glEnd();
	
	FLOAT64 z = (origin-target).len();
	FLOAT64 z_tan_fov = z*tan(fov);
	vec3d corner0 = target-0.5*side*z_tan_fov-0.5*up*z_tan_fov/aspect;
	vec3d corner1 = target+0.5*side*z_tan_fov-0.5*up*z_tan_fov/aspect;
	vec3d corner2 = target+0.5*side*z_tan_fov+0.5*up*z_tan_fov/aspect;
	vec3d corner3 = target-0.5*side*z_tan_fov+0.5*up*z_tan_fov/aspect;
	glColor4d(1.0,1.0,1.0,0.8);
	glBegin(GL_LINE_LOOP);
	glVertex3dv(corner0.v);
	glVertex3dv(corner1.v);
	glVertex3dv(corner2.v);
	glVertex3dv(corner3.v);
	glEnd();
	
	glBegin(GL_LINES);
	glVertex3dv(origin.v);
	glVertex3dv(corner0.v);
	glVertex3dv(origin.v);
	glVertex3dv(corner1.v);
	glVertex3dv(origin.v);
	glVertex3dv(corner2.v);
	glVertex3dv(origin.v);
	glVertex3dv(corner3.v);
	glEnd();
}
