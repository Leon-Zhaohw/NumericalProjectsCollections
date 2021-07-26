/*
 *	opengl.h
 *	
 *	Created by Ryoichi Ando on 11/3/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

// OpenGL Header

#if defined(__APPLE__) || defined(MACOSX) // For Mac OS X
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#elif defined(WIN32) // For Windows
#include <windows.h>
#include "glut.h"
#else // Mostly for Linux
#include <GL/gl.h>
#include <GL/glut.h>
#endif

#ifndef _OPENGL_H
#define _OPENGL_H

static void glVertex2dv( const long double *v ) {
	double v1 = v[0];
	double v2 = v[1];
	glVertex2d(v1,v2);
}

static void glRasterPos2dv( const long double *v ) {
	double v1 = v[0];
	double v2 = v[1];
	glRasterPos2d(v1,v2);
}

#endif