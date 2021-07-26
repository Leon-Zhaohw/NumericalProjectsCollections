/*
 *	opengl.h
 *	
 *	Created by Ryoichi Ando on 11/20/11
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

static void glVertex3dv( const long double *v ) {
	double v1 = v[0];
	double v2 = v[1];
	double v3 = v[2];
	glVertex3d(v1,v2,v3);
}

static void glRasterPos3dv( const long double *v ) {
	double v1 = v[0];
	double v2 = v[1];
	double v3 = v[2];
	glRasterPos3d(v1,v2,v3);
}

#endif