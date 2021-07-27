#ifndef BUNNY_SETTINGS_H
#define BUNNY_SETTINGS_H

#include <cassert>
#include <cmath>
#include <float.h>

#define Real double
#define REAL_MAX DBL_MAX
#define REAL_MIN FLT_MIN

#ifndef NDEBUG
#define NDEBUG
#endif

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/freeglut.h>
#include <GL/glu.h>
#endif

#endif
