// SETTINGS.h: Project-wide options set in one place
//
//////////////////////////////////////////////////////////////////////

#ifndef SETTINGS_H
#define SETTINGS_H

#include <cassert>

#define USING_OSX __APPLE__
#define USING_PETSC 1
#define USING_GLVU 1

// uses OpenMP in SUBSPACE_INTEGRATOR
//#define USING_OPENMP 1 
//#define USING_RENDERMAN 1
//#define USING_MKL 1

// select single or double precision
#define DOUBLE_PRECISION
#define Real double
//#define SINGLE_PRECISION
//#define Real float

// turn off asserts?
#define NDEBUG

#endif
