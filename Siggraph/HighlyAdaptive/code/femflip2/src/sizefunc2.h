/*
 *	sizefunc2.h
 *
 *	Created by Ryoichi Ando on Dec 25 2012
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "vec2.h"
#include "util.h"
#include "kernel.h"

#ifndef _SIZEFUNC2_H
#define _SIZEFUNC2_H

// Size functions with embedded magic numbers...
//
inline FLOAT64 curvature_func( FLOAT64 curvature ) {
	return 0.1 / fmax(fabs(curvature),1e-18);
}

inline FLOAT64 diff_curvature_func( FLOAT64 diff_curvature ) {
	return 3.0 / fmax(fabs(diff_curvature),1e-18);
}

inline FLOAT64 strain_func( FLOAT64 strain ) {
	return 0.1 / fmax(strain,1e-18);
}

inline FLOAT64 solid_func( FLOAT64 solid_levelset, FLOAT64 solid_curvature, FLOAT64 min_dx, FLOAT64 max_dx ) {
	return 0.025 / fmax(fabs(solid_curvature),1e-18) / fmax(1e-18,kernel::smooth_kernel(sqr(solid_levelset),max_dx));
}

inline FLOAT64 frustum_func( vec2d camera_coord, FLOAT64 dx ) {
	if( fabs(camera_coord[1]) < 1.0 ) {
		dx = fabs(camera_coord[0]) * dx;
	}
	return dx;
}

inline FLOAT64 depth_func( FLOAT64 dist ) {
	return dist;
}

// Global size function (returns desired BCC cell size dx)
// ---------------- List of arguments ------------------
// pos:					Position to evaluate
// levelset:			Fluid levelset ( negative = inside, positive = outside )
// curvature:			Curvature of fluid levelset
// diff_curvature:		Time derivative of curvature
// strain:				Norm of strain tensor (always positive)
// screen:				Camera screen coordinate x: depth (0~1) y: screen (-1~1,-1~1)
// solid_levelset:		Solid levelset ( negative = inside, positive = outside )
// solid_curvature:		Curvature of solid levelset
// min_dx:				Minimal BCC grid cell that can be expressed
// max_dx:				Maximum BCC grid cell size allowed on surfaces
// startup:				Only true at the initialization stage
// -----------------------------------------------------
static FLOAT64 size_func( vec2d pos, FLOAT64 levelset, FLOAT64 curvature, FLOAT64 diff_curvature, FLOAT64 strain,
						  vec2d screen, FLOAT64 solid_levelset, FLOAT64 solid_curvature, FLOAT64 min_dx, FLOAT64 max_dx, bool startup ) {
	// Ensure that dx is smaller than max_dx
	FLOAT64 dx = max_dx;
	
	// Curvature
	dx = fmin(dx,curvature_func(curvature));
	
	// Time derivative curvature
	// dx = fmin(dx,diff_curvature_func(diff_curvature));
	
	// Strain
	dx = fmin(dx,strain_func(strain));
	
	// View frustum
	// dx = frustum_func(screen,dx);
	
	// Solid function
	// dx = fmin(dx,solid_func(solid_levelset,solid_curvature,min_dx,max_dx));
	
	// Depth (allowed to exceed max_dx)
	dx = fmax(dx,depth_func(fabs(levelset)));
	return dx;
}

#endif