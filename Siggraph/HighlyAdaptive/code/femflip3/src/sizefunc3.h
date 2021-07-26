/*
 *	sizefunc3.h
 *
 *	Created by Ryoichi Ando on Dec 25 2012
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "vec3.h"
#include "kernel.h"
#include "util.h"

#ifndef _SIZEFUNC3_H
#define _SIZEFUNC3_H

// Size functions with embedded magic numbers...
//
inline FLOAT64 curvature_func( FLOAT64 curvature ) {
	return 0.8 / fmax(fabs(curvature),1e-18);
}

inline FLOAT64 diff_curvature_func( FLOAT64 diff_curvature ) {
	return 10.0 / fmax(fabs(diff_curvature),1e-18);
}

inline FLOAT64 strain_func( FLOAT64 strain ) {
	return 30.0 / fmax(strain,1e-18);
}

inline FLOAT64 solid_func( FLOAT64 solid_levelset, FLOAT64 solid_curvature, FLOAT64 min_dx, FLOAT64 max_dx ) {
	return 1.6 / fmax(fabs(solid_curvature),1e-18) / fmax(1e-18,kernel::smooth_kernel(sqr(solid_levelset),max_dx));
}

inline FLOAT64 frustum_func( vec3d camera_coord, FLOAT64 dx, FLOAT64 min_dx, FLOAT64 max_dx,
							 vec3d pos, vec3d origin, vec3d target, bool visible ) {
	FLOAT64 dist = camera_coord[2];
	if( dist < 0.0 ) return max_dx;
	if( fabs(camera_coord[0]) > 1.0 || fabs(camera_coord[1]) > 1.0 ) { // for stanford party example
		return max_dx;
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
// origin:				View origin
// target:				View target
// screen:				Camera screen coordinate (x,y): screen coord (-1~1,-1~1) z: depth (0~1)
// visible:				Is visible from a camera
// solid_levelset:		Solid levelset ( negative = inside, positive = outside )
// solid_curvature:		Curvature of solid levelset
// min_dx:				Minimal BCC grid cell that can be expressed
// max_dx:				Maximum BCC grid cell size allowed on surfaces
// startup:				Only true at the initialization stage
// -----------------------------------------------------
static FLOAT64 size_func( vec3d pos, FLOAT64 levelset, FLOAT64 curvature, FLOAT64 diff_curvature, FLOAT64 strain,
						  vec3d origin, vec3d target, vec3d screen,  bool visible, FLOAT64 solid_levelset, FLOAT64 solid_curvature,
						  FLOAT64 min_dx, FLOAT64 max_dx, bool startup ) {
	// Ensure that dx is smaller than max_dx
	FLOAT64 dx = max_dx;
	
	// Curvature
	dx = fmin(dx,curvature_func(curvature));

	// Strain
	dx = fmin(dx,strain_func(strain));
	
	// For floating pool test
	dx = fmin(dx,0.5*(pos-vec3d(0.5,0.5,0.5)).len());
	
	// Solid function
	// dx = fmin(dx,solid_func(solid_levelset,solid_curvature,min_dx,max_dx));

	// View frustum
	// dx = frustum_func(screen,dx,min_dx,max_dx,pos,origin,target,visible);
	
	// Depth (allowed to exceed max_dx)
	dx = fmax(dx,depth_func(fabs(levelset)));
	
	return dx;
}

#endif