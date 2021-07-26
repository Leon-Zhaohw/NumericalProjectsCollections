/*
 *
 *	Created by Ryoichi Ando on Dec 15 2012
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include <vector>
#include "octree3.h"
#include "mesher3.h"
#include "bcc3.h"
#include "util3.h"
#include "levelset3.h"
#include "ann3.h"
#include "svd3.h"

#ifndef _EXTSURF3_H
#define _EXTSURF3_H

class extsurf3 : public levelset3 {
public:
	extsurf3();
	void writeFile( const char *path, const std::vector<particle3 *> &particles, FLOAT64 dpx, const octree3 *octree=NULL );
	void loadFile( const char *path );
	void loadParticles( const std::vector<particle3 *> &particles, FLOAT64 dpx, const octree3 *octree=NULL );
	void setMethod( int method );
	void writeMesh( uint frame, levelset3 *solid );
	virtual FLOAT64 evalLevelset(vec3d p) const;
	virtual bool getClosestSurfacePos(vec3d &pos) const;
	void draw();
	
	std::vector<vec3d> vertices;
	std::vector<vec3d> normals;
	std::vector<vec3i> faces;
	
	typedef struct {
		vec3d p;
		FLOAT64 r;
		FLOAT64 remesh_r;
		uint n;
		FLOAT64 levelset;
		vec3d gradient;
		bool isolated;
		int depth;
	} sphere3;
	std::vector<sphere3> spheres;
	std::vector<uint> surface_spheres;
	
	// For anisotropic method...
	std::vector<svd3> anisotropy;
	std::vector<FLOAT64> density;
	
	uint numSample;
	int method;
	mesher3 mesher;
	octree3 octree;
	FLOAT64 dpx;
	FLOAT64 kernel_size;
	levelset3 *solid;
	ann3 sorter;
	ann3 surf_sorter;
private:
	particle3 genParticle( const sphere3 &sphere ) const;
};

#endif
