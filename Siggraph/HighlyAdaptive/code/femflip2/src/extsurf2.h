/*
 *	;
 *
 *	Created by Ryoichi Ando on Dec 15 2012
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include <vector>
#include "octree2.h"
#include "bcc2.h"
#include "mesher2.h"
#include "util2.h"
#include "levelset2.h"
#include "ann2.h"
#include "svd2.h"

#ifndef _EXTSURF2_H
#define _EXTSURF2_H

class extsurf2 : public levelset2 {
public:
	extsurf2();
	void loadParticles( const std::vector<particle2 *> &particles, FLOAT64 dpx, const octree2 *octree=NULL );
	void setMethod( int method );
	void generateMesh( levelset2 *solid );
	virtual FLOAT64 evalLevelset(vec2d p) const;
	virtual bool getClosestSurfacePos(vec2d &pos) const;
	void draw();
	
	std::vector<vec2d> vertices;
	std::vector<vec2d> normals;
	std::vector<vec2i> faces;
	
	typedef struct {
		vec2d p;
		FLOAT64 r;
		FLOAT64 remesh_r;
		uint n;
		FLOAT64 levelset;
		vec2d gradient;
		bool isolated;
		int depth;
	} sphere2;
protected:
	std::vector<sphere2> spheres;
	std::vector<uint> surface_spheres;
	// For anisotropic method...
	std::vector<svd2> anisotropy;
	std::vector<FLOAT64> density;
	
	uint numSample;
	int method;
	mesher2 mesher;
	octree2 octree;
	bcc2 bcc;
	FLOAT64 dpx;
	FLOAT64 kernel_size;
	levelset2 *solid;
	ann2 sorter;
	ann2 surf_sorter;
private:
	particle2 genParticle( const sphere2 &sphere ) const;
};

#endif
