/*
 *	surf2.h
 *
 *	Created by Ryoichi Ando on 8/19/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "levelset2.h"
#include "meshsurf2.h"
#include "bccsurf2.h"
#include "cellsurf2.h"
#include "ann2.h"

#ifndef _SURF2_H
#define _SURF2_H

class fluid2;
class surf2 : public levelset2 {
public:
	surf2();
	surf2( const surf2 &surf );
	void setResolution( uint gm );
	void setExtrapolationDist( FLOAT64 dist );
	void buildSurface(const levelset2 *fluid, const levelset2 *solid, fluid2 *fluidSolver, bool enclose=false );
	void drawGL( const levelset2 *solid, uint kind ) const;
	virtual FLOAT64 evalLevelset(vec2d p) const;
	virtual vec2d evalGradient(vec2d p) const;
	virtual FLOAT64 evalCurvature(vec2d p) const;
	virtual void operator=( const surf2 &surf );
	
	levelset2 *surf;
	// Basic mesh information
	std::vector<vec2d> vertices;
	std::vector<vec2d> normals;
	std::vector<vec2i> faces;
	std::vector<FLOAT64> curvature;
	std::vector<FLOAT64> areas;
	std::vector<std::vector<uint> > v2f;
	
	void computeCurvature( const levelset2 *solid );
protected:
	mesher2 mesh;						// Mesh for mesh surface
	meshsurf2 meshsurf;					// Mesh surface
	bccsurf2 bccsurf;					// BCC surface
	cellsurf2 cellsurf;					// Cell surface
	uint gm;
	void fitSurface( const levelset2 *levelset );
	ann2 ann;
};

#endif