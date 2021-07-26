/*
 *	surf3.h
 *
 *	Created by Ryoichi Ando on 8/21/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "levelset3.h"
#include "meshsurf3.h"
#include "bccsurf3.h"
#include "cellsurf3.h"
#include "ann3.h"

#ifndef _SURF3_H
#define _SURF3_H

class fluid3;
class surf3 : public levelset3 {
public:
	surf3();
	surf3( const surf3 &surf );
	void setResolution( uint gm );
	void setExtrapolationDist( FLOAT64 dist );
	void buildSurface(const fluid3 *fluidSolver, const levelset3 *fluid, const levelset3 *solid, bool enclose=false );
	void drawGL( const levelset3 *solid, uint kind ) const;
	virtual FLOAT64 evalLevelset(vec3d p) const;
	virtual vec3d evalGradient(vec3d p) const;
	virtual FLOAT64 evalCurvature(vec3d p) const;
	virtual void operator=( const surf3 &surf );
	
	levelset3 *surf;
	// Basic surface information
	std::vector<vec3d> vertices;
	std::vector<vec3d> normals;
	std::vector<vec3i> faces;
	std::vector<vec3d> curvature;
	std::vector<FLOAT64> ringArea;
	std::vector<std::vector<uint> > v2f;
	
	void computeCurvature( const levelset3 *solid );
protected:
	mesher3 mesh;						// Mesh for mesh surface
	meshsurf3 meshsurf;					// Mesh surface
	bccsurf3 bccsurf;					// BCC surface
	cellsurf3 cellsurf;					// Cell surface
	uint gm;
	ann3 ann;
};

#endif