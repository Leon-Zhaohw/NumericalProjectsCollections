/*
 *	cellsurf2.h
 *	
 *	Created by Ryoichi Ando on 2/11/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "levelset2.h"
#include "macros.h"
#include "vec2.h"
#include "array2.h"
#include <vector>

#ifndef _CELLSURF2_H
#define _CELLSURF2_H

class cellsurf2 : public levelset2 {
public:
	cellsurf2();
	cellsurf2( const cellsurf2 &cellsurf );
	void setResolution( uint gn );
	void setExtrapolationDist( FLOAT64 dist );
	virtual void buildSurface( std::vector<vec2d> &vertices, std::vector<vec2d> &normals, std::vector<vec2i> &faces, const levelset2 *fluid, const levelset2 *solid, bool enclose );
	virtual void smoothMesh( std::vector<vec2d> &vertices, std::vector<vec2d> &normals, std::vector<vec2i> &faces, uint iterations );
	virtual FLOAT64 evalLevelset(vec2d p) const;
	virtual vec2d evalGradient(vec2d p) const;
	virtual void drawGL( const levelset2 *solid, uint kind ) const;
protected:
	std::vector<std::vector<uint> > v2v;
	array2<std::vector<uint> > hash;
	array2<std::vector<uint> > vhash;
	array2<FLOAT64> nodalLevelset;
	array2<vec2d> referencePos;
	uint gn;
	FLOAT64 extrapolate_dist;
	FLOAT64 dx;
	void buildLevelsetInfo( std::vector<vec2d> &vertices, std::vector<vec2d> &normals, std::vector<vec2i> &faces );
	FLOAT64 findClosestSurfPosition( std::vector<vec2d> &vertices, std::vector<vec2i> &faces, int nodal_i, int nodal_j, vec2d &p ); // nodal pos
};

#endif