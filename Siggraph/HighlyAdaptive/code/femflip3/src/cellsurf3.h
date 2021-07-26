/*
 *	cellsurf3.h
 *	
 *	Created by Ryoichi Ando on 1/9/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "levelset3.h"
#include "macros.h"
#include "vec3.h"
#include "array3.h"
#include <vector>

#ifndef _CELLSURF3_H
#define _CELLSURF3_H

class cellsurf3 : public levelset3 {
public:
	cellsurf3();
	cellsurf3( const cellsurf3& cellsurf );
	void setResolution( uint gn );
	void setExtrapolationDist( FLOAT64 dist );
	virtual void buildSurface( std::vector<vec3d> &vertices, std::vector<vec3d> &normals, std::vector<vec3i> &faces, const levelset3 *fluid, const levelset3 *solid, bool enclose );
	virtual FLOAT64 evalLevelset(vec3d p) const;
	virtual vec3d evalGradient(vec3d p) const;
	virtual void drawGL( const levelset3 *solid, uint kind ) const;
protected:
	array3<std::vector<uint> > hash;
	array3<std::vector<uint> > vhash;
	array3<FLOAT64> nodalLevelset;
	array3<vec3d>	referencePos;
	uint gn;
	FLOAT64 extrapolate_dist;
	FLOAT64 dx;
	void buildLevelsetInfo( std::vector<vec3d> &vertices, std::vector<vec3d> &normals, std::vector<vec3i> &faces );
	FLOAT64 findClosestSurfPosition( std::vector<vec3d> &vertices, std::vector<vec3i> &faces, int nodal_i, int nodal_j, int nodal_k, vec3d &p );
	void removeNan( std::vector<vec3d> &vertices, std::vector<vec3d> &normals, std::vector<vec3i> &faces, const std::vector<std::vector<uint> > &connections, const std::vector<FLOAT64> &areas );
	void smoothMesh( std::vector<vec3d> &vertices, std::vector<vec3d> &normals, std::vector<vec3i> &faces, const levelset3 *solid, uint iterations );
};

#endif