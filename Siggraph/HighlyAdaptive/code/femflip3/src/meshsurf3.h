/*
 *	meshsurf3.h
 *
 *	Created by Ryoichi Ando on 8/21/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "macros.h"
#include "mesher3.h"
#include "levelset3.h"
#include <vector>

#ifndef _MESHSURF3_H
#define _MESHSURF3_H

class meshsurf3 : public levelset3 {
public:
	meshsurf3();
	meshsurf3( const meshsurf3 &meshsurf );
	void setReference(const mesher3 *g);
	void setExtrapolationDist( FLOAT64 dist );
	virtual void buildSurface( std::vector<vec3d> &vertices, std::vector<vec3d> &normals, std::vector<vec3i> &faces, const levelset3 *fluid, const levelset3 *solid, bool enclose=false, FLOAT64 dpx=0.0, uint iteration=1, bool doFit=false);
	virtual void drawGL( const levelset3 *solid, uint kind ) const;
	virtual FLOAT64 evalLevelset(vec3d p) const;
	virtual vec3d evalGradient(vec3d p) const;
	virtual void smoothMesh( std::vector<vec3d> &vertices, std::vector<vec3d> &normals, std::vector<vec3i> &faces, const levelset3 *solid, FLOAT64 dpx, uint iterations );
	static void fillHoles( std::vector<FLOAT64> &values, const levelset3 *solid, const std::vector<vec3d> &nodes, const std::vector<std::vector<uint> > &elements, const std::vector<std::vector<uint> > &node2node, FLOAT64 dx );
	const mesher3 *g;
protected:
	FLOAT64 extrapolate_dist;
	std::vector<FLOAT64> values;
	std::vector<std::vector<vec3d> > surfacePos;
	//
	bool isOpEdge( uint idx0, uint idx1 );
private:
	vec3d getElementGradient(uint elm) const;
	void computeNormals( std::vector<vec3d> &vertices, std::vector<vec3d> &normals, std::vector<vec3i> &faces, const std::vector<int> &cutIndices );
	void computeCurvature( const levelset3 *solid );
	void removeNan( std::vector<vec3d> &vertices, std::vector<vec3d> &normals, std::vector<vec3i> &faces, const std::vector<std::vector<uint> > &connections, const std::vector<FLOAT64> &areas );
	void flipFacet( std::vector<vec3d> &vertices, std::vector<vec3d> &normals, std::vector<vec3i> &faces, FLOAT64 min_dx );
};

#endif