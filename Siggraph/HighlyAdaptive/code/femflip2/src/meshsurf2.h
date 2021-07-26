/*
 *	meshsurf2.h
 *
 *	Created by Ryoichi Ando on 8/13/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "macros.h"
#include "mesher2.h"
#include "levelset2.h"
#include <vector>

#ifndef _MESHSURF2_H
#define _MESHSURF2_H

class meshsurf2 : public levelset2 {
public:
	meshsurf2();
	meshsurf2( const meshsurf2 &meshsurf );
	void setReference(const mesher2 *g);
	void setExtrapolationDist( FLOAT64 dist );
	virtual void buildSurface( std::vector<vec2d> &vertices, std::vector<vec2d> &normals, std::vector<vec2i> &faces, const levelset2 *fluid, const levelset2 *solid, bool enclose=false, FLOAT64 dpx=0.0, uint iteration=1, bool doFit=true );
	virtual void smoothMesh( std::vector<vec2d> &vertices, std::vector<vec2d> &normals, std::vector<vec2i> &faces, const levelset2 *solid, FLOAT64 dpx, uint iterations );
	virtual void drawGL( const levelset2 *solid, uint kind ) const;
	virtual FLOAT64 evalLevelset(vec2d p) const;
	virtual vec2d evalGradient(vec2d p) const;

	const mesher2 *g;
	std::vector<FLOAT64> values;
	std::vector<std::vector<vec2d> > surfacePos;
	static void fillHoles( std::vector<FLOAT64> &values, const levelset2 *solid, const std::vector<vec2d> &nodes, const std::vector<std::vector<uint> > &elements, const std::vector<std::vector<uint> > &node2node, FLOAT64 dx );
protected:
	FLOAT64 extrapolate_dist;
private:
	vec2d getElementGradient(uint elm) const;
};

#endif