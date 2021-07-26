/*
 *	bccsurf2.h
 *
 *	Created by Ryoichi Ando on 7/17/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "meshsurf2.h"
#include "bccmesher2.h"

#ifndef _BCCSURF2_H
#define _BCCSURF2_H

class bccsurf2 : public levelset2 {
public:
	bccsurf2();
	bccsurf2( const bccsurf2 &bccsurf );
	void operator=( const bccsurf2 &bccsurf );
	void setExtrapolationDist( FLOAT64 dist );
	void setBCC( const bccmesher2 &g );
	virtual void buildSurface( std::vector<vec2d> &vertices, std::vector<vec2d> &normals, std::vector<vec2i> &faces, const levelset2 *fluid, const levelset2 *solid, bool enclose );
	virtual void drawGL( const levelset2 *solid, uint kind ) const;
	virtual FLOAT64 evalLevelset(vec2d p) const;
	virtual vec2d evalGradient(vec2d p) const;
	virtual FLOAT64 evalCurvature(vec2d p) const;

	meshsurf2 meshsurf;
	bccmesher2 g;
};

#endif