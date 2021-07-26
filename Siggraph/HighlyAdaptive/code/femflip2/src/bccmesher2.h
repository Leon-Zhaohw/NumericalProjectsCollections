/*
 *	bccmesher2.h
 *
 *	Created by Ryoichi Ando on 5/28/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "mesher2.h"
#include "octree2.h"
#include "bcc2.h"
#ifndef _BCCMESHER2_H
#define _BCCMESHER2_H

class bccmesher2 : public mesher2 {
public:
	virtual void generateElements( std::vector<vec2d> &nodes, std::vector<std::vector<uint> > &elements, const levelset2 *hint, uint gn );
	virtual int	 hitElements(const vec2d p) const;
	octree2 octree;
	bcc2 bcc;
};

#endif