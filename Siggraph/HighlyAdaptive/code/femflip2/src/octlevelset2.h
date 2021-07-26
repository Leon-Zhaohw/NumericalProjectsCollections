/*
 *	octlevelset2.h
 *
 *	Created by Ryoichi Ando on 9/5/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "levelset2.h"
#include "octgrid2.h"

#ifndef _OCTLEVELSET2_H
#define _OCTLEVELSET2_H

class octlevelset2 : public levelset2 {
public:
	octlevelset2() {}
	virtual void setOctree( const octree2 &octree ) {
		this->octree = octree;
		octgrid.setOctree(this->octree);
	}
	virtual void setNodalLevelset( const std::vector<FLOAT64> &q ) {
		octgrid.setValue(q);
	}
	virtual FLOAT64 evalLevelset(vec2d p) const {
		return octgrid.getValue(p);
	}
	octree2 octree;
	octgrid2<FLOAT64> octgrid;
};

#endif