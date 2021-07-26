/*
 *	corrector3.h
 *	
 *	Created by Ryoichi Ando on 1/10/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "macros.h"
#include "levelset3.h"
#include <vector>

#ifndef _CORRECTOR3_H
#define _CORRECTOR3_H

class particle3;
class levelset3;
class sorter3;
class corrector3 {
public:
	corrector3();
	void correct( sorter3 &sorter, int num, const levelset3 *level, std::vector<particle3 *> &particles, FLOAT64 stiffness, FLOAT64 sample_visc,
				  FLOAT64 dt, uint gsize, FLOAT64 dpx, const levelset3* fluid, const levelset3* solid, uint numSamplem );
};

#endif