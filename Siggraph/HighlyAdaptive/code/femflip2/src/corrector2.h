/*
 *	corrector2.h
 *	
 *	Created by Ryoichi Ando on 11/6/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "macros.h"
#include "levelset2.h"
#include <vector>

#ifndef _CORRECTOR2_H
#define _CORRECTOR2_H

class particle2;
class levelset2;
class sorter2;
class surf2;
class corrector2 {
public:
	corrector2();
	void correct( sorter2 &sorter, int num, const levelset2 *level, std::vector<particle2 *> &particles, FLOAT64 stiffness, FLOAT64 sample_visc,
				  FLOAT64 dt, uint gsize, FLOAT64 dpx, const levelset2* fluid, const levelset2* solid, uint numSample );
};

#endif