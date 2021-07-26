/*
 *	adaptive3.h
 *	
 *	Created by Ryoichi Ando on 1/10/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "macros.h"
#include "sorter3.h"
#include "array3.h"
#include <vector>

#ifndef _ADAPTIVE3_H
#define _ADAPTIVE3_H

class particle3;
class levelset3;
class adaptive3 {
public:
	void merge( sorter3& sorter, const levelset3 *level, const levelset3 *surface, std::vector<particle3 *> &particles, FLOAT64 dpx, FLOAT64 dx );
    void split( sorter3& sorter, const levelset3 *level, const levelset3 *surface, std::vector<particle3 *> &particles, FLOAT64 dpx, FLOAT64 dx, bool splitAll=false );
private:
    particle3 *findClosestParticle( const sorter3& sorter, const particle3 &p, FLOAT64 r );
};

#endif