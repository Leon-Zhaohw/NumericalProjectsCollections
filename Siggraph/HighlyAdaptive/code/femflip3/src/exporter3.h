/*
 *	exporter3.h
 *	
 *	Created by Ryoichi Ando on 1/9/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "macros.h"
#include "object3.h"
#include "surf3.h"
#include "flip3.h"
#include <vector>
#include <string>

#ifndef _exporter3_H_
#define _exporter3_H_

class particle3;
class bcc3;
class octree3;
class camera3;
class object3;
class exporter3 {
public:
	void writeMesh(uint frame, const surf3 &surf, const std::vector<particle3 *> &particles, FLOAT64 dpx,
				   const mesher3 *mesher, const octree3 *octree, const camera3 *camera, const std::vector<object3 *> *objects,
				   bool writeParticleData=false );
protected:
	void write_mitsuba(uint frame, const surf3 &surf, const std::vector<particle3 *> &particles, FLOAT64 dpx, bool writeParticleData,
					   const camera3 *camera, const std::vector<object3 *> *objects, const mesher3 *tet );
};

#endif