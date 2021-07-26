/*
 *	testcase.h
 *	
 *	Created by Ryoichi Ando on 11/9/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include <math.h>
#include "levelset2.h"

#ifndef _TESTCASE_H
#define _TESTCASE_H

// Feel free to add your test case here...
namespace testcase {
	static FLOAT64 s = 0.9;
	
	//////////////////////////// Dambreak test case (default) ////////////////////
	
	class dambreakSolid : public levelset2 {
	public:
		virtual FLOAT64 evalLevelset(vec2d p) const {
			FLOAT64 sgndist = 1.0;
			// Wall levelset
			sgndist = fmin(sgndist,-levelset2::box(p,vec2d(s*dx,s*dx),vec2d(1.-s*dx,1.-s*dx)));
			// Sphere Levelset
			sgndist = fmin(sgndist,levelset2::circle(p,vec2d(1.0,0.0),0.4));
			return sgndist;
		}
	};
	
	class dambreakFluid : public levelset2 {
	public:
		virtual FLOAT64 evalLevelset(vec2d p) const {
			return levelset2::box(p,vec2d(-1.0,-1.0),vec2d(0.5,0.75));
		}
	};
	
	/////////////////////////// Waterdrop test case ///////////////////////////
	
	class waterdropSolid : public levelset2 {
	public:
		virtual FLOAT64 evalLevelset(vec2d p) const {
			return -levelset2::box(p,vec2d(s*dx,s*dx),vec2d(1.-s*dx,1.-s*dx));
		}
	};
	
	class waterdropFluid : public levelset2 {
	public:
		virtual FLOAT64 evalLevelset(vec2d p) const {
			return fmin(levelset2::circle(p,vec2d(0.5,0.75),0.15), levelset2::box(p,vec2d(-1.0,-1.0),vec2d(2.0,0.5)) );
		}
	};
	
	/////////////////////////// Ocean wave test case ///////////////////////////
	
	// Don't forget add test cases in this list...
	static levelset2 *scenes[][2] = {
		{ new dambreakFluid(), new dambreakSolid() },
		{ new waterdropFluid(), new waterdropSolid() },
	};
	
	static const uint NUMSCENE = sizeof(scenes)/(2*sizeof(levelset2 *));
}

#endif