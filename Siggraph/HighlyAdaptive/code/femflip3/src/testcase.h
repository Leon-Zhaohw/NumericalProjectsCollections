/*
 *	testcase.h
 *	
 *	Created by Ryoichi Ando on 2013/01/05
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "object3.h"
#include "flycamera3.h"
#include <vector>

#ifndef _TESTCASE_H
#define _TESTCASE_H

// Feel free to add your test case here...
namespace testcase {
	// Test case
	static std::vector<object3 *> buildTestcase() {
		std::vector<object3 *> objects;
		objects.push_back(new containerObject());
		objects.push_back(new boxObject(object3::FLUID,vec3d(-1.0,-1.0,-1.0),vec3d(2.0,0.375,2.0)));
		objects.push_back(new boxObject(object3::FLUID,vec3d(0.2,-1.0,0.2),vec3d(0.4,0.55,0.8)));
		objects.push_back(new sphereObject(object3::FLUID,vec3d(0.8,0.4,0.6),0.1));
		objects.push_back(new sphereObject(object3::FLUID,vec3d(0.6,0.4,0.3),0.1));
		return objects;
	}	
	// Change this line to switch testcase
	static flycamera3 *flycamera = new static_camera3();
	const static FLOAT64 gravity = 9.8;
}

#endif