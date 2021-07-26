/*
 *	object3.h
 *	
 *	Created by Ryoichi Ando on 1/8/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "opengl.h"
#include "macros.h"
#include "levelset3.h"
#include "vec3.h"

#ifndef _OBJECT3_H
#define _OBJECT3_H

static FLOAT64 s = 0.9;

class object3 : public levelset3 {
public:
	// type
	enum { FLUID, SOLID };
	object3( int type ) : type(type) {}
	virtual ~object3() {}
	virtual void draw() const {
		if( type == SOLID ) {
			std::vector<vec3d> vertices;
			std::vector<vec3i> faces;
			getPolygon(vertices,faces);
			glBegin(GL_LINES);
			glColor4f(1.0,1.0,1.0,0.4);
			for( uint n=0; n<faces.size(); n++ ) {
				for( uint m=0; m<3; m++ ) {
					uint idx0 = faces[n][m];
					uint idx1 = faces[n][(m+1)%3];
					if( idx0 < idx1 ) {
						glVertex3dv(vertices[idx0].v);
						glVertex3dv(vertices[idx1].v);
					}
				}
			}
			glEnd();
		}
	}
	virtual bool isVisible() const { return true; }
	virtual void getPolygon( std::vector<vec3d> &vertices, std::vector<vec3i> &faces ) const = 0;
	virtual bool intersect( vec3d p0, vec3d p1 ) const { return false; }
	int type;
};

class containerObject : public object3 {
public:
	containerObject() : object3(SOLID) {}
	virtual FLOAT64 evalLevelset(vec3d p) const {
		return -levelset3::box(p,vec3d(s*dx,s*dx,s*dx),vec3d(1.-s*dx,1.-s*dx,1.-s*dx));
	}
	virtual void getPolygon( std::vector<vec3d> &vertices, std::vector<vec3i> &faces ) const {
	}
};

class rectangleContainerObject : public object3 {
public:
	rectangleContainerObject(vec3d p0, vec3d p1) : object3(SOLID), p0(p0), p1(p1) {}
	virtual FLOAT64 evalLevelset(vec3d p) const {
		vec3d p0_f = p0;
		vec3d p1_f = p1;
		for( uint dim=0; dim<DIM; dim++ ) {
			p0_f[dim] = fmin(1.-s*dx,fmax(s*dx,p0_f[dim]));
			p1_f[dim] = fmin(1.-s*dx,fmax(s*dx,p1_f[dim]));
		}
		return -levelset3::box(p,p0_f,p1_f);
	}
	virtual void getPolygon( std::vector<vec3d> &vertices, std::vector<vec3i> &faces ) const {
	}
	vec3d p0;
	vec3d p1;
};

class boxObject : public object3 {
public:
	boxObject( int type, vec3d p0, vec3d p1 ) : object3(type), p0(p0), p1(p1) {
		vertices.resize(8);
		vertices[0] = vec3d(p0[0],p0[1],p0[2]);
		vertices[1] = vec3d(p0[0],p1[1],p0[2]);
		vertices[2] = vec3d(p1[0],p1[1],p0[2]);
		vertices[3] = vec3d(p1[0],p0[1],p0[2]);
		vertices[4] = vec3d(p0[0],p0[1],p1[2]);
		vertices[5] = vec3d(p0[0],p1[1],p1[2]);
		vertices[6] = vec3d(p1[0],p1[1],p1[2]);
		vertices[7] = vec3d(p1[0],p0[1],p1[2]);		
		faces.resize(12);
		faces[0] = vec3i(0,1,2);
		faces[1] = vec3i(2,3,0);
		faces[2] = vec3i(6,5,4);
		faces[3] = vec3i(4,7,6);
		faces[4] = vec3i(5,1,0);
		faces[5] = vec3i(0,4,5);
		faces[6] = vec3i(3,2,6);
		faces[7] = vec3i(6,7,3);
		faces[8] = vec3i(6,2,1);
		faces[9] = vec3i(1,5,6);
		faces[10] = vec3i(4,0,3);
		faces[11] = vec3i(3,7,4);
	}
	virtual FLOAT64 evalLevelset(vec3d p) const {
		return levelset3::box(p,p0,p1);
	}
	virtual void getPolygon( std::vector<vec3d> &vertices, std::vector<vec3i> &faces ) const {
		vertices = this->vertices;
		faces = this->faces;
	}
protected:
	std::vector<vec3d> vertices;
	std::vector<vec3i> faces;
	vec3d p0, p1;
};

class sphereObject : public object3 {
public:
	sphereObject( int type, vec3d center, FLOAT64 radius ) : object3(type), cp(center), r(radius) {}
	virtual FLOAT64 evalLevelset(vec3d p) const {
		return (cp-p).len()-r;
	}
	virtual void getPolygon( std::vector<vec3d> &vertices, std::vector<vec3i> &faces ) const {
	}
protected:
	vec3d cp;
	FLOAT64 r;
};

class objLevelset3 : public levelset3 {
public:
	objLevelset3( int type, std::vector<object3 *> objects ) : objects(objects), type(type){
	}
	virtual ~objLevelset3() {
		for( uint n=0; n<objects.size(); n++ ) {
			delete objects[n];
		}
	}
	virtual void resize( uint gsize ) {
		levelset3::resize(gsize);
		for( uint n=0; n<objects.size(); n++ ) {
			objects[n]->resize(gsize);
		}
	}
	virtual FLOAT64 evalCurvature(vec3d p) const {
		const object3 *obj = getDominated(p);
		if( obj ) {
			return obj->evalCurvature(p);
		} else {
			return levelset3::evalCurvature(p);
		}
	}
	virtual vec3d evalGradient(vec3d p) const {
		const object3 *obj = getDominated(p);
		if( obj ) {
			return obj->evalGradient(p);
		} else {
			return levelset3::evalGradient(p);
		}
	}
	virtual FLOAT64 evalLevelset(vec3d p) const {
		FLOAT64 sgndist = 1e9;
		for( uint n=0; n<objects.size(); n++ ) {
			if( objects[n]->type == type ) {
				sgndist = fmin(sgndist,objects[n]->evalLevelset(p));
			}
		}
		return sgndist;
	}
	virtual void draw() const {
		for( uint n=0; n<objects.size(); n++ ) {
			if( objects[n]->type == object3::SOLID ) objects[n]->draw();
		}
	}
	int type;
	std::vector<object3 *> objects;
private:
	const object3 * getDominated(vec3d p) const {
		FLOAT64 min_dist[2] = { 1e9, 1e9 };
		object3 *obj = NULL;
		for( uint n=0; n<objects.size(); n++ ) {
			if( objects[n]->type == type ) {
				FLOAT64 d = objects[n]->evalLevelset(p);
				if( d < min_dist[0] ) {
					min_dist[1] = min_dist[0];
					min_dist[0] = d;
					obj = objects[n];
				}
			}
		}
		if( obj && fabs(min_dist[0]-min_dist[1]) > 2.0*dx ) {
			return obj;
		} else {
			return NULL;
		}
	}
};

#endif
