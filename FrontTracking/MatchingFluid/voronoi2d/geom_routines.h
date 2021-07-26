#ifndef GEOM_ROUTINES_H
#define GEOM_ROUTINES_H

#include "vec.h"

bool point_in_box(Vec2f& pt, Vec2f& boxMin, float width, float height);

bool point_in_tri(const Vec2f& point, const Vec2f& p0, const Vec2f& p1, const Vec2f& p2);

bool seg_hits_box(const Vec2f&a, const Vec2f&b, const Vec2f&x00, float width, float height);

bool tri_intersects_square(Vec2f& a, Vec2f& b, Vec2f& c, Vec2f& x00, float width, float height);

void get_bound_box(const Vec2f& a, const Vec2f& b, const Vec2f& c,  Vec2f& minPt, Vec2f& maxPt);
#endif

