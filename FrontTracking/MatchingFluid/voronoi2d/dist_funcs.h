#ifndef DISTFUNCS_H
#define DISTFUNCS_H

#include "vec.h"

float circle(const Vec2f& point, const Vec2f& centre, float radius);
float capsule(const Vec2f& point, const Vec2f& start, const Vec2f& end, float radius);
float box(const Vec2f& point, const Vec2f& centre, float width, float height);
float rotated_box(const Vec2f& point, const Vec2f& centre, float width, float height, float angle);

void make_circle(Vec2f centre, float radius, int segments, std::vector<Vec2f>& verts, std::vector<Vec2ui>& edges);
void make_box(Vec2f bottom_left, float height, float width, std::vector<Vec2f>& verts, std::vector<Vec2ui>& edges);

Vec2f box_gradient( const Vec2f& point, const Vec2f& centre, float width, float height );

#endif