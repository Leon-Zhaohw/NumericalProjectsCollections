/*
**	facecutter.h
**
**	Created by Anonymouse Authors
**
**	Permission is hereby granted, free of charge, to any person obtaining a copy of
**	this software and associated documentation files (the "Software"), to deal in
**	the Software without restriction, including without limitation the rights to use,
**	copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the
**	Software, and to permit persons to whom the Software is furnished to do so,
**	subject to the following conditions:
**
**	The above copyright notice and this permission notice shall be included in all copies
**	or substantial portions of the Software.
**
**	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
**	INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
**	PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
**	HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
**	CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
**	OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
//
#ifndef FACECUTTER_H
#define FACECUTTER_H
//
#include <functional>
#include <algorithm>
#include <limits>
#include <cmath>
//
struct point2 {
	float x[2];
	point2() { x[0]=x[1]=0.0; }
	point2 ( float x, float y ) { this->x[0]=x; this->x[1]=y; }
	point2 operator*(float s) const { return point2(s*x[0],s*x[1]); }
	point2 operator+(const point2 &a) const { return point2(a.x[0]+x[0],a.x[1]+x[1]); }
};
static inline point2 operator*(float s, const point2 &p) { return p*s; }
static inline char sgn( float val ) { return val < 0.0 ? -1 : 1; }
//
struct vertex2 {
	point2 xi;	// (xi,eta) coordinate
	point2 x;	// (x,z) coordinate
	vertex2 () {}
	vertex2 ( const point2 &xi ) { this->xi=xi; }
	vertex2 ( const point2 &xi, const point2 &x ) { this->xi=xi; this->x=x; }
};
//
static void subdivide( 
	const vertex2 *polygon,
	unsigned char count,
	std::function<void(const vertex2 *polygon, unsigned char count)> output_func,
	unsigned char dim=0
) {
	const unsigned const_max_tmp_points (32);
	float min_x (std::numeric_limits<float>::max());
	float max_x (std::numeric_limits<float>::lowest());
	for( int n=0; n<count; ++n ) {
		const float c (polygon[n].x.x[dim]);
		min_x = std::min(min_x,c);
		max_x = std::max(max_x,c);
	}
	//
	int start (std::ceil(min_x));
	int end (std::floor(max_x));
	if( min_x - start == 0.0 ) start += 1;
	if( max_x - end == 0.0 ) end -= 1;
	//
	bool cut_produced (false);
	if( start <= end ) {
		//
		const float iso_value (start);
		vertex2 vertex_array[2][const_max_tmp_points];
		int vertex_array_head[2] = { 0, 0 };
		int vertex_array_slot (0);
		//
		vertex_array[0][vertex_array_head[0]++] = polygon[0];
		for( int n=0; n<count; ++n ) {
			//
			const int m = (n+1) % count;
			const float v0 = polygon[n].x.x[dim] - iso_value;
			const float v1 = polygon[m].x.x[dim] - iso_value;
			//
			if( (v0 or v1) and sgn(v0) * sgn(v1) < 0 ) {
				const float inv_len(1.0/(v1-v0));
				const float t0 (v1 * inv_len), t1 (-v0 * inv_len);
				cut_produced = true;
				//
				vertex2 vertex_add(t0 * polygon[n].xi + t1 * polygon[m].xi);
				vertex_add.x.x[dim] = iso_value;
				vertex_add.x.x[1-dim] = t0 * polygon[n].x.x[1-dim] + t1 * polygon[m].x.x[1-dim];
				//
				if( vertex_array_slot == 0 ) {
					vertex_array[0][vertex_array_head[0]++] = vertex_add;
					vertex_array_slot = 1;
					vertex_array_head[1] = 0;
					if( t1 < 1.0 ) vertex_array[vertex_array_slot][vertex_array_head[vertex_array_slot]++] = vertex_add;
				} else {
					vertex_array[1][vertex_array_head[1]++] = vertex_add;
					subdivide(vertex_array[1],vertex_array_head[1],output_func,dim);
					vertex_array_slot = 0;
					if( t1 < 1.0 ) vertex_array[vertex_array_slot][vertex_array_head[vertex_array_slot]++] = vertex_add;
				}
			}
			//
			if( n < count-1 ) vertex_array[vertex_array_slot][vertex_array_head[vertex_array_slot]++] = polygon[m];
			//
			assert( vertex_array_head[0] < const_max_tmp_points );
			assert( vertex_array_head[1] < const_max_tmp_points );
		}
		if( cut_produced ) {
			subdivide(vertex_array[0],vertex_array_head[0],output_func,dim);
		}
	}
	if( ! cut_produced ) {
		if( dim == 0 ) {
			subdivide(polygon,count,output_func,1);
		} else {
			output_func(polygon,count);
		}
	}
}
//
#endif
