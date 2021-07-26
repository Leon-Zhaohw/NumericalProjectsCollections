/*
 *	util2.h
 *	
 *	Created by Ryoichi Ando on 11/4/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include "macros.h"
#include "vec2.h"
#include "array2.h"
#include "svd2.h"
#include "email.h"

#if defined(WIN32)
#include <windows.h>
#else
#include <sys/time.h>
#endif

#ifndef _UTILITY_H_
#define _UTILITY_H_

extern int debugLevel;
static void dump( const char *format, ...) {
	if( debugLevel ) {
		va_list args;
		if( 0 ) {
			FILE *console = fopen( "console.out", "a" );
			va_start(args, format);
			vfprintf(console, format, args);
			va_end(args);
			fclose(console);
		}
		va_start(args, format);
		vprintf(format, args);
		va_end(args);
	}
}

static void run( const char *format, ...) {
	char command[512];
	va_list args;
	va_start(args, format);
	vsprintf(command,format,args);
	va_end(args);
	system(command);
}

static const char *format_str( const char *format, ...) {
	static char temporary_string[256][512];
	static int slot = 0;
	if( slot >= 256 ) slot = 0;
	va_list args;
	va_start(args, format);
	vsprintf(temporary_string[slot],format,args);
	va_end(args);
	return temporary_string[slot++];
}

static bool is_nan(const double v) {
    return (v != v);
}

inline FLOAT64 nrand() {
	return 2.0 * ((rand() % 100001) / 100000.0 - 0.5);
}

namespace util2 {		
	// 2-Dimentional liner interpolation of quantity q qt position p
	template <class T> static T interp( vec2d p, const array2<T> &q ) {
		uint w, h;
		w = q.size().w;
		h = q.size().h;
		FLOAT64 x = fmax(0.0,fmin(w,p[0]));
		FLOAT64 y = fmax(0.0,fmin(h,p[1]));
		int i = imin(x,w-2);
		int j = imin(y,h-2);
		return ((i+1-x)*q[i][j]+(x-i)*q[i+1][j])*(j+1-y) + ((i+1-x)*q[i][j+1]+(x-i)*q[i+1][j+1])*(y-j);
	}
	
	static FLOAT64 square(FLOAT64 a) {
		return a*a;
	}
	
	template <class T> static void blur( array2<T> &q, uint r ) {
		array2<T> old = q;
		FOR_EACH_SIZE(q.size()) {
			T sum = T();
			int n=0;
			FOR_NEIGHBORS(i,j,r,r,q.size().w,q.size().h) {
				sum += old[ni][nj];
				n++;
			} END_FOR
			q[i][j] = sum/(FLOAT64)n;
		} END_FOR
	}
	
	static unsigned long getMicroseconds() {
#if defined(_WIN32)	// I haven't tested this (Windows version code) works...hopefully it does
		LARGE_INTEGER nFreq, Time;
		QueryPerformanceFrequency(&nFreq);
		QueryPerformanceCounter(&Time);
		return (double)Time.QuadPart / nFreq.QuadPart * 1000000;
#else
		struct timeval tv;
		gettimeofday(&tv, NULL);
		return tv.tv_sec*1000000 + tv.tv_usec;
#endif
	}
	
	static FLOAT64 getMilliseconds() {
		return getMicroseconds()/1000.0;
	}
	
	static FLOAT64 getSeconds() {
		return getMilliseconds()/1000.0;
	}
	
	template <class T> static T arrayMax( const array2<T> &q, bool abs=false ) {
		T maxv = -99999.0;
		FOR_EACH_SIZE(q.size()) {
			T v = abs ? fabs(q[i][j]) : q[i][j];
			if( maxv < v ) maxv = v;
		} END_FOR
		return maxv;
	}
	
	template <class T> static T arrayMin( const array2<T> &q, bool abs=false ) {
		T minv = 999999.0;
		FOR_EACH_SIZE(q.size()) {
			T v = abs ? fabs(q[i][j]) : q[i][j];
			if( minv > v ) minv = v;
		} END_FOR
		return minv;
	}
	
	static vec2d stretchPosition( const svd2 &svd, vec2d p, FLOAT64 min, FLOAT64 max, bool inverse=true ) {
		// Compute Dot Product
		FLOAT64 dot[3];
		for( int k=0; k<DIM; k++ ) {
			dot[k] = p[0]*svd.R[0][k] + p[1]*svd.R[1][k];
		}
		// Scale By Eigen Value
		for( int k=0; k<DIM; k++ ) {
			FLOAT64 r = svd.eig[k];
			r = fmax(min,fmin(max,r));
			if( inverse ) dot[k] /= r;
			else dot[k] *= r;
		}
		// Compute Final Position
		for( int k=0; k<DIM; k++ ) {
			p[k] = dot[0]*(svd.R[k][0]) + dot[1]*svd.R[k][1];
		}
		return p;
	}
	
	static std::vector<vec2d> marchPoints( const std::vector<vec2d> &nodes, const std::vector<FLOAT64> &levelsets, bool fill=true ) {
		std::vector<vec2d> points;
		for( uint n=0; n<nodes.size(); n++ ) {
			const vec2d &p0 = nodes[n];
			const vec2d &p1 = nodes[(n+1)%nodes.size()];
			FLOAT64 ls0 = levelsets[n];
			FLOAT64 ls1 = levelsets[(n+1)%nodes.size()];
			if( fill && ls0 < 0.0 ) {
				points.push_back(p0);
			}
			if( copysign(1.0,ls0)*copysign(1.0,ls1) < 0.0 ) {
				FLOAT64 a = ls0/(ls0-ls1);
				points.push_back((1.0-a)*p0+a*p1);
			}
		}
		return points;
	}

	static FLOAT64 computeVolume( const std::vector<vec2d> points ) {
		FLOAT64 v = 0.0;
		for( uint i=0; i<points.size(); i++ ) {
			FLOAT64 x0 = points[i][0];
			FLOAT64 y0 = points[i][1];
			FLOAT64 x1 = points[(i+1)%points.size()][0];
			FLOAT64 y1 = points[(i+1)%points.size()][1];
			v += x0*y1-y0*x1;
		}
		return 0.5*fabs(v);
	}
	
	static FLOAT64 distance( const vec2d &p0, const vec2d &p1, vec2d &p ) {
		vec2d a = p1-p0;
		vec2d b = p-p0;
		FLOAT64 alen = a.len();
		if( alen ) {
			FLOAT64 dot = a/alen * b;
			dot = fmin(alen,fmax(0.0,dot));
			vec2d o = p;
			p = p0 + (a/alen) * dot;
			return (p-o).len();
		} else {
			p = p0;
			return a.len();
		}
	}
	
	static bool centerCircumcircle(const vec2d& p0, const vec2d& p1, const vec2d& p2, vec2d& center ){
		FLOAT64 p0x = p0[0]; FLOAT64 p1x = p1[0]; FLOAT64 p2x = p2[0];
		FLOAT64 p0y = p0[1]; FLOAT64 p1y = p1[1]; FLOAT64 p2y = p2[1];
		FLOAT64 k1 = (p1x*p1x+p1y*p1y-p0x*p0x-p0y*p0y)/2.0;
		FLOAT64 k2 = (p2x*p2x+p2y*p2y-p1x*p1x-p1y*p1y)/2.0;
		FLOAT64 u[2] = {p1[0]-p0[0],p1[1]-p0[1]};
		FLOAT64 v[2] = {p2[0]-p1[0],p2[1]-p1[1]};
		FLOAT64 det = u[0]*v[1] - u[1]*v[0];
		if( ! det ) return false;
		center[0] = (1.0/det)*(k1*v[1] - k2*u[1] );
		center[1] = (1.0/det)*(k2*u[0] - k1*v[0] );
		return true;
	}
	
	static FLOAT64 detDelaunay( vec2d p0, vec2d p1, vec2d p2, vec2d p3 ) {
		vec2d center;
		if( ! centerCircumcircle(p0, p1, p2, center )) return -1.0;
		return (center-p3).len2() - (center-p0).len2();
	}
}

#endif