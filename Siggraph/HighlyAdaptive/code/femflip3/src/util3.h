/*
 *	util3.h
 *	
 *	Created by Ryoichi Ando on 1/8/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include <stack>
#include <math.h>
#include <map>
#include <string>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include "macros.h"
#include "vec3.h"
#include "array3.h"
#include "svd3.h"
#include "email.h"
#include <sys/stat.h>

#if defined(WIN32)
#include <windows.h>
#else
#include <sys/time.h>
#endif

#if defined(__APPLE__) || defined(MACOSX) // For Mac OS X
#define FLUSH		1
#else
#define FLUSH		0
#endif

#ifndef UTIL3_H
#define UTIL3_H

extern int nestLevel;
extern int stepNumber;
extern char *root_path;
extern uint console_num;
extern std::map<std::string,FLOAT64> status_table;

static void inc_nest() {
	nestLevel ++;
	if( nestLevel > 5 ) {
		printf( "Too nested dump !\n" );
		exit(0);
	}
}
static void dec_nest() {
	nestLevel --;
	if( nestLevel < 0 ) {
		printf( "Nest level broken !\n" );
		exit(0);
	}
}

static const char * tstr(FLOAT64 msec) {
	static char dummy[256];
	if( msec < 1000.0 ) {
		// Milli seconds
		sprintf(dummy,"%.2f msec",msec);
	} else if( msec < 1000.0 * 60 * 3 ) {
		// Seconds
		sprintf(dummy,"%.3f sec",msec/1000.0);
	} else if( msec < 1000.0 * 60 * 60 * 3 ) {
		// Minutes
		sprintf(dummy,"%.3f minutes",msec/(60.0*1000.0));
	} else if( msec < 1000.0 * 60 * 60 * 60 * 3 ) {
		// Hours
		sprintf(dummy,"%.3f hours",msec/(60*60.0*1000.0));
	} else {
		// Days
		sprintf(dummy,"%.3f days",msec/(60*60.0*24*1000.0));
	} 
	return dummy;
}

static const char * nth(int num) {
	static char dummy[256];
	switch(num) {
		case 1:
			sprintf(dummy,"1st");
			break;
		case 2:
			sprintf(dummy,"2nd");
			break;
		case 3:
			sprintf(dummy,"3rd");
			break;
		default:
			sprintf(dummy,"%dth",num);
			break;
	}
	return dummy;
}

static int is_exist( const char *format, ... ) {
	char path[512];
	va_list args;
	va_start(args, format);
	vsprintf(path,format,args);
	va_end(args);
	
	struct stat buffer;
	return stat( path, &buffer ) == 0;
}

extern bool ended_with_return;

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

static void dump(const char *format, ...) {
#if DUMP
	uint len = strlen(format);
	va_list args;
	if( len > 2 ) {
		if( format[0]=='<' && format[1]=='<' && format[2]=='<' ) {
			dec_nest();
		}
	}
	const char *console_path;
	if( console_num ) console_path = format_str("%s/console_%d.out",root_path,console_num);
	else console_path = format_str("%s/console.out",root_path);
	
	FILE *console = fopen( console_path, "a" );
	if( ! console ) {
		if( ! is_exist(root_path)) run( "mkdir %s", root_path);
		console = fopen( console_path, "a" );
	}
	if( console ) {
		if( ended_with_return ) {
			for( uint i=0; i<nestLevel; i++ ) fprintf(console,"   ");
		}
		va_start(args, format);
		vfprintf(console, format, args);
		va_end(args);
		fclose(console);
		if( FLUSH ) fflush(console);
	}
	
	if( ended_with_return ) {
		for( uint i=0; i<nestLevel; i++ ) printf("   ");
	}
	va_start(args, format);
	vprintf(format, args);
	va_end(args);
	ended_with_return = format[len-1] == '\n';
	if( len > 2 ) {
		if( format[0]=='>' && format[1]=='>' && format[2]=='>' ) {
			inc_nest();
		}
	}
	if( FLUSH ) fflush(stdout);
#endif
}

static bool is_nan(const double v) {
    return (v != v);
}

inline FLOAT64 nrand() {
	return 2.0 * ((rand() % 100001) / 100000.0 - 0.5);
}

namespace util3 {
	// 3-Dimentional liner interpolation of quantity q qt position p
	template <class T> static T interp( vec3d p, const array3<T> &q ) {
		uint w, h, d;
		w = q.size().w;
		h = q.size().h;
		d = q.size().d;
		FLOAT64 x = fmax(0.0,fmin(w,p[0]));
		FLOAT64 y = fmax(0.0,fmin(h,p[1]));
		FLOAT64 z = fmax(0.0,fmin(d,p[2]));
		int i = imin(x,w-2);
		int j = imin(y,h-2);
		int k = imin(z,d-2);
		return	(k+1-z)*(((i+1-x)*q[i][j][k]+(x-i)*q[i+1][j][k])*(j+1-y)+((i+1-x)*q[i][j+1][k]+(x-i)*q[i+1][j+1][k])*(y-j))+
				(z-k)*(((i+1-x)*q[i][j][k+1]+(x-i)*q[i+1][j][k+1])*(j+1-y)+((i+1-x)*q[i][j+1][k+1]+(x-i)*q[i+1][j+1][k+1])*(y-j));
	}
	
	static FLOAT64 square(FLOAT64 a) {
		return a*a;
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
	
	template <class T> static T arrayMax( const array3<T> &q, bool abs=false ) {
		T maxv = -99999.0;
		FOR_EACH_SIZE(q.size()) {
			T v = abs ? fabs(q[i][j][k]) : q[i][j][k];
			if( maxv < v ) maxv = v;
		} END_FOR
		return maxv;
	}
	
	template <class T> static T arrayMin( const array3<T> &q, bool abs=false ) {
		T minv = 999999.0;
		FOR_EACH_SIZE(q.size()) {
			T v = abs ? fabs(q[i][j][k]) : q[i][j][k];
			if( minv > v ) minv = v;
		} END_FOR
		return minv;
	}
	
	static vec3d stretchPosition( const svd3 &svd, vec3d p, FLOAT64 min, FLOAT64 max, bool inverse=true ) {
		// Compute Dot Product
		FLOAT64 dot[3];
		for( int k=0; k<DIM; k++ ) {
			dot[k] = p[0]*svd.R[0][k] + p[1]*svd.R[1][k] + p[2]*svd.R[2][k];
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
			p[k] = dot[0]*(svd.R[k][0]) + dot[1]*svd.R[k][1] + dot[2]*svd.R[k][2];
		}
		return p;
	}
	
	static std::vector<vec3d> march2DPoints( const std::vector<vec3d> &nodes, const std::vector<FLOAT64> &levelsets, bool fill=true ) {
		std::vector<vec3d> points;
		for( uint n=0; n<nodes.size(); n++ ) {
			const vec3d &p0 = nodes[n];
			const vec3d &p1 = nodes[(n+1)%nodes.size()];
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
	
	static bool projectTriangle( vec3d *points, uint num=3 ) {
		// Compute normal vector of this triangle
		for( uint i=0; i<3; i++ ) for( uint j=i+1; j<3; j++ ) {
			if( ! (points[i]-points[j]).len2() ) return false;
		}
		
		vec3d n = ((points[2]-points[0])^(points[1]-points[0])).normal();
		vec3d e[2];
		e[0] = (points[1]-points[0]).normal();
		e[1] = n ^ e[0];
		
		// Take a dot product
		std::vector<vec3d> old_points(num);
		for( uint i=0; i<num; i++ ) old_points[i] = points[i];
		for( uint i=0; i<num; i++ ) {
			points[i][0] = e[0] * (old_points[i]-old_points[0]);
			points[i][1] = e[1] * (old_points[i]-old_points[0]);
			points[i][2] = 0.0;
		}
		return true;
	}
	
	static FLOAT64 compute2DArea( const std::vector<vec3d> points ) {
		FLOAT64 v = 0.0;
		for( uint i=0; i<points.size(); i++ ) {
			if( points[i][2] ) {
				printf( "compute2DArea can only handle 2D points !\n" );
				exit(0);
			}
			FLOAT64 x0 = points[i][0];
			FLOAT64 y0 = points[i][1];
			FLOAT64 x1 = points[(i+1)%points.size()][0];
			FLOAT64 y1 = points[(i+1)%points.size()][1];
			v += x0*y1-y0*x1;
		}
		return 0.5*fabs(v);
	}
	
	static void writeContainer( FILE *fp ) {
		fprintf(fp,
				"v 0 0 0\n"
				"v 1 0 0\n"
				"v 0 1 0\n"
				"v 1 1 0\n"
				"v 0 0 1\n"
				"v 1 0 1\n"
				"v 0 1 1\n"
				"v 1 1 1\n"
				"vn 0 0 -1\n"
				"vn -1 0 0\n"
				"vn 1 0 0\n"
				"vn 0 -1 0\n"
				"vn 0 1 0\n"
				"vn 0 0 1\n"
				"f 1//1 3//1 4//1 2//1\n"
				"f 1//2 5//2 7//2 3//2\n"
				"f 2//3 4//3 8//3 6//3\n"
				"f 1//4 2//4 6//4 5//4\n"
				"f 3//5 7//5 8//5 4//5\n"
				"f 5//6 6//6 8//6 7//6\n");
	}
}

extern std::stack<FLOAT64> global_timestack;
static void tick() {
	global_timestack.push(util3::getMilliseconds());
	if( global_timestack.size() > 10 ) {
		printf( "Too nested tick tock stack ! (%d)\n", (int)global_timestack.size() );
		exit(0);
	}
}

static void setTimestepNumber( int number ) {
	stepNumber = number;
}

static void setRootPath( const char *path ) {
	if( ! is_exist(path)) {
		run("mkdir %s",path);
	}
	strcpy(root_path,path);
}

static void setConsoleNumber( uint num ) {
	console_num = num;
}

static bool writeNumber( const char *name, FLOAT64 number ) {
	static bool firstTime = true;
	const char *record_path;
	if( console_num ) record_path = format_str("%s/record_%d",root_path,console_num);
	else record_path = format_str("%s/record",root_path);
	
	if( firstTime ) {
		if( ! is_exist(record_path)) run( "mkdir %s", record_path);
		run( "cp src/plot.sh %s/plot.sh", record_path );
		firstTime = false;
	}
	FILE *console = fopen( format_str("%s/%s.out", record_path, name), "a" );
	if( console ) {
		fprintf( console, "%d %f\n", stepNumber, number );
		fclose(console);
	}
	status_table[std::string(name)] = number;
	return false;
}

static FLOAT64 getNumber( const char *name ) {
	if( status_table.find(std::string(name)) != status_table.end()) return status_table[std::string(name)];
	return 0.0;
}

static FLOAT64 tock( const char *name=NULL ) {
	if( global_timestack.empty() ) {
		printf( "Tick tock stack broken !\n" );
		exit(0);
	}
	FLOAT64 elapsed = util3::getMilliseconds()-global_timestack.top();
	global_timestack.pop();
	if( name ) {
		writeNumber(name,elapsed);
	}
	return elapsed;
}

static const char * stock( const char *name=NULL ) {
	return tstr(tock(name));
}

static const char *memstr( uint bytes ) {
	static char buf[256];
	buf[0]=0;
	if( bytes < 1024 ) {
		sprintf(buf,"%d bytes",bytes);
	} else {
		FLOAT64 kbytes = bytes/1024;
		if( kbytes < 1024 ) {
			sprintf(buf,"%.2f KB",kbytes);
		} else {
			FLOAT64 mbytes = kbytes/1024.0;
			if( mbytes < 1024 ) {
				sprintf(buf,"%.2f MB",mbytes);
			} else {
				FLOAT64 gbytes = mbytes/1024.0;
				sprintf(buf,"%.2f GB",gbytes);
			}
		}
	}
	return buf;
}
#endif
