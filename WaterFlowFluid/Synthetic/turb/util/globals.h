/******************************************************************************
 *  DDF Fluid solver with Turbulence extensions
 *
 *  copyright 2009 Nils Thuerey, Tobias Pfaff
 * 
 *  DDF is free software, distributed under the GNU General Public License (GPL v2).
 *  See the file COPYING for more information.
 *
 * Declarations for global functions and variables
 *
 *****************************************************************************/


#ifndef UTILITIES_H
#define UTILITIES_H


// "brute force" global defines
#ifndef DDF_DEBUG
#define DDF_DEBUG 0
#endif // DDF_DEBUG

// dimension of solver 
#ifndef DDF_DIMENSION
#define DDF_DIMENSION 3
#endif 

// openmp
#ifndef DDF_OPENMP
#define DDF_OPENMP 0
#endif 

#ifdef WIN32
// If windows.h is included, it defines a macro min() and max(),
// which can lead to all sorts of weird problems, e.g. conflicts
// with the STL algorithms min() and max(). By #defining NOMINMAX,
// we prevent windows.h from defining these macros.
#define NOMINMAX
#endif

// this serves as the main include file
// for all kinds of stuff that might be required.
// under windos there seem to be strange 
// errors when including the STL header too
// late...
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <math.h>
#include <string.h>
#include <stdio.h>
#if DDF_OPENMP==1
#include "omp.h"
#endif

// hack for MSVC6.0 compiler
#ifdef _MSC_VER
#	if _MSC_VER < 1300
#		define for     if(false); else for
#		define map     std::map
#		define vector  std::vector
#		define string  std::string
		// use this define for MSVC6 stuff hereafter
		//#		define USE_MSVC6FIXES
#	else 
		// _MSC_VER >= 1300 , 7.0 or higher
		using std::map;
		using std::vector;
		using std::string;
#	endif
#else // not MSVC6
	// for proper compilers...
	using std::map;
	using std::vector;
	using std::string;
#endif // MSVC6

#ifdef __APPLE_CC__
	// apple
#else

#	ifdef WIN32
		// windows values missing, see below
#		ifndef snprintf
#		define snprintf _snprintf
#		endif
	// rounding function
#		ifndef round
		inline double round(double x) { return floor(x + 0.5); }
#		endif
#		ifndef false
#		define false 0
#		endif
#		ifndef true
#		define true 1
#		endif
#	else 
		// not WIN32
		// floating point limits for linux,*bsd etc...
		#include <float.h>
#	endif // WIN32

#endif // __APPLE_CC__

// standard includes
#include "vectorbase.h"

namespace DDF {

// dimension variable
const int gDim = DDF_DIMENSION;

// debugging outputs , debug level 0 (off) to 10 (max) 
extern int gDebugLevel;

// get string of global settings
string getSettingsString();

// time measurements
typedef unsigned long myTime_t;


// state of the simulation world
// default state
#define SIMWORLD_OK       0
// error during init
#define SIMWORLD_INITERROR    -1
// error during simulation
#define SIMWORLD_PANIC        -2
// general error 
#define SIMWORLD_GENERICERROR -3
// error in grid handling
#define SIMWORLD_GRIDERROR    -4
// error in grid handling
#define SIMWORLD_ERRPARSE     -5
// error with plugins
#define SIMWORLD_PLUGINERROR  -6

// access global state of elbeem simulator
void setSimState(int set);
int  getSimState(void);
int  isSimworldOk(void);

// access elbeem simulator error string
void setSimErrorString(char* set);
char* getSimErrorString(void);

// debug messages on/off
#define USE_DEBUG_MESSAGES 1
#if USE_DEBUG_MESSAGES==1

/* debug output function */
#define DM_MSG        1
#define DM_WARNING    2
#define DM_ERROR      3
#define DM_FATAL      4
void messageOutputFunc(string from, int id, string msg, myTime_t interval);

// debugging messages defines 
//#define MSGSTREAM std::ostringstream msg; msg.precision(15); msg.width(17);
#define MSGSTREAM std::ostringstream msg; msg.precision(7); msg.width(9);

#define errMsg(from,mStr)      if(DDF::gDebugLevel>0){ MSGSTREAM; msg << mStr <<"\n"; DDF::messageOutputFunc(from, DM_ERROR,   msg.str(), 0); }
#define warnMsg(from,mStr)     if(DDF::gDebugLevel>1){ MSGSTREAM; msg << mStr <<"\n"; DDF::messageOutputFunc(from, DM_WARNING, msg.str(), 0); }
#define debMsg(from,mStr)      if(DDF::gDebugLevel>2){ MSGSTREAM; msg << mStr <<"\n"; DDF::messageOutputFunc(from, DM_MSG, msg.str(), 0); }

#else
// no messages at all...
#	define debMsg(mStr) 
#	define debMsgId(from,id,mStr,level)
#	define errMsg(from,mStr)
#	define warnMsg(from,mStr)
#endif


// fatal errors - have to be handled 
#define errFatal(from,mStr,errCode) { \
	DDF::setSimState(errCode); \
	MSGSTREAM; msg << mStr; \
	DDF::messageOutputFunc(from, DM_FATAL, msg.str(), 0); \
}


//! helper function that converts a string to integer
int convertString2Int(const char *str, int alt);

//! helper function that converts a flag field to a readable integer
string convertFlags2String(int flags);

//! get the current system time
myTime_t getTime();
//! convert time to readable string
string getTimeString(myTime_t usecs);

//! helper to check if a bounding box was specified in the right way
bool checkBoundingBox(Vec3 s, Vec3 e, string checker);

//! reset color output for elbeem init
void resetGlobalColorSetting();


/*! print some vector from 3 values e.g. for ux,uy,uz */
#define PRINT_VEC(x,y,z) " ["<<(x)<<","<<(y)<<","<<(z)<<"] "

/*! print i,j,k as a vector, as we need ijk all the time */
#define PRINT_IJK PRINT_VEC(i,j,k)


// write png image
int writeImage(const char *fileName, unsigned char **rowsp, int w, int h);


// template helper funcs

// minimum 
template < class T >
inline T
MIN( T a, T b )
{ return (a < b) ? a : b ; }

// maximum 
template < class T >
inline T
MAX( T a, T b )
{ return (a < b) ? b : a ; }

// maximum of a vector 
template < class T >
inline T
VMAX( ntlVector3Dim<T> a )
{ 
	if(a[0]<a[1]) {
		if(a[1]<a[2]) return a[2];
		else          return a[1];
	} else {
		if(a[0]<a[2]) return a[2];
		else          return a[0];
	}
}
// minimum of a vector 
template < class T >
inline T
VMIN( ntlVector3Dim<T> a )
{ 
	if(a[0]>a[1]) {
		if(a[1]>a[2]) return a[2];
		else          return a[1];
	} else {
		if(a[0]>a[2]) return a[2];
		else          return a[0];
	}
}

// swap two values
template < class T >
inline void
TRISWAP( T& a, T& b )
{ 
	T tmp = a;
	a = b;
	b = tmp;
}

// absolute value 
template < class T >
inline T
ABS( T a )
{ return (0 < a) ? a : -a ; }

// sign of the value 
template < class T >
inline T
SIGNUM( T a )
{ return (0 < a) ? 1 : -1 ; }

// sign, returns -1,0,1 depending on sign/value=0 
template < class T >
inline T
SIGNUM0( T a )
{ return (0 < a) ? 1 : ( a < 0 ? -1 : 0 ) ; }

// round to nearest integer 
inline int
ROUND(double d)
{ return int(d + 0.5); }

// square function 
template < class T >
inline T
SQUARE( T a )
{ return a*a; }

// cube function 
template < class T >
inline T
CUBEFUNC( T a )
{ return a*a*a; }


// helper class, track min,max and avg of some scalar type T
template<class T> class MinMaxTracker {
public:
	MinMaxTracker() {
		mMin = 1e20;
		mMax = -1e20; 
		mAccumulation = 0.0;
		mCnt = 0;
	};
	~MinMaxTracker() {};

	// take into account value
	inline void check(T& val) {
		if (mMin>val) mMin=val;
		if (mMax<val) mMax=val;
		mAccumulation += val;
		mCnt++;
	}

	// in out put 

	T getMin() {return mMin;}
	T getMax() {return mMax;}
	T getAvg() {return mAccumulation/(double)(mCnt);}

	std::string str() {
		std::ostringstream out;
		out <<"min="<<mMin<<", max="<<mMax<<", avg="<<getAvg();
		return out.str();
	}

protected:
	T mMin,mMax, mAccumulation; 
	int mCnt;
}; // MinMaxTracker

template<class T>
std::string helperPrintMap(std::map<std::string, T*>& pmap, std::string caller) {
	std::ostringstream out;
	//debMsg("print map","for '"<<caller<<"', "<<pmap.size()<<" entries ");
	out << "print map for '"<<caller<<"', "<<pmap.size()<<" entries ";
	int cnt=0;
	for(typename std::map<std::string, T*>::iterator iter=pmap.begin(); iter != pmap.end(); iter++) {
		long int pntval = (long int)( (*iter).second );
		//debMsg("entry","i "<<cnt<<"/"<<pmap.size()<<" name '"<<(*iter).first <<"' is " << pntval );
		out << "entry "<<cnt<<"/"<<pmap.size()<<" name '"<<(*iter).first <<"' is " << pntval<<", ";
		cnt++;
	}
	return out.str();
} // helperPrintMap


//*****************************************************************************
// utility defines

#define FOR_IJK(kmin,kmax,jmin,jmax,imin,imax) \
	for (int k=(kmin); k<(kmax); k++) \
		for (int j=(jmin); j<(jmax); j++) \
			for (int i=(imin); i<(imax); i++)  /* loop body */
#define FOR_IJK_GRID(grid) \
	FOR_IJK( (grid)->getMinZLoopValue(), (grid)->getMaxZLoopValue(),   0, (grid)->getSizeY(),   0, (grid)->getSizeX() ) 
#define FOR_IJK_GRID_BND(grid,bnd) \
	FOR_IJK( (grid)->getMinZLoopValue(bnd), (grid)->getMaxZLoopValue(bnd),   0+(bnd), (grid)->getSizeY()-(bnd),   0+(bnd), (grid)->getSizeX()-(bnd) ) 
#define FOR_IJK_VEC(vecmin, vecmax) \
	FOR_IJK(vecmin[2],vecmax[2],vecmin[1],vecmax[1],vecmin[0],vecmax[0]) 

// reverse loops
#define FOR_IJKREV(kmin,kmax,jmin,jmax,imin,imax) \
	for (int k=(kmin); k>=(kmax); k--) \
		for (int j=(jmin); j>=(jmax); j--) \
			for (int i=(imin); i>=(imax); i--)  /* loop body */
#define FOR_IJKREV_GRID(grid) \
	FOR_IJKREV( (grid)->getMaxZLoopValue()-1,(grid)->getMinZLoopValue(),   (grid)->getSizeY()-1, 0,   (grid)->getSizeX()-1, 0  ) 

//*****************************************************************************
// clamping helper functions

// clamp a in range b<= a <= c
// write to a
template < class T >
inline void
CLAMP( T &a, T b , T c )
{ 
	if(a < b) a=b;
	else if(a > c) a=c;
}

// clamp , and return value
template < class T >
inline T
CLAMP_RET( T a, T b , T c ) { 
	if(a < b) return b;
	else if(a > c) return c;
	return a;
}

template < class T >
inline T
CLAMP_DOWN_RET( T a, T b ) { 
	if(a < b) return b;
	return a;
}

template < class T >
inline T
CLAMP_UP_RET( T a, T b ) { 
	if(a > b) return b;
	return a;
}

// clamp and return whether clamping occured
template < class T >
inline bool
CLAMP_BOOL( T &a, T b , T c )
{ 
	if(a < b) {
		a=b; 
		return true;
	}

	if(a > c) {
		a=c;
		return true;
	}
	return false;
}

// clamp to grid boundaries
template < class T >
inline void CLAMP_TO_GRID(int &i, int &j, int &k, T grid) {
	CLAMP(i, 0, grid->getSizeX()-2 );
	CLAMP(j, 0, grid->getSizeY()-2 );
	CLAMP(k, 0, grid->getSizeZ()-2 );
}

// clamp to grid and return whether clamping occured
template < class T >
inline bool CLAMP_TO_GRID_BOOL(int &i, int &j, int &k, T grid) {
	bool ret = false;
	ret |= CLAMP_BOOL(i, 0, grid->getSizeX()-2 );
	ret |= CLAMP_BOOL(j, 0, grid->getSizeY()-2 );
	ret |= CLAMP_BOOL(k, 0, grid->getSizeZ()-2 );
	return ret;
}


// constants for easier neighbor traversal
static const int nbX[6] = {-1,1,   0,0,0,0};
static const int nbY[6] = { 0,0, -1,1, 0,0};
static const int nbZ[6] = { 0,0,0,0,  -1,1,};


} // DDF 

#endif

