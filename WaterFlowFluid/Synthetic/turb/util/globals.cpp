/******************************************************************************
 *  DDF Fluid solver with Turbulence extensions
 *
 *  copyright 2009 Nils Thuerey, Tobias Pfaff
 * 
 *  DDF is free software, distributed under the GNU General Public License (GPL v2).
 *  See the file COPYING for more information.
 *
 * Global functions for logging etc.
 *
 *****************************************************************************/

#include "globals.h"
#include "grid.h"

#include <stdlib.h>
#include <limits.h>
#include <iostream>
#include <sstream>
#include "waveletnoise.h"
#ifdef WIN32
// for timing
#include <windows.h>
#else
#include <time.h>
#include <sys/time.h>
#include <sys/times.h>
#endif

// debugging, global id Counter
long int idCounter = 0;

namespace DDF {

// global debug level 
// 0=all off, 1=err only, 2=err+warn, >=3=all on
int gDebugLevel = 10;

// globals for grids.h
int gPatchIdCounter = 0;
int gGridIdCounter = 0;
class PatchManager;
PatchManager *gpPatchManager = NULL;

// globals for operators.h, misc. debugging profiling
int gGridPassCounter = 0;
double gPatchCount[10] = {0,0,0,0,0, 0,0,0,0,0}; 
double gSerialPatchCnt = 0.;
unsigned int gCountApplyOp = 0;
unsigned int gCountApplyOpNoFlag = 0;
unsigned int gCountApplyOpSimple = 0;
unsigned int gCountPatchLocks[10] = {0,0,0,0,0, 0,0,0,0,0}; 

// globals for vectorbase.h
const Real gVecEpsilon = VECTOR_EPSILON;
// instantiate zero vectors
template<> const ntlVector3Dim<float > ntlVector3Dim<float >::ZERO( 0.f, 0.f, 0.f );
template<> const ntlVector3Dim<double> ntlVector3Dim<double>::ZERO( 0., 0., 0. );
template<> const ntlVector3Dim<int   > ntlVector3Dim<int   >::ZERO( 0, 0, 0 );

// globals for waveletnoise.h
float* stdTileData = NULL;

//! for interval debugging output
myTime_t globalIntervalTime = 0;
//! color output setting for messages (0==off, else on)
#ifdef WIN32
// switch off first call (=-1)
#define DEF_globalColorSetting -1 
#else // WIN32
// linux etc., on by default (=1)
#define DEF_globalColorSetting -1 
#endif // WIN32
int globalColorSetting = DEF_globalColorSetting; // linux etc., on by default
int globalFirstEnvCheck = 0;
void resetGlobalColorSetting() { globalColorSetting = DEF_globalColorSetting; }

// global string for formatting vector output, TODO test!? 
#if DDF_DIMENSION==3
//const char *globVecFormatStr = "[%4.2f,%4.2f,%4.2f]";
const char *globVecFormatStr = "[%6.4f,%6.4f,%6.4f]";
#else // DDF_DIMENSION==3
const char *globVecFormatStr = "[%4.2f,%4.2f]";
#endif

// global state & error
int gSimState = 0;
char gSimErrorString[256];

void setSimState(int set) { gSimState = set; }
int  getSimState(void) { return gSimState; }
int  isSimworldOk(void) { return gSimState>=0; }

void setSimErrorString(char* set) { strncpy(gSimErrorString,set,256); }
char* getSimErrorString(void) { return gSimErrorString; }

// get string of global settings
string getSettingsString() {
	std::ostringstream ret;
	ret << 
		"Settings"<<
			//", Bits:"<<((sizeof(int*)==8) ? string("64") : string("32") )<<
			", B:"<<sizeof(int*)<<
			", Dim:"<<DDF::gDim<<
			", Patches: OFF"<<
			", PatchSize:"<<DDF::gPatchSize<<
			", fp.Prec.:"<<FLOATINGPOINT_PRECISION<<
			", Omp:"<<DDF_OPENMP<<
			" "
		;
#if DDF_OPENMP==1
	//omp_set_num_threads(2); // test
#pragma omp parallel default(shared) 
	{ // omp region
		const int id = omp_get_thread_num();
		const int Nthrds = omp_get_num_threads(); 
		if(id==0) ret << ", #threads:"<<Nthrds;
	}
#endif // DDF_OPENMP==1
	return ret.str();
}

//-----------------------------------------------------------------------------
// helper function that converts a string to integer, 
// and returns an alternative value if the conversion fails
int convertString2Int(const char *str, int alt)
{
	int val;
	char *endptr;
	bool success=true;

	val = strtol(str, &endptr, 10);
	if( (str==endptr) ||
			((str!=endptr) && (*endptr != '\0')) ) success = false;

	if(!success) {
		return alt;
	}
	return val;
}

//-----------------------------------------------------------------------------
//! helper function that converts a flag field to a readable integer
string convertFlags2String(int flags) {
	std::ostringstream ret;
	ret <<"(";
	int max = sizeof(int)*8;
	for(int i=0; i<max; i++) {
		if(flags & (1<<31)) ret <<"1";
		else ret<<"0";
		if(i<max-1) {
			//ret << ",";
			if((i%8)==7) ret << " ";
		}
		flags = flags << 1;
	}	
	ret <<")";
	return ret.str();
}

//-----------------------------------------------------------------------------
// helper function to determine current time
myTime_t getTime()
{
	myTime_t ret = 0;
#ifdef WIN32
	LARGE_INTEGER liTimerFrequency;
	QueryPerformanceFrequency(&liTimerFrequency);
	LARGE_INTEGER liLastTime;
	QueryPerformanceCounter(&liLastTime);
	ret = (INT)( ((double)liLastTime.QuadPart / liTimerFrequency.QuadPart)*1000 ); // - mFirstTime;
#else
	struct timeval tv;
	struct timezone tz;
	tz.tz_minuteswest = 0;
	tz.tz_dsttime = 0;
	gettimeofday(&tv,&tz);
 	ret = (tv.tv_sec*1000)+(tv.tv_usec/1000); //-mFirstTime;
#endif
	return (myTime_t)ret;
}

//-----------------------------------------------------------------------------
// convert time to readable string
string getTimeString(myTime_t usecs) {
	std::ostringstream ret;
	//myTime_t us = usecs % 1000;
	myTime_t ms = (myTime_t)(   (double)usecs / (60.0*1000.0)  );
	myTime_t ss = (myTime_t)(  ((double)usecs / 1000.0) - ((double)ms*60.0)  );
	int      ps = (int)(       ((double)usecs - (double)ss*1000.0)/10.0 );

 	//ret.setf(ios::showpoint|ios::fixed);
 	//ret.precision(5); ret.width(7);

	if(ms>0) {
		ret << ms<<"m"<< ss<<"s" ;
	} else {
		if(ps>0) {
			ret << ss<<".";
			if(ps<10) { ret <<"0"; }
			ret <<ps<<"s" ;
		} else {
			ret << ss<<"s" ;
		}
	}
	return ret.str();
}

//! helper to check if a bounding box was specified in the right way
bool checkBoundingBox(Vec3 s, Vec3 e, string checker) {
	if( (s[0]>e[0]) ||
			(s[1]>e[1]) ||
			(s[2]>e[2]) ) {
		errFatal("checkBoundingBox","Check by '"<<checker<<"' for BB "<<s<<":"<<e<<" failed! Aborting...",SIMWORLD_INITERROR);
		return 1;
	}
	return 0;
}



//-----------------------------------------------------------------------------
// debug message output

static string col_black ( "\033[0;30m");
static string col_dark_gray ( "\033[1;30m");
static string col_bright_gray ( "\033[0;37m");
static string col_red ( "\033[0;31m");
static string col_bright_red ( "\033[1;31m");
static string col_green ( "\033[0;32m");
static string col_bright_green ( "\033[1;32m");
static string col_bright_yellow ( "\033[1;33m");
static string col_yellow ( "\033[0;33m");
static string col_cyan ( "\033[0;36m");
static string col_bright_cyan ( "\033[1;36m");
static string col_purple ( "\033[0;35m");
static string col_bright_purple ( "\033[1;35m");
static string col_neutral ( "\033[0m");
static string col_std = col_bright_gray;

std::ostringstream globOutstr;
bool               globOutstrForce=false;
#define DM_NONE      100
void messageOutputForce(string from) {
	bool org = globOutstrForce;
	globOutstrForce = true;
	messageOutputFunc(from, DM_NONE, "\n", 0);
	globOutstrForce = org;
}

void messageOutputFunc(string from, int id, string msg, myTime_t interval) {
	// fast skip
	if((id!=DM_FATAL)&&(gDebugLevel<=0)) return;

	if(interval>0) {
		myTime_t currTime = getTime();
		if((currTime - globalIntervalTime)>interval) {
			globalIntervalTime = getTime();
		} else {
			return;
		}
	}

	// colors off?
	if( (globalColorSetting == -1) || // off for e.g. win32 
		  ((globalColorSetting==1) && ((id==DM_FATAL)||( getenv("DDF_NOCOLOROUT") )) )
		) {
		// only reset once
		col_std = col_black = col_dark_gray = col_bright_gray =  
		col_red =  col_bright_red =  col_green =  
		col_bright_green =  col_bright_yellow =  
		col_yellow =  col_cyan =  col_bright_cyan =  
		col_purple =  col_bright_purple =  col_neutral =  "";
		globalColorSetting = 0;
	}

	std::ostringstream sout;

	if(id==DM_FATAL) {
		// dont print?
		if(gDebugLevel==0) return;
		sout << "\n\n\n"; // add newline for output
		sout << " !!!!!!!!!!!!!!!!!!!!!!!! \n"; // add newline for output
		sout << " !!!!!!!!!!!!!!!!!!!!!!!! \n"; // add newline for output
		sout << " !!!!!!!!!!!!!!!!!!!!!!!! \n"; // add newline for output
	}

	sout << col_cyan<< from;
	switch(id) {
		case DM_MSG:
			sout << col_std << ":";
			break;
		case DM_WARNING:
			sout << col_bright_red << " warning:" << col_std;
			break;
		case DM_ERROR:
			sout << col_red << " error:" << col_red;
			break;
		case DM_FATAL:
			sout << col_red << " fatal("<<gSimState<<"):" << col_red;
			break;
		default:
			// this shouldnt happen...
			sout << col_red << " --- messageOutputFunc error: invalid id ("<<id<<") --- aborting... \n\n" << col_std;
			break;
	}
	sout <<" "<< msg << col_std;

	if(id==DM_FATAL) {
		strncpy(gSimErrorString, sout.str().c_str(), 256);
		// dont print?
		if(gDebugLevel==0) return;
		sout << " ! \n"; // add newline for output
		sout << " !!!!!!!!!!!!!!!!!!!!!!!! \n"; // add newline for output
		sout << " !!!!!!!!!!!!!!!!!!!!!!!! \n"; // add newline for output
		sout << " !!!!!!!!!!!!!!!!!!!!!!!! \n\n\n\n\n\n\n"; // add newline for output
		sout << "\n"; // add newline for output
	}

	// determine output - file==1/stdout==0 / globstr==2
	char filen[256];
	strcpy(filen,"debug_unini.txt");
	int fileout = 0;
	/*std::ostringstream mpin;
	if(glob_mpindex>=0) {
		mpin << "ddf_log_"<< glob_mpindex <<".txt";
	} else {
		mpin << "ddf_log_ini.txt";
	}
	fileout = 1;
	strncpy(filen, mpin.str().c_str(),255); filen[255]='\0';// */
	strncpy(filen, "ddf_debug_log.txt",255);

	//#ifdef WIN32
	// windows causes trouble with direct output
	//fileout = 1;
	//#endif // WIN32

#if WITH_OPENMP==1
	fileout = 2;// buffer out, switch off again...
	if(globOutstrForce) fileout=1;
#endif
	if(getenv("DDF_FORCESTDOUT")) {
		fileout = 0;// always direct out
	}
	//fprintf(stdout,"out deb %d, %d, '%s',l%d \n",globOutstrForce,fileout, filen, globOutstr.str().size() );

#if WITH_OPENMP==1
#pragma omp critical MESSAGE_DEBUG_OUT
#endif // WITH_OPENMP==1
	{
	if(fileout==1) {
		// debug level is >0 anyway, so write to file...
		FILE *logf = fopen(filen,"a+");
		// dont complain anymore here...
		if(logf) {
			if(globOutstrForce) {
				fprintf(logf, "%s",globOutstr.str().c_str() );
				globOutstr.str(""); // reset
			}
			fprintf(logf, "%s",sout.str().c_str() );
			fclose(logf);
		}
	} else if(fileout==2) {
		globOutstr << sout.str();
	} else {
		// normal stdout output
		fprintf(stdout, "%s",sout.str().c_str() );
		fflush(stdout); 
	}
	} // omp crit

	// don't segfault in gui
	if(id==DM_FATAL) {
		// DEBUG force segfault
		*(int *)(0) = 0;
		exit(999);
	}
#if DDF_GLUTGUI!=1
#endif
}

// grids.h implementations
template<> void Grid<int>::setDummyZero() { mZero = 0; } 
template<> void Grid<float>::setDummyZero() { mZero = 0.f; }
template<> void Grid<double>::setDummyZero() { mZero = 0.; }
template<> void Grid<Vec3>::setDummyZero() { mZero = Vec3(0.); }

} // DDF 

// global list
map<string, WAVELETNOISE::WaveletNoiseField*> gNoiseFields;



