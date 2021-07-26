// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#ifdef WIN32
#pragma once
#endif


//#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers
// Windows Header Files:

#ifdef _DEBUG
#define CRTDBG_MAP_ALLOC
#define _CRTDBG_MAP_ALLOC
#include <cstdlib>
#include <crtdbg.h>

#define DEBUG_NEW new(_NORMAL_BLOCK ,__FILE__, __LINE__)
#else
#define DEBUG_NEW new
#endif

#ifdef WIN32
#include <windows.h>
#include <tchar.h>
#endif


// TODO: reference additional headers your program requires here
#include "twigg/util.h"
#include "twigg/exception.h"

#ifndef CONDOR
#include <GL/glew.h>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
#endif

#ifdef WIN32
#pragma warning (disable:4786)
#endif

#include <algorithm>
#include <string>
#include <vector>
#include <deque>
#include <memory>

#include <boost/smart_ptr.hpp>
#include <boost/numeric/conversion/bounds.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/limits.hpp>
#include <boost/mem_fn.hpp>

// threading support
#ifndef CONDOR
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#endif

#include <boost/array.hpp>


