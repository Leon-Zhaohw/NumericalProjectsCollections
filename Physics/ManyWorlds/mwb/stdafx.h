// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#ifdef WIN32
#pragma once

#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers
// Windows Header Files:

#include <windows.h>
#include <winsock2.h>
#include <tchar.h>
#include <wmsdk.h>
#endif

#include <boost/pool/object_pool.hpp>

#ifdef _DEBUG
#define CRTDBG_MAP_ALLOC
#define _CRTDBG_MAP_ALLOC
#include <cstdlib>
#include <crtdbg.h>
#define DEBUG_NEW new(_NORMAL_BLOCK ,__FILE__, __LINE__)
//#define DEBUG_NEW new
#else
#define DEBUG_NEW new
#endif

#ifdef GUI
#include <GL/glew.h>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
#endif // GUI

#ifdef GUI
#include "wx/wxprec.h"               // wxWidgets precompiled / standard headers
#endif

// TODO: reference additional headers your program requires here
#include "twigg/util.h"
#include "twigg/stlext.h"
#include "twigg/exception.h"


#pragma warning (disable:4786)
#include <algorithm>
#include <string>
#include <vector>
#include <deque>
#include <queue>
#include <memory>
#include <list>
#include <sstream>

#include <boost/smart_ptr.hpp>
#include <boost/numeric/conversion/bounds.hpp>
#include <boost/limits.hpp>
#include <boost/random.hpp>
#include <boost/bind.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/crc.hpp>

#include <boost/iterator_adaptors.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <boost/ptr_container/ptr_deque.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

// threading support
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/recursive_mutex.hpp>
#include <boost/thread/condition.hpp>
#include <boost/thread/once.hpp>
#include <boost/array.hpp>
#include <boost/mem_fn.hpp>
#include <boost/lexical_cast.hpp>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>

namespace planning {

template <class _Pair>
struct compare2nd
{
	bool operator()(const _Pair& left, const _Pair& right)
	{
		std::less<typename _Pair::second_type> comparator;
		return comparator(left.second, right.second);
	}
};

template <typename T>
T toDegrees( T value )
{
	return value * (180.0 / M_PI);
}

template <typename T>
T toRadians( T value )
{
	return value * (M_PI / 180.0);
}


#ifdef WIN32
template <typename T>
struct hash_ptr
	: stdext::hash_compare<T>
{
	size_t operator()( const T value ) const
	{
		BOOST_STATIC_ASSERT( sizeof(T) == sizeof(size_t) );
		const size_t* intVal = reinterpret_cast<const size_t*>(&value);
		return *intVal;
	}

	bool operator()( const T left, const T right ) const
	{
		return left < right;
	}
};
#else
template <typename T>
struct hash_ptr
{
	size_t operator()( const T value ) const
	{
		const size_t* intVal = reinterpret_cast<const size_t*>(&value);
		return *intVal;
	}
};
#endif

template <typename T>
unsigned int checksum( const T& value )
{
	boost::crc_32_type  result;
	const size_t size = sizeof(T);
	const char* bytes = reinterpret_cast<const char*>( &value );
	result.process_bytes( bytes, size );
	return result.checksum();
}

template <typename T1>
unsigned int checksum( const std::vector<T1>& values )
{
	boost::crc_32_type  result;
	if( values.empty() )
		return result.checksum();

	const size_t size = values.size() * sizeof(T1);
	const char* bytes = reinterpret_cast<const char*>( &values[0] );
	result.process_bytes( bytes, size );
	return result.checksum();
}

template <>
unsigned int checksum< std::string >( const std::string& value )
{
	boost::crc_32_type  result;
	if( value.empty() )
		return result.checksum();

	const size_t size = value.size();
	result.process_bytes( &value[0], size );
	return result.checksum();
}


} // namespace planning


