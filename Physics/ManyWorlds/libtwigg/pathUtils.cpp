#include "stdafx.h"
#include "twigg/pathUtils.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

std::string toNativePath( const std::string& path )
{
	// FLTK tends to use '/' always but this confuses some other libraries;
	// this routine converts everything to the proper format.
	std::string result;
	result.reserve( path.size() );

#ifdef _WIN32
	char pathSeparator = '\\';
#else
	char pathSeparator = '/';
#endif

	for( std::string::const_iterator itr = path.begin();
		itr != path.end();
		++itr )
	{
		if( *itr == '/' || *itr == '\\' )
			result.push_back( pathSeparator );
		else
			result.push_back( *itr );
	}

	return result;
}

std::string toUniversalPath( const std::string& path )
{
	// switch to '/'s everywhere.
	std::string result;
	result.reserve( path.size() );
	for( std::string::const_iterator itr = path.begin();
		itr != path.end();
		++itr )
	{
		if( *itr == '\\' )
			result.push_back( '/' );
		else
			result.push_back( *itr );
	}

	return result;
}

std::string pathForFile( const std::string& fn )
{
	std::string filename = toUniversalPath(fn);

	std::string path;

	// for now assume Unix-style paths since this is what FLTK returns.
	std::string::size_type pos = filename.rfind( '/' );
	if( pos != std::string::npos )
		path = std::string( filename.begin(), filename.begin() + pos );

	return path;
}

std::string fileInPath( const std::string& filename, const std::string& path )
{
	if( path.empty() )
		return toNativePath(filename);

	if( (*path.rbegin()) == '/' )
		return toNativePath(path + filename);

	return toNativePath(path + std::string("/") + filename);
}
