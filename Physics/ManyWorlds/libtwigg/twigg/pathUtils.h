#ifdef WIN32
#pragma once
#endif

std::string toNativePath( const std::string& path );
std::string toUniversalPath( const std::string& path );
std::string pathForFile( const std::string& fn );
std::string fileInPath( const std::string& filename, const std::string& path );
