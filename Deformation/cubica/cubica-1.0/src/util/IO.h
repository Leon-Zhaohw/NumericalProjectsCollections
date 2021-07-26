/*
This file is part of Cubica.
 
Cubica is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Cubica is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Cubica.  If not, see <http://www.gnu.org/licenses/>.
*/
//////////////////////////////////////////////////////////////////////
// IO.h: Some helpful I/O routines
//
//////////////////////////////////////////////////////////////////////

#ifndef IO_H
#define IO_H

#include <MATRIX.h>
#include <MATRIX3.h>
#include <VECTOR.h>

#include <SETTINGS.h>

#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>

#include <ctime>
#include <cstring>
#include <time.h>
#include <stdlib.h>// for _MAX_PATH
#include <stdio.h>
#include <sstream>

#if defined(__unix__) || defined (__LINUX__)
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <unistd.h>
#endif

#ifdef WIN32
#include <windows.h>
#endif

using namespace std;

#define SDUMP(x)	" " << #x << "=[ " << x << " ] "

#if 0
class IO {
	public:
		// Save an array of vectors
		static void saveVec3Data( const Vector3Array &data, const char *fileName )
		{
			// Use the existing read/write functionality in VECTOR
			VECTOR v_data( data.size() * 3 );

			for ( int i = 0; i < data.size(); i++ )
			{
				v_data( i * 3 + 0 ) = data[i][0];
				v_data( i * 3 + 1 ) = data[i][1];
				v_data( i * 3 + 2 ) = data[i][2];
			}

			v_data.write( fileName );
		}

		// Load an array of vectors
		static void loadVec3Data( Vector3Array &data, const char *fileName )
		{
			// Use the existing read/write functionality in VECTOR
			VECTOR v_data;

			v_data.read( fileName );

			data.resize( v_data.size() / 3 );

			for ( int i = 0; i < v_data.size() / 3; i++ )
			{
				data[i][0] = v_data( i * 3 + 0 );
				data[i][1] = v_data( i * 3 + 1 );
				data[i][2] = v_data( i * 3 + 2 );
			}
		}

		// Save an array of scalar data
		static void saveScalarData( const FloatArray &data, const char *fileName )
		{
			// Use the existing read/write functionality in VECTOR
			VECTOR v_data( data.size() );

			for ( int i = 0; i < data.size(); i++ )
			{
				v_data(i) = data[i];
			}

			v_data.write( fileName );
		}

		// Load an array of scalar data
		static void loadScalarData( FloatArray &data, const char *fileName )
		{
			// Use the existing read/write functionality in VECTOR
			VECTOR v_data;

			v_data.read( fileName );

			data.resize( v_data.size() );

			for ( int i = 0; i < v_data.size(); i++ )
			{
				data[i] = v_data(i);
			}
		}

		static vector<string> split( const string& s, const string& splitters )
		{
			vector<string> splat;
			int start = 0;
			while( 1 )
			{
				int occur = s.find_first_of( splitters, start );

				if( occur == string::npos ) {
					// we're done. add the last string
					splat.push_back( s.substr( start, string::npos ) );
					break;
				}
				else {
					splat.push_back( s.substr( start, occur-start ) );
					start = occur + 1;
				}	
			}

			return splat;
		}

		//----------------------------------------
		//  Tries to create a directory of the given path.
		//  Returns false if the directory already existed.
		//----------------------------------------
		enum CREATE_DIR_RESULT { OK, EXISTS, FAIL };

#if defined( WIN32 )

		static enum CREATE_DIR_RESULT create_dir( const string& name )
		{
			wchar_t* wideStr = NULL;
			size_t numWideChars = 0;
			BOOL ok = false;

			wideStr = new wchar_t[ name.size() ];
			mbstowcs_s( &numWideChars, wideStr, name.size()+1,
									name.c_str(), _TRUNCATE );

			string temp = wchar2string( wideStr );
			ok = CreateDirectory( temp.c_str(), NULL );
			delete wideStr;

			if( ok )
			{
				return OK;
			}
			else
			{
				cerr << "** Creating directory " << SDUMP(name) << " failed!" << endl;
				// TODO - use GetLastError to get more detailed errror info
				return FAIL;
			}
		}

#elif defined( __unix__ )

		static enum CREATE_DIR_RESULT create_dir( const string& name )
		{
			int status = mkdir( name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );

			if( status == 0 )
			{
				return OK;
			}
			else
			{
				switch( errno )
				{
					case EEXIST: return EXISTS;
					default:
					{
						cerr << "** Creating directory " << SDUMP(name);
						cerr << " failed!" << endl;
						return FAIL;
					}
				}
			}
		}

#endif

};
#endif

#endif
