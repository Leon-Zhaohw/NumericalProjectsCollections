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
// ThreadSpecificData.cpp: Implementation
//
//////////////////////////////////////////////////////////////////////

#include "ThreadSpecificData.h"

#include <omp.h>

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
template <class T>
ThreadSpecificData<T>::ThreadSpecificData()
	: _useDefault( false )
{
}

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
template <class T>
ThreadSpecificData<T>::ThreadSpecificData( T defaultValue )
	: _useDefault( true ),
		_defaultValue( defaultValue )
{
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
template <class T>
ThreadSpecificData<T>::~ThreadSpecificData()
{
}

//////////////////////////////////////////////////////////////////////
// Returns the data associated with the current thread
//////////////////////////////////////////////////////////////////////
template <class T>
T & ThreadSpecificData<T>::get()
{
#if 0
#ifdef USING_OMP
	int threadNum = omp_get_thread_num();
	int numThreads = omp_get_max_threads();
#else
	int threadNum = 0;
	int numThreads = 1;
#endif
#endif
	int threadNum = omp_get_thread_num();

#if 0
	// Enlarge the data array if we need to.  This
	// needs to be done atomically
#pragma omp critical
	{
		if ( numThreads > _data.size() )
		{
			_data.resize( numThreads );
		}
	}

	return _data[ threadNum ];
#endif
#pragma omp critical
	{
		if ( _data.find( threadNum ) == _data.end() )
		{
			if ( _useDefault )
			{
				_data[ threadNum ] = _defaultValue;
			}
			else
			{
				_data[ threadNum ] = T();
			}
		}
	}

	return _data[ threadNum ];
}
