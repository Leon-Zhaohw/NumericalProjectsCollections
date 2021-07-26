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
// ThreadSpecificData.h: Interface for the ThreadSpecificData class
//
//////////////////////////////////////////////////////////////////////

#ifndef THREAD_SPECIFIC_DATA_H
#define THREAD_SPECIFIC_DATA_H

#include <MATRIX.h>
#include <MATRIX3.h>
#include <VECTOR.h>

#include <SETTINGS.h>

#include <vector>
#include <set>

//////////////////////////////////////////////////////////////////////
// ThreadSpecificData class
//
// A generic container for data which should be replicated
// across threads.
// The template class T must provide a default constructor
// with no arguments
//////////////////////////////////////////////////////////////////////
template <class T>
class ThreadSpecificData {
	public:
		ThreadSpecificData();
		ThreadSpecificData( T defaultValue );

		// Destructor
		virtual ~ThreadSpecificData();

		// Returns the data associated with the current thread
		T &get();

		// Returns all data
		std::map<int, T> &getAll() { return _data; }

	private:
#if 0
		std::vector<T>				_data;
#endif
		std::map<int, T>			_data;

		bool									_useDefault;
		T											_defaultValue;
};

#include "ThreadSpecificData.cpp"

#endif
