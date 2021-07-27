/*
 This file is part of SSFR (Zephyr).
 
 Zephyr is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 Zephyr is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Zephyr.  If not, see <http://www.gnu.org/licenses/>.
 
 Copyright 2013 Theodore Kim
 */
// 
// Simple config file parser
//
// Original source code was provided courtesy of 
// Nils Theurey (www.ntoken.com)
//
// Minor modifications started on 6/27/2008 by Ted Kim
//
//////////////////////////////////////////////////////////////////////

#ifndef SIMPLE_PARSER_H
#define SIMPLE_PARSER_H

#include <iostream>
#include <string>
#include <vector>
#include <map>

using namespace std;

// read the whole test file of pairs
// a = b
// and look up values with conversion
class SIMPLE_PARSER  
{
	public:
		SIMPLE_PARSER(std::string file);
		virtual ~SIMPLE_PARSER();

		// get values from file
		int      getInt(string name,    int defaultValue,    bool needed=false);
		bool     getBool(string name,   bool defaultValue,   bool needed=false);
		double   getFloat(string name,  double defaultValue, bool needed=false);
    string   getString(string name, string defaultValue, bool needed=false);

		// check if there were any unused pairs
		bool haveUnusedValues();
		string printAllUnused();

    // check if the a parameters was specified in the config file
    bool defined(string name);

	protected:
		// value pair storage
		map<string, string> mVals;
		// value pair check if used...
		map<string, bool> mUsed;

		template<class T> T getScalarValue(string name, T defaultValue, bool needed=false);
};

#endif

