//////////////////////////////////////////////////////////////////////
// This file is part of Closest Point Turbulence.
// 
// Closest Point Turbulence is free software: you can redistribute it 
// and/or modify it under the terms of the GNU General Public License 
// as published by the Free Software Foundation, either version 3 of 
// the License, or (at your option) any later version.
// 
// Closest Point Turbulence is distributed in the hope that it will 
// be useful, but WITHOUT ANY WARRANTY; without even the implied 
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Closest Point Turbulence. 
// If not, see <http://www.gnu.org/licenses/>.
// 
// Copyright 2013 Theodore Kim and Nils Thuerey
//////////////////////////////////////////////////////////////////////

#ifndef TIMER_H
#define TIMER_H

#include <iostream>
#include <string>
#include <map>
#include <stack>

#include <sys/time.h>

// To time a function, just put:
//
//  TIMER functionTimer(__FUNCTION__);
//
// at the beginning of the function

class TIMER
{
public:
  // start the timer by default -- if a tick is called later,
  // it will just stomp it
  TIMER(std::string blockName); 
  ~TIMER();

  static double timing(timeval& begin = _tick, timeval& end  = _tock) {
    double beginTime = (double)begin.tv_sec + 1e-6 * begin.tv_usec;
    double endTime = (double)end.tv_sec + 1e-6 * end.tv_usec;
    return endTime - beginTime;
  };
  static int hours(int seconds) { return seconds / (60 * 60); };
  static int minutes(int seconds) {
   int mod = seconds % (60 * 60);
   return mod / 60;
  };
  static int seconds(int seconds) {
    int mod = seconds % (60 * 60);
    return mod % 60;
  };

  static void printTimings();
  static void printTimingsPerFrame(const int frames);

private:
  // begin and end of current block being timed
  static timeval _tick;
  static timeval _tock;

  // hash table of all timings
  static std::map<std::string, double> _timings;

  // call stack
  static std::stack<std::string> _callStack;
};

#endif
