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
#include "TIMER.h"
#include <cassert>
#include <omp.h>
#include <cstdio>

using namespace std;

timeval TIMER::_tick;
timeval TIMER::_tock;
std::map<std::string, double> TIMER::_timings;
std::stack<std::string> TIMER::_callStack;

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
TIMER::TIMER(string blockName) : _stopped(false) 
{
  if (omp_in_parallel())
    return;

  // look at the back of the call stack,
  // if there's something there, then store its timing
  //
  // else, it's the first call, so set the global start
  if (_callStack.size() > 0)
  {
    string function = _callStack.top();
    gettimeofday(&_tock, 0);

    _timings[function] += timing();
  }
  
  _callStack.push(blockName);
  gettimeofday(&_tick, 0);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
TIMER::~TIMER()
{
  if (!_stopped) stop();

  /*
  if (omp_in_parallel())
    return;
  assert(_callStack.size() > 0);
  
  string function = _callStack.top();
  _callStack.pop();
  gettimeofday(&_tock, 0);

  _timings[function] += timing();
  gettimeofday(&_tick, 0);
  */
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TIMER::stop()
{
  if (omp_in_parallel())
    return;
  assert(_callStack.size() > 0);
  
  // remove this function from the internall tracked call stack
  string function = _callStack.top();
  _callStack.pop();

  // get the stopping time
  gettimeofday(&_tock, 0);

  // add the timing to the current global totals
  _timings[function] += timing();
  
  // store the current timer's timing, just in case
  _elapsed = timing(_tick, _tock);
  
  gettimeofday(&_tick, 0);

  // record we stopped so we don't double-stop
  _stopped = true;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TIMER::printTimings()
{
  timeval now;
  gettimeofday(&now, 0);

  double currentTimer = timing(_tick, now);
  string currentName;
  if (_callStack.size() > 0) 
    currentName = _callStack.top();

  // create an inverse map so that it will sort by time
  map<double, string> inverseMap;
  map<string, double>::iterator forwardIter;
  double totalTime = 0.0;
  for (forwardIter = _timings.begin(); forwardIter != _timings.end(); forwardIter++)
  {
    string name = forwardIter->first;
    double time = forwardIter->second;

    // if we are dinside this function, add the timing
    if (name.compare(currentName) == 0)
      time += currentTimer;

    inverseMap[time] = name;
    totalTime += time;
  }

  // print the map out backwards since it sorts from least to greatest
  cout << "===========================================================================" << endl;
  cout << " TIMING BREAKDOWN: " << endl;
  cout << "===========================================================================" << endl;
  map<double,string>::reverse_iterator backwardIter;
  for (backwardIter = inverseMap.rbegin(); backwardIter != inverseMap.rend(); backwardIter++)
  {
    string name = (*backwardIter).second;
    double time = (*backwardIter).first;

    name = name + string("                                   ");
    name = name.substr(0,30);

    cout << "[" << time / totalTime * 100.0 << "%\t]: "
         << name.c_str() << " " << hours(time) << ":" << minutes(time) << ":" << seconds(time) << "s" << endl;
  }
  cout << "===========================================================================" << endl;
  cout << " total time: " << hours(totalTime) << ":" << minutes(totalTime) << ":" << seconds(totalTime)  << " (" << totalTime << " seconds) " << endl;
  cout << "===========================================================================" << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TIMER::printTimingsPerFrame(const int frames)
{
  timeval now;
  gettimeofday(&now, 0);

  double currentTimer = timing(_tick, now);
  string currentName;
  if (_callStack.size() > 0) 
    currentName = _callStack.top();

  // create an inverse map so that it will sort by time
  map<double, string> inverseMap;
  map<string, double>::iterator forwardIter;
  double totalTime = 0.0;
  for (forwardIter = _timings.begin(); forwardIter != _timings.end(); forwardIter++)
  {
    string name = forwardIter->first;
    double time = forwardIter->second;

    // if we are dinside this function, add the timing
    if (name.compare(currentName) == 0)
      time += currentTimer;

    inverseMap[time] = name;
    totalTime += time;
  }

  // print the map out backwards since it sorts from least to greatest
  cout << "====================================================================================" << endl;
  cout << " TIMING BREAKDOWN, FRAME " << frames << ": " << endl;
  cout << "====================================================================================" << endl;
  map<double,string>::reverse_iterator backwardIter;
  char buffer[256];
  string timeString;
  string hoursString;
  string minutesString;
  string secondsString;
  for (backwardIter = inverseMap.rbegin(); backwardIter != inverseMap.rend(); backwardIter++)
  {
    string name = (*backwardIter).second;
    double time = (*backwardIter).first;

    string padding("                                   ");

    name = name + padding;
    name = name.substr(0,30);

    sprintf(buffer, "%f", time / totalTime * 100.0);
    timeString = string(buffer);
    timeString = timeString + padding;
    timeString = timeString.substr(0,10);

    sprintf(buffer, "%02i", hours(time));
    hoursString = string(buffer);
    
    sprintf(buffer, "%02i", minutes(time));
    minutesString = string(buffer);

    sprintf(buffer, "%02i", seconds(time));
    secondsString = string(buffer);

    cout << "[" << timeString.c_str() << "%]: "
         << name.c_str() << " " << hoursString.c_str() << ":" 
                         << minutesString << ":"
                         << secondsString << " s total ";

    double perFrame = time / frames;
    sprintf(buffer, "%02i", hours(perFrame));
    hoursString = string(buffer);
    
    sprintf(buffer, "%02i", minutes(perFrame));
    minutesString = string(buffer);

    sprintf(buffer, "%02i", seconds(perFrame));
    secondsString = string(buffer);
    cout << " " << hoursString.c_str() << ":" << minutesString.c_str() << ":" << secondsString.c_str() << " s/frame" << endl;
  }
  sprintf(buffer, "%02i", hours(totalTime));
  hoursString = string(buffer);
    
  sprintf(buffer, "%02i", minutes(totalTime));
  minutesString = string(buffer);

  sprintf(buffer, "%02i", seconds(totalTime));
  secondsString = string(buffer);
  cout << "====================================================================================" << endl;
  cout << " Total time: " << hoursString.c_str() << ":" << minutesString.c_str() << ":" << secondsString.c_str() << " (" << totalTime << "s) " << endl;
  cout << " Time per frame: " << totalTime / frames << endl;
  cout << "====================================================================================" << endl;
}

