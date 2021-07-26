/* 
 * File: Times.hpp
 *
 * Copyright (C) 2019  Camille Schreck
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef TIMES_HPP
#define TIMES_HPP

#include <boost/date_time/posix_time/posix_time.hpp>

class Times {
  
public :
  static Times *TIMES;
  static Times *TIMES_UP;

  enum timing_t {simu_time_ = 0, 
		 display_time_,
		 solve_time_,
		 sum_up_time_,
		 total_time_,
		 nTimes};

private :
  boost::posix_time::ptime init_time[nTimes];
  double loop_time[nTimes];
  double time_sum[nTimes];

  //count
  unsigned int frame;
  unsigned int step;

public :
  Times();
  
  void init();

  void tick(unsigned int i);
  void tock(unsigned int i);
  double getTime(unsigned int i);
  double getAverageTime(unsigned int i);
  double getAverageTimeByStep(unsigned int i);
  double getAverageNbStepsByFrame();
  
  void next_loop();
  void next_step();
};

#endif //TIMES_HPP
