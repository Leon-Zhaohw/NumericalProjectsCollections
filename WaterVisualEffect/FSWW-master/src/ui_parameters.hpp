/* 
 * File: ui_parameters.hpp
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

#ifndef UI_PARAMETERS_HPP
#define UI_PARAMETERS_HPP

#include "definitions.hpp"

namespace ui_parameters {
  extern bool show_in_field;
  extern bool show_scattered_field;

  extern bool circle_;
  extern bool square_;
  extern bool island_;
  extern bool harbour_;
  extern bool line_;
  extern bool test_;
  extern bool wavy_;

  extern bool point_source1_;
  extern bool point_source2_;
  extern bool point_source3_;
  extern bool linear_wave1_;
  extern bool linear_wave2_;
  extern bool linear_wave3_;
  extern bool user_def_source_;

  extern VEC2 user_def_pos_;

  extern bool running_;
  extern FLOAT size_drop_;
};

#endif
