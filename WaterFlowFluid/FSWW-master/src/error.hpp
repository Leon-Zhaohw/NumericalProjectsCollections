/* 
 * File: error.hpp
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

#ifndef ERROR_HPP
#define ERROR_HPP

#include <assert.h>
#include <iostream>
#include "ExceptionSimu.hpp"

//****************************************************

#ifndef __MODE_DEBUG  
#define __MODE_DEBUG 3
#endif

#if __MODE_DEBUG > 0

#include <assert.h>

void close();

#define ERROR(cond, msg, debug_info) if (!(cond)) {std::cerr<<"ERROR: "<<msg<<"\n"<<debug_info<<std::endl;throw ExceptionSimu();}

#define WARNING(cond, msg, debug_info) if (!(cond)) {std::cout<<"WARNING: "<<msg<<"\n"<<debug_info<<std::endl;}

#define VERBOSE(level, debug_info) if (__MODE_DEBUG >= level) {std::cout<<"INFO: \n"<<debug_info<<std::endl;}

#define INFO(debug_info) if (__MODE_DEBUG >= 3) {std::cout<<"INFO: \n"<<debug_info<<std::endl;}

#define TEST(cond) if (!(cond)) {close(); assert(cond);}

#define IS_DEF(nombre) assert(!std::isnan(nombre) && !std::isinf(nombre))

//*********************

#else

#define ERROR(cond, msg, debug_info) do {if (!(cond)) {std::cerr<<"\n***ERREUR*** : \n"<<msg<<"\n"<<std::endl; std::exit(-1)}}while(0)

#define WARNING(cond, msg, debug_info) do {if (!(cond)) {std::cerr<<"\n****WARNING*** : \n"<<msg<<"\n"<<std::endl;}while(0)

#define VERBOSE(level, debug_info) if (__MODE_DEBUG >= level) {std::cout<<"INFO: \n"<<debug_info<<std::endl;}

#define INFO(debug_info) 

#define TEST(cond)

#define IS_DEF(nombre)


#endif

//****************************************************


#endif
