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
#include "SETTINGS.h"
#include "VECTOR.h"

#if USING_MKL
#include "VECTOR_FAST.cpp"
#elif USING_OSX
#include "VECTOR_FAST.cpp"
#elif __linux__
#include "VECTOR_FAST.cpp"
#else
#include "VECTOR_DEBUG.cpp"
#endif

//////////////////////////////////////////////////////////////////////
// clamp entries smaller than a threshold to zero
//////////////////////////////////////////////////////////////////////
void VECTOR::clampToZero(const Real threshold)
{
  for (int x = 0; x < _size; x++)
    if (std::fabs(_vector[x]) < threshold)
      _vector[x] = 0.0;
}
