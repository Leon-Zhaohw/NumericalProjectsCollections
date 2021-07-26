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
// PRECONDITIONER.h: interface for the PRECONDITIONER class.
//
//////////////////////////////////////////////////////////////////////

#ifndef PRECONDITIONER_H
#define PRECONDITIONER_H

#include <SPARSE_MATRIX.h>

//////////////////////////////////////////////////////////////////////
// Preconditioner for SPARSE_PCG_MATRIX
//////////////////////////////////////////////////////////////////////
class PRECONDITIONER {

public:
  PRECONDITIONER()  {};
  virtual ~PRECONDITIONER() {};

  virtual void init() = 0;
  virtual void solve(VECTOR& x, VECTOR& b) = 0;
};

#endif
