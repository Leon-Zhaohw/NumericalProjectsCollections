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
 
 Copyright 2018 Qiaodong Cui (qiaodong@ucsb.edu)
 */

#include "2D/FIELD2D.h"
#include <stdio.h>

void FIELD2D::clear() {
  for(int i = 0; i < xRes*yRes; i++)
    data[i] = 0;
}

void FIELD2D::swapPointer(FIELD2D& field) {
  assert(xRes == field.getxRes());
  assert(yRes == field.getyRes());
  Real* tmp = NULL;
  tmp = data;
  data = field.data;
  field.data = tmp;
}

void FIELD2D::setZeroBorder() {
  int totals = xRes*yRes;

  for(int i = 0; i < xRes; i++) {
    data[i] = 0;
    data[i+ totals  - xRes ] = 0;
  }
  for(int j = 0; j < yRes; j++) {
    data[j*xRes] = 0;
    data[xRes - 1 + j*xRes] = 0.;
  }
}

void FIELD2D::setZeroBorder2Layers() {
  int totals = xRes*yRes;
  
  for(int i = 0; i < xRes; i++) {
    // Bottom.
    data[i] = 0;
    data[i + xRes] = 0;
    // Top.
    data[i + totals - xRes] = 0;
    data[i + totals - 2*xRes] = 0.;
  }
  for(int j = 0; j < yRes; j++) {
    // Left.
    data[j*xRes] = 0;
    data[j*xRes + 1] = 0;
    // Right.
    data[xRes - 1 + j*xRes - 1] = 0;
    data[xRes - 1 + j*xRes] = 0;
  }
}

void FIELD2D::copyBorderAll() {

  for(int i = 0; i < xRes; i++) {
    data[i] = data[i + xRes];
    data[i + totalsize -xRes] = data[i+totalsize -2*xRes];
  }
  
  for(int j = 0; j < yRes; j++) {
    data[j*xRes] = data[j*xRes + 1];
    data[j*xRes + xRes - 1] = data[j*xRes + xRes -2];
  }
}

Real FIELD2D::dot(const FIELD2D& in) {
  assert(xRes == in.getxRes());
  assert(yRes == in.getyRes());
  Real result = 0.;
  for(int i = 0; i < totalsize; i++)
  {
  	result += data[i]*in[i];
  }

  return result;
}
void FIELD2D::axpy(const Real& alpha, const FIELD2D& input) {
  assert(xRes == input.getxRes());
  assert(yRes == input.getyRes());
  for(int i = 0; i < totalsize; i++) {
    data[i] += alpha*input[i];
  }
}
// FIELD2D& FIELD2D::operator=(const Real& alpha)
// {
// 	for(int i=0;i<totalsize;i++)
// 		data[i] = alpha;
// 	return *this;
// }
// FIELD2D& FIELD2D::operator=(const FIELD2D& A)
// {
// 	for(int i=0;i<totalsize;i++)
// 		data[i] = A[i];
// 	return *this;
//}
FIELD2D& FIELD2D::operator *=(const Real alpha) {
  for (int x = 0; x < totalsize; x++)
    data[x] *= alpha;
  return *this;
}

FIELD2D& FIELD2D::operator +=(const Real alpha) {
  for(int x = 0; x < totalsize; x++)
    data[x] += alpha;

  return *this;
}

void FIELD2D::CheckNegative() {
  for (int i = 0; i < totalsize; i++) {
    if (data[i] < 0) {
      printf("Negative");
    }
  }
}

FIELD2D& FIELD2D::operator +=(const FIELD2D& input) {
  assert(xRes == input.getxRes());
  assert(yRes == input.getyRes());

  for(int x = 0; x < totalsize; x++)
    data[x] += input[x];

  return *this;
}

FIELD2D& FIELD2D::operator *=(const FIELD2D& input) {
  assert(xRes == input.getxRes());
  assert(yRes == input.getyRes());

  for(int x = 0; x < totalsize; x++)
    data[x] *= input[x];

  return *this;
}

void FIELD2D::Field2D_Save(FILE* out) {
  fwrite(data, sizeof(Real), totalsize, out);
}
