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

#include <math.h>

#include <GL/glut.h>
#include "2D/VFIELD2D.h"

//set every thing on border to zero
void VFIELD2D::set_zero_border() {
  int totals = xRes*yRes;

  for(int i = 0; i < xRes; i++) {
    // Bottom.
    data[i] = 0;
    // Top.
    data[i+ totals  - xRes ] = 0;
  }
  
  for(int j = 0; j < yRes; j++) {
    // Left.
    data[j*xRes] = 0;
    // Right.
    data[xRes - 1 + j*xRes] = 0.;
  }
}
// Set X to zero.
void VFIELD2D::SetZeroLeft() {
  for(int j = 0; j < yRes; j++) {
    // Left.
    data[j*xRes][0] = 0;
  }
}
// Set X to zero.
void VFIELD2D::SetZeroRight() {
  for(int j = 0; j < yRes; j++) {
    // Right.
    data[j*xRes - 1 + xRes][0] = 0;
  }
}
// Set Y to zero.
void VFIELD2D::SetZeroTop() {
  int total = xRes*yRes;
  for(int i = 0;i < xRes; i++) {
    // Top.
    data[i+total  -xRes][1] = 0.;
  }
}
// Set Y to zero.
void VFIELD2D::SetZeroBottom() {
  for(int i = 0; i < xRes; i++) {
    // Bottom.
    data[i][1] = 0;
  }
}

// Set neumann on the left.
void VFIELD2D::SetNeumannLeft() {
  for(int j = 0; j < yRes; j++) {
    // Left.
    data[j*xRes][0] = data[j*xRes + 2][0];
   }
}

void VFIELD2D::SetNeumannRight() {
  for(int j = 0; j < yRes; j++) {
    // Right.
    data[j*xRes - 1 + xRes][0] = data[j*xRes -1 + xRes - 2][0];
   }
}

void VFIELD2D::SetNeumannTop() {
  int total = xRes*yRes;
  for(int i = 0; i < xRes; i++) {
    // Top.
    data[i + total -xRes ][1] = data[i + total -xRes - 2*xRes][1];
  }
}

void VFIELD2D::SetNeumannBottom() {
  for(int i = 0; i < xRes; i++) {
    // Bottom.
    data[i][1] = data[i + 2*xRes][1];
  }
}

// Copy border to the left.
void VFIELD2D::CopyBorderLeft() {
  for(int j = 0; j < yRes; j++) {
    //Left.
    data[j*xRes][0] = data[j*xRes + 1][0];
  }
}

void VFIELD2D::CopyBorderRight() {
  for(int j = 0; j < yRes; j++) {
    //Right.
    data[j*xRes - 1 + xRes][0] = data[j*xRes -1 + xRes - 1][0];
  }
}

void VFIELD2D::CopyBorderTop() {
  int total = xRes*yRes;
  for(int i = 0; i < xRes; i++) {
    //Top.
    data[i + total -xRes ][1] = data[i + total -xRes - xRes][1];
  }
}

void VFIELD2D::CopyBorderBottom() {
  for(int i = 0; i < xRes; i++) {
    //Bottom.
    data[i][1] = data[i + xRes][1];
  }
}

void VFIELD2D::axpy(const Real alpha, const VFIELD2D& field) {
  for(int i=0;i<totalcell;i++)
    data[i] += alpha *  field[i];
}

void VFIELD2D::advect(const Real dt, const VFIELD2D& velocity_field,
                      const FIELD2D& oldField, FIELD2D& newField) {
  int yRes = velocity_field.getyRes();
  int xRes = velocity_field.getxRes();
  #pragma omp parallel for
  for(int y=0; y<yRes;y++) {
    for(int x=0; x<xRes;x++) {
      int index = x + y*xRes;
      bool outside_bc = false;
      const VEC2F velo = velocity_field[index];
      Real xTrace = x - velo[0]*dt;
      Real yTrace = y - velo[1]*dt;
      if (std::isnan(xTrace)) {
        xTrace = x;
      }
      if (std::isnan(yTrace)) {
        yTrace = y;
      }
      if (xTrace < 1.5 || xTrace > xRes - 2.5 || yTrace < 1.5 || yTrace > yRes - 2.5) {
        outside_bc = true;
      }
      xTrace = (xTrace < 1.5) ? 1.5 : xTrace;
      xTrace = (xTrace > xRes - 2.5) ? xRes - 2.5 : xTrace;
      yTrace = (yTrace < 1.5) ? 1.5 : yTrace;
      yTrace = (yTrace > yRes - 2.5) ? yRes - 2.5 : yTrace;
      const int x0 = (int)xTrace;
      const int x1 = x0 + 1;
      const int y0 = (int)yTrace;
      const int y1 = y0 + 1;
      //interpolate
      const Real s1 = xTrace - x0;
      const Real s0 = 1.0f - s1;
      const Real t1 = yTrace - y0;
      const Real t0 = 1.0f - t1;

      const int i000 = x0 + y0 * xRes;
      const int i010 = x0 + y1 * xRes;
      const int i100 = x1 + y0 * xRes;
      const int i110 = x1 + y1 * xRes;
//      if (!outside_bc) {
        newField[index] =  s0 * (t0 * oldField[i000] +
          t1 * oldField[i010]) +
          s1 * (t0 * oldField[i100] +
          t1 * oldField[i110]);
//      } else {
//        newField[index] = 0.;
//      }
    }
  }
}

void VFIELD2D::advect(const Real dt, const VFIELD2D& velocity_field,
                      const VFIELD2D& oldField, VFIELD2D& newField) {
  int yRes = velocity_field.getyRes();
  int xRes = velocity_field.getxRes();
  for(int y=0; y<yRes;y++) {
    for(int x=0; x<xRes;x++) {
      int index = x + y*xRes;
      const VEC2F velo = velocity_field[index];
      Real xTrace = x - velo[0]*dt;
      Real yTrace = y - velo[1]*dt;
      xTrace = (xTrace < 1.5) ? 1.5 : xTrace;
      xTrace = (xTrace > xRes - 2.5) ? xRes - 2.5 : xTrace;
      yTrace = (yTrace < 1.5) ? 1.5 : yTrace;
      yTrace = (yTrace > yRes - 2.5) ? yRes - 2.5 : yTrace;
      const int x0 = (int)xTrace;
      const int x1 = x0 + 1;
      const int y0 = (int)yTrace;
      const int y1 = y0 + 1;
      //interpolate
      const Real s1 = xTrace - x0;
      const Real s0 = 1.0f - s1;
      const Real t1 = yTrace - y0;
      const Real t0 = 1.0f - t1;

      const int i000 = x0 + y0 * xRes;
      const int i010 = x0 + y1 * xRes;
      const int i100 = x1 + y0 * xRes;
      const int i110 = x1 + y1 * xRes;

      newField[index] =  s0 * (t0 * oldField[i000] +
        t1 * oldField[i010]) +
        s1 * (t0 * oldField[i100] +
        t1 * oldField[i110]);
      }
  }
}
   
void VFIELD2D::advectMacCormack(const Real dt, const VFIELD2D& velocityGrid, FIELD2D& oldField, 
                              FIELD2D& newField, FIELD2D& temp1, FIELD2D& temp2) {
  FIELD2D& phiHatN  = temp1;
	FIELD2D& phiHatN1 = temp2;

	const int sx = oldField.getxRes();
	const int sy = oldField.getxRes();

	for (int x = 0; x < sx * sy; x++)
		phiHatN[x] = phiHatN1[x] = oldField[x];

	FIELD2D& phiN    = oldField;
	FIELD2D& phiN1   = newField;
  // phiHatN1 = A(phiN)
	advect(dt, velocityGrid, phiN, phiHatN1);

	// phiHatN = A^R(phiHatN1)
	advect(-1.0 * dt, velocityGrid, phiHatN1, phiHatN);

  // phiN1 = phiHatN1 + (phiN - phiHatN) / 2
	const int border = 0; 
		for (int y = border; y < sy - border; y++)
			for (int x = border; x < sx - border; x++) {
				int index = x + y * sx;
				phiN1[index] = phiHatN1[index] + (phiN[index] - phiHatN[index]) * 0.50f;
			}

  phiN1.copyBorderAll();
  // clamp any newly created extrema
	clampExtrema(dt, velocityGrid, oldField, newField);

	// if the error estimate was bad, revert to first order
	clampOutsideRays(dt, velocityGrid, oldField, phiHatN1, newField);
}

void VFIELD2D::clampExtrema(const Real dt, const VFIELD2D& velocityField, const FIELD2D& 
                  oldField, FIELD2D& newField) {
  const int xRes = velocityField.getxRes();
	const int yRes = velocityField.getyRes();
	// const int zRes = velocityField.zRes();
	// const int slabSize = velocityField.slabSize();
  #pragma omp parallel
  #pragma omp for schedule(static)
  for (int y = 1; y < yRes-1; y++)
    for (int x = 1; x < xRes-1; x++)
    {
      const int index = x + y * xRes;
				// backtrace
        const VEC2F velocity = velocityField[index];
        Real xTrace = x - dt * velocity[0];
        Real yTrace = y - dt * velocity[1];
    
				// clamp backtrace to grid boundaries
        xTrace = (xTrace < 0.5) ? 0.5 : xTrace;
        xTrace = (xTrace > xRes - 1.5) ? xRes - 1.5 : xTrace;
        yTrace = (yTrace < 0.5) ? 0.5 : yTrace;
        yTrace = (yTrace > yRes - 1.5) ? yRes - 1.5 : yTrace;
        
				// locate neighbors to interpolate
				const int x0 = (int)xTrace;
				const int x1 = x0 + 1;
				const int y0 = (int)yTrace;
				const int y1 = y0 + 1;

				const int i000 = x0 + y0 * xRes;
				const int i010 = x0 + y1 * xRes;
				const int i100 = x1 + y0 * xRes;
				const int i110 = x1 + y1 * xRes;
//				const int i001 = x0 + y0 * xRes;
//				const int i011 = x0 + y1 * xRes;
//				const int i101 = x1 + y0 * xRes;
//				const int i111 = x1 + y1 * xRes;

				Real minField = oldField[i000];
				Real maxField = oldField[i000];

				minField = (oldField[i010] < minField) ? oldField[i010] : minField;
				maxField = (oldField[i010] > maxField) ? oldField[i010] : maxField;

				minField = (oldField[i100] < minField) ? oldField[i100] : minField;
				maxField = (oldField[i100] > maxField) ? oldField[i100] : maxField;

				minField = (oldField[i110] < minField) ? oldField[i110] : minField;
				maxField = (oldField[i110] > maxField) ? oldField[i110] : maxField;

				/*minField = (oldField[i001] < minField) ? oldField[i001] : minField;
				maxField = (oldField[i001] > maxField) ? oldField[i001] : maxField;

				minField = (oldField[i011] < minField) ? oldField[i011] : minField;
				maxField = (oldField[i011] > maxField) ? oldField[i011] : maxField;

				minField = (oldField[i101] < minField) ? oldField[i101] : minField;
				maxField = (oldField[i101] > maxField) ? oldField[i101] : maxField;

				minField = (oldField[i111] < minField) ? oldField[i111] : minField;
				maxField = (oldField[i111] > maxField) ? oldField[i111] : maxField;
        */

				newField[index] = (newField[index] > maxField) ? maxField : newField[index];
				newField[index] = (newField[index] < minField) ? minField : newField[index];
    }
}

void VFIELD2D::clampOutsideRays(const Real dt, const VFIELD2D& velocityField, const FIELD2D&  
                             oldField, const FIELD2D& oldAdvection, FIELD2D& newField) {
  const int sx = velocityField.getxRes();
	const int sy = velocityField.getyRes();
//	const int sz = velocityField.zRes();
//	const int slabSize = velocityField.slabSize();
  #pragma omp parallel
  #pragma omp for schedule(static)
  for (int y = 1; y < sy - 1; y++)
    for (int x = 1; x < sx - 1; x++)
    {
        const int index = x + y * sx;
				// backtrace
        VEC2F velocity = velocityField[index];
        velocity *= dt;
				float xBackward = x + velocity[0];
				float yBackward = y + velocity[1];
				
				float xTrace    = x - velocity[0];
				float yTrace    = x - velocity[1];
				
				// see if it goes outside the boundaries
				bool hasObstacle = 
					(yTrace < 1.0f)    || (yTrace > sy - 2.0f) ||
					(xTrace < 1.0f)    || (xTrace > sx - 2.0f) ||
					(yBackward < 1.0f) || (yBackward > sy - 2.0f) ||
					(xBackward < 1.0f) || (xBackward > sx - 2.0f);

				// reuse old advection instead of doing another one...
				if(hasObstacle) { newField[index] = oldAdvection[index]; }
      
    }
}

void VFIELD2D::swapPointer( VFIELD2D& field) {
  VEC2F* temp = data;
  data = field.data;
  field.data = temp;
}

void VFIELD2D::clear() {
  for(int i=0;i<totalcell;i++)
    data[i] = 0;
}

void VFIELD2D::SaveVelocity(FILE* out) {
  fwrite(data, sizeof(VEC2F),totalcell,out); 
}

void VFIELD2D::Draw(float magnitute) {
  float x, y, h;
  int maxres = std::max(xRes, yRes);
  h = 1.0/(Real)maxres;;

  glColor3f ( 1.0f, 1.0f, 1.0f );
  glLineWidth ( 1.0f );

  glBegin ( GL_LINES );

  for (int j=0 ; j<yRes-1; j++ ) {
    y = (j-0.5f)*h;
    for (int i=0 ; i<xRes-1; i++ ) {
      x = (i-0.5f)*h;
      int index = i + j*xRes;
      {glVertex2f ( x, y );
       glVertex2f ( x+data[index][0]*magnitute, y+data[index][1]*magnitute);
      }
    }
  }
  glEnd ();
}

void VFIELD2D::Draw_X(float magnitute, int mode) {
  float x, y, h, d00, d01, d10, d11;
  h = 1.0 / std::max(xRes, yRes);
  glBegin ( GL_QUADS );
  for (int j=0 ; j<yRes-1; j++ )
  {
    y = (j-0.5f)*h;
    for (int i=0 ; i<xRes-1; i++ ) 
    {
      x = (i-0.5f)*h;
      int index = i + j*xRes;
      d00 = (data[index][mode]) * magnitute;
      d01 = (data[index + xRes][mode]) * magnitute;
      d10 = (data[index + 1][mode]) * magnitute;
      d11 = (data[index + 1 + xRes][mode]) * magnitute;
      if (d00 > 0.0)  
        glColor3f ( d00, 0.0, 0.0 ); 
      else  glColor3f ( 0.0, -d00, 0.0 );
      glVertex2f ( x, y );
      
      if (d10 > 0)
        glColor3f ( d10, 0.0, 0.0 );
      else glColor3f ( 0.0, -d10, 0.0 );
      glVertex2f ( x+h, y );
      
      if (d11 > 0)
        glColor3f ( d11, 0.0, 0.0 );
      else glColor3f (0.0, -d11, 0.0 );
        glVertex2f ( x+h, y+h );
      
      if (d01 > 0) 
        glColor3f ( d01, 0.0, 0.0 );
      else glColor3f ( 0.0, -d01, 0.0 );
        glVertex2f ( x, y+h );
    }
  }
  glEnd ();
  glColor3f(1.0, 0,0);glPointSize(3.0);
  glBegin(GL_POINTS);
  for (int i = 0; i < 6; i++) {
    glVertex2f(static_cast<float>(i)*0.2, 0);
  }
  glEnd();
}

double VFIELD2D::GetEnergy() {
  double result = 0;
  int index = 0;
  for (int j = 0; j < yRes; j++) {
    for (int i = 0; i < xRes; i++) {
      result += data[index][0] * data[index][0];
      result += data[index][1] * data[index][1];
      index++;
    }
  }
  return result;
}

double VFIELD2D::GetMaximumVelocity() {
  double result = -1.0;
  int index = 0;
  for (int j = 0; j < yRes; j++) {
    for (int i = 0; i < xRes; i++) {
      double temp = data[index][0] * data[index][0] + data[index][1] * data[index][1];
      if (temp > result) {
        result = temp;
      }
      index++;
    }
  }
  return sqrt(result);
}

VEC2 VFIELD2D::GetVelocity(const float pos_x, const float pos_y) const {


  int x0 = (int)pos_x;
  int x1 = x0 + 1;
  int y0 = (int)pos_y;
  int y1 = y0 + 1;
  
  // clamp everything
  x0 = (x0 < 0) ? 0 : x0;
  y0 = (y0 < 0) ? 0 : y0;
  x1 = (x1 < 0) ? 0 : x1;
  y1 = (y1 < 0) ? 0 : y1;
  
  x0 = (x0 > xRes - 1) ? xRes - 1 : x0;
  y0 = (y0 > yRes - 1) ? yRes - 1 : y0;
  x1 = (x1 > xRes - 1) ? xRes - 1 : x1;
  y1 = (y1 > yRes - 1) ? yRes - 1 : y1;
  
  //interpolate
  const Real s1 = pos_x - x0;
  const Real s0 = 1.0f - s1;
  const Real t1 = pos_y - y0;
  const Real t0 = 1.0f - t1;

  const int i000 = x0 + y0 * xRes;
  const int i010 = x0 + y1 * xRes;
  const int i100 = x1 + y0 * xRes;
  const int i110 = x1 + y1 * xRes;
  
  return s0 * (t0 * data[i000] +
    t1 * data[i010]) +
    s1 * (t0 * data[i100] +
    t1 * data[i110]);
}

double VFIELD2D::DotProduct(const VFIELD2D& b) {
  double result = 0;
  for (int i = 0; i < totalcell; i++) {
    result += data[i][0] * b[i][0] + data[i][1] * b[i][1];
  }
  
  return result;
}


