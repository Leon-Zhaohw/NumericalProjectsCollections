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

#include "BOX.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


BOX::BOX(const VEC3F& center, const VEC3F& lengths, double period, double translationPeriod) :
  _center(center), _lengths(lengths), _period(period), _translationPeriod(translationPeriod)
{
  set_halfLengths();
  set_angularVelocity();
  initialize_rotationMatrix();
  _theta = 0.0;
  _phi = 0.0;
  _currentTime = 0.0;
  _dt = 0.0;
  _rotationAxis = VEC3F(0.0, 0.0, 1.0);
  _displacement = VEC3F(0.0, 0.0, 0.0);
  _displacementVelocity = VEC3F(0.0, 0.0, 0.0);
  _original_center = _center;
}

BOX::BOX() 
{
}

BOX::~BOX()
{
}

void BOX::initialize_rotationMatrix()
{
  MATRIX3 eye = MATRIX3::I();
  _rotation = eye;
}

void BOX::update_rotationMatrix()
{
  // we must rotate *clockwise*, hence the minus sign
  _rotation = MATRIX3::rotation(_rotationAxis, -1 * _theta);
}

bool BOX::inside(const VEC3F& point) const
{
  // translate so that the center of rotation is the origin
  VEC3F pointTransformed = point - _center;
  // rotate *clockwise* back to the original position
  pointTransformed = _rotation * pointTransformed;

  return ( fabs(pointTransformed[0]) <= _halfLengths[0]  &&
           fabs(pointTransformed[1]) <= _halfLengths[1]  &&
           fabs(pointTransformed[2]) <= _halfLengths[2] );

}

void BOX::update_r(const VEC3F& point, VEC3F* r)
{
  VEC3F u = point - _center;
  // we want to project u onto the plane z = z0
  // this plane has unit normal n = (0, 0, 1)
  // so we project u onto this unit normal and then subtract off that component
  VEC3F normalComponent = project_onto(u, _rotationAxis);
  *r = u - normalComponent;
}

void BOX::calculate_displacementVelocity()
{
  // The position equation for the displacement in the x coordinate is:
  // x(t) = 1/4 * sin(2pi t / T), where T is the translation period.
  // Hence, the velocity v(t) is given by x'(t) as follows:
  // v(t) = 1/4 * 2pi / T * cos(2pi t / T)
  // Since phi will be updated according to phi_t = 2pi t / T, we simply use phi.
  _displacementVelocity[0] = 0.25 * 2 * M_PI / _translationPeriod * cos(_phi);
}


void BOX::spin() 
{
  double d_theta = 2 * M_PI * _dt / _period;
  _theta += d_theta;
  // printf("Theta is updated to: %f\n", _theta);
  // work modulo 2pi
  if (_theta >= 2 * M_PI) { _theta = 0.0; };
  this->update_rotationMatrix();
}

void BOX::translate_center()
{
  double d_phi= 2 * M_PI * _dt / _translationPeriod;
  _phi += d_phi;
  if (_phi >= 2 * M_PI) { _phi = 0.0; };
  // move along a sinusoidally controlled horizontal path of length 0.5,
  // centered at _original_center[0]
  _displacement[0] = 0.25 * sin(_phi);
  _center[0] = _original_center[0] + _displacement[0];
}


void BOX::draw() const
{
  glPushMatrix();
    glTranslatef(_center[0], _center[1], _center[2]);
    // OpenGL uses degrees, so we need to convert
    glRotatef(this->thetaDegrees(), _rotationAxis[0], _rotationAxis[1], _rotationAxis[2]);
    glScalef(_lengths[0], _lengths[1], _lengths[2]);
    glColor4b(30, 10, 5, 40);
    glutSolidCube(1.0);
  glPopMatrix();  
}


