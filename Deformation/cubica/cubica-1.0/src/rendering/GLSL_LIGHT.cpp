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
// GLSL_LIGHT.cpp
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include <GLSL_LIGHT.h>
#include <iostream>

//////////////////////////////////////////////////////////////////////
// Constructor for GLSL_LIGHT
//////////////////////////////////////////////////////////////////////
GLSL_LIGHT::GLSL_LIGHT() :
  _type(DIRECTIONAL_LIGHT), _active(false),
  _att0(1.0f), _att1(0.0f), _att2(0.0f),  
  _spotExponent(1.0f), _spotCutoff(180.0f),
  _spotFov(45.0), _spotNear(1.0), _spotFar(100.0),
  _shadow(false), _shadowMapID(0), _shadowWidth(512), _shadowHeight(512)
{
  _position[0] = _position[2] = _position[3] = 0.0f;
  _position[1] = -1.0f;

  _diffuse[0] = _diffuse[1] = _diffuse[2] = 0.7f;
  _ambient[0] = _ambient[1] = _ambient[2] = 0.3f;
  _specular[0] = _specular[1] = _specular[2] = 0.0f;
  _diffuse[3] = _ambient[3] = _specular[3] = 1.0f;

  computeMatrices();
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
GLSL_LIGHT::~GLSL_LIGHT()
{
  freeShadowMap();
}

//////////////////////////////////////////////////////////////////////
// Set light's position/direction
//////////////////////////////////////////////////////////////////////
void GLSL_LIGHT::setPosition(GLfloat *position)
{
  this->setPosition(position[0], position[1], position[2]);
}

void GLSL_LIGHT::setPosition(GLfloat x, GLfloat y, GLfloat z)
{
  _position[0] = x;
  _position[1] = y;
  _position[2] = z;
  _position[3] = 1.0f;

  // Directional lights use _position to store normalized direction
  if (_type == DIRECTIONAL_LIGHT)
  {
    GLfloat len = _position[0]*_position[0] 
                + _position[1]*_position[1] 
                + _position[2]*_position[2];
    if (len > 0.0f)
    {
      // Normalize direction if vector isn't 0
      len = sqrt(len);
      _position[0] = _position[0] / len;
      _position[1] = _position[1] / len;
      _position[2] = _position[2] / len;
    }
    else 
    {
      // If vector is 0, set to look down y axis
      _position[0] = 0.0f;
      _position[1] =-1.0f;
      _position[2] = 0.0f;
    }
    _position[3] = 0.0f;
  }

  computeMatrices();
}

//////////////////////////////////////////////////////////////////////
// Set light's diffuse color
//////////////////////////////////////////////////////////////////////
void GLSL_LIGHT::setDiffuse(GLfloat *diffuse)
{ this->setDiffuse(diffuse[0], diffuse[1], diffuse[2]);   }

void GLSL_LIGHT::setDiffuse(GLfloat r, GLfloat g, GLfloat b)
{
  _diffuse[0] = r;
  _diffuse[1] = g;
  _diffuse[2] = b;
  _diffuse[3] = 1.0f;
}

//////////////////////////////////////////////////////////////////////
// Set light's ambient color
void GLSL_LIGHT::setAmbient(GLfloat *ambient)
{ this->setAmbient(ambient[0], ambient[1], ambient[2]);   }

void GLSL_LIGHT::setAmbient(GLfloat r, GLfloat g, GLfloat b)
{
  _ambient[0] = r;
  _ambient[1] = g;
  _ambient[2] = b;
  _ambient[3] = 1.0f;
}

//////////////////////////////////////////////////////////////////////
// Set light's specular color
//////////////////////////////////////////////////////////////////////
void GLSL_LIGHT::setSpecular(GLfloat *specular)
{ this->setSpecular(specular[0], specular[1], specular[2]);   }

void GLSL_LIGHT::setSpecular(GLfloat r, GLfloat g, GLfloat b)
{
  _specular[0] = r;
  _specular[1] = g;
  _specular[2] = b;
  _specular[3] = 1.0f;
}

//////////////////////////////////////////////////////////////////////
// Set light's attentuation coefficients
//////////////////////////////////////////////////////////////////////

// Set light's constant distance attenuation coefficient
void GLSL_LIGHT::setAttConst(GLfloat attenuation)
{
  _att0 = attenuation;
}
// Set light's linear distance attenuation coefficient
void GLSL_LIGHT::setAttLin(GLfloat attenuation)
{
  _att1 = attenuation;
}
// Set light's quadratic distance attenuation coefficient
void GLSL_LIGHT::setAttQuad(GLfloat attenuation)
{
  _att2 = attenuation;
}

//////////////////////////////////////////////////////////////////////
// Set light's spot direction given a vector
//////////////////////////////////////////////////////////////////////
void GLSL_LIGHT::setSpotDirection(GLfloat *direction)
{ this->setSpotDirection(direction[0], direction[1], direction[2]);   }

void GLSL_LIGHT::setSpotDirection(GLfloat x, GLfloat y, GLfloat z)
{
  _spotDirection[0] = x;
  _spotDirection[1] = y;
  _spotDirection[2] = z;

  GLfloat len = _spotDirection[0]*_spotDirection[0] 
              + _spotDirection[1]*_spotDirection[1] 
              + _spotDirection[2]*_spotDirection[2];
  if (len > 0.0f)
  {
    // If direction vector is not 0, normalize it
    len = sqrt(len);
    _spotDirection[0] = _spotDirection[0] / len;
    _spotDirection[1] = _spotDirection[1] / len;
    _spotDirection[2] = _spotDirection[2] / len;
  }
  else
  {
    // If direction vector provided is 0, set direction to look down y axis
    _spotDirection[0] = 0.0f;
    _spotDirection[1] =-1.0f;
    _spotDirection[2] = 0.0f;
  }

  computeMatrices();
}

//////////////////////////////////////////////////////////////////////
// Set light's spot direction given a lookat position
//////////////////////////////////////////////////////////////////////
void GLSL_LIGHT::setSpotLookat(GLfloat *look)
{ this->setSpotLookat(look[0], look[1], look[2]);   }

void GLSL_LIGHT::setSpotLookat(GLfloat x, GLfloat y, GLfloat z)
{
  GLfloat dx = x - _position[0];
  GLfloat dy = y - _position[1];
  GLfloat dz = z - _position[2];

  setSpotDirection(dx, dy, dz);
}

//////////////////////////////////////////////////////////////////////
// Set light's spot exponent
//////////////////////////////////////////////////////////////////////
void GLSL_LIGHT::setSpotExponent(GLfloat exponent)
{
  _spotExponent = exponent;
}

//////////////////////////////////////////////////////////////////////
// Set light's spot field of view (zoom level), in degrees
//////////////////////////////////////////////////////////////////////
void GLSL_LIGHT::setSpotFov(GLdouble fov)
{
  if (fov == _spotFov)
    return;

  if (fov < 5.0)
  {
    _spotFov = 5.0;
  }
  else if (fov > 180.0)
  {
    _spotFov = 180.0;
  }
  else
  {
    _spotFov = fov;
  }

  computeMatrices();
}

//////////////////////////////////////////////////////////////////////
// Set light's spot cutoff (in degrees)
//////////////////////////////////////////////////////////////////////
void GLSL_LIGHT::setSpotCutoff(GLfloat cutoff)
{
  _spotCutoff = cutoff;
}

//////////////////////////////////////////////////////////////////////
// Set light's type
//////////////////////////////////////////////////////////////////////
void GLSL_LIGHT::setType(GLSL_LIGHT_TYPE type)
{
  _type = type;

  if (_type == DIRECTIONAL_LIGHT)
  {
    // Use setPosition to normalize _position for use as a direction
    setPosition(_position[0], _position[1], _position[2]);
  }
  else if (_type == POINT_LIGHT)
  {
    _shadow = false;
  }
  else if (_type == SPOT_LIGHT)
  {
    if (_spotDirection[0] == _spotDirection[1] == _spotDirection[2] == 0.0f)
    {
      // Make spot light look at the origin
      setSpotLookat(0.0f, 0.0f, 0.0f);
    }
  }
  else
  {
    // Incorrect type passed in. Set to default DIRECTIONAL_LIGHT.
    _type = DIRECTIONAL_LIGHT;

    // Use setPosition to normalize _position for use as a direction
    setPosition(_position[0], _position[1], _position[2]);
  }
}

//////////////////////////////////////////////////////////////////////
// Set whether light is active or not
//////////////////////////////////////////////////////////////////////
void GLSL_LIGHT::setActive(bool active)
{
  _active = active;
}

//////////////////////////////////////////////////////////////////////
// Set whether light uses shadow mapping or not
//////////////////////////////////////////////////////////////////////
void GLSL_LIGHT::setShadowed(bool shadow)
{
  // Point light would need environment shadow map, and I don't want to do that right now
  if (_type == POINT_LIGHT)
  {
    // Point light cannot use shadow mapping
    _shadow = false;
    return;
  }

  // Set shadow mapping flag
  _shadow = shadow;
  if (_shadow)
  {
    // If we want to use shadow mapping, check if a shadow map has been allocated
    if (_shadowMapID == 0)
    {
      // If it hasn't, make a new one
      allocateShadowMap();
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Set shadow map dimensions. They default to 512x512.
//////////////////////////////////////////////////////////////////////
void GLSL_LIGHT::setShadowDimensions(GLsizei width, GLsizei height)
{
  // If either dimension is set to 0, don't do anything
  if (width == 0 || height == 0)
  {
    return;
  }
  // If new dimensions are identical to old dimensions, don't do anything
  if (_shadowWidth == width && _shadowHeight == height)
  {
    return;
  }

  // Set new dimensions
  _shadowWidth = width;
  _shadowHeight = height;

  if (_shadowMapID != 0)
  {
    // Free old shadowmap
    freeShadowMap();
    // allocate new shadowmap
    allocateShadowMap();
  }
}

//////////////////////////////////////////////////////////////////////
// Calculate and store view matrices of the light
//////////////////////////////////////////////////////////////////////
void GLSL_LIGHT::computeMatrices()
{
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
  if (_type == SPOT_LIGHT)
  {
	  gluPerspective(_spotFov, 1.0, _spotNear, _spotFar);
	  glGetDoublev(GL_PROJECTION_MATRIX, _projectMatrix);
  }
  else if (_type == DIRECTIONAL_LIGHT)
  {
    glOrtho(-5.0, 5.0, -5.0, 5.0, -5.0, 5.0);
    glGetDoublev(GL_PROJECTION_MATRIX, _projectMatrix);
  }
	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
  if (_type == SPOT_LIGHT)
  {
  	gluLookAt(_position[0], _position[1], _position[2], 
              _position[0] + _spotDirection[0],
              _position[1] + _spotDirection[1],
              _position[2] + _spotDirection[2], 
              0.0, 1.0, 0.0);
    glGetDoublev(GL_MODELVIEW_MATRIX, _modelviewMatrix);
  }
  else if (_type == DIRECTIONAL_LIGHT)
  {
    gluLookAt(_position[0], _position[1], _position[2], 
      				0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    glGetDoublev(GL_MODELVIEW_MATRIX, _modelviewMatrix);
  }

	glPopMatrix();
}

//////////////////////////////////////////////////////////////////////
// Transform the view to the light's location
//////////////////////////////////////////////////////////////////////
void GLSL_LIGHT::viewTransformation()
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glMultMatrixd(_projectMatrix);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMultMatrixd(_modelviewMatrix);//*/
}

//////////////////////////////////////////////////////////////////////
// Store the ModelViewProjection matrix in the current GL_TEXTURE matrix
//////////////////////////////////////////////////////////////////////
void GLSL_LIGHT::textureTransformation(GLdouble *inverseModelview)
{
	glMatrixMode(GL_TEXTURE);
	glLoadIdentity();
	glTranslated(0.5, 0.5, 0.5);
	glScaled(0.5, 0.5, 0.5);		
	glMultMatrixd(_projectMatrix);
	glMultMatrixd(_modelviewMatrix);
	glMultMatrixd(inverseModelview);
	glMatrixMode(GL_MODELVIEW);
}

//////////////////////////////////////////////////////////////////////
// Store the ModelViewProjection matrix in the current GL_TEXTURE matrix
//////////////////////////////////////////////////////////////////////
void GLSL_LIGHT::textureTransformation()
{
	glMatrixMode(GL_TEXTURE);
	glLoadIdentity();
	glTranslated(0.5, 0.5, 0.5);
	glScaled(0.5, 0.5, 0.5);		
	glMultMatrixd(_projectMatrix);
	glMultMatrixd(_modelviewMatrix);
	glMatrixMode(GL_MODELVIEW);
}


//////////////////////////////////////////////////////////////////////
// Generate a shadow map texture object
//////////////////////////////////////////////////////////////////////
void GLSL_LIGHT::allocateShadowMap()
{
  if (_shadowMapID != 0)
  {
    cout << " GLSL_LIGHT trying to allocate shadow map texture when one already exists" << endl;
    freeShadowMap();
  }

  glGenTextures(1, &_shadowMapID);
  if (_shadowMapID == 0)
  {
    cout << " GLSL_LIGHT could not generate a shadow map texture" << endl;
    return;
  }
  
  glBindTexture(GL_TEXTURE_2D, _shadowMapID);
   
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F_ARB, _shadowWidth, _shadowHeight, 0, GL_RGB, GL_FLOAT, 0);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  glBindTexture(GL_TEXTURE_2D, 0);
}

//////////////////////////////////////////////////////////////////////
// Delete the shadow map texture object
//////////////////////////////////////////////////////////////////////
void GLSL_LIGHT::freeShadowMap()
{
  if (_shadowMapID != 0)
  {
    glDeleteTextures(1, &_shadowMapID);
  }
  _shadowMapID = 0;
}


