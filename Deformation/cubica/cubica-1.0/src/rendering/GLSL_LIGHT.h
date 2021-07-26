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
// GLSL_LIGHT.h : 
/////////////////////////////

#ifndef GLSL_LIGHT_H
#define GLSL_LIGHT_H

#if _WIN32
#include <gl/glew.h>
//#include <gl/glut.h>
#elif USING_OSX
#include <GL/glew.h>
//#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#endif

using namespace std;

//////////////////////////////////////////////////////////////////////
// Light object to hold OpenGL/GLSL light information 
//////////////////////////////////////////////////////////////////////
class GLSL_LIGHT {

public:
  enum GLSL_LIGHT_TYPE 
  {
    DIRECTIONAL_LIGHT,
    POINT_LIGHT, 
    SPOT_LIGHT
  };

  GLSL_LIGHT();
  virtual ~GLSL_LIGHT();

  // accessors
  GLfloat* position()       { return _position; };
  GLfloat* diffuse()        { return _diffuse;  };
  GLfloat* ambient()        { return _ambient;  };
  GLfloat* specular()       { return _specular; };

  GLfloat attConst()        { return _att0; };
  GLfloat attLin()          { return _att1; };
  GLfloat attQuad()         { return _att2; };

  GLfloat* spotDirection()  { return _spotDirection; };
  GLfloat spotExponent()    { return _spotExponent;  };
  GLfloat spotCutoff()      { return _spotCutoff;    };

  GLSL_LIGHT_TYPE type()    { return _type; };

  GLuint shadowMapID()      { return _shadowMapID;  };
  GLsizei shadowWidth()     { return _shadowWidth;  };
  GLsizei shadowHeight()    { return _shadowHeight; };

  // queries
  bool isActive()           { return _active; };
  bool isShadowed()         { return _shadow; };

  // modifiers
  void setPosition(GLfloat *position);
  void setPosition(GLfloat x, GLfloat y, GLfloat z);

  void setDiffuse(GLfloat *diffuse);
  void setDiffuse(GLfloat r, GLfloat g, GLfloat b);

  void setAmbient(GLfloat *ambient);
  void setAmbient(GLfloat r, GLfloat g, GLfloat b);

  void setSpecular(GLfloat *specular);
  void setSpecular(GLfloat r, GLfloat g, GLfloat b);

  void setAttConst(GLfloat attenuation);
  void setAttLin(GLfloat attenuation);
  void setAttQuad(GLfloat attenuation);

  void setSpotDirection(GLfloat *direction);
  void setSpotDirection(GLfloat x, GLfloat y, GLfloat z);

  // Alternate method of setting spot direction, by giving a lookat position
  void setSpotLookat(GLfloat *look);
  void setSpotLookat(GLfloat x, GLfloat y, GLfloat z);

  void setSpotExponent(GLfloat exponent);
  void setSpotCutoff(GLfloat cutoff);
  void setSpotFov(GLdouble fov);
  
  void setType(GLSL_LIGHT_TYPE type);

  void setActive(bool active = true);
  void setShadowed(bool shadow);

  void setShadowDimensions(GLsizei width, GLsizei height);

  void computeMatrices();
  void viewTransformation();
  void textureTransformation(GLdouble *inverse_modelview);
  void textureTransformation();

protected:

  // Member variables:

  // What type of light is this
  GLSL_LIGHT_TYPE _type;

  // Is this light used to illuminate the scene
  bool _active;

  // Position of light in world space if _type == POINT_LIGHT|SPOT_LIGHT, 
  // or direction vector if _type == DIRECTIONAL_LIGHT (should be normalized)
  GLfloat _position[4];
  
  // RGB color triples for different illumination types
  GLfloat _diffuse[4];
  GLfloat _ambient[4];
  GLfloat _specular[4];

  // Constant, linear, and quadratic attentuation factors
  GLfloat _att0;
  GLfloat _att1;
  GLfloat _att2;

  // Spot light parameters. Ignored if _type != SPOT_LIGHT
  GLfloat _spotDirection[3];
  GLfloat _spotExponent;
  GLfloat _spotCutoff;

  GLdouble _spotFov;
  GLdouble _spotNear;
  GLdouble _spotFar;

  // Shadow map texture allocation is lazy, and won't create the texture until shadow mapping is
  // enabled for this light.

  // Does this light cast shadows? If _type == POINT_LIGHT, _shadow is set to be false
  bool _shadow;

  // OpenGL Shadow map texture object ID
  GLuint _shadowMapID;
  // Dimensions of the shadow map (should be power-of-two, and equal)
  GLsizei _shadowWidth;
  GLsizei _shadowHeight;

  // Viewing matrices, used for shadow mapping
  GLdouble _projectMatrix[16];
  GLdouble _modelviewMatrix[16];

  // Protected methods:

  // Allocates memory on the GPU for a shadow map using the given dimensions, 
  // storing the texture object ID in _shadowMapID
  void allocateShadowMap();
  void freeShadowMap();
};

#endif

