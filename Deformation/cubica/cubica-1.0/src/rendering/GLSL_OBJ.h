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
// GLSL_OBJ.h: Interface for the GLSL_OBJ class.
//
///////////////////////////////////////

#ifndef __GLSL_OBJ_H__
#define __GLSL_OBJ_H__

#if _WIN32
#include <gl/glew.h>
#include <gl/glut.h>
#elif USING_OSX
#include <GL/glew.h>
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#endif

#include <OBJ.h>
#include <GLSL_SHADER.h>

class GLSL_OBJ
{
  public:
    GLSL_OBJ(OBJ *obj);
    ~GLSL_OBJ();

    void update();

    void renderWithShadows();
    void renderNoShadows();
    void renderSM();
    void renderBasic();
    void renderDeferred();
    void renderAOV(float delta=0.1f, float width=2.0f, float height=2.0f);

    void setColor(GLfloat r, GLfloat g, GLfloat b);
    void setPosition(GLfloat x, GLfloat y, GLfloat z);

  protected:
    OBJ *_obj;

    int _numVertices;
    GLfloat *_vertex;
    GLfloat *_normal;

    GLfloat _color[3];
    GLfloat _position[3];

    GLSL_SHADER *_renderWithShadowShader;
    GLSL_SHADER *_renderNoShadowShader;
    GLSL_SHADER *_renderSMShader;
    GLSL_SHADER *_renderBasicShader;
    GLSL_SHADER *_renderDeferredShader;
    GLSL_SHADER *_renderAOVShader;

    void render(bool normals);   
};

#endif

