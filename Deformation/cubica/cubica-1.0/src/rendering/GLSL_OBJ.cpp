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
// GLSL_OBJ.cpp
//
/////////////////////////////////////////////////////////////////////

#include <GLSL_OBJ.h>
#include <stdio.h>

/////////////////////////////////////////////////////////////////////
// Constructor
/////////////////////////////////////////////////////////////////////
GLSL_OBJ::GLSL_OBJ(OBJ *obj)
{
  _obj = obj;
  if (_obj == NULL)
  {
    cout << "GLSL_OBJ was given a NULL obj pointer! Bad things will happen..." << endl;
    return;
  }

  _numVertices = _obj->faces.size() * 3;
  
  _vertex = new GLfloat[_numVertices * 3];
  if (_vertex == NULL) {
    cout << "GLSL_OBJ could not allocate memory for vertices." << endl;
  }

  _normal = new GLfloat[_numVertices * 3];
  if (_normal == NULL) {
    cout << "GLSL_OBJ could not allocate memory for normals." << endl;
  }

  update();

  _color[0] = _color[1] = _color[2] = 1.0f;
  _position[0] = _position[1] = _position[2] = 0.0f;

  _renderWithShadowShader = new GLSL_SHADER();
  if (_renderWithShadowShader != NULL)
  {
    _renderWithShadowShader->attachVert("src/rendering/glsl/drawobj00.vert");
    _renderWithShadowShader->attachFrag("src/rendering/glsl/drawobj00.frag");
    _renderWithShadowShader->compile();
  }

  _renderNoShadowShader = new GLSL_SHADER();
  if (_renderNoShadowShader != NULL)
  {
    _renderNoShadowShader->attachVert("src/rendering/glsl/drawobj00.vert");
    _renderNoShadowShader->attachFrag("src/rendering/glsl/drawobj01.frag");
    _renderNoShadowShader->compile();
  }

  _renderSMShader = new GLSL_SHADER();
  if (_renderSMShader != NULL)
  {
    _renderSMShader->attachVert("src/rendering/glsl/drawobj02.vert");
    _renderSMShader->attachFrag("src/rendering/glsl/drawobj02.frag");
    _renderSMShader->compile();
  }

  _renderBasicShader = new GLSL_SHADER();
  if (_renderBasicShader != NULL)
  {
    _renderBasicShader->attachVert("src/rendering/glsl/drawobj05.vert");
    _renderBasicShader->attachFrag("src/rendering/glsl/drawobj05.frag");
    _renderBasicShader->compile();
  }

  _renderDeferredShader = new GLSL_SHADER();
  if (_renderDeferredShader != NULL)
  {
    _renderDeferredShader->attachVert("src/rendering/glsl/drawobj03.vert");
    _renderDeferredShader->attachFrag("src/rendering/glsl/drawobj03.frag");
    _renderDeferredShader->compile();
  }

  _renderAOVShader = new GLSL_SHADER();
  if (_renderAOVShader != NULL)
  {
    _renderAOVShader->attachVert("src/rendering/glsl/drawobj04.vert");
    _renderAOVShader->attachGeom("src/rendering/glsl/drawobj04.geom");
    _renderAOVShader->attachFrag("src/rendering/glsl/drawobj04.frag");
    _renderAOVShader->setGeomInput(GL_TRIANGLES);
    _renderAOVShader->setGeomOutput(GL_TRIANGLE_STRIP);
    _renderAOVShader->setGeomVerticesOut(16);
    _renderAOVShader->compile();
  }
}  
  
/////////////////////////////////////////////////////////////////////
// Destructor
/////////////////////////////////////////////////////////////////////
GLSL_OBJ::~GLSL_OBJ()
{
  // Clean up allocated vertex arrays
  if (_vertex != NULL)   {
    delete [] _vertex;
    _vertex = NULL;
  }
  if (_normal != NULL)   {
    delete [] _normal;
    _normal = NULL;
  }

  // Clean up all the shader objects
  if (_renderWithShadowShader != NULL)   {
    delete _renderWithShadowShader;
    _renderWithShadowShader = NULL;
  }
  if (_renderNoShadowShader != NULL)   {
    delete _renderNoShadowShader;
    _renderNoShadowShader = NULL;
  }
  if (_renderSMShader != NULL)   {
    delete _renderSMShader;
    _renderSMShader = NULL;
  }
  if (_renderBasicShader != NULL) {
    delete _renderBasicShader;
    _renderBasicShader = NULL;
  }
  if (_renderDeferredShader != NULL)   {
    delete _renderDeferredShader;
    _renderDeferredShader = NULL;
  }
  if (_renderAOVShader != NULL)   {
    delete _renderAOVShader;
    _renderAOVShader = NULL;
  }
}

/////////////////////////////////////////////////////////////////////
// Set OBJ origin's position
/////////////////////////////////////////////////////////////////////
void GLSL_OBJ::setPosition(GLfloat x, GLfloat y, GLfloat z)
{
  _position[0] = x;
  _position[1] = y;
  _position[2] = z;
}

/////////////////////////////////////////////////////////////////////
// Set OBJ color
/////////////////////////////////////////////////////////////////////
void GLSL_OBJ::setColor(GLfloat r, GLfloat g, GLfloat b)
{
  _color[0] = r;
  _color[1] = g;
  _color[2] = b;
}

/////////////////////////////////////////////////////////////////////
// Update function to copy latest obj geometry to local vertex arrays
/////////////////////////////////////////////////////////////////////
void GLSL_OBJ::update()
{
  vector<OBJ::Face>& faces = _obj->faces;
  int numFaces = faces.size();

  // May want to recompute normals here.
  _obj->ComputeVertexNormals();

  VEC3F v;
  VEC3F n;
  for (int i=0; i<numFaces; i++)
  {
    v = _obj->vertices[faces[i].vertices[0]];
    _vertex[i*9]   = (GLfloat)v[0];
    _vertex[i*9+1] = (GLfloat)v[1];
    _vertex[i*9+2] = (GLfloat)v[2];
    //cout << " (" << v[0] << "," << v[1] << "," << v[2] << ")   ";
    v = _obj->vertices[faces[i].vertices[1]];
    _vertex[i*9+3] = (GLfloat)v[0];
    _vertex[i*9+4] = (GLfloat)v[1];
    _vertex[i*9+5] = (GLfloat)v[2];
    //cout << " (" << v[0] << "," << v[1] << "," << v[2] << ")   ";
    v = _obj->vertices[faces[i].vertices[2]];
    _vertex[i*9+6] = (GLfloat)v[0];
    _vertex[i*9+7] = (GLfloat)v[1];
    _vertex[i*9+8] = (GLfloat)v[2];
    //cout << " (" << v[0] << "," << v[1] << "," << v[2] << ")" << endl;
    
    n = _obj->normals[faces[i].vertices[0]];
    _normal[i*9]   = (GLfloat)n[0];
    _normal[i*9+1] = (GLfloat)n[1];
    _normal[i*9+2] = (GLfloat)n[2];
    n = _obj->normals[faces[i].vertices[1]];
    _normal[i*9+3] = (GLfloat)n[0];
    _normal[i*9+4] = (GLfloat)n[1];
    _normal[i*9+5] = (GLfloat)n[2];
    n = _obj->normals[faces[i].vertices[2]];
    _normal[i*9+6] = (GLfloat)n[0];
    _normal[i*9+7] = (GLfloat)n[1];
    _normal[i*9+8] = (GLfloat)n[2];
  }
}

/////////////////////////////////////////////////////////////////////
// Render OBJ with shadow mapping
/////////////////////////////////////////////////////////////////////
void GLSL_OBJ::renderWithShadows()
{
  glUseProgram(_renderWithShadowShader->getHandle());
  glUniform1i(_renderWithShadowShader->uniformLoc("shadowTex"), 4);
  // Also has Uniform1f for minVariance and reduceBleed.

  glColor3f(_color[0], _color[1], _color[2]);
  render(true);

  glUseProgram(0);
}

/////////////////////////////////////////////////////////////////////
// Render OBJ without shadows
/////////////////////////////////////////////////////////////////////
void GLSL_OBJ::renderNoShadows()
{
  glUseProgram(_renderNoShadowShader->getHandle());

  glColor3f(_color[0], _color[1], _color[2]);
  render(true);

  glUseProgram(0);
}

/////////////////////////////////////////////////////////////////////
// Render OBJ into a shadow map
/////////////////////////////////////////////////////////////////////
void GLSL_OBJ::renderSM()
{
  glUseProgram(_renderSMShader->getHandle());

  render(false);

  glUseProgram(0);
}

/////////////////////////////////////////////////////////////////////
// Render OBJ with current color
/////////////////////////////////////////////////////////////////////
void GLSL_OBJ::renderBasic()
{
  glUseProgram(_renderBasicShader->getHandle());

  glColor3f(_color[0], _color[1], _color[2]);
  render(false);

  glUseProgram(0);
}

/////////////////////////////////////////////////////////////////////
// Render OBJ to deferred render targets
/////////////////////////////////////////////////////////////////////
void GLSL_OBJ::renderDeferred()
{
  glUseProgram(_renderDeferredShader->getHandle());

  render(true);

  glUseProgram(0);
}

/////////////////////////////////////////////////////////////////////
// Render Ambient Occlusion Volume contribution from this OBJ.
/////////////////////////////////////////////////////////////////////
void GLSL_OBJ::renderAOV(float delta, float width, float height)
{
  glUseProgram(_renderAOVShader->getHandle());
  glUniform1f(_renderAOVShader->uniformLoc("delta"), delta);
  glUniform1i(_renderAOVShader->uniformLoc("normalTex"), 1);
  glUniform1i(_renderAOVShader->uniformLoc("positionTex"), 2);
  glUniform2f(_renderAOVShader->uniformLoc("screenSize"), width, height);

  render(false);

  glUseProgram(0);
}

/////////////////////////////////////////////////////////////////////
// Basic render function, binds vertex array(s) and draws.
/////////////////////////////////////////////////////////////////////
void GLSL_OBJ::render(bool normals)
{
  glTranslatef(_position[0], _position[1], _position[2]);
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(3, GL_FLOAT, 0, _vertex);

  if (normals)
  {
    glEnableClientState(GL_NORMAL_ARRAY);
    glNormalPointer(GL_FLOAT, 0, _normal);
  }

  glDrawArrays(GL_TRIANGLES, 0, _numVertices);

  if (normals)
  {
    glDisableClientState(GL_NORMAL_ARRAY);
  }
  glDisableClientState(GL_VERTEX_ARRAY);
}

