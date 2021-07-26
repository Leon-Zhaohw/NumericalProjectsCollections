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
// GLSL_TET_MESH_INSTANCE.cpp: Implementation of the GLSL_TET_MESH_INSTANCE class.
//
///////////

#include <GLSL_TET_MESH_INSTANCE.h>

#include <SIMPLE_PARSER.h>

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
GLSL_TET_MESH_INSTANCE::GLSL_TET_MESH_INSTANCE(GLSL_TET_MESH *mesh,
                                               string cfgFilename) :
  _mesh(mesh), _integrator(NULL), _integratorOptions(INTEGRATOR_STEP_DEFAULT),
  _meshRenderWithShadowShader(NULL), _meshRenderNoShadowShader(NULL),
  _meshRenderSMShader(NULL), _tetVertexShader(NULL), 
  _tetNormTanTex0Shader(NULL), _tetNormTanTex1Shader(NULL), _tetNormTanTex2Shader(NULL),
  _q(NULL)
{
  GLSL_SHADER *shader = _mesh->tetVertexShader();
  constructVertexTextures();

  SIMPLE_PARSER configFile(cfgFilename); 

    // If the uniforms don't match what we're expecting, don't store it.
  if ((shader->uniformLoc("Umatrix") != -1) && (shader->uniformLoc("q") != -1))
  {
    _tetVertexShader = shader;
  }
  else 
  {
    cout << "GLSL_TET_MESH_INSTANCE could not grab tetVertexShader from GLSL_TET_MESH" << endl;
  }

  _position[0] = _position[1] = _position[2] = 0.0f;
  _rotation[0] = _rotation[1] = _rotation[2] = 0.0f;
  startPoint = 0;

  Real timestep = 1.0f / 60.0f;
  timestep = configFile.getFloat("timestep", timestep);
  cout << " Using timestep: " << timestep << endl;

  Real rayleighAlpha = 0.01;
  Real rayleighBeta = 0.01;

  rayleighAlpha = configFile.getFloat( "rayleigh alpha", rayleighAlpha );
  rayleighBeta = configFile.getFloat( "rayleigh beta", rayleighBeta );

  _integrator = new SUBSPACE_INTEGRATOR(_mesh->tetMesh(), timestep,
                                        rayleighAlpha, rayleighBeta);

  // factor to increase the mouse-input force by
  double forceMultiplier = 100.0;
  forceMultiplier = configFile.getFloat("force multiplier", forceMultiplier);
  _integrator->forceMultiplier() = forceMultiplier;
  cout << " Using force multiplier: " << forceMultiplier << endl;

  _integrator->maxNewtonSteps() = 20;
  
  // Copy the current viewports dimensions to we can resize it during render-to-texture operations
  glGetIntegerv(GL_VIEWPORT, _windowViewport);
  //cout << "Viewport: " << _windowViewport[2] << "x" << _windowViewport[3] << endl;
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
GLSL_TET_MESH_INSTANCE::~GLSL_TET_MESH_INSTANCE()
{
  // Delete SUBSPACE_INTEGRATOR
  if (_integrator != NULL)  {
    delete _integrator;
    _integrator = NULL;
  }

  // Delete framebuffer objects
  if (_tetVertexFboID != 0)  {
    glDeleteFramebuffersEXT(1, &_tetVertexFboID);
    _tetVertexFboID = 0;
  }
  if (_tetRotationFboID != 0)  {
    glDeleteFramebuffersEXT(1, &_tetRotationFboID);
    _tetRotationFboID = 0;
  }
  if (_embeddedVertexFboID != 0)  {
    glDeleteFramebuffersEXT(1, &_embeddedVertexFboID);
    _embeddedVertexFboID = 0;
  }

  // Delete texture objects
  if (_tetVertexTexID != 0)
  {
    glDeleteTextures(1, &_tetVertexTexID);
    _tetVertexTexID = 0;
  }
  if (_tetRotationTexID[0] != 0)
  {
    glDeleteTextures(3, _tetRotationTexID);
    _tetRotationTexID[0] = 0;
  }
  if (_embeddedVertexTexID != 0)
  {
    glDeleteTextures(1, &_embeddedVertexTexID);
    _embeddedVertexTexID = 0;
  }
  if (_embeddedNormalTexID != 0)
  {
    glDeleteTextures(1, &_embeddedNormalTexID);
    _embeddedNormalTexID = 0;
  }


  if (_q != NULL)
    delete _q;
}

//////////////////////////////////////////////////////////////////////
// Draw function to display the tetrahedral mesh deformed by q
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::drawEmbeddedMeshDeformed(RenderStyleType renderStyle)
{

  // Depending on parameter, use a different shader to render
  GLSL_SHADER *activeShader;
  switch (renderStyle) 
  {
    case RENDER_STYLE_SHADOW:
      activeShader = _meshRenderWithShadowShader;
      break;
    case RENDER_STYLE_NO_SHADOW:
      activeShader = _meshRenderNoShadowShader;
      break;
    case RENDER_STYLE_SM:
      activeShader = _meshRenderSMShader;
      break;
    case RENDER_STYLE_DEFERRED:
      activeShader = _meshRenderDeferredShader;
      break;
    case RENDER_STYLE_NORMALS:
      activeShader = _meshRenderNormalsShader;
      break;
    case RENDER_STYLE_BASIC:
      activeShader = _meshRenderBasicShader;
      break;
    default:
      activeShader = _meshRenderNoShadowShader;
      break;
  }

  glUseProgram(activeShader->getHandle());

  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, _embeddedVertexTexID);
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D, _embeddedNormalTexID);

  glBindBuffer(GL_ARRAY_BUFFER, _mesh->embeddedVertLoc2VboID());
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(2, GL_FLOAT, 0, 0);

  glUniform1i(activeShader->uniformLoc("embeddedVertexTex"), 0);

  if (renderStyle != RENDER_STYLE_SM && renderStyle != RENDER_STYLE_BASIC)
    glUniform1i(activeShader->uniformLoc("embeddedNormalTex"), 1);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _mesh->embeddedMeshIboID());

  if (renderStyle == RENDER_STYLE_SHADOW)
    glUniform1i(activeShader->uniformLoc("shadowTex"), 4);
  
  if (renderStyle == RENDER_STYLE_NORMALS)
  {
    glDrawElements(GL_POINTS, _mesh->embeddedMeshIboSize(), GL_UNSIGNED_INT, 0);
  } else {
    glDrawElements(GL_TRIANGLES, _mesh->embeddedMeshIboSize(), GL_UNSIGNED_INT, 0);
  }

  glUseProgram(0);

  glBindTexture(GL_TEXTURE_2D, 0);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, 0);
  glDisableClientState(GL_VERTEX_ARRAY);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  
}

//////////////////////////////////////////////////////////////////////
// Renders the embedded mesh using the pregenerated deformed vertex positions 
// in the embedded mesh textures.
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::drawEmbeddedMeshFromTexture()
{
  /*glEnable(GL_TEXTURE_2D);
  glDisable(GL_CULL_FACE);

  glBindTexture(GL_TEXTURE_2D, _tetVertexTexID);
  
  glBegin(GL_QUADS);
    glTexCoord2f(0.0f,0.0f);    glVertex3f(-0.75f,  1.0f, 0.0f);
    glTexCoord2f(1.0f,0.0f);    glVertex3f( 0.25f,  1.0f, 0.0f);
    glTexCoord2f(1.0f,1.0f);    glVertex3f( 0.25f,  2.0f, 0.0f);
    glTexCoord2f(0.0f,1.0f);    glVertex3f(-0.75f,  2.0f, 0.0f);
  glEnd();

  glBindTexture(GL_TEXTURE_2D, _embeddedVertexTexID);
  
  glBegin(GL_QUADS);
    glTexCoord2f(0.0f,0.0f);    glVertex3f( 0.75f,  1.0f, 0.0f);
    glTexCoord2f(1.0f,0.0f);    glVertex3f( 1.75f,  1.0f, 0.0f);
    glTexCoord2f(1.0f,1.0f);    glVertex3f( 1.75f,  2.0f, 0.0f);
    glTexCoord2f(0.0f,1.0f);    glVertex3f( 0.75f,  2.0f, 0.0f);
  glEnd();

  glBindTexture(GL_TEXTURE_2D, 0);
  glDisable(GL_TEXTURE_2D);
  glEnable(GL_CULL_FACE);// */

  cout << " Should not be called." << endl;


  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, _embeddedVertexTexID);
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D, _embeddedNormalTexID);

  glBindBuffer(GL_ARRAY_BUFFER, _mesh->embeddedVertLoc2VboID());
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(2, GL_FLOAT, 0, 0);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _mesh->embeddedMeshIboID());

  GLSL_SHADER *shader = _mesh->renderFromTextureTest();
  glUseProgram(shader->getHandle());
  glUniform1i(shader->uniformLoc("embeddedVertexTex"), 0);
  glUniform1i(shader->uniformLoc("embeddedNormalTex"), 1);

  glDrawElements(GL_TRIANGLES, _mesh->embeddedMeshIboSize(), GL_UNSIGNED_INT, 0);

  glBindTexture(GL_TEXTURE_2D, 0);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, 0);
  glDisableClientState(GL_VERTEX_ARRAY);

  glUseProgram(0);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
}

//////////////////////////////////////////////////////////////////////
// Generates deformed tet mesh texture, then the deformed embedded mesh textures.
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::genTetVertexTexture()
{
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _tetVertexFboID);

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  //glClearColor(0.0,0.0,0.0,0.0);

  glOrtho(0.0, 1.0, 0.0, 1.0, -1.0, 1.0);

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  glViewport(0, 0, _tetVertexTexWidth, _tetVertexTexHeight);

  glClear(GL_COLOR_BUFFER_BIT);

  glPointSize(1.0);

  copyQForUniform();

  glBindBuffer(GL_ARRAY_BUFFER, _mesh->tetMeshVboID());
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(3, GL_FLOAT, 0, 0);

  glClientActiveTexture(GL_TEXTURE0);
  glBindBuffer(GL_ARRAY_BUFFER, _mesh->UcoordVboID());
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);
  glTexCoordPointer(2, GL_FLOAT, 0, 0);

  glClientActiveTexture(GL_TEXTURE1);
  glBindBuffer(GL_ARRAY_BUFFER, _mesh->tetMeshInfoVboID());
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);
  glTexCoordPointer(3, GL_FLOAT, 0, 0);

  // Don't need to do a glEnable(GL_TEXTURE_2D) when using shaders
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, _mesh->UbasisTexID());

  if (_tetVertexShader != NULL)
  {
    glUseProgram(_tetVertexShader->getHandle());

    glUniform1i(_tetVertexShader->uniformLoc("Umatrix"), 0);
    glUniform4fv(_tetVertexShader->uniformLoc("q"), _qVecSize, _q);
  }

  glDrawArrays(GL_POINTS, 0, _mesh->tetMeshSize());

  glViewport(0, 0, _embeddedVertexTexWidth, _embeddedVertexTexHeight);
  //glClearColor(1.0,1.0,1.0,1.0);
  
  glBlendEquation(GL_FUNC_ADD);
  glBlendFunc(GL_ONE, GL_ONE);
  glEnable(GL_BLEND);
  glDisable(GL_DEPTH_TEST);

  glBindBuffer(GL_ARRAY_BUFFER, _mesh->embeddedMeshVboID());
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(3, GL_FLOAT, 0, 0);

  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, _tetVertexTexID);

  glClientActiveTexture(GL_TEXTURE0);
  glBindBuffer(GL_ARRAY_BUFFER, _mesh->embeddedVertLoc0VboID());
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);
  glTexCoordPointer(4, GL_FLOAT, 0, 0);

  glClientActiveTexture(GL_TEXTURE1);
  glBindBuffer(GL_ARRAY_BUFFER, _mesh->embeddedVertLoc1VboID());
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);
  glTexCoordPointer(4, GL_FLOAT, 0, 0);

  glClientActiveTexture(GL_TEXTURE2);
  glBindBuffer(GL_ARRAY_BUFFER, _mesh->embeddedVertLoc2VboID());
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);
  glTexCoordPointer(2, GL_FLOAT, 0, 0);

  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _embeddedVertexFboID);

  glClear(GL_COLOR_BUFFER_BIT);

  if (_meshRenderToTextureShader != NULL)
  {
    glUseProgram(_meshRenderToTextureShader->getHandle());
    glUniform1i(_meshRenderToTextureShader->uniformLoc("tetMeshTex"), 0);
  }

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _mesh->embeddedMeshIboID());

  glDrawElements(GL_TRIANGLES, _mesh->embeddedMeshIboSize(), GL_UNSIGNED_INT, 0);

  ///////////
 
  glDisableClientState(GL_TEXTURE_COORD_ARRAY);
  glClientActiveTexture(GL_TEXTURE1);
  glDisableClientState(GL_TEXTURE_COORD_ARRAY);
  glClientActiveTexture(GL_TEXTURE0);
  glDisableClientState(GL_TEXTURE_COORD_ARRAY);
  glDisableClientState(GL_VERTEX_ARRAY);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0); 

  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
  glUseProgram(0);

  glBindTexture(GL_TEXTURE_2D, 0);

  glDisable(GL_BLEND);
  glEnable(GL_DEPTH_TEST);
  
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();

  glViewport(_windowViewport[0], _windowViewport[1], _windowViewport[2], _windowViewport[3]);
}

//////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::renderBases()
{
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, _tetVertexTexID);
  /*glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D, _tetRotationTexID[0]);
  glActiveTexture(GL_TEXTURE2);
  glBindTexture(GL_TEXTURE_2D, _tetRotationTexID[1]);
  glActiveTexture(GL_TEXTURE3);
  glBindTexture(GL_TEXTURE_2D, _tetRotationTexID[2]);// */
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D, _mesh->n1RestTexID());
  glActiveTexture(GL_TEXTURE2);
  glBindTexture(GL_TEXTURE_2D, _mesh->n2RestTexID());// */

  glDisable(GL_CULL_FACE);


  glBindBuffer(GL_ARRAY_BUFFER, _mesh->surfaceVertVboID());
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(2, GL_FLOAT, 0, 0);

  glClientActiveTexture(GL_TEXTURE0);
  glBindBuffer(GL_ARRAY_BUFFER, _mesh->surfaceVertNeighbor0VboID());
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);
  glTexCoordPointer(4, GL_FLOAT, 0, 0);

  glClientActiveTexture(GL_TEXTURE1);
  glBindBuffer(GL_ARRAY_BUFFER, _mesh->surfaceVertNeighbor1VboID());
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);
  glTexCoordPointer(4, GL_FLOAT, 0, 0);// */

  //glClear(GL_COLOR_BUFFER_BIT);

  GLfloat offset[2];
  offset[0] = (0.5 / (double)_tetVertexTexWidth);
  offset[1] = (0.5 / (double)_tetVertexTexHeight);

  glLineWidth(1.0);

  assert(_mesh->drawBasesShader() != NULL);

  glUseProgram(_mesh->drawBasesShader()->getHandle());
  glUniform1i(_mesh->drawBasesShader()->uniformLoc("tetMeshTex"), 0);
  /*glUniform1i(_mesh->drawBasesShader()->uniformLoc("Rotation1"), 1);
  glUniform1i(_mesh->drawBasesShader()->uniformLoc("Rotation2"), 2);
  glUniform1i(_mesh->drawBasesShader()->uniformLoc("Rotation3"), 3); // */
  glUniform1i(_mesh->drawBasesShader()->uniformLoc("n1RestFrame"), 1);
  glUniform1i(_mesh->drawBasesShader()->uniformLoc("n2RestFrame"), 2); // */

  glUniform2f(_mesh->drawBasesShader()->uniformLoc("offset"), offset[0], offset[1]);

  glDrawArrays(GL_POINTS, 0, _mesh->surfaceVertSize());

  glUseProgram(_mesh->drawNeighborsShader()->getHandle());
  glUniform1i(_mesh->drawNeighborsShader()->uniformLoc("tetMeshTex"), 0);
  glUniform2f(_mesh->drawNeighborsShader()->uniformLoc("offset"), offset[0], offset[1]);

  glLineWidth(4.0);

  glDrawArrays(GL_POINTS, 0, _mesh->surfaceVertSize());

  glDisableClientState(GL_TEXTURE_COORD_ARRAY);
  glClientActiveTexture(GL_TEXTURE0);
  glDisableClientState(GL_TEXTURE_COORD_ARRAY);
  glDisableClientState(GL_VERTEX_ARRAY);// */

  //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0); 

  glUseProgram(0);

  glEnable(GL_CULL_FACE);

  /*glBindTexture(GL_TEXTURE_2D, 0);
  glActiveTexture(GL_TEXTURE2);// */
  glBindTexture(GL_TEXTURE_2D, 0);
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D, 0);
  glActiveTexture(GL_TEXTURE0);// */
  glBindTexture(GL_TEXTURE_2D, 0);
}

//////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::renderFaces()
{
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, _tetVertexTexID);

  glBindBuffer(GL_ARRAY_BUFFER, _mesh->tetFaceInfo0VboID());
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(4, GL_FLOAT, 0, 0);

  glClientActiveTexture(GL_TEXTURE0);
  glBindBuffer(GL_ARRAY_BUFFER, _mesh->tetFaceInfo1VboID());
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);
  glTexCoordPointer(2, GL_FLOAT, 0, 0);
  
  //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  GLfloat offset[2];
  offset[0] = (0.5 / (double)_tetVertexTexWidth);
  offset[1] = (0.5 / (double)_tetVertexTexHeight);

  assert(_mesh->drawFacesShader() != NULL);
  
  glUseProgram(_mesh->drawFacesShader()->getHandle());
  glUniform1i(_mesh->drawFacesShader()->uniformLoc("tetMeshTex"), 0);
  glUniform2f(_mesh->drawFacesShader()->uniformLoc("offset"), offset[0], offset[1]);
  
  glDrawArrays(GL_POINTS, 0, _mesh->tetFaceSize());

  glDisableClientState(GL_TEXTURE_COORD_ARRAY);
  glDisableClientState(GL_VERTEX_ARRAY);// */

  //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0); 

  glUseProgram(0);

  glBindTexture(GL_TEXTURE_2D, 0);
}

//////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::renderAOVEmbedded(float delta, float width, float height)
{
  // Bind required vertex buffers
  glBindBuffer(GL_ARRAY_BUFFER, _mesh->embeddedVertLoc2VboID());
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(2, GL_FLOAT, 0, 0);

   // Don't need to do a glEnable(GL_TEXTURE_2D) when using shaders
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, _embeddedVertexTexID);

  assert(_mesh->drawAOV2Shader() != NULL);
  
  glUseProgram(_mesh->drawAOV2Shader()->getHandle());
  glUniform1i(_mesh->drawAOV2Shader()->uniformLoc("embeddedVertexTex"), 0);
  glUniform1f(_mesh->drawAOV2Shader()->uniformLoc("delta"), delta);
  
  //cout << "screen size: " << _windowViewport[2] << " x " << _windowViewport[3] << endl;
  glUniform1i(_mesh->drawAOV2Shader()->uniformLoc("normalTex"), 1);
  glUniform1i(_mesh->drawAOV2Shader()->uniformLoc("positionTex"), 2);
  glUniform2f(_mesh->drawAOV2Shader()->uniformLoc("screenSize"), width, height);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _mesh->embeddedMeshIboID());
  

  if (startPoint < 0)
    startPoint += _mesh->tetFaceSize();
  if (startPoint >= _mesh->tetFaceSize())
    startPoint -= _mesh->tetFaceSize();

//  glDrawArrays(GL_POINTS, 0, _mesh->tetFaceSize());
  //glDrawArrays(GL_POINTS, startPoint, 20);
  //glDrawElements(GL_TRIANGLES, 27, GL_UNSIGNED_INT, (GLvoid*)(3*startPoint*sizeof(GLint)));
  glDrawElements(GL_TRIANGLES, _mesh->embeddedMeshIboSize(), GL_UNSIGNED_INT, 0);

  
  glDisableClientState(GL_VERTEX_ARRAY);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0); 

  glUseProgram(0);

  glBindTexture(GL_TEXTURE_2D, 0);
}

//////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::renderAOVWireframe(float delta)
{
  // Bind required vertex buffers
  glBindBuffer(GL_ARRAY_BUFFER, _mesh->embeddedMeshVboID());
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(3, GL_FLOAT, 0, 0);

   // Don't need to do a glEnable(GL_TEXTURE_2D) when using shaders
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, _tetVertexTexID);

  glClientActiveTexture(GL_TEXTURE0);
  glBindBuffer(GL_ARRAY_BUFFER, _mesh->embeddedVertLoc0VboID());
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);
  glTexCoordPointer(4, GL_FLOAT, 0, 0);

  glClientActiveTexture(GL_TEXTURE1);
  glBindBuffer(GL_ARRAY_BUFFER, _mesh->embeddedVertLoc1VboID());
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);
  glTexCoordPointer(4, GL_FLOAT, 0, 0);
  
  glUseProgram(_mesh->drawAOV3Shader()->getHandle());
  glUniform1i(_mesh->drawAOV3Shader()->uniformLoc("tetMeshTex"), 0);
  glUniform1f(_mesh->drawAOV3Shader()->uniformLoc("delta"), delta);
  
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _mesh->embeddedMeshIboID());
  
  if (startPoint < 0)
    startPoint += _mesh->tetFaceSize();
  if (startPoint >= _mesh->tetFaceSize())
    startPoint -= _mesh->tetFaceSize();

  //glDrawElements(GL_TRIANGLES, 27, GL_UNSIGNED_INT, (GLvoid*)(3*startPoint*sizeof(GLint)));
  glDrawElements(GL_TRIANGLES, _mesh->embeddedMeshIboSize(), GL_UNSIGNED_INT, 0);

  glDisableClientState(GL_TEXTURE_COORD_ARRAY);
  glClientActiveTexture(GL_TEXTURE0);
  glDisableClientState(GL_TEXTURE_COORD_ARRAY);

  glDisableClientState(GL_VERTEX_ARRAY);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0); 

  glUseProgram(0);

  glBindTexture(GL_TEXTURE_2D, 0);
}


//////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::renderAOV(float delta)
{
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, _tetVertexTexID);

  glBindBuffer(GL_ARRAY_BUFFER, _mesh->tetFaceInfo0VboID());
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(4, GL_FLOAT, 0, 0);

  glClientActiveTexture(GL_TEXTURE0);
  glBindBuffer(GL_ARRAY_BUFFER, _mesh->tetFaceInfo1VboID());
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);
  glTexCoordPointer(2, GL_FLOAT, 0, 0);
  
  //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  GLfloat offset[2];
  offset[0] = (0.5 / (double)_tetVertexTexWidth);
  offset[1] = (0.5 / (double)_tetVertexTexHeight);

  assert(_mesh->drawAOVShader() != NULL);
  
  glUseProgram(_mesh->drawAOVShader()->getHandle());
  glUniform1i(_mesh->drawAOVShader()->uniformLoc("tetMeshTex"), 0);
  glUniform2f(_mesh->drawAOVShader()->uniformLoc("offset"), offset[0], offset[1]);
  glUniform1f(_mesh->drawAOVShader()->uniformLoc("delta"), delta);

  //cout << "screen size: " << _windowViewport[2] << " x " << _windowViewport[3] << endl;
  glUniform1i(_mesh->drawAOVShader()->uniformLoc("normalTex"), 1);
  glUniform1i(_mesh->drawAOVShader()->uniformLoc("positionTex"), 2);
  glUniform2f(_mesh->drawAOVShader()->uniformLoc("screenSize"), (GLfloat)_windowViewport[2], (GLfloat)_windowViewport[3]);
  
  if (startPoint < 0)
    startPoint += _mesh->tetFaceSize();
  if (startPoint >= _mesh->tetFaceSize())
    startPoint -= _mesh->tetFaceSize();

  glDrawArrays(GL_POINTS, 0, _mesh->tetFaceSize());
  //glDrawArrays(GL_POINTS, startPoint, 20);

  glDisableClientState(GL_TEXTURE_COORD_ARRAY);
  glDisableClientState(GL_VERTEX_ARRAY);// */

  //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0); 

  glUseProgram(0);

  glBindTexture(GL_TEXTURE_2D, 0);
}


//////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::setIntegratorOptions(IntegratorOptionType options)
{
  _integratorOptions = options;
}

//////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::integratorStep( bool debug )
{
  if (_integrator != NULL)
  {
    switch (_integratorOptions)
    {
      case INTEGRATOR_STEP_QUASISTATIC:
        _integrator->stepQuasistatic();
        break;
      case INTEGRATOR_STEP_IMPLICIT:
        _integrator->stepImplicit(true);
        break;
      case INTEGRATOR_STEP_INVERTIBLE_IMPLICIT:
        _integrator->stepInvertibleImplicit(true, debug);
        break;
      case INTEGRATOR_STEP_INVERTIBLE_QUASISTATIC:
        _integrator->stepInvertibleQuasistatic(true);
        break;
      case INTEGRATOR_STEP_EXPLICIT:
        _integrator->stepExplicit();
        break;
      case INTEGRATOR_STEP_INVERTIBLE_SEMI_IMPLICIT:
        _integrator->stepInvertibleSemiImplicit(true);
        break;
      case INTEGRATOR_STEP_DEFAULT:
      default:
        _integrator->stepInvertibleQuasistatic(true);
        break;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////
bool GLSL_TET_MESH_INSTANCE::integratorClick(VEC3F& point, Real maxRadius)
{
  if (_integrator != NULL)
  {
    return _integrator->click(point, maxRadius);
  }
  else
  {
    return false;
  }
}

//////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::integratorDrag(VEC3F& point)
{
  if (_integrator != NULL)
  {
    _integrator->drag(point);
  }
}

//////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::integratorUnclick()
{
  if (_integrator != NULL)
  {
    _integrator->unclick();
  }
}

//////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::setMeshRenderToTextureShader(GLSL_SHADER *shader)
{
  // If the uniforms don't match what we're expecting, don't store it.
  if (shader->uniformLoc("tetMeshTex") == -1)
    return;
  _meshRenderToTextureShader = shader;
}

//////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::setMeshRenderWithShadowShader(GLSL_SHADER *shader)
{
  // If the uniforms don't match what we're expecting, don't store it.
  //if (shader->uniformLoc("tetMeshTex") == -1)
  //  return;
  //if (shader->uniformLoc("tetRotation1Tex") == -1)
  //  return;
  //if (shader->uniformLoc("tetRotation2Tex") == -1)
  //  return;
  //if (shader->uniformLoc("tetRotation3Tex") == -1)
  //  return;
  _meshRenderWithShadowShader = shader;
}

//////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::setMeshRenderNoShadowShader(GLSL_SHADER *shader)
{
  // If the uniforms don't match what we're expecting, don't store it.
  //if (shader->uniformLoc("tetMeshTex") == -1)
  //  return;
  //if (shader->uniformLoc("tetRotation1Tex") == -1)
  //  return;
  //if (shader->uniformLoc("tetRotation2Tex") == -1)
  //  return;
  //if (shader->uniformLoc("tetRotation3Tex") == -1)
  //  return;
  _meshRenderNoShadowShader = shader;
}


//////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::setMeshRenderSMShader(GLSL_SHADER *shader)
{
  // If the uniforms don't match what we're expecting, don't store it.
  //if (shader->uniformLoc("embeddedVertexTex") == -1)
  //  return;
  _meshRenderSMShader = shader;
}

//////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::setMeshRenderDeferredShader(GLSL_SHADER *shader)
{
  // If the uniforms don't match what we're expecting, don't store it.
  //if (shader->uniformLoc("tetMeshTex") == -1)
  //  return;
  _meshRenderDeferredShader = shader;
}

//////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::setMeshRenderBasicShader(GLSL_SHADER *shader)
{
  _meshRenderBasicShader = shader;
}

//////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::setMeshRenderNormalsShader(GLSL_SHADER *shader)
{
   // If the uniforms don't match what we're expecting, don't store it.
  if (shader->uniformLoc("tetMeshTex") == -1)
    return;
  if (shader->uniformLoc("tetRotation1Tex") == -1)
    return;
  if (shader->uniformLoc("tetRotation2Tex") == -1)
    return;
  if (shader->uniformLoc("tetRotation3Tex") == -1)
    return;
  _meshRenderNormalsShader = shader;
}

//////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::setTetVertexShader(GLSL_SHADER *shader)
{
  if (shader->uniformLoc("Umatrix") == -1)
    return;
  if (shader->uniformLoc("q") == -1)
    return;
  _tetVertexShader = shader;
}

//////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::setNormTanTexShader(GLSL_SHADER *shader, int passnum)
{
  // If the uniforms don't match what we're expecting, don't store it.
  if (shader->uniformLoc("tetMeshTex") == -1)
    return;
  if (shader->uniformLoc("offset") == -1)
    return;
  switch (passnum)
  {
    case 0:
      _tetNormTanTex0Shader = shader;
      break;
    case 1:
      _tetNormTanTex1Shader = shader;
      break;
    case 2:
      _tetNormTanTex2Shader = shader;
      break;
    default:
      break;
  }
}


//////////////////////////////////////////////////////////////////////
// Helper function to copy q from tetMesh to our local _q, convert
// it and pack it correctly for use as a uniform in the shaders.
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::copyQForUniform()
{
  int meshQSize = _mesh->tetMesh()->q().size();
  //cout << "CopyQForUniform: q().size() is " << meshQSize << endl;

  if (_q == NULL)
  {
    _qVecSize = (GLuint)(((meshQSize % 4) > 0) 
                  ? ((meshQSize / 4) + 1) 
                  :  (meshQSize / 4));
    //cout << "q was found to be NULL, allocating array of " << _qVecSize << " vec4 variables." << endl;
    _q = new GLfloat[_qVecSize*4];
    if (_q == NULL)
    {
      // Something really went wrong so bail
      cout << " Could not allocate memory for SUBSPACE_TET_MESH_RENDERER::_q" << endl;
      exit(1);
    }

    int vecQSize = _qVecSize*4;
    for (int x = 0; x < vecQSize; x++)
    {
      _q[x] = 0.0f;      
    }

  }

  /*float max, min, sum, mean, stdev;
  max = -1000000.0f;
  min = 1000000.0f;
  sum = 0.0f;
  mean = 0.0f;
  stdev = 0.0f;// */

  for (int x = 0; x < meshQSize; x++)
  {
    _q[x] = (GLfloat)(_mesh->tetMesh()->q()[x]);
    /*if (_q[x] > max)
      max = _q[x];
    if (_q[x] < min)
      min = _q[x];
    sum += _q[x];// */
  }
  /*mean = sum / (float)meshQSize;
  for (int x = 0; x < meshQSize; x++)
  {
    stdev += ((_q[x] - mean)*(_q[x] - mean));
  }
  stdev = stdev / (float)(meshQSize-1);
  stdev = sqrt(stdev);

  cout << "-[ q ]-" << endl;
  cout << "mean : " << mean << endl;
  cout << "stdev: " << stdev << endl;
  cout << "range: " << min << ", " << max << endl;// */

}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::constructVertexTextures()
{
  glGenTextures(1, &_tetVertexTexID);
  if (_tetVertexTexID == 0)
  {
    cout << "Could not allocate memory for texture object" << endl;
    return;
  }

  glBindTexture(GL_TEXTURE_2D, _tetVertexTexID);

  // Grab texture dimensions from the mesh
  _mesh->getVertexTexDims(&_tetVertexTexWidth, &_tetVertexTexHeight,
                          &_embeddedVertexTexWidth, &_embeddedVertexTexHeight);

  // _tetVertexTexID will hold deformed positions of tet vertices
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, _tetVertexTexWidth,
               _tetVertexTexHeight, 0, GL_RGBA, GL_FLOAT, 0);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

  glGenTextures(3, _tetRotationTexID);
  if (_tetRotationTexID == 0)
  {
    cout << "Could not allocate memory for texture object" << endl;
    return;
  }

  for (int x = 0; x < 3; x++)
  {
    glBindTexture(GL_TEXTURE_2D, _tetRotationTexID[x]);
   
    // _tetVertexTexID will hold deformed positions of tet vertices
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, _tetVertexTexWidth,
                 _tetVertexTexHeight, 0, GL_RGBA, GL_FLOAT, 0);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  }

  glGenFramebuffersEXT(1, &_tetVertexFboID);
  if (_tetVertexFboID == 0)
  {
    cout << " Could not generate tet vertex framebuffer object." << endl;
    return;
  }

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _tetVertexFboID);
  GLenum vbuf[1] = { GL_COLOR_ATTACHMENT0_EXT };
  glDrawBuffers(1, vbuf);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT,
                            GL_TEXTURE_2D, _tetVertexTexID, 0);
 
  GLenum fboStatus = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
  if (fboStatus != GL_FRAMEBUFFER_COMPLETE_EXT)
  {
   cout << " Framebuffer status is not complete; an error occurred at ";
   cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
   return; 
  }

  glGenFramebuffersEXT(1, &_tetRotationFboID);
  if (_tetRotationFboID == 0)
  {
    cout << " Could not generate tet vertex framebuffer object." << endl;
    return;
  }

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _tetRotationFboID);
  GLenum rbuf[3] = { GL_COLOR_ATTACHMENT0_EXT, GL_COLOR_ATTACHMENT1_EXT,
                     GL_COLOR_ATTACHMENT2_EXT };
  glDrawBuffers(3, rbuf);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT,
                            GL_TEXTURE_2D, _tetRotationTexID[0], 0);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT1_EXT,
                            GL_TEXTURE_2D, _tetRotationTexID[1], 0);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT2_EXT,
                            GL_TEXTURE_2D, _tetRotationTexID[2], 0);

  fboStatus = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
  if (fboStatus != GL_FRAMEBUFFER_COMPLETE_EXT)
  {
   cout << " Framebuffer status is not complete; an error occurred at "  << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
   return; 
  }

  glGenTextures(1, &_embeddedVertexTexID);
  if (_embeddedVertexTexID == 0)
  {
    cout << "Could not allocate memory for texture object" << endl;
    return;
  }
  glBindTexture(GL_TEXTURE_2D, _embeddedVertexTexID);

   // _embeddedVertexTexID will hold deformed positions of embedded vertices
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, _embeddedVertexTexWidth,
               _embeddedVertexTexHeight, 0, GL_RGBA, GL_FLOAT, 0);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

  glGenTextures(1, &_embeddedNormalTexID);
  if (_embeddedNormalTexID == 0)
  {
    cout << "Could not allocate memory for texture object" << endl;
    return;
  }
  glBindTexture(GL_TEXTURE_2D, _embeddedNormalTexID);

   // _embeddedNormalTexID will hold deformed normals of embedded vertices
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, _embeddedVertexTexWidth,
               _embeddedVertexTexHeight, 0, GL_RGBA, GL_FLOAT, 0);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

  glGenFramebuffersEXT(1, &_embeddedVertexFboID);
  if (_embeddedVertexFboID == 0)
  {
    cout << " Could not generate embedded vertex framebuffer object." << endl;
    return;
  }

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _embeddedVertexFboID);
  GLenum ebuf[2] = { GL_COLOR_ATTACHMENT0_EXT, GL_COLOR_ATTACHMENT1_EXT };
  glDrawBuffers(2, ebuf);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT,
                            GL_TEXTURE_2D, _embeddedVertexTexID, 0);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT1_EXT,
                            GL_TEXTURE_2D, _embeddedNormalTexID, 0);

  fboStatus = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
  if (fboStatus != GL_FRAMEBUFFER_COMPLETE_EXT)
  {
   cout << " Framebuffer status is not complete; an error occurred at "  << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
   return; 
  }

  //glClearColor(0.2, 0.2, 0.2, 1.0);

  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
}


//////////////////////////////////////////////////////////////////////
// Draw function to display the tetrahedral mesh deformed by q
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::oldDrawEmbeddedMeshDeformed(int renderStyle)
{
  // Draw the tetVertexTex
  /*glEnable(GL_TEXTURE_2D);
  glDisable(GL_CULL_FACE);

  glBindTexture(GL_TEXTURE_2D, _tetRotationTexID[0]);
  
  glBegin(GL_QUADS);
    glTexCoord2f(0.0f,0.0f);    glVertex3f(-1.0f,  0.0f, 0.0f);
    glTexCoord2f(1.0f,0.0f);    glVertex3f( 0.0f,  0.0f, 0.0f);
    glTexCoord2f(1.0f,1.0f);    glVertex3f( 0.0f,  1.0f, 0.0f);
    glTexCoord2f(0.0f,1.0f);    glVertex3f(-1.0f,  1.0f, 0.0f);
  glEnd();

  glBindTexture(GL_TEXTURE_2D, _tetRotationTexID[1]);
  
  glBegin(GL_QUADS);
    glTexCoord2f(0.0f,0.0f);    glVertex3f( 0.0f,  0.0f, 0.0f);
    glTexCoord2f(1.0f,0.0f);    glVertex3f( 1.0f,  0.0f, 0.0f);
    glTexCoord2f(1.0f,1.0f);    glVertex3f( 1.0f,  1.0f, 0.0f);
    glTexCoord2f(0.0f,1.0f);    glVertex3f( 0.0f,  1.0f, 0.0f);
  glEnd();

  glBindTexture(GL_TEXTURE_2D, _tetRotationTexID[2]);
  
  glBegin(GL_QUADS);
    glTexCoord2f(0.0f,0.0f);    glVertex3f( 1.0f,  0.0f, 0.0f);
    glTexCoord2f(1.0f,0.0f);    glVertex3f( 2.0f,  0.0f, 0.0f);
    glTexCoord2f(1.0f,1.0f);    glVertex3f( 2.0f,  1.0f, 0.0f);
    glTexCoord2f(0.0f,1.0f);    glVertex3f( 1.0f,  1.0f, 0.0f);
  glEnd();

  glBindTexture(GL_TEXTURE_2D, _mesh->n1RestTexID());
  
  glBegin(GL_QUADS);
    glTexCoord2f(0.0f,0.0f);    glVertex3f(-0.75f,  1.0f, 0.0f);
    glTexCoord2f(1.0f,0.0f);    glVertex3f( 0.25f,  1.0f, 0.0f);
    glTexCoord2f(1.0f,1.0f);    glVertex3f( 0.25f,  2.0f, 0.0f);
    glTexCoord2f(0.0f,1.0f);    glVertex3f(-0.75f,  2.0f, 0.0f);
  glEnd();

  glBindTexture(GL_TEXTURE_2D, _mesh->n2RestTexID());
  
  glBegin(GL_QUADS);
    glTexCoord2f(0.0f,0.0f);    glVertex3f( 0.75f,  1.0f, 0.0f);
    glTexCoord2f(1.0f,0.0f);    glVertex3f( 1.75f,  1.0f, 0.0f);
    glTexCoord2f(1.0f,1.0f);    glVertex3f( 1.75f,  2.0f, 0.0f);
    glTexCoord2f(0.0f,1.0f);    glVertex3f( 0.75f,  2.0f, 0.0f);
  glEnd();


  glBindTexture(GL_TEXTURE_2D, 0);
  glDisable(GL_TEXTURE_2D);
  glEnable(GL_CULL_FACE);// */

  // Depending on parameter, use a different shader to render
  GLSL_SHADER *activeShader;
  switch (renderStyle) 
  {
    case RENDER_STYLE_SHADOW:
      activeShader = _meshRenderWithShadowShader;
      break;
    case RENDER_STYLE_NO_SHADOW:
      activeShader = _meshRenderNoShadowShader;
      break;
    case RENDER_STYLE_SM:
      activeShader = _meshRenderSMShader;
      break;
    case RENDER_STYLE_DEFERRED:
      activeShader = _meshRenderDeferredShader;
      break;
    case RENDER_STYLE_NORMALS:
      activeShader = _meshRenderNormalsShader;
      break;
    case RENDER_STYLE_BASIC:
      activeShader = _meshRenderBasicShader;
      break;
    default:
      activeShader = _meshRenderNoShadowShader;
      break;
  }

  glUseProgram(activeShader->getHandle());

  // Bind required vertex buffers
  glBindBuffer(GL_ARRAY_BUFFER, _mesh->embeddedMeshVboID());
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(3, GL_FLOAT, 0, 0);

   // Don't need to do a glEnable(GL_TEXTURE_2D) when using shaders
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, _tetVertexTexID);

  glClientActiveTexture(GL_TEXTURE0);
  glBindBuffer(GL_ARRAY_BUFFER, _mesh->embeddedVertLoc0VboID());
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);
  glTexCoordPointer(4, GL_FLOAT, 0, 0);

  glClientActiveTexture(GL_TEXTURE1);
  glBindBuffer(GL_ARRAY_BUFFER, _mesh->embeddedVertLoc1VboID());
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);
  glTexCoordPointer(4, GL_FLOAT, 0, 0);

  glUniform1i(activeShader->uniformLoc("tetMeshTex"), 0);

  if (renderStyle != RENDER_STYLE_SM)
  {
    glBindBuffer(GL_ARRAY_BUFFER, _mesh->embeddedNormalVboID());
    glEnableClientState(GL_NORMAL_ARRAY);
    glNormalPointer(GL_FLOAT, 0, 0);
    
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, _tetRotationTexID[0]);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, _tetRotationTexID[1]);
    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, _tetRotationTexID[2]);

    glUniform1i(activeShader->uniformLoc("tetRotation1Tex"), 1);
    glUniform1i(activeShader->uniformLoc("tetRotation2Tex"), 2);
    glUniform1i(activeShader->uniformLoc("tetRotation3Tex"), 3);
  }
    
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _mesh->embeddedMeshIboID());
  
  if (renderStyle == RENDER_STYLE_SHADOW)
    glUniform1i(activeShader->uniformLoc("shadowTex"), 4);
  
  if (renderStyle == RENDER_STYLE_NORMALS)
  {
    glDrawElements(GL_POINTS, _mesh->embeddedMeshIboSize(), GL_UNSIGNED_INT, 0);
  } else {
    glDrawElements(GL_TRIANGLES, _mesh->embeddedMeshIboSize(), GL_UNSIGNED_INT, 0);
  }

  glUseProgram(0);

  if (renderStyle != RENDER_STYLE_SM)
  {
    glDisableClientState(GL_NORMAL_ARRAY);

    glBindTexture(GL_TEXTURE_2D, 0);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, 0);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, 0);
    glActiveTexture(GL_TEXTURE0);
  }

  glBindTexture(GL_TEXTURE_2D, 0);

  glDisableClientState(GL_TEXTURE_COORD_ARRAY);
  glClientActiveTexture(GL_TEXTURE0);
  glDisableClientState(GL_TEXTURE_COORD_ARRAY);

  glDisableClientState(GL_VERTEX_ARRAY);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  
}


//////////////////////////////////////////////////////////////////////
// Draw function to display the tetrahedral mesh deformed by q
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::oldGenTetVertexTexture()
{
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _tetVertexFboID);

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();

  glOrtho(0.0, 1.0, 0.0, 1.0, -1.0, 1.0);

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  glViewport(0, 0, _tetVertexTexWidth, _tetVertexTexHeight);

  glClear(GL_COLOR_BUFFER_BIT);

  glPointSize(1.0);

  copyQForUniform();

  glBindBuffer(GL_ARRAY_BUFFER, _mesh->tetMeshVboID());
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(3, GL_FLOAT, 0, 0);

  glClientActiveTexture(GL_TEXTURE0);
  glBindBuffer(GL_ARRAY_BUFFER, _mesh->UcoordVboID());
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);
  glTexCoordPointer(2, GL_FLOAT, 0, 0);

  glClientActiveTexture(GL_TEXTURE1);
  glBindBuffer(GL_ARRAY_BUFFER, _mesh->tetMeshInfoVboID());
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);
  glTexCoordPointer(3, GL_FLOAT, 0, 0);

  // Don't need to do a glEnable(GL_TEXTURE_2D) when using shaders
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, _mesh->UbasisTexID());

  if (_tetVertexShader != NULL)
  {
    glUseProgram(_tetVertexShader->getHandle());

    glUniform1i(_tetVertexShader->uniformLoc("Umatrix"), 0);
    glUniform4fv(_tetVertexShader->uniformLoc("q"), _qVecSize, _q);
  }

  glDrawArrays(GL_POINTS, 0, _mesh->tetMeshSize());

  glBindTexture(GL_TEXTURE_2D, _tetVertexTexID);
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D, _mesh->n1RestTexID());
  glActiveTexture(GL_TEXTURE2);
  glBindTexture(GL_TEXTURE_2D, _mesh->n2RestTexID());


  //glBindBuffer(GL_ARRAY_BUFFER, _mesh->tetFaceInfo0VboID());
  glBindBuffer(GL_ARRAY_BUFFER, _mesh->surfaceVertVboID());
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(2, GL_FLOAT, 0, 0);

  //glDisableClientState(GL_TEXTURE_COORD_ARRAY);
  glClientActiveTexture(GL_TEXTURE0);
  glBindBuffer(GL_ARRAY_BUFFER, _mesh->surfaceVertNeighbor0VboID());
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);
  glTexCoordPointer(4, GL_FLOAT, 0, 0);

  glClientActiveTexture(GL_TEXTURE1);
  glBindBuffer(GL_ARRAY_BUFFER, _mesh->surfaceVertNeighbor1VboID());
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);
  glTexCoordPointer(4, GL_FLOAT, 0, 0);

  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _tetRotationFboID);

  glClear(GL_COLOR_BUFFER_BIT);

  //glBlendFunc(GL_ONE, GL_ONE);

  //glEnable(GL_BLEND);
  //glDisable(GL_DEPTH_TEST);

  GLfloat offset[2];
  offset[0] = (0.5 / (double)_tetVertexTexWidth);
  offset[1] = (0.5 / (double)_tetVertexTexHeight);

  // Calculate normals and tangents for each tet vertex. 
  // Done in three passes, since each pass draws to one vertex position of each face
  /*if (_tetNormTanTex0Shader != NULL)
  {
    glUseProgram(_tetNormTanTex0Shader->getHandle());
    glUniform1i(_tetNormTanTex0Shader->uniformLoc("tetMeshTex"), 0);
    glUniform2f(_tetNormTanTex0Shader->uniformLoc("offset"), offset[0], offset[1]);
  }
  glDrawArrays(GL_POINTS, 0, _mesh->tetFaceSize());

  if (_tetNormTanTex1Shader != NULL)
  {
    glUseProgram(_tetNormTanTex1Shader->getHandle());
    glUniform1i(_tetNormTanTex1Shader->uniformLoc("tetMeshTex"), 0);
    glUniform2f(_tetNormTanTex1Shader->uniformLoc("offset"), offset[0], offset[1]);
  }
  glDrawArrays(GL_POINTS, 0, _mesh->tetFaceSize());

  if (_tetNormTanTex2Shader != NULL)
  {
    glUseProgram(_tetNormTanTex2Shader->getHandle());
    glUniform1i(_tetNormTanTex2Shader->uniformLoc("tetMeshTex"), 0);
    glUniform2f(_tetNormTanTex2Shader->uniformLoc("offset"), offset[0], offset[1]);
  }
  glDrawArrays(GL_POINTS, 0, _mesh->tetFaceSize());*/

  if (_mesh->genRotationShader() != NULL)
  {
    glUseProgram(_mesh->genRotationShader()->getHandle());
    glUniform1i(_mesh->genRotationShader()->uniformLoc("tetMeshTex"), 0);
    glUniform1i(_mesh->genRotationShader()->uniformLoc("n1RestFrame"), 1);
    glUniform1i(_mesh->genRotationShader()->uniformLoc("n2RestFrame"), 2);
    glUniform2f(_mesh->genRotationShader()->uniformLoc("offset"),
                offset[0], offset[1]);
  }
  glDrawArrays(GL_POINTS, 0, _mesh->surfaceVertSize());

  //glEnable(GL_DEPTH_TEST);
  //glDisable(GL_BLEND);
  
  glDisableClientState(GL_TEXTURE_COORD_ARRAY);
  glClientActiveTexture(GL_TEXTURE0);
  glDisableClientState(GL_TEXTURE_COORD_ARRAY);
  glDisableClientState(GL_VERTEX_ARRAY);

  //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0); 

  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
  if (_mesh->genRotationShader() != NULL)
  {
    glUseProgram(0);
  }

  glBindTexture(GL_TEXTURE_2D, 0);
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D, 0);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, 0);
  
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();

  glViewport(_windowViewport[0], _windowViewport[1],
             _windowViewport[2], _windowViewport[3]);
}

//////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::oldRenderAOVEmbedded(float delta, float width, float height)
{
  // Bind required vertex buffers
  glBindBuffer(GL_ARRAY_BUFFER, _mesh->embeddedMeshVboID());
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(3, GL_FLOAT, 0, 0);

   // Don't need to do a glEnable(GL_TEXTURE_2D) when using shaders
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, _tetVertexTexID);

  glClientActiveTexture(GL_TEXTURE0);
  glBindBuffer(GL_ARRAY_BUFFER, _mesh->embeddedVertLoc0VboID());
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);
  glTexCoordPointer(4, GL_FLOAT, 0, 0);

  glClientActiveTexture(GL_TEXTURE1);
  glBindBuffer(GL_ARRAY_BUFFER, _mesh->embeddedVertLoc1VboID());
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);
  glTexCoordPointer(4, GL_FLOAT, 0, 0);

  assert(_mesh->drawAOV2Shader() != NULL);
  
  glUseProgram(_mesh->drawAOV2Shader()->getHandle());
  glUniform1i(_mesh->drawAOV2Shader()->uniformLoc("tetMeshTex"), 0);
  glUniform1f(_mesh->drawAOV2Shader()->uniformLoc("delta"), delta);
  
  //cout << "screen size: " << _windowViewport[2] << " x " << _windowViewport[3] << endl;
  glUniform1i(_mesh->drawAOV2Shader()->uniformLoc("normalTex"), 1);
  glUniform1i(_mesh->drawAOV2Shader()->uniformLoc("positionTex"), 2);
  glUniform2f(_mesh->drawAOV2Shader()->uniformLoc("screenSize"), width, height);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _mesh->embeddedMeshIboID());
  

  if (startPoint < 0)
    startPoint += _mesh->tetFaceSize();
  if (startPoint >= _mesh->tetFaceSize())
    startPoint -= _mesh->tetFaceSize();

//  glDrawArrays(GL_POINTS, 0, _mesh->tetFaceSize());
  //glDrawArrays(GL_POINTS, startPoint, 20);
  //glDrawElements(GL_TRIANGLES, 27, GL_UNSIGNED_INT, (GLvoid*)(3*startPoint*sizeof(GLint)));
  glDrawElements(GL_TRIANGLES, _mesh->embeddedMeshIboSize(), GL_UNSIGNED_INT, 0);

  
  glDisableClientState(GL_TEXTURE_COORD_ARRAY);
  glClientActiveTexture(GL_TEXTURE0);
  glDisableClientState(GL_TEXTURE_COORD_ARRAY);

  glDisableClientState(GL_VERTEX_ARRAY);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0); 

  glUseProgram(0);

  glBindTexture(GL_TEXTURE_2D, 0);
}



//////////////////////////////////////////////////////////////////////
// Set the world position of the instance
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::setPosition(GLfloat *position)
{
  this->setPosition(position[0], position[1], position[2]);
}

void GLSL_TET_MESH_INSTANCE::setPosition(GLfloat x, GLfloat y, GLfloat z)
{
  _position[0] = x;
  _position[1] = y;
  _position[2] = z;
}

//////////////////////////////////////////////////////////////////////
// Translate the world position by the vector
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::addPosition(GLfloat *delta)
{
  this->addPosition(delta[0], delta[1], delta[2]);
}

void GLSL_TET_MESH_INSTANCE::addPosition(GLfloat dx, GLfloat dy, GLfloat dz)
{
  _position[0] += dx;
  _position[1] += dy;
  _position[2] += dz;
}


//////////////////////////////////////////////////////////////////////
// Set three rotation angles
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::setRotation(GLfloat *angles)
{ this->setRotation(angles[0], angles[1], angles[2]);  }

void GLSL_TET_MESH_INSTANCE::setRotation(GLfloat rx, GLfloat ry, GLfloat rz)
{
  _rotation[0] = rx;
  _rotation[1] = ry;
  _rotation[2] = rz;
}

//////////////////////////////////////////////////////////////////////
// Change three rotation angles, given delta angles
//////////////////////////////////////////////////////////////////////
void GLSL_TET_MESH_INSTANCE::addRotation(GLfloat *delta)
{ this->addRotation(delta[0], delta[1], delta[0]);    }

void GLSL_TET_MESH_INSTANCE::addRotation(GLfloat drx, GLfloat dry, GLfloat drz)
{
  _rotation[0] += drx;
  _rotation[1] += dry;
  _rotation[2] += drz;
}


