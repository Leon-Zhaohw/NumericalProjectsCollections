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
// GLSL_SCENE_VIEWER.h: Interface for the GLSL_SCENE_VIEWER class.
//
///////////

#ifndef GLSL_SCENE_VIEWER_H
#define GLSL_SCENE_VIEWER_H

#include <GLSL_LIGHT.h>
#include <GLSL_TET_MESH.h>
#include <GLSL_TET_MESH_INSTANCE.h>
#include <GLSL_OBJ.h>
#include <VEC3.h>
#include <SIMPLE_PARSER.h>
#include <SUBSPACE_TET_MESH.h>
#include <OBJ.h>

#include <string>
#include <vector>
#include <map>

#define MAX_LIGHTS 8

class GLSL_SCENE_VIEWER
{
public:
  GLSL_SCENE_VIEWER();
  virtual ~GLSL_SCENE_VIEWER();

  void init();
  void update( bool debug = false );
  void stepSystem( bool debug = false );
  void display();

  int addMesh(string cfgFilename);
  GLSL_TET_MESH_INSTANCE* getMesh(int meshid);

  GLSL_LIGHT *getLight(int index);

  void addObj(OBJ *obj);
  void addObj(OBJ *obj, GLfloat r, GLfloat g, GLfloat b);

  bool mouseClick(VEC3F& point, Real maxRadius = -1.0);
  void mouseDrag(VEC3F& point);
  void mouseUnclick();

  void setDrawScene(bool drawScene)               { _drawScene = drawScene; };
  void setUseAmbientOcclusion(bool useAO)         { _useAO = useAO;         };
  void setClearColor(GLfloat r, GLfloat g, GLfloat b, GLfloat a)
        { _clearColor[0] = r;   _clearColor[1] = g;   _clearColor[2] = b;   _clearColor[3] = a;   };

  // Methods to modify uniform variables used by shadowing
  void changeFactor(GLfloat factor);
  void changeUnits(GLfloat units); 

  // Used for debugging
  void changeStartPoint(int delta); 

  void setGroundHeight(GLfloat height);
  void setAOVDelta(GLfloat delta);

protected:
  map<string, int> _cfgIndexMap;
  vector<GLSL_TET_MESH*> _meshList;
  vector<GLSL_TET_MESH_INSTANCE*> _instanceList;
  vector<GLSL_OBJ*> _objList;

  bool _drawScene;
  bool _useAO;
  GLfloat _clearColor[4];

  GLuint _shadowMapFboID;
  GLuint _shadowMapRboID;
  GLsizei _shadowMapDimensions;

  GLuint _deferredTex0ID;
  GLuint _deferredTex1ID;
  GLuint _deferredFboID;
  GLuint _deferredRboID;

  GLuint _aovAccumTexID;
  GLuint _aovAccumFboID;
  GLuint _aovAccumRboID;
  GLuint _aovDeferredTex0ID;
  GLuint _aovDeferredTex1ID;
  GLuint _aovDeferredFboID;
  GLint  _aovAccumWidth;
  GLint  _aovAccumHeight;
  GLfloat _aovDelta;

  GLSL_SHADER *_bilateralShader;
  
  GLSL_TET_MESH_INSTANCE* _animatedMesh;

  GLSL_SHADER *_embeddedMeshRenderShader;

  GLSL_SHADER *_embeddedMeshToTextureShader;

  GLSL_SHADER *_meshRenderWithShadowShader;
  GLSL_SHADER *_meshRenderNoShadowShader;
  GLSL_SHADER *_meshRenderSMShader;
  GLSL_SHADER *_meshRenderDeferredShader;
  GLSL_SHADER *_meshRenderBasicShader;

  GLSL_SHADER *_meshRenderNormalsShader;

  GLSL_SHADER *_sceneRenderWithShadowShader;
  GLSL_SHADER *_sceneRenderNoShadowShader;
  GLSL_SHADER *_sceneRenderSMShader;
  GLSL_SHADER *_sceneRenderDeferredShader;
  GLSL_SHADER *_sceneRenderBasicShader;
  GLSL_SHADER *_sceneRenderAOVShader;

  GLSL_SHADER *_tetNormTanTex0Shader;
  GLSL_SHADER *_tetNormTanTex1Shader;
  GLSL_SHADER *_tetNormTanTex2Shader;

  GLSL_SHADER *_shadowBlurShader;
  GLuint _shadowBlurTempTexID;

  GLSL_SHADER *_ssoaRenderShader;
  GLuint _ssaoTexID;
  GLuint _ssBlurTexID;

  GLSL_SHADER *_drawTextureShader;

  GLSL_LIGHT *_lights[MAX_LIGHTS];

  // Vertex buffer of the ground plane
  GLuint _groundPlaneVboID;
  GLfloat _groundHeight;
    
  // Shadow map generation uniform variable values
  GLfloat _minVariance;
  GLfloat _reduceBleed;

  // window viewport size
  GLint _windowViewport[4];

  float _lastForceMultiplier;

  MATERIAL* readMaterial(SIMPLE_PARSER& parser);
  bool loadTetMesh(string cfgFilename, SUBSPACE_TET_MESH **tetMesh, OBJ **embeddedMesh);

  bool reportFboStatus(GLenum fboStatus);
  bool reportGLError(GLenum error);

  void drawGroundPlane(int style);
  void blurShadowMaps();
  void blurTex(GLuint src, GLuint blur, GLsizei dim);
  
  void drawToDeferred();
  void drawSSAO();
};

#endif // GLSL_SCENE_VIEWER_H

