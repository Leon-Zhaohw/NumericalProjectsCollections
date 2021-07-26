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
// GLSL_SCENE_VIEWER.cpp: Interface for the GLSL_SCENE_VIEWER class.
//
//////////////////////////////////////////////////////////////////////

#include <GLSL_SCENE_VIEWER.h>
#include <STVK.h>
#include <MOONEY_RIVLIN.h>
#include <ARRUDA_BOYCE.h>
#include <NEO_HOOKEAN.h>
#include <INVERTIBLE.h>

#define GROUND_DRAW_STYLE_SHADOW    1
#define GROUND_DRAW_STYLE_NO_SHADOW 2
#define GROUND_DRAW_STYLE_SM        3
#define GROUND_DRAW_STYLE_DEFERRED  4
#define GROUND_DRAW_STYLE_AOV       5
#define GROUND_DRAW_STYLE_BASIC     6

#define OBJ_MOVE_DOWN (-0.15f)

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
GLSL_SCENE_VIEWER::GLSL_SCENE_VIEWER()
{
  _cfgIndexMap.clear();
  _meshList.clear();
  _instanceList.clear();
  _objList.clear();

  _drawScene = false;
  _useAO = false;

  _animatedMesh = NULL;
  _embeddedMeshRenderShader = NULL;
  _embeddedMeshToTextureShader = NULL;
  _meshRenderWithShadowShader = NULL;
  _meshRenderNoShadowShader = NULL;
  _meshRenderSMShader = NULL;
  _meshRenderDeferredShader = NULL;
  _meshRenderBasicShader = NULL;

  _meshRenderNormalsShader = NULL;

  _sceneRenderWithShadowShader = NULL;
  _sceneRenderNoShadowShader = NULL;
  _sceneRenderSMShader = NULL;
  _sceneRenderDeferredShader = NULL;
  _sceneRenderBasicShader = NULL;
  _sceneRenderAOVShader = NULL;

  _tetNormTanTex0Shader = NULL;
  _tetNormTanTex1Shader = NULL;
  _tetNormTanTex2Shader = NULL;

  _shadowBlurShader = NULL;
  _shadowBlurTempTexID = 0;

  _ssoaRenderShader = NULL;
  _ssaoTexID = 0;
  _ssBlurTexID = 0;

  for (int x = 0; x < MAX_LIGHTS; x++)
    _lights[x] = NULL;

  _groundPlaneVboID = 0;

  _aovDelta = 0.1f;
  
  _minVariance = 0.001f;
  _reduceBleed = 0.65f;

  _shadowMapDimensions = 1024;
}


//////////////////////////////////////////////////////////////////////
void GLSL_SCENE_VIEWER::changeFactor(GLfloat factor)
{ 
  _minVariance += (0.01f * factor); 
  if (_minVariance < 0.0f)
    _minVariance = 0.0f;
  cout << " Min Variance: " << _minVariance << endl;
}

void GLSL_SCENE_VIEWER::changeUnits(GLfloat units)
{ 
  _reduceBleed += units;
  if (_reduceBleed < 0.0f)
    _reduceBleed = 0.0f;
  cout << " Reduce Bleed: " << _reduceBleed << endl;
}

void GLSL_SCENE_VIEWER::changeStartPoint(int delta)
{
  for (int x = 0; x < (int)_instanceList.size(); x++)
  {
    _instanceList[x]->startPoint += delta;
    cout << " startpoint " << _instanceList[x]->startPoint << endl;
  }
}


//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
GLSL_SCENE_VIEWER::~GLSL_SCENE_VIEWER()
{
  // Clean up mesh instance objects
  for (int x = 0; x < (int)_instanceList.size(); x++)
  {
    delete _instanceList[x];
  }
  _instanceList.clear();
  _cfgIndexMap.clear();

  // Clean up mesh objects
  for (int x = 0; x < (int)_meshList.size(); x++)
  {
    delete _meshList[x];
  }
  _meshList.clear();

  // Clean up GLSL_OBJ objects
  for (int x = 0; x < (int)_objList.size(); x++)
  {
    delete _objList[x];
  }
  _objList.clear();

  if (_shadowMapFboID != 0)  {
    glDeleteFramebuffersEXT(1, &_shadowMapFboID);
    _shadowMapFboID = 0;
  }
  if (_shadowMapRboID != 0)  {
    glDeleteRenderbuffersEXT(1, &_shadowMapRboID);
    _shadowMapRboID = 0;
  }

  if (_deferredFboID != 0)
  {
    glDeleteFramebuffersEXT(1, &_deferredFboID);
    _deferredFboID = 0;
  }
  if (_deferredRboID != 0)  {
    glDeleteRenderbuffersEXT(1, &_deferredRboID);
    _deferredRboID = 0;
  }
  if (_deferredTex0ID != 0)  {
    glDeleteTextures(1, &_deferredTex0ID);
    _deferredTex0ID = 0;
  }
  if (_deferredTex1ID != 0)  {
    glDeleteTextures(1, &_deferredTex1ID);
    _deferredTex1ID = 0;
  }

  if (_aovAccumRboID != 0)  {
    glDeleteRenderbuffersEXT(1, &_aovAccumRboID);
    _aovAccumRboID = 0;
  }
  if (_aovAccumFboID != 0)
  {
    glDeleteFramebuffersEXT(1, &_aovAccumFboID);
    _aovAccumFboID = 0;
  }
  if (_aovDeferredFboID != 0)
  {
    glDeleteFramebuffersEXT(1, &_aovDeferredFboID);
    _aovDeferredFboID = 0;
  }

  if (_aovAccumTexID != 0)  {
    glDeleteTextures(1, &_aovAccumTexID);
    _aovAccumTexID = 0;
  }
  if (_aovAccumTexID != 0)  {
    glDeleteTextures(1, &_aovDeferredTex0ID);
    _aovAccumTexID = 0;
  }
  if (_aovAccumTexID != 0)  {
    glDeleteTextures(1, &_aovDeferredTex1ID);
    _aovAccumTexID = 0;
  }

  if (_bilateralShader != NULL)  {
    delete _bilateralShader;
    _bilateralShader = NULL;
  }

  _animatedMesh = NULL;

  // Clean up shader objects
  if (_embeddedMeshRenderShader != NULL)  {
    delete _embeddedMeshRenderShader;
    _embeddedMeshRenderShader = NULL;
  }
  if (_embeddedMeshToTextureShader != NULL)  {
    delete _embeddedMeshToTextureShader;
    _embeddedMeshToTextureShader = NULL;
  }
  if (_meshRenderWithShadowShader != NULL)  {
    delete _meshRenderWithShadowShader;
    _meshRenderWithShadowShader = NULL;
  }
  if (_meshRenderNoShadowShader != NULL)  {
    delete _meshRenderNoShadowShader;
    _meshRenderNoShadowShader = NULL;
  }
  if (_meshRenderSMShader != NULL)  {
    delete _meshRenderSMShader;
    _meshRenderSMShader = NULL;
  }
  if (_meshRenderDeferredShader != NULL)  {
    delete _meshRenderDeferredShader;
    _meshRenderDeferredShader = NULL;
  }
  if (_meshRenderBasicShader != NULL) {
    delete _meshRenderBasicShader;
    _meshRenderBasicShader = NULL;
  }
  if (_meshRenderNormalsShader != NULL)  {
    delete _meshRenderNormalsShader;
    _meshRenderNormalsShader = NULL;
  }

  if (_sceneRenderWithShadowShader != NULL)  {
    delete _sceneRenderWithShadowShader;
    _sceneRenderWithShadowShader = NULL;
  }
  if (_sceneRenderNoShadowShader != NULL)  {
    delete _sceneRenderNoShadowShader;
    _sceneRenderNoShadowShader = NULL;
  }
  if (_sceneRenderSMShader != NULL)  {
    delete _sceneRenderSMShader;
    _sceneRenderSMShader = NULL;
  }
  if (_sceneRenderDeferredShader != NULL)  {
    delete _sceneRenderDeferredShader;
    _sceneRenderDeferredShader = NULL;
  }
  if (_sceneRenderBasicShader != NULL) {
    delete _sceneRenderBasicShader;
    _sceneRenderBasicShader = NULL;
  }
  if (_sceneRenderAOVShader != NULL)  {
    delete _sceneRenderAOVShader;
    _sceneRenderAOVShader = NULL;
  }
  
  if (_tetNormTanTex0Shader != NULL)  {
    delete _tetNormTanTex0Shader;
    _tetNormTanTex0Shader = NULL;
  }
  if (_tetNormTanTex1Shader != NULL)  {
    delete _tetNormTanTex1Shader;
    _tetNormTanTex1Shader = NULL;
  }
  if (_tetNormTanTex2Shader != NULL)  {
    delete _tetNormTanTex2Shader;
    _tetNormTanTex2Shader = NULL;
  }

  if (_shadowBlurShader != NULL)  {
    delete _shadowBlurShader;
    _shadowBlurShader = NULL;
  }
  if (_shadowBlurTempTexID != 0)  {
    glDeleteTextures(1, &_shadowBlurTempTexID);
    _shadowBlurTempTexID = 0;
  }
  
  if (_ssoaRenderShader != NULL)  {
    delete _ssoaRenderShader;
    _ssoaRenderShader = NULL;
  }
  if (_ssaoTexID != 0)  {
    glDeleteTextures(1, &_ssaoTexID);
    _ssaoTexID = 0;
  }
  if (_ssBlurTexID != 0)  {
    glDeleteTextures(1, &_ssBlurTexID);
    _ssBlurTexID = 0;
  }

  if (_drawTextureShader != NULL) {
    delete _drawTextureShader;
    _drawTextureShader = NULL;
  }

  if (_groundPlaneVboID != 0)  {
    glDeleteBuffers(1, &_groundPlaneVboID);
    _groundPlaneVboID = 0;
  }

  for (int x = 0; x < MAX_LIGHTS; x++)
  {
    if (_lights[x] != NULL)
    {
      delete _lights[x];
      _lights[x] = NULL;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Initialize the viewer
//////////////////////////////////////////////////////////////////////
void GLSL_SCENE_VIEWER::init()
{
  // Initialize OpenGL extension wrangler, so we can use fancy GL features
  GLenum err = glewInit();
  if (GLEW_OK != err)
  {
    // Problem: glewInit failed, something is seriously wrong. */
    //cout << " Error occurred when initializing GLEW" << endl;
    fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
  } else {
    cout << " GLEW was initialized correctly" << endl;
  }
  //_clearColor[0] = _clearColor[1] = _clearColor[2] = 0.0f;
  _clearColor[0] = 1.0f;
  _clearColor[1] = 1.0f;
  _clearColor[2] = 1.0f;
  _clearColor[3] = 0.0f;

  // Init OpenGL default states
  glClearColor(_clearColor[0],_clearColor[1],_clearColor[2],_clearColor[3]);
	glClearDepth(1.0f);
	//glClearDepth(0.0f);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_DEPTH_TEST);
	glFrontFace(GL_CCW);
	glEnable(GL_CULL_FACE);
	glShadeModel(GL_SMOOTH);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	glHint(GL_GENERATE_MIPMAP_HINT, GL_NICEST);

  glColorMaterial	(GL_FRONT_AND_BACK , GL_AMBIENT_AND_DIFFUSE );
  glEnable(GL_COLOR_MATERIAL);

  GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
  
  int resolutionFraction = 1;
  //_aovDelta = 0.25f;
  
  // Start loading up all GLSL shaders
  
  // Upsample that takes a deferred framebuffer texture to perform a geometry aware blur
  _bilateralShader = new GLSL_SHADER();
  if (_bilateralShader != NULL)
  {
    _bilateralShader->attachVert("src/rendering/glsl/aovup.vert");
    _bilateralShader->attachFrag("src/rendering/glsl/aovup2.frag");
    if (resolutionFraction > 1)
    {
      _bilateralShader->setProgramConst("SAMPLE_RADIUS", (float)resolutionFraction, 2);
      _bilateralShader->setProgramConst("_USE_FILTER_", 1, 2);
    }
    _bilateralShader->compile();
  }

  _embeddedMeshToTextureShader = new GLSL_SHADER();
  if (_embeddedMeshToTextureShader != NULL)
  {
    _embeddedMeshToTextureShader->attachVert("src/rendering/glsl/embedded16.vert");
    _embeddedMeshToTextureShader->attachFrag("src/rendering/glsl/embedded16.frag");
    _embeddedMeshToTextureShader->attachGeom("src/rendering/glsl/embedded16.geom");
    _embeddedMeshToTextureShader->setGeomInput(GL_TRIANGLES);
    _embeddedMeshToTextureShader->setGeomOutput(GL_POINTS);
    _embeddedMeshToTextureShader->setGeomVerticesOut(3);
    _embeddedMeshToTextureShader->compile();
  }

  // Render an embedded mesh of a subspace tetmesh with variance shadow mapping
  _meshRenderWithShadowShader = new GLSL_SHADER();
  if (_meshRenderWithShadowShader != NULL)
  {
    _meshRenderWithShadowShader->attachVert("src/rendering/glsl/embedded17.vert");
    _meshRenderWithShadowShader->attachFrag("src/rendering/glsl/embedded18.frag");
    _meshRenderWithShadowShader->compile();
  }

  // Render an embedded mesh of a subspace tetmesh with basic NdotL diffuse lighting
  _meshRenderNoShadowShader = new GLSL_SHADER();
  if (_meshRenderNoShadowShader != NULL)
  {
    _meshRenderNoShadowShader->attachVert("src/rendering/glsl/embedded17.vert");
    _meshRenderNoShadowShader->attachFrag("src/rendering/glsl/embedded17.frag");
    //_meshRenderNoShadowShader->attachVert("src/rendering/glsl/embedded17_3light.vert");
    //_meshRenderNoShadowShader->attachFrag("src/rendering/glsl/embedded17_3light.frag");
    // setProgramConst(constName, constValue, whichShader)
    //_meshRenderNoShadowShader->setProgramConst("TOTAL_LIGHTS", 1, 3);
    _meshRenderNoShadowShader->compile();
  }

  // Render an embedded mesh of a subspace tetmesh with the glColor set
  _meshRenderSMShader = new GLSL_SHADER();
  if (_meshRenderSMShader != NULL)
  {
    _meshRenderSMShader->attachVert("src/rendering/glsl/embedded19.vert");
    _meshRenderSMShader->attachFrag("src/rendering/glsl/embedded19.frag");
    _meshRenderSMShader->compile();

    //glUseProgram(_meshRenderSMShader->getHandle());
    //glUniform1f(_meshRenderSMShader->uniformLoc("factor"), _smFactor);
    //glUniform1f(_meshRenderSMShader->uniformLoc("units"), _smUnits);
    //glUseProgram(0);
  }

  // Render an embedded mesh of a subspace tetmesh into a deferred
  // framebuffer (i.e. frag[0] = normal,depth, frag[1] = eye space pos)
  _meshRenderDeferredShader = new GLSL_SHADER();
  if (_meshRenderDeferredShader != NULL)
  {
    _meshRenderDeferredShader->attachVert("src/rendering/glsl/deferred03.vert");
    _meshRenderDeferredShader->attachFrag("src/rendering/glsl/deferred03.frag");
    _meshRenderDeferredShader->compile();
  }

  _meshRenderBasicShader = new GLSL_SHADER();
  if (_meshRenderBasicShader != NULL)
  {
    _meshRenderBasicShader->attachVert("src/rendering/glsl/embedded20.vert");
    _meshRenderBasicShader->attachFrag("src/rendering/glsl/embedded20.frag");
    _meshRenderBasicShader->compile();
  }

  // Render the vertex normals of an embedded mesh
  _meshRenderNormalsShader = new GLSL_SHADER();
  if (_meshRenderNormalsShader != NULL)
  {
    _meshRenderNormalsShader->attachVert("src/rendering/glsl/embedded13.vert");
    _meshRenderNormalsShader->attachFrag("src/rendering/glsl/embedded13.frag");
    _meshRenderNormalsShader->attachGeom("src/rendering/glsl/embedded13.geom");
    _meshRenderNormalsShader->setGeomInput(GL_POINTS);
    _meshRenderNormalsShader->setGeomOutput(GL_LINE_STRIP);
    _meshRenderNormalsShader->setGeomVerticesOut(2);
    _meshRenderNormalsShader->compile();
  }

  _sceneRenderWithShadowShader = new GLSL_SHADER();
  if (_sceneRenderWithShadowShader != NULL)
  {
    _sceneRenderWithShadowShader->attachVert("src/rendering/glsl/ground01.vert");
    _sceneRenderWithShadowShader->attachFrag("src/rendering/glsl/embedded18.frag");
    _sceneRenderWithShadowShader->compile();
  }

  _sceneRenderNoShadowShader = new GLSL_SHADER();
  if (_sceneRenderNoShadowShader != NULL)
  {
    _sceneRenderNoShadowShader->attachVert("src/rendering/glsl/ground01.vert");
    _sceneRenderNoShadowShader->attachFrag("src/rendering/glsl/embedded17.frag");
    _sceneRenderNoShadowShader->compile();
  }

  _sceneRenderSMShader = new GLSL_SHADER();
  if (_sceneRenderSMShader != NULL)
  {
    _sceneRenderSMShader->attachVert("src/rendering/glsl/ground06.vert");
    _sceneRenderSMShader->attachFrag("src/rendering/glsl/embedded19.frag");
    _sceneRenderSMShader->compile();
  }

  _sceneRenderDeferredShader = new GLSL_SHADER();
  if (_sceneRenderDeferredShader != NULL)
  {
    _sceneRenderDeferredShader->attachVert("src/rendering/glsl/ground03.vert");
    _sceneRenderDeferredShader->attachGeom("src/rendering/glsl/deferred03.geom");
    _sceneRenderDeferredShader->attachFrag("src/rendering/glsl/deferred03.frag");
    _sceneRenderDeferredShader->setGeomInput(GL_TRIANGLES);
    _sceneRenderDeferredShader->setGeomOutput(GL_TRIANGLE_STRIP);
    _sceneRenderDeferredShader->setGeomVerticesOut(3);
    _sceneRenderDeferredShader->compile();
  }

  _sceneRenderBasicShader = new GLSL_SHADER();
  if (_sceneRenderBasicShader != NULL)
  {
    _sceneRenderBasicShader->attachVert("src/rendering/glsl/ground07.vert");
    _sceneRenderBasicShader->attachFrag("src/rendering/glsl/embedded20.frag");
    _sceneRenderBasicShader->compile();
  }

  _sceneRenderAOVShader = new GLSL_SHADER();
  if (_sceneRenderAOVShader != NULL)
  {
    _sceneRenderAOVShader->attachVert("src/rendering/glsl/ground04.vert");
    _sceneRenderAOVShader->attachFrag("src/rendering/glsl/ground04.frag");
    _sceneRenderAOVShader->compile();
  }// */
  /*_sceneRenderAOVShader = new GLSL_SHADER();
  if (_sceneRenderAOVShader != NULL)
  {
    _sceneRenderAOVShader->attachVert("src/rendering/glsl/ground05.vert");
    _sceneRenderAOVShader->attachGeom("src/rendering/glsl/aov2.geom");
    _sceneRenderAOVShader->attachFrag("src/rendering/glsl/aov1.frag");
    _sceneRenderAOVShader->setGeomInput(GL_TRIANGLES);
    _sceneRenderAOVShader->setGeomOutput(GL_TRIANGLE_STRIP);
    _sceneRenderAOVShader->setGeomVerticesOut(16);
    _sceneRenderAOVShader->compile();
  }// */

  /*_sceneRenderAOVShader = new GLSL_SHADER();
  if (_sceneRenderAOVShader != NULL)
  {
    _sceneRenderAOVShader->attachVert("src/rendering/glsl/ground05.vert");
    _sceneRenderAOVShader->attachGeom("src/rendering/glsl/aov3.geom");
    _sceneRenderAOVShader->attachFrag("src/rendering/glsl/aov3.frag");
    _sceneRenderAOVShader->setGeomInput(GL_TRIANGLES);
    _sceneRenderAOVShader->setGeomOutput(GL_LINE_STRIP);
    _sceneRenderAOVShader->setGeomVerticesOut(26);
    _sceneRenderAOVShader->compile();
  }// */

  _tetNormTanTex0Shader = new GLSL_SHADER();
  if (_tetNormTanTex0Shader != NULL)
  {
    _tetNormTanTex0Shader->attachVert("src/rendering/glsl/tetnormtan.vert");
    _tetNormTanTex0Shader->attachFrag("src/rendering/glsl/tetnormtan.frag");
    _tetNormTanTex0Shader->setProgramConst("TARGET", 0, 1);
    _tetNormTanTex0Shader->compile();
  }
  _tetNormTanTex1Shader = new GLSL_SHADER();
  if (_tetNormTanTex1Shader != NULL)
  {
    _tetNormTanTex1Shader->attachVert("src/rendering/glsl/tetnormtan.vert");
    _tetNormTanTex1Shader->attachFrag("src/rendering/glsl/tetnormtan.frag");
    _tetNormTanTex1Shader->setProgramConst("TARGET", 1, 1);
    _tetNormTanTex1Shader->compile();
  }
  _tetNormTanTex2Shader = new GLSL_SHADER();
  if (_tetNormTanTex2Shader != NULL)
  {
    _tetNormTanTex2Shader->attachVert("src/rendering/glsl/tetnormtan.vert");
    _tetNormTanTex2Shader->attachFrag("src/rendering/glsl/tetnormtan.frag");
    _tetNormTanTex2Shader->setProgramConst("TARGET", 2, 1);
    _tetNormTanTex2Shader->compile();
  }

  _ssoaRenderShader = new GLSL_SHADER();
  if (_ssoaRenderShader != NULL)
  {
    _ssoaRenderShader->attachVert("src/rendering/glsl/ssao01.vert");
    _ssoaRenderShader->attachFrag("src/rendering/glsl/ssao01.frag");
    _ssoaRenderShader->compile();
  }

  _shadowBlurShader = new GLSL_SHADER();
  if (_shadowBlurShader != NULL)
  {
    _shadowBlurShader->attachVert("src/rendering/glsl/shadowblur.vert");
    _shadowBlurShader->attachFrag("src/rendering/glsl/shadowblur5.frag");
    _shadowBlurShader->compile();
  }

  _drawTextureShader = new GLSL_SHADER();
  if (_drawTextureShader != NULL)
  {
    _drawTextureShader->attachVert("src/rendering/glsl/drawtexture.vert");
    _drawTextureShader->attachFrag("src/rendering/glsl/drawtexture.frag");
    _drawTextureShader->compile();
  }

  // Ground vertices:     Position            Normal              Color
  GLfloat ground[36] = { -5.0f, 5.0f, 0.0f,   0.0f, 0.0f, 1.0f,   0.7f, 0.7f, 0.7f,
                         -5.0f,-5.0f, 0.0f,   0.0f, 0.0f, 1.0f,   0.7f, 0.7f, 0.7f,
                          5.0f, 5.0f, 0.0f,   0.0f, 0.0f, 1.0f,   0.7f, 0.7f, 0.7f,
                          5.0f,-5.0f, 0.0f,   0.0f, 0.0f, 1.0f,   0.7f, 0.7f, 0.7f  };
  /*
  // Ground vertices:     Position            Normal              Color
  GLfloat ground[36] = { -4.5f, 0.0f, 5.5f,   0.0f, 1.0f, 0.0f,   0.7f, 0.7f, 0.7f,
                         -4.5f, 0.0f, -4.5f,   0.0f, 1.0f, 0.0f,   0.7f, 0.7f, 0.7f,
                          5.5f, 0.0f, 5.5f,   0.0f, 1.0f, 0.0f,   0.7f, 0.7f, 0.7f,
                          5.5f, 0.0f,-4.5f,   0.0f, 1.0f, 0.0f,   0.7f, 0.7f, 0.7f  };
                          */

  glGenBuffers(1, &_groundPlaneVboID);
  if (_groundPlaneVboID == 0)
  {
    cout << "Could not generate a GL buffer object" << endl;
  } else {
	  // Attempt to copy to GPU.
    glBindBuffer(GL_ARRAY_BUFFER, _groundPlaneVboID);
    glBufferData(GL_ARRAY_BUFFER, 4 * 9 * sizeof(GLfloat),
                 (GLvoid*)ground, GL_STATIC_DRAW);
		
    glBindBuffer(GL_ARRAY_BUFFER, 0);  
  }
  _groundHeight = 0.0f;

  glGenTextures(1, &_shadowBlurTempTexID);
  if (_shadowBlurTempTexID == 0)
    cout << " Error occurred generating shadow blur temp texture object for shadow map fbo" << endl;
  glBindTexture(GL_TEXTURE_2D, _shadowBlurTempTexID);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F_ARB, _shadowMapDimensions,
               _shadowMapDimensions, 0, GL_RGB, GL_FLOAT, 0);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  
  glBindTexture(GL_TEXTURE_2D, 0);

  glGenFramebuffersEXT(1, &_shadowMapFboID);
  if (_shadowMapFboID == 0)  {
    cout << " Could not generate shadow mapping framebuffer object." << endl;
  }
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _shadowMapFboID);
  //GLenum smbuf[1] = { GL_COLOR_ATTACHMENT0_EXT };
  //glDrawBuffers(1, smbuf);

  glGenRenderbuffersEXT(1, &_shadowMapRboID);
  if (_shadowMapRboID == 0)  {
    cout << " Could not generate shadow mapping renderbuffer object." << endl;
  }
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, _shadowMapRboID);
  glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT32,
                           _shadowMapDimensions, _shadowMapDimensions);
  glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT,
                               GL_RENDERBUFFER_EXT, _shadowMapRboID);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT,
                            GL_TEXTURE_2D, _shadowBlurTempTexID, 0);
 
  GLenum fboStatus = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
  reportFboStatus(fboStatus);

  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, 0);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);

  // Grab the viewport dimensions
  glGetIntegerv(GL_VIEWPORT, _windowViewport);
 
  glGenTextures(1, &_deferredTex0ID);
  if (_deferredTex0ID == 0)
    cout << " Error occurred generating deferred render target" << endl;
  glBindTexture(GL_TEXTURE_2D, _deferredTex0ID);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, _windowViewport[2],
               _windowViewport[3], 0, GL_RGBA, GL_FLOAT, 0);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  glGenTextures(1, &_deferredTex1ID);
  if (_deferredTex1ID == 0)
    cout << " Error occurred generating deferred render target" << endl;
  glBindTexture(GL_TEXTURE_2D, _deferredTex1ID);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, _windowViewport[2],
               _windowViewport[3], 0, GL_RGBA, GL_FLOAT, 0);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);


  glGenTextures(1, &_ssaoTexID);
  if (_ssaoTexID == 0)
    cout << " Error occurred generating deferred render target" << endl;
  glBindTexture(GL_TEXTURE_2D, _ssaoTexID);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, _windowViewport[2],
               _windowViewport[3], 0, GL_RGBA, GL_FLOAT, 0);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  glGenTextures(1, &_ssBlurTexID);
  if (_ssBlurTexID == 0)
    cout << " Error occurred generating deferred render target" << endl;
  glBindTexture(GL_TEXTURE_2D, _ssBlurTexID);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, _windowViewport[2],
               _windowViewport[3], 0, GL_RGBA, GL_FLOAT, 0);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  glBindTexture(GL_TEXTURE_2D, 0);

  glGenFramebuffersEXT(1, &_deferredFboID);
  if (_deferredFboID == 0)  {
    cout << " Could not generate shadow mapping framebuffer object." << endl;
  }
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _deferredFboID);
  GLenum defBuffer[2] = { GL_COLOR_ATTACHMENT0_EXT, GL_COLOR_ATTACHMENT1_EXT };
  glDrawBuffers(2, defBuffer);

  glGenRenderbuffersEXT(1, &_deferredRboID);
  if (_deferredRboID == 0)  {
    cout << " Could not generate shadow mapping renderbuffer object." << endl;
  }
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, _deferredRboID);
  glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT32,
                           _windowViewport[2], _windowViewport[3]);
  glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT,
                               GL_RENDERBUFFER_EXT, _deferredRboID);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT,
                            GL_TEXTURE_2D, _deferredTex0ID, 0);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT1_EXT,
                            GL_TEXTURE_2D, _deferredTex1ID, 0);
 
  fboStatus = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
  reportFboStatus(fboStatus);

  _aovAccumWidth = _windowViewport[2] / resolutionFraction;
  _aovAccumHeight = _windowViewport[3] / resolutionFraction;
  
  glGenTextures(1, &_aovAccumTexID);
  if (_aovAccumTexID == 0)
    cout << " Error occurred generating deferred render target" << endl;
  glBindTexture(GL_TEXTURE_2D, _aovAccumTexID);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB,
               _aovAccumWidth, _aovAccumHeight,
               0, GL_RGBA, GL_FLOAT, 0);
  //glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8,
                 //_aovAccumWidth, _aovAccumHeight, 0,
                 //GL_RGBA, GL_UNSIGNED_BYTE, 0);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  glGenTextures(1, &_aovDeferredTex0ID);
  if (_aovDeferredTex0ID == 0)
    cout << " Error occurred generating AOV deferred render target" << endl;
  glBindTexture(GL_TEXTURE_2D, _aovDeferredTex0ID);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB,
               _aovAccumWidth, _aovAccumHeight, 0, GL_RGBA, GL_FLOAT, 0);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  glGenTextures(1, &_aovDeferredTex1ID);
  if (_aovDeferredTex1ID == 0)
    cout << " Error occurred generating AOV deferred render target" << endl;
  glBindTexture(GL_TEXTURE_2D, _aovDeferredTex1ID);
  glTexImage2D(GL_TEXTURE_2D, 0,
               GL_RGBA32F_ARB, _aovAccumWidth,
               _aovAccumHeight, 0, GL_RGBA, GL_FLOAT, 0);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  
  glGenFramebuffersEXT(1, &_aovAccumFboID);
  if (_aovAccumFboID == 0)  {
    cout << " Could not generate AOV accumulation framebuffer object." << endl;
  }
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _aovAccumFboID);
  //GLenum accBuffer[1] = { GL_COLOR_ATTACHMENT0_EXT };
  //glDrawBuffers(1, accBuffer);
  
  glGenRenderbuffersEXT(1, &_aovAccumRboID);
  if (_aovAccumRboID == 0)  {
    cout << " Could not generate AOV accumulation renderbuffer object." << endl;
  }
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, _aovAccumRboID);
  glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT,
                           GL_DEPTH_COMPONENT32,
                           _aovAccumWidth,
                           _aovAccumHeight);
  glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT,
                               GL_DEPTH_ATTACHMENT_EXT,
                               GL_RENDERBUFFER_EXT,
                               _aovAccumRboID);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,
                            GL_COLOR_ATTACHMENT0_EXT,
                            GL_TEXTURE_2D,
                            _aovAccumTexID, 0);
   
  fboStatus = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
  reportFboStatus(fboStatus);

  glGenFramebuffersEXT(1, &_aovDeferredFboID);
  if (_aovDeferredFboID == 0)  {
    cout << " Could not generate AOV Deferred framebuffer object." << endl;
  }
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _aovDeferredFboID);
  GLenum accBuffer[2] = { GL_COLOR_ATTACHMENT0_EXT, GL_COLOR_ATTACHMENT1_EXT };
  glDrawBuffers(2, accBuffer);
  
  glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT,
                               GL_DEPTH_ATTACHMENT_EXT,
                               GL_RENDERBUFFER_EXT,
                               _aovAccumRboID);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,
                            GL_COLOR_ATTACHMENT0_EXT,
                            GL_TEXTURE_2D,
                            _aovDeferredTex0ID, 0);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,
                            GL_COLOR_ATTACHMENT1_EXT,
                            GL_TEXTURE_2D,
                            _aovDeferredTex1ID, 0);
   
  fboStatus = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
  reportFboStatus(fboStatus);
}

//////////////////////////////////////////////////////////////////////
// This will be where all mesh instances are updated, and shadow maps are drawn
//////////////////////////////////////////////////////////////////////
void GLSL_SCENE_VIEWER::update( bool debug )
{
  // Update/deform all subspace temtmeshes, and generate the tet vertex position textures
  for (int x = 0; x < (int)_instanceList.size(); x++)
  {
    _instanceList[x]->integratorStep( debug );
    _instanceList[x]->genTetVertexTexture();
  }
  //_instanceList[0]->setPosition(0.0f, 0.0f, 0.0f);
  
  for (int x = 0; x < (int)_objList.size(); x++)
  {
    // If OBJ object has changed, this needs to be called to
    // update the stored vertex arrays.
    _objList[x]->update();
  }

  //glUseProgram(_meshRenderSMShader->getHandle());
  //glUniform1f(_meshRenderSMShader->uniformLoc("factor"), _smFactor);
  //glUniform1f(_meshRenderSMShader->uniformLoc("units"), _smUnits);
  //glUseProgram(0);

  GLenum fboStatus;
	glShadeModel(GL_FLAT);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _shadowMapFboID);
  glViewport(0, 0, _shadowMapDimensions, _shadowMapDimensions);
  glPushMatrix();

  glClearColor(1.0,1.0,1.0,1.0);

  // Loop through all lights, and if shadowing is enabled, render the shadow map for it
  for (int l = 0; l < MAX_LIGHTS; l++)
  {    
    if (_lights[l] != NULL && _lights[l]->isActive() == true
     && _lights[l]->isShadowed() == true)
    {
      glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,
                                GL_COLOR_ATTACHMENT0_EXT,
                                GL_TEXTURE_2D,
                                _lights[l]->shadowMapID(), 0);
      fboStatus = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
      if (fboStatus != GL_FRAMEBUFFER_COMPLETE_EXT) {
        cout << " Framebuffer status is not complete; an error occurred." << endl;
      }
      
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      _lights[l]->viewTransformation();

      //drawGroundPlane(GROUND_DRAW_STYLE_SM);
       
      for (int x = 0; x < (int)_instanceList.size(); x++)
      {
        glPushMatrix();
#if 0
          glTranslatef(_instanceList[x]->position()[0],
                       _instanceList[x]->position()[1],
                       _instanceList[x]->position()[2]);
          glRotatef(_instanceList[x]->rotation()[0], 1.0f, 0.0f, 0.0f);
          glRotatef(_instanceList[x]->rotation()[1], 0.0f, 1.0f, 0.0f);
          glRotatef(_instanceList[x]->rotation()[2], 0.0f, 0.0f, 1.0f);
#endif
#if 0
          cout << "Calling instanceList[" << x << "]->renderSM()" << endl;
#endif
          _instanceList[x]->renderSM();
        glPopMatrix();
      }
    }
  }

  // Since shadow mapping is variance shadowmapping, blur the maps with separable gaussian blur
  blurShadowMaps();

  // Reset OpenGL state variables
  glClearColor(_clearColor[0],_clearColor[1],_clearColor[2],_clearColor[3]);

  glPopMatrix();
  glShadeModel(GL_SMOOTH);
  glViewport(_windowViewport[0], _windowViewport[1],
             _windowViewport[2], _windowViewport[3]);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
  
}

//////////////////////////////////////////////////////////////////////
// Steps the integrator systems
//////////////////////////////////////////////////////////////////////
void GLSL_SCENE_VIEWER::stepSystem( bool debug )
{
  // Update/deform all subspace temtmeshes, and generate the tet vertex position textures
  for (int x = 0; x < (int)_instanceList.size(); x++)
  {
    _instanceList[x]->integratorStep( debug );
  }
}

int MAX_NUM_LIGHTS = 3;

//////////////////////////////////////////////////////////////////////
// This will be where all mesh instances are drawn, on a background scene (ground plane)
//////////////////////////////////////////////////////////////////////
void GLSL_SCENE_VIEWER::display()
{
////////////////////////////////////////////////////////////////

/*  glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
  gluPerspective(20.0, 1.0, 0.1, 10.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0.0, 1.0, 2.5, 0.45, 0.45, 0.45, 0.0, 1.0, 0.0);
*/ // _lights[0]->viewTransformation();

////////////////////////////////////////////////////////////////
//
  
  // Render the scene with a deferred pass, and generate the ambient occlusion mask
  drawToDeferred();
 
  // FIXME - manually hack in a fixed number of lights here
  /*glLightfv(GL_LIGHT0, GL_POSITION, _lights[0]->position());
  glLightfv(GL_LIGHT0, GL_DIFFUSE, _lights[0]->diffuse());
  glLightfv(GL_LIGHT0, GL_AMBIENT, _lights[0]->ambient());
  glLightfv(GL_LIGHT0, GL_SPECULAR, _lights[0]->specular());

  glLightfv(GL_LIGHT1, GL_POSITION, _lights[1]->position());
  glLightfv(GL_LIGHT1, GL_DIFFUSE, _lights[1]->diffuse());
  glLightfv(GL_LIGHT1, GL_AMBIENT, _lights[1]->ambient());
  glLightfv(GL_LIGHT1, GL_SPECULAR, _lights[1]->specular());

  glLightfv(GL_LIGHT2, GL_POSITION, _lights[2]->position());
  glLightfv(GL_LIGHT2, GL_DIFFUSE, _lights[2]->diffuse());
  glLightfv(GL_LIGHT2, GL_AMBIENT, _lights[2]->ambient());
  glLightfv(GL_LIGHT2, GL_SPECULAR, _lights[2]->specular());

  glLightfv(GL_LIGHT3, GL_POSITION, _lights[3]->position());
  glLightfv(GL_LIGHT3, GL_DIFFUSE, _lights[3]->diffuse());
  glLightfv(GL_LIGHT3, GL_AMBIENT, _lights[3]->ambient());
  glLightfv(GL_LIGHT3, GL_SPECULAR, _lights[3]->specular());

  glLightfv(GL_LIGHT4, GL_POSITION, _lights[4]->position());
  glLightfv(GL_LIGHT4, GL_DIFFUSE, _lights[4]->diffuse());
  glLightfv(GL_LIGHT4, GL_AMBIENT, _lights[4]->ambient());
  glLightfv(GL_LIGHT4, GL_SPECULAR, _lights[4]->specular());*/
 
  // Render scene without shadows
  drawGroundPlane(GROUND_DRAW_STYLE_BASIC);
  //drawGroundPlane(GROUND_DRAW_STYLE_NO_SHADOW);

  //glLineWidth(2.0);
  
  for (int x = 0; x < (int)_instanceList.size(); x++)
  {
    glPushMatrix();
#if 0
      glTranslatef(_instanceList[x]->position()[0],
                   _instanceList[x]->position()[1],
                   _instanceList[x]->position()[2]);
      glRotatef(_instanceList[x]->rotation()[0], 1.0f, 0.0f, 0.0f);
      glRotatef(_instanceList[x]->rotation()[1], 0.0f, 1.0f, 0.0f);
      glRotatef(_instanceList[x]->rotation()[2], 0.0f, 0.0f, 1.0f);
#endif
#if 0
      cout << "Calling instanceList[" << x << "]->renderNoShadows()" << endl;
      cout << "instanceList[" << x << "]->position() = [";
      cout << _instanceList[x]->position()[0] << ", ";
      cout << _instanceList[x]->position()[1] << ", ";
      cout << _instanceList[x]->position()[2] << "]" << endl;
#endif
      //_instanceList[x]->renderNoShadows();
      _instanceList[x]->renderBasic();
      //_instanceList[x]->drawEmbeddedMeshFromTexture();

      //_instanceList[x]->renderFaces();
      //_instanceList[x]->renderBases();
     
    glPopMatrix();
  }// */

  for (int x = 0; x < (int)_objList.size(); x++)
  {
    glPushMatrix();
    //glTranslatef(0.0f,(GLfloat)(x-1),0.0f);
    _objList[x]->renderBasic();
    glPopMatrix();
  }
  // objList->renderBasic sets color, so reset to white
  glColor3f(1.0f,1.0f,1.0f);

  if (_useAO)
  {
    /*glMatrixMode(GL_PROJECTION);  glPushMatrix(); glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);   glPushMatrix(); glLoadIdentity();
    glDisable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, _aovAccumTexID);
   
    glUseProgram(_drawTextureShader->getHandle());
    glUniform1i(_drawTextureShader->uniformLoc("tex"), 0);

    glBegin(GL_QUADS);
      glTexCoord2f(0.0f, 0.0f);   glVertex2f(-1.0f, -1.0f);
      glTexCoord2f(1.0f, 0.0f);   glVertex2f( 1.0f, -1.0f);
      glTexCoord2f(1.0f, 1.0f);   glVertex2f( 1.0f,  1.0f);
      glTexCoord2f(0.0f, 1.0f);   glVertex2f(-1.0f,  1.0f);
    glEnd();

    glUseProgram(0);
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
    glMatrixMode(GL_PROJECTION);  glPopMatrix();
    glMatrixMode(GL_MODELVIEW);   glPopMatrix(); // */

    // Render AOV contribution by the ground plane (very fast) 
    glBlendEquation(GL_FUNC_REVERSE_SUBTRACT);
    glBlendFunc(GL_ONE, GL_ONE);
    glEnable(GL_BLEND);
    glDisable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);
    //glDisable(GL_CULL_FACE);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, _deferredTex0ID);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, _deferredTex1ID);

    drawGroundPlane(GROUND_DRAW_STYLE_AOV);

    // Blend the AOV buffer for the models over the current color buffer
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    //glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, _aovDeferredTex0ID);
    //glActiveTexture(GL_TEXTURE1);
    //glBindTexture(GL_TEXTURE_2D, _deferredTex0ID);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, _aovAccumTexID);

    glUseProgram(_bilateralShader->getHandle());
    glUniform1i(_bilateralShader->uniformLoc("defLowTex"), 0);
    glUniform1i(_bilateralShader->uniformLoc("defHiTex"), 1);
    glUniform1i(_bilateralShader->uniformLoc("aovTex"), 2);
    glUniform2f(_bilateralShader->uniformLoc("screenSize"),
                (GLfloat)_windowViewport[2],
                (GLfloat)_windowViewport[3] );
    //glEnable(GL_TEXTURE_2D);
            
    //glColor3f(1.0f,1.0f,1.0f);
    glBegin(GL_QUADS);
      glVertex3f(-1.0f,-1.0f, 0.0f);
      glVertex3f( 1.0f,-1.0f, 0.0f);
      glVertex3f( 1.0f, 1.0f, 0.0f);
      glVertex3f(-1.0f, 1.0f, 0.0f);
    glEnd();
  
    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    //glDisable(GL_TEXTURE_2D);
    glDepthMask(GL_TRUE);
    glBlendEquation(GL_FUNC_ADD);

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

    glEnable(GL_DEPTH_TEST);
    //glEnable(GL_CULL_FACE);
    //glDepthMask(GL_TRUE);
    glDisable(GL_BLEND);// */
  }
   
  glBindTexture(GL_TEXTURE_2D, 0);
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D, 0);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, 0);
  //glLineWidth(1.0);

  int lightCount = 0;
#if 0
  for (int l = 0; l < MAX_LIGHTS; l++)
  {
    // Render the scene once for every active light (this isn't good, should maybe change)
    if (_lights[l] != NULL && _lights[l]->isActive() == true)
    {
      lightCount++;

      //cout << "P: " << _lights[l]->position()[0] << ", " << _lights[l]->position()[1] << ", " << _lights[l]->position()[2] << endl;
      glLightfv(GL_LIGHT0, GL_POSITION, _lights[l]->position());
      glLightfv(GL_LIGHT0, GL_DIFFUSE, _lights[l]->diffuse());
      glLightfv(GL_LIGHT0, GL_AMBIENT, _lights[l]->ambient());
      glLightfv(GL_LIGHT0, GL_SPECULAR, _lights[l]->specular());

      if (_lights[l]->isShadowed() == true)
      {
        
//// Debug /////////////////////////////////////////////////////
  /*glPushMatrix();
  glEnable(GL_TEXTURE_2D);
  glDisable(GL_CULL_FACE);

  glBindTexture(GL_TEXTURE_2D, _lights[l]->shadowMapID());
  
  glBegin(GL_QUADS);
    glTexCoord2f(0.0f,0.0f);    glVertex3f( 0.0f,  0.0f, 0.1f);
    glTexCoord2f(1.0f,0.0f);    glVertex3f( 1.0f,  0.0f, 0.1f);
    glTexCoord2f(1.0f,1.0f);    glVertex3f( 1.0f,  1.0f, 0.1f);
    glTexCoord2f(0.0f,1.0f);    glVertex3f( 0.0f,  1.0f, 0.1f);
  glEnd();

  glBindTexture(GL_TEXTURE_2D, 0);
  glDisable(GL_TEXTURE_2D);
  glEnable(GL_CULL_FACE);
  glPopMatrix();  // */
////////////////////////////////////////////////////////////////

        // Render with shadow mapping. .renderWithShadows() expects shadowmap in texture 4
        glActiveTexture(GL_TEXTURE4);
        glBindTexture(GL_TEXTURE_2D, _lights[l]->shadowMapID());
        _lights[l]->textureTransformation();

        //drawGroundPlane(GROUND_DRAW_STYLE_SHADOW);

        glUseProgram(_meshRenderWithShadowShader->getHandle());
        glUniform1f(_meshRenderWithShadowShader->uniformLoc("minVariance"), _minVariance);
        glUniform1f(_meshRenderWithShadowShader->uniformLoc("reduceBleed"), _reduceBleed);

        // Render all instances with shadows
        for (int x = 0; x < (int)_instanceList.size(); x++)
        {
#if 0
          glPushMatrix();
            glTranslatef(_instanceList[x]->position()[0],
                         _instanceList[x]->position()[1],
                         _instanceList[x]->position()[2]);
            glRotatef(_instanceList[x]->rotation()[0], 1.0f, 0.0f, 0.0f);
            glRotatef(_instanceList[x]->rotation()[1], 0.0f, 1.0f, 0.0f);
            glRotatef(_instanceList[x]->rotation()[2], 0.0f, 0.0f, 1.0f);
#endif
#if 0
            cout << "Calling instanceList[" << x << "]->renderWithShadows()" << endl;
#endif
            _instanceList[x]->renderWithShadows();
            //_instanceList[x]->renderNoShadows();
            
          glPopMatrix();
        }
        glActiveTexture(GL_TEXTURE4);
        glBindTexture(GL_TEXTURE_2D, 0);
        glMatrixMode(GL_TEXTURE);
        glLoadIdentity();
        glMatrixMode(GL_MODELVIEW);
      }
      else // _lights[l]->isShadowed() == false
      {
        // Render scene without shadows
        //drawGroundPlane(GROUND_DRAW_STYLE_NO_SHADOW);

        //glLineWidth(2.0);
        
        for (int x = 0; x < (int)_instanceList.size(); x++)
        {
          glPushMatrix();
#if 0
            glTranslatef(_instanceList[x]->position()[0],
                         _instanceList[x]->position()[1],
                         _instanceList[x]->position()[2]);
            glRotatef(_instanceList[x]->rotation()[0], 1.0f, 0.0f, 0.0f);
            glRotatef(_instanceList[x]->rotation()[1], 0.0f, 1.0f, 0.0f);
            glRotatef(_instanceList[x]->rotation()[2], 0.0f, 0.0f, 1.0f);
#endif
#if 0
            cout << "Calling instanceList[" << x << "]->renderNoShadows()" << endl;
            cout << "instanceList[" << x << "]->position() = [";
            cout << _instanceList[x]->position()[0] << ", ";
            cout << _instanceList[x]->position()[1] << ", ";
            cout << _instanceList[x]->position()[2] << "]" << endl;
#endif
            _instanceList[x]->renderNoShadows();
            //_instanceList[x]->drawEmbeddedMeshFromTexture();

            //_instanceList[x]->renderFaces();
            //_instanceList[x]->renderBases();
           
          glPopMatrix();
        }// */

        if (_useAO)
        {
          // Render AOV contribution by the ground plane (very fast) 
          glBlendEquation(GL_FUNC_REVERSE_SUBTRACT);
          glBlendFunc(GL_ONE, GL_ONE);
          glEnable(GL_BLEND);
		      glDisable(GL_DEPTH_TEST);
		      glDepthMask(GL_FALSE);
          //glDisable(GL_CULL_FACE);

          glActiveTexture(GL_TEXTURE1);
          glBindTexture(GL_TEXTURE_2D, _deferredTex0ID);
          glActiveTexture(GL_TEXTURE2);
          glBindTexture(GL_TEXTURE_2D, _deferredTex1ID);

          //drawGroundPlane(GROUND_DRAW_STYLE_AOV);

          // Blend the AOV buffer for the models over the current color buffer
          glMatrixMode(GL_PROJECTION);
		      glPushMatrix();
		      glLoadIdentity();
		      glMatrixMode(GL_MODELVIEW);
		      glPushMatrix();
		      glLoadIdentity();
		      //glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
        
          glActiveTexture(GL_TEXTURE0);
		      glBindTexture(GL_TEXTURE_2D, _aovDeferredTex0ID);
          //glActiveTexture(GL_TEXTURE1);
		      //glBindTexture(GL_TEXTURE_2D, _deferredTex0ID);
          glActiveTexture(GL_TEXTURE2);
		      glBindTexture(GL_TEXTURE_2D, _aovAccumTexID);

          glUseProgram(_bilateralShader->getHandle());
          glUniform1i(_bilateralShader->uniformLoc("defLowTex"), 0);
          glUniform1i(_bilateralShader->uniformLoc("defHiTex"), 1);
          glUniform1i(_bilateralShader->uniformLoc("aovTex"), 2);
          glUniform2f(_bilateralShader->uniformLoc("screenSize"),
                      (GLfloat)_windowViewport[2],
                      (GLfloat)_windowViewport[3] );
          //glEnable(GL_TEXTURE_2D);
		    		   		
		      //glColor3f(1.0f,1.0f,1.0f);
		      glBegin(GL_QUADS);
			      glVertex3f(-1.0f,-1.0f, 0.0f);
			      glVertex3f( 1.0f,-1.0f, 0.0f);
			      glVertex3f( 1.0f, 1.0f, 0.0f);
			      glVertex3f(-1.0f, 1.0f, 0.0f);
		      glEnd();
        
		      glDisable(GL_BLEND);
		      glEnable(GL_DEPTH_TEST);
		      //glDisable(GL_TEXTURE_2D);
		      glDepthMask(GL_TRUE);
          glBlendEquation(GL_FUNC_ADD);

		      glMatrixMode(GL_PROJECTION);
		      glPopMatrix();
		      glMatrixMode(GL_MODELVIEW);
		      glPopMatrix();

          glEnable(GL_DEPTH_TEST);
          //glEnable(GL_CULL_FACE);
          //glDepthMask(GL_TRUE);
          glDisable(GL_BLEND);// */
        }
         
        glBindTexture(GL_TEXTURE_2D, 0);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, 0);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, 0);
        //glLineWidth(1.0);
      }
    }
  }
#endif
  
  reportGLError(glGetError());

  //drawToDeferred();
  //drawSSAO();

  if (lightCount == 0)
  {    
    /*for (int x = 0; x < (int)_instanceList.size(); x++)
    {
      //glPushMatrix();
        //glTranslatef(_instanceList[x]->position()[0], _instanceList[x]->position()[1], _instanceList[x]->position()[2]);
        //glRotatef(_instanceList[x]->rotation()[0], 1.0f, 0.0f, 0.0f);
        //glRotatef(_instanceList[x]->rotation()[1], 0.0f, 1.0f, 0.0f);
        //glRotatef(_instanceList[x]->rotation()[2], 0.0f, 0.0f, 1.0f);
        //_instanceList[x]->renderNoShadows();
      //glPushMatrix();
    }*/
  }
}

//////////////////////////////////////////////////////////////////////
// Draw the scene to a deferred render target (rgb = normal, a = depth)
//////////////////////////////////////////////////////////////////////
void GLSL_SCENE_VIEWER::drawToDeferred()
{
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _deferredFboID);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,
                            GL_COLOR_ATTACHMENT0_EXT,
                            GL_TEXTURE_2D,
                            _deferredTex0ID, 0);
  GLenum fboStatus = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
  reportFboStatus(fboStatus);
  glClearColor(0.0,0.0,0.0,0.0);
 
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

  drawGroundPlane(GROUND_DRAW_STYLE_DEFERRED);

  for (int x = 0; x < (int)_instanceList.size(); x++)
  {
    glPushMatrix();
#if 0
      glTranslatef(_instanceList[x]->position()[0],
                   _instanceList[x]->position()[1],
                   _instanceList[x]->position()[2]);
      glRotatef(_instanceList[x]->rotation()[0], 1.0f, 0.0f, 0.0f);
      glRotatef(_instanceList[x]->rotation()[1], 0.0f, 1.0f, 0.0f);
      glRotatef(_instanceList[x]->rotation()[2], 0.0f, 0.0f, 1.0f);
#endif
#if 0
      cout << "Calling instanceList[" << x << "]->renderDeferred()" << endl;
#endif
      _instanceList[x]->renderDeferred();
    glPopMatrix();
  }

  for (int x = 0; x < (int)_objList.size(); x++)
  {
    glPushMatrix();
    //glTranslatef(0.0f,(GLfloat)(x-1),0.0f);
    _objList[x]->renderDeferred();
    glPopMatrix();
  }
 
  if (_useAO)
  {
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _aovDeferredFboID);
    //fboStatus = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
    //reportFboStatus(fboStatus);

    // Render lower resolution deferred buffers for AOV
    glViewport(0.0,0.0,_aovAccumWidth,_aovAccumHeight);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    drawGroundPlane(GROUND_DRAW_STYLE_DEFERRED);

    for (int x = 0; x < (int)_instanceList.size(); x++)
    {
      glPushMatrix();
#if 0
        glTranslatef(_instanceList[x]->position()[0],
                     _instanceList[x]->position()[1],
                     _instanceList[x]->position()[2]);
        glRotatef(_instanceList[x]->rotation()[0], 1.0f, 0.0f, 0.0f);
        glRotatef(_instanceList[x]->rotation()[1], 0.0f, 1.0f, 0.0f);
        glRotatef(_instanceList[x]->rotation()[2], 0.0f, 0.0f, 1.0f);
#endif
#if 0
        cout << "Calling instanceList[" << x << "]->renderDeferred()" << endl;
#endif
        _instanceList[x]->renderDeferred();
      glPopMatrix();
    }

    for (int x = 0; x < (int)_objList.size(); x++)
    {
      glPushMatrix();
      //glTranslatef(0.0f,(GLfloat)(x-1),0.0f);
      _objList[x]->renderDeferred();
      glPopMatrix();
    }
 
    // Using low res buffers, render AOV (accum and deferred share a depth buffer)
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _aovAccumFboID);
    //fboStatus = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
    //reportFboStatus(fboStatus);
    //glClearColor(0.0,0.0,0.0,0.0);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, _aovDeferredTex0ID);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, _aovDeferredTex1ID);

    glClear(GL_COLOR_BUFFER_BIT);

    glBlendEquation(GL_FUNC_ADD);
    //glBlendEquation(GL_FUNC_REVERSE_SUBTRACT);
    glBlendFunc(GL_ONE, GL_ONE);
    glEnable(GL_BLEND);
    glDepthMask(GL_FALSE);
    //glDisable(GL_CULL_FACE);
    //glDisable(GL_DEPTH_TEST);

    //_aovDelta = 0.05f;
       
    drawGroundPlane(GROUND_DRAW_STYLE_AOV);
       
    for (int x = 0; x < (int)_instanceList.size(); x++)
    {
      glPushMatrix();
#if 0
        glTranslatef(_instanceList[x]->position()[0],
                     _instanceList[x]->position()[1],
                     _instanceList[x]->position()[2]);
        glRotatef(_instanceList[x]->rotation()[0], 1.0f, 0.0f, 0.0f);
        glRotatef(_instanceList[x]->rotation()[1], 0.0f, 1.0f, 0.0f);
        glRotatef(_instanceList[x]->rotation()[2], 0.0f, 0.0f, 1.0f);
#endif
        //_instanceList[x]->renderFaces();
        //_instanceList[x]->renderBases();
#if 0
        cout << "Calling instanceList[" << x << "]->renderAOVEmbedded()" << endl;
#endif
        _instanceList[x]->renderAOVEmbedded(_aovDelta, _aovAccumWidth,
                                            _aovAccumHeight);
        //_instanceList[x]->renderAOVWireframe(_aovDelta);

      glPopMatrix();
    }// */

    for (int x = 0; x < (int)_objList.size(); x++)
    {
      glPushMatrix();
      //glTranslatef(0.0f,(GLfloat)(x-1),0.0f);
      _objList[x]->renderAOV(_aovDelta, _aovAccumWidth, _aovAccumHeight);
      glPopMatrix();
    }

    glDisable(GL_BLEND);
    glDepthMask(GL_TRUE);
  }

  glClearColor(_clearColor[0],_clearColor[1],_clearColor[2],_clearColor[3]);
  glViewport(_windowViewport[0], _windowViewport[1],
             _windowViewport[2], _windowViewport[3]);

  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
}

//////////////////////////////////////////////////////////////////////
// Using the deferred render target, generate a screen-space ambient
// occlusion modulation texture
//////////////////////////////////////////////////////////////////////
void GLSL_SCENE_VIEWER::drawSSAO()
{
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glDisable(GL_DEPTH_TEST);

  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, _deferredFboID);
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT,
                            GL_TEXTURE_2D, _ssaoTexID, 0);
  
  glViewport(0, 0, _windowViewport[2], _windowViewport[3]);
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

  glUseProgram(_ssoaRenderShader->getHandle());
  glUniform1i(_ssoaRenderShader->uniformLoc("texture"), 0);
				
  glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, _deferredTex0ID);
	glBegin(GL_QUADS);
    glTexCoord2f(0.0f, 0.0f);        glVertex2f(-1.0f, -1.0f);
    glTexCoord2f(1.0f, 0.0f);        glVertex2f( 1.0f, -1.0f);
    glTexCoord2f(1.0f, 1.0f);        glVertex2f( 1.0f,  1.0f);
    glTexCoord2f(0.0f, 1.0f);        glVertex2f(-1.0f,  1.0f);
  glEnd();

  blurTex(_ssaoTexID, _ssBlurTexID, _windowViewport[3]);
  
  glUseProgram(0);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);

  glViewport(0, 0, _windowViewport[2], _windowViewport[3]);

  glBlendFunc(GL_ZERO, GL_SRC_ALPHA);
  glEnable(GL_TEXTURE_2D);
  glEnable(GL_BLEND);
  glDepthMask(GL_FALSE);
  
  glBindTexture(GL_TEXTURE_2D, _ssaoTexID);
	glBegin(GL_QUADS);
    glTexCoord2f(0.0f, 0.0f);        glVertex2f(-1.0f, -1.0f);
    glTexCoord2f(1.0f, 0.0f);        glVertex2f( 1.0f, -1.0f);
    glTexCoord2f(1.0f, 1.0f);        glVertex2f( 1.0f,  1.0f);
    glTexCoord2f(0.0f, 1.0f);        glVertex2f(-1.0f,  1.0f);
  glEnd();

  glDepthMask(GL_TRUE);
  glDisable(GL_BLEND);
  glDisable(GL_TEXTURE_2D);

  glBindTexture(GL_TEXTURE_2D, 0);
  
  glEnable(GL_DEPTH_TEST);
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}

//////////////////////////////////////////////////////////////////////
// Given a texture 'src', and another temp texture 'blur', the same 
// square dimensions 'dim', blur the texture
//////////////////////////////////////////////////////////////////////
void GLSL_SCENE_VIEWER::blurTex(GLuint src, GLuint blur, GLsizei dim)
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

	// Change to the 2-pass separable gaussian blurring shader for shadowmaps,
  // and pass it some variables
  glUseProgram(_shadowBlurShader->getHandle());

  glUniform1i(_shadowBlurShader->uniformLoc("texture"), 0);
  glUniform1f(_shadowBlurShader->uniformLoc("texDim"), (float)dim);

  glDisable(GL_DEPTH_TEST);
  glActiveTexture(GL_TEXTURE0);

  // Set the shadow blur texture as the render target
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT,
                            GL_TEXTURE_2D, blur, 0);
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
				
  // Bind the recently created shadowmap, and draw it.
  glBindTexture(GL_TEXTURE_2D, src);
  glBegin(GL_QUADS);
    // Note: this gets drawn as if rotated 90 degrees to the right
    glTexCoord2f(1.0f, 0.0f);        glVertex2f(-1.0f, -1.0f);
    glTexCoord2f(1.0f, 1.0f);        glVertex2f( 1.0f, -1.0f);
    glTexCoord2f(0.0f, 1.0f);        glVertex2f( 1.0f,  1.0f);
    glTexCoord2f(0.0f, 0.0f);        glVertex2f(-1.0f,  1.0f);
  glEnd();

  // Switch render target and bound texture
  glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT,
                            GL_TEXTURE_2D, src, 0);
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
			
  glBindTexture(GL_TEXTURE_2D, blur);
  glBegin(GL_QUADS);
    // Note: this gets drawn as if rotated 90 degrees to the left
    glTexCoord2f(0.0f, 1.0f);        glVertex2f(-1.0f, -1.0f);
    glTexCoord2f(0.0f, 0.0f);        glVertex2f( 1.0f, -1.0f);
    glTexCoord2f(1.0f, 0.0f);        glVertex2f( 1.0f,  1.0f);
    glTexCoord2f(1.0f, 1.0f);        glVertex2f(-1.0f,  1.0f);
  glEnd();

  // Bind the shadowmap, and generate some mipmaps (used in variance
  // and exponential shadow mapping)
  glBindTexture(GL_TEXTURE_2D, src);
  glGenerateMipmapEXT(GL_TEXTURE_2D);	

  glBindTexture(GL_TEXTURE_2D, 0);

  glEnable(GL_DEPTH_TEST);
  glUseProgram(0);

}

//////////////////////////////////////////////////////////////////////
// Draw simple ground plane
//////////////////////////////////////////////////////////////////////
void GLSL_SCENE_VIEWER::drawGroundPlane(int style)
{
  if (_drawScene == false)
    return;

  if (style == GROUND_DRAW_STYLE_SHADOW)
  {
    glUseProgram(_sceneRenderWithShadowShader->getHandle());
    glUniform1i(_sceneRenderWithShadowShader->uniformLoc("shadowTex"), 4);
    glUniform1f(_sceneRenderWithShadowShader->uniformLoc("minVariance"),
                _minVariance);
    glUniform1f(_sceneRenderWithShadowShader->uniformLoc("reduceBleed"),
                _reduceBleed);

  }
  else if (style == GROUND_DRAW_STYLE_NO_SHADOW)
  {
    glUseProgram(_sceneRenderNoShadowShader->getHandle());
  }
  else if (style == GROUND_DRAW_STYLE_SM)
  {
    glUseProgram(_sceneRenderSMShader->getHandle());
    //glUniform1f(_sceneRenderSMShader->uniformLoc("factor"), _smFactor);
    //glUniform1f(_sceneRenderSMShader->uniformLoc("units"), _smUnits);
  }
  else if (style == GROUND_DRAW_STYLE_DEFERRED)
  {
    glUseProgram(_sceneRenderDeferredShader->getHandle());
  }
  else if (style == GROUND_DRAW_STYLE_BASIC)
  {
    glUseProgram(_sceneRenderBasicShader->getHandle());
  }
  else if (style == GROUND_DRAW_STYLE_AOV)
  {
    glUseProgram(_sceneRenderAOVShader->getHandle());
      
    glUniform1i(_sceneRenderAOVShader->uniformLoc("normalTex"), 1);
    glUniform1i(_sceneRenderAOVShader->uniformLoc("positionTex"), 2);
    //glUniform2f(_sceneRenderAOVShader->uniformLoc("screenSize"),
                  //(GLfloat)_aovAccumWidth, (GLfloat)_aovAccumHeight);
    glUniform2f(_sceneRenderAOVShader->uniformLoc("screenSize"),
                (GLfloat)_windowViewport[2], (GLfloat)_windowViewport[3]);

    //glUniform1f(_sceneRenderAOVShader->uniformLoc("delta"), _aovDelta);

    /*glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glTranslatef(0.0f, 0.0f,-0.999f);// */

  }

  glBindBuffer(GL_ARRAY_BUFFER, _groundPlaneVboID);
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(3, GL_FLOAT, 9*sizeof(GLfloat), (void*)0);

  if (style == GROUND_DRAW_STYLE_SHADOW || style == GROUND_DRAW_STYLE_NO_SHADOW)
  {
    glBindBuffer(GL_ARRAY_BUFFER, _groundPlaneVboID);
    glEnableClientState(GL_NORMAL_ARRAY);
    glNormalPointer(GL_FLOAT, 9*sizeof(GLfloat), (void*)(3*sizeof(GLfloat)));

    glBindBuffer(GL_ARRAY_BUFFER, _groundPlaneVboID);
    glEnableClientState(GL_COLOR_ARRAY);
    glColorPointer(3, GL_FLOAT, 9*sizeof(GLfloat), (void*)(6*sizeof(GLfloat)));
  }

  glPushMatrix();
  glTranslatef(0.0f,0.0f,_groundHeight);

  glRotatef(-90.0, 1,0,0);

  glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

  glPopMatrix();

  if (style == GROUND_DRAW_STYLE_SHADOW || style == GROUND_DRAW_STYLE_NO_SHADOW)
  {
    glDisableClientState(GL_COLOR_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
  }
  else if (style == GROUND_DRAW_STYLE_AOV)
  {
    /*glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();// */
  }


  glDisableClientState(GL_VERTEX_ARRAY);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
}


//////////////////////////////////////////////////////////////////////
// Conduct a 5x5 (or 7x7) separable gaussian blur on all active shadow maps
// Assumes that _shadowMapFboID is currently bound FBO
//////////////////////////////////////////////////////////////////////
void GLSL_SCENE_VIEWER::blurShadowMaps()
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

	// Change to the 2-pass separable gaussian blurring shader for shadowmaps,
  // and pass it some variables
  glUseProgram(_shadowBlurShader->getHandle());

  glUniform1i(_shadowBlurShader->uniformLoc("texture"), 0);
  glUniform1f(_shadowBlurShader->uniformLoc("texDim"),
              (float)_shadowMapDimensions);

  glDisable(GL_DEPTH_TEST);
  glActiveTexture(GL_TEXTURE0);

	for (int x = 0; x < MAX_LIGHTS; x++)
  {
    if (_lights[x] != NULL && _lights[x]->isActive() == true
     && _lights[x]->isShadowed() == true)
    {
	    // Set the shadow blur texture as the render target
    	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,
                                GL_COLOR_ATTACHMENT0_EXT,
                                GL_TEXTURE_2D,
                                _shadowBlurTempTexID, 0);
	    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
				
	    // Bind the recently created shadowmap, and draw it.
	    glBindTexture(GL_TEXTURE_2D, _lights[x]->shadowMapID());
	    glBegin(GL_QUADS);
        // Note: this gets drawn as if rotated 90 degrees to the right
        glTexCoord2f(1.0f, 0.0f);        glVertex2f(-1.0f, -1.0f);
        glTexCoord2f(1.0f, 1.0f);        glVertex2f( 1.0f, -1.0f);
        glTexCoord2f(0.0f, 1.0f);        glVertex2f( 1.0f,  1.0f);
        glTexCoord2f(0.0f, 0.0f);        glVertex2f(-1.0f,  1.0f);
      glEnd();

      // Switch render target and bound texture
      glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,
                                GL_COLOR_ATTACHMENT0_EXT,
                                GL_TEXTURE_2D, _lights[x]->shadowMapID(), 0);
      glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
			
      glBindTexture(GL_TEXTURE_2D, _shadowBlurTempTexID);
      glBegin(GL_QUADS);
        // Note: this gets drawn as if rotated 90 degrees to the left
        glTexCoord2f(0.0f, 1.0f);        glVertex2f(-1.0f, -1.0f);
        glTexCoord2f(0.0f, 0.0f);        glVertex2f( 1.0f, -1.0f);
        glTexCoord2f(1.0f, 0.0f);        glVertex2f( 1.0f,  1.0f);
        glTexCoord2f(1.0f, 1.0f);        glVertex2f(-1.0f,  1.0f);
      glEnd();

      // Bind the shadowmap, and generate some mipmaps (used in variance
      // and exponential shadow mapping)
      glBindTexture(GL_TEXTURE_2D, _lights[x]->shadowMapID());
      glGenerateMipmapEXT(GL_TEXTURE_2D);	
    }
  }

  glBindTexture(GL_TEXTURE_2D, 0);

  glEnable(GL_DEPTH_TEST);
  glUseProgram(0);
}


//////////////////////////////////////////////////////////////////////
// Add an OBJ pointer to list of OBJs.
//////////////////////////////////////////////////////////////////////
void GLSL_SCENE_VIEWER::addObj(OBJ *obj)
{
  if (obj != NULL)
  {
    GLSL_OBJ *glslObj = new GLSL_OBJ(obj);
    _objList.push_back(glslObj);
  }
}

//////////////////////////////////////////////////////////////////////
// Add an OBJ pointer to list of OBJs.
//////////////////////////////////////////////////////////////////////
void GLSL_SCENE_VIEWER::addObj(OBJ *obj, GLfloat r, GLfloat g, GLfloat b)
{
  if (obj != NULL)
  {
    GLSL_OBJ *glslObj = new GLSL_OBJ(obj);
    glslObj->setColor(r,g,b);
    _objList.push_back(glslObj);
  }
}

//////////////////////////////////////////////////////////////////////
// Add a mesh instance given a config file
//////////////////////////////////////////////////////////////////////
int GLSL_SCENE_VIEWER::addMesh(string cfgFilename)
{
  GLSL_TET_MESH* glslTetMesh = NULL;

  // Check to see if this cfg file has been opened before
  if (_cfgIndexMap.find(cfgFilename) == _cfgIndexMap.end())
  {
    // If not, then load up the SUBSPACE_TET_MESH and OBJ given by the cfg:
    SUBSPACE_TET_MESH *tetMesh;
    OBJ *embeddedMesh;
    if (loadTetMesh(cfgFilename, &tetMesh, &embeddedMesh) == true)
    {
      glslTetMesh = new GLSL_TET_MESH(tetMesh, embeddedMesh);
      if (glslTetMesh == NULL)
      {
        cout << "Add Mesh failed to create a new GLSL_TET_MESH using";
        cout << " config file \'" << cfgFilename << "\'" << endl;
        return -1;
      }
    }
    else 
    {
      cout << "Add Mesh failed to open the file \'" << cfgFilename;
      cout << "\'" << endl;
      return -1;
    }
    
    // put the new GLSL_TET_MESH object into the array
    _meshList.push_back(glslTetMesh);
    int oldmeshsize = _meshList.size() - 1;
    
    // Store the _meshList index of the GLSL_TET_MESH in the cfg index map
    _cfgIndexMap[cfgFilename] = oldmeshsize;
  }
  else
  {
    glslTetMesh = _meshList[_cfgIndexMap[cfgFilename]];
    if (glslTetMesh == NULL)
    {
      cout << "Stored GLSL_TET_MESH from config file \'" << cfgFilename;
      cout << "\' was NULL" << endl;
      return -1;
    }
  }

  GLSL_TET_MESH_INSTANCE *meshInstance = new GLSL_TET_MESH_INSTANCE(glslTetMesh,
                                                                    cfgFilename);

  meshInstance->setIntegratorOptions( INTEGRATOR_STEP_INVERTIBLE_IMPLICIT );

  if (meshInstance == NULL)
  {
    cout << "Add Mesh failed to create a new GLSL_TET_MESH_INSTANCE" << endl;
    return -1;
  }

  meshInstance->setMeshRenderToTextureShader(_embeddedMeshToTextureShader);
  meshInstance->setMeshRenderWithShadowShader(_meshRenderWithShadowShader);
  meshInstance->setMeshRenderNoShadowShader(_meshRenderNoShadowShader);
  meshInstance->setMeshRenderSMShader(_meshRenderSMShader);
  meshInstance->setMeshRenderDeferredShader(_meshRenderDeferredShader);
  meshInstance->setMeshRenderBasicShader(_meshRenderBasicShader);
  meshInstance->setMeshRenderNormalsShader(_meshRenderNormalsShader);

  meshInstance->setNormTanTexShader(_tetNormTanTex0Shader, 0);
  meshInstance->setNormTanTexShader(_tetNormTanTex1Shader, 1);
  meshInstance->setNormTanTexShader(_tetNormTanTex2Shader, 2);

  glClearColor(0.0,0.0,0.0,0.0);
  meshInstance->genTetVertexTexture();
  
  _instanceList.push_back(meshInstance);
  int oldinstsize = _instanceList.size() - 1;
  
  // Return the meshid: the index of the mesh instance in _instanceList
  return oldinstsize;
}

//////////////////////////////////////////////////////////////////////////////
// Return a mesh instance pointer given a meshid
//////////////////////////////////////////////////////////////////////////////
GLSL_TET_MESH_INSTANCE* GLSL_SCENE_VIEWER::getMesh(int meshid)
{
  if (meshid < (int)_instanceList.size())
  {
    return _instanceList[meshid];
  }

  cout << "Get Mesh was passed a bad meshid" << endl;
  return NULL;
}

//////////////////////////////////////////////////////////////////////////////
// Return a mesh instance pointer given a meshid
//////////////////////////////////////////////////////////////////////////////
GLSL_LIGHT* GLSL_SCENE_VIEWER::getLight(int index)
{
  // Check index parameter
  if (index < 0 || index >= MAX_LIGHTS)
  {
    cout << " getLight passed an invalid index parameter" << endl;
    return NULL;
  }
  // Lazy allocation, so only allocate new light at this index if it is requested
  if (_lights[index] == NULL)
  {
    _lights[index] = new GLSL_LIGHT();
    _lights[index]->setShadowDimensions(_shadowMapDimensions,
                                        _shadowMapDimensions);
  }

  // Return light object pointer
  return _lights[index];
}

//////////////////////////////////////////////////////////////////////
// Sets the height of the ground plane
//////////////////////////////////////////////////////////////////////
void GLSL_SCENE_VIEWER::setGroundHeight(GLfloat height)
{
  _groundHeight = height;
}

//////////////////////////////////////////////////////////////////////
// Set the scene wide delta value for AOV
//////////////////////////////////////////////////////////////////////
void GLSL_SCENE_VIEWER::setAOVDelta(GLfloat delta)
{
  _aovDelta = (delta < 0.0) ? 0.0f : delta;
}

//////////////////////////////////////////////////////////////////////
// Given a point in eye space, determine affected mesh, and pass point along
//////////////////////////////////////////////////////////////////////
bool GLSL_SCENE_VIEWER::mouseClick(VEC3F& point, Real maxRadius)
{
  // determine closest instance mesh
  // pass point to that instance's integrator
  // store that instance in _animatedMesh
  int closest = -1;
  float closestDist = 10000.0f;
  VEC3F displace;
  float displaceDist;

  for (int x = 0; x < (int)_instanceList.size(); x++)
  {
    displace[0] = _instanceList[x]->position()[0] - point[0];
    displace[1] = _instanceList[x]->position()[1] - point[1];
    displace[2] = _instanceList[x]->position()[2] - point[2];
    displaceDist = displace[0]*displace[0] + displace[1]*displace[1]
                 + displace[2]*displace[2];
    if (displaceDist < closestDist)
    {
      closest = x;
      closestDist = displaceDist;
    }
  }

  if (closest != -1)
  {
    _animatedMesh = _instanceList[closest];
    return _animatedMesh->integratorClick(point, maxRadius);
  }
  else
  {
    return false;
  }
}

//////////////////////////////////////////////////////////////////////
// Pass point to currently affected mesh
//////////////////////////////////////////////////////////////////////
void GLSL_SCENE_VIEWER::mouseDrag(VEC3F& point)
{
  // pass point as a drag point to integrator of _animatedMesh
  if (_animatedMesh != NULL)
  {
    _animatedMesh->integratorDrag(point);
  }
}

//////////////////////////////////////////////////////////////////////
// Tell affected mesh that the mouse is unclicked
//////////////////////////////////////////////////////////////////////
void GLSL_SCENE_VIEWER::mouseUnclick()
{
  // tell _animatedMesh's integrator to unclick
  // set _animatedMesh to NULL
  if (_animatedMesh != NULL)
  {
    _animatedMesh->integratorUnclick();
    _animatedMesh = NULL;
  }
}

//////////////////////////////////////////////////////////////////////////////
// Read a material file and return a pointer to it
//////////////////////////////////////////////////////////////////////////////
MATERIAL* GLSL_SCENE_VIEWER::readMaterial(SIMPLE_PARSER& parser)
{
  return SIMPLE_PARSER::READ_MATERIAL( parser );
}

//////////////////////////////////////////////////////////////////////////////
// Read a config file, load in the tet mesh and embedded mesh specified.
//////////////////////////////////////////////////////////////////////////////
bool GLSL_SCENE_VIEWER::loadTetMesh(string cfgFilename,
                                    SUBSPACE_TET_MESH **mesh,
                                    OBJ **obj)
{
  // output pointers
  SUBSPACE_TET_MESH* tetMesh;
  OBJ* originalMesh;

  // input parameters
  int bccRes = 32;
  string triangleMeshPath;
  string triangleMeshName;
  string outputPath;
  string tetMeshName;

  // read in different parameters
  string configName(cfgFilename);
  SIMPLE_PARSER configFile(configName); 
  bccRes           = configFile.getInt("bccRes", bccRes);
  triangleMeshPath = configFile.getString("triangle mesh path", triangleMeshPath);
  triangleMeshName = configFile.getString("triangle mesh name", triangleMeshName);
  outputPath       = configFile.getString("output path", outputPath);
  string originalMeshFile = triangleMeshPath + triangleMeshName;

  if (configFile.defined("tet mesh name"))
    tetMeshName = configFile.getString("tet mesh name", tetMeshName);
  else
    tetMeshName = outputPath + triangleMeshName + string(".tetmesh");

  // read in how many materials there are
  int totalMaterials = 0;
  totalMaterials = configFile.getInt("total materials", totalMaterials);
  if (totalMaterials == 0)
  {
    cout << " NO MATERIALS SPECIFIED!!!!" << endl;
    return false;
  }
  
  // read in the actual materials
  MATERIAL** materials = new MATERIAL*[totalMaterials];
  for (int x = 0; x < totalMaterials; x++)
  {
    // read in the config file name for the material
    char buffer[256];
    sprintf(buffer, "material %i", x);
    string materialString(buffer);
    string materialFilename;
    materialFilename = configFile.getString(materialString.c_str(),
                                            materialFilename);

    // open the config file
    SIMPLE_PARSER materialFile(materialFilename);
    
    // get the material
    MATERIAL* material = readMaterial(materialFile);
    
    // set the invertible wrapper if necessary
    bool invertible = false;
    invertible = configFile.getBool("invertible", invertible);
    if (invertible)
    {
      cout << " Setting material to invertible" << endl;
      material = new INVERTIBLE(material);
    }

    materials[x] = material;
  }

  Real meshScale           = configFile.getFloat("mesh scale", 1.0);
  
  cout << "Got mesh scale " << meshScale << endl;

  // output precision
  cout << " Precision: " << sizeof(Real) * 8 << " bit" << endl;
  
  cout << "==================================================" << endl;
  cout << " Reading mesh file " << tetMeshName.c_str() << endl;
  cout << "==================================================" << endl;
  tetMesh = new SUBSPACE_TET_MESH(tetMeshName.c_str(), materials, totalMaterials,
                                  false, NULL, NULL, meshScale);

  // load the head file, since we need the normalization
  originalMesh = new OBJ();
  originalMesh->Load(originalMeshFile.c_str(), meshScale);
  originalMesh->normalize(bccRes);

  tetMesh->normalizeEmbedding(originalMesh);
  tetMesh->embedMesh(originalMesh);

  // Write out the mesh and obj
  (*mesh) = tetMesh;
  (*obj) = originalMesh;

  Real tetBBox[6];
  Real objBBox[6];

  tetMesh->boundingBox( tetBBox );
  originalMesh->BoundingBox( objBBox );

  return true;
}

// Prints framebuffer object status, and returns true if an error occurred, false otherwise
bool GLSL_SCENE_VIEWER::reportFboStatus(GLenum fboStatus)
{
  switch (fboStatus)
  {
    case GL_FRAMEBUFFER_COMPLETE_EXT:
      //cout << " Framebuffer status is complete." << endl;
      return false;
      break;
    case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT_EXT:
      cout << " Framebuffer status is incomplete attachment!" << endl;
      break;
    case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT_EXT:
      cout << " Framebuffer status is imcomplete missing attachment!" << endl;
      break;
    case GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS_EXT:
      cout << " Framebuffer status is incomplete dimensions!" << endl;
      break;
    case GL_FRAMEBUFFER_INCOMPLETE_FORMATS_EXT:
      cout << " Framebuffer status is incomplete formats!" << endl;
      break;
    case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER_EXT:
      cout << " Framebuffer status is incomplete draw buffer!" << endl;
      break;
    case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER_EXT:
      cout << " Framebuffer status is incomplete read buffer!" << endl;
      break;
    case GL_FRAMEBUFFER_UNSUPPORTED_EXT:
      cout << " Framebuffer status is unsupported!" << endl;
      break;
    default:
      cout << " Framebuffer status is not complete; an error occurred: ";
      cout << (int)fboStatus << endl;
      break;
  }
  return true;
}

// Prints OpenGL error code, and returns true if an error occurred, false otherwise
bool GLSL_SCENE_VIEWER::reportGLError(GLenum error)
{
  switch (error)
  {
    case GL_NO_ERROR:
      //cout << " OpenGL did not report any error." << endl;
      return false;
      break;
    case GL_INVALID_ENUM:
      cout << " OpenGL Error: Invalid Enum!" << endl;
      break;
    case GL_INVALID_VALUE:
      cout << " OpenGL Error: Invalid Value!" << endl;
      break;
    case GL_INVALID_OPERATION:
      cout << " OpenGL Error: Invalid Operation!" << endl;
      break;
    case GL_STACK_OVERFLOW:
      cout << " OpenGL Error: Stack Overflow!" << endl;
      break;
    case GL_STACK_UNDERFLOW:
      cout << " OpenGL Error: Stack Underflow!" << endl;
      break;
    case GL_OUT_OF_MEMORY:
      cout << " OpenGL Error: Out of Memory!" << endl;
      break;
    case GL_TABLE_TOO_LARGE:
      cout << " OpenGL Error: Table Too Large!" << endl;
      break;
    default:
      cout << " OpenGL Error: Unknown Error!" << endl;
      break;
  }
  return true;
}

