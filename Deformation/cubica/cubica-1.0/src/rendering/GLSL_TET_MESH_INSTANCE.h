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
// GLSL_TET_MESH_INSTANCE.h: Interface for the GLSL_TET_MESH_INSTANCE class.
//
///////////

#ifndef GLSL_TET_MESH_INSTANCE_H
#define GLSL_TET_MESH_INSTANCE_H

#include <GLSL_TET_MESH.h>
#include <SUBSPACE_INTEGRATOR.h>

enum RenderStyleType { 
  RENDER_STYLE_SHADOW = 0,
  RENDER_STYLE_NO_SHADOW = 1,
  RENDER_STYLE_SM = 2,
  RENDER_STYLE_DEFERRED = 3,
  RENDER_STYLE_NORMALS = 4,
  RENDER_STYLE_BASIC = 5
};

enum IntegratorOptionType {
  INTEGRATOR_STEP_DEFAULT = 0,
  INTEGRATOR_STEP_QUASISTATIC = 1,
  INTEGRATOR_STEP_IMPLICIT = 2,
  INTEGRATOR_STEP_INVERTIBLE_IMPLICIT = 3,
  INTEGRATOR_STEP_INVERTIBLE_QUASISTATIC = 4,
  INTEGRATOR_STEP_EXPLICIT = 5,
  INTEGRATOR_STEP_INVERTIBLE_SEMI_IMPLICIT = 6
};

/*#define RENDER_STYLE_SHADOW     0
#define RENDER_STYLE_NO_SHADOW  1
#define RENDER_STYLE_BASIC      2
#define RENDER_STYLE_DEFERRED   3
#define RENDER_STYLE_NORMALS   4*/


class GLSL_TET_MESH_INSTANCE
{
public:
  GLSL_TET_MESH_INSTANCE(GLSL_TET_MESH *mesh, string configFileName);
  virtual ~GLSL_TET_MESH_INSTANCE();

  int startPoint;

  // Don't have surfaceIbo anymore
  //void drawRestPose();
  //void drawMeshDeformed();
  //void drawEmbeddedMeshDeformed();
  void renderWithShadows()    { drawEmbeddedMeshDeformed(RENDER_STYLE_SHADOW);     };
  void renderNoShadows()      { drawEmbeddedMeshDeformed(RENDER_STYLE_NO_SHADOW);  };
  void renderSM()             { drawEmbeddedMeshDeformed(RENDER_STYLE_SM);         };
  void renderDeferred()       { drawEmbeddedMeshDeformed(RENDER_STYLE_DEFERRED);   };
  void renderBasic()          { drawEmbeddedMeshDeformed(RENDER_STYLE_BASIC);      };
  //void renderNormals()        { drawEmbeddedMeshDeformed(RENDER_STYLE_NORMALS);    };

  void renderAOV(float delta=0.1f);
  void renderAOVEmbedded(float delta=0.1f, float width=2.0f, float height=2.0f);
  void renderAOVWireframe(float delta=0.1f);
  
  void genTetVertexTexture();
  void integratorStep( bool debug = false );

  GLfloat* position()     { return _position;   };
  GLfloat* rotation()     { return _rotation;   };

  void setPosition(GLfloat *position);
  void setPosition(GLfloat x, GLfloat y, GLfloat z);
  void addPosition(GLfloat *delta);
  void addPosition(GLfloat dx, GLfloat dy, GLfloat dz);

  void setRotation(GLfloat *angles);
  void setRotation(GLfloat rx, GLfloat ry, GLfloat rz);
  void addRotation(GLfloat *delta);
  void addRotation(GLfloat drx, GLfloat dry, GLfloat drz);

  void setMeshRenderToTextureShader(GLSL_SHADER *shader);

  void setMeshRenderWithShadowShader(GLSL_SHADER *shader);
  void setMeshRenderNoShadowShader(GLSL_SHADER *shader);
  void setMeshRenderSMShader(GLSL_SHADER *shader);
  void setMeshRenderDeferredShader(GLSL_SHADER *shader);
  void setMeshRenderBasicShader(GLSL_SHADER *shader);
  void setMeshRenderNormalsShader(GLSL_SHADER *shader);
  
  void setTetVertexShader(GLSL_SHADER *shader);
  void setNormTanTexShader(GLSL_SHADER *shader, int passnum);

  void drawEmbeddedMeshFromTexture();

  bool integratorClick(VEC3F& point, Real maxRadius = -1.0);
  void integratorDrag(VEC3F& point);
  void integratorUnclick(); 

  IntegratorOptionType integratorOptions()      { return _integratorOptions; };
  void setIntegratorOptions(IntegratorOptionType options);

  GLSL_TET_MESH *mesh() { return _mesh; }

  SUBSPACE_INTEGRATOR *integrator() { return _integrator; }

protected:
  // Shared tetMesh object
  GLSL_TET_MESH *_mesh;
  // This instance's integrator
  SUBSPACE_INTEGRATOR *_integrator;

  // Bit vector set of options for the SUBSPACE_INTEGRATOR
  IntegratorOptionType _integratorOptions;
  
  // Some basic location/orientation info of this instance
  GLfloat _position[3];
  GLfloat _rotation[3];

  // Axis aligned bounding box of the unrotated, rest pose mesh (min/max x, min/max y, min/max z)
  GLfloat _aabb[6];

  // Shader Pointers. Not to be deleted from here
  GLSL_SHADER *_meshRenderToTextureShader;

  GLSL_SHADER *_meshRenderWithShadowShader;
  GLSL_SHADER *_meshRenderNoShadowShader;
  GLSL_SHADER *_meshRenderSMShader;
  GLSL_SHADER *_meshRenderDeferredShader;
  GLSL_SHADER *_meshRenderBasicShader;
  GLSL_SHADER *_meshRenderNormalsShader;

  GLSL_SHADER *_tetVertexShader;
  GLSL_SHADER *_tetNormTanTex0Shader;
  GLSL_SHADER *_tetNormTanTex1Shader;
  GLSL_SHADER *_tetNormTanTex2Shader;

  // tet vertex framebuffer object to hold the RTT _tetVertexTexID 
  GLuint _tetVertexFboID;  
  // texture object ID which holds deformed tet vertices
  GLuint _tetVertexTexID;
  // dimensions of _tetVertexTexID
  GLsizei _tetVertexTexWidth;
  GLsizei _tetVertexTexHeight;

  // embedded vertex framebuffer object to hold the RTT _embeddedVertexTexID and _embeddedNormalTexID
  GLuint _embeddedVertexFboID;  
  // texture object IDs which holds deformed embedded vertices and normals
  GLuint _embeddedVertexTexID;
  GLuint _embeddedNormalTexID;                                               
  // dimensions of _tetVertexTexID
  GLsizei _embeddedVertexTexWidth;
  GLsizei _embeddedVertexTexHeight;

  // framebuffer to hold the rotation textures
  GLuint _tetRotationFboID;  
  // Three textures to hold rotation matrix at every vertex. Same dimensions as _tetVertexTexID
  GLuint _tetRotationTexID[3];

  // window viewport size
  GLint _windowViewport[4];

  // local array to hold q, for passing it to the shader
  GLfloat *_q;

  // qVecSize holds the number of vec4's it expects in q (i.e. ceil(q.size() / 4))
  GLuint _qVecSize;

  void copyQForUniform();
  void constructVertexTextures();
  void oldGenTetVertexTexture();
  void oldDrawEmbeddedMeshDeformed(int renderStyle);
  void oldRenderAOVEmbedded(float delta, float width, float height);

  void drawEmbeddedMeshDeformed(RenderStyleType renderStyle);

  void renderBases();
  void renderFaces();


};


#endif // GLSL_TET_MESH_INSTANCE_H

