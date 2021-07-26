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
// GLSL_TET_MESH.h: Interface for the GLSL_TET_MESH class.
//
//////////////////////////////////////////////////////////////////////

#ifndef GLSL_TET_MESH_H
#define GLSL_TET_MESH_H

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
#include <SUBSPACE_TET_MESH.h>
#include <GLSL_SHADER.h>

using namespace std;

//////////////////////////////////////////////////////////////////////
// OpenGL GLSL tet mesh (from SUBSPACE_TET_MESH)
//////////////////////////////////////////////////////////////////////
class GLSL_TET_MESH {

public:
  GLSL_TET_MESH(SUBSPACE_TET_MESH *tetMesh, OBJ *embeddedMesh);
  virtual ~GLSL_TET_MESH();

  // accessors
  SUBSPACE_TET_MESH *tetMesh();
  GLSL_SHADER *tetVertexShader();
  GLSL_SHADER *genRotationShader();
  GLSL_SHADER *drawBasesShader();
  GLSL_SHADER *drawFacesShader();
  GLSL_SHADER *drawNeighborsShader();
  GLSL_SHADER *drawAOVShader();
  GLSL_SHADER *drawAOV2Shader();
  GLSL_SHADER *drawAOV3Shader();

  GLSL_SHADER *renderFromTextureTest();

  GLuint tetMeshVboID();
  GLuint tetMeshInfoVboID();
  GLsizei tetMeshSize();

  GLuint tetFaceInfo0VboID();
  GLuint tetFaceInfo1VboID();
  GLsizei tetFaceSize();

  GLuint embeddedMeshVboID();
  GLuint embeddedNormalVboID();
  GLuint embeddedVertLoc0VboID();
  GLuint embeddedVertLoc1VboID();
  GLuint embeddedVertLoc2VboID();
  GLuint embeddedMeshIboID();
  GLsizei embeddedMeshIboSize();
  GLsizei embeddedVertSize();

  GLuint surfaceVertVboID();
  GLuint surfaceVertNeighbor0VboID();
  GLuint surfaceVertNeighbor1VboID();
  GLsizei surfaceVertSize();

  GLuint n1RestTexID();
  GLuint n2RestTexID();
  
  GLuint UcoordVboID();
  GLuint UbasisTexID();

  void getVertexTexDims(GLsizei *width, GLsizei *height,
                        GLsizei *embeddedwidth, GLsizei *embeddedheight);
  void constructDeltaVector(OBJ *embeddedMesh, Real maxDelta,
                            Real threshold, int ra_max);
  
  void debugDraw(int pose, bool drawdeltas = false);
  bool drawEmbedded();
  void debugTestDeltaFunctions();

protected:
   // pointer to the mesh
  SUBSPACE_TET_MESH *_tetMesh;

  // mesh in its resting pose
  GLuint _tetMeshVboID;
  // holds the vertex info for constructing _tetVertexTexID
  GLuint _tetMeshInfoVboID;
  GLsizei _tetMeshSize;

  GLsizei _tetVertexTexWidth;
  GLsizei _tetVertexTexHeight;
  GLsizei _embeddedVertexTexWidth;
  GLsizei _embeddedVertexTexHeight;
      
  // holds the mesh face info needed for constructing _tetNormalTexID and _tetTangentTexID
  GLuint _tetFaceInfo0VboID;      // holds tex coords for v0, v1
  GLuint _tetFaceInfo1VboID;      // holds tex coords for v2
  GLsizei _tetFaceSize;

  // vertex buffer to hold the baricentric coords of the embedded mesh
  GLuint _embeddedMeshVboID;
   // vertex buffer to hold vertex normal 
  GLuint _embeddedNormalVboID;
  // Two VBO's to hold four sets of tex-coords into _tetVertexTexID:
  GLuint _embeddedVertLoc0VboID;  // Holds v0 (.st) and v1 (.pq)
  GLuint _embeddedVertLoc1VboID;  // Holds v2 (.st) and v3 (.pq)
  GLuint _embeddedVertLoc2VboID;  // Holds .st of embedded vertex
   // Index buffer to hold the correct index info for rendering embedded mesh
  GLuint _embeddedMeshIboID;
  // Size of embedded mesh index buffer
  GLsizei _embeddedMeshIboSize;
  GLsizei _embeddedVertSize;
  
  // texcoords to index into _UbasisTexID
  GLuint _UcoordVboID;

  // texture object ID which holds the U basis matrix
  GLuint _UbasisTexID;
  GLsizei _UbasisTexWidth;
  GLsizei _UbasisTexHeight;
  // Holds the number of columns of U from tetMesh, padded to be multiple of 4
  GLsizei _UmatrixPaddedWidth;

  // qVecSize holds the number of vec4's it expects in q (i.e. ceil(q.size() / 4))
  GLuint _qVecSize;

  // Store surface vertex texcoords for rotation matrix generation
  GLuint _surfaceVertVboID;
  GLuint _surfaceVertNeighbor0VboID;
  GLuint _surfaceVertNeighbor1VboID;
  GLsizei _surfaceVertSize;

  // Rest frame textures (used for creating rotation matrix)
  GLuint _n1RestTexID;
  GLuint _n2RestTexID; // */

  GLSL_SHADER *_tetVertexShader;
  GLSL_SHADER *_genRotationShader;
  GLSL_SHADER *_drawBasesShader;
  GLSL_SHADER *_drawFacesShader;
  GLSL_SHADER *_drawNeighborsShader;
  GLSL_SHADER *_drawAOVShader;
  GLSL_SHADER *_drawAOV2Shader;
  GLSL_SHADER *_drawAOV3Shader;

  GLSL_SHADER *_renderFromTextureTest;
 
  void orderTet(TET& tet, VEC3F& normal, int *returnOrder);
  void orderNeighbors(VEC3F* center, vector<int>& neighbors, int *returnOrder);

  // Helper methods called in the constructor to initialize the many variables and states of the class
  void constructTetMeshVbo();
  void constructUBasisTex();
  void constructEmbeddedMeshBuffers(OBJ *embeddedMesh);
  void constructTetFaceInfo();
  
  void genRestFrame();

  // All these variables are TEMPORARY, storing info for drawing
  // the embedded model and calculated deltas
  MATRIX *_Dv;
  MATRIX *_Df;
  MATRIX *_Q;
  OBJ *_embeddedMesh;
  int _prevPose;

  Real calculateDelta(TRIANGLE &T, Real delta_max, Real threshold, OBJ &obj);
  bool isPointInVolume(VEC3F& x, TRIANGLE& t, vector<VEC3F>& m);
  Real calculateAO(VEC3F& x, VEC3F& n, TRIANGLE& t, bool debugmode=false);
  Real calculatePointTriDistance(VEC3F& x, TRIANGLE& t, vector<VEC3F>& m, Real delta);
};

#endif // GLSL_TET_MESH_H

