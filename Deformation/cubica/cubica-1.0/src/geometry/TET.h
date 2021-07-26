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
// TET.h: interface for the TET class.
//
//////////////////////////////////////////////////////////////////////

#ifndef TET_H
#define TET_H

#include <SETTINGS.h>
#include <VEC3.h>
#include <MATRIX3.h>
#include <TRIANGLE.h>
#include <VECTOR.h>
#include <MATRIX.h>
#if WIN32
#include <gl/glut.h>
#elif USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

using namespace std;

//////////////////////////////////////////////////////////////////////
// Tetrahedron class
//
// Note these use a LEFT handed coordinate system (like in screenspace)
//////////////////////////////////////////////////////////////////////
class TET {
public:
  TET(VEC3F* v0, VEC3F* v1, VEC3F* v2, VEC3F* v3);
  TET(VEC3F* v[4]);
  TET();
  TET(const TET& tet);
  TET& operator=(const TET& tet);

  // the only public variables, since everybody uses them so much
  // -- this is the biggest deviation from the "no public variables"
  // rule
  VEC3F* vertices[4];
  //vector<TRIANGLE> faces;
  TRIANGLE face(int x);
  
  int partition;

  int& materialIndex() { return _materialIndex; };

  // GL drawing functions
  void drawFaces();
  void drawTriangles();
  void drawLines();

  // tet volume
  Real volume();

  // check that the winding order of the tet is valid
  bool invalid();

  // dihedral angle statistics
  void dihedrals(Real* angles);
  Real minDihedral();
  Real maxDihedral();

  // deformation gradient - F seems standard, both 
  // Bonet and Wood and Teran use it
  MATRIX3 F();

  // area vectors
  const VEC3F* b() { return _b; };

  // materal inverse
  const MATRIX3 DmInv() { return _DmInv; };

  // see if the tet is inverted
  bool inverted() { return det(F()) < 0.0f; };

  // flatten the deformation gradient out into a correctly ordered vector
  static VECTOR flattenF(MATRIX3& F);
  static void flattenF(MATRIX3& F, VECTOR &output);
  
  // repack the deformation gradient into a correctly ordered matrix
  static MATRIX3 repackF(VECTOR& F);
  static MATRIX3 repackF(const Real* F);

  // Converts a block column vector of 3x3 matrices to a block row
  // vector of the same 3x3 matrices.
  static void repackBlockColumnToRow( const Real *input, Real *output,
                                      int numBlocks );

  // Same as the above, but row to column format
  static void repackBlockRowToColumn( const Real *input, Real *output,
                                      int numBlocks );

  // Flattens a block row vector of numBlock 3x3 matrices in to a
  // 9 by numBlock matrix
  static void flattenBlockRow( const Real *input, Real *output,
                               int numBlocks );

  // Repacks a 9 x numBlock matrix in to a block column vector
  // of 3x3 matrices
  static void repackBlockColumn( const Real *input, Real *output,
                                 int numBlocks );

  // accessor for skinning matrices
  MATRIX*& skinningMatrix(int index) { return _skinningMatrix[index]; };
  
  // average of all the vertices
  VEC3F center();

  // get the rest version of the deformed vertex 
  map<VEC3F*, VEC3F*>* deformedToRest;

  // initialize various variables
  void init();

  // flatten state out into a vector
  VECTOR position();

  // do these two tets share a face?
  // -1 if they do not,
  // face # in this if they do
  int sharedFace(TET& otherTet);

  // DEBUG
  static bool verbose;

private:
  // material inverse - following Teran's notation
  MATRIX3 _DmInv;

  // area vectors - following Teran's notation
  VEC3F _b[4];
  
  // GL drawing functions
  void drawNormal(int v0, int v1, int v2);

  // skinning support
  MATRIX* _skinningMatrix[4];
  
  // index into the TET_MESH::_materials array
  int _materialIndex;
};

//////////////////////////////////////////////////////////////////////
// dump tet vertices to iostream
//////////////////////////////////////////////////////////////////////
inline ostream &operator<<(ostream &out, TET& tet)
{ 
  for (int x = 0; x < 4; x++)
    out << *(tet.vertices[x]) << endl;
  out << endl;
  return out;
}

#endif
