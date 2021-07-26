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
// TET.cpp: implementation of the TET class.
//
//////////////////////////////////////////////////////////////////////

#include "TET.h"

// DEBUG
bool TET::verbose = false;

//////////////////////////////////////////////////////////////////////
// Constructor for tetrahedra
//////////////////////////////////////////////////////////////////////
TET::TET()
{
  vertices[0] = NULL;
  vertices[1] = NULL;
  vertices[2] = NULL;
  vertices[3] = NULL;
  partition = 0;
  _materialIndex = 0;

  _skinningMatrix[0] = NULL;
  _skinningMatrix[1] = NULL;
  _skinningMatrix[2] = NULL;
  _skinningMatrix[3] = NULL;
  deformedToRest = NULL;
}

TET::TET(VEC3F* v0, VEC3F* v1, VEC3F* v2, VEC3F* v3) :
  _materialIndex(0)
{
  vertices[0] = v0;
  vertices[1] = v1;
  vertices[2] = v2;
  vertices[3] = v3;
  partition = 0;

  _skinningMatrix[0] = NULL;
  _skinningMatrix[1] = NULL;
  _skinningMatrix[2] = NULL;
  _skinningMatrix[3] = NULL;
  deformedToRest = NULL;

  init();
}

TET::TET(VEC3F* v[4]) :
  _materialIndex(0)
{
  vertices[0] = v[0];
  vertices[1] = v[1];
  vertices[2] = v[2];
  vertices[3] = v[3];
  partition = 0;

  _skinningMatrix[0] = NULL;
  _skinningMatrix[1] = NULL;
  _skinningMatrix[2] = NULL;
  _skinningMatrix[3] = NULL;
  deformedToRest = NULL;

  init();
}

TET::TET(const TET& tet)
{
  vertices[0] = tet.vertices[0];
  vertices[1] = tet.vertices[1];
  vertices[2] = tet.vertices[2];
  vertices[3] = tet.vertices[3];
  partition = 0;
  deformedToRest = tet.deformedToRest;

  _materialIndex = tet._materialIndex;

  _skinningMatrix[0] = tet._skinningMatrix[0];
  _skinningMatrix[1] = tet._skinningMatrix[1];
  _skinningMatrix[2] = tet._skinningMatrix[2];
  _skinningMatrix[3] = tet._skinningMatrix[3];

  init();
}

TET& TET::operator=(const TET& tet)
{
  vertices[0] = tet.vertices[0];
  vertices[1] = tet.vertices[1];
  vertices[2] = tet.vertices[2];
  vertices[3] = tet.vertices[3];
  partition = 0;
  deformedToRest = tet.deformedToRest;

  _materialIndex = tet._materialIndex;

  _skinningMatrix[0] = tet._skinningMatrix[0];
  _skinningMatrix[1] = tet._skinningMatrix[1];
  _skinningMatrix[2] = tet._skinningMatrix[2];
  _skinningMatrix[3] = tet._skinningMatrix[3];

  init();

  return *this;
}

void TET::init()
{
  // if tet is uninitialized, bail
  if (vertices[0] == NULL) return;

  vector<TRIANGLE> faces;
  faces.push_back(TRIANGLE(vertices[0], vertices[1], vertices[3]));
  faces.push_back(TRIANGLE(vertices[0], vertices[2], vertices[1]));
  faces.push_back(TRIANGLE(vertices[3], vertices[2], vertices[0]));
  faces.push_back(TRIANGLE(vertices[1], vertices[2], vertices[3]));

  // cache material inverse - following Teran's notation
	_DmInv = MATRIX3(*vertices[1] - *vertices[0],
               		 *vertices[2] - *vertices[0],
               		 *vertices[3] - *vertices[0]);
  _DmInv = _DmInv.transpose();
  _DmInv = _DmInv.inverse();

  // calculate the area vectors
  // v0 is incident on faces (0,1,2)
  _b[0] = faces[0].normal() * faces[0].area() +
          faces[1].normal() * faces[1].area() +
          faces[2].normal() * faces[2].area();
  
  // v1 is incident on faces (0,1,3)
  _b[1] = faces[0].normal() * faces[0].area() +
          faces[1].normal() * faces[1].area() +
          faces[3].normal() * faces[3].area();
  
  // v2 is incident on faces (1,2,3)
  _b[2] = faces[1].normal() * faces[1].area() +
          faces[2].normal() * faces[2].area() +
          faces[3].normal() * faces[3].area();
  
  // v3 is incident on faces (0,2,3)
  _b[3] = faces[0].normal() * faces[0].area() +
          faces[2].normal() * faces[2].area() +
          faces[3].normal() * faces[3].area();
  _b[0] *= -1.0f / 3.0f;
  _b[1] *= -1.0f / 3.0f;
  _b[2] *= -1.0f / 3.0f;
  _b[3] *= -1.0f / 3.0f;
}

//////////////////////////////////////////////////////////////////////
// set the GL normal
//////////////////////////////////////////////////////////////////////
void TET::drawNormal(int v0, int v1, int v2)
{
  VEC3F normal = cross(*vertices[v1] - *vertices[v0],
                       *vertices[v2] - *vertices[v0]);
  normal.normalize();

  // PetSc doesn't like glNormal3fv
#ifdef SINGLE_PRECISION
  glNormal3f(normal[0], normal[1], normal[2]);
#else
  glNormal3d(normal[0], normal[1], normal[2]);
#endif
}

//////////////////////////////////////////////////////////////////////
// compute the tet volume
//////////////////////////////////////////////////////////////////////
Real TET::volume()
{
  // formula for a tet volume with vertices (a,b,c,d) is:
  // |(a - d) dot ((b - d) cross (c - d))| / 6
  VEC3F a = (*vertices[1]) - (*vertices[0]);
  VEC3F b = (*vertices[2]) - (*vertices[0]);
  VEC3F c = (*vertices[3]) - (*vertices[0]);

  VEC3F crossProduct = cross(b,c);

  return fabs(crossProduct * a) / 6.0f;
}

//////////////////////////////////////////////////////////////////////
// see if any of the normal point in the wrong direction
//////////////////////////////////////////////////////////////////////
bool TET::invalid()
{
  VEC3F edges[6];
  edges[0] = *(vertices[1]) - *(vertices[0]);
  edges[1] = *(vertices[2]) - *(vertices[0]);
  edges[2] = *(vertices[3]) - *(vertices[0]);
  edges[3] = *(vertices[1]) - *(vertices[2]);
  edges[4] = *(vertices[2]) - *(vertices[3]);
  edges[5] = *(vertices[1]) - *(vertices[3]);

  // LEFT handed coordinate system
  VEC3F normals[4];
  normals[0] = cross(edges[0], edges[2]);
  normals[1] = cross(edges[5], edges[4]);
  normals[2] = cross(edges[1], edges[0]);
  normals[3] = cross(edges[2], edges[1]);

  // indices of a point on the face
  int origin[] = {0, 3, 0, 3};

  // indices of the one point not on the face
  int opposing[] = {2,0,3,1};

  for (int x = 0; x < 4; x++)
  {
    // calc a point in direction of normal and opposite
    VEC3F positive = *(vertices[origin[x]]) + normals[x];
    VEC3F negative = *(vertices[origin[x]]) - normals[x];

    // "negative" should now be closer to the opposing vertex
    Real distancePositive = norm(*(vertices[opposing[x]]) - positive);
    Real distanceNegative = norm(*(vertices[opposing[x]]) - negative);

    // if it is not, return that badness is true
    if (distanceNegative > distancePositive)
      return true;
  }
  return false;
}

//////////////////////////////////////////////////////////////////////
// get the dihedral angles
//////////////////////////////////////////////////////////////////////
void TET::dihedrals(Real* angles) {
  VEC3F edges[6];
  edges[0] = *(vertices[1]) - *(vertices[0]);
  edges[1] = *(vertices[2]) - *(vertices[0]);
  edges[2] = *(vertices[3]) - *(vertices[0]);
  edges[3] = *(vertices[1]) - *(vertices[2]);
  edges[4] = *(vertices[2]) - *(vertices[3]);
  edges[5] = *(vertices[1]) - *(vertices[3]);

  VEC3F normals[4];
  normals[0] = cross(edges[0], edges[2]);
  normals[1] = cross(edges[5], edges[4]);
  normals[2] = cross(edges[1], edges[0]);
  normals[3] = cross(edges[2], edges[1]);
  unitize(normals[0]);
  unitize(normals[1]);
  unitize(normals[2]);
  unitize(normals[3]);

  Real dots[6];
  dots[0] = normals[0] * -normals[2];
  dots[1] = normals[2] * -normals[3];
  dots[2] = normals[0] * -normals[3];
  dots[3] = normals[1] * -normals[2];
  dots[4] = normals[1] * -normals[3];
  dots[5] = normals[0] * -normals[1];

  angles[0] = acos(dots[0]);
  angles[1] = acos(dots[1]);
  angles[2] = acos(dots[2]);
  angles[3] = acos(dots[3]);
  angles[4] = acos(dots[4]);
  angles[5] = acos(dots[5]);
}

//////////////////////////////////////////////////////////////////////
// get the minimum of the dihedral angles
//////////////////////////////////////////////////////////////////////
Real TET::minDihedral()
{
  Real angles[6];
  dihedrals(&angles[0]);
  Real minFound = angles[0];
  for (int x = 1; x < 6; x++)
    minFound = (minFound < angles[x]) ? minFound : angles[x];
  return minFound;
}

//////////////////////////////////////////////////////////////////////
// get the maxmimum of the dihedral angles
//////////////////////////////////////////////////////////////////////
Real TET::maxDihedral()
{
  Real angles[6];
  dihedrals(&angles[0]);
  Real maxFound = angles[0];
  for (int x = 1; x < 6; x++)
    maxFound = (maxFound > angles[x]) ? maxFound : angles[x];
  return maxFound;
}

//////////////////////////////////////////////////////////////////////
// tetrahedral triangle drawing function
//////////////////////////////////////////////////////////////////////
void TET::drawTriangles() {
  glBegin(GL_TRIANGLES);
#ifdef SINGLE_PRECISION
    drawNormal(0,2,1);
    glVertex3fv(*vertices[0]);
    glVertex3fv(*vertices[2]);
    glVertex3fv(*vertices[1]);

    drawNormal(0,1,3);
    glVertex3fv(*vertices[0]);
    glVertex3fv(*vertices[1]);
    glVertex3fv(*vertices[3]);

    drawNormal(0,3,2);
    glVertex3fv(*vertices[0]);
    glVertex3fv(*vertices[3]);
    glVertex3fv(*vertices[2]);

    drawNormal(1,2,3);
    glVertex3fv(*vertices[1]);
    glVertex3fv(*vertices[2]);
    glVertex3fv(*vertices[3]);
#else
    drawNormal(0,2,1);
    glVertex3dv(*vertices[0]);
    glVertex3dv(*vertices[2]);
    glVertex3dv(*vertices[1]);

    drawNormal(0,1,3);
    glVertex3dv(*vertices[0]);
    glVertex3dv(*vertices[1]);
    glVertex3dv(*vertices[3]);

    drawNormal(0,3,2);
    glVertex3dv(*vertices[0]);
    glVertex3dv(*vertices[3]);
    glVertex3dv(*vertices[2]);

    drawNormal(1,2,3);
    glVertex3dv(*vertices[1]);
    glVertex3dv(*vertices[2]);
    glVertex3dv(*vertices[3]);
#endif
  glEnd();
}

//////////////////////////////////////////////////////////////////////
// return the triangle corresponding to a face
//////////////////////////////////////////////////////////////////////
TRIANGLE TET::face(int x) {
  if (x == 0) return TRIANGLE(vertices[0], vertices[1], vertices[3]);
  if (x == 1) return TRIANGLE(vertices[0], vertices[2], vertices[1]);
  if (x == 2) return TRIANGLE(vertices[3], vertices[2], vertices[0]);
  return TRIANGLE(vertices[1], vertices[2], vertices[3]);
}

//////////////////////////////////////////////////////////////////////
// draw the stored faces with precomputed normals
//////////////////////////////////////////////////////////////////////
void TET::drawFaces() {
  vector<TRIANGLE> faces;
  faces.push_back(TRIANGLE(vertices[0], vertices[1], vertices[3]));
  faces.push_back(TRIANGLE(vertices[0], vertices[2], vertices[1]));
  faces.push_back(TRIANGLE(vertices[3], vertices[2], vertices[0]));
  faces.push_back(TRIANGLE(vertices[1], vertices[2], vertices[3]));
  
  for (int x = 0; x < 4; x++)
    faces[x].draw();
}

//////////////////////////////////////////////////////////////////////
// tetrahedral line drawing function
//////////////////////////////////////////////////////////////////////
void TET::drawLines() {
  glBegin(GL_LINE_STRIP);
#ifdef SINGLE_PRECISION
    glVertex3fv(*vertices[0]);
    glVertex3fv(*vertices[2]);
    glVertex3fv(*vertices[1]);
    glVertex3fv(*vertices[0]);

    glVertex3fv(*vertices[0]);
    glVertex3fv(*vertices[1]);
    glVertex3fv(*vertices[3]);
    glVertex3fv(*vertices[0]);

    glVertex3fv(*vertices[0]);
    glVertex3fv(*vertices[3]);
    glVertex3fv(*vertices[2]);
    glVertex3fv(*vertices[0]);

    glVertex3fv(*vertices[1]);
    glVertex3fv(*vertices[2]);
    glVertex3fv(*vertices[3]);
    glVertex3fv(*vertices[1]);
#else
    glVertex3dv(*vertices[0]);
    glVertex3dv(*vertices[2]);
    glVertex3dv(*vertices[1]);
    glVertex3dv(*vertices[0]);

    glVertex3dv(*vertices[0]);
    glVertex3dv(*vertices[1]);
    glVertex3dv(*vertices[3]);
    glVertex3dv(*vertices[0]);

    glVertex3dv(*vertices[0]);
    glVertex3dv(*vertices[3]);
    glVertex3dv(*vertices[2]);
    glVertex3dv(*vertices[0]);

    glVertex3dv(*vertices[1]);
    glVertex3dv(*vertices[2]);
    glVertex3dv(*vertices[3]);
    glVertex3dv(*vertices[1]);
#endif
  glEnd();
}

//////////////////////////////////////////////////////////////////////
// Deformation gradient - F seems standard, both Bonet and Wood 
// and Teran use it
//////////////////////////////////////////////////////////////////////
MATRIX3 TET::F() {
  // if there's not skinning, do the computation as usual
  if (_skinningMatrix[0] == NULL)
  {
    if (verbose)
    {
      cout << __FILE__ << " " << __LINE__ << " : " << endl;
      cout << " No skinning F" << endl;
    }

    // Following Teran's notation - spatial tensor
    MATRIX3 Ds = MATRIX3(*(vertices[1]) - *(vertices[0]),
                         *(vertices[2]) - *(vertices[0]),
                         *(vertices[3]) - *(vertices[0]));

    return Ds.transpose() * _DmInv;
  }

  if (verbose)
  {
    cout << __FILE__ << " " << __LINE__ << " : " << endl;
    cout << " Skinning F" << endl;
  }

  // get all the rest versions of the current vertices
  map<VEC3F*, VEC3F*>& translate = *deformedToRest;
  VEC3F* restVertices[4];
  for (int x = 0; x < 4; x++)
    restVertices[x] = translate[vertices[x]];

  // repack everything into a homogeneous vector
  VECTOR homogeneous(4);
  homogeneous[3] = 1.0;

  // get the initial positions
  VEC3F v0 = *(restVertices[0]);
  VEC3F v1 = *(restVertices[1]);
  VEC3F v2 = *(restVertices[2]);
  VEC3F v3 = *(restVertices[3]);

  VEC3F diff0 = (*(vertices[0])) - v0;
  VEC3F diff1 = (*(vertices[1])) - v1;
  VEC3F diff2 = (*(vertices[2])) - v2;
  VEC3F diff3 = (*(vertices[3])) - v3;
  
  homogeneous[0] = v0[0];
  homogeneous[1] = v0[1];
  homogeneous[2] = v0[2];
  homogeneous = (*(_skinningMatrix[0])) * homogeneous;
  v0[0] = homogeneous[0];
  v0[1] = homogeneous[1];
  v0[2] = homogeneous[2];

  homogeneous[0] = v1[0];
  homogeneous[1] = v1[1];
  homogeneous[2] = v1[2];
  homogeneous = (*(_skinningMatrix[1])) * homogeneous;
  v1[0] = homogeneous[0];
  v1[1] = homogeneous[1];
  v1[2] = homogeneous[2];

  homogeneous[0] = v2[0];
  homogeneous[1] = v2[1];
  homogeneous[2] = v2[2];
  homogeneous = (*(_skinningMatrix[2])) * homogeneous;
  v2[0] = homogeneous[0];
  v2[1] = homogeneous[1];
  v2[2] = homogeneous[2];

  homogeneous[0] = v3[0];
  homogeneous[1] = v3[1];
  homogeneous[2] = v3[2];
  homogeneous = (*(_skinningMatrix[3])) * homogeneous;
  v3[0] = homogeneous[0];
  v3[1] = homogeneous[1];
  v3[2] = homogeneous[2];

  v0 += diff0;
  v1 += diff1;
  v2 += diff2;
  v3 += diff3;

  // Following Teran's notation - spatial tensor
  MATRIX3 Ds = MATRIX3(v1 - v0, v2 - v0, v3 - v0);
  return Ds.transpose() * _DmInv;

  /*
  // repack everything into a homogeneous vector
  VECTOR homogeneous(4);
  homogeneous[3] = 1.0;

  // get the initial positions
  VEC3F v0 = *(vertices[0]);
  VEC3F v1 = *(vertices[1]);
  VEC3F v2 = *(vertices[2]);
  VEC3F v3 = *(vertices[3]);

  homogeneous[0] = v0[0];
  homogeneous[1] = v0[1];
  homogeneous[2] = v0[2];
  homogeneous = (*(_skinningMatrix[0])) * homogeneous;
  v0[0] = homogeneous[0];
  v0[1] = homogeneous[1];
  v0[2] = homogeneous[2];

  homogeneous[0] = v1[0];
  homogeneous[1] = v1[1];
  homogeneous[2] = v1[2];
  homogeneous = (*(_skinningMatrix[1])) * homogeneous;
  v1[0] = homogeneous[0];
  v1[1] = homogeneous[1];
  v1[2] = homogeneous[2];

  homogeneous[0] = v2[0];
  homogeneous[1] = v2[1];
  homogeneous[2] = v2[2];
  homogeneous = (*(_skinningMatrix[2])) * homogeneous;
  v2[0] = homogeneous[0];
  v2[1] = homogeneous[1];
  v2[2] = homogeneous[2];

  homogeneous[0] = v3[0];
  homogeneous[1] = v3[1];
  homogeneous[2] = v3[2];
  homogeneous = (*(_skinningMatrix[3])) * homogeneous;
  v3[0] = homogeneous[0];
  v3[1] = homogeneous[1];
  v3[2] = homogeneous[2];

  // Following Teran's notation - spatial tensor
  MATRIX3 Ds = MATRIX3(v1 - v0, v2 - v0, v3 - v0);
  return Ds.transpose() * _DmInv;
  */
}

//////////////////////////////////////////////////////////////////////
// flatten the deformation gradient out into a correctly ordered 
// vector
//////////////////////////////////////////////////////////////////////
VECTOR TET::flattenF(MATRIX3& F)
{
  VECTOR final(9);
  final(0) = F(0,0);
  final(1) = F(1,0);
  final(2) = F(2,0);
  final(3) = F(0,1);
  final(4) = F(1,1);
  final(5) = F(2,1);
  final(6) = F(0,2);
  final(7) = F(1,2);
  final(8) = F(2,2);

  return final;
}

//////////////////////////////////////////////////////////////////////
// flatten the deformation gradient out into a correctly ordered 
// vector
//////////////////////////////////////////////////////////////////////
void TET::flattenF(MATRIX3& F, VECTOR &output)
{
  if ( output.size() != 9 )
    output.resizeAndWipe( 9 );

  output(0) = F(0,0);
  output(1) = F(1,0);
  output(2) = F(2,0);
  output(3) = F(0,1);
  output(4) = F(1,1);
  output(5) = F(2,1);
  output(6) = F(0,2);
  output(7) = F(1,2);
  output(8) = F(2,2);
}

//////////////////////////////////////////////////////////////////////
// repack the deformation gradient into a correctly ordered matrix
//////////////////////////////////////////////////////////////////////
MATRIX3 TET::repackF(VECTOR& F)
{
  MATRIX3 final;
  final(0,0) = F(0);
  final(1,0) = F(1);
  final(2,0) = F(2);
  final(0,1) = F(3);
  final(1,1) = F(4);
  final(2,1) = F(5);
  final(0,2) = F(6);
  final(1,2) = F(7);
  final(2,2) = F(8);

  return final;
}

//////////////////////////////////////////////////////////////////////
// repack the deformation gradient into a correctly ordered matrix
//////////////////////////////////////////////////////////////////////
MATRIX3 TET::repackF(const Real* F)
{
  MATRIX3 final;
  final(0,0) = F[0];
  final(1,0) = F[1];
  final(2,0) = F[2];
  final(0,1) = F[3];
  final(1,1) = F[4];
  final(2,1) = F[5];
  final(0,2) = F[6];
  final(1,2) = F[7];
  final(2,2) = F[8];

  return final;
}

//////////////////////////////////////////////////////////////////////
// Converts a block column vector of 3x3 matrices to a block row
// vector of the same 3x3 matrices.
//////////////////////////////////////////////////////////////////////
void TET::repackBlockColumnToRow( const Real *input, Real *output,
                                  int numBlocks )
{
  int inputOffset, outputOffset;

  for ( int i = 0; i < numBlocks; i++ )
  {
    // Copy each row - this assumes row-major storage
    for ( int j = 0; j < 3; j++ )
    {
      // Figure out row-major format offsets
      inputOffset = i * 9 + j * 3;
      outputOffset = j * numBlocks * 3 + i * 3;

      memcpy( (void *)( output + outputOffset ), input + inputOffset,
              3 * sizeof( Real ) );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Same as the above, but row to column format
//////////////////////////////////////////////////////////////////////
void TET::repackBlockRowToColumn( const Real *input, Real *output,
                                  int numBlocks )
{
  int inputOffset, outputOffset;

  for ( int i = 0; i < numBlocks; i++ )
  {
    // Copy each row - this assumes row-major storage
    for ( int j = 0; j < 3; j++ )
    {
      // Figure out row-major format offsets
      inputOffset = j * numBlocks * 3 + i * 3;
      outputOffset = i * 9 + j * 3;

      memcpy( (void *)( output + outputOffset ), input + inputOffset,
              3 * sizeof( Real ) );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Flattens a block row vector of numBlock 3x3 matrices in to a
// 9 by numBlock matrix
//////////////////////////////////////////////////////////////////////
void TET::flattenBlockRow( const Real *input, Real *output,
                           int numBlocks )
{
  int inputCols = 3 * numBlocks;
  int colIndex, blockIndex;

  for ( int i = 0; i < 3; i++ )
  for ( int j = 0; j < inputCols; j++ )
  {
    colIndex = j % 3;
    blockIndex = j / 3;

    MATRIX::access( output, 9, numBlocks, colIndex * 3 + i, blockIndex )
      = MATRIX::access( input, 3, inputCols, i, j );
  }
}

//////////////////////////////////////////////////////////////////////
// Repacks a 9 x numBlock matrix in to a block column vector
// of 3x3 matrices
//////////////////////////////////////////////////////////////////////
void TET::repackBlockColumn( const Real *input, Real *output,
                             int numBlocks )
{
  int outputRows = 3 * numBlocks;

  Real *blockData;

  for ( int i = 0; i < 9; i++ )
  for ( int j = 0; j < numBlocks; j++ )
  {
    blockData = output + j * 9;

    MATRIX::access( blockData, 3, 3, i % 3, i / 3 )
      = MATRIX::access( input, 9, numBlocks, i, j );
  }
}

//////////////////////////////////////////////////////////////////////
// Center of the tet
//////////////////////////////////////////////////////////////////////
VEC3F TET::center()
{
  VEC3F sum;
  sum += *(vertices[0]);
  sum += *(vertices[1]);
  sum += *(vertices[2]);
  sum += *(vertices[3]);
  sum *= 0.25;
  return sum;
}

//////////////////////////////////////////////////////////////////////
// Flatten state out into a vector
//////////////////////////////////////////////////////////////////////
VECTOR TET::position()
{
  VECTOR state(12);

  int i = 0;
  for (int x = 0; x < 4; x++)
    for (int y = 0; y < 3; y++, i++)
      state[i] =(*(vertices[x]))[y];

  return state;
}

//////////////////////////////////////////////////////////////////////
// Do these two tets share a face?
//////////////////////////////////////////////////////////////////////
int TET::sharedFace(TET& otherTet)
{
  // generate the triangles
  TRIANGLE thisFaces[4];
  TRIANGLE otherFaces[4];
  for (int x = 0; x < 4; x++)
  {
    thisFaces[x] = this->face(x);
    otherFaces[x] = otherTet.face(x);
  }

  for (int x = 0; x < 4; x++)
    for (int y = 0; y < 4; y++)
    {
      if (thisFaces[x] == otherFaces[y])
      {
        //cout << " this tet:   " << this << endl;
        //cout << " other tet:  " << &otherTet << endl;
        //cout << " this face:  " << thisFaces[x] << endl;
        //cout << " other face: " << otherFaces[y] << endl;
        return x;
      }
    }

  //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  return -1;
}
