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
#include "ISO_STUFFER.h"
#include <float.h>

//////////////////////////////////////////////////////////////////////
// Constructor - dimensions are of center mesh, not corner mesh
//////////////////////////////////////////////////////////////////////
ISO_STUFFER::ISO_STUFFER(int xRes, int yRes, int zRes) :
  _xRes(xRes), _yRes(yRes), _zRes(zRes), _useExistingCaches(true), _totalTime(0)
{
  _slabSize = _xRes * _yRes;
  _maxRes = (_xRes > _yRes) ? _xRes : _yRes;
  _maxRes = (_maxRes > _zRes) ? _maxRes : _zRes;
  _outside = -2.0f * _maxRes;

  _blackLength = 1.0f / _maxRes;
  _redLength = sqrtf(3.0f * 0.25f * _blackLength * _blackLength);

  _totalCorners = (_xRes + 1) * (_yRes + 1) * (_zRes + 1);
  _totalCenters = _xRes * _yRes * _zRes;

  // row 6 from Table 1 of the paper
  //_alphaLong = 0.21509f;
  //_alphaShort = 0.35900f;
  //_proofMin = 6.4917;
  //_proofMax = 164.1013;

  // row 7 from Table 1 of the paper
  _alphaLong = 0.22383;
  _alphaShort = 0.39700;
  _proofMin = 7.6872;
  _proofMax = 168.0481;
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
ISO_STUFFER::~ISO_STUFFER()
{
  map<EDGE, VEC3F*>::iterator i;
  for (i = _cutPoints.begin(); i != _cutPoints.end(); ++i)
    delete i->second;

  map<int, VEC3F*>::iterator j;
  for (j = _vertices.begin(); j != _vertices.end(); ++j)
    delete i->second;

  for (unsigned int x = 0; x < _tets.size(); x++)
    delete _tets[x];
}

//////////////////////////////////////////////////////////////////////
// Draw all the tets as triangles
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::drawTets()
{
  for (unsigned int x = 0; x < _tets.size(); x++)
    _tets[x]->drawTriangles();
}

//////////////////////////////////////////////////////////////////////
// Draw all the surface tets as triangles
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::drawSurfaceTets()
{
  for (unsigned int x = 0; x < _surfaceTets.size(); x++)
    _surfaceTets[x]->drawTriangles();
}

//////////////////////////////////////////////////////////////////////
// Draw a z slice of a tet set
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::drawZSlice(Real zSlice, vector<TET*>& tets)
{
  srand(12345);

  for (unsigned int x = 0; x < tets.size(); x++)
  {
    Real zMin = (*tets[x]->vertices[0])[2];
    for (int y = 1; y < 4; y++)
      zMin = (zMin < (*tets[x]->vertices[y])[2]) ? zMin : (*tets[x]->vertices[y])[2];
    if (zMin < zSlice)
    {
      // color them randomly
      Real red = (Real)rand() / RAND_MAX;
      Real green = (Real)rand() / RAND_MAX;
      Real blue = (Real)rand() / RAND_MAX;

      glColor4f(red, green, blue, 1.0f);
      
      tets[x]->drawTriangles();
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Draw a z slice of a tet set
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::drawXSlice(Real xSlice, vector<TET*>& tets)
{
  srand(12345);

  for (unsigned int x = 0; x < tets.size(); x++)
  {
    Real xMax = (*tets[x]->vertices[0])[0];
    for (int y = 1; y < 4; y++)
      xMax = (xMax > (*tets[x]->vertices[y])[0]) ? xMax : (*tets[x]->vertices[y])[0];

    if (xMax < xSlice)
    {
      // color them randomly
      Real red = (Real)rand() / RAND_MAX;
      Real green = (Real)rand() / RAND_MAX;
      Real blue = (Real)rand() / RAND_MAX;

      glColor4f(red, green, blue, 1.0f);
      
      tets[x]->drawTriangles();
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Draw all the inside tets as triangles
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::drawInsideTets()
{
  for (unsigned int x = 0; x < _insideTets.size(); x++)
    _insideTets[x]->drawTriangles();
}

//////////////////////////////////////////////////////////////////////
// Draw all the final tets as triangles
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::drawFinalTets()
{
  srand(12345);

  for (unsigned int x = 0; x < _finalTets.size(); x++)
  {
    // color them randomly
    Real red = (Real)rand() / RAND_MAX;
    Real green = (Real)rand() / RAND_MAX;
    Real blue = (Real)rand() / RAND_MAX;

    glColor4f(red, green, blue, 1.0f);
    _finalTets[x]->drawTriangles();
  }
}

//////////////////////////////////////////////////////////////////////
// Draw all the outside tets as triangles
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::drawOutsideTets()
{
  for (unsigned int x = 0; x < _outsideTets.size(); x++)
    _outsideTets[x]->drawTriangles();
}

//////////////////////////////////////////////////////////////////////
// Draw all the tets as lines
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::drawLines()
{
  for (unsigned int x = 0; x < _tets.size(); x++)
    _tets[x]->drawLines();
}

//////////////////////////////////////////////////////////////////////
// Draw all the surface tets as lines
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::drawSurfaceLines()
{
  for (unsigned int x = 0; x < _surfaceTets.size(); x++)
    _surfaceTets[x]->drawLines();
}

//////////////////////////////////////////////////////////////////////
// Draw all the inside tets as lines
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::drawInsideLines()
{
  for (unsigned int x = 0; x < _insideTets.size(); x++)
    _insideTets[x]->drawLines();
}

//////////////////////////////////////////////////////////////////////
// Draw all the final tets as lines
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::drawFinalLines()
{
  for (unsigned int x = 0; x < _finalTets.size(); x++)
    _finalTets[x]->drawLines();
}

//////////////////////////////////////////////////////////////////////
// Draw all the outside tets as lines
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::drawOutsideLines()
{
  for (unsigned int x = 0; x < _outsideTets.size(); x++)
    _outsideTets[x]->drawLines();
}

//////////////////////////////////////////////////////////////////////
// Draw a point at each grid center
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::drawCenters()
{
  glBegin(GL_POINTS);
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        VEC3F* center = centerVec(x,y,z); 
#ifdef SINGLE_PRECISION      
        glVertex3fv(*center);
#else
        glVertex3dv(*center);
#endif
      }
  glEnd();
}

//////////////////////////////////////////////////////////////////////
// Draw the nodes that got sliced
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::drawSlicedNodes()
{
  glBegin(GL_POINTS);
    for (unsigned int x = 0; x < _sliced.size(); x++)
#ifdef SINGLE_PRECISION      
      glVertex3fv(*_sliced[x]);
#else
      glVertex3dv(*_sliced[x]);
#endif
  glEnd();
}

//////////////////////////////////////////////////////////////////////
// Draw the cut points
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::drawCutPoints()
{
  glBegin(GL_POINTS);
    map<EDGE, VEC3F*>::iterator i;
    for (i = _cutPoints.begin(); i != _cutPoints.end(); ++i)
#ifdef SINGLE_PRECISION      
      glVertex3fv(*(i->second));
#else
      glVertex3dv(*(i->second));
#endif
  glEnd();
}

//////////////////////////////////////////////////////////////////////
// Draw the non-face-shared points
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::drawNonSharedFaceNodes()
{
  glBegin(GL_POINTS);
    for (unsigned int x = 0; x < _nonSharedFaceVertices.size(); x++)
#ifdef SINGLE_PRECISION      
      glVertex3fv((*_nonSharedFaceVertices[x]));
#else
      glVertex3dv((*_nonSharedFaceVertices[x]));
#endif
  glEnd();
}

//////////////////////////////////////////////////////////////////////
// Draw the unconstrained points
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::drawUnconstrainedNodes()
{
  glBegin(GL_POINTS);
    for (unsigned int x = 0; x < _unconstrainedNodes.size(); x++)
#ifdef SINGLE_PRECISION      
      glVertex3fv((*_unconstrainedNodes[x]));
#else
      glVertex3dv((*_unconstrainedNodes[x]));
#endif
  glEnd();
}

//////////////////////////////////////////////////////////////////////
// Draw the constrained points
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::drawConstrainedNodes()
{
  glBegin(GL_POINTS);
    for (unsigned int x = 0; x < _constrainedNodes.size(); x++)
#ifdef SINGLE_PRECISION      
      glVertex3fv((*_constrainedNodes[x]));
#else
      glVertex3dv((*_constrainedNodes[x]));
#endif
  glEnd();
}

//////////////////////////////////////////////////////////////////////
// get VEC3 corresponding to the center of this grid cell
//////////////////////////////////////////////////////////////////////
VEC3F* ISO_STUFFER::centerVec(int x, int y, int z)
{
  //int index = (2 * x + 1) + (2 * y + 1) * _xRes + (2 * z + 1) * _slabSize;
  int index = (2 * x + 1) + 
              (2 * y + 1) * 2 * _xRes + 
              (2 * z + 1) * 4 * _slabSize;

  // first see if this vertex was created before
  map<int, VEC3F*>::iterator i = _vertices.find(index);
  if (i != _vertices.end())
    return (i->second);

  VEC3F* center = new VEC3F();
  (*center)[0] = (Real)(x + 0.5f) / _xRes;
  (*center)[1] = (Real)(y + 0.5f) / _yRes;
  (*center)[2] = (Real)(z + 0.5f) / _zRes;
  _vertices[index] = center;
  _vertexIDs[_vertices[index]] = index;

  return _vertices[index];
}

//////////////////////////////////////////////////////////////////////
// get VEC3 corresponding to a corner of this cell
//////////////////////////////////////////////////////////////////////
bool ISO_STUFFER::cornerVecExists(int x, int y, int z, 
                                  int cornerX, int cornerY, int cornerZ)
{
  int index = (2 * (x + cornerX)) + 
              (2 * (y + cornerY)) * 2 * _xRes + 
              (2 * (z + cornerZ)) * 4 * _slabSize;

  // first see if this vertex was created before
  map<int, VEC3F*>::iterator i = _vertices.find(index);
  if (i != _vertices.end())
    return true;

  return false;
}

//////////////////////////////////////////////////////////////////////
// get VEC3 corresponding to a corner of this cell
//////////////////////////////////////////////////////////////////////
VEC3F* ISO_STUFFER::cornerVec(int x, int y, int z, 
                             int cornerX, int cornerY, int cornerZ)
{
  int index = (2 * (x + cornerX)) + 
              (2 * (y + cornerY)) * 2 * _xRes + 
              (2 * (z + cornerZ)) * 4 * _slabSize;

  // first see if this vertex was created before
  map<int, VEC3F*>::iterator i = _vertices.find(index);
  if (i != _vertices.end())
    return i->second;

  VEC3F* corner = new VEC3F();
  (*corner)[0] = (Real)(x + cornerX) / _xRes;
  (*corner)[1] = (Real)(y + cornerY) / _yRes;
  (*corner)[2] = (Real)(z + cornerZ) / _zRes;
  _vertices[index] = corner;
  _vertexIDs[_vertices[index]] = index;

  return _vertices[index];
}

//////////////////////////////////////////////////////////////////////
// add a vector to the list, unguarded
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::addVecUnguarded(VEC3F* vec)
{
  // stack them all at the end of map
  static int index = 32 * _xRes * _yRes * _zRes;
  index++;

  map<int, VEC3F*>::iterator i = _vertices.find(index);
  assert(i == _vertices.end());

  _vertices[index] = vec;
  _vertexIDs[_vertices[index]] = index;
}

//////////////////////////////////////////////////////////////////////
// get VEC3 corresponding to a corner of this cell
//////////////////////////////////////////////////////////////////////
VEC3F* ISO_STUFFER::faceCenterVecUnguarded(int x, int y, int z, 
                                  int faceX, int faceY, int faceZ)
{
  // stack them all at the end of map
  static int index = 8 * _xRes * _yRes * _zRes;
  index++;

  VEC3F* face = new VEC3F();
  (*face)[0] = (Real)(x + 0.5 + 0.5 * faceX) / _xRes;
  (*face)[1] = (Real)(y + 0.5 + 0.5 * faceY) / _yRes;
  (*face)[2] = (Real)(z + 0.5 + 0.5 * faceZ) / _zRes;

  map<int, VEC3F*>::iterator i = _vertices.find(index);
  assert(i == _vertices.end());

  _vertices[index] = face;
  _vertexIDs[_vertices[index]] = index;

  return _vertices[index];
}

//////////////////////////////////////////////////////////////////////
// get VEC3 corresponding to a corner of this cell
//////////////////////////////////////////////////////////////////////
VEC3F* ISO_STUFFER::faceCenterVec(int x, int y, int z, 
                                  int faceX, int faceY, int faceZ)
{
  int index = (2 * x + 1 + faceX) + 
              (2 * y + 1 + faceY) * 2 * _xRes + 
              (2 * z + 1 + faceZ) * 4 * _slabSize;

  // first see if this vertex was created before
  map<int, VEC3F*>::iterator i = _vertices.find(index);
  if (i != _vertices.end())
    return i->second;

  VEC3F* face = new VEC3F();
  (*face)[0] = (Real)(x + 0.5 + 0.5 * faceX) / _xRes;
  (*face)[1] = (Real)(y + 0.5 + 0.5 * faceY) / _yRes;
  (*face)[2] = (Real)(z + 0.5 + 0.5 * faceZ) / _zRes;
  _vertices[index] = face;
  _vertexIDs[_vertices[index]] = index;

  return _vertices[index];
}

//////////////////////////////////////////////////////////////////////
// Record the connectivity information of the passed in tet
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::recordConnectivity(TET* tet)
{
  // for each vertex
  for (int x = 0; x < 4; x++)
  {
    // retrieve its connectivity vector
    vector<VEC3F*>& edges = _connectivity[tet->vertices[x]];

    // insert all the incident vertices
    for (int y = x + 1; y % 4 != x; y++)
    {
      int index = y % 4;
      VEC3F* incidentVertex = tet->vertices[index];

      // make sure not to make any duplicates
      // (is there an STL type that will do this for you?)
      bool found = false;
      for (unsigned int z = 0; z < edges.size(); z++)
        if (edges[z] == incidentVertex) found = true;

      if (!found)
        edges.push_back(incidentVertex);
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Record the edge colors of the passed in tet
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::recordColors(TET* tet)
{
  recordColors(tet->vertices[0], tet->vertices[1]);
  recordColors(tet->vertices[0], tet->vertices[2]);
  recordColors(tet->vertices[0], tet->vertices[3]);
  recordColors(tet->vertices[1], tet->vertices[2]);
  recordColors(tet->vertices[1], tet->vertices[3]);
  recordColors(tet->vertices[2], tet->vertices[3]);
}

//////////////////////////////////////////////////////////////////////
// Record the edge colors of the passed in edge
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::recordColors(VEC3F* v0, VEC3F* v1)
{
  EDGE edge = createEdge(v0, v1);
  CENTER_CORNER_TYPE type0 = _centerCornerType[v0];
  CENTER_CORNER_TYPE type1 = _centerCornerType[v1];
  if (type0 == type1)
    _edgeTypes[edge] = BLACK;
  else
    _edgeTypes[edge] = RED;
}

//////////////////////////////////////////////////////////////////////
// Generate all possible tets in the grid (for debugging)
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateAllTets()
{
  cout << "Generating exhaustive tets ... "; flush(cout);
  // stomp all old tets;
  _tets.clear();

  // stomp old connectivity info
  _connectivity.clear();

  for (int z = 1; z < _zRes - 1; z++)
  {
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++)
      {
        VEC3F* vert000 = corner000(x,y,z);
        VEC3F* vert100 = corner100(x,y,z);
        VEC3F* vert010 = corner010(x,y,z);
        VEC3F* vert110 = corner110(x,y,z);
        VEC3F* vert001 = corner001(x,y,z);
        VEC3F* vert101 = corner101(x,y,z);
        VEC3F* vert011 = corner011(x,y,z);

        VEC3F* center = centerVec(x,y,z);
        VEC3F* left   = centerVec(x-1,y,z);
        VEC3F* below  = centerVec(x,y-1,z);
        VEC3F* neary   = centerVec(x,y,z-1);

        _centerCornerType[vert000] = CORNER;
        _centerCornerType[vert100] = CORNER;
        _centerCornerType[vert010] = CORNER;
        _centerCornerType[vert110] = CORNER;
        _centerCornerType[vert001] = CORNER;
        _centerCornerType[vert101] = CORNER;
        _centerCornerType[vert011] = CORNER;
        _centerCornerType[center] = CENTER;
        _centerCornerType[left]   = CENTER;
        _centerCornerType[below]  = CENTER;
        _centerCornerType[neary]   = CENTER;

        // type 0 tet
        // the first from the left in Fig. 2 of Isosurface Stuffing
        TET* type0 = new TET(center, left, vert011, vert010);
        _tets.push_back(type0);
        recordConnectivity(_tets.back());
        recordColors(_tets.back());

        // type 0 reflected about y
        TET* type0Reflected = new TET(center, left, vert000, vert001);
        _tets.push_back(type0Reflected);
        recordConnectivity(_tets.back());
        recordColors(_tets.back());

        // type 1 tet
        // the second from the left in Fig. 2 of Isosurface Stuffing
        TET* type1 = new TET(center, left, vert001, vert011);
        _tets.push_back(type1);
        recordConnectivity(_tets.back());
        recordColors(_tets.back());

        // type 1 reflected about y
        TET* type1Reflected = new TET(center, left, vert010, vert000);
        _tets.push_back(type1Reflected);
        recordConnectivity(_tets.back());
        recordColors(_tets.back());

        // type 2 tet
        // rightmost in Fig. 2 of Isosurface Stuffing
        TET* type2 = new TET(center, neary, vert000, vert010);
        _tets.push_back(type2);
        recordConnectivity(_tets.back());
        recordColors(_tets.back());

        // type 2 tet reflecte about x
        // rightmost in Fig. 2 of Isosurface Stuffing
        TET* type2Reflected = new TET(center, neary, vert110, vert100);
        _tets.push_back(type2Reflected);
        recordConnectivity(_tets.back());
        recordColors(_tets.back());

        // type 3 tet
        TET* type3 = new TET(center, below, vert001, vert000);
        _tets.push_back(type3);
        recordConnectivity(_tets.back());
        recordColors(_tets.back());

        // type 3 tet reflected
        TET* type3Reflected = new TET(center, below, vert100, vert101);
        _tets.push_back(type3Reflected);
        recordConnectivity(_tets.back());
        recordColors(_tets.back());

        // type 4 tet
        TET* type4 = new TET(center, below, vert000, vert100);
        _tets.push_back(type4);
        recordConnectivity(_tets.back());
        recordColors(_tets.back());

        // type 4 tet reflected
        TET* type4Reflected = new TET(center, below, vert101, vert001);
        _tets.push_back(type4Reflected);
        recordConnectivity(_tets.back());
        recordColors(_tets.back());

        // type 5 tet
        TET* type5 = new TET(center, neary, vert010, vert110);
        _tets.push_back(type5);
        recordConnectivity(_tets.back());
        recordColors(_tets.back());

        // type 5 tet reflected
        TET* type5Reflected = new TET(center, neary, vert100, vert000);
        _tets.push_back(type5Reflected);
        recordConnectivity(_tets.back());
        recordColors(_tets.back());
      }

    // output status every 10% (making sure not to mod by 0)
    if (_zRes / 10 > 0)
    {
      if (z % (int)(_zRes / 10) == 0)
      {
        cout << 100 * ((Real)z / _zRes) << "% ";
        flush(cout);
      }
    }
  }
  cout << "done." << endl;
}

//////////////////////////////////////////////////////////////////////
// Check if the tet corresponding to v0,v1,v2,v3 could potentially
// be in the final mesh. If it cannot (ie it is entirely outside)
// then do not create a tet for it.
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::emitLimitedTet(VEC3F* v0, VEC3F* v1, VEC3F* v2, VEC3F* v3, 
                                 bool& v0Used, bool& v1Used, bool& v2Used, bool& v3Used)
{
  bool inside = _inside[_vertexIDs[v0]] ||
                _inside[_vertexIDs[v1]] ||
                _inside[_vertexIDs[v2]] ||
                _inside[_vertexIDs[v3]];

  if (inside)
  {
    TET* tet = new TET(v0, v1, v2, v3);
    _tets.push_back(tet);
    recordConnectivity(_tets.back());
    recordColors(_tets.back());
    v0Used = true;
    v1Used = true;
    v2Used = true;
    v3Used = true;
  }
}

//////////////////////////////////////////////////////////////////////
// Delete all the auxiliary information if this vertex was not used
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::cleanupVertex(bool used, VEC3F* vertex)
{
  if (!used)
  {
    int index = _vertexIDs[vertex];
    _vertexIDs.erase(vertex);
    delete _vertices[index];
    _vertices.erase(index);
    _distances.erase(index);
    _inside.erase(index);
  }
  else
    _vertexUsed[vertex] = true;
}

//////////////////////////////////////////////////////////////////////
// Generate a limited set of tets on the grid -- only those that are
// entirely inside the surface, or could potentially intersect the
// surface. Tets entirely outside are not instantiated.
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateLimitedTets(SURFACE& surface)
{

  cout << "Generating limited tets ... "; flush(cout);
  // stomp all old tets;
  _tets.clear();

  // stomp old connectivity info
  _connectivity.clear();

  for (int z = 1; z < _zRes - 1; z++)
  {
    TIMER total;
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++)
      {
        VEC3F* vert000 = corner000(x,y,z);
        VEC3F* vert100 = corner100(x,y,z);
        VEC3F* vert010 = corner010(x,y,z);
        VEC3F* vert110 = corner110(x,y,z);
        VEC3F* vert001 = corner001(x,y,z);
        VEC3F* vert101 = corner101(x,y,z);
        VEC3F* vert011 = corner011(x,y,z);

        VEC3F* center = centerVec(x,y,z);
        VEC3F* left   = centerVec(x-1,y,z);
        VEC3F* below  = centerVec(x,y-1,z);
        VEC3F* neary   = centerVec(x,y,z-1);

        _centerCornerType[vert000] = CORNER;
        _centerCornerType[vert100] = CORNER;
        _centerCornerType[vert010] = CORNER;
        _centerCornerType[vert110] = CORNER;
        _centerCornerType[vert001] = CORNER;
        _centerCornerType[vert101] = CORNER;
        _centerCornerType[vert011] = CORNER;
        _centerCornerType[center] = CENTER;
        _centerCornerType[left]   = CENTER;
        _centerCornerType[below]  = CENTER;
        _centerCornerType[neary]   = CENTER;

        TIMER distance;
        _distances[_vertexIDs[vert000]] = surface.distance(*vert000);
        _distances[_vertexIDs[vert100]] = surface.distance(*vert100);
        _distances[_vertexIDs[vert010]] = surface.distance(*vert010);
        _distances[_vertexIDs[vert110]] = surface.distance(*vert110);
        _distances[_vertexIDs[vert001]] = surface.distance(*vert001);
        _distances[_vertexIDs[vert101]] = surface.distance(*vert101);
        _distances[_vertexIDs[vert011]] = surface.distance(*vert011);
        _distances[_vertexIDs[center]]  = surface.distance(*center);
        _distances[_vertexIDs[left]]    = surface.distance(*left);
        _distances[_vertexIDs[below]]   = surface.distance(*below);
        _distances[_vertexIDs[neary]]    = surface.distance(*neary);
        _timingBreakdown["Distance Test"] += distance.timing();
     
        TIMER inside;
        if (_inside.find(_vertexIDs[vert000]) == _inside.end())
          _inside[_vertexIDs[vert000]] = surface.inside(*vert000);
        if (_inside.find(_vertexIDs[vert100]) == _inside.end())
          _inside[_vertexIDs[vert100]] = surface.inside(*vert100);
        if (_inside.find(_vertexIDs[vert010]) == _inside.end())
          _inside[_vertexIDs[vert010]] = surface.inside(*vert010);
        if (_inside.find(_vertexIDs[vert110]) == _inside.end())
          _inside[_vertexIDs[vert110]] = surface.inside(*vert110);
        if (_inside.find(_vertexIDs[vert001]) == _inside.end())
          _inside[_vertexIDs[vert001]] = surface.inside(*vert001);
        if (_inside.find(_vertexIDs[vert101]) == _inside.end())
          _inside[_vertexIDs[vert101]] = surface.inside(*vert101);
        if (_inside.find(_vertexIDs[vert011]) == _inside.end())
          _inside[_vertexIDs[vert011]] = surface.inside(*vert011);
        if (_inside.find(_vertexIDs[center]) == _inside.end())
          _inside[_vertexIDs[center]]  = surface.inside(*center);
        if (_inside.find(_vertexIDs[left]) == _inside.end())
          _inside[_vertexIDs[left]]    = surface.inside(*left);
        if (_inside.find(_vertexIDs[below]) == _inside.end())
          _inside[_vertexIDs[below]]   = surface.inside(*below);
        if (_inside.find(_vertexIDs[neary]) == _inside.end())
          _inside[_vertexIDs[neary]]    = surface.inside(*neary);
        _timingBreakdown["Inside Test"] += inside.timing();

        // track if a vertex was used
        bool vert000Used = _vertexUsed[vert000];
        bool vert100Used = _vertexUsed[vert100];
        bool vert010Used = _vertexUsed[vert010];
        bool vert110Used = _vertexUsed[vert110];
        bool vert001Used = _vertexUsed[vert001];
        bool vert101Used = _vertexUsed[vert101];
        bool vert011Used = _vertexUsed[vert011];
        bool centerUsed  = _vertexUsed[center];
        bool leftUsed    = _vertexUsed[left];
        bool belowUsed   = _vertexUsed[below];
        bool nearUsed    = _vertexUsed[neary];

        emitLimitedTet(center, left,  vert011, vert010, centerUsed, leftUsed,  vert011Used, vert010Used);
        emitLimitedTet(center, left,  vert000, vert001, centerUsed, leftUsed,  vert000Used, vert001Used);
        emitLimitedTet(center, left,  vert001, vert011, centerUsed, leftUsed,  vert001Used, vert011Used);
        emitLimitedTet(center, left,  vert010, vert000, centerUsed, leftUsed,  vert010Used, vert000Used);
        emitLimitedTet(center, neary,  vert000, vert010, centerUsed, nearUsed,  vert000Used, vert010Used);
        emitLimitedTet(center, neary,  vert110, vert100, centerUsed, nearUsed,  vert110Used, vert100Used);
        emitLimitedTet(center, below, vert001, vert000, centerUsed, belowUsed, vert001Used, vert000Used);
        emitLimitedTet(center, below, vert100, vert101, centerUsed, belowUsed, vert100Used, vert101Used);
        emitLimitedTet(center, below, vert000, vert100, centerUsed, belowUsed, vert000Used, vert100Used);
        emitLimitedTet(center, below, vert101, vert001, centerUsed, belowUsed, vert101Used, vert001Used);
        emitLimitedTet(center, neary,  vert010, vert110, centerUsed, nearUsed,  vert010Used, vert110Used);
        emitLimitedTet(center, neary,  vert100, vert000, centerUsed, nearUsed,  vert100Used, vert000Used);

        cleanupVertex(vert000Used, vert000);
        cleanupVertex(vert100Used, vert100);
        cleanupVertex(vert010Used, vert010);
        cleanupVertex(vert110Used, vert110);
        cleanupVertex(vert001Used, vert001);
        cleanupVertex(vert101Used, vert101);
        cleanupVertex(vert011Used, vert011);
      }
    _totalTime += total.timing();

    // output status every 10% if we have enough work to do (extra
    // check here to avoid division by zero for small meshes).
    if (_zRes / 10 > 0)
    {
      if (z % (int)(_zRes / 10) == 0)
      {
        cout << 100 * ((Real)z / _zRes) << "% ";
        flush(cout);
      }
    }
  }
  cout << "done." << endl;

}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
bool ISO_STUFFER::isOccupied(Real x, Real y, Real z, SURFACE& surface)
{
  Real point[3];
  point[0] = x * (1.0 / _xRes);
  point[1] = y * (1.0 / _yRes);
  point[2] = z * (1.0 / _zRes);
  return surface.isOccupied(point);
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateCubeTetsWithClonedVertices(OBJ& surface)
{
  map<VEC3F*, bool> xCuts;
  map<VEC3F*, bool> yCuts;
  map<VEC3F*, bool> zCuts;

  cout << "Generating Cube tets with vertex cloning ... "; flush(cout);
  // stomp all old tets;
  _finalTets.clear();

  for (int z = 0; z < _zRes; z++)
  {
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        float xCell = x + 0.5;
        float yCell = y + 0.5;
        float zCell = z + 0.5;

        // Before creating the centerVec, see if we should even bother
        if (!isOccupied(xCell, yCell, zCell, surface)) continue;

        // if this is going to create a single "hinge cube", don't bother
        int neighbors = 0;
        if (isOccupied(xCell + 1, yCell, zCell, surface)) neighbors++;
        if (isOccupied(xCell - 1, yCell, zCell, surface)) neighbors++;
        if (isOccupied(xCell, yCell + 1, zCell, surface)) neighbors++;
        if (isOccupied(xCell, yCell - 1, zCell, surface)) neighbors++;
        if (isOccupied(xCell, yCell, zCell + 1, surface)) neighbors++;
        if (isOccupied(xCell, yCell, zCell - 1, surface)) neighbors++;

        if (neighbors < 1) continue;

        VEC3F* center = centerVec(x,y,z);

        VEC3F* vert000 = corner000(x,y,z);
        VEC3F* vert100 = corner100(x,y,z);
        VEC3F* vert010 = corner010(x,y,z);
        VEC3F* vert110 = corner110(x,y,z);
        VEC3F* vert001 = corner001(x,y,z);
        VEC3F* vert101 = corner101(x,y,z);
        VEC3F* vert011 = corner011(x,y,z);
        VEC3F* vert111 = corner111(x,y,z);

        VEC3F* left  = NULL;
        if (isOccupied(xCell - 1, yCell, zCell, surface))
        {
          if (isOccupied(xCell - 0.5, yCell, zCell, surface))
            left = centerVec(x-1,y,z);
          else
          {
            left = faceCenterVecUnguarded(x,y,z,-1,0,0);
            // record that these vertices later need to be cut
            xCuts[vert000] = true;
            xCuts[vert001] = true;
            xCuts[vert010] = true;
            xCuts[vert011] = true;
          }
        }
        else
          left = faceCenterVec(x,y,z,-1,0,0);

        VEC3F* below = NULL;
        if (isOccupied(xCell, yCell - 1, zCell, surface))
        {
          if (isOccupied(xCell, yCell - 0.5, zCell, surface))
            below = centerVec(x,y-1,z);
          else
          {
            below = faceCenterVecUnguarded(x,y,z,0,-1,0);
            yCuts[vert000] = true;
            yCuts[vert001] = true;
            yCuts[vert100] = true;
            yCuts[vert101] = true;
          }
        }
        else
          below = faceCenterVec(x,y,z,0,-1,0);

        VEC3F* near = NULL;
        if  (isOccupied(xCell, yCell, zCell - 1, surface))
        {
          if (isOccupied(xCell, yCell, zCell - 0.5, surface))
            near = centerVec(x,y,z-1);
          else
          { 
            near = faceCenterVecUnguarded(x,y,z,0,0,-1);
            zCuts[vert000] = true;
            zCuts[vert010] = true;
            zCuts[vert100] = true;
            zCuts[vert110] = true;
          }
        }
        else
          near = faceCenterVec(x,y,z,0,0,-1);

        // do all the ones to the left
        TET* newTet0 = new TET(center, vert010, left, vert011);
        _finalTets.push_back(newTet0);
        TET* newTet1 = new TET(center, vert000, left, vert010);
        _finalTets.push_back(newTet1);
        TET* newTet2 = new TET(center, vert000, vert001, left);
        _finalTets.push_back(newTet2);
        TET* newTet3 = new TET(center, left, vert001, vert011);
        _finalTets.push_back(newTet3);

        // do all the ones that are near
        TET* newTet4 = new TET(center, vert110, near, vert010);
        _finalTets.push_back(newTet4);
        TET* newTet5 = new TET(center, vert100, near, vert110);
        _finalTets.push_back(newTet5);
        TET* newTet6 = new TET(center, vert100, vert000, near);
        _finalTets.push_back(newTet6);
        TET* newTet7 = new TET(center, near, vert000, vert010);
        _finalTets.push_back(newTet7);

        // do all the ones that are below
        TET* newTet8 = new TET(center, vert100, below, vert000);
        _finalTets.push_back(newTet8);
        TET* newTet9 = new TET(center, vert101, below, vert100);
        _finalTets.push_back(newTet9);
        TET* newTet10 = new TET(center, vert101, vert001, below);
        _finalTets.push_back(newTet10);
        TET* newTet11 = new TET(center, below, vert001, vert000);
        _finalTets.push_back(newTet11);

        // check the top, right, and far ones to see if they
        // won't be creating the other sides of the cube
        if (!isOccupied(xCell + 1, yCell, zCell,surface))
        {
          VEC3F* right = faceCenterVec(x,y,z,1,0,0);
          TET* newTet12 = new TET(right, vert110, center, vert111);
          TET* newTet13 = new TET(right, vert100, center, vert110);
          TET* newTet14 = new TET(right, vert100, vert101, center);
          TET* newTet15 = new TET(right, vert101, vert111, center);
          _finalTets.push_back(newTet12);
          _finalTets.push_back(newTet13);
          _finalTets.push_back(newTet14);
          _finalTets.push_back(newTet15);
        }
        else if (!isOccupied(xCell + 0.5, yCell, zCell,surface))
        {
          VEC3F* right = faceCenterVecUnguarded(x,y,z,1,0,0);
          TET* newTet12 = new TET(right, vert110, center, vert111);
          TET* newTet13 = new TET(right, vert100, center, vert110);
          TET* newTet14 = new TET(right, vert100, vert101, center);
          TET* newTet15 = new TET(right, vert101, vert111, center);
          _finalTets.push_back(newTet12);
          _finalTets.push_back(newTet13);
          _finalTets.push_back(newTet14);
          _finalTets.push_back(newTet15);
          xCuts[vert100] = true;
          xCuts[vert101] = true;
          xCuts[vert110] = true;
          xCuts[vert111] = true;
        }

        if (!isOccupied(xCell, yCell + 1, zCell,surface))
        {
          VEC3F* up = faceCenterVec(x,y,z,0,1,0);
          TET* newTet12 = new TET(up, vert010, center, vert011);
          TET* newTet13 = new TET(up, vert110, center, vert010);
          TET* newTet14 = new TET(up, vert110, vert111, center);
          TET* newTet15 = new TET(up, vert111, vert011, center);
          _finalTets.push_back(newTet12);
          _finalTets.push_back(newTet13);
          _finalTets.push_back(newTet14);
          _finalTets.push_back(newTet15);
        }
        else if (!isOccupied(xCell, yCell + 0.5, zCell,surface))
        {
          VEC3F* up = faceCenterVecUnguarded(x,y,z,0,1,0);
          TET* newTet12 = new TET(up, vert010, center, vert011);
          TET* newTet13 = new TET(up, vert110, center, vert010);
          TET* newTet14 = new TET(up, vert110, vert111, center);
          TET* newTet15 = new TET(up, vert111, vert011, center);
          _finalTets.push_back(newTet12);
          _finalTets.push_back(newTet13);
          _finalTets.push_back(newTet14);
          _finalTets.push_back(newTet15);
          yCuts[vert010] = true;
          yCuts[vert011] = true;
          yCuts[vert110] = true;
          yCuts[vert111] = true;
        }

        if (!isOccupied(xCell, yCell, zCell + 1,surface))
        {
          VEC3F* far = faceCenterVec(x,y,z,0,0,1);
          TET* newTet12 = new TET(far, vert111, center, vert011);
          TET* newTet13 = new TET(far, vert101, center, vert111);
          TET* newTet14 = new TET(far, vert101, vert001, center);
          TET* newTet15 = new TET(far, vert001, vert011, center);
          _finalTets.push_back(newTet12);
          _finalTets.push_back(newTet13);
          _finalTets.push_back(newTet14);
          _finalTets.push_back(newTet15);
        }
        else if (!isOccupied(xCell, yCell, zCell + 0.5,surface))
        {
          VEC3F* far = faceCenterVecUnguarded(x,y,z,0,0,1);
          TET* newTet12 = new TET(far, vert111, center, vert011);
          TET* newTet13 = new TET(far, vert101, center, vert111);
          TET* newTet14 = new TET(far, vert101, vert001, center);
          TET* newTet15 = new TET(far, vert001, vert011, center);
          _finalTets.push_back(newTet12);
          _finalTets.push_back(newTet13);
          _finalTets.push_back(newTet14);
          _finalTets.push_back(newTet15);
          zCuts[vert001] = true;
          zCuts[vert011] = true;
          zCuts[vert101] = true;
          zCuts[vert111] = true;
        }
      }

    // output status every 10% (making sure not to mod by 0)
    if (_zRes / 10 > 0)
    {
      if (z % (int)(_zRes / 10) == 0)
      {
        cout << 100 * ((Real)z / _zRes) << "% ";
        flush(cout);
      }
    }
  }
  cout << "done." << endl;

  cout << " Slicing vertices ... "; flush(cout);
  int totalSliced = sliceVertices(xCuts, yCuts, zCuts);
  cout << "sliced " << totalSliced << " vertices." << endl;

  for (unsigned int x = 0; x < _finalTets.size(); x++)
    if (_finalTets[x]->invalid())
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " Generated an invalid tet ordering!!!! " << endl;
    }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateCubeTets(SURFACE& surface)
{
  cout << "Generating Cube tets ... "; flush(cout);
  // stomp all old tets;
  _finalTets.clear();

  for (int z = 0; z < _zRes; z++)
  {
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        float xCell = x + 0.5;
        float yCell = y + 0.5;
        float zCell = z + 0.5;

        // Before creating the centerVec, see if we should even bother
        if (!isOccupied(xCell, yCell, zCell, surface)) continue;

        // if this is going to create a single "hinge cube", don't bother
        int neighbors = 0;
        if (isOccupied(xCell + 1, yCell, zCell, surface)) neighbors++;
        if (isOccupied(xCell - 1, yCell, zCell, surface)) neighbors++;
        if (isOccupied(xCell, yCell + 1, zCell, surface)) neighbors++;
        if (isOccupied(xCell, yCell - 1, zCell, surface)) neighbors++;
        if (isOccupied(xCell, yCell, zCell + 1, surface)) neighbors++;
        if (isOccupied(xCell, yCell, zCell - 1, surface)) neighbors++;

        if (neighbors < 1) continue;

        VEC3F* center = centerVec(x,y,z);

        VEC3F* vert000 = corner000(x,y,z);
        VEC3F* vert100 = corner100(x,y,z);
        VEC3F* vert010 = corner010(x,y,z);
        VEC3F* vert110 = corner110(x,y,z);
        VEC3F* vert001 = corner001(x,y,z);
        VEC3F* vert101 = corner101(x,y,z);
        VEC3F* vert011 = corner011(x,y,z);
        VEC3F* vert111 = corner111(x,y,z);

        VEC3F* left  = isOccupied(xCell - 1, yCell, zCell, surface) ? centerVec(x-1,y,z) : faceCenterVec(x,y,z,-1,0,0);
        VEC3F* below = isOccupied(xCell, yCell - 1, zCell, surface) ? centerVec(x,y-1,z) : faceCenterVec(x,y,z,0,-1,0);
        VEC3F* near  = isOccupied(xCell, yCell, zCell - 1, surface) ? centerVec(x,y,z-1) : faceCenterVec(x,y,z,0,0,-1);
      
        // do all the ones to the left
        TET* newTet0 = new TET(center, vert010, left, vert011);
        _finalTets.push_back(newTet0);
        TET* newTet1 = new TET(center, vert000, left, vert010);
        _finalTets.push_back(newTet1);
        TET* newTet2 = new TET(center, vert000, vert001, left);
        _finalTets.push_back(newTet2);
        TET* newTet3 = new TET(center, left, vert001, vert011);
        _finalTets.push_back(newTet3);

        // do all the ones that are near
        TET* newTet4 = new TET(center, vert110, near, vert010);
        _finalTets.push_back(newTet4);
        TET* newTet5 = new TET(center, vert100, near, vert110);
        _finalTets.push_back(newTet5);
        TET* newTet6 = new TET(center, vert100, vert000, near);
        _finalTets.push_back(newTet6);
        TET* newTet7 = new TET(center, near, vert000, vert010);
        _finalTets.push_back(newTet7);

        // do all the ones that are below
        TET* newTet8 = new TET(center, vert100, below, vert000);
        _finalTets.push_back(newTet8);
        TET* newTet9 = new TET(center, vert101, below, vert100);
        _finalTets.push_back(newTet9);
        TET* newTet10 = new TET(center, vert101, vert001, below);
        _finalTets.push_back(newTet10);
        TET* newTet11 = new TET(center, below, vert001, vert000);
        _finalTets.push_back(newTet11);

        // check the top, right, and far ones to see if they
        // won't be creating the other sides of the cube
        if (!isOccupied(xCell + 1, yCell, zCell,surface))
        {
          VEC3F* right = faceCenterVec(x,y,z,1,0,0);
          TET* newTet12 = new TET(right, vert110, center, vert111);
          TET* newTet13 = new TET(right, vert100, center, vert110);
          TET* newTet14 = new TET(right, vert100, vert101, center);
          TET* newTet15 = new TET(right, vert101, vert111, center);
          _finalTets.push_back(newTet12);
          _finalTets.push_back(newTet13);
          _finalTets.push_back(newTet14);
          _finalTets.push_back(newTet15);
        }
        if (!isOccupied(xCell, yCell + 1, zCell,surface))
        {
          VEC3F* up = faceCenterVec(x,y,z,0,1,0);
          TET* newTet12 = new TET(up, vert010, center, vert011);
          TET* newTet13 = new TET(up, vert110, center, vert010);
          TET* newTet14 = new TET(up, vert110, vert111, center);
          TET* newTet15 = new TET(up, vert111, vert011, center);
          _finalTets.push_back(newTet12);
          _finalTets.push_back(newTet13);
          _finalTets.push_back(newTet14);
          _finalTets.push_back(newTet15);
        }
        if (!isOccupied(xCell, yCell, zCell + 1,surface))
        {
          VEC3F* far = faceCenterVec(x,y,z,0,0,1);
          TET* newTet12 = new TET(far, vert111, center, vert011);
          TET* newTet13 = new TET(far, vert101, center, vert111);
          TET* newTet14 = new TET(far, vert101, vert001, center);
          TET* newTet15 = new TET(far, vert001, vert011, center);
          _finalTets.push_back(newTet12);
          _finalTets.push_back(newTet13);
          _finalTets.push_back(newTet14);
          _finalTets.push_back(newTet15);
        }
      }

    // output status every 10% (making sure not to mod by 0)
    if (_zRes / 10 > 0)
    {
      if (z % (int)(_zRes / 10) == 0)
      {
        cout << 100 * ((Real)z / _zRes) << "% ";
        flush(cout);
      }
    }
  }
  cout << "done." << endl;
  cout << " Generated " << _finalTets.size() << " tets " << endl;

  for (unsigned int x = 0; x < _finalTets.size(); x++)
    if (_finalTets[x]->invalid())
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " Generated an invalid tet ordering!!!! " << endl;
    }
}


//////////////////////////////////////////////////////////////////////
// Partition the tets into inside and outside
// 
// Use an arbitrary implicit function, not just an triangle mesh
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateInsideTets(SURFACE& surface)
{
  cout << "Generating inside tets ... "; flush(cout);
  map<int, VEC3F*>::iterator i;

  // generate the inverse mapping from vertex address to 
  // grid index -- this is needed to index into _distances in the next
  // loop
  _vertexIDs.clear();
  for(i = _vertices.begin(); i != _vertices.end(); ++i)
    _vertexIDs[i->second] = i->first;

  _inside.clear();
  string cacheName;
  char buffer[256];
  sprintf(buffer, "%i", _xRes);
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    // do the inside/outside test for each vertex
    for (int y = 0; y < 4; y++)
    {
      // check first if the test has already been done for this vertex
      // (this does appear to speed it up quite a bit)
      int vertexID = _vertexIDs[_tets[x]->vertices[y]];
      if (_inside.find(vertexID) != _inside.end())
        continue;
   
      // do the test 
      VEC3F vertex = *(_tets[x]->vertices[y]);
      _inside[vertexID] = surface.inside(vertex); 
    }
  }

  // flood the rest of the vertices with the inside/outside info
  floodInsideOutside();

  // generate the inside/outside tet lists
  _insideTets.clear();
  _outsideTets.clear();
  _surfaceTets.clear();
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    int insideHits = 0;

    // do the inside/outside test for each vertex
    for (int y = 0; y < 4; y++)
    {
      int vertexID = _vertexIDs[_tets[x]->vertices[y]];
      if (_inside[vertexID])
      {
        if (insideHits == 0)
          _insideTets.push_back(_tets[x]);
        insideHits++;
      }
    }

    // if none of the vertices were inside, add the tet to the outside list
    if (insideHits == 0)
      _outsideTets.push_back(_tets[x]);

    // if there is a sign flip across the vertices, add to the surface list
    if (insideHits > 0 && insideHits < 4)
      _surfaceTets.push_back(_tets[x]);
  }
  cout << "done." << endl;

  cout << "  Inside tets: " << _insideTets.size() << endl;
  cout << "  Outside tets: " << _outsideTets.size() << endl;
  cout << "  Surface tets: " << _surfaceTets.size() << endl;

  cout << "  Generating surface distances " << endl;
  generateSurfaceDistances(surface);

  cout << "  Generating surface edges" << endl;
  // compute a list of all the edges along the surface
  generateSurfaceEdges();

  cout << "  Generating cut points" << endl;
  // generate cut points
  generateCutPoints();

  cout << "  Snapping to cut points " << endl;
  snapToCutPoints();

  cout << "  Generating final tets " << endl;
  generateFinalTets();
}

//////////////////////////////////////////////////////////////////////
// propagate sign information to all the vertices
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::floodInsideOutside()
{
  // push all the vertices with valid 
  // inside/outside values onto the stack
  map<int,bool>::iterator i;
  vector<int> floodFront; 
  for(i = _inside.begin(); i != _inside.end(); ++i)
    floodFront.push_back(i->first);

  // for each item on the stack
  while (floodFront.size() != 0)
  {
    // get the index of the next vertex to flood from
    int vertexID = floodFront.back();
    floodFront.pop_back();

    // get the actual vertex
    VEC3F* vertex = _vertices[vertexID];

    // get its neighbors
    vector<VEC3F*> neighbors = _connectivity[vertex];

    // check each neighbor
    for (unsigned int x = 0; x < neighbors.size(); x++)
    {
      // get the neighbor ID
      int neighborID = _vertexIDs[neighbors[x]];

      // if this neighbor has never been flooded before
      map<int, bool>::iterator visited = _inside.find(neighborID);
      if (visited == _inside.end())
      {
        // assign it the current sign
        _inside[neighborID] = _inside[vertexID];

        // look at this node's neighbors
        vector<VEC3F*> neighborNeighbors = _connectivity[_vertices[neighborID]];

        for (unsigned int y = 0; y < neighborNeighbors.size(); y++)
        {
          // if it has never been looked at, push to the back
          // of the floodFront.
          int neighborNeighborID = _vertexIDs[neighborNeighbors[x]];
          map<int, bool>::iterator visited = _inside.find(neighborNeighborID);
          if (visited == _inside.end())
            floodFront.push_back(neighborNeighborID);
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// generate distance values for nodes that are at the surface,
// ie experience a sign flip
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateSurfaceDistances(SURFACE& surface)
{
  // clear any previous distance values
  _distances.clear();

  // for each tet on the surface
  for (unsigned int x = 0; x < _surfaceTets.size(); x++)
  {
    // for each vertex in the tet
    for (int y = 0; y < 4; y++)
    {
      VEC3F* vertex = _surfaceTets[x]->vertices[y];

      // see if we computed a distance for this one already
      int vertexID = _vertexIDs[vertex];
      map<int,Real>::iterator finder = _distances.find(vertexID);
      if (finder != _distances.end())
        continue;

      // compute the distance
      _distances[vertexID] = surface.distance(*vertex);
    }
  }
}

//////////////////////////////////////////////////////////////////////
// generate the cut points where the edges intersect the isosurface
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateCutPoints()
{
  // clear any previous cut points
  _cutPoints.clear();

  int totalCuts = 0;

  // for each edge
  map<EDGE,EDGE>::iterator i;
  for (i = _surfaceEdges.begin(); i != _surfaceEdges.end(); ++i)
  {
    // get the vertex IDs
    int v0 = _vertexIDs[i->second.first];
    int v1 = _vertexIDs[i->second.second];

    // check for a sign flip
    Real distance0 = _distances[v0];
    Real distance1 = _distances[v1];

    Real sign0 = (_inside[v0]) ? 1.0f : -1.0f;
    Real sign1 = (_inside[v1]) ? 1.0f : -1.0f;

    // if the signs are the same, skip this edge
    if (sign0 * sign1 > 0.0f) continue;
   
    distance0 *= sign0;
    distance1 *= sign1;

    // get the intersection point
    Real cutDistance = -distance0 / (distance1 - distance0);
    if (fabs(distance1 - distance0) < 1e-7)
      cutDistance = 0.0f;
    VEC3F edge = (*_vertices[v1]) - (*_vertices[v0]);
    Real magnitude = norm(edge);
    unitize(edge);
    VEC3F* cutPoint = new VEC3F();
    *cutPoint = (*_vertices[v0]) + cutDistance * edge * magnitude;

    // add it to the list
    _cutPoints[i->first] = cutPoint;

    totalCuts++;
  }
}

//////////////////////////////////////////////////////////////////////
// snap vertices to the cut points
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::snapToCutPoints()
{
  // reset snapped map
  _snapped.clear();

  // for each tet vertex
  map<int, VEC3F*>::iterator i;
  for(i = _vertices.begin(); i != _vertices.end(); ++i)
  {
    VEC3F* vertex = i->second;
    int vertexID = _vertexIDs[vertex];

    // for each adjoining edge
    map<VEC3F*, vector<VEC3F*> >::iterator finder = _connectivity.find(vertex);

    // if there are no neighbors, skip vertex
    if (finder == _connectivity.end()) continue;

    // if there are neighbors, retrieve them
    vector<VEC3F*>& neighbors = finder->second;

    // for each neighboring vertex
    bool snapped = false;
    unsigned int x = 0;
    while (x < neighbors.size() && snapped == false)
    {
      // construct the edge
      EDGE edge = createEdge(vertex, neighbors[x]);

      // get the cut points corresponding to the edge
      map<EDGE, VEC3F*>::iterator cutFinder = _cutPoints.find(edge);

      // if there is a cut, check for a snap
      if (cutFinder != _cutPoints.end()) 
      {
        // get the edge threshold
        EDGE edge = createEdge(vertex, neighbors[x]);
        map<EDGE, EDGE_TYPE>::iterator colorFinder;
        colorFinder = _edgeTypes.find(edge);
        if (colorFinder == _edgeTypes.end())
          cout << __FILE__ << " " << __LINE__ << " NO COLOR FOUND " << endl;

        EDGE_TYPE edgeColor = _edgeTypes[edge];
        Real threshold = (edgeColor == BLACK) ? _blackLength * _alphaLong : _redLength * _alphaShort;

        // calc the distance to the cut point
        VEC3F& cutPoint = *(cutFinder->second);
        Real distance = norm(cutPoint - *vertex);

        // if it is below a threshold, snap it
        if (distance < threshold)
        {
          *vertex = cutPoint;
          _distances[vertexID] = 0.0f;
          _snapped[vertexID] = true;
          snapped = true;

          // delete it so that the vertex isn't double-counted
          _cutPoints.erase(cutFinder);
        }
      }
      x++;
    }

    // if there was a snap, delete the other cut points
    if (snapped)
      for (unsigned int x = 0; x < neighbors.size(); x++)
      {
        // construct the edge
        EDGE edge = createEdge(vertex, neighbors[x]);

        // get the cut points corresponding to the edge
        map<EDGE, VEC3F*>::iterator cutFinder = _cutPoints.find(edge);

        // if there are no cut points, skip this edge
        if (cutFinder != _cutPoints.end()) 
          _cutPoints.erase(cutFinder);
      }
  }
}

//////////////////////////////////////////////////////////////////////
// generate a list of all the edges belonging to tets that intersect 
// the surface
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateSurfaceEdges()
{
  // clear any old surface edges
  _surfaceEdges.clear();

  // for each tet on the surface
  for (unsigned int x = 0; x < _surfaceTets.size(); x++)
  {
    insertEdge(_surfaceTets[x]->vertices[0], _surfaceTets[x]->vertices[1]);
    insertEdge(_surfaceTets[x]->vertices[0], _surfaceTets[x]->vertices[2]);
    insertEdge(_surfaceTets[x]->vertices[0], _surfaceTets[x]->vertices[3]);
    insertEdge(_surfaceTets[x]->vertices[1], _surfaceTets[x]->vertices[2]);
    insertEdge(_surfaceTets[x]->vertices[1], _surfaceTets[x]->vertices[3]);
    insertEdge(_surfaceTets[x]->vertices[2], _surfaceTets[x]->vertices[3]);
  }
}

//////////////////////////////////////////////////////////////////////
// insert an edge into the edge list, enforcing the ordering and
// making sure not to make duplicates (the map takes care of this)
//////////////////////////////////////////////////////////////////////
EDGE ISO_STUFFER::createEdge(VEC3F* v0, VEC3F* v1) {
  VEC3F* first  = v0 < v1 ? v0 : v1;
  VEC3F* second = v0 < v1 ? v1 : v0;
  EDGE edge(first, second);
  return edge;
}

//////////////////////////////////////////////////////////////////////
// insert an edge into the edge list, enforcing the ordering and
// making sure not to make duplicates (the map takes care of this)
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::insertEdge(VEC3F* v0, VEC3F* v1) {
  EDGE edge = createEdge(v0, v1);
  _surfaceEdges[edge] = edge;
}

//////////////////////////////////////////////////////////////////////
// generate the final tets
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateFinalTetsDebug()
{
  // clear previous final tets
  _finalTets.clear();

  // for each inside tet
  for (unsigned int x = 0; x < _insideTets.size(); x++)
  {
    // count number of inside/outside/snapped points
    int insidePoints[] = {-1,-1,-1,-1};
    int outsidePoints[] = {-1,-1,-1,-1};
    int snappedPoints[] = {-1,-1,-1,-1};
    int insideTotal = 0;
    int outsideTotal = 0;
    int snappedTotal = 0;
    for (int y = 0; y < 4; y++)
    {
      int vertexID = _vertexIDs[_insideTets[x]->vertices[y]];
      map<int,bool>::iterator finder = _snapped.find(vertexID);
      if (finder != _snapped.end())
      {
        snappedPoints[y] = vertexID;
        snappedTotal++;
      }
      else if (_inside[vertexID])
      {
        insidePoints[y] = vertexID;
        insideTotal++;
      }
      else
      {
        outsidePoints[y] = vertexID;
        outsideTotal++;
      }
    }

    // count number of edge cuts
    VEC3F* cuts[] = {NULL, NULL, NULL, NULL, NULL, NULL};
    EDGE edges[6];
    edges[0] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[1]);
    edges[1] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[2]);
    edges[2] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[3]);
    edges[3] = createEdge(_insideTets[x]->vertices[1],
                          _insideTets[x]->vertices[2]);
    edges[4] = createEdge(_insideTets[x]->vertices[2],
                          _insideTets[x]->vertices[3]);
    edges[5] = createEdge(_insideTets[x]->vertices[3],
                          _insideTets[x]->vertices[1]);

    // look for a cut along each edge
    int cutsTotal = 0;
    for (int y = 0; y < 6; y++)
    {
      map<EDGE, VEC3F*>::iterator cutFinder;
      cutFinder = _cutPoints.find(edges[y]);
      if (cutFinder != _cutPoints.end())
      {
        cuts[y] = cutFinder->second;
        cutsTotal++;
      }
    }
    
    // case 0-3 in Figure 3 of paper
    if (outsideTotal == 0 && insideTotal > 0)
    {
      _finalTets.push_back(_insideTets[x]);
      continue;
    }
  }
  cout << " After case 0-3 ===============================" << endl;
  collectTetStats();

  // for each inside tet
  for (unsigned int x = 0; x < _insideTets.size(); x++)
  {
    // count number of inside/outside/snapped points
    int insidePoints[] = {-1,-1,-1,-1};
    int outsidePoints[] = {-1,-1,-1,-1};
    int snappedPoints[] = {-1,-1,-1,-1};
    int insideTotal = 0;
    int outsideTotal = 0;
    int snappedTotal = 0;
    for (int y = 0; y < 4; y++)
    {
      int vertexID = _vertexIDs[_insideTets[x]->vertices[y]];
      map<int,bool>::iterator finder = _snapped.find(vertexID);
      if (finder != _snapped.end())
      {
        snappedPoints[y] = vertexID;
        snappedTotal++;
      }
      else if (_inside[vertexID])
      {
        insidePoints[y] = vertexID;
        insideTotal++;
      }
      else
      {
        outsidePoints[y] = vertexID;
        outsideTotal++;
      }
    }

    // count number of edge cuts
    VEC3F* cuts[] = {NULL, NULL, NULL, NULL, NULL, NULL};
    EDGE edges[6];
    edges[0] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[1]);
    edges[1] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[2]);
    edges[2] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[3]);
    edges[3] = createEdge(_insideTets[x]->vertices[1],
                          _insideTets[x]->vertices[2]);
    edges[4] = createEdge(_insideTets[x]->vertices[2],
                          _insideTets[x]->vertices[3]);
    edges[5] = createEdge(_insideTets[x]->vertices[3],
                          _insideTets[x]->vertices[1]);

    // look for a cut along each edge
    int cutsTotal = 0;
    for (int y = 0; y < 6; y++)
    {
      map<EDGE, VEC3F*>::iterator cutFinder;
      cutFinder = _cutPoints.find(edges[y]);
      if (cutFinder != _cutPoints.end())
      {
        cuts[y] = cutFinder->second;
        cutsTotal++;
      }
    }
    // case 4
    if (outsideTotal == 1 && cutsTotal == 1)
    {
      generateCase4(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      continue;
    }
  }
  cout << " After case 4 ===============================" << endl;
  collectTetStats();

  // for each inside tet
  for (unsigned int x = 0; x < _insideTets.size(); x++)
  {
    // count number of inside/outside/snapped points
    int insidePoints[] = {-1,-1,-1,-1};
    int outsidePoints[] = {-1,-1,-1,-1};
    int snappedPoints[] = {-1,-1,-1,-1};
    int insideTotal = 0;
    int outsideTotal = 0;
    int snappedTotal = 0;
    for (int y = 0; y < 4; y++)
    {
      int vertexID = _vertexIDs[_insideTets[x]->vertices[y]];
      map<int,bool>::iterator finder = _snapped.find(vertexID);
      if (finder != _snapped.end())
      {
        snappedPoints[y] = vertexID;
        snappedTotal++;
      }
      else if (_inside[vertexID])
      {
        insidePoints[y] = vertexID;
        insideTotal++;
      }
      else
      {
        outsidePoints[y] = vertexID;
        outsideTotal++;
      }
    }

    // count number of edge cuts
    VEC3F* cuts[] = {NULL, NULL, NULL, NULL, NULL, NULL};
    EDGE edges[6];
    edges[0] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[1]);
    edges[1] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[2]);
    edges[2] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[3]);
    edges[3] = createEdge(_insideTets[x]->vertices[1],
                          _insideTets[x]->vertices[2]);
    edges[4] = createEdge(_insideTets[x]->vertices[2],
                          _insideTets[x]->vertices[3]);
    edges[5] = createEdge(_insideTets[x]->vertices[3],
                          _insideTets[x]->vertices[1]);

    // look for a cut along each edge
    int cutsTotal = 0;
    for (int y = 0; y < 6; y++)
    {
      map<EDGE, VEC3F*>::iterator cutFinder;
      cutFinder = _cutPoints.find(edges[y]);
      if (cutFinder != _cutPoints.end())
      {
        cuts[y] = cutFinder->second;
        cutsTotal++;
      }
    }
    // case 5
    if (outsideTotal == 2 && cutsTotal == 2)
    {
      generateCase5(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      continue;
    }
  }
  cout << " After case 5 ===============================" << endl;
  collectTetStats();

  // for each inside tet
  for (unsigned int x = 0; x < _insideTets.size(); x++)
  {
    // count number of inside/outside/snapped points
    int insidePoints[] = {-1,-1,-1,-1};
    int outsidePoints[] = {-1,-1,-1,-1};
    int snappedPoints[] = {-1,-1,-1,-1};
    int insideTotal = 0;
    int outsideTotal = 0;
    int snappedTotal = 0;
    for (int y = 0; y < 4; y++)
    {
      int vertexID = _vertexIDs[_insideTets[x]->vertices[y]];
      map<int,bool>::iterator finder = _snapped.find(vertexID);
      if (finder != _snapped.end())
      {
        snappedPoints[y] = vertexID;
        snappedTotal++;
      }
      else if (_inside[vertexID])
      {
        insidePoints[y] = vertexID;
        insideTotal++;
      }
      else
      {
        outsidePoints[y] = vertexID;
        outsideTotal++;
      }
    }

    // count number of edge cuts
    VEC3F* cuts[] = {NULL, NULL, NULL, NULL, NULL, NULL};
    EDGE edges[6];
    edges[0] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[1]);
    edges[1] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[2]);
    edges[2] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[3]);
    edges[3] = createEdge(_insideTets[x]->vertices[1],
                          _insideTets[x]->vertices[2]);
    edges[4] = createEdge(_insideTets[x]->vertices[2],
                          _insideTets[x]->vertices[3]);
    edges[5] = createEdge(_insideTets[x]->vertices[3],
                          _insideTets[x]->vertices[1]);

    // look for a cut along each edge
    int cutsTotal = 0;
    for (int y = 0; y < 6; y++)
    {
      map<EDGE, VEC3F*>::iterator cutFinder;
      cutFinder = _cutPoints.find(edges[y]);
      if (cutFinder != _cutPoints.end())
      {
        cuts[y] = cutFinder->second;
        cutsTotal++;
      }
    }
    // case 6
    if (insideTotal == 1 && cutsTotal == 3)
    {
      generateCase6(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      continue;
    }
  }
  cout << " After case 6 ===============================" << endl;
  collectTetStats();

  // for each inside tet
  for (unsigned int x = 0; x < _insideTets.size(); x++)
  {
    // count number of inside/outside/snapped points
    int insidePoints[] = {-1,-1,-1,-1};
    int outsidePoints[] = {-1,-1,-1,-1};
    int snappedPoints[] = {-1,-1,-1,-1};
    int insideTotal = 0;
    int outsideTotal = 0;
    int snappedTotal = 0;
    for (int y = 0; y < 4; y++)
    {
      int vertexID = _vertexIDs[_insideTets[x]->vertices[y]];
      map<int,bool>::iterator finder = _snapped.find(vertexID);
      if (finder != _snapped.end())
      {
        snappedPoints[y] = vertexID;
        snappedTotal++;
      }
      else if (_inside[vertexID])
      {
        insidePoints[y] = vertexID;
        insideTotal++;
      }
      else
      {
        outsidePoints[y] = vertexID;
        outsideTotal++;
      }
    }

    // count number of edge cuts
    VEC3F* cuts[] = {NULL, NULL, NULL, NULL, NULL, NULL};
    EDGE edges[6];
    edges[0] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[1]);
    edges[1] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[2]);
    edges[2] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[3]);
    edges[3] = createEdge(_insideTets[x]->vertices[1],
                          _insideTets[x]->vertices[2]);
    edges[4] = createEdge(_insideTets[x]->vertices[2],
                          _insideTets[x]->vertices[3]);
    edges[5] = createEdge(_insideTets[x]->vertices[3],
                          _insideTets[x]->vertices[1]);

    // look for a cut along each edge
    int cutsTotal = 0;
    for (int y = 0; y < 6; y++)
    {
      map<EDGE, VEC3F*>::iterator cutFinder;
      cutFinder = _cutPoints.find(edges[y]);
      if (cutFinder != _cutPoints.end())
      {
        cuts[y] = cutFinder->second;
        cutsTotal++;
      }
    }
    // cases 7 and 9
    if (outsideTotal == 1 && snappedTotal == 1)
    {
      // get the snapped vertex
      VEC3F* snappedVertex = NULL;
      for (int y = 0; y < 4; y++)
        if (snappedPoints[y] >= 0)
          snappedVertex = _vertices[snappedPoints[y]];

      // get the outside vertex
      VEC3F* outsideVertex = NULL;
      for (int y = 0; y < 4; y++)
        if (outsidePoints[y] >= 0)
          outsideVertex = _vertices[outsidePoints[y]];

      // make an edge
      EDGE edge = createEdge(snappedVertex, outsideVertex);

      // get the edge color
      map<EDGE, EDGE_TYPE>::iterator finder;
      finder = _edgeTypes.find(edge);

      // pick the case based on the edge color
      if (finder->second == RED)
        generateCase7(_insideTets[x], snappedPoints, insidePoints, outsidePoints, cuts, edges);
      else
        generateCase9(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      continue;
    }
  }
  cout << " After case 7 and 9 ===============================" << endl;
  collectTetStats();

  // for each inside tet
  for (unsigned int x = 0; x < _insideTets.size(); x++)
  {
    // count number of inside/outside/snapped points
    int insidePoints[] = {-1,-1,-1,-1};
    int outsidePoints[] = {-1,-1,-1,-1};
    int snappedPoints[] = {-1,-1,-1,-1};
    int insideTotal = 0;
    int outsideTotal = 0;
    int snappedTotal = 0;
    for (int y = 0; y < 4; y++)
    {
      int vertexID = _vertexIDs[_insideTets[x]->vertices[y]];
      map<int,bool>::iterator finder = _snapped.find(vertexID);
      if (finder != _snapped.end())
      {
        snappedPoints[y] = vertexID;
        snappedTotal++;
      }
      else if (_inside[vertexID])
      {
        insidePoints[y] = vertexID;
        insideTotal++;
      }
      else
      {
        outsidePoints[y] = vertexID;
        outsideTotal++;
      }
    }

    // count number of edge cuts
    VEC3F* cuts[] = {NULL, NULL, NULL, NULL, NULL, NULL};
    EDGE edges[6];
    edges[0] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[1]);
    edges[1] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[2]);
    edges[2] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[3]);
    edges[3] = createEdge(_insideTets[x]->vertices[1],
                          _insideTets[x]->vertices[2]);
    edges[4] = createEdge(_insideTets[x]->vertices[2],
                          _insideTets[x]->vertices[3]);
    edges[5] = createEdge(_insideTets[x]->vertices[3],
                          _insideTets[x]->vertices[1]);

    // look for a cut along each edge
    int cutsTotal = 0;
    for (int y = 0; y < 6; y++)
    {
      map<EDGE, VEC3F*>::iterator cutFinder;
      cutFinder = _cutPoints.find(edges[y]);
      if (cutFinder != _cutPoints.end())
      {
        cuts[y] = cutFinder->second;
        cutsTotal++;
      }
    }
    // case 8 and 11
    if (outsideTotal == 2 && insideTotal == 2)
    {
      // check the color of the edge between the two inside points
      vector<VEC3F*> whichInside;
      for (int y = 0; y < 4; y++)
        if (insidePoints[y] != -1)
          whichInside.push_back(_vertices[insidePoints[y]]);
      EDGE edge = createEdge(whichInside[0],
                             whichInside[1]);
      map<EDGE, EDGE_TYPE>::iterator finder;
      finder = _edgeTypes.find(edge);

      if (finder->second == RED)
        generateCase8(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      else
        generateCase11(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      continue;
    }
  }
  cout << " After case 8 and 11 ===============================" << endl;
  collectTetStats();

  // for each inside tet
  for (unsigned int x = 0; x < _insideTets.size(); x++)
  {
    // count number of inside/outside/snapped points
    int insidePoints[] = {-1,-1,-1,-1};
    int outsidePoints[] = {-1,-1,-1,-1};
    int snappedPoints[] = {-1,-1,-1,-1};
    int insideTotal = 0;
    int outsideTotal = 0;
    int snappedTotal = 0;
    for (int y = 0; y < 4; y++)
    {
      int vertexID = _vertexIDs[_insideTets[x]->vertices[y]];
      map<int,bool>::iterator finder = _snapped.find(vertexID);
      if (finder != _snapped.end())
      {
        snappedPoints[y] = vertexID;
        snappedTotal++;
      }
      else if (_inside[vertexID])
      {
        insidePoints[y] = vertexID;
        insideTotal++;
      }
      else
      {
        outsidePoints[y] = vertexID;
        outsideTotal++;
      }
    }

    // count number of edge cuts
    VEC3F* cuts[] = {NULL, NULL, NULL, NULL, NULL, NULL};
    EDGE edges[6];
    edges[0] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[1]);
    edges[1] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[2]);
    edges[2] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[3]);
    edges[3] = createEdge(_insideTets[x]->vertices[1],
                          _insideTets[x]->vertices[2]);
    edges[4] = createEdge(_insideTets[x]->vertices[2],
                          _insideTets[x]->vertices[3]);
    edges[5] = createEdge(_insideTets[x]->vertices[3],
                          _insideTets[x]->vertices[1]);

    // look for a cut along each edge
    int cutsTotal = 0;
    for (int y = 0; y < 6; y++)
    {
      map<EDGE, VEC3F*>::iterator cutFinder;
      cutFinder = _cutPoints.find(edges[y]);
      if (cutFinder != _cutPoints.end())
      {
        cuts[y] = (cutFinder->second);
        cutsTotal++;
      }
    }
    // case 10
    if (outsideTotal == 1 && insideTotal == 3)
    {
      generateCase10(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      continue;
    }
  }
  cout << " After case 10 ===============================" << endl;
  collectTetStats();
}

//////////////////////////////////////////////////////////////////////
// generate the final tets
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateFinalTets()
{
  // clear previous final tets
  _finalTets.clear();

  // for each inside tet
  for (unsigned int x = 0; x < _insideTets.size(); x++)
  {
    // count number of inside/outside/snapped points
    int insidePoints[] = {-1,-1,-1,-1};
    int outsidePoints[] = {-1,-1,-1,-1};
    int snappedPoints[] = {-1,-1,-1,-1};
    int insideTotal = 0;
    int outsideTotal = 0;
    int snappedTotal = 0;

    for (int y = 0; y < 4; y++)
    {
      int vertexID = _vertexIDs[_insideTets[x]->vertices[y]];
      map<int,bool>::iterator finder = _snapped.find(vertexID);
      if (finder != _snapped.end())
      {
        snappedPoints[y] = vertexID;
        snappedTotal++;
      }
      else if (_inside[vertexID])
      {
        insidePoints[y] = vertexID;
        insideTotal++;
      }
      else
      {
        outsidePoints[y] = vertexID;
        outsideTotal++;
      }
    }

    // count number of edge cuts
    VEC3F* cuts[] = {NULL, NULL, NULL, NULL, NULL, NULL};
    EDGE edges[6];
    edges[0] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[1]);
    edges[1] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[2]);
    edges[2] = createEdge(_insideTets[x]->vertices[0],
                          _insideTets[x]->vertices[3]);
    edges[3] = createEdge(_insideTets[x]->vertices[1],
                          _insideTets[x]->vertices[2]);
    edges[4] = createEdge(_insideTets[x]->vertices[2],
                          _insideTets[x]->vertices[3]);
    edges[5] = createEdge(_insideTets[x]->vertices[3],
                          _insideTets[x]->vertices[1]);

    // look for a cut along each edge
    int cutsTotal = 0;
    for (int y = 0; y < 6; y++)
    {
      map<EDGE, VEC3F*>::iterator cutFinder;
      cutFinder = _cutPoints.find(edges[y]);
      if (cutFinder != _cutPoints.end())
      {
        cuts[y] = (cutFinder->second);
        cutsTotal++;
      }
    }
    
    // case 0-3 in Figure 3 of paper
    if (outsideTotal == 0 && insideTotal > 0)
    {
      _finalTets.push_back(_insideTets[x]);
      continue;
    }

    // case 4
    if (outsideTotal == 1 && cutsTotal == 1)
    {
      generateCase4(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      continue;
    }

    // case 5
    if (outsideTotal == 2 && cutsTotal == 2)
    {
      generateCase5(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      continue;
    }

    // case 6
    if (insideTotal == 1 && cutsTotal == 3)
    {
      generateCase6(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      continue;
    }

    // cases 7 and 9
    if (outsideTotal == 1 && snappedTotal == 1)
    {
      // get the snapped vertex
      VEC3F* snappedVertex = NULL;
      for (int y = 0; y < 4; y++)
        if (snappedPoints[y] >= 0)
          snappedVertex = _vertices[snappedPoints[y]];

      // get the outside vertex
      VEC3F* outsideVertex = NULL;
      for (int y = 0; y < 4; y++)
        if (outsidePoints[y] >= 0)
          outsideVertex = _vertices[outsidePoints[y]];

      // make an edge
      EDGE edge = createEdge(snappedVertex, outsideVertex);

      // get the edge color
      map<EDGE, EDGE_TYPE>::iterator finder;
      finder = _edgeTypes.find(edge);

      // pick the case based on the edge color
      if (finder->second == RED)
        generateCase7(_insideTets[x], snappedPoints, insidePoints, outsidePoints, cuts, edges);
      else
        generateCase9(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      continue;
    }

    // case 8 and 11
    if (outsideTotal == 2 && insideTotal == 2)
    {
      // check the color of the edge between the two inside points
      vector<VEC3F*> whichInside;
      for (int y = 0; y < 4; y++)
        if (insidePoints[y] != -1)
          whichInside.push_back(_vertices[insidePoints[y]]);
      EDGE edge = createEdge(whichInside[0],
                             whichInside[1]);
      map<EDGE, EDGE_TYPE>::iterator finder;
      finder = _edgeTypes.find(edge);

      if (finder->second == RED)
        generateCase8(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      else
        generateCase11(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      continue;
    }
    // case 10
    if (outsideTotal == 1 && insideTotal == 3)
    {
      generateCase10(_insideTets[x], insidePoints, outsidePoints, cuts, edges);
      continue;
    }
    if (outsideTotal + snappedTotal != 4)
    {
      cout << __FILE__ << " " << __LINE__ << " CASE MISSED: " << endl;
      cout << __FILE__ << " " << __LINE__ << " insideTotal: " << insideTotal << endl;
      cout << __FILE__ << " " << __LINE__ << " outsideTotal: " << outsideTotal << endl;
      cout << __FILE__ << " " << __LINE__ << " snappedTotal: " << snappedTotal << endl;
      cout << __FILE__ << " " << __LINE__ << " cutsTotal: " << cutsTotal << endl;
    }
  }
  generateNodeMasses();
  cullUnusedNodes();
}

//////////////////////////////////////////////////////////////////////
// generate case 4 from Figure 3 of the paper
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateCase4(TET* tet, int insidePoints[4], int outsidePoints[4], 
                                VEC3F* cuts[6], EDGE edges[6])
{
  // if the inside point is v0
  if (insidePoints[0] >= 0)
  {
    TET* newTet = new TET(tet->vertices[0],
               tet->vertices[1],
               tet->vertices[2],
               tet->vertices[3]);
    if (cuts[0] != NULL) newTet->vertices[1] = cuts[0];
    if (cuts[1] != NULL) newTet->vertices[2] = cuts[1];
    if (cuts[2] != NULL) newTet->vertices[3] = cuts[2];
    
    _tets.push_back(newTet);
    _finalTets.push_back((_tets.back()));
    return;
  }

  // if the inside point is v1
  if (insidePoints[1] >= 0)
  {
    TET* newTet = new TET(tet->vertices[1],
               tet->vertices[2],
               tet->vertices[0],
               tet->vertices[3]);
    if (cuts[0] != NULL) newTet->vertices[2] = cuts[0];
    if (cuts[3] != NULL) newTet->vertices[1] = cuts[3];
    if (cuts[5] != NULL) newTet->vertices[3] = cuts[5];

    _tets.push_back(newTet);
    _finalTets.push_back((_tets.back()));
    return;
  }

  // if the inside point is v2
  if (insidePoints[2] >= 0)
  {
    TET* newTet = new TET(tet->vertices[2],
               tet->vertices[0],
               tet->vertices[1],
               tet->vertices[3]);
    if (cuts[1] != NULL) newTet->vertices[1] = cuts[1];
    if (cuts[3] != NULL) newTet->vertices[2] = cuts[3];
    if (cuts[4] != NULL) newTet->vertices[3] = cuts[4];

    _tets.push_back(newTet);
    _finalTets.push_back((_tets.back()));
    return;
  }

  // if the inside point is v3
  if (insidePoints[3] >= 0)
  {
    TET* newTet = new TET(tet->vertices[3],
               tet->vertices[0],
               tet->vertices[2],
               tet->vertices[1]);
    if (cuts[2] != NULL) newTet->vertices[1] = cuts[2];
    if (cuts[4] != NULL) newTet->vertices[2] = cuts[4];
    if (cuts[5] != NULL) newTet->vertices[3] = cuts[5];

    _tets.push_back(newTet);
    _finalTets.push_back((_tets.back()));
    return;
  }
}

//////////////////////////////////////////////////////////////////////
// generate case 5 from Figure 3 of the paper
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateCase5(TET* tet, int insidePoints[4], int outsidePoints[4], 
                                VEC3F* cuts[6], EDGE edges[6])
{
  // if the inside point is v0
  if (insidePoints[0] >= 0)
  {
    TET* newTet = new TET(tet->vertices[0],
               tet->vertices[1],
               tet->vertices[2],
               tet->vertices[3]);
    if (cuts[0] != NULL && cuts[2] != NULL)
    { 
      newTet->vertices[1] = cuts[0];
      newTet->vertices[3] = cuts[2];
    }
    if (cuts[1] != NULL && cuts[2] != NULL)
    { 
      newTet->vertices[2] = cuts[1];
      newTet->vertices[3] = cuts[2];
    }
    if (cuts[0] != NULL && cuts[1] != NULL)
    { 
      newTet->vertices[1] = cuts[0];
      newTet->vertices[2] = cuts[1];
    }

    _tets.push_back(newTet);
    _finalTets.push_back((_tets.back()));
    return;
  }

  // if the inside point is v1
  if (insidePoints[1] >= 0)
  {
    TET* newTet = new TET(tet->vertices[1],
               tet->vertices[2],
               tet->vertices[0],
               tet->vertices[3]);
    if (cuts[0] != NULL && cuts[3] != NULL)
    { 
      newTet->vertices[1] = cuts[3];
      newTet->vertices[2] = cuts[0];
    }
    if (cuts[0] != NULL && cuts[5] != NULL)
    { 
      newTet->vertices[2] = cuts[0];
      newTet->vertices[3] = cuts[5];
    }
    if (cuts[3] != NULL && cuts[5] != NULL)
    { 
      newTet->vertices[1] = cuts[3];
      newTet->vertices[3] = cuts[5];
    }

    _tets.push_back(newTet);
    _finalTets.push_back(_tets.back());
    return;
  }

  // if the inside point is v2
  if (insidePoints[2] >= 0)
  {
    TET* newTet = new TET(tet->vertices[2],
               tet->vertices[0],
               tet->vertices[1],
               tet->vertices[3]);
    if (cuts[1] != NULL && cuts[3] != NULL)
    { 
      newTet->vertices[1] = cuts[1];
      newTet->vertices[2] = cuts[3];
    }
    if (cuts[1] != NULL && cuts[4] != NULL)
    { 
      newTet->vertices[1] = cuts[1];
      newTet->vertices[3] = cuts[4];
    }
    if (cuts[3] != NULL && cuts[4] != NULL)
    { 
      newTet->vertices[2] = cuts[3];
      newTet->vertices[3] = cuts[4];
    }

    _tets.push_back(newTet);
    _finalTets.push_back(_tets.back());
    return;
  }

  // if the inside point is v3
  if (insidePoints[3] >= 0)
  {
    TET* newTet = new TET(tet->vertices[3],
               tet->vertices[0],
               tet->vertices[2],
               tet->vertices[1]);
    if (cuts[2] != NULL && cuts[4] != NULL)
    { 
      newTet->vertices[1] = cuts[2];
      newTet->vertices[2] = cuts[4];
    }
    if (cuts[2] != NULL && cuts[5] != NULL)
    { 
      newTet->vertices[1] = cuts[2];
      newTet->vertices[3] = cuts[5];
    }
    if (cuts[4] != NULL && cuts[5] != NULL)
    { 
      newTet->vertices[2] = cuts[4];
      newTet->vertices[3] = cuts[5];
    }

    _tets.push_back(newTet);
    _finalTets.push_back(_tets.back());
    return;
  }
}

//////////////////////////////////////////////////////////////////////
// generate case 6 from Figure 3 of the paper
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateCase6(TET* tet, int insidePoints[4], int outsidePoints[4], 
                                VEC3F* cuts[6], EDGE edges[6])
{
  // if the inside point is v0
  if (insidePoints[0] >= 0)
  {
    TET* newTet = new TET(tet->vertices[0],
               cuts[0],
               cuts[1],
               cuts[2]);

    _tets.push_back(newTet);
  }

  // if the inside point is v1
  if (insidePoints[1] >= 0)
  {
    TET* newTet = new TET(tet->vertices[1],
               cuts[3],
               cuts[0],
               cuts[5]);
    
    _tets.push_back(newTet);
  }

  // if the inside point is v2
  if (insidePoints[2] >= 0)
  {
    TET* newTet = new TET(tet->vertices[2],
               cuts[1],
               cuts[3],
               cuts[4]);

    _tets.push_back(newTet);
  }

  // if the inside point is v3
  if (insidePoints[3] >= 0)
  {
    TET* newTet = new TET(tet->vertices[3],
               cuts[2],
               cuts[4],
               cuts[5]);

    _tets.push_back(newTet);
  }
  _finalTets.push_back((_tets.back()));
  return;
}

//////////////////////////////////////////////////////////////////////
// generate case 7 from Figure 3 of the paper
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateCase7(TET* tet, int snappedPoints[4], int insidePoints[4], 
                                int outsidePoints[4], VEC3F* cuts[6], EDGE edges[6])
{
  EDGE edge;
  map<EDGE, EDGE_TYPE>::iterator finder;
  if (snappedPoints[0] >= 0)
  {
    edge = createEdge(tet->vertices[0], tet->vertices[1]);
    finder = _edgeTypes.find(edge);
    if (finder->second == BLACK)
    {
      if (outsidePoints[2] >= 0)
      {
        TET* newTet0 = new TET(tet->vertices[0],
                    cuts[3],
                    cuts[4],
                    tet->vertices[1]);
        TET* newTet1 = new TET(tet->vertices[0],
                    cuts[4],
                    tet->vertices[3],
                    tet->vertices[1]);
        _tets.push_back(newTet0);
        _finalTets.push_back((_tets.back()));
        _tets.push_back(newTet1);
        _finalTets.push_back((_tets.back()));
        return;
      }
      if (outsidePoints[3] >= 0)
      {
        TET* newTet0 = new TET(tet->vertices[0],
                    tet->vertices[1],
                    tet->vertices[2],
                    cuts[4]);
        TET* newTet1 = new TET(tet->vertices[0],
                    tet->vertices[1],            
                    cuts[4],
                    cuts[5]);
        _tets.push_back(newTet0);
        _finalTets.push_back((_tets.back()));
        _tets.push_back(newTet1);
        _finalTets.push_back((_tets.back()));
        return;
      }
    }
  }

  if (snappedPoints[1] >= 0)
  {
    // seems to like this first one
    edge = createEdge(tet->vertices[1], tet->vertices[0]);
    finder = _edgeTypes.find(edge);
    if (finder->second == BLACK)
    {
      if (outsidePoints[3] >= 0)
      {
        TET* newTet0 = new TET(tet->vertices[1],
                    cuts[4],
                    tet->vertices[0],
                    cuts[2]);
        TET* newTet1 = new TET(tet->vertices[1],
                    tet->vertices[2],
                    tet->vertices[0],
                    cuts[4]);
        _tets.push_back(newTet0);
        _finalTets.push_back((_tets.back()));
        _tets.push_back(newTet1);
        _finalTets.push_back((_tets.back()));
        return;
      }
      if (outsidePoints[2] >= 0)
      {
        TET* newTet0 = new TET(tet->vertices[1],
                    tet->vertices[0],
                    tet->vertices[3],
                    cuts[4]);
        TET* newTet1 = new TET(tet->vertices[1],
                    tet->vertices[0],
                    cuts[4],
                    cuts[1]);
        _tets.push_back(newTet0);
        _finalTets.push_back((_tets.back()));
        _tets.push_back(newTet1);
        _finalTets.push_back((_tets.back()));
        return;
      }
    }
  }
  if (snappedPoints[2] >= 0)
  {
    // seems to like this one
    edge = createEdge(tet->vertices[2], tet->vertices[3]);
    finder = _edgeTypes.find(edge);
    if (finder->second == BLACK)
    {
      if (outsidePoints[0] >= 0)
      {
        TET* newTet0 = new TET(tet->vertices[2],
                    cuts[0],
                    tet->vertices[3],
                    cuts[2]);
        TET* newTet1 = new TET(tet->vertices[2],
                    cuts[0],
                    tet->vertices[1],
                    tet->vertices[3]);
        _tets.push_back(newTet0);
        _finalTets.push_back((_tets.back()));
        _tets.push_back(newTet1);
        _finalTets.push_back((_tets.back()));
        return;
      }
      if (outsidePoints[1] >= 0)
      {
        TET* newTet0 = new TET(tet->vertices[2],
                    tet->vertices[0],
                    cuts[0],
                    tet->vertices[3]);
        TET* newTet1 = new TET(tet->vertices[2],
                    cuts[0],
                    cuts[5],
                    tet->vertices[3]);
        _tets.push_back(newTet0);
        _finalTets.push_back((_tets.back()));
        _tets.push_back(newTet1);
        _finalTets.push_back((_tets.back()));
        return;
      }
    }
  }
  if (snappedPoints[3] >= 0)
  {
    // seems to like this one
    edge = createEdge(tet->vertices[3], tet->vertices[2]);
    finder = _edgeTypes.find(edge);
    if (finder->second == BLACK)
    {
      if (outsidePoints[1] >= 0)
      {
        TET* newTet0 = new TET(tet->vertices[3],
                    cuts[0],
                    tet->vertices[2],
                    cuts[3]);
        TET* newTet1 = new TET(tet->vertices[3],
                    cuts[0],
                    tet->vertices[0],
                    tet->vertices[2]);
        _tets.push_back(newTet0);
        _finalTets.push_back((_tets.back()));
        _tets.push_back(newTet1);
        _finalTets.push_back((_tets.back()));
        return;
      }
      if (outsidePoints[0] >= 0)
      {
        TET* newTet0 = new TET(tet->vertices[3],
                    tet->vertices[1],
                    cuts[0],
                    tet->vertices[2]);
        TET* newTet1 = new TET(tet->vertices[3],
                    cuts[0],
                    cuts[1],
                    tet->vertices[2]);
        _tets.push_back(newTet0);
        _finalTets.push_back((_tets.back()));
        _tets.push_back(newTet1);
        _finalTets.push_back((_tets.back()));
        return;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// generate case 8 from Figure 3 of the paper
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateCase8(TET* tet, int insidePoints[4], int outsidePoints[4], 
                                VEC3F* cuts[6], EDGE edges[6])
{
  map<EDGE, VEC3F>::iterator cutFinder;
  map<EDGE, EDGE_TYPE>::iterator finder;
  EDGE edge;
  VEC3F* vertices[4];
  if (insidePoints[0] >= 0)
  {
    edge = createEdge(tet->vertices[0], tet->vertices[1]);
    finder = _edgeTypes.find(edge);
    if (finder->second == BLACK)
    {
      vertices[0] = tet->vertices[2];
      vertices[1] = tet->vertices[3];
      vertices[2] = tet->vertices[0];
      vertices[3] = tet->vertices[1];
      emitCase8(vertices);
      return;
    }
  }
  if (insidePoints[1] >= 0)
  {
    edge = createEdge(tet->vertices[1], tet->vertices[0]);
    finder = _edgeTypes.find(edge);
    if (finder->second == BLACK)
    {
      vertices[0] = tet->vertices[3];
      vertices[1] = tet->vertices[2];
      vertices[2] = tet->vertices[1];
      vertices[3] = tet->vertices[0];
      emitCase8(vertices);
      return;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// auxiliary function for case 8
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::emitCase8(VEC3F* vertices[4])
{
  // see which rotation it is
  map<EDGE, VEC3F*>::iterator cutFinder;
  EDGE edge = createEdge(vertices[0], vertices[2]);
  cutFinder = _cutPoints.find(edge);
  VEC3F* cuts[4];

  // rotation 0 case
  if (cutFinder != _cutPoints.end())
  {
    edge = createEdge(vertices[0], vertices[2]);
    cuts[0] = (_cutPoints[edge]);
    edge = createEdge(vertices[0], vertices[1]);
    cuts[1] = (_cutPoints[edge]);
    edge = createEdge(vertices[1], vertices[3]);
    cuts[2] = (_cutPoints[edge]);
    edge = createEdge(vertices[2], vertices[3]);
    cuts[3] = (_cutPoints[edge]);

    TET* newTet0 = new TET(cuts[3],
                cuts[1],
                vertices[2],
                vertices[1]);
    TET* newTet1 = new TET(cuts[3],
                cuts[1],
                vertices[1],
                cuts[2]);
    TET* newTet2 = new TET(cuts[0],
                cuts[1],
                vertices[2],
                cuts[3]);
    _tets.push_back(newTet0);
    _finalTets.push_back((_tets.back()));
    _tets.push_back(newTet1);
    _finalTets.push_back((_tets.back()));
    _tets.push_back(newTet2);
    _finalTets.push_back((_tets.back()));
  }
  // rotation 1 case
  else
  {
    edge = createEdge(vertices[0], vertices[3]);
    cuts[0] = (_cutPoints[edge]);
    edge = createEdge(vertices[0], vertices[1]);
    cuts[1] = (_cutPoints[edge]);
    edge = createEdge(vertices[1], vertices[2]);
    cuts[2] = (_cutPoints[edge]);
    edge = createEdge(vertices[2], vertices[3]);
    cuts[3] = (_cutPoints[edge]);

    TET* newTet0 = new TET(vertices[0],
                cuts[1],
                cuts[2],
                cuts[0]);
    TET* newTet1 = new TET(vertices[0],
                vertices[2],
                cuts[0],
                cuts[2]);
    TET* newTet2 = new TET(vertices[2],
                cuts[0],
                cuts[2],
                cuts[3]);
    _tets.push_back(newTet0);
    _finalTets.push_back((_tets.back()));
    _tets.push_back(newTet1);
    _finalTets.push_back((_tets.back()));
    _tets.push_back(newTet2);
    _finalTets.push_back((_tets.back()));
  }
}

//////////////////////////////////////////////////////////////////////
// generate case 9 from Figure 3 of the paper
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateCase9(TET* tet, int insidePoints[4], int outsidePoints[4], 
                                VEC3F* cuts[6], EDGE edges[6])
{
  map<EDGE, EDGE_TYPE>::iterator finder;
  if (insidePoints[0] >= 0 && insidePoints[1] >= 0)
  {
    // sanity check
    EDGE edge = createEdge(tet->vertices[0], tet->vertices[1]);
    finder = _edgeTypes.find(edge);
    if (finder->second != BLACK)
      cout << __FILE__ << " " << __LINE__ << " :  SANITY CHECK FAILED" << endl;

    if (outsidePoints[2] >= 0)
    {
      TET* newTet0 = new TET(tet->vertices[0],
                  tet->vertices[1],
                  cuts[1],
                  tet->vertices[3]);
      TET* newTet1 = new TET(cuts[1],
                  tet->vertices[1],
                  cuts[3],
                  tet->vertices[3]);
      _tets.push_back(newTet0);
      _finalTets.push_back((_tets.back()));
      _tets.push_back(newTet1);
      _finalTets.push_back((_tets.back()));
      return;
    }
    if (outsidePoints[3] >= 0)
    {
      TET* newTet0 = new TET(tet->vertices[0],
                  tet->vertices[1],
                  tet->vertices[2],
                  cuts[2]);
      TET* newTet1 = new TET(cuts[2],
                  tet->vertices[1],
                  tet->vertices[2],
                  cuts[5]);
      _tets.push_back(newTet0);
      _finalTets.push_back((_tets.back()));
      _tets.push_back(newTet1);
      _finalTets.push_back((_tets.back()));
      return;
    }
  }
  if (insidePoints[2] >= 0 && insidePoints[3] >= 0)
  {
    // sanity check
    EDGE edge = createEdge(tet->vertices[2], tet->vertices[3]);
    finder = _edgeTypes.find(edge);
    if (finder->second != BLACK)
      cout << __FILE__ << " " << __LINE__ << " :  SANITY CHECK FAILED" << endl;

    if (outsidePoints[0] >= 0)
    {
      TET* newTet0 = new TET(tet->vertices[2],
                  cuts[1],
                  tet->vertices[1],
                  tet->vertices[3]);
      TET* newTet1 = new TET(cuts[1],
                  cuts[2],
                  tet->vertices[1],
                  tet->vertices[3]);
      _tets.push_back(newTet0);
      _finalTets.push_back((_tets.back()));
      _tets.push_back(newTet1);
      _finalTets.push_back((_tets.back()));
      return;
    }
    if (outsidePoints[1] >= 0)
    {
      TET* newTet0 = new TET(tet->vertices[2],
                  tet->vertices[0],
                  cuts[3],
                  tet->vertices[3]);
      TET* newTet1 = new TET(cuts[3],
                  tet->vertices[0],
                  cuts[5],
                  tet->vertices[3]);
      _tets.push_back(newTet0);
      _finalTets.push_back((_tets.back()));
      _tets.push_back(newTet1);
      _finalTets.push_back((_tets.back()));
      return;
    }
  }
  if (insidePoints[0] >= 0 && insidePoints[2] >= 0)
  {
    if (outsidePoints[3] >= 0)
    {
      TET* newTet0 = new TET(tet->vertices[0],
                  tet->vertices[1],
                  tet->vertices[2],
                  cuts[2]);
      TET* newTet1 = new TET(cuts[2],
                  tet->vertices[1],
                  tet->vertices[2],
                  cuts[4]);
      _tets.push_back(newTet0);
      _finalTets.push_back((_tets.back()));
      _tets.push_back(newTet1);
      _finalTets.push_back((_tets.back()));
      return;
    }
    if (outsidePoints[1] >= 0)
    {
      TET* newTet0 = new TET(tet->vertices[0],
                  cuts[0],
                  tet->vertices[2],
                  tet->vertices[3]);
      TET* newTet1 = new TET(cuts[0],
                  cuts[3],
                  tet->vertices[2],
                  tet->vertices[3]);
      _tets.push_back(newTet0);
      _finalTets.push_back((_tets.back()));
      _tets.push_back(newTet1);
      _finalTets.push_back((_tets.back()));
      return;
    }
  }
  if (insidePoints[0] >= 0 && insidePoints[3] >= 0)
  {
    if (outsidePoints[1] >= 0)
    {
      TET* newTet0 = new TET(tet->vertices[0],
                  tet->vertices[2],
                  tet->vertices[3],
                  cuts[0]);
      TET* newTet1 = new TET(cuts[0],
                  tet->vertices[2],
                  tet->vertices[3],
                  cuts[5]);
      _tets.push_back(newTet0);
      _finalTets.push_back((_tets.back()));
      _tets.push_back(newTet1);
      _finalTets.push_back((_tets.back()));
      return;
    }
    if (outsidePoints[2] >= 0)
    {
      TET* newTet0 = new TET(tet->vertices[0],
                  cuts[1],
                  tet->vertices[3],
                  tet->vertices[1]);
      TET* newTet1 = new TET(cuts[1],
                  cuts[4],
                  tet->vertices[3],
                  tet->vertices[1]);
      _tets.push_back(newTet0);
      _finalTets.push_back((_tets.back()));
      _tets.push_back(newTet1);
      _finalTets.push_back((_tets.back()));
      return;
    }
  }
  if (insidePoints[1] >= 0 && insidePoints[2] >= 0)
  {
    if (outsidePoints[0] >= 0)
    {
      TET* newTet0 = new TET(tet->vertices[1],
                  tet->vertices[3],
                  tet->vertices[2],
                  cuts[0]);
      TET* newTet1 = new TET(cuts[0],
                  tet->vertices[3],
                  tet->vertices[2],
                  cuts[1]);
      _tets.push_back(newTet0);
      _finalTets.push_back((_tets.back()));
      _tets.push_back(newTet1);
      _finalTets.push_back((_tets.back()));
      return;
    }
    if (outsidePoints[3] >= 0)
    {
      TET* newTet0 = new TET(tet->vertices[1],
                  cuts[5],
                  tet->vertices[2],
                  tet->vertices[0]);
      TET* newTet1 = new TET(cuts[5],
                  cuts[4],
                  tet->vertices[2],
                  tet->vertices[0]);
      _tets.push_back(newTet0);
      _finalTets.push_back((_tets.back()));
      _tets.push_back(newTet1);
      _finalTets.push_back((_tets.back()));
      return;
    }
  }
  if (insidePoints[1] >= 0 && insidePoints[3] >= 0)
  {
    if (outsidePoints[2] >= 0)
    {
      TET* newTet0 = new TET(tet->vertices[1],
                  tet->vertices[0],
                  tet->vertices[3],
                  cuts[3]);
      TET* newTet1 = new TET(cuts[3],
                  tet->vertices[0],
                  tet->vertices[3],
                  cuts[4]);
      _tets.push_back(newTet0);
      _finalTets.push_back((_tets.back()));
      _tets.push_back(newTet1);
      _finalTets.push_back((_tets.back()));
      return;
    }
    if (outsidePoints[0] >= 0)
    {
      TET* newTet0 = new TET(tet->vertices[1],
                  cuts[0],
                  tet->vertices[3],
                  tet->vertices[2]);
      TET* newTet1 = new TET(cuts[0],
                  cuts[2],
                  tet->vertices[3],
                  tet->vertices[2]);
      _tets.push_back(newTet0);
      _finalTets.push_back((_tets.back()));
      _tets.push_back(newTet1);
      _finalTets.push_back((_tets.back()));
      return;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// generate case 10 from Figure 3 of the paper
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateCase10(TET* tet, int insidePoints[4], int outsidePoints[4], 
                                 VEC3F* cuts[6], EDGE edges[6])
{
  EDGE edge;
  VEC3F* vertices[4];
  VEC3F* cutsReordered[3];
  map<EDGE, EDGE_TYPE>::iterator finder;
  map<EDGE, VEC3F*>::iterator cutFinder;
  if (outsidePoints[0] >= 0)
  {
    edge = createEdge(tet->vertices[0], tet->vertices[1]);
    finder = _edgeTypes.find(edge);
    if (finder->second == BLACK)
    {
      vertices[0] = tet->vertices[2];
      vertices[1] = tet->vertices[3];
      vertices[2] = tet->vertices[0];
      vertices[3] = tet->vertices[1];

      edge = createEdge(tet->vertices[0], tet->vertices[2]);
      cutsReordered[0] = (_cutPoints[edge]);
      edge = createEdge(tet->vertices[0], tet->vertices[3]);
      cutsReordered[1] = (_cutPoints[edge]);
      edge = createEdge(tet->vertices[0], tet->vertices[1]);
      cutsReordered[2] = (_cutPoints[edge]);
      emitCase10(vertices, cutsReordered);
      return;
    }
    edge = createEdge(tet->vertices[0], tet->vertices[2]);
    finder = _edgeTypes.find(edge);
    if (finder->second == BLACK)
    {
      vertices[0] = tet->vertices[3];
      vertices[1] = tet->vertices[1];
      vertices[2] = tet->vertices[0];
      vertices[3] = tet->vertices[2];

      edge = createEdge(tet->vertices[0], tet->vertices[3]);
      cutsReordered[0] = (_cutPoints[edge]);
      edge = createEdge(tet->vertices[0], tet->vertices[1]);
      cutsReordered[1] = (_cutPoints[edge]);
      edge = createEdge(tet->vertices[0], tet->vertices[2]);
      cutsReordered[2] = (_cutPoints[edge]);
      emitCase10(vertices, cutsReordered);
      return;
    }
    edge = createEdge(tet->vertices[0], tet->vertices[3]);
    finder = _edgeTypes.find(edge);
    if (finder->second == BLACK)
    {
      vertices[0] = tet->vertices[1];
      vertices[1] = tet->vertices[2];
      vertices[2] = tet->vertices[0];
      vertices[3] = tet->vertices[3];

      edge = createEdge(tet->vertices[0], tet->vertices[1]);
      cutsReordered[0] = (_cutPoints[edge]);
      edge = createEdge(tet->vertices[0], tet->vertices[2]);
      cutsReordered[1] = (_cutPoints[edge]);
      edge = createEdge(tet->vertices[0], tet->vertices[3]);
      cutsReordered[2] = (_cutPoints[edge]);
      emitCase10(vertices, cutsReordered);
      return;
    }
  }
  if (outsidePoints[1] >= 0)
  {
    edge = createEdge(tet->vertices[1], tet->vertices[0]);
    finder = _edgeTypes.find(edge);
    if (finder->second == BLACK)
    {
      vertices[0] = tet->vertices[3];
      vertices[1] = tet->vertices[2];
      vertices[2] = tet->vertices[1];
      vertices[3] = tet->vertices[0];

      edge = createEdge(tet->vertices[1], tet->vertices[3]);
      cutsReordered[0] = (_cutPoints[edge]);
      edge = createEdge(tet->vertices[1], tet->vertices[2]);
      cutsReordered[1] = (_cutPoints[edge]);
      edge = createEdge(tet->vertices[1], tet->vertices[0]);
      cutsReordered[2] = (_cutPoints[edge]);
      emitCase10(vertices, cutsReordered);
      return;
    }
    edge = createEdge(tet->vertices[1], tet->vertices[2]);
    finder = _edgeTypes.find(edge);
    if (finder->second == BLACK)
    {
      vertices[0] = tet->vertices[0];
      vertices[1] = tet->vertices[3];
      vertices[2] = tet->vertices[1];
      vertices[3] = tet->vertices[2];

      edge = createEdge(tet->vertices[1], tet->vertices[0]);
      cutsReordered[0] = (_cutPoints[edge]);
      edge = createEdge(tet->vertices[1], tet->vertices[3]);
      cutsReordered[1] = (_cutPoints[edge]);
      edge = createEdge(tet->vertices[1], tet->vertices[2]);
      cutsReordered[2] = (_cutPoints[edge]);
      emitCase10(vertices, cutsReordered);
      return;
    }
    edge = createEdge(tet->vertices[1], tet->vertices[3]);
    finder = _edgeTypes.find(edge);
    if (finder->second == BLACK)
    {
      vertices[0] = tet->vertices[2];
      vertices[1] = tet->vertices[0];
      vertices[2] = tet->vertices[1];
      vertices[3] = tet->vertices[3];

      edge = createEdge(tet->vertices[1], tet->vertices[2]);
      cutsReordered[0] = (_cutPoints[edge]);
      edge = createEdge(tet->vertices[1], tet->vertices[0]);
      cutsReordered[1] = (_cutPoints[edge]);
      edge = createEdge(tet->vertices[1], tet->vertices[3]);
      cutsReordered[2] = (_cutPoints[edge]);
      emitCase10(vertices, cutsReordered);
      return;
    }
  }
  if (outsidePoints[2] >= 0)
  {
    edge = createEdge(tet->vertices[2], tet->vertices[0]);
    finder = _edgeTypes.find(edge);
    if (finder->second == BLACK)
    {
      vertices[0] = tet->vertices[1];
      vertices[1] = tet->vertices[3];
      vertices[2] = tet->vertices[2];
      vertices[3] = tet->vertices[0];

      edge = createEdge(tet->vertices[2], tet->vertices[1]);
      cutsReordered[0] = (_cutPoints[edge]);
      edge = createEdge(tet->vertices[2], tet->vertices[3]);
      cutsReordered[1] = (_cutPoints[edge]);
      edge = createEdge(tet->vertices[2], tet->vertices[0]);
      cutsReordered[2] = (_cutPoints[edge]);
      emitCase10(vertices, cutsReordered);
      return;
    }
    edge = createEdge(tet->vertices[2], tet->vertices[1]);
    finder = _edgeTypes.find(edge);
    if (finder->second == BLACK)
    {
      vertices[0] = tet->vertices[3];
      vertices[1] = tet->vertices[0];
      vertices[2] = tet->vertices[2];
      vertices[3] = tet->vertices[1];

      edge = createEdge(tet->vertices[2], tet->vertices[3]);
      cutsReordered[0] = (_cutPoints[edge]);
      edge = createEdge(tet->vertices[2], tet->vertices[0]);
      cutsReordered[1] = (_cutPoints[edge]);
      edge = createEdge(tet->vertices[2], tet->vertices[1]);
      cutsReordered[2] = (_cutPoints[edge]);
      emitCase10(vertices, cutsReordered);
      return;
    }
    edge = createEdge(tet->vertices[2], tet->vertices[3]);
    finder = _edgeTypes.find(edge);
    if (finder->second == BLACK)
    {
      vertices[0] = tet->vertices[0];
      vertices[1] = tet->vertices[1];
      vertices[2] = tet->vertices[2];
      vertices[3] = tet->vertices[3];

      edge = createEdge(tet->vertices[2], tet->vertices[0]);
      cutsReordered[0] = (_cutPoints[edge]);
      edge = createEdge(tet->vertices[2], tet->vertices[1]);
      cutsReordered[1] = (_cutPoints[edge]);
      edge = createEdge(tet->vertices[2], tet->vertices[3]);
      cutsReordered[2] = (_cutPoints[edge]);
      emitCase10(vertices, cutsReordered);
      return;
    }
  }
  if (outsidePoints[3] >= 0)
  {
    edge = createEdge(tet->vertices[3], tet->vertices[0]);
    finder = _edgeTypes.find(edge);
    if (finder->second == BLACK)
    {
      vertices[0] = tet->vertices[2];
      vertices[1] = tet->vertices[1];
      vertices[2] = tet->vertices[3];
      vertices[3] = tet->vertices[0];

      edge = createEdge(tet->vertices[3], tet->vertices[2]);
      cutsReordered[0] = (_cutPoints[edge]);
      edge = createEdge(tet->vertices[3], tet->vertices[1]);
      cutsReordered[1] = (_cutPoints[edge]);
      edge = createEdge(tet->vertices[3], tet->vertices[0]);
      cutsReordered[2] = (_cutPoints[edge]);
      emitCase10(vertices, cutsReordered);
      return;
    }
    edge = createEdge(tet->vertices[3], tet->vertices[1]);
    finder = _edgeTypes.find(edge);
    if (finder->second == BLACK)
    {
      vertices[0] = tet->vertices[0];
      vertices[1] = tet->vertices[2];
      vertices[2] = tet->vertices[3];
      vertices[3] = tet->vertices[1];

      edge = createEdge(tet->vertices[3], tet->vertices[0]);
      cutsReordered[0] = (_cutPoints[edge]);
      edge = createEdge(tet->vertices[3], tet->vertices[2]);
      cutsReordered[1] = (_cutPoints[edge]);
      edge = createEdge(tet->vertices[3], tet->vertices[1]);
      cutsReordered[2] = (_cutPoints[edge]);
      emitCase10(vertices, cutsReordered);
      return;
    }
    edge = createEdge(tet->vertices[3], tet->vertices[2]);
    finder = _edgeTypes.find(edge);
    if (finder->second == BLACK)
    {
      vertices[0] = tet->vertices[1];
      vertices[1] = tet->vertices[0];
      vertices[2] = tet->vertices[3];
      vertices[3] = tet->vertices[2];

      edge = createEdge(tet->vertices[3], tet->vertices[1]);
      cutsReordered[0] = (_cutPoints[edge]);
      edge = createEdge(tet->vertices[3], tet->vertices[0]);
      cutsReordered[1] = (_cutPoints[edge]);
      edge = createEdge(tet->vertices[3], tet->vertices[2]);
      cutsReordered[2] = (_cutPoints[edge]);
      emitCase10(vertices, cutsReordered);
      return;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// auxiliary function for case 10
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::emitCase10(VEC3F* vertices[4], VEC3F* cuts[3])
{
  VEC3F* parityVertex = vertices[0];
  VEC3F* parityCut = cuts[1];

  // pick the parity rule based on whether the vertex is a corner
  // or cell center -- paragraph 4 and 5 of Section 3.3
  CENTER_CORNER_TYPE centerCornerType = _centerCornerType[parityVertex];
  if (centerCornerType == NEITHER)
    cout << __FILE__ << " " << __LINE__ << " DUD CENTER/CORNER TYPE: " << endl;
  bool parityCase = (centerCornerType == CORNER) ? false : true;

  if (parityRule(parityVertex, parityCut) == parityCase)
  {
    TET* newTet0 = new TET(cuts[1],
                cuts[0],
                vertices[1],
                cuts[2]);
    TET* newTet1 = new TET(cuts[0],
                vertices[0],
                vertices[1],
                cuts[2]);
    TET* newTet2 = new TET(vertices[0],
                vertices[1],
                cuts[2],
                vertices[3]);
    _tets.push_back(newTet0);
    _finalTets.push_back((_tets.back()));
    _tets.push_back(newTet1);
    _finalTets.push_back((_tets.back()));
    _tets.push_back(newTet2);
    _finalTets.push_back((_tets.back()));
  }
  else
  {
    TET* newTet0 = new TET(cuts[0],
                cuts[1],
                cuts[2],
                vertices[3]);
    TET* newTet1 = new TET(cuts[0],
                cuts[1],
                vertices[3],
                vertices[0]);
    TET* newTet2 = new TET(vertices[0],
                vertices[1],
                cuts[1],
                vertices[3]);
    _tets.push_back(newTet0);
    _finalTets.push_back((_tets.back()));
    _tets.push_back(newTet1);
    _finalTets.push_back((_tets.back()));
    _tets.push_back(newTet2);
    _finalTets.push_back((_tets.back()));
  }
}

//////////////////////////////////////////////////////////////////////
// generate case 11 from Figure 3 of the paper
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateCase11(TET* tet, int insidePoints[4], int outsidePoints[4], 
                                 VEC3F* cuts[6], EDGE edges[6])
{
  EDGE edge;
  VEC3F* vertices[4];
  map<EDGE, EDGE_TYPE>::iterator finder;
  if (outsidePoints[0] >= 0)
  {
    edge = createEdge(tet->vertices[0], tet->vertices[1]);
    finder = _edgeTypes.find(edge);
    if (finder->second == BLACK)
    {
      vertices[0] = tet->vertices[2];
      vertices[1] = tet->vertices[3];
      vertices[2] = tet->vertices[0];
      vertices[3] = tet->vertices[1];
      emitCase11(vertices);
      return;
    }
  }
  if (outsidePoints[2] >= 0)
  {
    edge = createEdge(tet->vertices[2], tet->vertices[1]);
    finder = _edgeTypes.find(edge);
    if (finder->second == BLACK)
    {
      vertices[0] = tet->vertices[3];
      vertices[1] = tet->vertices[0];
      vertices[2] = tet->vertices[2];
      vertices[3] = tet->vertices[1];
      emitCase11(vertices);
      return;
    }
    edge = createEdge(tet->vertices[2], tet->vertices[3]);
    finder = _edgeTypes.find(edge);
    if (finder->second == BLACK)
    {
      vertices[0] = tet->vertices[0];
      vertices[1] = tet->vertices[1];
      vertices[2] = tet->vertices[2];
      vertices[3] = tet->vertices[3];
      emitCase11(vertices);
      return;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// auxiliary function for case 11
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::emitCase11(VEC3F* vertices[4])
{
  VEC3F* cuts[4];
  map<EDGE, VEC3F*>::iterator cutFinder;
  EDGE edges[4];

  edges[0] = createEdge(vertices[0], vertices[2]);
  edges[1] = createEdge(vertices[1], vertices[2]);
  edges[2] = createEdge(vertices[1], vertices[3]);
  edges[3] = createEdge(vertices[0], vertices[3]);

  bool dud = false;
  for (int x = 0; x < 4; x++)
  {
    cutFinder = _cutPoints.find(edges[x]);
    if (cutFinder == _cutPoints.end())
      dud = true;
  }

  if (dud == true)
  {
    cout << __FILE__ << " " << __LINE__ << " DUD CUT " << endl;
    EDGE dudEdges[6];
    dudEdges[0] = createEdge(vertices[0], vertices[1]);
    dudEdges[1] = createEdge(vertices[0], vertices[2]);
    dudEdges[2] = createEdge(vertices[0], vertices[3]);
    dudEdges[3] = createEdge(vertices[1], vertices[2]);
    dudEdges[4] = createEdge(vertices[2], vertices[3]);
    dudEdges[5] = createEdge(vertices[1], vertices[3]);

    cout << " CUT POINTS: " << endl;
    for (int x = 0; x < 6; x++)
    {
      cutFinder = _cutPoints.find(dudEdges[x]);
      if (cutFinder == _cutPoints.end())
        cout << "NULL" << endl;
      else
        cout << "CUT" << endl;
    }
    cout << " EDGE COLORS: " << endl;
    map<EDGE, EDGE_TYPE>::iterator colorFinder;
    for (int x = 0; x < 6; x++)
      cout << _edgeTypes[dudEdges[x]] << endl;

    cout << " INSIDE/OUTSIDE: " << endl;
    for (int x = 0; x < 4; x++)
      cout << _inside[_vertexIDs[vertices[x]]] << endl;
    return;
  }

  cuts[0] = (_cutPoints[edges[0]]);
  cuts[1] = (_cutPoints[edges[1]]);
  cuts[2] = (_cutPoints[edges[2]]);
  cuts[3] = (_cutPoints[edges[3]]);

  VEC3F* parityVertex = vertices[0];
  VEC3F* parityCut = cuts[1];

  // pick the parity rule based on whether the vertex is a corner
  // or cell center -- paragraph 4 and 5 of Section 3.3
  CENTER_CORNER_TYPE centerCornerType = _centerCornerType[vertices[0]];
  if (centerCornerType == NEITHER)
    cout << __FILE__ << " " << __LINE__ << " DUD CENTER/CORNER TYPE: " << endl;
  bool parityCase = (centerCornerType == CORNER) ? true : false;

  if (parityRule(parityVertex, parityCut) == parityCase)
  {
    TET* newTet0 = new TET(cuts[0],
                cuts[1],
                cuts[2],
                vertices[1]);
    TET* newTet1 = new TET(cuts[2],
                cuts[0],
                vertices[1],
                vertices[0]);
    TET* newTet2 = new TET(cuts[3],
                cuts[0],
                cuts[2],
                vertices[0]);
    _tets.push_back(newTet0);
    _finalTets.push_back((_tets.back()));
    _tets.push_back(newTet1);
    _finalTets.push_back((_tets.back()));
    _tets.push_back(newTet2);
    _finalTets.push_back((_tets.back()));
  }
  else
  {
    TET* newTet0 = new TET(cuts[3],
                cuts[1],
                cuts[2],
                vertices[1]);
    TET* newTet1 = new TET(cuts[3],
                cuts[1],
                vertices[1],
                vertices[0]);
    TET* newTet2 = new TET(cuts[3],
                cuts[0],
                cuts[1],
                vertices[0]);
    _tets.push_back(newTet0);
    _finalTets.push_back((_tets.back()));
    _tets.push_back(newTet1);
    _finalTets.push_back((_tets.back()));
    _tets.push_back(newTet2);
    _finalTets.push_back((_tets.back()));
  }
}


//////////////////////////////////////////////////////////////////////
// implements the "parity rule"
// if "vertex" has an even number of components greater than "cut",
// returns TRUE
//////////////////////////////////////////////////////////////////////
bool ISO_STUFFER::parityRule(VEC3F* vertex, VEC3F* cut)
{
  VEC3F diff = (*vertex) - (*cut);
  int greater = 0;
  if (diff[0] > 0) greater++;
  if (diff[1] > 0) greater++;
  if (diff[2] > 0) greater++;

  return greater % 2;
}

//////////////////////////////////////////////////////////////////////
// Run some tests on the final tets
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::collectTetStats()
{
  cout << "==================================================" << endl;
  cout << "FINAL MESH STATISTICS" << endl;
  cout << "==================================================" << endl;

  // check for invalid tets
  int invalid = 0;
  for (unsigned int x = 0; x < _finalTets.size(); x++)
  {
    if (_finalTets[x]->invalid())
      invalid++;
  }
  cout << " Invalid tets: " << invalid << " of " << _finalTets.size() << endl;

  // compute volume
  Real minVolume = 0;
  Real averageVolume = 0.0f;
  for (unsigned int x = 0; x < _finalTets.size(); x++)
  {
    averageVolume += _finalTets[x]->volume();
    if (x == 0)
      minVolume = _finalTets[0]->volume();

    if (_finalTets[x]->volume() < minVolume)
      minVolume = _finalTets[x]->volume();
  }
  averageVolume /= _finalTets.size();

  Real regularVolume = _insideTets[0]->volume();
  cout << " Smallest tet volume found: " << minVolume << "\t(" << 100.0f * minVolume / regularVolume << "\% of regular)" << endl;
  cout << " Average tet volume: " << averageVolume << "\t(" << 100.0f * averageVolume / regularVolume << "\% of regular)" << endl;

  Real minEdgeLength = 0;
  Real maxEdgeLength = 0;
  for (unsigned int x = 0; x < _finalTets.size(); x++)
  {
    TET& tet = *(_finalTets[x]);

    // compute normals first
    VEC3F edges[6];

    edges[0] = *(tet.vertices[0]) - *(tet.vertices[1]);
    edges[1] = *(tet.vertices[0]) - *(tet.vertices[2]);
    edges[2] = *(tet.vertices[0]) - *(tet.vertices[3]);
    edges[3] = *(tet.vertices[1]) - *(tet.vertices[2]);
    edges[4] = *(tet.vertices[1]) - *(tet.vertices[3]);
    edges[5] = *(tet.vertices[2]) - *(tet.vertices[3]);

    Real thisMin = 0;
    Real thisMax = 0;
    for (int y = 0; y < 6; y++)
    {
      Real length = norm(edges[y]);
      if (y == 0)
      {
        thisMin = length;
        thisMax = length;
      }
      thisMin = (length < thisMin) ? length : thisMin;
      thisMax = (length > thisMax) ? length : thisMax;
    }

    // initialize
    if (x == 0)
    {
      minEdgeLength = thisMin;
      maxEdgeLength = thisMax;
    }
    minEdgeLength = (thisMin < minEdgeLength) ? thisMin : minEdgeLength;
    maxEdgeLength = (thisMax > minEdgeLength) ? thisMax : maxEdgeLength;
  }
  cout << " Minimum edge length found: " << minEdgeLength;
  cout << "\t(" << 100.0f * minEdgeLength / (_redLength * _alphaShort)  << "\% of allowable)" << endl;
  cout << " Maxmimum edge length found: " << maxEdgeLength << endl;

  Real minDihedral = 0;
  Real maxDihedral = 0;
  Real averageDihedral = 0.0f;
  for (unsigned int x = 0; x < _finalTets.size(); x++)
  {
    TET& tet = *(_finalTets[x]);

    Real dihedrals[6];
    tet.dihedrals(&dihedrals[0]);

    // initialize
    if (x == 0)
    {
      minDihedral = dihedrals[0];
      maxDihedral = dihedrals[0];
    }

    for (int y = 0; y < 6; y++)
    {
      minDihedral = (minDihedral < dihedrals[y]) ? minDihedral : dihedrals[y];
      maxDihedral = (maxDihedral > dihedrals[y]) ? maxDihedral : dihedrals[y];
      averageDihedral += dihedrals[y];
    }
  }
  cout << " Min dihedral angle: " << 360.0f * minDihedral / (2.0f * M_PI) << " (should not be less than " << _proofMin << ") " << endl;
  cout << " Max dihedral angle: " << 360.0f * maxDihedral / (2.0f * M_PI) << " (should not be greater than " << _proofMax << ") " << endl;
  cout << " Average dihedral angle: " << 360.0f * (averageDihedral / (_finalTets.size() * 6)) / (2.0f * M_PI) << endl;
}

//////////////////////////////////////////////////////////////////////
// Run some tests on the final tets
//////////////////////////////////////////////////////////////////////
bool ISO_STUFFER::satisfiesProof(TET& tet, int caseNumber)
{
  if (tet.minDihedral() * 180.0f / M_PI < _proofMin)
  {
    cout << " MIN ANGLE VIOLATED BY CASE " << caseNumber << ": " << tet.minDihedral() * 180.0f / M_PI << endl;
    return false;
  }
  if (tet.maxDihedral() * 180.0f / M_PI > _proofMax)
  {
    cout << " MAX ANGLE VIOLATED BY CASE " << caseNumber << ": " << tet.maxDihedral() * 180.0f / M_PI << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////
// Constrain the nodes with minimum axis value
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::constrainMin(int axis)
{
  // find the min axis value
  VEC3F* vertex = _vertices.begin()->second;
  Real leastFound = (*vertex)[axis];
  map<int,VEC3F*>::iterator i;
  for (i = _vertices.begin(); i != _vertices.end(); i++)
  {
    vertex = i->second;
    if ((*vertex)[axis] < leastFound)
      leastFound = (*vertex)[axis];
  }

  Real boxBounds[] = {1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f};
  //boxBounds[2 * axis] = leastFound + _blackLength * 2.0;
  boxBounds[2 * axis] = leastFound + _blackLength;
  //boxBounds[2 * axis] = leastFound + _blackLength * 0.25;
  boxBounds[2 * axis + 1] = leastFound;
  BOX box(boxBounds[0], boxBounds[1], 
          boxBounds[2], boxBounds[3], 
          boxBounds[4], boxBounds[5]);
  generateConstrainedNodes(box);
}

//////////////////////////////////////////////////////////////////////
// Create an unconstrained mesh
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::constrainNone()
{
  if (_vertices.size() == 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << "No vertices found!!! " << endl;
  }

  // find the min x axis value
  VEC3F* vertex = _vertices.begin()->second;
  Real leastFound = (*vertex)[0];
  map<int,VEC3F*>::iterator i;
  for (i = _vertices.begin(); i != _vertices.end(); i++)
  {
    vertex = i->second;
    if ((*vertex)[0] < leastFound)
      leastFound = (*vertex)[0];
  }

  // pass a constraint box that misses the mesh entirely
  Real boxBounds[] = {1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f};
  boxBounds[0] = leastFound - 9.0f * _blackLength;
  boxBounds[1] = leastFound - 10.0f * _blackLength;
  BOX box(boxBounds[0], boxBounds[1], 
          boxBounds[2], boxBounds[3], 
          boxBounds[4], boxBounds[5]);
  generateConstrainedNodes(box);
}

//////////////////////////////////////////////////////////////////////
// Constrain the nodes with maximum axis value
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::constrainMax(int axis)
{
  // find the min axis value
  VEC3F* vertex = _vertices.begin()->second;
  Real greatestFound = (*vertex)[axis];
  map<int,VEC3F*>::iterator i;
  for (i = _vertices.begin(); i != _vertices.end(); i++)
  {
    vertex = i->second;
    if ((*vertex)[axis] > greatestFound)
      greatestFound = (*vertex)[axis];
  }

  Real boxBounds[] = {1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f};
  boxBounds[2 * axis] = greatestFound;
  boxBounds[2 * axis + 1] = greatestFound - _blackLength * 0.25;
  BOX box(boxBounds[0], boxBounds[1], 
          boxBounds[2], boxBounds[3], 
          boxBounds[4], boxBounds[5]);
  generateConstrainedNodes(box);
}

//////////////////////////////////////////////////////////////////////
// Constrain the nodes along both the minimum and maximum axis value
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::constrainMinMax(int axis)
{
  // find the min and max axis values
  VEC3F* vertex = _vertices.begin()->second;
  Real greatestFound = (*vertex)[axis];
  Real leastFound = (*vertex)[axis];
  map<int,VEC3F*>::iterator i;
  for (i = _vertices.begin(); i != _vertices.end(); i++)
  {
    vertex = i->second;
    if ((*vertex)[axis] > greatestFound)
      greatestFound = (*vertex)[axis];
    if ((*vertex)[axis] < leastFound)
      leastFound = (*vertex)[axis];
  }

  Real maxBoxBounds[] = {1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f};
  maxBoxBounds[2 * axis] = greatestFound;
  maxBoxBounds[2 * axis + 1] = greatestFound - _blackLength * 0.25;

  Real minBoxBounds[] = {1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f};
  minBoxBounds[2 * axis] = leastFound + _blackLength * 0.25;
  minBoxBounds[2 * axis + 1] = leastFound;

  // join the two boxes into the same object
  COMPOUND compound;
  compound.addSurface(new BOX(maxBoxBounds[0], maxBoxBounds[1], 
                              maxBoxBounds[2], maxBoxBounds[3], 
                              maxBoxBounds[4], maxBoxBounds[5]));
  compound.addSurface(new BOX(minBoxBounds[0], minBoxBounds[1], 
                              minBoxBounds[2], minBoxBounds[3], 
                              minBoxBounds[4], minBoxBounds[5]));

  generateConstrainedNodes(compound);
}

//////////////////////////////////////////////////////////////////////
// generate the constrained and unconstrained vertex lists
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateConstrainedNodes(SURFACE& surface)
{
  // clear the previous entries
  _constrainedNodes.clear();
  _unconstrainedNodes.clear();
  _constrainedIDs.clear();
  _unconstrainedIDs.clear();

  // for each vertex, see if it's inside the constraint
  // and then add it to the right list
  map<int,VEC3F*>::iterator i;
  for (i = _vertices.begin(); i != _vertices.end(); i++)
    if (surface.inside(*(i->second)))
      _constrainedNodes.push_back(i->second);
    else
      _unconstrainedNodes.push_back(i->second);

  map<EDGE,VEC3F*>::iterator cut;
  for (cut = _cutPoints.begin(); cut != _cutPoints.end(); cut++)
  {
    if (surface.inside(*(cut->second)))
      _constrainedNodes.push_back(cut->second);
    else
      _unconstrainedNodes.push_back(cut->second);
  }

  // generate the inverse mapping to that the tets can later
  // figure out which indices to output
  for (unsigned int x = 0; x < _constrainedNodes.size(); x++)
    _constrainedIDs[_constrainedNodes[x]] = x;

  for (unsigned int x = 0; x < _unconstrainedNodes.size(); x++)
    _unconstrainedIDs[_unconstrainedNodes[x]] = x;

  generateNodeMasses();
}

//////////////////////////////////////////////////////////////////////
// write out the final file - legacy support for files created for
// the deadline
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::writeFileLegacy(const char* filename)
{
  if (_unconstrainedNodes.size() == 0 &&
      _constrainedNodes.size() == 0)
  {
    cout << " No nodes to output! Did you run generateConstrainedNodes()?" << endl;
    return;
  }

  ofstream in(filename);
  if (in.fail () )
  {
    cerr << "Can't open tetmesh file " << filename << endl;
    return;
  }

  // calculate the minimum bounding box
  map<int,VEC3F*>::iterator i;
  VEC3F minBox;
  for (i = _vertices.begin(); i != _vertices.end(); i++)
  {
    if (i == _vertices.begin())
      minBox = *(i->second);
    minBox[0] = (minBox[0] < (*(i->second))[0]) ? minBox[0] : (*(i->second))[0];
    minBox[1] = (minBox[1] < (*(i->second))[1]) ? minBox[1] : (*(i->second))[1];
    minBox[2] = (minBox[2] < (*(i->second))[2]) ? minBox[2] : (*(i->second))[2];
  }

  // output the min box
  in << minBox << endl;
  Real radius = 1.0f;
  in << radius << endl;

  // output the unconstrained nodes -- must write to 32 bit int so that
  // 64 bit windows doesn't think it's an int_64
  unsigned int size = _unconstrainedNodes.size();
  in << size << endl;
  for (unsigned int x = 0; x < _unconstrainedNodes.size(); x++)
  {
    VEC3F* node = _unconstrainedNodes[x];
    in << (*node) << endl;

    // divide by 12 since each node has at most 12
    // incident tets
    Real mass = _masses[node];
    in << mass << endl;
  }

  // output the constrained nodes
  size = _constrainedNodes.size();
  in << size << endl;
  for (unsigned int x = 0; x < _constrainedNodes.size(); x++)
  {
    VEC3F* node = _constrainedNodes[x];
    in << (*node) << endl;
  }

  // output the tets
  size = _finalTets.size();
  in << size << endl;
  // for each tet
  for (unsigned int x = 0; x < _finalTets.size(); x++)
  {
    // for each node in the tet
    for (int y = 0; y < 4; y++)
    {
      // see if it is constrained
      map<VEC3F*, int>::iterator finder;
      finder = _constrainedIDs.find(_finalTets[x]->vertices[y]);
      bool constrained = false;
      
      if (finder != _constrainedIDs.end())
      {
        constrained = true;
        in << constrained << "\t";
        in << finder->second << endl;
      }
      else
      {
        // sanity check
        finder = _unconstrainedIDs.find(_finalTets[x]->vertices[y]);
        if (finder == _unconstrainedIDs.end())
          cout << __FILE__ << " " << __LINE__ << " SANITY CHECK FAILED: " << endl;

        in << constrained << "\t";
        in << _unconstrainedIDs[_finalTets[x]->vertices[y]] << endl;
      }
    }
    VEC3F v0 = *(_finalTets[x]->vertices[0]);
    VEC3F v1 = *(_finalTets[x]->vertices[1]);
    VEC3F v2 = *(_finalTets[x]->vertices[2]);
    VEC3F v3 = *(_finalTets[x]->vertices[3]);
    in << "1 2 3" << endl;
    in << cross(v2 - v1, v3 - v1) << endl;
    in << "0 3 2" << endl;
    in << cross(v2 - v0, v3 - v0) << endl;
    in << "0 1 3" << endl;
    in << cross(v3 - v0, v1 - v0) << endl;
    in << "0 2 1" << endl;
    in << cross(v1 - v0, v2 - v0) << endl;
  }

  // still need to write the face stuff here
  in.close();
}

//////////////////////////////////////////////////////////////////////
// write out the final file
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::writeFile(const char* filename)
{
  if (_unconstrainedNodes.size() == 0 &&
      _constrainedNodes.size() == 0)
  {
    cout << " No nodes to output! Did you run generateConstrainedNodes()?" << endl;
    return;
  }

  FILE* file = fopen(filename, "wb");

  cout << " Total nodes: " << _unconstrainedNodes.size() + _constrainedNodes.size() << endl;
  cout << " Total tets:  " << _finalTets.size() << endl;

  // output vertex array sizes
  int size = _unconstrainedNodes.size();
  fwrite((void*)&size, sizeof(int), 1, file);
  size = _constrainedNodes.size();
  fwrite((void*)&size, sizeof(int), 1, file);

  // output vertex positions
  for (unsigned int x = 0; x < _unconstrainedNodes.size(); x++)
  {
    VEC3F* node = _unconstrainedNodes[x];
    for (int y = 0; y < 3; y++)
    {
      double position = (*node)[y];
      fwrite((void*)&position, sizeof(double), 1, file);
    }
  }
  for (unsigned int x = 0; x < _constrainedNodes.size(); x++)
  {
    VEC3F* node = _constrainedNodes[x];
    for (int y = 0; y < 3; y++)
    {
      double position = (*node)[y];
      fwrite((void*)&position, sizeof(double), 1, file);
    }
  }

  // output tet vertex lists
  int unconstrainedSize = _unconstrainedNodes.size();
  int totalTets = _finalTets.size();
  fwrite((void*)&(totalTets), sizeof(int), 1, file);
  for (unsigned int x = 0; x < _finalTets.size(); x++)
  {
    // for each node in the tet
    for (int y = 0; y < 4; y++)
    {
      // see if it is constrained
      map<VEC3F*, int>::iterator finder;
      finder = _constrainedIDs.find(_finalTets[x]->vertices[y]);
     
      int index; 
      if (finder != _constrainedIDs.end())
        index = unconstrainedSize + finder->second;
      else
        index = _unconstrainedIDs[_finalTets[x]->vertices[y]];
      fwrite((void*)&(index), sizeof(int), 1, file);
    }
  }

  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// generate a mass for each node
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateNodeMasses()
{
  _masses.clear();

  for (unsigned int x = 0; x < _finalTets.size(); x++)
    for (int y = 0; y < 4; y++)
      _masses[_finalTets[x]->vertices[y]] =  1.0f / _unconstrainedNodes.size();
}

//////////////////////////////////////////////////////////////////////
// cull away the tets not connected to the rest
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::cullOrphanTets()
{
  if (_constrainedNodes.size() == 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Need a constrained starting point in order to find orphans!" << endl;
    return;
  }
  cout << "Culling orphan tets and vertices ... ";
  flush(cout);

  // record the connectivity of all the vertices
  _connectivity.clear();
  for (unsigned int x = 0; x < _finalTets.size(); x++)
    recordConnectivity(_finalTets[x]);

  // crawl the connected vertices
  _crawled.clear();
  crawlConnectedVertices(_constrainedNodes[0]);
  cout << _crawled.size() << " of " << _constrainedNodes.size() + _unconstrainedNodes.size()
                          << " nodes crawled " << endl;

  cout << " Crawling connected faces ... "; flush(cout);
  crawlConnectedFaces();
  cout << _crawledTets.size() << " of " << _finalTets.size() << " tets found " << endl;

  // find the orphans, ie the ones that contain an unconnected vertex
  map<TET*, bool> orphanTets;
  for (unsigned int x = 0; x < _finalTets.size(); x++)
  {
    // see if the vertices are orphans
    for (int y = 0; y < 4; y++)
      if (_crawled.find(_finalTets[x]->vertices[y]) == _crawled.end())
      {
        orphanTets[_finalTets[x]] = true;
        y = 4;
      }

    // see if it's not face connected
    if (_crawledTets.find(x) == _crawledTets.end())
      orphanTets[_finalTets[x]] = true;
  }
  cout << orphanTets.size() << " orphan tets found " << endl;

  // find the orphan vertices that only belong to orphan tets
  map<VEC3F*, bool> orphanVertices;
  for (unsigned int x = 0; x < _unconstrainedNodes.size(); x++)
  {
    // get its tet membership list
    vector<int>& membership = _tetMembership[_unconstrainedNodes[x]];
  
    // if it is a member of a non-orphan tet 
    bool orphan = true;
    for (unsigned int y = 0; y < membership.size(); y++)
      if (orphanTets.find(_finalTets[membership[y]]) == orphanTets.end())
      {
        orphan = false;
        break;
      }

    // else, tag it as an orphan
    if (orphan)
      orphanVertices[_unconstrainedNodes[x]] = true;
  }
  cout << orphanVertices.size() << " orphan vertices found " << endl;

  // delete the orphan vertices
  vector<VEC3F*> unconstrainedOld = _unconstrainedNodes;
  _unconstrainedNodes.clear();
  for (unsigned int x = 0; x < unconstrainedOld.size(); x++)
    if (_crawled.find(unconstrainedOld[x]) != _crawled.end() &&
        orphanVertices.find(unconstrainedOld[x]) == orphanVertices.end())
      _unconstrainedNodes.push_back(unconstrainedOld[x]);
  vector<VEC3F*> constrainedOld = _constrainedNodes;
  _constrainedNodes.clear();
  for (unsigned int x = 0; x < constrainedOld.size(); x++)
    if (_crawled.find(constrainedOld[x]) != _crawled.end() &&
        orphanVertices.find(constrainedOld[x]) == orphanVertices.end())
      _constrainedNodes.push_back(constrainedOld[x]);

  // delete the orphan tets
  vector<TET*> finalTetsOld = _finalTets;
  _finalTets.clear();
  for (unsigned int x = 0; x < finalTetsOld.size(); x++)
    if (orphanTets.find(finalTetsOld[x]) == orphanTets.end())
      _finalTets.push_back(finalTetsOld[x]);
  
  // recompute the constrained and unconstrained IDs
  _constrainedIDs.clear();
  for (unsigned int x = 0; x < _constrainedNodes.size(); x++)
    _constrainedIDs[_constrainedNodes[x]] = x;
  _unconstrainedIDs.clear();
  for (unsigned int x = 0; x < _unconstrainedNodes.size(); x++)
    _unconstrainedIDs[_unconstrainedNodes[x]] = x;

  // build a list of the non-face-shared nodes
  _nonSharedFaceVertices.clear();
  for (unsigned int x = 0; x < _unconstrainedNodes.size(); x++)
    if (_sharedFaceMembership.find(_unconstrainedNodes[x]) ==
        _sharedFaceMembership.end())
      _nonSharedFaceVertices.push_back(_unconstrainedNodes[x]);
  for (unsigned int x = 0; x < _constrainedNodes.size(); x++)
    if (_sharedFaceMembership.find(_constrainedNodes[x]) ==
        _sharedFaceMembership.end())
      _nonSharedFaceVertices.push_back(_constrainedNodes[x]);

  cout << " Total face unshared nodes: " << _nonSharedFaceVertices.size() << endl;

  generateNodeMasses();
}

//////////////////////////////////////////////////////////////////////
// crawl the connected vertices to sniff out the orphans
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::crawlConnectedVertices(VEC3F* vertex)
{
  vector<VEC3F*>& connections = _connectivity[vertex];
  vector<VEC3F*> toCrawl;

  // prime the stack
  for (unsigned int x = 0; x < connections.size(); x++)
    toCrawl.push_back(connections[x]);

  while (toCrawl.size() > 0)
  {
    VEC3F* connected = toCrawl.back();
    toCrawl.pop_back();
    _crawled[connected] = true;

    connections = _connectivity[connected];
    for (unsigned int x = 0; x < connections.size(); x++)
      if (_crawled.find(connections[x]) == _crawled.end())
        toCrawl.push_back(connections[x]);
  }
}

//////////////////////////////////////////////////////////////////////
// get the face-based one ring of a tet
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::getTetOneRing(int tetIndex, vector<int>& oneRing)
{
  TET* tet = _finalTets[tetIndex];

  // get the one rings of its vertices
  map<int, bool> allTets;
  for (int x = 0; x < 4; x++)
  {
    VEC3F* vertex = tet->vertices[x];
    vector<int>& membership = _tetMembership[vertex];

    for (unsigned int y = 0; y < membership.size(); y++)
      if (membership[y] != tetIndex)
        allTets[membership[y]] = true;
  }

  // see which tets share a face
  map<int, bool>::iterator i;
  for (i = allTets.begin(); i != allTets.end(); i++)
  {
    // get the index of the tet to compare to
    int otherIndex = i->first;
    TET& otherTet = *_finalTets[otherIndex];

    // if they do in fact share a face, store it
    int sharedFace = tet->sharedFace(otherTet);
    if (sharedFace >= 0)
    {
      TRIANGLE face = tet->face(sharedFace);
      for (int x = 0; x < 3; x++)
        _sharedFaceMembership[face.vertex(x)] = true;
      oneRing.push_back(otherIndex);
    }
  }
}

//////////////////////////////////////////////////////////////////////
// crawl the connected faces to sniff out the orphans
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::crawlConnectedFaces()
{
  // stomp any previous results
  _crawledTets.clear();

  // track what vertices are members of shared faces
  _sharedFaceMembership.clear();

  // build tet membership list
  _tetMembership.clear();
  for (unsigned int x = 0; x < _finalTets.size(); x++)
    for (int y = 0; y < 4; y++)
      _tetMembership[_finalTets[x]->vertices[y]].push_back(x);
 
  if (_constrainedNodes.size() == 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " No constrained nodes! Can't crawl faces! " << endl;
  }

  // get a constrained node's one ring
  vector<int>& membership = _tetMembership[_constrainedNodes[0]];
  assert(membership.size() > 0);
 
  // seed the crawl with the first tet in its one ring
  vector<int> tetStack;
  tetStack.push_back(membership[0]);
  
  // for each tet in the queue
  while (tetStack.size() > 0)
  {
    // pop from the stack
    int currentTet = tetStack.back();
    tetStack.pop_back();

    // if the tet has never been crawled before
    if (_crawledTets.find(currentTet) == _crawledTets.end())
    {
      // set it to crawled
      _crawledTets[currentTet] = true;

      // add the tets in the one ring to the queue, as long as they
      // have not been crawled either
      vector<int> oneRing;
      getTetOneRing(currentTet, oneRing);
      for (unsigned int x = 0; x < oneRing.size(); x++)
        if (_crawledTets.find(oneRing[x]) == _crawledTets.end())
          tetStack.push_back(oneRing[x]);
    }
  }
}

//////////////////////////////////////////////////////////////////////
// cull away the nodes that are not used by any tet
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::cullUnusedNodes()
{
  // track the ones to erase
  // can't erase them in place since that would invalidate
  // the iterator
  vector<map<int, VEC3F*>::iterator> toErase;

  map<int, VEC3F*>::iterator i;
  map<VEC3F*,Real>::iterator finder;
  for (i = _vertices.begin(); i != _vertices.end(); i++)
  {
    finder = _masses.find(i->second);
    if (finder == _masses.end())
      toErase.push_back(i);
  }
  for (unsigned int x = 0; x < toErase.size(); x++)
  {
    delete toErase[x]->second;
    _vertices.erase(toErase[x]);
  }
}

//////////////////////////////////////////////////////////////////////
// emit appropriate vertex of tet
//////////////////////////////////////////////////////////////////////
int ISO_STUFFER::getVertexID(TET* tet, int y)
{
  // see if it is constrained
  map<VEC3F*, int>::iterator finder;
  finder = _constrainedIDs.find(tet->vertices[y]);
  
  if (finder != _constrainedIDs.end())
    return finder->second + _unconstrainedNodes.size();
  return _unconstrainedIDs[tet->vertices[y]];
}

//////////////////////////////////////////////////////////////////////
// write an OBJ file
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::writeOBJ(const char* filename)
{
  if (_unconstrainedNodes.size() == 0 &&
      _constrainedNodes.size() == 0)
  {
    cout << " No nodes to output! Did you run generateConstrainedNodes()?" << endl;
    return;
  }

  ofstream in(filename);
  if (in.fail () )
  {
    cerr << "Can't open tetmesh file " << filename << endl;
    return;
  }

  // output the unconstrained nodes -- must write to 32 bit int so that
  // 64 bit windows doesn't think it's an int_64
  int size = _unconstrainedNodes.size();
  for (unsigned int x = 0; x < _unconstrainedNodes.size(); x++)
  {
    VEC3F* node = _unconstrainedNodes[x];
    in << "v " << (*node) << endl;
  }

  // output the constrained nodes
  size = _constrainedNodes.size();
  in << size << endl;
  for (unsigned int x = 0; x < _constrainedNodes.size(); x++)
  {
    VEC3F* node = _constrainedNodes[x];
    in << "v " << (*node) << endl;
  }

  // output the tets
  size = _finalTets.size();
  // for each tet
  for (unsigned int x = 0; x < _finalTets.size(); x++)
  {
    in << "f " << getVertexID(_finalTets[x], 0) + 1 << " "
               << getVertexID(_finalTets[x], 1) + 1 << " "
               << getVertexID(_finalTets[x], 3) + 1 << endl;

    in << "f " << getVertexID(_finalTets[x], 1) + 1 << " "
               << getVertexID(_finalTets[x], 2) + 1 << " "
               << getVertexID(_finalTets[x], 3) + 1 << endl;

    in << "f " << getVertexID(_finalTets[x], 0) + 1 << " "
               << getVertexID(_finalTets[x], 2) + 1 << " "
               << getVertexID(_finalTets[x], 1) + 1 << endl;

    in << "f " << getVertexID(_finalTets[x], 3) + 1 << " "
               << getVertexID(_finalTets[x], 2) + 1 << " "
               << getVertexID(_finalTets[x], 0) + 1 << endl;
  }

  // still need to write the face stuff here
  in.close();
}

//////////////////////////////////////////////////////////////////////
// center the mesh at 0,0,0
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::centerMesh()
{
  VEC3F sum(0,0,0);
  for (unsigned int x = 0; x < _constrainedNodes.size(); x++)
    sum = sum + *(_constrainedNodes[x]);
  for (unsigned int x = 0; x < _unconstrainedNodes.size(); x++)
    sum = sum + *(_unconstrainedNodes[x]);

  Real inv = 1.0f / (_constrainedNodes.size() + _unconstrainedNodes.size());
  sum = sum * inv;
  cout << "centering: " << sum << endl;

  for (unsigned int x = 0; x < _constrainedNodes.size(); x++)
    *(_constrainedNodes[x]) -= sum;
  for (unsigned int x = 0; x < _unconstrainedNodes.size(); x++)
    *(_unconstrainedNodes[x]) -= sum;
}

//////////////////////////////////////////////////////////////////////
// scale the mesh by a constant factor
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::scaleMesh(Real factor)
{
  for (unsigned int x = 0; x < _constrainedNodes.size(); x++)
    *(_constrainedNodes[x]) *= factor;
  for (unsigned int x = 0; x < _unconstrainedNodes.size(); x++)
    *(_unconstrainedNodes[x]) *= factor;
}

//////////////////////////////////////////////////////////////////////
// Print out a detailed timing breakdown
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::printTimingBreakdown()
{
  // create an inverse map so that it will sort by time
  map<double, string> inverseMap;
  map<string, double>::iterator forwardIter;
  for (forwardIter = _timingBreakdown.begin(); forwardIter != _timingBreakdown.end(); forwardIter++)
    inverseMap[forwardIter->second] = forwardIter->first;

  // print the map out backwards since it sorts from least to greatest
  cout << " ISO_STUFFER TIMING BREAKDOWN: " << endl;
  cout << "===============================================================================" << endl;
  map<double,string>::reverse_iterator backwardIter;
  double totalSeen = 0.0;
  for (backwardIter = inverseMap.rbegin(); backwardIter != inverseMap.rend(); backwardIter++)
  {
    string name = (*backwardIter).second + string("               ");
    name = name.substr(0,20);

    cout << "[" << (*backwardIter).first / _totalTime * 100.0 << "%\t]: "
         << name.c_str() << "\t" << (*backwardIter).first << "s total" << endl;
    totalSeen += (*backwardIter).first;
  }
  cout << "[" << (_totalTime - totalSeen) / _totalTime * 100.0 << "%\t]: "
       << "Misc. " << endl;
  cout << "===============================================================================" << endl;
  cout << " Total running time: " << _totalTime << endl;
  cout << "===============================================================================" << endl;
}

//////////////////////////////////////////////////////////////////////
// Find the biggest tet
//////////////////////////////////////////////////////////////////////
float ISO_STUFFER::maxTetVolume()
{
  if (_finalTets.size() == 0)
    return -1.0;
  float maxVolume = _finalTets[0]->volume();

  for (unsigned int x = 0; x < _finalTets.size(); x++)
    if (_finalTets[x]->volume() > maxVolume)
      maxVolume = _finalTets[x]->volume();

  return maxVolume;
}

//////////////////////////////////////////////////////////////////////
// Find the smallest tet
//////////////////////////////////////////////////////////////////////
float ISO_STUFFER::minTetVolume()
{
  if (_finalTets.size() == 0)
    return -1.0;
  float minVolume = _finalTets[0]->volume();

  for (unsigned int x = 0; x < _finalTets.size(); x++)
    if (_finalTets[x]->volume() < minVolume)
      minVolume = _finalTets[x]->volume();

  return minVolume;
}

//////////////////////////////////////////////////////////////////////
// Slice up the vertices slated for slicing
//////////////////////////////////////////////////////////////////////
int ISO_STUFFER::sliceVertices(map<VEC3F*, bool>& xCuts, map<VEC3F*, bool>& yCuts, map<VEC3F*, bool>& zCuts)
{
  int totalSliced = 0;
  map<VEC3F*, bool>::iterator i;

  // rebuild tet membership list
  _tetMembership.clear();
  for (unsigned int x = 0; x < _finalTets.size(); x++)
    for (int y = 0; y < 4; y++)
      _tetMembership[_finalTets[x]->vertices[y]].push_back(x);

  // do the y vertices first
  cout << " y cuts: " << yCuts.size() << endl;
  for (i = yCuts.begin(); i != yCuts.end(); i++)
  {
    VEC3F* original = i->first;

    // make the clone
    VEC3F* clone = new VEC3F();
    *clone = *original;
    totalSliced++;
    addVecUnguarded(clone);
    _sliced.push_back(clone);

    if (xCuts.find(original) != xCuts.end())
      xCuts[clone] = true;
    if (zCuts.find(original) != zCuts.end())
      zCuts[clone] = true;

    // get the one ring
    vector<int>& membership = _tetMembership[original];
    assert(membership.size() != 0);
    
    // divide the one ring tets into above and below
    vector<int> aboveTets;
    vector<int> belowTets;
    Real originalY = (*original)[1];
    for (unsigned int x = 0; x < membership.size(); x++)
    {
      TET& tet = *_finalTets[membership[x]];

      bool above = false;
      for (int y = 0; y < 4; y++)
      {
        VEC3F& vertex = *(tet.vertices[y]);

        // if it's not coplanar
        if (fabs(vertex[1] - originalY) > 1e-4)
        {
          if (vertex[1] - originalY > 0.0)
            above = true;
          else
            above = false;

          // found the non co-planar one, so break out
          break;
        }
      }

      if (above)
        aboveTets.push_back(membership[x]);
      else
        belowTets.push_back(membership[x]);
    }

    // reset the above ones to the clone
    for (unsigned int x = 0; x < aboveTets.size(); x++)
    {
      TET& tet = *_finalTets[aboveTets[x]];

      for (int y = 0; y < 4; y++)
        if (tet.vertices[y] == original)
          tet.vertices[y] = clone;
    }
  }

  // rebuild tet membership list
  _tetMembership.clear();
  for (unsigned int x = 0; x < _finalTets.size(); x++)
    for (int y = 0; y < 4; y++)
      _tetMembership[_finalTets[x]->vertices[y]].push_back(x);

  // do the x vertices next
  cout << " x cuts: " << xCuts.size() << endl;
  for (i = xCuts.begin(); i != xCuts.end(); i++)
  {
    VEC3F* original = i->first;

    // make the clone
    VEC3F* clone = new VEC3F();
    *clone = *original;
    totalSliced++;
    addVecUnguarded(clone);
    _sliced.push_back(clone);

    if (zCuts.find(original) != zCuts.end())
      zCuts[clone] = true;

    // get the one ring
    vector<int>& membership = _tetMembership[original];
    assert(membership.size() != 0);
    
    // divide the one ring tets into right and left
    vector<int> rightTets;
    vector<int> leftTets;
    Real originalX = (*original)[0];
    for (unsigned int x = 0; x < membership.size(); x++)
    {
      TET& tet = *_finalTets[membership[x]];

      bool right = false;
      for (int y = 0; y < 4; y++)
      {
        VEC3F& vertex = *(tet.vertices[y]);

        // if it's not coplanar
        if (fabs(vertex[0] - originalX) > 1e-4)
        {
          if (vertex[0] - originalX > 0.0)
            right = true;
          else
            right = false;

          // found the non co-planar one, so break out
          break;
        }
      }

      if (right)
        rightTets.push_back(membership[x]);
      else
        leftTets.push_back(membership[x]);
    }

    // reset the right ones to the clone
    for (unsigned int x = 0; x < rightTets.size(); x++)
    {
      TET& tet = *_finalTets[rightTets[x]];

      for (int y = 0; y < 4; y++)
        if (tet.vertices[y] == original)
          tet.vertices[y] = clone;
    }
  }

  // rebuild tet membership list
  _tetMembership.clear();
  for (unsigned int x = 0; x < _finalTets.size(); x++)
    for (int y = 0; y < 4; y++)
      _tetMembership[_finalTets[x]->vertices[y]].push_back(x);

  // do the z vertices last
  cout << " z cuts: " << zCuts.size() << endl;
  for (i = zCuts.begin(); i != zCuts.end(); i++)
  {
    VEC3F* original = i->first;

    // make the clone
    VEC3F* clone = new VEC3F();
    *clone = *original;
    totalSliced++;
    addVecUnguarded(clone);
    _sliced.push_back(clone);
    cout << " z clone: " << *clone << endl;

    // get the one ring
    vector<int>& membership = _tetMembership[original];
    
    // divide the one ring tets into far and near
    vector<int> farTets;
    vector<int> nearTets;
    Real originalZ = (*original)[2];
    for (unsigned int x = 0; x < membership.size(); x++)
    {
      TET& tet = *_finalTets[membership[x]];

      bool far = false;
      for (int y = 0; y < 4; y++)
      {
        VEC3F& vertex = *(tet.vertices[y]);

        // if it's not coplanar
        if (fabs(vertex[2] - originalZ) > 1e-4)
        {
          if (vertex[2] - originalZ > 0.0)
            far = true;
          else
            far = false;

          // found the non co-planar one, so break out
          break;
        }
      }

      if (far)
        farTets.push_back(membership[x]);
      else
        nearTets.push_back(membership[x]);
    }

    // reset the far ones to the clone
    for (unsigned int x = 0; x < farTets.size(); x++)
    {
      TET& tet = *_finalTets[farTets[x]];

      for (int y = 0; y < 4; y++)
        if (tet.vertices[y] == original)
          tet.vertices[y] = clone;
    }
  }

  return totalSliced;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void ISO_STUFFER::generateCubeTetsNonBCC(SURFACE& surface)
{
  cout << "Generating Non-BCC Cube tets ... "; flush(cout);
  // stomp all old tets;
  _finalTets.clear();

  for (int z = 0; z < _zRes; z++)
  {
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        float xCell = x + 0.5;
        float yCell = y + 0.5;
        float zCell = z + 0.5;

        // Before creating the centerVec, see if we should even bother
        if (!isOccupied(xCell, yCell, zCell, surface)) continue;

        // if this is going to create a single "hinge cube", don't bother
        int neighbors = 0;
        if (isOccupied(xCell + 1, yCell, zCell, surface)) neighbors++;
        if (isOccupied(xCell - 1, yCell, zCell, surface)) neighbors++;
        if (isOccupied(xCell, yCell + 1, zCell, surface)) neighbors++;
        if (isOccupied(xCell, yCell - 1, zCell, surface)) neighbors++;
        if (isOccupied(xCell, yCell, zCell + 1, surface)) neighbors++;
        if (isOccupied(xCell, yCell, zCell - 1, surface)) neighbors++;

        if (neighbors < 1) continue;

        VEC3F* center = centerVec(x,y,z);

        VEC3F* vert000 = corner000(x,y,z);
        VEC3F* vert100 = corner100(x,y,z);
        VEC3F* vert010 = corner010(x,y,z);
        VEC3F* vert110 = corner110(x,y,z);
        VEC3F* vert001 = corner001(x,y,z);
        VEC3F* vert101 = corner101(x,y,z);
        VEC3F* vert011 = corner011(x,y,z);
        VEC3F* vert111 = corner111(x,y,z);

        VEC3F* left  = faceCenterVec(x,y,z,-1,0,0);
        VEC3F* below = faceCenterVec(x,y,z,0,-1,0);
        VEC3F* near  = faceCenterVec(x,y,z,0,0,-1);
     
        // do all the ones to the left
        TET* newTet0 = new TET(center, vert010, left, vert011);
        _finalTets.push_back(newTet0);
        TET* newTet1 = new TET(center, vert000, left, vert010);
        _finalTets.push_back(newTet1);
        TET* newTet2 = new TET(center, vert000, vert001, left);
        _finalTets.push_back(newTet2);
        TET* newTet3 = new TET(center, left, vert001, vert011);
        _finalTets.push_back(newTet3);

        // do all the ones that are near
        TET* newTet4 = new TET(center, vert110, near, vert010);
        _finalTets.push_back(newTet4);
        TET* newTet5 = new TET(center, vert100, near, vert110);
        _finalTets.push_back(newTet5);
        TET* newTet6 = new TET(center, vert100, vert000, near);
        _finalTets.push_back(newTet6);
        TET* newTet7 = new TET(center, near, vert000, vert010);
        _finalTets.push_back(newTet7);

        // do all the ones that are below
        TET* newTet8 = new TET(center, vert100, below, vert000);
        _finalTets.push_back(newTet8);
        TET* newTet9 = new TET(center, vert101, below, vert100);
        _finalTets.push_back(newTet9);
        TET* newTet10 = new TET(center, vert101, vert001, below);
        _finalTets.push_back(newTet10);
        TET* newTet11 = new TET(center, below, vert001, vert000);
        _finalTets.push_back(newTet11);

        // check the top, right, and far ones to see if they
        // won't be creating the other sides of the cube
        {
          VEC3F* right = faceCenterVec(x,y,z,1,0,0);
          TET* newTet12 = new TET(right, vert110, center, vert111);
          TET* newTet13 = new TET(right, vert100, center, vert110);
          TET* newTet14 = new TET(right, vert100, vert101, center);
          TET* newTet15 = new TET(right, vert101, vert111, center);
          _finalTets.push_back(newTet12);
          _finalTets.push_back(newTet13);
          _finalTets.push_back(newTet14);
          _finalTets.push_back(newTet15);
        }
        {
          VEC3F* up = faceCenterVec(x,y,z,0,1,0);
          TET* newTet12 = new TET(up, vert010, center, vert011);
          TET* newTet13 = new TET(up, vert110, center, vert010);
          TET* newTet14 = new TET(up, vert110, vert111, center);
          TET* newTet15 = new TET(up, vert111, vert011, center);
          _finalTets.push_back(newTet12);
          _finalTets.push_back(newTet13);
          _finalTets.push_back(newTet14);
          _finalTets.push_back(newTet15);
        }
        {
          VEC3F* far = faceCenterVec(x,y,z,0,0,1);
          TET* newTet12 = new TET(far, vert111, center, vert011);
          TET* newTet13 = new TET(far, vert101, center, vert111);
          TET* newTet14 = new TET(far, vert101, vert001, center);
          TET* newTet15 = new TET(far, vert001, vert011, center);
          _finalTets.push_back(newTet12);
          _finalTets.push_back(newTet13);
          _finalTets.push_back(newTet14);
          _finalTets.push_back(newTet15);
        }
      }

    // output status every 10% (making sure not to mod by 0)
    if (_zRes / 10 > 0)
    {
      if (z % (int)(_zRes / 10) == 0)
      {
        cout << 100 * ((Real)z / _zRes) << "% ";
        flush(cout);
      }
    }
  }
  cout << "done." << endl;
  cout << " Generated " << _finalTets.size() << " tets " << endl;

  for (unsigned int x = 0; x < _finalTets.size(); x++)
    if (_finalTets[x]->invalid())
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " Generated an invalid tet ordering!!!! " << endl;
    }
}
