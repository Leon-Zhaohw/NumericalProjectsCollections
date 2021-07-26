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
// ISO_STUFFER.h: interface for the ISO_STUFFER class.
//
//////////////////////////////////////////////////////////////////////

#ifndef ISO_STUFFER_H
#define ISO_STUFFER_H

#include <SETTINGS.h>
#include <map>
#include <vector>
#include <fstream>
#include <VEC3.h>
#include <TET.h>
#include "OBJ.h"
#include <cstdlib>
#include "BOX.h"
#include "COMPOUND.h"
#include <TIMER.h>

#if _WIN32
#include <GL/glut.h>
#elif USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

using namespace std;

typedef pair<VEC3F*,VEC3F*> EDGE;

class ISO_STUFFER {
public:
  ISO_STUFFER(int xRes, int yRes, int zRes);
  ~ISO_STUFFER();

  enum EDGE_TYPE { UNKNOWN, RED, BLACK };
  enum CENTER_CORNER_TYPE { NEITHER, CENTER, CORNER };

  // lots of different drawing options
  void drawTets();
  void drawSurfaceTets();
  void drawSurfaceTetsZSlice(Real zSlice) { drawZSlice(zSlice, _surfaceTets); };
  void drawInsideTetsZSlice(Real zSlice) { drawZSlice(zSlice, _insideTets); };
  void drawFinalTetsZSlice(Real zSlice) { drawZSlice(zSlice, _finalTets); };
  void drawFinalTetsXSlice(Real xSlice) { drawXSlice(xSlice, _finalTets); };
  void drawInsideTets();
  void drawFinalTets();
  void drawOutsideTets();
  void drawLines();
  void drawSurfaceLines();
  void drawInsideLines();
  void drawFinalLines();
  void drawOutsideLines();
  void drawCenters();
  void drawCutPoints();
  void drawConstrainedNodes();
  void drawUnconstrainedNodes();
  void drawNonSharedFaceNodes();
  void drawSlicedNodes();

  // generate all the tets in the grid (memory efficient version)
  //
  // Seems to have problems when compiled in Visual Studio. Compile in Cygwin
  // works just fine. I'm going to go ahead and say this is Microsoft's problem,
  // not mine.
  //
  // This method is much more memory efficient than generateAllTets but is also
  // slower because it does lots of inside/outside tets to cull tets. If you
  // have tons of memory and/or absolutely have to compile in VC++, you can use
  // generateAllTets() instead.
  void generateLimitedTets(SURFACE& object);

  // insert an SURFACE instance into the grid --
  // partition the mesh into inside and outside tets
  //
  // Note this uses the "inside" and "distance" functions to determine
  // inside/outside -- not the more specialized functions used for
  // the triangle meshes
  void generateInsideTets(SURFACE& object);
  
  // generate a coarse grid of tets just based on where the triangles hashed
  // to the accel grid in the OBJ
  void generateCubeTets(SURFACE& object);

  // generate a coarse grid of tets just based on where the triangles hashed
  // to the accel grid in the OBJ
  //
  // Don't stick to the BCC lattice -- make sure even on the inside everything
  // is a grid
  void generateCubeTetsNonBCC(SURFACE& object);

  // generate a coarse grid of tets just based on where the triangles hashed
  // to the accel grid in the OBJ -- additionally, generate cloned
  // vertices where the grid should be cut
  void generateCubeTetsWithClonedVertices(OBJ& object);

  // generate the final isosurface stuffing tets
  void generateFinalTets();

  // generate the constrained and unconstrained vertex lists
  void generateConstrainedNodes(SURFACE& object);

  // generate the mass of each node
  void generateNodeMasses();

  // delete the nodes that are not used
  void cullUnusedNodes();

  // gather some statistics about the tets
  void collectTetStats();

  // output the final file
  void writeFile(const char* filename);
  void writeFileLegacy(const char* filename);

  // output an OBJ file
  void writeOBJ(const char* filename);
  int getVertexID(TET* tet, int y);

  // use an existing cache?
  bool& useExistingCaches() { return _useExistingCaches; };
  int xRes() { return _xRes; };
  int yRes() { return _yRes; };
  int zRes() { return _zRes; };

  // center the mesh at 0,0,0
  void centerMesh();

  // scale the whole mesh by a factor
  void scaleMesh(Real factor);

  // constrain extremal nodes along a Cartesian axis
  void constrainMinX() { constrainMin(0); };
  void constrainMaxX() { constrainMax(0); };
  void constrainMinY() { constrainMin(1); };
  void constrainMaxY() { constrainMax(1); };
  void constrainMinZ() { constrainMin(2); };
  void constrainMaxZ() { constrainMax(2); };
  void constrainMinMaxX() { constrainMinMax(0); };
  void constrainMinMaxY() { constrainMinMax(1); };
  void constrainMinMaxZ() { constrainMinMax(2); };
  void constrainNone();

  // print out some arbitrary timing breakdown
  void printTimingBreakdown();

  // cull away tets that are not connected to the structure as a whole
  void cullOrphanTets();

  // find the smallest tet
  float minTetVolume();
  float maxTetVolume();

  // accessors
  vector<TET*>& finalTets() { return _finalTets; };
  int vertexID(VEC3F* vertex) { return _vertexIDs[vertex]; };
  vector<VEC3F*>& unconstrainedNodes() { return _unconstrainedNodes; };
  vector<VEC3F*>& constrainedNodes() { return _constrainedNodes; };

private:
  // centered grid dimensions
  // for corner grid, add one to each dimension
  int _xRes;
  int _yRes;
  int _zRes;
  int _slabSize;

  // alpha long and short from the paper
  Real _alphaLong;
  Real _alphaShort;

  // min and max possible angles given by proof
  Real _proofMin;
  Real _proofMax;

  // maximum of all the dimensions
  int _maxRes;

  // red and black edge lengths
  Real _blackLength;
  Real _redLength;

  // distance value considered outside the grid, 
  // -2 * _maxRes
  Real _outside;

  // total corner grid cells
  int _totalCorners;

  // total body-centered grid cells
  int _totalCenters;

  // corners of the BCC gridcells
  Real* _corners;

  // centers of the BCC gridcells 
  Real* _centers;

  // all the BCC tetrahedra
  vector<TET*> _tets;

  // tets that intersect the isosurface
  vector<TET*> _surfaceTets;
  
  // edges that intersect the isosurface
  // this is a map instead of a vector to ensure there are no duplicate
  // edges. The convention is that the smaller address comes first
  // in the pair.
  map<EDGE,EDGE> _surfaceEdges;
  map<EDGE,EDGE_TYPE> _edgeTypes;
  map<VEC3F*,CENTER_CORNER_TYPE> _centerCornerType;

  // tets on the inside
  vector<TET*> _insideTets;

  // tets on the outside 
  vector<TET*> _outsideTets;

  // final isosurface stuffed tets
  vector<TET*> _finalTets;

  // cut points where edges intersect the isosurface
  map<EDGE, VEC3F*> _cutPoints;

  // the tet vertices
  //
  // They are indexed by their regular grid index. In order to
  // take into account the center grid, the "regular grid" being indexed
  // into is of size (2 * _xRes, 2 * _yRes, 2 * _zRes)
  map<int, VEC3F*> _vertices;

  // inverse mapping for unique vertex IDs, based on object addresses
  map<VEC3F*, int> _vertexIDs;

  // distances from the mesh to each vertex 
  map<int, Real> _distances;

  // stores if a vertex is inside or outside the mesh
  map<int, bool> _inside;

  // stores whether a vertex was snapped to the zero isosurface
  map<int, bool> _snapped;

  // store if a vertex is used by any tet
  map<VEC3F*, bool> _vertexUsed;

  // Connectivity information
  //
  // Given a vertex (the key), returns a vector of the other vertices
  // that neighbor that vertex
  map<VEC3F*, vector<VEC3F*> > _connectivity;

  // use existing caches?
  bool _useExistingCaches;

  // lists of constrained and unconstrained nodes
  vector<VEC3F*> _constrainedNodes;
  vector<VEC3F*> _unconstrainedNodes;
  map<VEC3F*, int> _constrainedIDs;
  map<VEC3F*, int> _unconstrainedIDs;

  // node masses
  map<VEC3F*, Real> _masses;

  // timing info
  map<string, double> _timingBreakdown;
  double _totalTime;
  
  // record connectivity information of a tet
  void recordConnectivity(TET* tet);

  // record color information of a tet
  void recordColors(TET* tet);
  void recordColors(VEC3F* v0, VEC3F* v1);

  // given the inside/outside values of the vertices near the surface,
  // flood fill the rest of the vertices with the information
  void floodInsideOutside();

  // generate distance values for nodes that are at the surface,
  // ie experience a sign flip
  void generateSurfaceDistances(SURFACE& object);

  // generate a list of all the edges belonging to tets that intersect the surface
  void generateSurfaceEdges();

  // insert an edge into the edge list, enforcing the ordering and
  // making sure not to make duplicates
  void insertEdge(VEC3F* v0, VEC3F* v1);

  // generate the cut points where the edges intersect the isosurface
  void generateCutPoints();

  // snap vertices to the cut points
  void snapToCutPoints();

  // final generation cases - these are the cases from Figure 3 of the Shewchuk 2007 paper
  void generateCase4(TET* tet, int insidePoints[4], int outsidePoints[4], VEC3F* cuts[6], EDGE edges[6]);
  void generateCase5(TET* tet, int insidePoints[4], int outsidePoints[4], VEC3F* cuts[6], EDGE edges[6]);
  void generateCase6(TET* tet, int insidePoints[4], int outsidePoints[4], VEC3F* cuts[6], EDGE edges[6]);
  void generateCase7(TET* tet, int snappedPoints[4], int insidePoints[4], int outsidePoints[4], VEC3F* cuts[6], EDGE edges[6]);
  void generateCase8(TET* tet, int insidePoints[4], int outsidePoints[4], VEC3F* cuts[6], EDGE edges[6]);
  void generateCase9(TET* tet, int insidePoints[4], int outsidePoints[4], VEC3F* cuts[6], EDGE edges[6]);
  void generateCase10(TET* tet, int insidePoints[4], int outsidePoints[4], VEC3F* cuts[6], EDGE edges[6]);
  void generateCase11(TET* tet, int insidePoints[4], int outsidePoints[4], VEC3F* cuts[6], EDGE edges[6]);

  // aux functions for the more complex stencils
  void emitCase8(VEC3F* vertices[4]);
  void emitCase10(VEC3F* vertices[4], VEC3F* cuts[3]);
  void emitCase11(VEC3F* vertices[4]);

  // aux function for lower memory footprint version
  void emitLimitedTet(VEC3F* v0, VEC3F* v1, VEC3F* v2, VEC3F* v3, 
                      bool& v0Used, bool& v1Used, bool& v2Used, bool& v3Used);

  // aux cleanup function for lower memory footprint version
  void cleanupVertex(bool used, VEC3F* vertex);

  // implements the partity rule
  // if "vertex" has an even number of components greater than "cut",
  // returns TRUE
  bool parityRule(VEC3F* vertex, VEC3F* cut);

  // generate an edge pair
  EDGE createEdge(VEC3F* v0, VEC3F* v1);

  // VEC3 creators
  VEC3F* centerVec(int x, int y, int z);
  VEC3F* cornerVec(int x, int y, int z, int cornerX, int cornerY, int cornerZ);
  VEC3F* faceCenterVec(int x, int y, int z, int cornerX, int cornerY, int cornerZ);

  // has corner vec already been created?
  bool cornerVecExists(int x, int y, int z, int cornerX, int cornerY, int cornerZ);
  
  // generate one even if an identical one exists
  VEC3F* faceCenterVecUnguarded(int x, int y, int z, int cornerX, int cornerY, int cornerZ);

  // add a general VEC3, unguarded
  void addVecUnguarded(VEC3F* vec);

  VEC3F* corner000(int x, int y, int z) { return cornerVec(x,y,z,0,0,0); };
  VEC3F* corner100(int x, int y, int z) { return cornerVec(x,y,z,1,0,0); };
  VEC3F* corner010(int x, int y, int z) { return cornerVec(x,y,z,0,1,0); };
  VEC3F* corner110(int x, int y, int z) { return cornerVec(x,y,z,1,1,0); };
  VEC3F* corner001(int x, int y, int z) { return cornerVec(x,y,z,0,0,1); };
  VEC3F* corner101(int x, int y, int z) { return cornerVec(x,y,z,1,0,1); };
  VEC3F* corner011(int x, int y, int z) { return cornerVec(x,y,z,0,1,1); };
  VEC3F* corner111(int x, int y, int z) { return cornerVec(x,y,z,1,1,1); };

  bool corner000Exists(int x, int y, int z) { return cornerVecExists(x,y,z,0,0,0); };
  bool corner100Exists(int x, int y, int z) { return cornerVecExists(x,y,z,1,0,0); };
  bool corner010Exists(int x, int y, int z) { return cornerVecExists(x,y,z,0,1,0); };
  bool corner110Exists(int x, int y, int z) { return cornerVecExists(x,y,z,1,1,0); };
  bool corner001Exists(int x, int y, int z) { return cornerVecExists(x,y,z,0,0,1); };
  bool corner101Exists(int x, int y, int z) { return cornerVecExists(x,y,z,1,0,1); };
  bool corner011Exists(int x, int y, int z) { return cornerVecExists(x,y,z,0,1,1); };
  bool corner111Exists(int x, int y, int z) { return cornerVecExists(x,y,z,1,1,1); };

  // draw the z slice of one of the tet sets
  void drawZSlice(Real zSlice, vector<TET*>& tets);

  // draw the z slice of one of the tet sets
  void drawXSlice(Real zSlice, vector<TET*>& tets);

  // generic coordinate axis constaint functions
  void constrainMin(int axis);
  void constrainMax(int axis);
  void constrainMinMax(int axis);

  map<VEC3F*, bool> _crawled;
  void crawlConnectedVertices(VEC3F* vertex);

  // test if the cube mesher should even look at this vertex
  bool isOccupied(Real x, Real y, Real z, SURFACE& surface);
  bool isOccupied(VEC3F& vec, SURFACE& surface) { return isOccupied(vec[0], vec[1], vec[2], surface); };
  bool isOccupied(VEC3F* vec, SURFACE& surface) { return isOccupied(*vec, surface); };

  // crawl the face connectivity of the mesh
  void crawlConnectedFaces();

  // needed to crawl the face connectivity
  map<VEC3F*, vector<int> > _tetMembership;
  map<int, bool> _crawledTets;
  map<VEC3F*, bool> _sharedFaceMembership;

  // which vertices are not a member of a shared face?
  vector<VEC3F*> _nonSharedFaceVertices;

  // get the face one-ring of a tet
  void getTetOneRing(int tetIndex, vector<int>& oneRing);

  // slice all the necessary vertices
  int sliceVertices(map<VEC3F*, bool>& xCuts, map<VEC3F*, bool>& yCuts, map<VEC3F*, bool>& zCuts);

  // list of vertices that got sliced
  vector<VEC3F*> _sliced;

// from here on the functions are for DEBUGGING purposes only
public:
  void drawTets0();
  void drawTets1();
  void drawTets2();
  void drawTets3();
  void drawTets4();
  void drawTets5();
  void generateType0Tets();
  void generateType1Tets();
  void generateType2Tets();
  void generateType3Tets();
  void generateType4Tets();
  void generateType5Tets();

  // generate the final isosurface stuffing tets -- reordered to be more debug-friendly
  void generateFinalTetsDebug();

  // generate all the tets in the grid (debugging version)
  void generateAllTets();

private:
  vector<TET> _tets0;
  vector<TET> _tets1;
  vector<TET> _tets2;
  vector<TET> _tets3;
  vector<TET> _tets4;
  vector<TET> _tets5;

  // cache distance values for faster debugging cycle
  void saveDistances(const char* filename);
  bool loadDistances(const char* filename);
  void saveInsideOutside(const char* filename);
  bool loadInsideOutside(const char* filename);

  bool satisfiesProof(TET& tet, int caseNumber);
};

#endif
