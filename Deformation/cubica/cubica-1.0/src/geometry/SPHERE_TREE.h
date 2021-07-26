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
// SPHERE_TREE.h: interface for the SPHERE_TREE class.
//
//////////////////////////////////////////////////////////////////////

#ifndef SPHERE_TREE_H
#define SPHERE_TREE_H

#include <TRIANGLE.h>
#include <OBJ.h>
#include <TET_MESH.h>
//#include <UNCONSTRAINED_SUBSPACE_TET_MESH.h>
#include <SPHERE.h>
#include <PLANE.h>
#include <MERSENNETWISTER.h>
#include <BOX.h>

#if _WIN32
#include <gl/glut.h>
#elif USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

//////////////////////////////////////////////////////////////////////
// actual nodes in the BD-Tree
//////////////////////////////////////////////////////////////////////
class SPHERE_TREE_NODE {
public:
  SPHERE_TREE_NODE(SPHERE& sphere, vector<int>& faceIDs, int depth);
  SPHERE& boundingSphere()          { return _boundingSphere; };
  vector<SPHERE_TREE_NODE>& children()  { return _children; };
  bool hasChildren()                { return (_children.size() != 0); };
  vector<int>& faceIDs()            { return _faceIDs; };
  vector<int>& vertexIDs()          { return _vertexIDs; };
  int& depth()                      { return _depth; };

private:
  SPHERE _boundingSphere;
  vector<SPHERE_TREE_NODE> _children;
  vector<int> _faceIDs;
  vector<int> _vertexIDs;
  int _depth;
};

//////////////////////////////////////////////////////////////////////
// high-level BD-Tree class
//////////////////////////////////////////////////////////////////////
class SPHERE_TREE {

public:
  SPHERE_TREE();
  SPHERE_TREE(OBJ* obj, int maxDepth = 12);
  SPHERE_TREE(TET_MESH* tetMesh, int maxDepth = 12);
  ~SPHERE_TREE();

  //void init(SUBSPACE_TET_MESH* tetMesh);

  void drawSpheres(int depth)  { drawSpheres(*_root, 0, depth); };
  void drawFaces(int depth)    { drawFaces(*_root, 0, depth); };
  const int maxDepth() const   { return _maxDepth; };
  int fundamentalTests()       { return _fundamentalTests; };
  int sphereTests()       { return _sphereTests; };
  bool& unconstrained()   { return _unconstrained; };
  TET_MESH* tetMesh() { return _tetMesh; };
  static MERSENNETWISTER& twister() { return _twister; };

  // intersect against the passed in tree and pass back a vector of
  // colliding faceID pairs
  void intersect(SPHERE_TREE& rightTree, vector<pair<int, int> >& collisionPairs);
 
  // intersect against a surface, and pass back a vector of surface, triangle ID
  // pairs 
  void intersect(SURFACE* surface, vector<pair<SURFACE*, int> >& collisionPairs);

  // do a brute force search of colliding triangles
  void exhaustiveIntersect(SPHERE_TREE& rightTree, vector<pair<int, int> >& collisionPairs);
  void exhaustiveIntersect(PLANE* plane, vector<pair<PLANE*, int> >& collisionPairs);

  // DEBUG: do some analysis on the collison pairs
  // pass back the vector of missed pairs
  vector<pair<int, int> > analyzeCollisions(SPHERE_TREE& rightTree, vector<pair<int, int> >& collisionPairs);

  // update the spheres based on new geometry
  void update() { update(_root); };

private:
  // actual geometry being collided
  OBJ* _obj;

  // tet mesh being collided
  TET_MESH* _tetMesh;

  // how deep to construct the tree hierarchy?
  int _maxDepth;

  // the actual tree
  SPHERE_TREE_NODE* _root;

  // keep a pointer to the vertices around
  vector<VEC3F>* _vertices;

  // vector of face centers
  vector<VEC3F> _faceCenters;

  // count the nodes at different levels
  vector<int> _spheresAtDepth;
  
  // hash table of what vertices belong to which triangles
  map<VEC3F*, vector<int> > _faceHash;

  // colors for the different partitions
  float _colors[8][4];

  // vector of triangles being intersected against --
  // store globally to avoid having the pass the right SPHERE_TREE down the recursion as well
  vector<TRIANGLE*>* _rightFaces;

  // tet mesh being intersected against --
  // store globally to avoid having the pass the right SPHERE_TREE down the recursion as well
  TET_MESH* _rightTetMesh;

  // track how many fundamental tests were done last time
  int _fundamentalTests;

  // total faces at each level
  vector<map<int, bool> > _facesAtDepth;

  // track how many sphere tests were done last time
  int _sphereTests;

  // mesh is unconstrained?
  bool _unconstrained;

  // random number generator for multi basis
  static MERSENNETWISTER _twister;

  // Compute the bounding sphere -- just use the layered version for now
  SPHERE computeBoundingSphere(vector<int>& faceIDs);

  // compute the best bounding sphere - very expensive
  SPHERE computeBestBoundingSphere(vector<int>& faceIDs);
 
  // topmost call for recursively constructing the BD-Tree
  void initSphereTree();

  // recursively build the BD-Tree
  void buildSphereTree(SPHERE_TREE_NODE* node, int depth);

  // draw a single sphere
  void drawSphere(SPHERE_TREE_NODE& node, int depth);
  
  // recursively draw the spheres at a desired depth
  void drawSpheres(SPHERE_TREE_NODE& node, int currentDepth, int depthToDraw);

  // recursively draw the faces at a desired depth
  void drawFaces(SPHERE_TREE_NODE& node, int currentDepth, int depthToDraw);

  // initialize colors
  void initColors();

  // recursive intersection check
  void intersect(SPHERE_TREE_NODE& leftRoot, SPHERE_TREE_NODE& rightRoot, vector<pair<int, int> >& collisionPairs);
  
  // recursive intersection check
  void intersect(SPHERE_TREE_NODE& leftRoot, SURFACE* surface, vector<pair<SURFACE*, int> >& collisionPairs);
  
  // build vertex IDs
  void buildVertexIDs(SPHERE_TREE_NODE& node);

  // do a plane-sphere overlap test
  bool overlap(SPHERE& sphere, PLANE* plane);

  // do a box-sphere overlap test
  bool overlap(SPHERE& sphere, BOX* box);

  // do a plane-triangle overlap test
  bool overlap(TRIANGLE& triangle, PLANE* plane);

  // do a sphere-triangle overlap test
  bool overlap(TRIANGLE& triangle, SPHERE* sphere);

  // recursively update the sphere tree
  void update(SPHERE_TREE_NODE* node);
};

#endif
