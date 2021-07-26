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
// SPHERE_TREE.cpp: Implementation of BD-Tree collision detection
// data structure
//
//////////////////////////////////////////////////////////////////////

#include "SPHERE_TREE.h"

MERSENNETWISTER SPHERE_TREE::_twister(123456);

SPHERE_TREE::SPHERE_TREE(OBJ* obj, int maxDepth) :
  _obj(obj), _tetMesh(NULL), _maxDepth(maxDepth), 
  _root(NULL), _rightFaces(NULL), _fundamentalTests(0),
  _unconstrained(false)
{
  _vertices = &(obj->vertices);
  _facesAtDepth.resize(maxDepth);
  initSphereTree();
  initColors();
}

SPHERE_TREE::SPHERE_TREE(TET_MESH* tetMesh, int maxDepth) :
  _obj(NULL), _tetMesh(tetMesh), _maxDepth(maxDepth), 
  _root(NULL), _rightFaces(NULL), _fundamentalTests(0),
  _unconstrained(false)
{
  _vertices = &(tetMesh->vertices());
  _facesAtDepth.resize(maxDepth);
  initSphereTree();
  initColors();
}

SPHERE_TREE::SPHERE_TREE() :
  _obj(NULL), _tetMesh(NULL), _maxDepth(12), _root(NULL), 
  _rightFaces(NULL), _fundamentalTests(0), _unconstrained(false)
{
}

SPHERE_TREE_NODE::SPHERE_TREE_NODE(SPHERE& sphere, vector<int>& faceIDs, int depth) :
  _boundingSphere(sphere), _faceIDs(faceIDs), _depth(depth)
{
}

SPHERE_TREE::~SPHERE_TREE()
{
  if (_root) delete _root;
}

/*
//////////////////////////////////////////////////////////////////////
// Initialize the tree
//////////////////////////////////////////////////////////////////////
void SPHERE_TREE::init(SUBSPACE_TET_MESH* tetMesh)
{
  _tetMesh = tetMesh;
  _vertices = &(tetMesh->vertices());
  _facesAtDepth.resize(_maxDepth);
  initSphereTree();
  initColors();
}
*/

//////////////////////////////////////////////////////////////////////
// Initialize colors for drawing
//////////////////////////////////////////////////////////////////////
void SPHERE_TREE::initColors()
{
  _colors[0][0] = 1.0f; _colors[0][1] = 0.0f; _colors[0][2] = 0.0f; _colors[0][3] = 1.0f;
  _colors[1][0] = 0.0f; _colors[1][1] = 1.0f; _colors[1][2] = 0.0f; _colors[1][3] = 1.0f;
  _colors[2][0] = 0.0f; _colors[2][1] = 0.0f; _colors[2][2] = 1.0f; _colors[2][3] = 1.0f;
  _colors[3][0] = 1.0f; _colors[3][1] = 0.0f; _colors[3][2] = 1.0f; _colors[3][3] = 1.0f;
  _colors[4][0] = 0.0f; _colors[4][1] = 1.0f; _colors[4][2] = 1.0f; _colors[4][3] = 1.0f;
  _colors[5][0] = 1.0f; _colors[5][1] = 1.0f; _colors[5][2] = 0.0f; _colors[5][3] = 1.0f;
  _colors[6][0] = 1.0f; _colors[6][1] = 1.0f; _colors[6][2] = 1.0f; _colors[6][3] = 1.0f;
  _colors[7][0] = 0.5f; _colors[7][1] = 0.5f; _colors[7][2] = 0.5f; _colors[7][3] = 1.0f;
}

//////////////////////////////////////////////////////////////////////
// compute the best bounding sphere - very expensive
//////////////////////////////////////////////////////////////////////
SPHERE SPHERE_TREE::computeBestBoundingSphere(vector<int>& faceIDs)
{
  assert(faceIDs.size() > 0);

  vector<TRIANGLE*>& faces = _tetMesh->explicitSurfaceFaces();
  VEC3F* leftVertex = faces[faceIDs[0]]->vertex(0);
  VEC3F* rightVertex = faces[faceIDs[0]]->vertex(1);
  VEC3F diff = *leftVertex - *rightVertex;
  Real maxDistance = norm(diff);

  // do an exhaustive pairwise search for the best center
  for (unsigned int x = 0; x < faceIDs.size(); x++)
  {    
    TRIANGLE& leftTriangle = *faces[faceIDs[x]];
    for (unsigned int y = 0; y < faceIDs.size(); y++)
    {
      TRIANGLE& rightTriangle = *faces[faceIDs[y]];
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
          diff = *leftTriangle.vertex(i) - *rightTriangle.vertex(j);
          if (norm(diff) > maxDistance)
          {
            leftVertex = leftTriangle.vertex(i);
            rightVertex = rightTriangle.vertex(j);
            maxDistance = norm(diff);
          }
        }
    }
  }

  // recompute the radius to be safe
  VEC3F center = (*leftVertex + *rightVertex) * 0.5;
  Real maxRadius = maxDistance * 0.5;
  for (unsigned int x = 0; x < faceIDs.size(); x++)
    for (int y = 0; y < 3; y++)
    {
      VEC3F* vertex = faces[faceIDs[x]]->vertex(y);
      VEC3F diff = *vertex - center;
      Real radius = norm(diff);
      if (radius > maxRadius)
        maxRadius = radius;
    }

  maxRadius *= 1.01;
  return SPHERE(maxRadius, center);
}

//////////////////////////////////////////////////////////////////////
// Compute the bounding sphere -- just use the layered version
// for now
//////////////////////////////////////////////////////////////////////
SPHERE SPHERE_TREE::computeBoundingSphere(vector<int>& faceIDs)
{
  // get the center
  VEC3F center;
  for (unsigned int x = 0; x < faceIDs.size(); x++)
    center += _faceCenters[faceIDs[x]];
  center *= 1.0 / faceIDs.size();

  // get the radius
  Real maxRadius = 0.0;
  vector<TRIANGLE*>& faces = _tetMesh->explicitSurfaceFaces();
  VEC3F* furthestVertex = NULL; 
  for (unsigned int x = 0; x < faceIDs.size(); x++)
    for (int y = 0; y < 3; y++)
    {
      VEC3F* vertex = faces[faceIDs[x]]->vertex(y);
      VEC3F diff = *vertex - center;
      Real radius = norm(diff);
      if (radius > maxRadius)
      {
        maxRadius = radius;
        furthestVertex = vertex;
      }
    }

  SPHERE bestSphere(maxRadius, center);
  bestSphere.radius() *= 1.01;

  if (faceIDs.size() == 1) return bestSphere;

  // try improving the sphere
  Real shrinkFactor = 0.5;
  for (int i = 0; i < 8; i++)
  {
    VEC3F newCenter = bestSphere.center() + shrinkFactor * (*furthestVertex - bestSphere.center());
    Real newRadius = 0.0;
    VEC3F* newFurthestVertex = NULL;
    for (unsigned int x = 0; x < faceIDs.size(); x++)
      for (int y = 0; y < 3; y++)
      {
        VEC3F* vertex = faces[faceIDs[x]]->vertex(y);
        VEC3F diff = *vertex - newCenter;
        Real radius = norm(diff);
        if (radius > newRadius)
        {
          newRadius = radius;
          newFurthestVertex = vertex;
        }
      }

    if (newRadius < bestSphere.radius())
    {
      bestSphere = SPHERE(newRadius, newCenter);
      shrinkFactor = 0.5;
    }
    else
      shrinkFactor *= 0.5;
  }

  // fatten it up a bit
  bestSphere.radius() *= 1.01;
  return bestSphere;
}

//////////////////////////////////////////////////////////////////////
// Build the sphere tree
//////////////////////////////////////////////////////////////////////
void SPHERE_TREE::initSphereTree()
{
  cout << "=====================================" << endl;
  cout << " Building Sphere-Tree for " 
       << _tetMesh->explicitSurfaceFaces().size() <<" faces " << endl;
  cout << "=====================================" << endl;
  _spheresAtDepth.resize(_maxDepth);

  // build a list of face centers
  vector<TRIANGLE*>& faces = _tetMesh->explicitSurfaceFaces(); 
  for (unsigned int x = 0; x < faces.size(); x++)
    _faceCenters.push_back(faces[x]->centroid());  

  // build an ID list that includes all the faces
  vector<int> faceIDs;
  for (unsigned int x = 0; x < _faceCenters.size(); x++)
    faceIDs.push_back(x);

  cout << " Building root node ... ";
  flush(cout); 

  // build the root node
  SPHERE boundingSphere = computeBoundingSphere(faceIDs);
  //SPHERE boundingSphere = computeBestBoundingSphere(faceIDs);
  _root = new SPHERE_TREE_NODE(boundingSphere, faceIDs, 0);
  buildVertexIDs(*_root);

  cout << "done." << endl;

  // build the tree recursively
  buildSphereTree(_root, 1);
  _spheresAtDepth[0]++;

  cout << " Total spheres at depths: " << endl;
  for (int x = 0; x < _maxDepth; x++)
    cout << " Depth: " << x << " Spheres: " << _spheresAtDepth[x] << " Faces: " << _facesAtDepth[x].size() << endl;

  cout << "=====================================" << endl;
}

//////////////////////////////////////////////////////////////////////
// recursively build the BD-Tree
//////////////////////////////////////////////////////////////////////
void SPHERE_TREE::buildSphereTree(SPHERE_TREE_NODE* node, int depth)
{
  if (depth >= _maxDepth) 
    return;

  // sort the vertices into octants
  vector<int> childIDs[8];
  vector<int>& faceIDs = node->faceIDs();
  VEC3F& center = node->boundingSphere().center();
  for (unsigned int x = 0; x < faceIDs.size(); x++)
  {
    VEC3F diff = _faceCenters[faceIDs[x]] - center;
    int index = 0;

    // do the slick octant-sorting bit shifting
    for (int y = 0; y < 3; y++)
	    index += (1 << y) * (diff[y] > 0);

    // store the vertex in the appropriate octant
    childIDs[index].push_back(faceIDs[x]);
  }

  // if there's only one child, don't bother -- the splitting heuristic has
  // broken down
  int totalChildren = 0;
  for (int x = 0; x < 8; x++)
    if (childIDs[x].size() > 0)
      totalChildren++;
  if (totalChildren <= 1)
    return;
  
  // build the child nodes
  vector<SPHERE_TREE_NODE>& children = node->children();
  for (int x = 0; x < 8; x++)
  {
    if (childIDs[x].size() > 0)
    {
      SPHERE boundingSphere = computeBoundingSphere(childIDs[x]);
      //SPHERE boundingSphere = computeBestBoundingSphere(childIDs[x]);
      children.push_back(SPHERE_TREE_NODE(boundingSphere, childIDs[x], depth));
      buildVertexIDs(children.back());

      _spheresAtDepth[depth]++;
      //_facesAtDepth[depth] += childIDs[x].size();
      for (unsigned int y = 0; y < childIDs[x].size(); y++)
        _facesAtDepth[depth][childIDs[x][y]] = true;
    }
  }

  // build the next level for nodes that contain more than one vertex
  for (unsigned int x = 0; x < children.size(); x++)
    if (children[x].faceIDs().size() > 1)
      buildSphereTree(&(children[x]), depth + 1);
}

//////////////////////////////////////////////////////////////////////
// draw a single sphere to GL
//////////////////////////////////////////////////////////////////////
void SPHERE_TREE::drawSphere(SPHERE_TREE_NODE& node, int depth)
{
  SPHERE& sphere = node.boundingSphere();
  VEC3F& center = sphere.center();

  glPushMatrix();
    glTranslatef(center[0], center[1], center[2]);
    glutSolidSphere(sphere.radius(), 100 / depth, 100 / depth);
  glPopMatrix();
}

//////////////////////////////////////////////////////////////////////
// draw the sphere tree to GL
//////////////////////////////////////////////////////////////////////
void SPHERE_TREE::drawSpheres(SPHERE_TREE_NODE& node, int currentDepth, int depthToDraw)
{
  // if it's this depth, draw it
  if (currentDepth == depthToDraw)
  {
    drawSphere(node, currentDepth + 1);
    return;
  }

  // else, try bubbling down a level
  vector<SPHERE_TREE_NODE>& children = node.children();
  for (unsigned int x = 0; x < children.size(); x++)
    drawSpheres(children[x], currentDepth + 1, depthToDraw);
}

//////////////////////////////////////////////////////////////////////
// draw the faces in the sphere tree to GL
//////////////////////////////////////////////////////////////////////
void SPHERE_TREE::drawFaces(SPHERE_TREE_NODE& node, int currentDepth, int depthToDraw)
{
  vector<TRIANGLE*>& faces = _tetMesh->explicitSurfaceFaces();

  // if it's this depth, draw it
  if (currentDepth == depthToDraw)
  {
    vector<int>& faceIDs = node.faceIDs();
    for (unsigned int x = 0; x < faceIDs.size(); x++)
      faces[faceIDs[x]]->draw();
    return;
  }

  // else, try bubbling down a level
  vector<SPHERE_TREE_NODE>& children = node.children();
  for (unsigned int x = 0; x < children.size(); x++)
  {
    int mod = x % 8;
    glColor4f(_colors[mod][0], _colors[mod][1], _colors[mod][2], _colors[mod][3]);
    drawFaces(children[x], currentDepth + 1, depthToDraw);
  }
}

//////////////////////////////////////////////////////////////////////
// intersect against the passed in tree and pass back a vector of
// colliding faceID pairs
//////////////////////////////////////////////////////////////////////
void SPHERE_TREE::intersect(SPHERE_TREE& rightTree, vector<pair<int, int> >& collisionPairs)
{
  // cache the other tree's surface faces
  _rightFaces = &(rightTree.tetMesh()->explicitSurfaceFaces());
  _rightTetMesh = rightTree.tetMesh();

  // clear last fundamental test result
  _fundamentalTests = 0;
  _sphereTests = 0;

  // traverse the hierarchy
  intersect(*_root, *rightTree._root, collisionPairs);

  // null out the cached faces
  _rightFaces = NULL;
  _rightTetMesh = NULL;
}

//////////////////////////////////////////////////////////////////////
// intersect against the passed in tree and pass back a vector of
// colliding faceID pairs
//////////////////////////////////////////////////////////////////////
void SPHERE_TREE::intersect(SURFACE* surface, vector<pair<SURFACE*, int> >& collisionPairs)
{
  // clear last fundamental test result
  _fundamentalTests = 0;
  _sphereTests = 0;

  intersect(*_root, surface, collisionPairs);
}

//////////////////////////////////////////////////////////////////////
// Update a sphere based on its current geometry
//////////////////////////////////////////////////////////////////////
void SPHERE_TREE::update(SPHERE_TREE_NODE* node)
{
  // if it has children, take the mean center and radius
  if (node->hasChildren())
  {
    // compute the new children first
    vector<SPHERE_TREE_NODE>& children = node->children();
    for (unsigned int x = 0; x < children.size(); x++)
      update(&(children[x]));

    // compute the new center
    VEC3F newCenter;
    for (unsigned int x = 0; x < children.size(); x++)
      newCenter += children[x].boundingSphere().center();

    newCenter *= 1.0 / children.size();

    // compute the largest radius
    Real maxRadius = 0.0;
    for (unsigned int x = 0; x < children.size(); x++)
    {
      VEC3F diff = children[x].boundingSphere().center() - newCenter;
      Real radius = norm(diff) + children[x].boundingSphere().radius();

      maxRadius = (radius > maxRadius) ? radius : maxRadius;
    }
    
    node->boundingSphere().center() = newCenter;
    node->boundingSphere().radius() = maxRadius;

    return;
  }

  // otherwise, compute the new center and radius directly

  // compute the new center
  VEC3F center;
  vector<int>& vertexIDs = node->vertexIDs();
  vector<VEC3F>& vertices = _tetMesh->vertices();
  for (unsigned int x = 0; x < vertexIDs.size(); x++)
    center += vertices[vertexIDs[x]];

  center *= 1.0 / vertexIDs.size();

  // compute the new radius
  Real maxRadius = 0.0;
  for (unsigned int x = 0; x < vertexIDs.size(); x++)
  {
    VEC3F diff = center - vertices[vertexIDs[x]];
    Real radius = norm(diff);
    maxRadius = (radius > maxRadius) ? radius : maxRadius;
  }

  node->boundingSphere().center() = center;
  node->boundingSphere().radius() = maxRadius;
}

//////////////////////////////////////////////////////////////////////
// do a brute force search of colliding triangles
//////////////////////////////////////////////////////////////////////
void SPHERE_TREE::exhaustiveIntersect(SPHERE_TREE& rightTree, vector<pair<int, int> >& collisionPairs)
{
  // cache the other tree's surface faces
  vector<TRIANGLE*> leftFaces = _tetMesh->explicitSurfaceFaces();
  vector<TRIANGLE*> rightFaces = rightTree.tetMesh()->explicitSurfaceFaces();

  // do the pairwise comparisons
  for (unsigned int x = 0; x < leftFaces.size(); x++)
    for (unsigned int y = 0; y < rightFaces.size(); y++)
      if (leftFaces[x]->intersects(*rightFaces[y]))
        collisionPairs.push_back(pair<int,int>(x,y));
}

//////////////////////////////////////////////////////////////////////
// do a brute force search of colliding triangles
//////////////////////////////////////////////////////////////////////
void SPHERE_TREE::exhaustiveIntersect(PLANE* plane, vector<pair<PLANE*, int> >& collisionPairs)
{
  // cache the other tree's surface faces
  vector<TRIANGLE*> leftFaces = _tetMesh->explicitSurfaceFaces();

  // do the pairwise comparisons
  for (unsigned int x = 0; x < leftFaces.size(); x++)
    if (overlap(*leftFaces[x], plane))
      collisionPairs.push_back(pair<PLANE*,int>(plane, x));
}

//////////////////////////////////////////////////////////////////////
// recursive intersection check
//////////////////////////////////////////////////////////////////////
void SPHERE_TREE::intersect(SPHERE_TREE_NODE& leftRoot, SPHERE_TREE_NODE& rightRoot, vector<pair<int, int> >& collisionPairs)
{
  SPHERE leftSphere = leftRoot.boundingSphere();
  SPHERE rightSphere = rightRoot.boundingSphere();

  // check for an overlap
  _sphereTests++;
  bool overlap = leftSphere.intersect(rightSphere);

  if (overlap == false) return;

  // if both have children, choose a recursion strategy
  if (leftRoot.hasChildren() && rightRoot.hasChildren())
  {
    // if left root is at higher or equal depth, split it
    if (leftRoot.depth() <= rightRoot.depth())
    {
      assert(leftRoot.hasChildren());
      vector<SPHERE_TREE_NODE>& children = leftRoot.children();
      for (unsigned int x = 0; x < children.size(); x++)
        intersect(children[x], rightRoot, collisionPairs);
    }
    // else split the right one
    else
    {
      assert(rightRoot.hasChildren());
      vector<SPHERE_TREE_NODE>& children = rightRoot.children();
      for (unsigned int x = 0; x < children.size(); x++)
        intersect(leftRoot, children[x], collisionPairs);
    }

    // all done -- skip the fundamental tests
    return;
  }

  // if just one has children, split it
  if (leftRoot.hasChildren())
  {
    vector<SPHERE_TREE_NODE>& children = leftRoot.children();
    for (unsigned int x = 0; x < children.size(); x++)
      intersect(children[x], rightRoot, collisionPairs);
    
    // all done -- skip the fundamental tests
    return;
  }
  if (rightRoot.hasChildren())
  {
    vector<SPHERE_TREE_NODE>& children = rightRoot.children();
    for (unsigned int x = 0; x < children.size(); x++)
      intersect(leftRoot, children[x], collisionPairs);
    
    // all done -- skip the fundamental tests
    return;
  }

  // make sure both are leaves
  assert(!leftRoot.hasChildren());
  assert(!rightRoot.hasChildren());

  // cache all faces and face IDs in preparation of fundamental tests
  vector<int>& leftFaceIDs = leftRoot.faceIDs();
  vector<int>& rightFaceIDs = rightRoot.faceIDs();
  vector<TRIANGLE*>& leftFaces = _tetMesh->explicitSurfaceFaces();
  assert(_rightFaces != NULL);
  vector<TRIANGLE*>& rightFaces = *_rightFaces;

  // do all pairwise fundamental tests
  for (unsigned int x = 0; x < leftFaceIDs.size(); x++)
  {
    TRIANGLE* leftTriangle = leftFaces[leftFaceIDs[x]];
    for (unsigned int y = 0; y < rightFaceIDs.size(); y++)
    {
      TRIANGLE* rightTriangle = rightFaces[rightFaceIDs[y]];

      if (leftTriangle->intersects(*rightTriangle))
        collisionPairs.push_back(pair<int, int>(leftFaceIDs[x], rightFaceIDs[y]));

      _fundamentalTests++;
    }
  }
  return;
}

//////////////////////////////////////////////////////////////////////
// recursive intersection check
//////////////////////////////////////////////////////////////////////
void SPHERE_TREE::intersect(SPHERE_TREE_NODE& leftRoot, SURFACE* surface, vector<pair<SURFACE*, int> >& collisionPairs)
{
  SPHERE leftSphere = leftRoot.boundingSphere();

  // check for an overlap
  _sphereTests++;
  bool overlapping = false;

  if (surface->type().compare("PLANE") == 0)
    overlapping = overlap(leftSphere, (PLANE*)surface);
  else if (surface->type().compare("SPHERE") == 0)
  {
    SPHERE* rightSphere = (SPHERE*)surface;
    overlapping = leftSphere.intersect(*rightSphere);
  }
  else if (surface->type().compare("BOX") == 0)
  {
    overlapping = overlap(leftSphere, (BOX*)surface);
  }
  else
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " SURFACE type unknown! " << endl;
  }

  if (overlapping == false) return;

  // if there are still children, recurse
  if (leftRoot.hasChildren())
  {
    vector<SPHERE_TREE_NODE>& children = leftRoot.children();
    for (unsigned int x = 0; x < children.size(); x++)
      intersect(children[x], surface, collisionPairs);
    
    // all done -- skip the fundamental tests
    return;
  }

  // make sure both are leaves
  assert(!leftRoot.hasChildren());

  // cache all faces and face IDs in preparation of fundamental tests
  vector<int>& leftFaceIDs = leftRoot.faceIDs();
  vector<TRIANGLE*>& leftFaces = _tetMesh->explicitSurfaceFaces();

  // do all pairwise fundamental tests
  for (unsigned int x = 0; x < leftFaceIDs.size(); x++)
  {
    TRIANGLE* leftTriangle = leftFaces[leftFaceIDs[x]];

    bool overlapping = false;
    if (surface->type().compare("PLANE") == 0)
      overlapping = overlap(*leftTriangle, (PLANE*)surface);
    else if (surface->type().compare("SPHERE") == 0)
      overlapping = overlap(*leftTriangle, (SPHERE*)surface);
    else if (surface->type().compare("BOX") == 0)
      //overlapping = overlap(*leftTriangle, (BOX*)surface);
      overlapping = ((BOX*)surface)->overlap(*leftTriangle);
    else
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " SURFACE type unknown! " << endl;
    }

    if (overlapping)
      collisionPairs.push_back(pair<SURFACE*, int>(surface, leftFaceIDs[x]));
    _fundamentalTests++;
  }
  return;
}

//////////////////////////////////////////////////////////////////////
// DEBUG: do some analysis on the collison pairs
//////////////////////////////////////////////////////////////////////
vector<pair<int, int> > SPHERE_TREE::analyzeCollisions(SPHERE_TREE& rightTree, vector<pair<int, int> >& collisionPairs)
{
  vector<pair<int, int> > finals;
  cout << " Original collision pairs: " << collisionPairs.size() << endl;

  for (unsigned int x = 0; x < collisionPairs.size(); x++)
  {
    vector<int> leftFace;
    leftFace.push_back(collisionPairs[x].first);
    vector<int> rightFace;
    rightFace.push_back(collisionPairs[x].second);

    SPHERE leftSphere = computeBoundingSphere(leftFace);
    SPHERE rightSphere = rightTree.computeBoundingSphere(rightFace);

    bool overlap = leftSphere.intersect(rightSphere);
    if (!overlap)
    {
      cout << " Pair missed! " << endl;

      VEC3F diff = leftSphere.center() - rightSphere.center();
      Real combined = leftSphere.radius() + rightSphere.radius();
      cout << " diff: " << norm(diff) << endl;
      cout << " combined: " << combined << endl;

      vector<TRIANGLE*> leftTriangles = this->tetMesh()->explicitSurfaceFaces();
      vector<TRIANGLE*> rightTriangles = rightTree.tetMesh()->explicitSurfaceFaces();

      TRIANGLE& leftTri = *leftTriangles[collisionPairs[x].first];
      TRIANGLE& rightTri = *rightTriangles[collisionPairs[x].second];

      cout << " left centroid: " << leftTri.centroid() << endl;
      cout << " left sphere center: " << leftSphere.center() << endl;

      cout << " right centroid: " << rightTri.centroid() << endl;
      cout << " right sphere center: " << rightSphere.center() << endl;

      finals.push_back(collisionPairs[x]);
    }
  }
  return finals;
}

//////////////////////////////////////////////////////////////////////
// build vertex IDs
//////////////////////////////////////////////////////////////////////
void SPHERE_TREE::buildVertexIDs(SPHERE_TREE_NODE& node)
{
  vector<TRIANGLE*>& faces = _tetMesh->explicitSurfaceFaces();
  vector<int>& faceIDs = node.faceIDs();

  // hash all the vertices to eliminate duplicates
  map<int, bool> vertexHash;
  for (unsigned int x = 0; x < faceIDs.size(); x++)
  {
    TRIANGLE& triangle = *faces[faceIDs[x]];
    for (int y = 0; y < 3; y++)
    {
      int vertexID = _tetMesh->vertexID(triangle.vertex(y));
      vertexHash[vertexID] = true;
    }
  }

  // get a list of the vertices
  vector<int>& vertexIDs = node.vertexIDs();
  vertexIDs.clear();
  map<int, bool>::iterator iter;
  for (iter = vertexHash.begin(); iter != vertexHash.end(); iter++)
    vertexIDs.push_back(iter->first);
}

//////////////////////////////////////////////////////////////////////
// do an AABB -sphere overlap test
//
// A Simple Method for Box-Sphere Intersection Testing
// by Jim Arvo
// from "Graphics Gems", Academic Press, 1990
//////////////////////////////////////////////////////////////////////
bool SPHERE_TREE::overlap(SPHERE& sphere, BOX* box)
{
  // solid sphere and box
  float dmin = 0;
  Real r2 = sphere.radius() * sphere.radius();
  VEC3F& C = sphere.center();
  VEC3F Bmin = box->boxMins();
  VEC3F Bmax = box->boxMaxs();
  for(int i = 0; i < 3; i++ ) {
    if( C[i] < Bmin[i] ) 
    {
      float diff = C[i] - Bmin[i]; 
      dmin += diff * diff;
    }
    else if( C[i] > Bmax[i] )
    { 
      float diff = C[i] - Bmax[i]; 
      dmin += diff * diff;
    }
  }
  if( dmin <= r2 ) return true;

  // hollow sphere and box
  /*
  Real dmin = 0;
  Real dmax = 0;
  Real r2 = sphere.radius() * sphere.radius();
  bool face = false;

  VEC3F& C = sphere.center();
  VEC3F Bmin = box->boxMins();
  VEC3F Bmax = box->boxMaxs();

  for (int i = 0; i < 3; i++) 
  {
    float a = C[i] - Bmin[i];
    float b = C[i] - Bmax[i];
    a = a * a;
    b = b * b;

    dmax += (a > b) ? a : b;
    if (C[i] < Bmin[i]) {
      face = true;
      dmin += a;
    }
    else if (C[i] > Bmax[i]) {
      face = true;
      dmin += b;
    }
    else {
      float pairMin = (a < b) ? a : b;
      if( pairMin <= r2 )
        face = true;
    }
  }
  if (face && ( dmin <= r2 ) && ( r2 <= dmax)) return true;
  */

  return false;
}

//////////////////////////////////////////////////////////////////////
// do a plane-sphere overlap test
//////////////////////////////////////////////////////////////////////
bool SPHERE_TREE::overlap(SPHERE& sphere, PLANE* plane)
{
  // get the distance to the plane
  VEC3F sphereCenter = sphere.center();
  VEC3F planePoint = plane->point();
  VEC3F planeNormal = plane->normal();

  VEC3F diff = sphereCenter - planePoint;

  Real distance = diff * planeNormal;
  return (distance <= sphere.radius());
}

//////////////////////////////////////////////////////////////////////
// do a plane-triangle overlap test
//////////////////////////////////////////////////////////////////////
bool SPHERE_TREE::overlap(TRIANGLE& triangle, PLANE* plane)
{
  VEC3F point = plane->point();
  Real dots[3];
  for (int x = 0; x < 3; x++)
  {
    VEC3F diff = (*triangle.vertex(x)) - point;
    dots[x] = plane->normal() * diff;
  }

  /*
  // if they're all on the same side of the plane, it does not overlap
  if (dots[0] * dots[1] > 0.0 &&
      dots[1] * dots[2] > 0.0 &&
      dots[0] * dots[2] > 0.0)
    return false;

  return true;
  */

  for (int x = 0; x < 3; x++)
    if (dots[x] < plane->stickiness())
      return true;

  return false;
}

//////////////////////////////////////////////////////////////////////
// do a plane-triangle overlap test
//////////////////////////////////////////////////////////////////////
bool SPHERE_TREE::overlap(TRIANGLE& triangle, SPHERE* sphere)
{
  VEC3F center = sphere->center();
  Real distances[3];
  for (int x = 0; x < 3; x++)
  {
    VEC3F diff = (*triangle.vertex(x)) - center;
    distances[x] = sphere->radius() - sqrt(diff * diff);

    if (distances[x] > 0.0) return true;
  }

  // if they're all on the same side of the sphere, it does not overlap
  if (distances[0] * distances[1] > 0.0 &&
      distances[1] * distances[2] > 0.0 &&
      distances[0] * distances[2] > 0.0)
    return false;

  return true;
}
