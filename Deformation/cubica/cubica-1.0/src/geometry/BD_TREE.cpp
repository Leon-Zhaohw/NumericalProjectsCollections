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
// BD_TREE.cpp: Implementation of BD-Tree collision detection
// data structure
//
//////////////////////////////////////////////////////////////////////

#include "BD_TREE.h"

MERSENNETWISTER BD_TREE::_twister(123456);

BD_TREE::BD_TREE(OBJ* obj, int maxDepth) :
  _obj(obj), _tetMesh(NULL), _maxDepth(maxDepth), 
  _root(NULL), _rightFaces(NULL), _fundamentalTests(0),
  _unconstrained(false)
{
  _vertices = &(obj->vertices);
  _facesAtDepth.resize(maxDepth);
  initSphereTree();
  initColors();
}

BD_TREE::BD_TREE(SUBSPACE_TET_MESH* tetMesh, int maxDepth) :
  _obj(NULL), _tetMesh(tetMesh), _maxDepth(maxDepth), 
  _root(NULL), _rightFaces(NULL), _fundamentalTests(0),
  _unconstrained(false)
{
  _vertices = &(tetMesh->vertices());
  _facesAtDepth.resize(maxDepth);
  initSphereTree();
  initColors();
}

BD_TREE::BD_TREE() :
  _obj(NULL), _tetMesh(NULL), _maxDepth(12), _root(NULL), 
  _rightFaces(NULL), _fundamentalTests(0), _unconstrained(false)
{
}

BD_TREE_NODE::BD_TREE_NODE(SPHERE& sphere, vector<int>& faceIDs, int depth) :
  _boundingSphere(sphere), _faceIDs(faceIDs), _depth(depth)
{
}

BD_TREE::~BD_TREE()
{
  if (_root) delete _root;
}

//////////////////////////////////////////////////////////////////////
// Initialize the tree
//////////////////////////////////////////////////////////////////////
void BD_TREE::init(SUBSPACE_TET_MESH* tetMesh)
{
  _tetMesh = tetMesh;
  _vertices = &(tetMesh->vertices());
  _facesAtDepth.resize(_maxDepth);
  initSphereTree();
  initColors();
}

//////////////////////////////////////////////////////////////////////
// Initialize colors for drawing
//////////////////////////////////////////////////////////////////////
void BD_TREE::initColors()
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
SPHERE BD_TREE::computeBestBoundingSphere(vector<int>& faceIDs)
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
SPHERE BD_TREE::computeBoundingSphere(vector<int>& faceIDs)
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
void BD_TREE::initSphereTree()
{
  cout << "=====================================" << endl;
  cout << " Building BD-Tree for " 
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
  _root = new BD_TREE_NODE(boundingSphere, faceIDs, 0);
  buildVertexIDs(*_root);
  buildVertexBasis(*_root);
  buildCenterBasis(*_root);
  buildRadiusBasis(*_root);
  buildRadiusMultiBasis(*_root);

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
void BD_TREE::buildSphereTree(BD_TREE_NODE* node, int depth)
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
  vector<BD_TREE_NODE>& children = node->children();
  for (int x = 0; x < 8; x++)
  {
    if (childIDs[x].size() > 0)
    {
      SPHERE boundingSphere = computeBoundingSphere(childIDs[x]);
      children.push_back(BD_TREE_NODE(boundingSphere, childIDs[x], depth));
      buildVertexIDs(children.back());
      buildVertexBasis(children.back());
      buildCenterBasis(children.back());
      buildRadiusBasis(children.back());
      buildRadiusMultiBasis(children.back());

      _spheresAtDepth[depth]++;
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
void BD_TREE::drawSphere(BD_TREE_NODE& node, int depth)
{
  SPHERE& sphere = node.boundingSphere();
  VEC3F& center = sphere.center();

  VECTOR& q = _tetMesh->q();
  VEC3F centerPrime = node.centerBasis() * q;
  VEC3F finalCenter = center + centerPrime;

  VECTOR qAbs = _tetMesh->q();
  qAbs.fabs();
  Real radiusDelta = node.radiusBasis() * qAbs;

  // if it's unconstrained, apply a rigid transform
  if (_unconstrained)
  {
    UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_tetMesh;
    VEC3F translation = mesh->rigidTranslation();
    MATRIX3 rotation = mesh->rotationQuaternion().toExplicitMatrix3x3();

    finalCenter = rotation * finalCenter + translation;
  }

  glPushMatrix();
    glTranslatef(finalCenter[0], finalCenter[1], finalCenter[2]);
    glutSolidSphere(sphere.radius() + radiusDelta, 100 / depth, 100 / depth);
  glPopMatrix();
}

//////////////////////////////////////////////////////////////////////
// draw the sphere tree to GL
//////////////////////////////////////////////////////////////////////
void BD_TREE::drawSpheres(BD_TREE_NODE& node, int currentDepth, int depthToDraw)
{
  // if it's this depth, draw it
  if (currentDepth == depthToDraw)
  {
    drawSphere(node, currentDepth + 1);
    return;
  }

  // else, try bubbling down a level
  vector<BD_TREE_NODE>& children = node.children();
  for (unsigned int x = 0; x < children.size(); x++)
    drawSpheres(children[x], currentDepth + 1, depthToDraw);
}

//////////////////////////////////////////////////////////////////////
// draw the faces in the sphere tree to GL
//////////////////////////////////////////////////////////////////////
void BD_TREE::drawFaces(BD_TREE_NODE& node, int currentDepth, int depthToDraw)
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
  vector<BD_TREE_NODE>& children = node.children();
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
void BD_TREE::intersect(BD_TREE& rightTree, vector<pair<int, int> >& collisionPairs)
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
void BD_TREE::intersect(SURFACE* surface, vector<pair<SURFACE*, int> >& collisionPairs)
{
  // clear last fundamental test result
  _fundamentalTests = 0;
  _sphereTests = 0;

  // traverse the hierarchy
  if (_unconstrained)
    intersectUnconstrained(*_root, surface, collisionPairs);
  else
    intersect(*_root, surface, collisionPairs);
}

//////////////////////////////////////////////////////////////////////
// intersect against the passed in tree and pass back a vector of
// colliding faceID pairs
//////////////////////////////////////////////////////////////////////
void BD_TREE::intersectUnconstrained(SURFACE* surface, vector<pair<SURFACE*, int> >& collisionPairs)
{
  // clear last fundamental test result
  _fundamentalTests = 0;
  _sphereTests = 0;

  // traverse the hierarchy
  intersectUnconstrained(*_root, surface, collisionPairs);
}

//////////////////////////////////////////////////////////////////////
// do a brute force search of colliding triangles
//////////////////////////////////////////////////////////////////////
void BD_TREE::exhaustiveIntersect(BD_TREE& rightTree, vector<pair<int, int> >& collisionPairs)
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
void BD_TREE::exhaustiveIntersect(PLANE* plane, vector<pair<PLANE*, int> >& collisionPairs)
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
void BD_TREE::intersect(BD_TREE_NODE& leftRoot, BD_TREE_NODE& rightRoot, vector<pair<int, int> >& collisionPairs)
{
  SPHERE leftSphere = leftRoot.boundingSphere();
  SPHERE rightSphere = rightRoot.boundingSphere();

  VECTOR qLeft = _tetMesh->q();
  VECTOR qRight = _rightTetMesh->q();

  // update the sphere
  VEC3F leftCenterPrime = leftRoot.centerBasis() * qLeft;
  VEC3F rightCenterPrime = rightRoot.centerBasis() * qRight;
  leftSphere.center() += leftCenterPrime;
  rightSphere.center() += rightCenterPrime;

  // take the absolute value for radius computation
  qLeft.fabs();
  qRight.fabs();
  Real leftRadiusPrime = leftRoot.radiusBasis() * qLeft;
  Real rightRadiusPrime = rightRoot.radiusBasis() * qRight;
  leftSphere.radius() += leftRadiusPrime;
  rightSphere.radius() += rightRadiusPrime;

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
      vector<BD_TREE_NODE>& children = leftRoot.children();
      for (unsigned int x = 0; x < children.size(); x++)
        intersect(children[x], rightRoot, collisionPairs);
    }
    // else split the right one
    else
    {
      assert(rightRoot.hasChildren());
      vector<BD_TREE_NODE>& children = rightRoot.children();
      for (unsigned int x = 0; x < children.size(); x++)
        intersect(leftRoot, children[x], collisionPairs);
    }

    // all done -- skip the fundamental tests
    return;
  }

  // if just one has children, split it
  if (leftRoot.hasChildren())
  {
    vector<BD_TREE_NODE>& children = leftRoot.children();
    for (unsigned int x = 0; x < children.size(); x++)
      intersect(children[x], rightRoot, collisionPairs);
    
    // all done -- skip the fundamental tests
    return;
  }
  if (rightRoot.hasChildren())
  {
    vector<BD_TREE_NODE>& children = rightRoot.children();
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
void BD_TREE::intersectUnconstrained(BD_TREE_NODE& leftRoot, SURFACE* surface, vector<pair<SURFACE*, int> >& collisionPairs)
{
  UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_tetMesh;
  VEC3F translation = mesh->rigidTranslation();
  MATRIX3 rotation = mesh->rotationQuaternion().toExplicitMatrix3x3();

  SPHERE leftSphere = leftRoot.boundingSphere();
  VECTOR qLeft = _tetMesh->q();

  // update the sphere
  VEC3F leftCenterPrime = leftRoot.centerBasis() * qLeft;
  leftSphere.center() += leftCenterPrime;

  // take the absolute value for radius computation
  qLeft.fabs();
  Real leftRadiusPrime = leftRoot.radiusBasis() * qLeft;
  leftSphere.radius() += leftRadiusPrime;

  // take into account the rigid transform
  leftSphere.center() = rotation * leftSphere.center() + translation;

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
    overlapping = overlap(leftSphere, (BOX*)surface);
  else
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " SURFACE type unknown! " << endl;
  }

  if (overlapping == false) return;

  // if there are still children, recurse
  if (leftRoot.hasChildren())
  {
    vector<BD_TREE_NODE>& children = leftRoot.children();
    for (unsigned int x = 0; x < children.size(); x++)
      intersectUnconstrained(children[x], surface, collisionPairs);
    
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

    // make a transformed version of the triangle -- applying it to the triangle
    // is better than to the surface, since it avoids the need to write a new
    // transform for each new surface type
    //
    // it will be slower however, since the transform has to be applied to each
    // vertex, as opposed to once to the surface
    vector<VEC3F> vertices;
    for (int y = 0; y < 3; y++)
    {
      vertices.push_back(*leftTriangle->vertex(y));
      vertices[y] = rotation * vertices[y] + translation;
    }
    TRIANGLE transformedTriangle(&vertices[0], &vertices[1], &vertices[2]);

    bool overlapping = false;
    if (surface->type().compare("PLANE") == 0)
      overlapping = overlap(transformedTriangle, (PLANE*)surface);
    else if (surface->type().compare("SPHERE") == 0)
      overlapping = overlap(transformedTriangle, (SPHERE*)surface);
    else if (surface->type().compare("BOX") == 0)
      overlapping = ((BOX*)surface)->overlap(transformedTriangle);
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
// recursive intersection check
//////////////////////////////////////////////////////////////////////
void BD_TREE::intersect(BD_TREE_NODE& leftRoot, SURFACE* surface, vector<pair<SURFACE*, int> >& collisionPairs)
{
  SPHERE leftSphere = leftRoot.boundingSphere();
  VECTOR qLeft = _tetMesh->q();

  // update the sphere
  VEC3F leftCenterPrime = leftRoot.centerBasis() * qLeft;
  leftSphere.center() += leftCenterPrime;

  // take the absolute value for radius computation
  qLeft.fabs();
  Real leftRadiusPrime = leftRoot.radiusBasis() * qLeft;
  leftSphere.radius() += leftRadiusPrime;

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
    overlapping = overlap(leftSphere, (BOX*)surface);
  else
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " SURFACE type unknown! " << endl;
  }

  if (overlapping == false) return;

  // if there are still children, recurse
  if (leftRoot.hasChildren())
  {
    vector<BD_TREE_NODE>& children = leftRoot.children();
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
vector<pair<int, int> > BD_TREE::analyzeCollisions(BD_TREE& rightTree, vector<pair<int, int> >& collisionPairs)
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
// build the face basis for a node
//////////////////////////////////////////////////////////////////////
void BD_TREE::buildVertexBasis(BD_TREE_NODE& node)
{
  MATRIX& U = _tetMesh->U();
  vector<int>& vertexIDs = node.vertexIDs();

  MATRIX& vertexBasis = node.vertexBasis();
  vertexBasis.resizeAndWipe(vertexIDs.size() * 3, _tetMesh->rank());

  for (unsigned int x = 0; x < vertexIDs.size(); x++)
  {
    int index = 3 * vertexIDs[x];

    // if the vertex is unconstrained, and therefore has a basis row
    if (index < U.rows())
    {
      SUBMATRIX submatrix(U, index, 3);
      submatrix.copiesInto(vertexBasis, 3 * x);
    }
    else
    {
      MATRIX zeros(3, U.cols());
      zeros.copiesInto(vertexBasis, 3 * x);
    }
  }
}

//////////////////////////////////////////////////////////////////////
// build the center basis for a node
//////////////////////////////////////////////////////////////////////
void BD_TREE::buildCenterBasis(BD_TREE_NODE& node)
{
  MATRIX& vertexBasis = node.vertexBasis();
  MATRIX& centerBasis = node.centerBasis();

  // resize the center basis
  int rank = _tetMesh->rank();
  centerBasis.resizeAndWipe(3, rank);

  // make the averaging basis
  MATRIX meanBasis(3, vertexBasis.rows());
  for (int x = 0; x < meanBasis.cols() / 3; x++)
  {
    meanBasis(0, 3 * x) = 1;
    meanBasis(1, 3 * x + 1) = 1;
    meanBasis(2, 3 * x + 2) = 1;
  }
  meanBasis *= 1.0 / node.vertexIDs().size();

  // compute the final center basis
  centerBasis = meanBasis * vertexBasis;
}

//////////////////////////////////////////////////////////////////////
// build vertex IDs
//////////////////////////////////////////////////////////////////////
void BD_TREE::buildVertexIDs(BD_TREE_NODE& node)
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
// build the radius basis for a node
//////////////////////////////////////////////////////////////////////
void BD_TREE::buildRadiusBasis(BD_TREE_NODE& node)
{
  MATRIX& vertexBasis = node.vertexBasis();
  VECTOR& radiusBasis = node.radiusBasis();
  MATRIX& radiusBasis3 = node.radiusBasis3();
  MATRIX& centerBasis = node.centerBasis();
  MATRIX& axisBasis = node.axisBasis();

  // resize the radius basis
  int rank = _tetMesh->rank();
  radiusBasis.resizeAndWipe(rank);
  radiusBasis3.resizeAndWipe(3, rank);

  axisBasis.resizeAndWipe(6, rank);

  // make the radius basis
  for (int x = 0; x < rank; x++)
  {
    // extract /bar/UU_j
    VEC3F centerColumn;
    centerColumn[0] = centerBasis(0, x);
    centerColumn[1] = centerBasis(1, x);
    centerColumn[2] = centerBasis(2, x);

    // for each contained vertex, subtract out the mean
    // and take the norm
    Real maxFound = 0.0;
    for (int y = 0; y < vertexBasis.rows() / 3; y++)
    {
      VEC3F basisColumn;
      basisColumn[0] = vertexBasis(3 * y, x);
      basisColumn[1] = vertexBasis(3 * y + 1, x);
      basisColumn[2] = vertexBasis(3 * y + 2, x);

      // component-wise search
      VEC3F diff = basisColumn - centerColumn;
      for (int z = 0; z < 3; z++)
        if (fabs(diff[z]) > maxFound)
          maxFound = fabs(diff[z]);

      if (diff[0] > radiusBasis3(0,x))
        radiusBasis3(0,x) = diff[0];
      if (diff[1] > radiusBasis3(1,x))
        radiusBasis3(1,x) = diff[1];
      if (diff[2] > radiusBasis3(2,x))
        radiusBasis3(2,x) = diff[2];

      if (diff[0] > 0.0 && diff[0] > axisBasis(0,x))
        axisBasis(0,x) = diff[0];

      if (diff[1] > 0.0 && diff[1] > axisBasis(1,x))
        axisBasis(1,x) = diff[1];

      if (diff[2] > 0.0 && diff[2] > axisBasis(2,x))
        axisBasis(2,x) = diff[2];

      if (diff[3] < 0.0 && diff[3] < axisBasis(3,x))
        axisBasis(3,x) = diff[3];

      if (diff[4] < 0.0 && diff[4] < axisBasis(4,x))
        axisBasis(4,x) = diff[4];

      if (diff[5] < 0.0 && diff[5] < axisBasis(5,x))
        axisBasis(5,x) = diff[5];
    }

    // store the largest found in the basis
    radiusBasis[x] = maxFound;
  }

  axisBasis.absoluteValue();
}

//////////////////////////////////////////////////////////////////////
// build the radius basis for a node
//////////////////////////////////////////////////////////////////////
void BD_TREE::buildRadiusMultiBasis(BD_TREE_NODE& node)
{
  MATRIX& vertexBasis = node.vertexBasis();
  MATRIX& centerBasis = node.centerBasis();
  MATRIX& multiBasis = node.radiusMultiBasis();

  // multibasis size for the time being
  int multiRows = vertexBasis.rows() / 3;

  // resize the center basis
  int rank = _tetMesh->rank();
  multiBasis.resizeAndWipe(multiRows, rank);

  // randomly distribute the vertices
  vector<vector<int> > vertexBins;
  vertexBins.resize(multiRows);
  for (int x = 0; x < vertexBasis.rows() / 3; x++)
    vertexBins[x].push_back(x);

  // make the radius basis for each row in the multibasis
  for (int i = 0; i < multiRows; i++)
    for (int x = 0; x < rank; x++)
    {
      // extract /bar/UU_j
      VEC3F centerColumn;
      centerColumn[0] = centerBasis(0, x);
      centerColumn[1] = centerBasis(1, x);
      centerColumn[2] = centerBasis(2, x);

      // for each contained vertex, subtract out the mean
      // and take the norm
      Real maxFound = 0.0;
      for (unsigned int y = 0; y < vertexBins[i].size(); y++)
      {
        int vertexIndex = vertexBins[i][y];

        VEC3F basisColumn;
        basisColumn[0] = vertexBasis(3 * vertexIndex, x);
        basisColumn[1] = vertexBasis(3 * vertexIndex + 1, x);
        basisColumn[2] = vertexBasis(3 * vertexIndex + 2, x);

        // vector-wise search
        VEC3F diff = basisColumn;
        Real diffNorm = norm(diff);
          if (diffNorm > maxFound)
            maxFound = diffNorm;
      }

      // store the largest found in the multi basis
      multiBasis(i,x) = maxFound;
    }
}

//////////////////////////////////////////////////////////////////////
// do a plane-sphere overlap test
//////////////////////////////////////////////////////////////////////
bool BD_TREE::overlap(SPHERE& sphere, PLANE* plane)
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
bool BD_TREE::overlap(TRIANGLE& triangle, PLANE* plane)
{
  VEC3F point = plane->point();
  Real dots[3];
  for (int x = 0; x < 3; x++)
  {
    VEC3F diff = (*triangle.vertex(x)) - point;
    dots[x] = plane->normal() * diff;
  }

  for (int x = 0; x < 3; x++)
    if (dots[x] < plane->stickiness())
      return true;

  return false;
}

//////////////////////////////////////////////////////////////////////
// do a plane-triangle overlap test
//////////////////////////////////////////////////////////////////////
bool BD_TREE::overlap(TRIANGLE& triangle, SPHERE* sphere)
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

//////////////////////////////////////////////////////////////////////
// print out stats on the tree, particularly the root node
//////////////////////////////////////////////////////////////////////
void BD_TREE::printStats()
{
  SPHERE sphere = _root->boundingSphere();
  VECTOR q = _tetMesh->q();
  VEC3F centerPrime = _root->centerBasis() * q;
  q.fabs();
  Real radiusPrime = _root->radiusBasis() * q;
  VECTOR multiRadius  = _root->radiusMultiBasis() * q;

  VECTOR radiusPrime3 = _root->radiusBasis3() * q;
  Real finalPrime3 = sqrt(radiusPrime3 * radiusPrime3);

  cout << "=====================================" << endl;
  cout << " BD-Tree stats" << endl;
  cout << "=====================================" << endl;
  cout << " Top radius: " << sphere.radius() << endl;
  cout << " Radius delta: " << radiusPrime << endl;
  cout << " Top center: " << sphere.center() << endl;
  cout << " Radius3: " << finalPrime3 << endl;
  BD_TREE_NODE& node = *_root;

  // find actual max displacement
  VECTOR displacement = node.vertexBasis() * q;
  vector<int>& vertexIDs = node.vertexIDs();
  vector<VEC3F>& restPose = _tetMesh->restPose();
  Real maxFoundTrue = 0.0;
  VEC3F maxDeltaTrue;
  Real maxFound = 0.0;
  VEC3F maxDelta;
  Real absFound = 0.0;
  VEC3F absDelta;

  MATRIX absU = node.vertexBasis();
  absU.absoluteValue();
  VECTOR qAbs = q;
  qAbs.fabs();
  VECTOR absDisplacement = absU * q;
  VEC3F newCenter = sphere.center() + centerPrime;

  for (unsigned int x = 0; x < vertexIDs.size(); x++)
  {
    VEC3F absDelta;
    absDelta[0] = absDisplacement[3 * x];
    absDelta[1] = absDisplacement[3 * x + 1];
    absDelta[2] = absDisplacement[3 * x + 2];
    VEC3F absDiff = absDelta - centerPrime;
    if (sqrt(absDiff * absDiff) > maxFound)
    {
      absFound = sqrt(absDiff * absDiff);
      absDelta = absDiff; 
    }

    VEC3F delta;
    delta[0] = displacement[3 * x];
    delta[1] = displacement[3 * x + 1];
    delta[2] = displacement[3 * x + 2];

    VEC3F diff = delta - centerPrime;
    if (sqrt(diff * diff) > maxFound)
    {
      maxFound = sqrt(diff * diff);
      maxDelta = diff; 
    }

    delta = (delta + restPose[vertexIDs[x]]) - newCenter;

    if (sqrt(delta * delta) > maxFoundTrue)
    {
      maxFoundTrue = sqrt(delta * delta);
      maxDeltaTrue = delta;
    }
  }

  // find max rest node
  VEC3F maxRest;
  for (unsigned int x = 0; x < vertexIDs.size(); x++)
  {
    VEC3F vertex = restPose[vertexIDs[x]];
    if (fabs(vertex[0]) > maxRest[0])
      maxRest[0] = fabs(vertex[0]);
    if (fabs(vertex[1]) > maxRest[1])
      maxRest[1] = fabs(vertex[1]);
    if (fabs(vertex[2]) > maxRest[2])
      maxRest[2] = fabs(vertex[2]);
  }

  VECTOR allMaxDelta = node.radiusBasis3() * qAbs;
  VECTOR inside2 = maxRest.toVector() + allMaxDelta - newCenter.toVector();

  VECTOR axialRadius = node.axisBasis() * qAbs;
  Real xComp = axialRadius[0] > axialRadius[3] ? axialRadius[0] : axialRadius[3];
  Real yComp = axialRadius[1] > axialRadius[4] ? axialRadius[1] : axialRadius[4];
  Real zComp = axialRadius[2] > axialRadius[5] ? axialRadius[2] : axialRadius[5];

  Real finalAxial =  sqrt(xComp * xComp + yComp * yComp + zComp * zComp);

  cout << " original radius: " << sphere.radius() << endl;
  cout << " max displacement: " << maxDeltaTrue << " " << maxFoundTrue << endl;
  cout << " true radius delta: " << maxFoundTrue - sphere.radius() << endl;

  // single triangle ineaulity
  cout << " displacement based delta: " << maxDelta << " " << maxFound << endl;

  // triangle inequality with absolute vals
  cout << " abs displacement based delta: " << absDelta << " " << absFound << endl;
  cout << " all max delta: " << allMaxDelta << " " << allMaxDelta.norm2() << endl;
  cout << " rank 1 delta: " << node.radiusBasis() * qAbs << endl;
  cout << " inside 2 delta: " << inside2.norm2()  << endl;
  cout << " axial delta: [" << xComp << ", " << yComp << ", " << zComp << "] " << finalAxial << endl;
}

//////////////////////////////////////////////////////////////////////
// do an AABB -sphere overlap test
//
// A Simple Method for Box-Sphere Intersection Testing
// by Jim Arvo
// from "Graphics Gems", Academic Press, 1990
//////////////////////////////////////////////////////////////////////
bool BD_TREE::overlap(SPHERE& sphere, BOX* box)
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

  return false;
}
