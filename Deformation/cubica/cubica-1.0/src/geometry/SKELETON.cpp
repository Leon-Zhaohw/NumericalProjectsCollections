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
// SKELETON.cpp: implementation of the SKELETON class.
//
//////////////////////////////////////////////////////////////////////

#include "SKELETON.h"
#include <MERSENNETWISTER.h>

//////////////////////////////////////////////////////////////////////
// Constructor for the skeleton
//////////////////////////////////////////////////////////////////////
SKELETON::SKELETON(TET_MESH* tetMesh, int bccRes) :
  _rootBone(NULL),
  _selected(0),
  _bccRes(bccRes),
  _tetMesh(tetMesh),
  _originalObj(NULL)
{
  // set the size of the skinning transforms vector
  vector<VEC3F>& vertices = _tetMesh->vertices();
  int size = vertices.size();
  _skinningTransforms.resize(vertices.size());
  _skinningTransformsOld.resize(vertices.size());
  VEC3F zero(0.0, 0.0, 0.0);
  MATRIX3 I = MATRIX3::I();
  for (int x = 0; x < size; x++)
  {
    _skinningTransforms[x] = MATRIX(I, zero);
    _skinningTransformsOld[x] = MATRIX(I, zero);
  }
 
  // tell the tets about the skinning transforms
  vector<TET>& tets = _tetMesh->tets();
  for (unsigned int x = 0; x < tets.size(); x++)
  {
    TET& tet = tets[x];
    for (int y = 0; y < 4; y++)
    {
      int vertexID = _tetMesh->vertexID(tet.vertices[y]);
      tet.skinningMatrix(y) = &(_skinningTransforms[vertexID]);
    }
  }

  // cache the location of all the skinning matrices
  for (int x = 0; x < size; x++)
  {
    _transformID[&(_skinningTransforms[x])] = x;
    _transformOldID[&(_skinningTransformsOld[x])] = x;
  }
  
  // cache the bone colors
  _colors[0][0] = 1.0f; _colors[0][1] = 0.0f; _colors[0][2] = 0.0f; _colors[0][3] = 1.0f;
  _colors[1][0] = 0.0f; _colors[1][1] = 1.0f; _colors[1][2] = 0.0f; _colors[1][3] = 1.0f;
  _colors[2][0] = 0.0f; _colors[2][1] = 0.0f; _colors[2][2] = 1.0f; _colors[2][3] = 1.0f;
  _colors[3][0] = 1.0f; _colors[3][1] = 0.0f; _colors[3][2] = 1.0f; _colors[3][3] = 1.0f;
  _colors[4][0] = 0.0f; _colors[4][1] = 1.0f; _colors[4][2] = 1.0f; _colors[4][3] = 1.0f;
  _colors[5][0] = 1.0f; _colors[5][1] = 1.0f; _colors[5][2] = 0.0f; _colors[5][3] = 1.0f;
  _colors[6][0] = 1.0f; _colors[6][1] = 1.0f; _colors[6][2] = 1.0f; _colors[6][3] = 1.0f;
  _colors[7][0] = 0.5f; _colors[7][1] = 0.5f; _colors[7][2] = 0.5f; _colors[7][3] = 1.0f;

  // hacks
  _first = true;
  _second = true;
  _third = true;
  _fourth = true;
}

SKELETON::~SKELETON()
{
  cleanup(_rootBone);
  if (_originalObj)
    delete _originalObj;
}

SKELETON::SKELETON(const string& filename, bool odeFile) :
  _rootBone(NULL),
  _selected(0),
  _bccRes(32),
  _tetMesh(NULL),
  _skinnedObj(NULL),
  _originalObj(NULL)
{
  // try to open the file
  FILE* file = fopen(filename.c_str(), "r");
  if (file == NULL)
  {
    cout << " No skeleton file " << filename.c_str() << " found!" << endl;
    return;
  }

  // read in each bone
  while (!feof(file))
  {
    int index;
    float trans[3];
    int previousIndex;
    fscanf(file,"%i %f %f %f %i\n", &index, &trans[0], &trans[1], &trans[2], &previousIndex);

    QUATERNION rotation;
    Real boneLength = 0.0;
    Real boneRadius = 0.0;
    if (odeFile)
    {
      float quat[4];
      fscanf(file,"%f %f %f %f\n", &quat[0], &quat[1], &quat[2], &quat[3]);
      rotation[0] = quat[0]; rotation[1] = quat[1]; 
      rotation[2] = quat[2]; rotation[3] = quat[3];

      float length, radius;
      fscanf(file,"%f %f\n", &length, &radius);
      boneLength = length;
      boneRadius = radius;
    }

    BONE* previousBone = NULL;
    if (previousIndex >= 0)
      previousBone = _bones[previousIndex];

    VEC3F translation(trans[0], trans[1], trans[2]);
    //MATRIX3 I = MATRIX3::I();
    MATRIX3 I = rotation.toExplicitMatrix3x3();
    BONE* newBone = new BONE(previousBone, translation, I);
    newBone->globalPosition() = translation;
    newBone->globalRestPosition() = translation;

    if (previousIndex < 0)
      _rootBone = newBone;

    _bones.push_back(newBone);
    _boneLengths.push_back(boneLength);
    _boneRadii.push_back(boneRadius);
  }

  fclose(file);

  cout << " Read in skeleton containing " << _bones.size() << " bones " << endl;

  // cache the bone colors
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
// recursively delete the skeleton
//////////////////////////////////////////////////////////////////////
void SKELETON::cleanup(BONE* bone)
{
  if (bone == NULL) return;

  vector<BONE*>& children = bone->children();
  for (unsigned int x = 0; x < children.size(); x++)
    cleanup(children[x]);

  delete bone;
}

//////////////////////////////////////////////////////////////////////
// build some skeleton
//////////////////////////////////////////////////////////////////////
void SKELETON::buildSkeleton()
{
  MATRIX3 I = MATRIX3::I();
  
  // initial vertex
  VEC3F first(0.0, 0.0, 0.0);
  _rootBone = new BONE(NULL, first, I);

  // middle joint
  VEC3F second(0.0, 0.333, 0.0);
  BONE* middle = new BONE(_rootBone, second, I);

  // top joint
  BONE* top = new BONE(middle, second, I);

  // top joint
  BONE* last = new BONE(top, second, I);

  // store all the bones for reference
  _bones.push_back(_rootBone);
  _bones.push_back(middle);
  _bones.push_back(top);
  _bones.push_back(last);
}

//////////////////////////////////////////////////////////////////////
// build skeleton for the dragon model
//////////////////////////////////////////////////////////////////////
void SKELETON::buildDragonSkeleton()
{
  MATRIX3 I = MATRIX3::I();
  
  // initial vertex
  VEC3F first(0.0, 0.0, 0.0);
  _rootBone = new BONE(NULL, first, I);

  // middle joint
  VEC3F second(-0.333, 0.0, 0.0);
  BONE* middle = new BONE(_rootBone, second, I);

  // top joint
  BONE* top = new BONE(middle, second, I);

  // top joint
  BONE* last = new BONE(top, second, I);

  // store all the bones for reference
  _bones.push_back(_rootBone);
  _bones.push_back(middle);
  _bones.push_back(top);
  _bones.push_back(last);
}

//////////////////////////////////////////////////////////////////////
// build skeleton for the head model
//////////////////////////////////////////////////////////////////////
void SKELETON::buildHeadSkeleton()
{
  MATRIX3 I = MATRIX3::I();
  
  // initial vertex
  VEC3F first(0.0, 0.0, 0.0);
  _rootBone = new BONE(NULL, first, I);

  // middle joint
  VEC3F second(-0.003, 0.504, 0.563);
  second = headNormalize(second, _bccRes) - _meshOrigin;
  BONE* middle = new BONE(_rootBone, second, I);

  // store all the bones for reference
  _bones.push_back(_rootBone);
  _bones.push_back(middle);
}

//////////////////////////////////////////////////////////////////////
// recursively update the skeleton
//////////////////////////////////////////////////////////////////////
void SKELETON::updateBones(BONE* bone, MATRIX transform)
{
  const MATRIX3 rotation = bone->rotation();
  VEC3F& translation = bone->translation();

  // load the position into a homogeneous vector
  VECTOR bonePosition(4);
  bonePosition[0] = translation[0];
  bonePosition[1] = translation[1];
  bonePosition[2] = translation[2];
  bonePosition[3] = 1.0;

  // transform the position
  bonePosition = transform * bonePosition;

  // repack the result into non-homogeneous coordinates
  VEC3F& globalPosition = bone->globalPosition();
  globalPosition[0] = bonePosition[0];
  globalPosition[1] = bonePosition[1];
  globalPosition[2] = bonePosition[2];

  // push new matrix on the stack
  MATRIX boneTransform(rotation, translation);
  transform = transform * boneTransform;

  // store the global transform
  bone->globalTransform() = transform;
  
  vector<BONE*>& children = bone->children();
  for (unsigned int x = 0; x < children.size(); x++)
    updateBones(children[x], transform);
}

//////////////////////////////////////////////////////////////////////
// Recreate the skinning matrix for a given time
//////////////////////////////////////////////////////////////////////
MATRIX SKELETON::recreateHeadSkinning(int vertexID, float time)
{
  ///////////////////////////////////////////////////////////////////
  // recreate "updateBones"
  ///////////////////////////////////////////////////////////////////
  bool dummy = false;
  MATRIX3 rotation = headShake(time, 1.0, dummy);
  VEC3F& translation = _bones[1]->translation();

  // load the position into a homogeneous vector
  VECTOR bonePosition(4);
  bonePosition[0] = translation[0];
  bonePosition[1] = translation[1];
  bonePosition[2] = translation[2];
  bonePosition[3] = 1.0;

  // get the base transform -- this should never have changed
  MATRIX transform = _bones[0]->globalTransform();

  // transform the position
  bonePosition = transform * bonePosition;

  // repack the result into non-homogeneous coordinates
  VEC3F& globalPosition = _bones[1]->globalPosition();
  globalPosition[0] = bonePosition[0];
  globalPosition[1] = bonePosition[1];
  globalPosition[2] = bonePosition[2];

  // push new matrix on the stack
  MATRIX boneTransform(rotation, translation);
  transform = transform * boneTransform;

  // store the global transform
  MATRIX globalTransform = transform;

  ///////////////////////////////////////////////////////////////////
  // recreate "updateSkinning"
  ///////////////////////////////////////////////////////////////////
  // get the skinning vector
  vector<pair<int, Real> > weights = _skinning[vertexID];

  // compute cumulative transform from each attached bone
  MATRIX skinningTransform(4,4);
  MATRIX3 I = MATRIX3::I();
  for (unsigned int y = 0; y < weights.size(); y++)
  {
    // get the bone
    BONE* bone = _bones[weights[y].first];
    
    // get the global rest position so we can compute
    // the transform to local coordinates
    VEC3F globalRestPosition = bone->globalRestPosition();
    globalRestPosition *= -1.0;
    MATRIX transformToLocal(I, globalRestPosition);

    // get the local bone transform
    MATRIX& boneTransform = bone->globalTransform();

    // if it's the neck bone, hijack it with the matrix we just computed
    if (weights[y].first == 1)
      boneTransform = globalTransform;

    // compose the final transform
    Real weight = weights[y].second;
    MATRIX finalMatrix = boneTransform * transformToLocal;
    finalMatrix = weight * finalMatrix;
    skinningTransform += finalMatrix;
  }

  // add the to and from mesh origin translation
  VEC3F negated = -1.0 * _meshOrigin;
  MATRIX toOrigin(I, negated);
  MATRIX fromOrigin(I, _meshOrigin);
  skinningTransform = skinningTransform * toOrigin;
  skinningTransform = fromOrigin * skinningTransform;

  return skinningTransform; 
}

//////////////////////////////////////////////////////////////////////
// draw the skeleton
//////////////////////////////////////////////////////////////////////
void SKELETON::drawBones(BONE* bone)
{
  if (bone == NULL) return;

  glPushMatrix();
    glTranslatef(_meshOrigin[0], _meshOrigin[1], _meshOrigin[2]);

    // draw the endpoint
    VEC3F& current = bone->globalPosition();
    glPointSize(10.0f);
    if (bone == _bones[_selected])
      glColor4f(10.0, 10.0, 10.0, 1.0);
    else
      glColor4f(10.0, 0.0, 0.0, 1.0);
    glBegin(GL_POINTS);
      glVertex3f(current[0], current[1], current[2]);
    glEnd();

    // if it has a parent, draw the bone
    if (bone->parent() != NULL)
    {
      VEC3F& parent = bone->parent()->globalPosition();
      /*
      glPointSize(10.0f);
        glColor4f(10.0, 0.0, 0.0, 1.0);
      glBegin(GL_POINTS);
        glVertex3f(parent[0], parent[1], parent[2]); 
      glEnd();
      */
      
      glLineWidth(3.0f);
      glColor4f(0.0, 10.0, 0.0, 1.0);
      glBegin(GL_LINES);
        glVertex3f(current[0], current[1], current[2]); 
        glVertex3f(parent[0], parent[1], parent[2]); 
      glEnd();
    }
  glPopMatrix();

  // draw the children
  vector<BONE*>& children = bone->children();
  for (unsigned int x = 0; x < children.size(); x++)
    drawBones(children[x]);
}

//////////////////////////////////////////////////////////////////////
// build some skinning
//////////////////////////////////////////////////////////////////////
void SKELETON::buildSkinning()
{
  // blow away any previous skinning
  _skinning.clear();
  
  // get rest vertices
  vector<VEC3F>& restPose = _tetMesh->restPose();

  // for each vertex
  for (unsigned int x = 0; x < restPose.size(); x++)
  {
    vector<pair<int, Real> > weights;

    if (restPose[x][1] < 1.0 / 6.0)
    {
      pair<int, Real> weight;
      weight.second = 1.0;
      weight.first = 0;
      weights.push_back(weight);
    }
    else if (restPose[x][1] < 3.0 / 6.0)
    {
      Real ramp = (restPose[x][1] - 1.0 / 6.0) / (1.0 / 3.0);

      pair<int, Real> weight0;
      weight0.first = 0;
      weight0.second = 1.0 - ramp;
      weights.push_back(weight0);

      pair<int, Real> weight1;
      weight1.first = 1;
      weight1.second = ramp;
      weights.push_back(weight1);
    }
    else if (restPose[x][1] < 5.0 / 6.0)
    {
      Real ramp = (restPose[x][1] - 3.0 / 6.0) / (1.0 / 3.0);

      pair<int, Real> weight0;
      weight0.first = 1;
      weight0.second = 1.0 - ramp;
      weights.push_back(weight0);

      pair<int, Real> weight1;
      weight1.first = 2;
      weight1.second = ramp;
      weights.push_back(weight1);
    }
    else 
    {
      pair<int, Real> weight;
      weight.second = 1.0;
      weight.first = 2;
      weights.push_back(weight);
    }

    _skinning.push_back(weights);
  }
}

//////////////////////////////////////////////////////////////////////
// build some skinning
//////////////////////////////////////////////////////////////////////
void SKELETON::buildDragonSkinning()
{
  // blow away any previous skinning
  _skinning.clear();
  
  // get rest vertices
  vector<VEC3F>& restPose = _tetMesh->restPose();

  // for each vertex
  for (unsigned int x = 0; x < restPose.size(); x++)
  {
    vector<pair<int, Real> > weights;

    /*
    if (restPose[x][0] < 1.0 / 6.0)
    {
      pair<int, Real> weight;
      weight.second = 1.0;
      weight.first = 0;
      weights.push_back(weight);
    }
    else if (restPose[x][0] < 3.0 / 6.0)
    {
      Real ramp = (restPose[x][0] - 1.0 / 6.0) / (1.0 / 3.0);

      pair<int, Real> weight0;
      weight0.first = 0;
      weight0.second = 1.0 - ramp;
      weights.push_back(weight0);

      pair<int, Real> weight1;
      weight1.first = 1;
      weight1.second = ramp;
      weights.push_back(weight1);
    }
    else if (restPose[x][0] < 5.0 / 6.0)
    {
      Real ramp = (restPose[x][1] - 3.0 / 6.0) / (1.0 / 3.0);

      pair<int, Real> weight0;
      weight0.first = 1;
      weight0.second = 1.0 - ramp;
      weights.push_back(weight0);

      pair<int, Real> weight1;
      weight1.first = 2;
      weight1.second = ramp;
      weights.push_back(weight1);
    }
    else 
    {
      pair<int, Real> weight;
      weight.second = 1.0;
      weight.first = 2;
      weights.push_back(weight);
    }
    */
    if (restPose[x][0] > 5.0 / 6.0)
    {
      pair<int, Real> weight;
      weight.second = 1.0;
      weight.first = 0;
      weights.push_back(weight);
    }
    else if (restPose[x][0] > 3.0 / 6.0)
    {
      Real ramp = (restPose[x][0] - 3.0 / 6.0) / (1.0 / 3.0);

      pair<int, Real> weight0;
      weight0.first = 0;
      //weight0.second = 1.0 - ramp;
      weight0.second = ramp;
      weights.push_back(weight0);

      pair<int, Real> weight1;
      weight1.first = 1;
      //weight1.second = ramp;
      weight1.second = 1.0 - ramp;
      weights.push_back(weight1);
    }
    else if (restPose[x][0] > 1.0 / 6.0)
    {
      Real ramp = (restPose[x][0] - 1.0 / 6.0) / (1.0 / 3.0);

      pair<int, Real> weight0;
      weight0.first = 1;
      //weight0.second = 1.0 - ramp;
      weight0.second = ramp;
      weights.push_back(weight0);

      pair<int, Real> weight1;
      weight1.first = 2;
      //weight1.second = ramp;
      weight1.second = 1.0 - ramp;
      weights.push_back(weight1);
    }
    else 
    {
      pair<int, Real> weight;
      weight.second = 1.0;
      weight.first = 2;
      weights.push_back(weight);
    }

    _skinning.push_back(weights);
  }
}

//////////////////////////////////////////////////////////////////////
// build some skinning
//////////////////////////////////////////////////////////////////////
void SKELETON::buildHeadSkinning()
{
  if (readSkinning()) return;
  
  // load the bone boxes
  OBJ neckBone;
  OBJ headBone;
  //neckBone.Load("./Skinner/neck.bone.cylinder.obj");
  neckBone.Load("./Skinner/neck.bone.cylinder.lower.obj");
  headBone.Load("./Skinner/head.bone.obj");

  // normalize the bone boxes
  vector<VEC3>& neckVertices = neckBone.vertices;
  vector<VEC3>& headVertices = headBone.vertices;
  for (unsigned int x = 0; x < neckVertices.size(); x++)
    neckVertices[x] = headNormalize(neckVertices[x], _bccRes);
  for (unsigned int x = 0; x < headVertices.size(); x++)
    headVertices[x] = headNormalize(headVertices[x], _bccRes);

  neckBone.setBCCRes(_bccRes);
  neckBone.createAccelGrid();
  //neckBone.createDistanceGrid(_bccRes);
  headBone.setBCCRes(_bccRes);
  headBone.createAccelGrid();
  //headBone.createDistanceGrid(_bccRes);

  // blow away any previous skinning
  _skinning.clear();
  
  // get rest vertices
  vector<VEC3F>& restPose = _tetMesh->restPose();

  // for each vertex
  _skinning.resize(restPose.size());
  for (unsigned int x = 0; x < restPose.size(); x++)
  {
    vector<pair<int, Real> > weights;

    float vert[3];
    vert[0] = restPose[x][0];
    vert[1] = restPose[x][1];
    vert[2] = restPose[x][2];
    if (headBone.inside(vert))
    {
      pair<int, Real> weight;
      weight.second = 1.0;
      weight.first = 1;
      weights.push_back(weight);
    }
    else if (neckBone.inside(vert))
    {
      pair<int, Real> weight;
      weight.second = 1.0;
      weight.first = 0;
      weights.push_back(weight);
    }
    else
    {
      Real w0 = neckBone.bruteForceDistance(vert);
      Real w1 = headBone.bruteForceDistance(vert);
      Real sum = w0 + w1;
      pair<int, Real> weight0;
      weight0.second = 1.0 - w0 / sum;
      weight0.first = 0;
      weights.push_back(weight0);

      pair<int, Real> weight1;
      weight1.second = 1.0 - w1 / sum;
      weight1.first = 1;
      weights.push_back(weight1);
    }
    _skinning[x] = weights;
  }
  writeSkinning();
}

//////////////////////////////////////////////////////////////////////
// draw the skinning weights
//////////////////////////////////////////////////////////////////////
void SKELETON::drawSkinning()
{
  // get vertices
  vector<VEC3F>& vertices = _tetMesh->vertices();

  // for each vertex
  glPointSize(1.0f);
  glBegin(GL_POINTS);
  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    // get the skinning vector
    vector<pair<int, Real> > weights = _skinning[x];

    VEC3F color;
    for (unsigned int y = 0; y < weights.size(); y++)
    {
      int bone = weights[y].first;
      Real weight = weights[y].second;

      color[0] += weight * _colors[bone % 8][0] * 10;
      color[1] += weight * _colors[bone % 8][1] * 10;
      color[2] += weight * _colors[bone % 8][2] * 10;
    }

    // draw the endpoint
    glColor4f(color[0], color[1], color[2], 1.0);
    glVertex3f(vertices[x][0], vertices[x][1], vertices[x][2]); 
  }
  glEnd();
}

//////////////////////////////////////////////////////////////////////
// draw partitioning weights
//////////////////////////////////////////////////////////////////////
void SKELETON::drawPartitioningWeights()
{
  // compute the one rings if needed
  vector<vector<int> > oneRings;
  vector<TET>& tets = _tetMesh->tets();
  if (!readOneRings(oneRings))
  {
    for (unsigned int x = 0; x < tets.size(); x++)
    {
      vector<int> oneRing;
      _tetMesh->tetOneRing(x, oneRing);
      oneRings.push_back(oneRing);
    }
    writeOneRings(oneRings);
  }

  MERSENNETWISTER twister(9876);
  vector<VEC3F> colors;
  for (unsigned int x = 0; x < _bones.size(); x++)
  {
    VEC3F color(twister.rand(), twister.rand(), twister.rand());
    colors.push_back(color);
  }

  // find the max weights for each bone
  /*
  vector<float> maxWeights(_objSkinningWeights[0].size());
  for (unsigned int x = 0; x < _objSkinningWeights.size(); x++)
    for (unsigned int y = 0; y < _objSkinningWeights[x].size(); y++)
    {
      if (_objSkinningWeights[x][y] > maxWeights[y])
        maxWeights[y] = _objSkinningWeights[x][y];
    }
    */
  vector<float> maxWeights(_bones.size());
  for (unsigned int x = 0; x < maxWeights.size(); x++)
    maxWeights[x] = 0;
  for (unsigned int x = 0; x < _skinning.size(); x++)
    for (unsigned int y = 0; y < _skinning[x].size(); y++)
    {
      int boneID = _skinning[x][y].first;
      float weight = _skinning[x][y].second;

      if (weight > maxWeights[boneID])
        maxWeights[boneID] = weight;
    }

  // draw an edge between each tet center and its one ring
  glLineWidth(0.1);
  for (unsigned int x = 0; x < tets.size(); x++)
  {
    vector<int>& ring = oneRings[x];
    VEC3F center = tets[x].center();

    VEC3F finalColor;

    TET& tet = tets[x];
    float maxFound = 0;
    VEC3F maxColor;
    for (int y = 0; y < 4; y++)
    {
      // get the weight list of the current vertex
      int vertexID = _tetMesh->vertexID(tet.vertices[y]);
      vector<pair<int, Real> >& weights = _skinning[vertexID];

      // see if there's a big one
      for (unsigned int z = 0; z < weights.size(); z++)
        if (weights[z].second / maxWeights[weights[z].first] > maxFound)
        {
          maxFound = weights[z].second / maxWeights[weights[z].first];
          maxColor = colors[weights[z].first];
        }
    }
    finalColor = maxColor;
    finalColor *= 3;
    //finalColor *= 0.25;

    glBegin(GL_LINES);
      glColor4f(finalColor[0], finalColor[1], finalColor[2], 1.0);
      for (unsigned int y = 0; y < ring.size(); y++)
      {
        VEC3F neighborCenter = tets[ring[y]].center();
        glVertex3dv(center);
        glVertex3dv(neighborCenter);
      }
    glEnd();
  }

  /*
  // draw an edge between each tet center and its one ring
  glLineWidth(0.1);
  for (unsigned int x = 0; x < tets.size(); x++)
  {
    vector<int>& ring = oneRings[x];
    VEC3F center = tets[x].center();

    bool bigFound = false;
    TET& tet = tets[x];
    for (int y = 0; y < 4; y++)
    {
      // get the weight list of the current vertex
      int vertexID = _tetMesh->vertexID(tet.vertices[y]);
      vector<pair<int, Real> >& weights = _skinning[vertexID];

      // see if there's a big one
      for (unsigned int z = 0; z < weights.size(); z++)
        if (weights[z].second > 0.95 * maxWeights[weights[z].first])
        //if (weights[z].second > 0.9 * maxWeights[weights[z].first])
        //if (weights[z].second > 0.9 * maxWeights[weights[z].first] || weights[z].second > 0.9)
          bigFound = true;
    }

    if (!bigFound) continue;

    glBegin(GL_LINES);
      glColor4f(10,10,10,1);
      for (unsigned int y = 0; y < ring.size(); y++)
      {
        VEC3F neighborCenter = tets[ring[y]].center();
        glVertex3dv(center);
        glVertex3dv(neighborCenter);
      }
    glEnd();
  }
  */
  /*
  // compute the largest weights in each tet
  vector<float> largestWeights;
  float largestOverall = 0;
  for (unsigned int x = 0; x < tets.size(); x++)
  {
    // get the largest weight of the vertices in this tet
    TET& tet = tets[x];
    float largestFound = 0.0;
    for (int y = 0; y < 4; y++)
    {
      // get the weight list of the current vertex
      int vertexID = _tetMesh->vertexID(tet.vertices[y]);
      vector<pair<int, Real> >& weights = _skinning[vertexID];

      // see if there's a big one
      for (unsigned int z = 0; z < weights.size(); z++)
        if (weights[z].second > largestFound)
          largestFound = weights[z].second;
    }
    largestWeights.push_back(largestFound);

    if (largestFound > largestOverall)
      largestOverall = largestFound;
  }

  // draw an edge between each tet center and its one ring
  glLineWidth(0.1);
  glBegin(GL_LINES);
  for (unsigned int x = 0; x < tets.size(); x++)
  {
    vector<int>& ring = oneRings[x];
    VEC3F center = tets[x].center();

    for (unsigned int y = 0; y < ring.size(); y++)
    {
      float larger = largestWeights[x];
      if (largestWeights[ring[y]] > larger)
        larger = largestWeights[ring[y]];
      Real edgeColor = larger / largestOverall;
      edgeColor *= 10;

      //glColor4f(1,1,1,1);
      glColor4f(edgeColor,1 - edgeColor,1 - edgeColor, 10);
      VEC3F neighborCenter = tets[ring[y]].center();
      glVertex3dv(center);
      glVertex3dv(neighborCenter);
    }
  }
  glEnd();
  */
}

//////////////////////////////////////////////////////////////////////
// Update the tet mesh based on the vertex skinnings
//////////////////////////////////////////////////////////////////////
void SKELETON::updateSkinning()
{
  if (_tetMesh == NULL) return;

  // get vertices
  //vector<VEC3F>& restVertices = _tetMesh->restPose();
  vector<VEC3F>& vertices = _tetMesh->vertices();

  // for each vertex
  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    // get the skinning vector
    vector<pair<int, Real> > weights = _skinning[x];

    // compute cumulative transform from each attached bone
    MATRIX skinningTransform(4,4);
    MATRIX3 I = MATRIX3::I();
    for (unsigned int y = 0; y < weights.size(); y++)
    {
      // get the bone
      BONE* bone = _bones[weights[y].first];
      
      // get the global rest position so we can compute
      // the transform to local coordinates
      VEC3F globalRestPosition = bone->globalRestPosition();
      globalRestPosition *= -1.0;
      MATRIX transformToLocal(I, globalRestPosition);

      // get the local bone transform
      MATRIX& boneTransform = bone->globalTransform();

      // compose the final transform
      Real weight = weights[y].second;
      MATRIX finalMatrix = boneTransform * transformToLocal;
      finalMatrix = weight * finalMatrix;
      skinningTransform += finalMatrix;
    }

    // add the to and from mesh origin translation
    VEC3F negated = -1.0 * _meshOrigin;
    MATRIX toOrigin(I, negated);
    MATRIX fromOrigin(I, _meshOrigin);
    skinningTransform = skinningTransform * toOrigin;
    skinningTransform = fromOrigin * skinningTransform;
   
    // store the old skinning matrix
    _skinningTransformsOld[x] = _skinningTransforms[x];

    // store the final transform matrix
    _skinningTransforms[x] = skinningTransform;
  }
}

//////////////////////////////////////////////////////////////////////
// add the skinning displacements to the tet mesh for display
//////////////////////////////////////////////////////////////////////
void SKELETON::addSkinningDisplacements()
{
  // get vertices
  vector<VEC3F>& vertices = _tetMesh->vertices();
  vector<VEC3F>& restPose = _tetMesh->restPose();
  VECTOR& displacement = _tetMesh->x();
  int unconstrainedNodes = _tetMesh->unconstrainedNodes();

  // for each vertex
  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    // pack the position into a homogeneous vector
    VECTOR homogeneous(4);
    homogeneous[0] = restPose[x][0];
    homogeneous[1] = restPose[x][1];
    homogeneous[2] = restPose[x][2];
    homogeneous[3] = 1.0;

    // add the displacement if one exists
    if (x < (unsigned int)unconstrainedNodes)
    {
      int index = 3 * x;
      homogeneous[0] += displacement[index];
      homogeneous[1] += displacement[index + 1];
      homogeneous[2] += displacement[index + 2];
    }
    homogeneous = _skinningTransforms[x] * homogeneous;

    // compute final skinned position
    VEC3F final;
    final[0] = homogeneous[0];
    final[1] = homogeneous[1];
    final[2] = homogeneous[2];

    // store the final position
    vertices[x] = final;
  }
}

//////////////////////////////////////////////////////////////////////
// undo the skinning displacements to the tet mesh for simulation
//////////////////////////////////////////////////////////////////////
void SKELETON::undoSkinningDisplacements()
{
  // get vertices
  vector<VEC3F>& vertices = _tetMesh->vertices();
  vector<VEC3F>& restPose = _tetMesh->restPose();
  VECTOR& displacement = _tetMesh->x();
  int unconstrainedNodes = _tetMesh->unconstrainedNodes();

  // for each vertex
  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    // add the displacement if one exists
    if (x < (unsigned int)unconstrainedNodes)
    {
      int index = 3 * x;
      vertices[x][0] = restPose[x][0] + displacement[index];
      vertices[x][1] = restPose[x][1] + displacement[index + 1];
      vertices[x][2] = restPose[x][2] + displacement[index + 2];
    }
    else
      vertices[x] = restPose[x];
  }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
int SKELETON::headShakeKeyframed(float time)
{
  float fracDuration = 1;
  //float fracDuration = 0.5;
  //float fracPose = 0.75;
  //float fracPose = 0.5;
  float fracPose = 0.6;
  
  // just so we're dealing with concrete frame times
  time *= 60.0;
  
  vector<MATRIX3> rotations;
  vector<int> durations;

  // push back all the rotations desired
  rotations.push_back(MATRIX3::I());

  // transition lean left
  rotations.push_back(multiRotation(0,-17,0, fracPose));
  durations.push_back(20 * fracDuration);

  // hold lean left
  rotations.push_back(multiRotation(0,-17,0, fracPose));
  durations.push_back(20 * fracDuration);
  
  // transition stretch
  rotations.push_back(multiRotation(-4,-21,-32, fracPose));
  durations.push_back(20 * fracDuration);

  // hold stretch
  rotations.push_back(multiRotation(-4,-21,-32, fracPose));
  durations.push_back(20 * fracDuration);
 
  for (int x = 0; x < 10; x++)
  {
    // Shake right:
    rotations.push_back(multiRotation(-12 + 4 * x, 21,-22, fracPose));
    durations.push_back(4 * fracDuration);

    // Shake left:
    rotations.push_back(multiRotation(-12 + 4 * x + 1,-21,22, fracPose));
    durations.push_back(4 * fracDuration);
  }

  int totalFrames = 90 * fracDuration;
  for (int x = 0; x < totalFrames; x++)
  {
    // y is head twist
    // z is side to side
    // x is up and down
    rotations.push_back(multiRotation(0, 
          20 * cos(2.0 * M_PI * 2.0 * x / (totalFrames / 3)),
          20 * cos(2.0 * M_PI * 3.0 * x / (totalFrames / 3)), fracPose));
    if (x == 0)
      durations.push_back(4 * fracDuration);
    else if (x == totalFrames - 1)
      durations.push_back(8 * fracDuration);
    else  
      durations.push_back(1 * fracDuration);
  }
  
  // back to rest
  rotations.push_back(MATRIX3::I());
  durations.push_back(20 * fracDuration);
  
  // find which slice we are in
  int durationsSeen = 0;
  int beginIndex = -1;
  for (unsigned int x = 0; x < durations.size(); x++)
  {
    durationsSeen += durations[x];
    if (durationsSeen >= time)
    {
      beginIndex = x;
      durationsSeen -= durations[x];
      x = durations.size();
    }
  }

  // get the matrix we want
  BONE* selected = _bones[1];
  const MATRIX3 rotation = selected->rotation();
  QUATERNION& quaternion = selected->quaternion();
  
  // if no frames are left, just return identity
  if (beginIndex == -1)
  {
    quaternion = QUATERNION(MATRIX3::I());
    //rotation = MATRIX3::I();
    //return MATRIX3::I();
  }

  // calc the lerp
  MATRIX3 oldMatrix;
  MATRIX3 newMatrix;
  if (beginIndex == -1)
  {
    oldMatrix = MATRIX3::I();
    newMatrix = MATRIX3::I();
  }
  else
  {
    oldMatrix = rotations[beginIndex];
    newMatrix = rotations[beginIndex + 1];
  }

  float lerp = (time - durationsSeen) / durations[beginIndex];
  if (beginIndex == -1)
    lerp = 1.0;
  float cubic = -2.0f * pow(lerp, 3.0f) + 3.0f * pow(lerp, 2.0f);
  //rotation = (1 - cubic) * oldMatrix + cubic * newMatrix;
  quaternion = QUATERNION((1 - cubic) * oldMatrix + cubic * newMatrix);

  int totalDurations = 0;
  for (unsigned int x = 0; x < durations.size(); x++)
    totalDurations += durations[x];
  
  return totalDurations;
}



//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//void SKELETON::headShake(float time, float speed)
MATRIX3 SKELETON::headShake(float time, float speed, bool& forceFull)
{
  time *=  2.0 * M_PI * speed;

  Real period = 2.0 * M_PI;

  Real twistAmp = -0.125;
  Real leftRightAmp = 0.125;
  Real upDownAmp = 0.125;

  forceFull = false;
  if (time > 4 * period && time < 6 * period)
  {
    upDownAmp *= 0.25;
    time *= 4.0;
    if (_first)
    {
      _first = false;
      forceFull = true;
    }
  }
  else if (time > 6 * period && time < 8 * period)
  {
    leftRightAmp *= 0.25;
    time *= 4.0;
    if (_second)
    {
      _second = false;
      forceFull = true;
    }
  }
  else if (time > 8 * period)
  {
    upDownAmp = 0.0;
    leftRightAmp = 0.0;
    twistAmp = 0.0;
    if (_third)
    {
      _third = false;
      forceFull = true;
    }
  }
  else
  {
    if (_fourth)
    {
      _fourth = false;
      forceFull = true;
    }
  }

  BONE* selected = _bones[1];
  const MATRIX3 rotation = selected->rotation();
  QUATERNION& quaternion = selected->quaternion();
  MATRIX3 zRotate = MATRIX3::rotation(VEC3F(0.0, 0.0, 1.0), twistAmp * sin(time));    // twist
  MATRIX3 yRotate = MATRIX3::rotation(VEC3F(0.0, 1.0, 0.0), leftRightAmp * sin(time));      // left right
  MATRIX3 xRotate = MATRIX3::rotation(VEC3F(1.0, 0.0, 0.0), upDownAmp * sin(2 * time));  // up down
  //rotation = zRotate * yRotate * xRotate;
  quaternion = QUATERNION(zRotate * yRotate * xRotate);

  return rotation;
}

//////////////////////////////////////////////////////////////////////
// cycle the selected bone on positive z axis
//////////////////////////////////////////////////////////////////////
void SKELETON::rotateSelectedZ(Real sign)
{
  BONE* selected = _bones[_selected];
  const MATRIX3 rotation = selected->rotation();
  MATRIX3 zRotate = MATRIX3::rotation(VEC3F(0.0, 0.0, 1.0), sign * M_PI / 300);
  selected->quaternion() = QUATERNION(rotation * zRotate);
  //rotation = rotation * zRotate;
  //updateSkeleton();
  //addSkinningDisplacements();
}

//////////////////////////////////////////////////////////////////////
// cycle the selected bone on positive y axis
//////////////////////////////////////////////////////////////////////
void SKELETON::rotateSelectedY(Real sign)
{
  BONE* selected = _bones[_selected];
  const MATRIX3 rotation = selected->rotation();
  MATRIX3 yRotate = MATRIX3::rotation(VEC3F(0.0, 1.0, 0.0), sign * M_PI / 300);
  //rotation = rotation * yRotate;
  selected->quaternion() = QUATERNION(rotation * yRotate);
  //updateSkeleton();
  //addSkinningDisplacements();
}

//////////////////////////////////////////////////////////////////////
// cycle the selected bone on positive x axis
//////////////////////////////////////////////////////////////////////
void SKELETON::rotateSelectedX(Real sign)
{
  BONE* selected = _bones[_selected];
  const MATRIX3 rotation = selected->rotation();
  MATRIX3 xRotate = MATRIX3::rotation(VEC3F(1.0, 0.0, 0.0), sign * M_PI / 300);
  //rotation = rotation * xRotate;
  selected->quaternion() = QUATERNION(rotation * xRotate);
  //updateSkeleton();
  //addSkinningDisplacements();
}

//////////////////////////////////////////////////////////////////////
// normalize according to the head mesh
//////////////////////////////////////////////////////////////////////
VEC3F SKELETON::headNormalize(VEC3F vertex, int res)
{
  VEC3F centerOfMass(0.0197512, 0.427389, -0.0937424);
  Real maxVal = 1.61086;

  vertex = vertex - centerOfMass;
  double scale = 0.5 - 4.0 / res;
  vertex *= scale / maxVal;
  VEC3F half(0.5, 0.5, 0.5);
  vertex += half;
  return vertex;
}

//////////////////////////////////////////////////////////////////////
// write skinning
//////////////////////////////////////////////////////////////////////
bool SKELETON::writeSkinning()
{
  string tetMeshFilename = _tetMesh->filename();
  string filename = tetMeshFilename;
  filename += string(".skinning");
  FILE* file = fopen(filename.c_str(), "wb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Skinning file " << filename.c_str() << " could not be opened!" << endl;
    return false;
  }

  int skinningSize = _skinning.size();
  fwrite((void*)&skinningSize, sizeof(int), 1, file);

  for (int x = 0; x < skinningSize; x++)
  {
    int totalPairs = _skinning[x].size();
    fwrite((void*)&totalPairs, sizeof(int), 1, file);
    for (int y = 0; y < totalPairs; y++)
    {
      pair<int, Real> currentPair = _skinning[x][y];
      int bone = currentPair.first;
      double weight = currentPair.second;
      fwrite((void*)&bone, sizeof(int), 1, file);
      fwrite((void*)&weight, sizeof(double), 1, file);
    }
  }

  fclose(file);
  return true;
}

//////////////////////////////////////////////////////////////////////
// read skinning
//////////////////////////////////////////////////////////////////////
bool SKELETON::readSkinning()
{
  string tetMeshFilename = _tetMesh->filename();
  string filename = tetMeshFilename;
  filename += string(".skinning");
  FILE* file = fopen(filename.c_str(), "rb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Skinning file " << filename.c_str() << " could not be opened!" << endl;
    return false;
  }
  cout << " Found skinning file " << filename.c_str() << endl;

  int skinningSize;
  fread((void*)&skinningSize, sizeof(int), 1, file);
  _skinning.resize(skinningSize);

  for (int x = 0; x < skinningSize; x++)
  {
    int totalPairs;
    fread((void*)&totalPairs, sizeof(int), 1, file);
    for (int y = 0; y < totalPairs; y++)
    {
      int bone;
      double weight;
      fread((void*)&bone, sizeof(int), 1, file);
      fread((void*)&weight, sizeof(double), 1, file);
      pair<int, Real> currentPair;
      currentPair.first = bone;
      currentPair.second = weight;
      _skinning[x].push_back(currentPair);
    }
  }

  fclose(file);
  return true;
}

//////////////////////////////////////////////////////////////////////
// write out state
//////////////////////////////////////////////////////////////////////
bool SKELETON::writeState(string filename)
{
  FILE* file = fopen(filename.c_str(), "wb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __LINE__ << " : " << endl;
    cout << "Failed to open file: " << filename.c_str() << endl;
    return false;
  }
  cout << " Dumping state to file " << filename.c_str() << endl;

  // write how many skinning matrices there are
  int matricesSize = _skinningTransforms.size();
  fwrite((void*)&matricesSize, sizeof(int), 1, file);

  for (unsigned int x = 0; x < _skinningTransforms.size(); x++)
    _skinningTransforms[x].write(file);
  for (unsigned int x = 0; x < _skinningTransforms.size(); x++)
    _skinningTransformsOld[x].write(file);

  fclose(file);
  return true;
}

//////////////////////////////////////////////////////////////////////
// read in state
//////////////////////////////////////////////////////////////////////
bool SKELETON::readState(string filename)
{
  FILE* file = fopen(filename.c_str(), "rb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __LINE__ << " : " << endl;
    cout << "Failed to open file: " << filename.c_str() << endl;
    return false;
  }
  cout << " Reading state from file " << filename.c_str() << endl;

  // write how many skinning matrices there are
  int matricesSize;
  fread((void*)&matricesSize, sizeof(int), 1, file);

  for (int x = 0; x < matricesSize; x++)
    _skinningTransforms[x].read(file);
  for (unsigned int x = 0; x < _skinningTransforms.size(); x++)
    _skinningTransformsOld[x].read(file);

  fclose(file);
  return true;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void SKELETON::reset()
{
  VEC3F zero(0.0, 0.0, 0.0);
  MATRIX3 I = MATRIX3::I();
  for (unsigned int x = 0; x < _skinningTransforms.size(); x++)
  {
    _skinningTransforms[x] = MATRIX(I, zero);
    _skinningTransformsOld[x] = MATRIX(I, zero);
  }

  for (unsigned int x = 0; x < _bones.size(); x++)
    _bones[x]->reset();

  updateSkeleton();

  _first = true;
  _second = true;
  _third = true;
  _fourth = true;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
MATRIX3 SKELETON::multiRotation(int xTics, int yTics, int zTics, float fracPose)
{
  MATRIX3 final = MATRIX3::I();
  float sign = 1.0;
  
  sign = xTics < 0 ? -1.0 : 1.0;
  for (int x = 0; x < abs(xTics) * fracPose; x++)
  {
    MATRIX3 xRotate = MATRIX3::rotation(VEC3F(1.0, 0.0, 0.0), sign * M_PI / 300);
    final = final * xRotate;
  }

  sign = yTics < 0 ? -1.0 : 1.0;
  for (int x = 0; x < abs(yTics) * fracPose; x++)
  {
    MATRIX3 yRotate = MATRIX3::rotation(VEC3F(0.0, 1.0, 0.0), sign * M_PI / 300);
    final = final * yRotate;
  }

  sign = zTics < 0 ? -1.0 : 1.0;
  for (int x = 0; x < abs(zTics) * fracPose; x++)
  {
    MATRIX3 zRotate = MATRIX3::rotation(VEC3F(0.0, 0.0, 1.0), sign * M_PI / 300);
    final = final * zRotate;
  }

  return final;
}

//////////////////////////////////////////////////////////////////////
// Read in motion capture data from Pinocchio
//////////////////////////////////////////////////////////////////////
bool SKELETON::readMocap(string filename)
{
  FILE* file = fopen(filename.c_str(), "rb");

  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Could not open motion capture file! " << endl;
    return false;
  }

  // get high level data
  _totalFrames = 0;
  fread((void*)&_totalFrames, sizeof(int), 1, file);
  int totalBones = -1;
  fread((void*)&totalBones, sizeof(int), 1, file);

  if ((unsigned int)totalBones != _bones.size() - 1)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Number of bones in motion capture file do not match the skeleton! " << endl;
    cout << " In mocap: " << totalBones << endl;
    cout << " In skeleton: " << _bones.size() <<  endl;
    fclose(file);
    return false;
  }

  cout << " Reading in mocap data ...";flush(cout);

  // read in the motion capture
  _mocapRotations.clear();
  _mocapTranslations.clear();
  _mocapScalings.clear();

  _mocapRotations.resize(totalBones);
  _mocapTranslations.resize(totalBones);
  _mocapScalings.resize(totalBones);
  for (int x = 0; x < _totalFrames; x++)
  {
    if (feof(file))
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " Motion capture file is corrupted! " << endl;
      fclose(file);
      return false;
    }

    for (unsigned int y = 0; y < _bones.size() - 1; y++)
    {
      float qFloats[4];
      float tFloats[3];
      float scale;
      fread((void*)&qFloats[0], sizeof(float), 1, file);
      fread((void*)&qFloats[1], sizeof(float), 1, file);
      fread((void*)&qFloats[2], sizeof(float), 1, file);
      fread((void*)&qFloats[3], sizeof(float), 1, file);
        
      fread((void*)&tFloats[0], sizeof(float), 1, file);
      fread((void*)&tFloats[1], sizeof(float), 1, file);
      fread((void*)&tFloats[2], sizeof(float), 1, file);
      fread((void*)&scale, sizeof(float), 1, file);

      // packing is slightly different from Pinocchio -- real component goes last here
      _mocapRotations[y].push_back(QUATERNION(qFloats[1], qFloats[2], qFloats[3], qFloats[0]));
      _mocapTranslations[y].push_back(VEC3F(tFloats[0], tFloats[1], tFloats[2]));
      _mocapScalings[y].push_back(scale);
    }
  }

  fclose(file);
  cout << " done. " << endl;

  return true;
}

//////////////////////////////////////////////////////////////////////
// load up a Pinocchio frame
//////////////////////////////////////////////////////////////////////
void SKELETON::loadMocapFrame(int frame)
{
  // scaling factor for tet mesh
  //Real scale = _tetMesh->meshScaling();
  Real scale = 1.0;
  //Real scale = 2.0;

  assert(frame >= 0);
  frame = frame % _totalFrames;

  // get the total body translation
  QUATERNION q = _mocapRotations[0][frame];
  VEC3F t = _mocapTranslations[0][frame];
  MATRIX3 rotation = q.toExplicitMatrix3x3();
  VEC3F totalTranslation = rotation * _bones[0]->globalRestPosition() + t;

  for (unsigned int x = 0; x < _bones.size(); x++)
  {
    int which = (x == 0) ? 0 : x - 1;

    QUATERNION q = _mocapRotations[which][frame];
    VEC3F t = _mocapTranslations[which][frame];
    MATRIX3 rotation = q.toExplicitMatrix3x3();

    _bones[x]->quaternion() = q;
    //_bones[x]->translation() = scale * t;
    //_bones[x]->globalPosition() = rotation * (scale * _bones[x]->globalRestPosition()) + t;
    _bones[x]->translation() = t;
    _bones[x]->globalPosition() = rotation * (_bones[x]->globalRestPosition()) + t;

    // apply the scaling
    _bones[x]->globalPosition() -= totalTranslation;
    _bones[x]->globalPosition() *= scale;
    _bones[x]->globalPosition() += scale * totalTranslation;
  }

  /*
  for (unsigned int x = 0; x < _bones.size(); x++)
  {
    int which = (x == 0) ? 0 : x - 1;

    QUATERNION q = _mocapRotations[which][frame];
    VEC3F t = _mocapTranslations[which][frame];

    MATRIX3 rotation = q.toExplicitMatrix3x3();
    _bones[x]->rotation() = rotation;
    _bones[x]->translation() = t;
    _bones[x]->globalPosition() = rotation * _bones[x]->globalRestPosition() + t;
  }
  */

  if (!_skinnedObj || !_originalObj)
    return;

  // apply to the OBJ vertices
  assert(_skinnedObj);
  assert(_originalObj);
  vector<VEC3F>& newVertices = _skinnedObj->vertices; 
  vector<VEC3F>& oldVertices = _originalObj->vertices;
  assert(newVertices.size() == oldVertices.size());

  for (unsigned int x = 0; x < oldVertices.size(); x++)
  {
    VEC3F original = oldVertices[x];
    newVertices[x] *= 0;

    for (unsigned int y = 1; y < _bones.size(); y++)
    {
      VEC3F transformed = _bones[y]->rotation() * original + _bones[y]->translation();

      // weight for bone 0 is always zero, so shift everything over one
      transformed *= _objSkinningWeights[x][y - 1];
      newVertices[x] += transformed;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// load up a Pinocchio frame
//////////////////////////////////////////////////////////////////////
void SKELETON::loadInterpolatedMocapFrame(float frame)
{
  // scaling factor for tet mesh
  //Real scale = _tetMesh->meshScaling();
  Real scale = 1.0;

  int firstFrame = (int)frame;
  Real blend = frame - firstFrame;

  // get the total body translation
  QUATERNION q0 = _mocapRotations[0][firstFrame];
  VEC3F t0 = _mocapTranslations[0][firstFrame];
  QUATERNION q1 = _mocapRotations[0][firstFrame+1];
  VEC3F t1 = _mocapTranslations[0][firstFrame+1];

  VEC3F tBlended = (1.0 - blend) * t0 + blend * t1;
  QUATERNION qBlended = (1.0 - blend) * q0 + blend * q1;
  qBlended.normalize();

  MATRIX3 rotation = qBlended.toExplicitMatrix3x3();
  VEC3F totalTranslation = rotation * _bones[0]->globalRestPosition() + tBlended;

  for (unsigned int x = 0; x < _bones.size(); x++)
  {
    int which = (x == 0) ? 0 : x - 1;

    QUATERNION q0 = _mocapRotations[which][firstFrame];
    VEC3F t0 = _mocapTranslations[which][firstFrame];
    QUATERNION q1 = _mocapRotations[which][firstFrame + 1];
    VEC3F t1 = _mocapTranslations[which][firstFrame + 1];
    VEC3F t = (1.0 - blend) * t0 + blend * t1;
    QUATERNION q = (1.0 - blend) * q0 + blend * q1;
    q.normalize();

    MATRIX3 rotation = q.toExplicitMatrix3x3();

    _bones[x]->quaternion() = q;
    _bones[x]->rotation() = rotation;
    _bones[x]->translation() = t;
    _bones[x]->globalPosition() = rotation * (_bones[x]->globalRestPosition()) + t;

    // apply the scaling
    _bones[x]->globalPosition() -= totalTranslation;
    _bones[x]->globalPosition() *= scale;
    _bones[x]->globalPosition() += scale * totalTranslation;
  }
}

//////////////////////////////////////////////////////////////////////
// load up a Pinocchio frame
//////////////////////////////////////////////////////////////////////
void SKELETON::loadMocapFrame(float time, float skinningDt)
{
  // scaling factor for tet mesh
  Real scale = _tetMesh->meshScaling();

  int frame = (int)(time / skinningDt) % _totalFrames;

  Real blend = (time / skinningDt) - frame;

  // get the total body translation
  QUATERNION q0 = _mocapRotations[0][frame];
  VEC3F t0 = _mocapTranslations[0][frame];
  QUATERNION q1 = _mocapRotations[0][frame+1];
  VEC3F t1 = _mocapTranslations[0][frame+1];

  VEC3F t = (1.0 - blend) * t0 + blend * t1;
  QUATERNION q = (1.0 - blend) * q0 + blend * q1;
  q.normalize();

  MATRIX3 rotation = q.toExplicitMatrix3x3();
  VEC3F totalTranslation = rotation * _bones[0]->globalRestPosition() + t;

  for (unsigned int x = 0; x < _bones.size(); x++)
  {
    int which = (x == 0) ? 0 : x - 1;

    QUATERNION q0 = _mocapRotations[which][frame];
    VEC3F t0 = _mocapTranslations[which][frame];
    QUATERNION q1 = _mocapRotations[which][frame + 1];
    VEC3F t1 = _mocapTranslations[which][frame + 1];
    VEC3F t = (1.0 - blend) * t0 + blend * t1;
    QUATERNION q = (1.0 - blend) * q0 + blend * q1;
    q.normalize();

    MATRIX3 rotation = q.toExplicitMatrix3x3();

    _bones[x]->rotation() = rotation;
    _bones[x]->translation() = t;
    _bones[x]->globalPosition() = rotation * (_bones[x]->globalRestPosition()) + t;

    // apply the scaling
    _bones[x]->globalPosition() -= totalTranslation;
    _bones[x]->globalPosition() *= scale;
    _bones[x]->globalPosition() += scale * totalTranslation;
  }
}

//////////////////////////////////////////////////////////////////////////////
// Load in the attachment weights from Pinocchio output
//////////////////////////////////////////////////////////////////////////////
void SKELETON::loadPinocchioSkinning(OBJ* obj, string attachmentFile)
{
  _skinnedObj = obj;
  _originalObj = new OBJ(*obj);

  // find out how many bones there are
  int totalBones = _bones.size() - 1;

  // stomp old weights
  _objSkinningWeights.clear();

  // read in the weights
  FILE* file = fopen(attachmentFile.c_str(), "r");
  if (file == NULL)
  {
    cout << " No skinning file " << attachmentFile.c_str() << " found!" << endl;
    return;
  }

  vector<VEC3F>& vertices = obj->vertices;
  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    if (feof(file))
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " Number of weights in skinning and in OBJ do not match! " << endl;
      exit(0);
    }

    vector<float> weights;
    for (int y = 0; y < totalBones; y++)
    {
      float weight = 0.0;
      fscanf(file, "%f", &weight);
      weights.push_back(weight);
    }
    _objSkinningWeights.push_back(weights);
    fscanf(file, "\n");

    if (x == 0)
    {
      cout << " first weights read in: " << endl;
      for (int y = 0; y < totalBones; y++)
        cout << _objSkinningWeights[0][y] << " ";
      cout << endl;
    }
  }
  fclose(file);
}

//////////////////////////////////////////////////////////////////////////////
// draw the Pinocchio skinning weights
//////////////////////////////////////////////////////////////////////////////
void SKELETON::drawObjWeights(int bone)
{
  if (_skinnedObj->normals.size() == 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Need to call ComputeVertexNormals before calling draw! " << endl;
    return;
  }

  glBegin(GL_TRIANGLES);
  vector<OBJ::Face>& faces = _skinnedObj->faces;
  for(unsigned int i = 0; i < faces.size(); i++)
  {
    vector<int>& verts = faces[i].vertices;

    float weights[3];
    weights[0] = _objSkinningWeights[verts[0]][bone];
    weights[1] = _objSkinningWeights[verts[1]][bone];
    weights[2] = _objSkinningWeights[verts[2]][bone];

    weights[0] *= 3;
    weights[1] *= 3;
    weights[2] *= 3;

    for (int y = 0; y < 3; y++)
      weights[y] = (weights[y] > 1.0) ? 1.0 : weights[y];

    glColor4f(weights[0], 0.0, 1.0 - weights[0],1);
    _skinnedObj->drawVert(verts[0]);
    glColor4f(weights[1], 0.0, 1.0 - weights[1],1);
    _skinnedObj->drawVert(verts[1]);
    glColor4f(weights[2], 0.0, 1.0 - weights[2],1);
    _skinnedObj->drawVert(verts[2]);
  }
  glEnd();
}

//////////////////////////////////////////////////////////////////////////////
// draw the Pinocchio skinning weights
//////////////////////////////////////////////////////////////////////////////
void SKELETON::drawTetWeights(int bone)
{
  vector<pair<int, int> >& surfaceFaces = _tetMesh->surfaceFaces();
  vector<TET>& tets = _tetMesh->tets();

  for (unsigned int x = 0; x < surfaceFaces.size(); x++)
  {
    TRIANGLE triangle = tets[surfaceFaces[x].first].face(surfaceFaces[x].second);
    VEC3F* vertices[3];
    int vertexIDs[3];
    for (int y = 0; y < 3; y++)
    {
      vertices[y] = triangle.vertex(y);
      vertexIDs[y] = _tetMesh->vertexID(vertices[y]);
    }

    // get the weights
    float weights[3];
    for (int y = 0; y < 3; y++)
    {
      weights[y] = 0;
      vector<pair<int, Real> >& weightVector = _skinning[vertexIDs[y]];
      for (unsigned int z = 0; z < weightVector.size(); z++)
      {
        if (weightVector[z].first == bone)
          weights[y] = weightVector[z].second;
      }

      weights[y] *= 3.0;

      weights[y] = (weights[y] > 1.0) ? 1.0 : weights[y];
    }

    VEC3F normal = cross(*vertices[1] - *vertices[0], 
                         *vertices[2] - *vertices[0]);
    normal.normalize();

    glBegin(GL_TRIANGLES);
      glNormal3d(normal[0], normal[1], normal[2]);
      //glColor4f(weights[0], 0.0, 1.0 - weights[0],1);
      glColor4f(weights[0], weights[0], weights[0], 1);
      glVertex3dv(*vertices[0]);
      //glColor4f(weights[1], 0.0, 1.0 - weights[1],1);
      glColor4f(weights[1], weights[1], weights[1], 1);
      glVertex3dv(*vertices[1]);
      //glColor4f(weights[2], 0.0, 1.0 - weights[2],1);
      glColor4f(weights[2], weights[2], weights[2], 1);
      glVertex3dv(*vertices[2]);
    glEnd();
  }
}

//////////////////////////////////////////////////////////////////////////////
// draw all the Pinocchio skinning weights
//////////////////////////////////////////////////////////////////////////////
void SKELETON::drawTetWeights()
{
  vector<pair<int, int> >& surfaceFaces = _tetMesh->surfaceFaces();
  vector<TET>& tets = _tetMesh->tets();

  MERSENNETWISTER twister(9876);

  vector<VEC3F> colors;
  for (unsigned int x = 0; x < _bones.size(); x++)
  {
    VEC3F color(twister.rand(), twister.rand(), twister.rand());
    colors.push_back(color);
  }

  for (unsigned int x = 0; x < surfaceFaces.size(); x++)
  {
    TRIANGLE triangle = tets[surfaceFaces[x].first].face(surfaceFaces[x].second);
    VEC3F* vertices[3];
    int vertexIDs[3];
    for (int y = 0; y < 3; y++)
    {
      vertices[y] = triangle.vertex(y);
      vertexIDs[y] = _tetMesh->vertexID(vertices[y]);
    }

    // get the weights
    VEC3F vertexColors[3];
    for (int y = 0; y < 3; y++)
    {
      vector<pair<int, Real> >& weightVector = _skinning[vertexIDs[y]];
      for (unsigned int z = 0; z < weightVector.size(); z++)
        vertexColors[y] += weightVector[z].second * colors[weightVector[z].first];
    }

    VEC3F normal = cross(*vertices[1] - *vertices[0], 
                         *vertices[2] - *vertices[0]);
    normal.normalize();

    glBegin(GL_TRIANGLES);
      glNormal3d(normal[0], normal[1], normal[2]);
      glColor4f(vertexColors[0][0], vertexColors[0][1], vertexColors[0][2], 1.0);
      glVertex3dv(*vertices[0]);
      glColor4f(vertexColors[1][0], vertexColors[1][1], vertexColors[1][2], 1.0);
      glVertex3dv(*vertices[1]);
      glColor4f(vertexColors[2][0], vertexColors[2][1], vertexColors[2][2], 1.0);
      glVertex3dv(*vertices[2]);
    glEnd();
  }
}

//////////////////////////////////////////////////////////////////////////////
// compute a Pinocchio skinning for the given tet mesh
//////////////////////////////////////////////////////////////////////////////
void SKELETON::buildOdeSkinning(TET_MESH* mesh)
{
  _tetMesh = mesh;

  // blow away any old skinning
  vector<VEC3F>& restPose = _tetMesh->restPose();
  _skinning.resize(restPose.size());

  // cache the bone line segments
  vector<pair<VEC3F, VEC3F> > boneSegments;
  for (unsigned int x = 0; x < _bones.size(); x++)
  {
    VEC3F& translation = _bones[x]->translation();
    MATRIX3 rotation = _bones[x]->quaternion().toExplicitMatrix3x3();
    Real length = _boneLengths[x];
    VEC3F beginVertex(0,0, length * 0.5);
    VEC3F endVertex(0,0, -length * 0.5);
    beginVertex = rotation * beginVertex + translation;
    endVertex = rotation * endVertex + translation;

    boneSegments.push_back(pair<VEC3F, VEC3F>(beginVertex, endVertex));
  }

  cout << " Building ODE skinning ..."; flush(cout);

  // each vertex just picks up the skinning of the nearest surface vertex
  for (unsigned int x = 0; x < restPose.size(); x++)
  {
    int closestBone = 0;
    Real closestDistance = 0;

    for (unsigned int y = 0; y < boneSegments.size(); y++)
    {
      VEC3F B = boneSegments[y].first;
      VEC3F P = restPose[x];
      VEC3F M = boneSegments[y].second - B;

      Real t0 = M * (P - B) / (M * M);

      Real distance = norm(P - (B + t0 * M));
      if (t0 <= 0)
        distance = norm(P - B);

      if (t0 >= 1)
        distance = norm(P - (B + M));

      if (distance < closestDistance || y == 0)
      {
        closestBone = y;
        closestDistance = distance;
      }
    }

    _skinning[x].clear();
    _skinning[x].push_back(pair<int, Real>(closestBone, 1.0));
  }
  cout << " done." << endl;
}

//////////////////////////////////////////////////////////////////////////////
// compute a Pinocchio skinning for the given tet mesh
//////////////////////////////////////////////////////////////////////////////
void SKELETON::buildPinocchioSkinning(TET_MESH* mesh)
{
  _tetMesh = mesh;

  // blow away any old skinning
  vector<VEC3F>& restPose = _tetMesh->restPose();
  _skinning.resize(restPose.size());

  if (_originalObj == NULL)
  {
    cout << " Need a OBJ to skin against! " << endl;
    return;
  }

  cout << " Building Pinocchio skinning ..."; flush(cout);
  vector<VEC3F>& objVertices = _originalObj->vertices;

  // each vertex just picks up the skinning of the nearest surface vertex
  for (unsigned int x = 0; x < restPose.size(); x++)
  {
    // doing a brute force search here -- can add a KD tree if this
    // gets to be very expensive
    //VEC3F diff = restPose[x] - objVertices[0];
    VEC3F diff = restPose[x] - objVertices[0];
    Real bestFound = norm2(diff);
    int bestIndex = 0;

    for (unsigned int y = 1; y < objVertices.size(); y++)
    {
      //diff = restPose[x] - objVertices[y];
      diff = restPose[x] - objVertices[y];
      Real current = norm2(diff);
      if (current < bestFound)
      {
        bestFound = current;
        bestIndex = y;
      }
    }

    // store the skinning of the closest OBJ vertex
    vector<float>& skinning = _objSkinningWeights[bestIndex];
    for (unsigned int y = 0; y < skinning.size(); y++)
    {
      // offset by one, because Pinocchio doesn't count the first bone
      if (skinning[y] > 0.0)
        _skinning[x].push_back(pair<int, Real>(y + 1, skinning[y]));
    }
  }
  cout << " done." << endl;
}

//////////////////////////////////////////////////////////////////////////////
// constrain vertices in tets that intersect the bones
//////////////////////////////////////////////////////////////////////////////
void SKELETON::constrainBoneTets(TET_MESH* tetMesh, string filename)
{
  if (tetMesh->constrained())
  {
    cout << " Tet mesh is already constrained! " << endl;
    cout << " Don't try to constrain along bones twice, it will create redundant, confusing meshes! " << endl;
    exit(0);
  }

  cout << " Constraining tets along bones ... "; flush(cout);
  _tetMesh = tetMesh;

  // vertices to constrain
  map<VEC3F*, bool> constrainedVertices;

  // scan all tets, intersect with bones
  vector<TET>& tets = _tetMesh->tets();
  for (unsigned int x = 0; x < tets.size(); x++)
  {
    bool hitTet = false;
    for (unsigned int y = 0; y < _bones.size(); y++)
    {
      BONE* bone = _bones[y];
      if (bone->parent() == NULL) continue;

      BONE* parent = bone->parent();

      // line segment to intersect against
      VEC3F start = bone->globalRestPosition();
      VEC3F end = parent->globalRestPosition();

      // intersect line segment with each face
      for (int z = 0; z < 4; z++)
      {
        TRIANGLE face = tets[x].face(z);
        if (face.intersects(start, end))
        {
          hitTet = true;
          z = 4;
          y = _bones.size();
        }
      }
    }

    // if the tet was hit, store its vertices
    if (hitTet)
      for (int y = 0; y < 4; y++)
        constrainedVertices[tets[x].vertices[y]] = true;
  }
  cout << " done." << endl;
  cout << " Found " << constrainedVertices.size() << " vertices to constrain. " << endl;

  // flatten out to a vector
  vector<VEC3F*> newConstraints;
  map<VEC3F*, bool>::iterator iter;
  for (iter = constrainedVertices.begin(); iter != constrainedVertices.end(); iter++)
    newConstraints.push_back(iter->first);

  _tetMesh->writeNewConstraints(newConstraints, filename.c_str());
}

//////////////////////////////////////////////////////////////////////////////
// constrain vertices in tets that intersect the bones
//////////////////////////////////////////////////////////////////////////////
void SKELETON::constrainOdeBoneTets(TET_MESH* tetMesh, string filename)
{
  if (tetMesh->constrained())
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Tet mesh is already constrained! " << endl;
    cout << " Don't try to constrain along bones twice, it will create redundant, confusing meshes! " << endl;
    exit(0);
  }

  cout << " Constraining tets along bones ... "; flush(cout);
  _tetMesh = tetMesh;

  // vertices to constrain
  map<VEC3F*, bool> constrainedVertices;

  // scan all tets, intersect with bones
  vector<TET>& tets = _tetMesh->tets();
  for (unsigned int x = 0; x < tets.size(); x++)
  {
    bool hitTet = false;
    for (unsigned int y = 0; y < _bones.size(); y++)
    {
      BONE* bone = _bones[y];

      VEC3F& translation = _bones[y]->translation();
      MATRIX3 rotation = _bones[y]->quaternion().toExplicitMatrix3x3();
      Real length = _boneLengths[y];

      VEC3F beginVertex(0,0, length * 0.5);
      VEC3F endVertex(0,0, -length * 0.5);

      VEC3F start = rotation * beginVertex + translation;
      VEC3F end = rotation * endVertex + translation;

      // intersect line segment with each face
      for (int z = 0; z < 4; z++)
      {
        TRIANGLE face = tets[x].face(z);
        if (face.intersects(start, end))
        {
          hitTet = true;
          z = 4;
          y = _bones.size();
        }
      }
    }

    // if the tet was hit, store its vertices
    if (hitTet)
      for (int y = 0; y < 4; y++)
        constrainedVertices[tets[x].vertices[y]] = true;
  }
  cout << " done." << endl;
  cout << " Found " << constrainedVertices.size() << " vertices to constrain. " << endl;

  // flatten out to a vector
  vector<VEC3F*> newConstraints;
  map<VEC3F*, bool>::iterator iter;
  for (iter = constrainedVertices.begin(); iter != constrainedVertices.end(); iter++)
    newConstraints.push_back(iter->first);

  cout << " Writing out new constrained mesh " << filename.c_str() << endl;
  _tetMesh->writeNewConstraints(newConstraints, filename.c_str());
}

//////////////////////////////////////////////////////////////////////////////
// apply the ODE skinning to the tet mesh
//////////////////////////////////////////////////////////////////////////////
void SKELETON::updateOdeSkinning(bool constrainedOnly)
{
  assert(_tetMesh);

  int start = (constrainedOnly) ? _tetMesh->unconstrainedNodes() : 0;
  vector<VEC3F>& vertices = _tetMesh->vertices();
  vector<VEC3F>& restPose = _tetMesh->restPose();

  // for each constrained vertex
  for (unsigned int x = start; x < vertices.size(); x++)
  {
    // get the skinning weights
    vector<pair<int, Real> >& weights = _skinning[x];

    VEC3F final;
    for (unsigned int y = 0; y < weights.size(); y++)
    {
      int bone = weights[y].first;
      Real weight = weights[y].second;

      MATRIX3 rotation = _bones[bone]->rotation() * _bones[bone]->rotationOriginal().transpose();
      VEC3F transformed = restPose[x] - _bones[bone]->translationOriginal();
      transformed = rotation * transformed;
      transformed += _bones[bone]->translation();
      transformed *= weight;
      final += transformed;
    }

    vertices[x] = final;
  }

  // compute the displacement vector
  _tetMesh->recoverX();
}

//////////////////////////////////////////////////////////////////////////////
// apply the Pinocchio skinning to the tet mesh
//////////////////////////////////////////////////////////////////////////////
void SKELETON::updatePinocchioSkinning(bool constrainedOnly)
{
  assert(_tetMesh);

  int start = (constrainedOnly) ? _tetMesh->unconstrainedNodes() : 0;
  vector<VEC3F>& vertices = _tetMesh->vertices();
  vector<VEC3F>& restPose = _tetMesh->restPose();

#if 1
  Real scale = _tetMesh->meshScaling();

  // for each constrained vertex
  for (unsigned int x = start; x < vertices.size(); x++)
  {
    // get the skinning weights
    vector<pair<int, Real> >& weights = _skinning[x];

    VEC3F final;
    for (unsigned int y = 0; y < weights.size(); y++)
    {
      int bone = weights[y].first;
      Real weight = weights[y].second;

      //VEC3F transformed = _bones[bone]->rotation() * restPose[x] + _bones[bone]->translation();
      VEC3F transformed = _bones[bone]->rotation() * restPose[x] + scale * _bones[bone]->translation();
      transformed *= weight;
      final += transformed;
    }

    vertices[x] = final;
  }
#else
  // for each vertex
  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    //VEC3F final = _bones[0]->rotation() * restPose[x] + _bones[0]->translation();

    VEC3F final = vertices[x];
    final[2] -= 0.1;

    vertices[x] = final;
  }
#endif

  // compute the displacement vector
  _tetMesh->recoverX();
}

//////////////////////////////////////////////////////////////////////////////
// write the current tet mesh to Metis weighted by the skinning
// weights
//////////////////////////////////////////////////////////////////////////////
void SKELETON::writeWeightedMetis(string filename)
{
  if (_tetMesh == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Need to load a tet mesh into the skeleton before it can do anything! " << endl;
    exit(0);
  }

  vector<TET>& tets = _tetMesh->tets();

  // compute the total number of arcs (i.e. 2x the number of edges)
  cout << " Building arcs ..."; flush(cout);
  vector<vector<int> > arcs;
  if (!readOneRings(arcs))
  {
    for (unsigned int x = 0; x < tets.size(); x++)
    {
      vector<int> oneRing;
      _tetMesh->tetOneRing(x, oneRing);
      arcs.push_back(oneRing);
    }
    writeOneRings(arcs);
  }
  cout << " done." << endl;

  // compute the unique edges
  map<pair<int, int>, int> uniqueEdges;
  for (unsigned int x = 0; x < tets.size(); x++)
  {
    vector<int>& ring = arcs[x];

    for (unsigned int y = 0; y < ring.size(); y++)
    {
      pair<int,int> toHash(x, ring[y]);
      if (x < ring[y])
        toHash = pair<int,int>(ring[y], x);
      uniqueEdges[toHash]++;
    }
  }
  int totalEdges = uniqueEdges.size();

  vector<float> maxWeights(_bones.size());
  for (unsigned int x = 0; x < maxWeights.size(); x++)
    maxWeights[x] = 0;
  for (unsigned int x = 0; x < _skinning.size(); x++)
    for (unsigned int y = 0; y < _skinning[x].size(); y++)
    {
      int boneID = _skinning[x][y].first;
      float weight = _skinning[x][y].second;

      if (weight > maxWeights[boneID])
        maxWeights[boneID] = weight;
    }

  FILE* file = fopen(filename.c_str(), "w");
  //fprintf(file, "%i %i 1\n", (int)(tets.size()), totalEdges);
  fprintf(file, "%i %i 0\n", (int)(tets.size()), totalEdges);

  cout << " Writing arcs ..."; flush(cout);
  for (unsigned int x = 0; x < tets.size(); x++)
  {
    vector<int>& oneRing = arcs[x];

    /*
    bool bigFound = false;
    for (int y = 0; y < 4; y++)
    {
      // get the weight list of the current vertex
      int vertexID = _tetMesh->vertexID(tets[x].vertices[y]);
      vector<pair<int, Real> >& weights = _skinning[vertexID];

      // see if there's a big one
      for (unsigned int z = 0; z < weights.size(); z++)
        if (weights[z].second > 0.9 * maxWeights[weights[z].first] ||
            weights[z].second > 0.9)
        //if (weights[z].second > 0.95 * maxWeights[weights[z].first])
          bigFound = true;
    }

    //if (oneRing.size() == 0) continue;
    for (unsigned int y = 0; y < oneRing.size(); y++)
    {
      //int weight = (bigFound) ? 100 : 1;
      int weight = (bigFound) ? 1 : 1;
      //fprintf(file, "%i %i ", oneRing[y], weight);
      //
      if (y != 0)
        fprintf(file, " ");

      fprintf(file, "%i", oneRing[y]);
    }
    fprintf(file, "\n");
    */
  }
  cout << " done." << endl;

  fclose(file);
}

//////////////////////////////////////////////////////////////////////////////
// write the current tet mesh to scotch weighted by the skinning
// weights
//////////////////////////////////////////////////////////////////////////////
void SKELETON::writeWeightedScotch(string filename)
{
  if (_tetMesh == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Need to load a tet mesh into the skeleton before it can do anything! " << endl;
    exit(0);
  }

  vector<TET>& tets = _tetMesh->tets();

  // compute the total number of arcs (i.e. 2x the number of edges)
  cout << " Building arcs ..."; flush(cout);
  int totalArcs = 0;
  vector<vector<int> > arcs;
  if (!readOneRings(arcs))
  {
    for (unsigned int x = 0; x < tets.size(); x++)
    {
      vector<int> oneRing;
      _tetMesh->tetOneRing(x, oneRing);
      totalArcs += oneRing.size();

      arcs.push_back(oneRing);
    }
    writeOneRings(arcs);
  }
  else
  {
    for (unsigned int x = 0; x < tets.size(); x++)
      totalArcs += arcs[x].size();
  }
  cout << " done." << endl;

  vector<float> maxWeights(_bones.size());
  for (unsigned int x = 0; x < maxWeights.size(); x++)
    maxWeights[x] = 0;
  for (unsigned int x = 0; x < _skinning.size(); x++)
    for (unsigned int y = 0; y < _skinning[x].size(); y++)
    {
      int boneID = _skinning[x][y].first;
      float weight = _skinning[x][y].second;

      if (weight > maxWeights[boneID])
        maxWeights[boneID] = weight;
    }

  FILE* file = fopen(filename.c_str(), "w");
  fprintf(file, "0\n");
  fprintf(file, "%i\t%i\n", (int)(tets.size()), totalArcs);
  fprintf(file, "0\t010\n");

  cout << " Writing arcs ..."; flush(cout);
  for (unsigned int x = 0; x < tets.size(); x++)
  {
    vector<int>& oneRing = arcs[x];
    fprintf(file, "%i", (int)(oneRing.size()));

    bool bigFound = false;
    for (int y = 0; y < 4; y++)
    {
      // get the weight list of the current vertex
      int vertexID = _tetMesh->vertexID(tets[x].vertices[y]);
      vector<pair<int, Real> >& weights = _skinning[vertexID];

      // see if there's a big one
      for (unsigned int z = 0; z < weights.size(); z++)
        //if (weights[z].second > 0.9 * maxWeights[weights[z].first] || weights[z].second > 0.9)
        if (weights[z].second > 0.95 * maxWeights[weights[z].first])
          bigFound = true;
    }

    for (unsigned int y = 0; y < oneRing.size(); y++)
    {
      //int weight = (bigFound) ? 100 : 1;
      int weight = (bigFound) ? 100 : 1;
      fprintf(file, "\t %i %i", weight, oneRing[y]);
    }
    fprintf(file, "\n");
  }
  cout << " done." << endl;

  fclose(file);

  /*
  cout << " Computing edge weights ... "; flush(cout);

  // compute the largest weights in each tet
  vector<float> largestWeights;
  float largestOverall = 0;
  for (unsigned int x = 0; x < tets.size(); x++)
  {
    // get the largest weight of the vertices in this tet
    TET& tet = tets[x];
    float largestFound = 0.0;
    for (int y = 0; y < 4; y++)
    {
      // get the weight list of the current vertex
      int vertexID = _tetMesh->vertexID(tet.vertices[y]);
      vector<pair<int, Real> >& weights = _skinning[vertexID];

      // see if there's a big one
      for (unsigned int z = 0; z < weights.size(); z++)
        if (weights[z].second > largestFound)
          largestFound = weights[z].second;
    }
    largestWeights.push_back(largestFound);

    if (largestFound > largestOverall)
      largestOverall = largestFound;
  }

  // compute the edge weights
  vector<vector<float> > edgeWeights;
  for (unsigned int x = 0; x < tets.size(); x++)
  {
    vector<float> maxWeights;
    vector<int>& oneRing = arcs[x];

    // take the larger max weight of the center tet and its neighbor
    for (unsigned int y = 0; y < oneRing.size(); y++)
    {
      float larger = largestWeights[x];
      if (largestWeights[oneRing[y]] > larger)
        larger = largestWeights[oneRing[y]];

      maxWeights.push_back(larger);
    }
    edgeWeights.push_back(maxWeights);
  }
  cout << " done." << endl;

  FILE* file = fopen(filename.c_str(), "w");
  fprintf(file, "0\n");
  fprintf(file, "%i\t%i\n", (int)(tets.size()), totalArcs);
  fprintf(file, "0\t010\n");
  //fprintf(file, "0\t000\n");

  cout << " Writing arcs ..."; flush(cout);
  for (unsigned int x = 0; x < tets.size(); x++)
  {
    vector<int>& oneRing = arcs[x];
    vector<float>& weights = edgeWeights[x];
    fprintf(file, "%i", (int)(oneRing.size()));
    for (unsigned int y = 0; y < oneRing.size(); y++)
    {
      float fraction = weights[y] / largestOverall;

      int weight = (1000 * (pow(fraction, 3.0f)));
      fprintf(file, "\t %i %i", weight, oneRing[y]);
    }
    fprintf(file, "\n");
  }
  cout << " done." << endl;

  fclose(file);
  */
}

//////////////////////////////////////////////////////////////////////////////
// write out the one rings for the tet mesh
//////////////////////////////////////////////////////////////////////////////
void SKELETON::writeOneRings(vector<vector<int> >& rings)
{
  assert(_tetMesh);

  string tetMeshFilename = _tetMesh->filename();
  string filename = tetMeshFilename;
  filename += string(".onerings");
  FILE* file = fopen(filename.c_str(), "wb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Couldn't write out one rings! " << endl;
    return;
  }

  cout << " Writing out one rings: " << filename.c_str() << endl;

  int totalRings = rings.size();
  fwrite((void*)&totalRings, sizeof(int), 1, file);

  for (unsigned int x = 0; x < rings.size(); x++)
  {
    vector<int>& ring = rings[x];
    int oneRingSize = ring.size();
    fwrite((void*)&oneRingSize, sizeof(int), 1, file);
    for (unsigned int y = 0; y < ring.size(); y++)
      fwrite((void*)&ring[y], sizeof(int), 1, file);
  }

  fclose(file);
}
  
//////////////////////////////////////////////////////////////////////////////
// write out the one rings for the tet mesh
//////////////////////////////////////////////////////////////////////////////
bool SKELETON::readOneRings(vector<vector<int> >& rings)
{
  assert(_tetMesh);

  //cout << " Looking for one ring file ..."; flush(cout);
  string tetMeshFilename = _tetMesh->filename();
  string filename = tetMeshFilename;
  filename += string(".onerings");
  FILE* file = fopen(filename.c_str(), "rb");
  if (file == NULL)
  {
    //cout << " none found." << endl;
    return false;
  }
  //cout << " found!" << endl;

  int totalRings;
  fread((void*)&totalRings, sizeof(int), 1, file);

  rings.clear();
  for (int x = 0; x < totalRings; x++)
  {
    vector<int> ring;
    int oneRingSize;
    fread((void*)&oneRingSize, sizeof(int), 1, file);
    for (int y = 0; y < oneRingSize; y++)
    {
      int index;
      fread((void*)&index, sizeof(int), 1, file);
      ring.push_back(index);
    }
    rings.push_back(ring);
  }
  fclose(file);

  return true;
}

//////////////////////////////////////////////////////////////////////
// partition based on the ODE bone weights -- 
// there should only be one per tet
//////////////////////////////////////////////////////////////////////
void SKELETON::odeSkinningPartition()
{
  vector<TET>& tets = _tetMesh->tets();

  for (unsigned int x = 0; x < tets.size(); x++)
  {
    TET& tet = tets[x];
    map<int, int> partitionVotes;

    for (int y = 0; y < 4; y++)
    {
      // get the weight list of the current vertex
      int vertexID = _tetMesh->vertexID(tet.vertices[y]);
      vector<pair<int, Real> >& weights = _skinning[vertexID];

      partitionVotes[weights[0].first]++;
    }

    int winningPartition = 0;
    int mostVotes = -1;
    map<int,int>::iterator iter;
    for (iter = partitionVotes.begin(); iter != partitionVotes.end(); iter++)
      if (iter->second > mostVotes)
      {
        mostVotes = iter->second;
        winningPartition = iter->first;
      } 

    tet.partition = winningPartition;
  }
}

//////////////////////////////////////////////////////////////////////
// partition based on the maximum skinning weights
//////////////////////////////////////////////////////////////////////
void SKELETON::maxSkinningPartition()
{
  vector<TET>& tets = _tetMesh->tets();
  vector<float> maxWeights(_bones.size());
  for (unsigned int x = 0; x < maxWeights.size(); x++)
    maxWeights[x] = 0;

  for (unsigned int x = 0; x < _skinning.size(); x++)
    for (unsigned int y = 0; y < _skinning[x].size(); y++)
    {
      int boneID = _skinning[x][y].first;
      float weight = _skinning[x][y].second;

      if (weight > maxWeights[boneID])
        maxWeights[boneID] = weight;
    }

  for (unsigned int x = 0; x < tets.size(); x++)
  {
    VEC3F center = tets[x].center();

    VEC3F finalColor;

    TET& tet = tets[x];
    float maxFound = 0;
    int maxPartition = 1;
    for (int y = 0; y < 4; y++)
    {
      // get the weight list of the current vertex
      int vertexID = _tetMesh->vertexID(tet.vertices[y]);
      vector<pair<int, Real> >& weights = _skinning[vertexID];

      for (unsigned int z = 0; z < weights.size(); z++)
        if (weights[z].second / maxWeights[weights[z].first] > maxFound)
        {
          maxFound = weights[z].second / maxWeights[weights[z].first];
          maxPartition = weights[z].first;
        }
    }

    tet.partition = maxPartition - 1;
  }
}

//////////////////////////////////////////////////////////////////////
// apply the Pinocchio skinning to a partitioned tet mesh
//////////////////////////////////////////////////////////////////////
void SKELETON::updatePinocchioSkinning(PARTITIONED_TET_MESH* partitionedMesh)
{
  // apply transform to the raw vertices
  //for (unsigned int x = 0; x < _bones.size(); x++)
  for (unsigned int x = 0; x < _bones.size() - 1; x++)
  {
    TET_MESH* mesh = partitionedMesh->mesh(x);
    vector<VEC3F>& vertices = mesh->vertices();
    vector<VEC3F>& restPose = mesh->restPose();

    int whichBone = x + 1;

    MATRIX3 rotation = _bones[whichBone]->rotation();
    VEC3F translation = _bones[whichBone]->translation();

    VEC3F& centerOfMass = mesh->centerOfMass();
    VEC3F& boneRest = _bones[whichBone]->globalRestPosition();

    VEC3F domainTranslation = rotation * centerOfMass + translation;

    for (unsigned int y = 0; y < vertices.size(); y++)
      vertices[y] = rotation * (restPose[y] - centerOfMass) + domainTranslation;
      //vertices[y] = rotation * restPose[y] + translation;
      //vertices[y] = restPose[y] - centerOfMass + translation;
      //vertices[y] = restPose[y] - boneRest + _bones[x]->globalPosition();
  }
}

//////////////////////////////////////////////////////////////////////
// load up the orientations of an ODE file
//
// same as constructor, but assumes all the bones have already been
// allocated, and that there will be rotations and bone dimensions
//////////////////////////////////////////////////////////////////////
void SKELETON::loadOdeFrame(const char* filename)
{
  // try to open the file
  FILE* file = fopen(filename, "r");
  if (file == NULL)
  {
    cout << " No skeleton file " << filename << " found!" << endl;
    return;
  }
  cout << " Loading skeleton file " << filename << endl;

  // read in each bone
  while (!feof(file))
  {
    int index;
    float trans[3];
    int previousIndex;
    fscanf(file,"%i %f %f %f %i\n", &index, &trans[0], &trans[1], &trans[2], &previousIndex);

    QUATERNION rotation;
    float quat[4];
    fscanf(file,"%f %f %f %f\n", &quat[0], &quat[1], &quat[2], &quat[3]);
    rotation[0] = quat[0]; rotation[1] = quat[1]; 
    rotation[2] = quat[2]; rotation[3] = quat[3];

    float length, radius;
    fscanf(file,"%f %f\n", &length, &radius);

    VEC3F translation(trans[0], trans[1], trans[2]);
    _bones[index]->quaternion() = rotation;
    _bones[index]->translation() = translation;
    _bones[index]->globalPosition() = translation;
  }

  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// draw the ODE version with cylinders
//////////////////////////////////////////////////////////////////////
void SKELETON::drawOdeCylinders()
{
  // polygon resolution for capsule bodies
  int CAPSULE_SLICES = 16;
  int CAPSULE_STACKS = 12;

  for (unsigned int x = 0; x < _bones.size(); x++)
  {
    VEC3F p = _bones[x]->translation();
    float length = _boneLengths[x];
    float radius = _boneRadii[x];

    glPushMatrix();
      VEC3F axis;
      Real angle;
      _bones[x]->quaternion().axisAngle(axis, angle);
      glTranslatef(p[0], p[1], p[2]);
      glRotatef(angle, axis[0], axis[1], axis[2]);

      Real cylHalfHeight = length / 2.0;
      glBegin(GL_QUAD_STRIP);
      for (int i = 0; i < CAPSULE_SLICES + 1; i++)
      {
        Real angle = i / float(CAPSULE_SLICES) * 2.0 * M_PI;
        Real ca = cos(angle);
        Real sa = sin(angle);
        glNormal3f(ca, sa, 0);
        glVertex3f(radius * ca, radius * sa, cylHalfHeight);
        glVertex3f(radius * ca, radius * sa, -cylHalfHeight);
      }
      glEnd();
      glTranslated(0, 0, cylHalfHeight);
      glutSolidSphere(radius, CAPSULE_SLICES, CAPSULE_STACKS);
      glTranslated(0, 0, -2.0 * cylHalfHeight);
      glutSolidSphere(radius, CAPSULE_SLICES, CAPSULE_STACKS);
    glPopMatrix();
  }
}

//////////////////////////////////////////////////////////////////////
// draw the ODE version of bones
//////////////////////////////////////////////////////////////////////
void SKELETON::drawOdeBones()
{
  for (unsigned int x = 0; x < _bones.size(); x++)
  {
    VEC3F& translation = _bones[x]->translation();
    MATRIX3 rotation = _bones[x]->quaternion().toExplicitMatrix3x3();
    Real length = _boneLengths[x];

    VEC3F beginVertex(0,0, length * 0.5);
    VEC3F endVertex(0,0, -length * 0.5);
    VEC3F centerVertex(0,0,0);

    beginVertex = rotation * beginVertex + translation;
    endVertex = rotation * endVertex + translation;
    centerVertex = rotation * centerVertex + translation;

    glPointSize(10.0f);
    glColor4f(10.0, 0.0, 0.0, 1.0);
    glBegin(GL_POINTS);
      glVertex3f(beginVertex[0], beginVertex[1], beginVertex[2]);
      glVertex3f(endVertex[0], endVertex[1], endVertex[2]);
    glEnd();

    glLineWidth(3.0f);
    glColor4f(0.0, 10.0, 0.0, 1.0);
    glBegin(GL_LINES);
      glVertex3f(beginVertex[0], beginVertex[1], beginVertex[2]);
      glVertex3f(endVertex[0], endVertex[1], endVertex[2]);
    glEnd();

    glPointSize(10.0f);
    glColor4f(10.0, 10.0, 10.0, 1.0);
    glBegin(GL_POINTS);
      glVertex3f(centerVertex[0], centerVertex[1], centerVertex[2]);
    glEnd();
  }
}

//////////////////////////////////////////////////////////////////////
// read in joint angle limits from ODE file
//////////////////////////////////////////////////////////////////////
void SKELETON::readOdeJointLimits(const char* filename)
{
  FILE* file = fopen(filename, "r");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Could not open joints file " << filename << "!!!" << endl;
    exit(0);
  }
  cout << " Reading joint file: " << filename << endl;

  int totalJoints = 0;
  fscanf(file, "%i\n", &totalJoints);
  for (int x = 0; x < totalJoints; x++)
  {
    // retrieve data to write out
    VEC3F axis1;
    VEC3F axis2;

    pair<Real, Real> limits1;
    pair<Real, Real> limits2;

    int bodyIndex = 0;

    VEC3F anchor;
    MATRIX3 rotation;
    VEC3F translation;

    char type;
    fscanf(file, "%c\n", &type);

    if (type == 'u')
    {
      fscanf(file, "%lf %lf %lf\n", &axis1[0], &axis1[1], &axis1[2]);
      fscanf(file, "%lf %lf %lf\n", &axis2[0], &axis2[1], &axis2[2]);
      fscanf(file, "%lf %lf\n", &limits1.first, &limits1.second);
      fscanf(file, "%lf %lf\n", &limits2.first, &limits2.second);
      fscanf(file, "%i\n", &bodyIndex);
      fscanf(file, "%lf %lf %lf\n", &anchor[0], &anchor[1], &anchor[2]);
      fscanf(file, "%lf %lf %lf\n", &rotation(0,0), &rotation(0,1), &rotation(0,2));
      fscanf(file, "%lf %lf %lf\n", &rotation(1,0), &rotation(1,1), &rotation(1,2));
      fscanf(file, "%lf %lf %lf\n", &rotation(2,0), &rotation(2,1), &rotation(2,2));
      fscanf(file, "%lf %lf %lf\n", &translation[0], &translation[1], &translation[2]);

      _jointAxis1.push_back(axis1);
      _jointAxis2.push_back(axis2);
      _jointLimits1.push_back(limits1);
      _jointLimits2.push_back(limits2);
      _jointBody2.push_back(bodyIndex);
      _originalRotations.push_back(rotation);
      _originalTranslations.push_back(translation);
      _jointAnchors.push_back(anchor);
    }
    else
    {
      cout << " Unsupported joint type found: " << type << endl;
      exit(0);
    }
  }

  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// randomize the ODE joints
//////////////////////////////////////////////////////////////////////
void SKELETON::randomizeOdeJoints()
{
  static MERSENNETWISTER twister;

  for (unsigned int x = 0; x < _jointBody2.size(); x++)
  {
    // retrieve all the joint information
    VEC3F axis1 = _jointAxis1[x];
    VEC3F axis2 = _jointAxis2[x];

    pair<Real, Real> limits1 = _jointLimits1[x];
    pair<Real, Real> limits2 = _jointLimits2[x];

    Real range1 = limits1.second - limits1.first;
    Real range2 = limits2.second - limits2.first;

    VEC3F translation = _originalTranslations[x];
    VEC3F anchor = _jointAnchors[x];

    // create the random matrix
    Real rand1 = limits1.first + twister.rand() * range1;
    Real rand2 = limits2.first + twister.rand() * range2;

    MATRIX3 original = _originalRotations[x];
    MATRIX3 rotation1 = MATRIX3::rotation(axis1, rand1);
    MATRIX3 rotation2 = MATRIX3::rotation(rotation1 * axis2, rand2);

    MATRIX3 rotation = rotation2 * rotation1 * _originalRotations[x];

    VEC3F diff = translation - anchor;
    translation += rotation2 * rotation1 *diff;
    translation -= diff;

    int boneIndex = _jointBody2[x];
    _bones[boneIndex]->translation() = translation;
    _bones[boneIndex]->quaternion() = QUATERNION(rotation);

    /*
    // hand it to ODE
    dMatrix3 bodyRotation;
    matrix3ToOde(rotation, bodyRotation);
    dBodySetRotation(bodyID, bodyRotation);

    dBodySetPosition(bodyID, translation[0], translation[1], translation[2]);
    */
  }
}

//////////////////////////////////////////////////////////////////////
// cycle an ODE joint according to the time
//////////////////////////////////////////////////////////////////////
void SKELETON::cycleJoint(int jointID, Real time, Real period, bool flip)
{
  // retrieve all the joint information
  VEC3F axis1 = _jointAxis1[jointID];
  VEC3F axis2 = _jointAxis2[jointID];

  pair<Real, Real> limits1 = _jointLimits1[jointID];
  pair<Real, Real> limits2 = _jointLimits2[jointID];

  Real range1 = limits1.second - limits1.first;
  Real range2 = limits2.second - limits2.first;

  VEC3F translation = _originalTranslations[jointID];
  VEC3F anchor = _jointAnchors[jointID];

  // create the random matrix
  //Real cycleTime = time * M_PI  / 2;
  Real cycleTime = time * M_PI * period;

  Real sineTime = sin(cycleTime);

  Real angle1 = 0;
  Real angle2 = 0;

  if (flip)
    sineTime *= -1;

  if (sineTime < 0.0)
    angle1 = -sineTime * limits1.first;
  else
    angle1 = sineTime * limits1.second;

  sineTime *= -1;
  if (sineTime < 0.0)
    angle2 = -sineTime * limits2.first;
  else
    angle2 = sineTime * limits2.second;

  /*
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " sine: " << sineTime << endl;
  cout << " angles: " << angle1 << " " << angle2 << endl;
  cout << " limits 1: " << limits1.first << " " << limits1.second << endl;
  cout << " limits 2: " << limits2.first << " " << limits2.second << endl;
  */

  //Real angle1 = limits1.first + 0.5 * (-cos(cycleTime) + 1) * range1;
  //Real angle2 = limits2.first + 0.5 * (-cos(cycleTime) + 1) * range2;
  //Real angle1 = 0;
  //Real angle2 = 0;

  MATRIX3 original = _originalRotations[jointID];
  MATRIX3 rotation1 = MATRIX3::rotation(axis1, angle1);
  MATRIX3 rotation2 = MATRIX3::rotation(rotation1 * axis2, angle2);

  MATRIX3 rotation = rotation2 * rotation1 * _originalRotations[jointID];

  VEC3F diff = translation - anchor;
  translation += rotation2 * rotation1 *diff;
  translation -= diff;

  int boneIndex = _jointBody2[jointID];
  _bones[boneIndex]->translation() = translation;
  _bones[boneIndex]->quaternion() = QUATERNION(rotation);
}
