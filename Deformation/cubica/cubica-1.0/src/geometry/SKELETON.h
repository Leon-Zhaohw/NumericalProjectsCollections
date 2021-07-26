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
// SKELETON.h: interface for the SKELETON class.
//
//////////////////////////////////////////////////////////////////////

#ifndef SKELETON_H
#define SKELETON_H

#include <VEC3.h>
#include <MATRIX3.h>
#if WIN32
#include <gl/glut.h>
#elif USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif
#include <MATRIX.h>
#include <TET_MESH.h>
#include <PARTITIONED_TET_MESH.h>
#include "BONE.h"
#include <OBJ.h>
#include <QUATERNION.h>

//////////////////////////////////////////////////////////////////////
// skeleton for skinning
//////////////////////////////////////////////////////////////////////
class SKELETON {
public:
  SKELETON(TET_MESH* tetMesh, int bccRes);
  SKELETON(const string& filename, bool odeFile = false);
  virtual ~SKELETON();

  // draw
  void draw() { 
    glDisable(GL_DEPTH_TEST);
    drawBones(_rootBone);
    glEnable(GL_DEPTH_TEST);
  };
  void drawSkinning();

  // draw the ODE version with cylinders
  void drawOdeCylinders(); 

  // draw the ODE version of bones
  void drawOdeBones(); 

  // update the skeleton with the most current transformations
  void updateSkeleton() {
    VEC3F zero(0.0, 0.0, 0.0);
    MATRIX3 I = MATRIX3::I();
    updateBones(_rootBone, MATRIX(I, zero));
    updateSkinning();
  };

  // build some skeleton
  void buildSkeleton();
  void buildDragonSkeleton();
  void buildHeadSkeleton();

  // build some skinning
  void buildSkinning();
  void buildDragonSkinning();
  void buildHeadSkinning(); // NAZIS!!!

  // build Pinocchio skinning for a tet mesh --
  // just look for the nearest surface vertex and use those weights
  void buildPinocchioSkinning(TET_MESH* tetMesh);

  // build Ode skinning for a tet mesh --
  // just look for the nearest bone and use those weights
  void buildOdeSkinning(TET_MESH* tetMesh);

  // add the skinning displacements to the tet mesh, just for display
  void addSkinningDisplacements();
  
  // undo the skinning displacements to the tet mesh, for simulation
  void undoSkinningDisplacements();

  // cycle the selected bone
  void cycleSelection() { 
    _selected++;
    cout << " Selected bone " << _selected << endl;
    _selected = _selected % _bones.size();
  };

  // return the selected bone
  BONE* selected() { return  _bones[_selected]; };

  // return a specific bone
  BONE* bone(int x) { return _bones[x]; };

  // access all bones
  vector<BONE*>& bones() { return _bones; };

  // bone dimensions
  Real boneLength(int x) { return _boneLengths[x]; };
  Real boneRadius(int x) { return _boneRadii[x]; };

  vector<Real>& boneLengths() { return _boneLengths; };

  // where do you want to the base of the mesh?
  VEC3F& meshOrigin() { return _meshOrigin; };
  vector<MATRIX>& skinningTransforms() { return _skinningTransforms; };
  vector<MATRIX>& skinningTransformsOld() { return _skinningTransformsOld; };
  vector<vector<pair<int, Real> > >& skinning() { return _skinning; };
  int totalBones() { return _bones.size(); };
  int totalMocapFrames() { return _totalFrames; };
  OBJ* originalObj() { return _originalObj; };
  OBJ* skinnedObj()  { return _skinnedObj; };
  TET_MESH*& tetMesh() { return _tetMesh; };
  vector<pair<Real, Real> >& jointLimits1() { return _jointLimits1; };
  vector<pair<Real, Real> >& jointLimits2() { return _jointLimits2; };

  // get the array index of the transform
  int transformID(MATRIX* transform) { return _transformID[transform]; };
  int transformOldID(MATRIX* transform) { return _transformOldID[transform]; };

  // cycle rotations
  void rotateSelectedX(Real sign);
  void rotateSelectedY(Real sign);
  void rotateSelectedZ(Real sign);

  MATRIX3 headShake(float time, float speed, bool& forceFull);
  int headShakeKeyframed(float time);

  VEC3F headNormalize(VEC3F vertex, int res);

  MATRIX recreateHeadSkinning(int vertexID, float time);

  void reset();

  // read/write skinning
  bool readSkinning();
  bool writeSkinning();

  bool readState(string filename);
  bool writeState(string filename);

  // read in motion capture data from Pinocchio
  bool readMocap(string filename);

  // load up a Pinocchio frame
  void loadMocapFrame(int frame);
 
  // load up an interpolated Pinocchio frame
  void loadInterpolatedMocapFrame(float frame);

  // load up a Pinocchio frame
  void loadMocapFrame(float time, float skinningDt);

  // load up the orientations of an ODE file
  void loadOdeFrame(const char* filename);

  // read in a Pinocchio skinning
  void loadPinocchioSkinning(OBJ* obj, string attachmentFile);

  // draw the Pinocchio skinning weights
  void drawObjWeights(int bone);

  // draw the Pinocchio skinning weights
  void drawTetWeights(int bone);
  void drawTetWeights();

  // draw partitioning weights
  void drawPartitioningWeights();

  // constrain vertices in tets that intersect the bones
  // and write a new file
  void constrainBoneTets(TET_MESH* tetMesh, string filename);

  // constrain vertices in tets that intersect the bones
  // and write a new file
  // same as above, except bones are represented slightly
  // differently in ODE
  void constrainOdeBoneTets(TET_MESH* tetMesh, string filename);

  // apply the Pinocchio skinning to the tet mesh
  void updatePinocchioSkinning(bool constrainedOnly);
 
  // apply the Pinocchio skinning to a partitioned tet mesh
  void updatePinocchioSkinning(PARTITIONED_TET_MESH* partitionedMesh);

  // apply the ODE skinning to the tet mesh
  void updateOdeSkinning(bool constrainedOnly);

  // write the current tet mesh to Scotch weighted by the skinning
  // weights
  void writeWeightedScotch(string filename);

  // write the current tet mesh to Metis weighted by the skinning
  // weights
  void writeWeightedMetis(string filename);

  // partition based on the maximum skinning weights
  void maxSkinningPartition();
  
  // partition based on the ODE bone weights -- there should only be one per tet
  void odeSkinningPartition();

  // read in joint angle limits from ODE file
  void readOdeJointLimits(const char* filename);

  // randomize the ODE joints
  void randomizeOdeJoints();

  // cycle an ODE joint according to the time
  void cycleJoint(int jointID, Real time, Real period, bool flip = false);

private:
  // pointer to the bone tree
  BONE* _rootBone;

  // list of the bones
  vector<BONE*> _bones;

  // ODE specific --
  // bone lengths and radii
  vector<Real> _boneLengths;
  vector<Real> _boneRadii;

  // bone currently selection by the GUI
  int _selected;

  int _bccRes;

  // tet mesh to skin
  TET_MESH* _tetMesh;

  // where is the mesh's (0,0,0)?
  VEC3F _meshOrigin;

  // skinning bones and weights
  //
  // _skinning is of size restPose.size(). Each pair is
  // a bone index, and the Real is the skinning weight.
  vector<vector<pair<int, Real> > > _skinning;

  // current skinning transform matrices
  vector<MATRIX> _skinningTransforms;
  vector<MATRIX> _skinningTransformsOld;

  // map the _skinningTransforms and _skinningTransformsOld entries to indices
  map<MATRIX*, int> _transformID;
  map<MATRIX*, int> _transformOldID;
  
  // colors for the different bones
  float _colors[8][4];

  // motion capture data --
  // _mocap[x][y] is for bone x at timestep y
  int _totalFrames;
  vector<vector<QUATERNION> > _mocapRotations;
  vector<vector<VEC3F> > _mocapTranslations;
  vector<vector<Real> > _mocapScalings;

  // Pinocchio skinning weights
  //
  // _skinningWeights[x][y] is the weight for vertex x, bone y
  vector<vector<float> > _objSkinningWeights;

  // skinned OBJ file
  OBJ* _skinnedObj;

  // copy of original OBJ file
  OBJ* _originalObj;

  // recursively delete the skeleton
  void cleanup(BONE* bone);

  // recursively draw the skeleton
  void drawBones(BONE* bone);

  // recursively update the skeleton
  void updateBones(BONE* bone, MATRIX transform);

  // update the mesh based on the skeleton
  void updateSkinning();

  MATRIX3 multiRotation(int xTics, int yTics, int zTics, float fracPose);

  // write out the one rings for the tet mesh
  void writeOneRings(vector<vector<int> >& rings);
  
  // write out the one rings for the tet mesh
  bool readOneRings(vector<vector<int> >& rings);

  // hacks
  bool _first;
  bool _second;
  bool _third;
  bool _fourth;

  // joint angle limit data for universal ODE joints
  vector<VEC3F> _jointAxis1;
  vector<VEC3F> _jointAxis2;
  vector<pair<Real, Real> > _jointLimits1;
  vector<pair<Real, Real> > _jointLimits2;
  vector<int> _jointBody2;
  vector<MATRIX3> _originalRotations;
  vector<VEC3F> _originalTranslations;
  vector<VEC3F> _jointAnchors;
};

#endif
