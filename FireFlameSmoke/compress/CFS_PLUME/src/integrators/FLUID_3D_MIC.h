/*
 This file is part of SSFR (Zephyr).
 
 Zephyr is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 Zephyr is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Zephyr.  If not, see <http://www.gnu.org/licenses/>.
 
 Copyright 2013 Theodore Kim
 */
#ifndef FLUID_3D_MIC_H
#define FLUID_3D_MIC_H

#include "EIGEN.h"

#include <cstdlib>
#include <cmath>
#include <iostream>
#include "OBSTACLE.h"
#include "BOX.h"
#include "VEC3.h"
#include "FIELD_3D.h"
#include "VECTOR3_FIELD_3D.h"
#include "MERSENNETWISTER.h"
#include "SPARSE_MATRIX.h"

using namespace std;

class FLUID_3D_MIC  
{
public:
  FLUID_3D_MIC(int xRes, int yRes, int zRes, int amplify, unsigned int* boundaries = NULL);
  FLUID_3D_MIC();
  virtual ~FLUID_3D_MIC();

  void addSmokeColumn();
  void addSmokeTestCase(Real* field, VEC3I res);

  void addSmokeSphere();
  void addSmokeSphereTestCase(Real* field, VEC3I res);

  void setInitialVelocity(BOX* box);

  void step();
  void stepWithMovingObstacle(BOX* box);

  int xRes() const { return _xRes; };
  int yRes() const { return _yRes; };
  int zRes() const { return _zRes; };
  int slabSize() const { return _slabSize; };
  int size() const { return _xRes*_yRes*_zRes; };
  Real dt() const { return _dt; };
  Real dx() const { return _dx; };
  Real& vorticityEps() { return _vorticityEps; };
  Real& buoyancy() { return _buoyancy; };
  Real& wind()     { return _wind; };
    
  VEC3F  center() const {return _center;}	
  VEC3F  lengths() const {return _lengths;}	
    
  VEC3F  velocity(int i) const      { return _velocity(i); };
  VECTOR3_FIELD_3D& velocity()      { return _velocity; };
  VECTOR3_FIELD_3D& velocityOld()   { return _velocityOld; };
  VectorXd& velocityNeumann()       { return _velocityNeumann; };
  VECTOR3_FIELD_3D& force()         { return _force; };
  VECTOR3_FIELD_3D& vorticity()     { return _vorticity; };
  FIELD_3D&         density()       { return _density; };
  FIELD_3D&         pressure()      { return _pressure; };
  FIELD_3D&         divergence()    { return _divergence; };
  VECTOR3_FIELD_3D& preprojection() { return _preprojection; };
  VECTOR3_FIELD_3D& prediffusion()  { return _prediffusion; };
  VECTOR3_FIELD_3D& prevorticity()  { return _prevorticity; };
  VECTOR3_FIELD_3D& preadvection()  { return _preadvection; };
  VECTOR3_FIELD_3D& postadvection() { return _postadvection; };
  VECTOR3_FIELD_3D& postIOP()       { return _postIOP; };
  VECTOR3_FIELD_3D& postIOPAndPressure() { return _postIOPAndPressure; };
  VECTOR3_FIELD_3D& cachedQuadratic()  { return _cachedQuadratic; };
  VECTOR3_FIELD_3D& midMaccormack()  { return _midMaccormack; };
  FIELD_3D& cachedPressure()    { return _cachedPressure; };
  FIELD_3D& cachedDivergence()  { return _cachedDivergence; };
  vector<OBSTACLE*>& dirichletObstacles() { return _dirichletObstacles; };
  vector<OBSTACLE*>& neumannObstacles() { return _neumannObstacles; };
  vector<VECTOR3_FIELD_3D>& cachedVelocitiesIOP() { return _cachedVelocitiesIOP; };
  vector<FIELD_3D>& cachedPressuresIOP() { return _cachedPressuresIOP; };
  vector<FIELD_3D>& cachedDivergencesIOP() { return _cachedDivergencesIOP; };
  unsigned int& domainBcFront()  { return _domainBcFront; };
  unsigned int& domainBcBack()   { return _domainBcBack; };
  unsigned int& domainBcLeft()   { return _domainBcLeft; };
  unsigned int& domainBcRight()  { return _domainBcRight; };
  unsigned int& domainBcTop()    { return _domainBcTop; };
  unsigned int& domainBcBottom() { return _domainBcBottom; };
  SPARSE_MATRIX& neumannIOPMatrix() { return _neumannIOP; };
  string& snapshotPath() { return _snapshotPath; };
  vector<VEC3F>& glyphPoints() { return _glyphPoints; };
  vector<VEC3F>& glyphNormals() { return _glyphNormals; };

  void addBottomForce();

  void writeGz(string filename) const;
  void readGz(string filename);

  void appendStreams() const;
  void appendStreamsIOP() const;

  // create a procedural velocity field
  void generateProceduralVelocity();

  // set the neumannIOP obstacle matrix
  void setPeeledSparseMovingIOP(BOX* box);

  // set the Neumann velocity by appending a 1 in the last entry
  void setVelocityNeumann();

protected:  
  // dimensions
  int _xRes, _yRes, _zRes, _maxRes;
	int _xs, _ys, _zs;
  VEC3I _res;
  int _totalCells;
  int _slabSize;
  Real _dx;
  Real _totalTime;
  int _totalSteps;
  int _totalImgDumps;
  int _totalVelDumps;

  // fields
  FIELD_3D _density;
  FIELD_3D _densityOld;
  FIELD_3D _heat;
  FIELD_3D _heatOld;
  FIELD_3D _pressure;
  FIELD_3D _divergence;


  VECTOR3_FIELD_3D _velocity;
  VECTOR3_FIELD_3D _velocityOld;
  VectorXd _velocityNeumann;

  VECTOR3_FIELD_3D _dudt0;
  VECTOR3_FIELD_3D _dudt1;

  VECTOR3_FIELD_3D _force;
  VECTOR3_FIELD_3D _vorticity;
  unsigned char*  _obstacles;
  bool* _dirichletObstacleField;

  // store the pre-projection velocity in case the reduced solver needs them
  VECTOR3_FIELD_3D _preprojection;
  
  // store the pre-diffusion velocity in case the reduced solver needs them
  VECTOR3_FIELD_3D _prediffusion;
  
  // store the pre-advection velocity in case the reduced solver needs them
  VECTOR3_FIELD_3D _preadvection;
  
  // cache the result of the Treuille-style quadratic in case the reduced solver needs it
  VECTOR3_FIELD_3D _cachedQuadratic;
  
  // store the post-advection velocity in case the reduced solver needs them
  VECTOR3_FIELD_3D _postadvection;

  // cache the velocity field from before vorticity confinement has been added
  VECTOR3_FIELD_3D _prevorticity;

  // cache the intermediate Maccormack field
  VECTOR3_FIELD_3D _midMaccormack;
  
  // cache the result of Neumann IOP
  VECTOR3_FIELD_3D _postIOP;

  // cache the result of a boundary stomp + pressure projection
  VECTOR3_FIELD_3D _postIOPAndPressure;

  FIELD_3D _cachedPressure;
  FIELD_3D _cachedDivergence;

  vector<FIELD_3D> _cachedPressuresIOP;
  vector<FIELD_3D> _cachedDivergencesIOP;
  vector<VECTOR3_FIELD_3D> _cachedVelocitiesIOP;

  VEC3F _center;
  VEC3F _lengths;

  // Preconditioner fields
  FIELD_3D _z;
  FIELD_3D _precon;
 
  // Matrix fields	
  FIELD_3D _Adiag;
  FIELD_3D _Aplusi;
  FIELD_3D _Aplusj;
  FIELD_3D _Aplusk;
  FIELD_3D _Aminui;
  FIELD_3D _Aminuj;
  FIELD_3D _Aminuk;
	
  // CG fields
  FIELD_3D _residual;
  FIELD_3D _direction;
  FIELD_3D _q;
  int _iterations;

  // simulation constants
  Real _dt;
  Real _dtOld;
  Real _buoyancy;
  Real _wind;
  Real _vorticityEps;
  Real _heatDiffusion;
  Real _solverEps;

  // cache the damping matrix
  SPARSE_MATRIX _dampingMatrix;
  SPARSE_MATRIX _peeledDampingFull;

  // array of Dirichlet obstacles
  vector<OBSTACLE*> _dirichletObstacles;
  
  // array of Neumann obstacles
  vector<OBSTACLE*> _neumannObstacles;

  VECTOR3_FIELD_3D _neumannVelocityField;
  bool* _neumannObstacleField;

  // file paths
  string _snapshotPath;

  // domain boundary conditions
  unsigned int _domainBcFront;
  unsigned int _domainBcBack;
  unsigned int _domainBcBottom;
  unsigned int _domainBcTop;
  unsigned int _domainBcLeft;
  unsigned int _domainBcRight;

  // Neumann matrix for IOP
  SPARSE_MATRIX _neumannIOP;
  SPARSE_MATRIX _dirichletIOP;

  // ADJ: added _fullIOP and peeledIOP
  SPARSE_MATRIX _fullIOP;
  SPARSE_MATRIX _peeledIOP;

  // procedurally generated velocity basis
  VECTOR3_FIELD_3D _proceduralVelocity;

  vector<MATRIX3> _glyphFrames;

  vector<VEC3F> _glyphPoints;
  vector<VEC3F> _glyphNormals;

  // timestepping functions
  void addVorticity();
  void addBuoyancy(Real *field);
  void addWind(Real *field);

  // preconditioner stuff
  void zeroBoundary(FIELD_3D& x, unsigned char* skip);
	
  void buildMICPreconditioner();
  void applyMICPreconditioner(FIELD_3D& zf, FIELD_3D& r, FIELD_3D& q);	
	
  void buildA(unsigned char* skip, bool heat);	
  void multA(const FIELD_3D& xf, FIELD_3D& bf);
  // solver stuff
  void project();
  void projectIOP();
  void solvePoisson(FIELD_3D& field, FIELD_3D& b, unsigned char* skip, bool heat = false);
	
  // handle obstacle boundaries
  void setObstacleBoundaries();

  // build the damping matrix
  void buildPeeledDampingMatrixFull();

  // build the peeled IOP stomping matrix
  void buildPeeledSparseIOP(SPARSE_MATRIX& A, const VEC3I& center, double radius);

 
  // advect using first order semi-Lagrangian
  void advectStam();
  VEC3F cellCenter(int x, int y, int z);
};

#endif

