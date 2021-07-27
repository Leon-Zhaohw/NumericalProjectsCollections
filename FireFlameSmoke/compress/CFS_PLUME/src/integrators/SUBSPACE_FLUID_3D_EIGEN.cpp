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
// SUBSPACE_FLUID_3D_EIGEN.cpp: implementation of the SUBSPACE_FLUID_3D_EIGEN class.
//
//////////////////////////////////////////////////////////////////////

#include "SUBSPACE_FLUID_3D_EIGEN.h"
#include "BIG_MATRIX.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

SUBSPACE_FLUID_3D_EIGEN::SUBSPACE_FLUID_3D_EIGEN(int xRes, int yRes, int zRes, const string& reducedPath, unsigned int* boundaries, bool usingIOP, bool loadNothing) :
  FLUID_3D_MIC(xRes, yRes, zRes, 4, boundaries), _reducedPath(reducedPath), _discardThreshold(-1), _usingIOP(usingIOP)
{
  _snapshotPath = string("");

  // set up the boundary conditions
  if (boundaries != NULL)
  {
    _domainBcFront  = boundaries[0];
    _domainBcBack   = boundaries[1];
    _domainBcLeft   = boundaries[2];
    _domainBcRight  = boundaries[3];
    _domainBcTop    = boundaries[4];
    _domainBcBottom = boundaries[5];
  }
  else
  {
    _domainBcFront  = 1; 
    _domainBcBack   = 1;
    _domainBcLeft   = 1;
    _domainBcRight  = 1;
    _domainBcTop    = 0;
    _domainBcBottom = 0;
  }

  if (!loadNothing)
  {
    if (usingIOP) {
      initOutOfCoreIOP();
    }
    else {
      cout << " Not using IOP in fluid constructor." << endl;
      initOutOfCore();
    }


    readAdvectionCubature();
  }
  else
  {
    // init the peeled dimensions
    _xPeeled = _xRes - 2;
    _yPeeled = _yRes - 2;
    _zPeeled = _zRes - 2;
    _slabPeeled = _xPeeled * _yPeeled;
  }
}

//////////////////////////////////////////////////////////////////////
// initialize the peeled version where there is a separate basis
// for each stage
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_EIGEN::initOutOfCore()
{

  // init the peeled dimensions
  _xPeeled = _xRes - 2;
  _yPeeled = _yRes - 2;
  _zPeeled = _zRes - 2;
  _slabPeeled = _xPeeled * _yPeeled;

  bool pcaBuilt = fileExists(_reducedPath + string("U.final.matrix")) &&
                  fileExists(_reducedPath + string("U.preproject.matrix")) &&
                  fileExists(_reducedPath + string("U.preadvect.matrix")) &&
                  fileExists(_reducedPath + string("U.prediffuse.matrix")) &&
                  fileExists(_reducedPath + string("U.pressure.matrix"));

  // check if pre-built matrices exist
  bool filesBuilt = fileExists(_reducedPath + string("projected.A.matrix")) &&
                    fileExists(_reducedPath + string("projected.ptof.matrix")) &&
                    fileExists(_reducedPath + string("projected.vtod.matrix")) &&
                    fileExists(_reducedPath + string("damping.peeled.matrix")) &&
                    fileExists(_reducedPath + string("projected.ptov.matrix")) &&
                    fileExists(_reducedPath + string("inverseProduct.matrix"));

  cout << " filesBuilt: " << filesBuilt << endl;
  if (!filesBuilt)
    buildOutOfCoreMatrices();
  else
  {
    string filename;
    
    filename = _reducedPath + string("projected.ptof.matrix");
    EIGEN::read(filename, _preprojectToFinal);

    filename = _reducedPath + string("projected.vtod.matrix");
    EIGEN::read(filename, _reducedVelocityToDivergence);

    filename = _reducedPath + string("projected.ptov.matrix");
    EIGEN::read(filename, _reducedPressureToVelocity);

    filename = _reducedPath + string("projected.A.matrix");
    EIGEN::read(filename, _reducedA);

    filename = _reducedPath + string("damping.peeled.matrix");
    EIGEN::read(filename, _dampingMatrixReduced);
    
    filename = _reducedPath + string("inverseProduct.matrix");
    bool success = EIGEN::read(filename, _inverseProduct);

    if (!success)
    {
      // Needs everything prior to be built already
      MatrixXd inverse = _reducedA.inverse();
      _inverseProduct = _reducedPressureToVelocity * inverse * _reducedVelocityToDivergence;
      filename = _reducedPath + string("inverseProduct.matrix");
      EIGEN::write(filename, _inverseProduct);
      TIMER::printTimings();
    }
  }
}
void SUBSPACE_FLUID_3D_EIGEN::initOutOfCoreIOP()
{
  cout << " Using initOutOfCoreIOP inside fluid constructor. " << endl;
  // init the peeled dimensions
  _xPeeled = _xRes - 2;
  _yPeeled = _yRes - 2;
  _zPeeled = _zRes - 2;
  _slabPeeled = _xPeeled * _yPeeled;

  bool pcaBuilt = fileExists(_reducedPath + string("U.final.matrix")) &&
                  fileExists(_reducedPath + string("U.iop.matrix"))   &&
                  fileExists(_reducedPath + string("U.preproject.matrix")) &&
                  fileExists(_reducedPath + string("U.preadvect.matrix")) &&
                  fileExists(_reducedPath + string("U.prediffuse.matrix")) &&
                  fileExists(_reducedPath + string("U.pressure.matrix"));

  // check if pre-built matrices exist
  bool filesBuilt = fileExists(_reducedPath + string("projected.A.matrix")) &&
                    fileExists(_reducedPath + string("projected.ptof.matrix")) &&
                    fileExists(_reducedPath + string("projected.vtod.matrix")) &&
                    fileExists(_reducedPath + string("damping.peeled.matrix")) &&
                    fileExists(_reducedPath + string("projected.ptov.matrix")) &&
                    // ADJ: switching this name to U'U.iop.subspace.matrix
                    // fileExists(_reducedPath + string("U'U.iop.subspace.matrix")) &&
                    fileExists(_reducedPath + string("inverseProduct.matrix"));
  
  cout << " filesBuilt: " << filesBuilt << endl;
  if (!filesBuilt)
    buildOutOfCoreMatricesIOP();
  else
  {
    string filename;
    
    filename = _reducedPath + string("projected.ptof.matrix");
    // actually using preproject to preadvect rather than
    // preproject to final here is necessary if
    // the splitting is permuted
    EIGEN::read(filename, _preprojectToFinal);

    filename = _reducedPath + string("projected.vtod.matrix");
    EIGEN::read(filename, _reducedVelocityToDivergence);

    filename = _reducedPath + string("projected.ptov.matrix");
    EIGEN::read(filename, _reducedPressureToVelocity);

    filename = _reducedPath + string("projected.A.matrix");
    EIGEN::read(filename, _reducedA);

    filename = _reducedPath + string("damping.peeled.matrix");
    EIGEN::read(filename, _dampingMatrixReduced);
    
    filename = _reducedPath + string("inverseProduct.matrix");
    bool success = EIGEN::read(filename, _inverseProduct);

    if (!success)
    {
      // Needs everything prior to be built already
      MatrixXd inverse = _reducedA.inverse();
      _inverseProduct = _reducedPressureToVelocity * inverse * _reducedVelocityToDivergence;
      filename = _reducedPath + string("inverseProduct.matrix");
      EIGEN::write(filename, _inverseProduct);
      TIMER::printTimings();
    }
  }
}
SUBSPACE_FLUID_3D_EIGEN::~SUBSPACE_FLUID_3D_EIGEN()
{
}

//////////////////////////////////////////////////////////////////////
// The reduced solver, with peeled boundaries, 
// with cubature enabled
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_EIGEN::stepReorderedCubatureStam()
{
  TIMER functionTimer(__FUNCTION__);
  Real goalTime = 0.1;
  Real currentTime = 0;

  // compute the CFL condition
  _dt = goalTime;

  // wipe forces
  _force.clear();

  // wipe boundaries
  _velocity.setZeroBorder();

  // compute the forces
  addBuoyancy(_heat.data());
  _velocity.axpy(_dt, _force);

  _force.clear();
  addVorticity();
  _velocity.axpy(_dt, _force);

  VECTOR::printVertical = false;

  advectHeatAndDensityStam();
  
  

  TIMER projectionTimer("Velocity projection");
  _qDot = _velocity.peeledProject(_preadvectU);
  

  projectionTimer.stop();

  reducedAdvectStagedStamFast();

  
  TIMER diffusionProjectionTimer("Diffusion projection");

  reducedPeeledDiffusion();

  diffusionProjectionTimer.stop();

  reducedStagedProject();

  // do the full space unprojection
  TIMER unprojectionTimer("Velocity unprojection");
  _velocity.peeledUnproject(_U, _qDot);
  unprojectionTimer.stop();


  currentTime += _dt;

  cout << " Simulation step " << _totalSteps << " done. " << endl;

	_totalTime += goalTime;
	_totalSteps++;

  // diff the current sim results against ground truth
  diffGroundTruth();
}
//////////////////////////////////////////////////////////////////////
// The reduced solver, with peeled boundaries, 
// with cubature enabled
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_EIGEN::stepReorderedCubatureStamTest()
{
  TIMER functionTimer(__FUNCTION__);
  VECTOR::printVertical = false;

  Real goalTime = 0.1;
  Real currentTime = 0;

  // compute the CFL condition
  _dt = goalTime;

  // benchmarks to test against
  VECTOR3_FIELD_3D velocityTest;
  FIELD_3D densityTest;

  // wipe forces
  _force.clear();
 
  // wipe boundaries
  _velocity.setZeroBorder();

  // compute the forces
  addBuoyancy(_heat.data());
  _velocity.axpy(_dt, _force);
  _force.clear();

  addVorticity();
  _velocity.axpy(_dt, _force);

 
  // grab the velocity preadvection
  velocityTest = _velocity; 

  // do a full-coordinates advection of just heat and density
  // advectHeatAndDensityStam();

  // for the test, we will do a full-coordinates advection of 
  // velocity, heat and density

  advectStam();

  // _velocity is now post-advected

  // grab the post-advection density
  densityTest = _density;
 
  // next, do a reduced-coordinates advection of pre-advect velocity 
  TIMER projectionTimer("Velocity projection");

  _qDot = velocityTest.peeledProject(_preadvectU);
  // _qDot = _velocity.peeledProject(_preadvectU);
  
  projectionTimer.stop();
 
  // this advects _qDot
  reducedAdvectStagedStamFast();

  // bring _qDot back into the full space
  // QUESTION: what basis shoud I pass to Unproject?
  velocityTest.peeledUnproject(_prediffuseU, _qDot);
  cout << "Advection test: \n";
  // compare velocityTest and densityTest to _velocity and _density
  diffTruth(velocityTest, densityTest);

  TIMER diffusionProjectionTimer("Diffusion projection");

  reducedPeeledDiffusion();

  diffusionProjectionTimer.stop();

  reducedStagedProject();

  // do the full space unprojection
  TIMER unprojectionTimer("Velocity unprojection");
  _velocity.peeledUnproject(_U, _qDot);
  unprojectionTimer.stop();

  currentTime += _dt;

  cout << " Simulation step " << _totalSteps << " done. " << endl;

	_totalTime += goalTime;
	_totalSteps++;

  // diff the current sim results against ground truth
  diffGroundTruth();
}

//////////////////////////////////////////////////////////////////////
// The reduced solver, with peeled boundaries, with a moving obstacle, 
// with cubature enabled.
// **No reordering of the splitting!**
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_EIGEN::stepMovingObstacle(BOX* box)
{
  TIMER functionTimer(__FUNCTION__);
  VECTOR::printVertical = false;

  Real goalTime = 0.1;
  Real currentTime = 0;

  // compute the CFL condition
  _dt = goalTime;

  // wipe forces
  _force.clear();

  // wipe boundaries
  _velocity.setZeroBorder();
  _density.setZeroBorder();

  // compute the forces
  _velocity.axpy(_dt, _force);
  _force.clear();

  addVorticity();
  _velocity.axpy(_dt, _force);

  // a debugging velocity
  // VECTOR3_FIELD_3D velocityTrue = _velocity;
  // FIELD_3D densityTrue = _density;

  TIMER projectionTimer("Initial velocity projection");
  _qDot = _velocity.peeledProject(_preadvectU);
  projectionTimer.stop();
  
  puts("Finished velocity projection.");

  // grab the preadvected values of _velocity and _density
  // VECTOR3_FIELD_3D preadvectVelocity = _velocity;
  // FIELD_3D preadvectDensity = _density;
  // FIELD_3D preadvectHeat = _heat;

  // grab the preadvected values of the 'old' fields

  // VECTOR3_FIELD_3D oldVelocity = _velocityOld;
  // FIELD_3D oldDensity = _densityOld;
  // FIELD_3D oldHeat = _heatOld;

  ////////////////////////////////////////////////////////////////////
  // full advection---modifies _velocity and _density and _heat

  // advectStam();
  // velocityTrue = _velocity;
  // densityTrue = _density;
  ////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////
  // reduced advection
  // _velocity = preadvectVelocity;
  // _density = preadvectDensity;
  // _heat = preadvectHeat;

  // _velocityOld = oldVelocity;
  // _densityOld = oldDensity;
  // _heatOld = oldHeat;

  TIMER advectionTimer("Advection timing");

  // advect just heat and density
  advectHeatAndDensityStam();
 
  puts("Finished full-space heat and density advection.");

  // subspace advection for velocity.
  // reads from _preadvectU
  reducedAdvectStagedStamFast();
  advectionTimer.stop();

  puts("Finished reduced advection.");
  ////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////
  // comparing the reduced advection to full advection
  // _velocity.peeledUnproject(_prediffuseU, _qDot);
  // cout << "Advection projection. \n";
  // diffTruth(velocityTrue, densityTrue);
  ///////////////////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////////////////////
  // Full-space diffusion
  // if (_peeledDampingFull.rows() <= 0) {
  //   buildPeeledDampingMatrixFull();
  // }
  // VectorXd after = _peeledDampingFull * _velocity.peelBoundary().flattenedEigen();
  // _velocity.setWithPeeled(after);
  // velocityTrue = _velocity;
  ////////////////////////////////////////////////////////////////////
  
  // subspace diffusion 
  TIMER diffusionProjectionTimer("Diffusion projection");
  reducedPeeledDiffusion();
  diffusionProjectionTimer.stop();

  puts("Finished reduced diffusion.");

  ////////////////////////////////////////////////////////////////////
  // full space comparison of diffusion
  // _velocity.peeledUnproject(_preprojectU, _qDot);
  // cout << "Diffusion projection.\n";
  // diffTruth(velocityTrue, densityTrue);
  ////////////////////////////////////////////////////////////////////
 
  ////////////////////////////////////////////////////////////////////
  // full space boundary stomping

  // update the full IOP matrix (for debugging)
  
  // setPeeledSparseMovingIOP(box);
  // setVelocityNeumann();
  // VectorXd afterIOP = _neumannIOP * _velocityNeumann; 
  // _velocity.setWithPeeled(afterIOP); 
  // velocityTrue = _velocity;

  // reduced IOP
  TIMER iopTiming("IOP timing");
  reducedSetMovingBox(box);

  puts("Finished reducedSetMovingBox.");

  ////////////////////////////////////////////////////////////////////
  // obstacle stomping comparison

  // _velocity.peeledUnproject(_projectionIOP, _qDot);
  // cout << "Stomping boundaries test. \n";
  // diffTruth(velocityTrue, densityTrue);

  ////////////////////////////////////////////////////////////////////
  // this will modify _velocity using a full space projection
  // project();
  // velocityTrue = _velocity;
  
  // reduced pressure project
  reducedStagedProject();

  puts("Finished reduced pressure project.");

  // _velocity.peeledUnproject(_U, _qDot);
  // cout << "Pressure projection. \n";
  // diffTruth(velocityTrue, densityTrue);
  iopTiming.stop();
  ////////////////////////////////////////////////////////////////////
  
  // do the full space unprojection
  TIMER unprojectionTimer("Velocity unprojection");
  _velocity.peeledUnproject(_U, _qDot);
  unprojectionTimer.stop();

  puts("Finished subspace unprojection.");

  currentTime += _dt;

  cout << " Simulation step " << _totalSteps << " done. " << endl;

	_totalTime += goalTime;
	_totalSteps++;

  // diff the current sim results against ground truth
  diffGroundTruth();
}

//////////////////////////////////////////////////////////////////////
// do a full-rank advection of heat and density
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_EIGEN::advectHeatAndDensityStam()
{
  TIMER functionTimer(__FUNCTION__);

	if(_domainBcLeft == 0) 
    _velocity.copyBorderX();
	else 
    _velocity.setZeroX();

	if(_domainBcTop == 0) 
    _velocity.copyBorderY();
	else 
    _velocity.setZeroY();

	if(_domainBcFront == 0) 
    _velocity.copyBorderZ();
	else 
    _velocity.setZeroZ();

	_density.swapPointers(_densityOld);
	_heat.swapPointers(_heatOld);

	const Real dt0 = _dt / _dx;
  VECTOR3_FIELD_3D::advect(dt0, _velocity, _densityOld, _density);
  VECTOR3_FIELD_3D::advect(dt0, _velocity, _heatOld, _heat);
	
  _density.setZeroBorder();
	_heat.setZeroBorder();
}

//////////////////////////////////////////////////////////////////////
// perform reduced order diffusion with separate boundary slabs
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_EIGEN::reducedPeeledDiffusion() 
{
  TIMER functionTimer(__FUNCTION__);

  // apply the reduced damping to the middle
  _qDot = _dampingMatrixReduced * _qDot;
}

//////////////////////////////////////////////////////////////////////
// get the projection error of a vector with respect to a basis
//////////////////////////////////////////////////////////////////////
Real SUBSPACE_FLUID_3D_EIGEN::projectionError(const MatrixXd& basis, const VectorXd& v)
{
  VectorXd after = basis * (basis.transpose() * v);

  return (v - after).norm();
}

//////////////////////////////////////////////////////////////////////
// compute pressure-to-velocity matrix
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_EIGEN::computePressureToVelocity()
{
  int xPeeled = _xRes - 2;
  int yPeeled = _yRes - 2;
  int zPeeled = _zRes - 2;
  int slabPeeled = xPeeled * yPeeled;
  int peeledDims = xPeeled * yPeeled * zPeeled;
  int totalCells = 3 * peeledDims;
  _pressureToVelocity.resize(totalCells, peeledDims);

	Real invDx = 1.0f / _dx;
	Real halfInvDx = 0.5f / _dx;
	int index = _slabSize + _xRes + 1;
	for (int z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
		for (int y = 1; y < _yRes - 1; y++, index += 2)
			for (int x = 1; x < _xRes - 1; x++, index++)
      {
        int peeledIndex = (x - 1) + (y - 1) * xPeeled + (z - 1) * slabPeeled;
        int peeledIndex3 = 3 * peeledIndex;

        if (x == 1)
        {
          _pressureToVelocity(peeledIndex3 + 0, peeledIndex + 1) -=  invDx;
          _pressureToVelocity(peeledIndex3 + 0, peeledIndex)     -= -invDx;
        }
        else if (x == _xRes - 2)
        {
          _pressureToVelocity(peeledIndex3 + 0, peeledIndex)     -=  invDx;
          _pressureToVelocity(peeledIndex3 + 0, peeledIndex - 1) -= -invDx;
        }
        else
        {
          _pressureToVelocity(peeledIndex3 + 0, peeledIndex + 1) -=  halfInvDx;
          _pressureToVelocity(peeledIndex3 + 0, peeledIndex - 1) -= -halfInvDx;
        }

        if (y == 1)
        {
          _pressureToVelocity(peeledIndex3 + 1, peeledIndex + xPeeled) -=  invDx;
          _pressureToVelocity(peeledIndex3 + 1, peeledIndex)           -= -invDx;
        }
        else if (y == _yRes - 2)
        {
          _pressureToVelocity(peeledIndex3 + 1, peeledIndex)           -=  invDx;
          _pressureToVelocity(peeledIndex3 + 1, peeledIndex - xPeeled) -= -invDx;
        }
        else
        {
          _pressureToVelocity(peeledIndex3 + 1, peeledIndex + xPeeled) -=  halfInvDx;
          _pressureToVelocity(peeledIndex3 + 1, peeledIndex - xPeeled) -= -halfInvDx;
        }

        if (z == 1)
        {
          _pressureToVelocity(peeledIndex3 + 2, peeledIndex + slabPeeled) -=  invDx;
          _pressureToVelocity(peeledIndex3 + 2, peeledIndex)              -= -invDx;
        }
        else if (z == _zRes - 2)
        {
          _pressureToVelocity(peeledIndex3 + 2, peeledIndex)              -=  invDx;
          _pressureToVelocity(peeledIndex3 + 2, peeledIndex - slabPeeled) -= -invDx;
        }
        else
        {
          _pressureToVelocity(peeledIndex3 + 2, peeledIndex + slabPeeled) -=  halfInvDx;
          _pressureToVelocity(peeledIndex3 + 2, peeledIndex - slabPeeled) -= -halfInvDx;
        }
			}
}

//////////////////////////////////////////////////////////////////////
// do a staged reduced order pressure projection
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_EIGEN::reducedStagedProject()
{
  TIMER functionTimer(__FUNCTION__);
  cout << "_preprojectToFinal.cols: " << _preprojectToFinal.cols() << endl;
  cout << "_qDot.size: " << _qDot.size() << endl;
  _qDot = _preprojectToFinal * _qDot + _inverseProduct * _qDot;
}

//////////////////////////////////////////////////////////////////////
// update the Neumann moving boundary condition complement matrix
// and the Neumann homogeneous vector
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_EIGEN::setPeeledSparseMovingIOPComplement(BOX* box)
{
  TIMER functionTimer(__FUNCTION__);

  int totalCells = 3 * (_xRes - 2) * (_yRes - 2) * (_zRes - 2);
  // if it's the first time this function is called, resize _neumannIOPcomplement
  if (_neumannIOPcomplement.cols() <= 0) {
     _neumannIOPcomplement = SPARSE_MATRIX(totalCells, totalCells);
  }
  
  _neumannIOPcomplement.clearAndStompSparsity();
  // if it's the first time, resize and zero out _neumannVector
  if (_neumannVector.size() <= 0) {
    _neumannVector = VectorXd::Zero(totalCells);
  }
  _neumannVector.setZero();

  int index = 0;
  VEC3F r(0.0, 0.0, 0.0);
  VEC3F point(0.0, 0.0, 0.0);
  VEC3F* velocityPointer = box->get_velocity();

  for (int z = 1; z < _zRes - 1; z++) {
    for (int y = 1; y < _yRes - 1; y++) {
      for (int x = 1; x < _xRes - 1; x++) {
        for (int i = 0; i < 3; i++, index++) {
          point = this->cellCenter(x, y, z);
          // if point is in the interior of the obstacle
          if ( box->inside(point) ) {
            // complement matrix, so put 1's instead of 0's
            _neumannIOPcomplement(index, index) = 1.0;
            // update the radial and linear velocities
            box->update_r(point, &r);
            box->update_linearVelocity(r);
            if (i == 0) {
              _neumannVector[index] = (*velocityPointer)[0];
            }
            else if (i == 1) {
              _neumannVector[index] = (*velocityPointer)[1];
            }
            else {
              _neumannVector[index] = (*velocityPointer)[2];
            }
          }
        }
      }
    }
  }        
}

//////////////////////////////////////////////////////////////////////
// do a reduced IOP Neumann boundary for the box
/////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_EIGEN::reducedSetMovingBox(BOX* box)
{
  TIMER functionTimer(__FUNCTION__);
  
  /*
  // update _neumannIOPcomplement and _neumannVector
  setPeeledSparseMovingIOPComplement(box);

  cout << " left rows: " << _projectionIOP.rows() << endl;
  cout << " right rows: " << _preprojectU.rows() << endl;
  cout << " left cols: " << _projectionIOP.cols() << endl;
  cout << " right cols: " << _preprojectU.cols() << endl;
  
  printf(" precomputed U'U dims: (%lu, %lu)\n",
      _projectionIOP_T_preprojectU.rows(), _projectionIOP_T_preprojectU.cols());

  _reducedIOP = _projectionIOP_T_preprojectU - 
                _neumannIOPcomplement.project(_projectionIOP, _preprojectU);

  puts(" Computed the IOP projection inside reducedSetMovingBox.");
  // project the homgeneous coordinate
  cout << " _neumannVector size: " << _neumannVector.size() << endl;
  VectorXd n_hat = _projectionIOP.transpose() * _neumannVector; 
  puts(" Computed the homogeneous coordinate projection inside reducedSetMovingBox.");
  */
  
  char matrixBuffer[256];
  char vectorBuffer[256];
  sprintf(matrixBuffer, "%sreducedIOP.matrix.%i", _reducedPath.c_str(), _totalSteps);
  sprintf(vectorBuffer, "%sreducedNeumann.vector.%i", _reducedPath.c_str(), _totalSteps);
  EIGEN::read(matrixBuffer, _reducedIOP);
  EIGEN::read(vectorBuffer, _reducedNeumannVector);
  _qDot = _reducedIOP * _qDot + _reducedNeumannVector;
}

//////////////////////////////////////////////////////////////////////
// diff the current sim results against ground truth
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_EIGEN::diffGroundTruth()
{
  TIMER functionTimer(__FUNCTION__);
  if (_fullRankPath.size() == 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " No snapshot path has been set to compare to! " << endl;
    exit(0);
  }

  char buffer[256];
  string path = _fullRankPath;
  sprintf(buffer, "%sfluid.%04i.fluid3d", path.c_str(), _totalSteps - 1);
  string filename(buffer);
  FLUID_3D_MIC ground;
  ground.readGz(filename);

  cout << "=====================================================" << endl;
  cout << " Ground truth difference for frame " << _totalSteps - 1 << endl;
  cout << "=====================================================" << endl;
  VectorXd diff = _velocity.peelBoundary().flattenedEigen() - ground.velocity().peelBoundary().flattenedEigen();
  _velocityErrorAbs.push_back(diff.norm());
  _velocityErrorRelative.push_back(diff.norm() / _velocity.peelBoundary().flattenedEigen().norm());
  cout << " velocity abs error:      " << _velocityErrorAbs.back() << endl;
  cout << " velocity relative error: " << _velocityErrorRelative.back() << endl;

  diff = _density.peelBoundary().flattenedEigen() - ground.density().peelBoundary().flattenedEigen();
  _densityErrorAbs.push_back(diff.norm());
  _densityErrorRelative.push_back(diff.norm() / _density.peelBoundary().flattenedEigen().norm());
  cout << " density abs error:      " << _densityErrorAbs.back() << endl;
  cout << " density relative error: " << _densityErrorRelative.back() << endl;
}

//////////////////////////////////////////////////////////////////////
// diff the passed in sim results against the true results
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_EIGEN::diffTruth(const VECTOR3_FIELD_3D& testVelocity, 
    const FIELD_3D& testDensity)
{
  cout << "Invoking diffTruth: \n";
  VectorXd diff = testVelocity.peelBoundary().flattenedEigen() - _velocity.peelBoundary().flattenedEigen();
  _velocityErrorAbs.push_back(diff.norm());
  _velocityErrorRelative.push_back(diff.norm() / _velocity.peelBoundary().flattened().norm2());
  cout << " velocity abs error:      " << _velocityErrorAbs.back() << endl;
  cout << " velocity relative error: " << _velocityErrorRelative.back() << endl;

  // diff = testDensity.peelBoundary().flattenedEigen() - _density.peelBoundary().flattenedEigen();
  // _densityErrorAbs.push_back(diff.norm());
  // _densityErrorRelative.push_back(diff.norm() / _density.peelBoundary().flattened().norm2());
  // cout << " density abs error:      " << _densityErrorAbs.back() << endl;
  // cout << " density relative error: " << _densityErrorRelative.back() << endl;
}

//////////////////////////////////////////////////////////////////////
// get the sub-basis associated with a cell
//////////////////////////////////////////////////////////////////////
MatrixXd SUBSPACE_FLUID_3D_EIGEN::cellBasisPeeled(const MatrixXd& U, const int index)
{
  TIMER functionTimer(__FUNCTION__);
  // decompose into x,y,z
  const int decompose = index;
  const int z = decompose / _slabPeeled;
  const int y = (decompose % _slabPeeled) / _xPeeled;
  const int x = (decompose % _slabPeeled) % _xPeeled;

  assert(x >= 0);
  assert(x < _xRes - 2);
  assert(y >= 0);
  assert(y < _yRes - 2);
  assert(z >= 0);
  assert(z < _zRes - 2);
  assert(3 * index < U.rows());

  return EIGEN::getRows(3 * index, 3, U);
}

//////////////////////////////////////////////////////////////////////
// get the sub-basis associated with a cell
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_EIGEN::cellBasisPeeledMemory(int prediffuseRows, int prediffuseCols, FILE* prediffuseFile, double* cellData,
    const int index)
{
  TIMER functionTimer(__FUNCTION__);
  // decompose into x,y,z
  const int decompose = index;
  const int z = decompose / _slabPeeled;
  const int y = (decompose % _slabPeeled) / _xPeeled;
  const int x = (decompose % _slabPeeled) % _xPeeled;

  assert(x >= 0);
  assert(x < _xRes - 2);
  assert(y >= 0);
  assert(y < _yRes - 2);
  assert(z >= 0);
  assert(z < _zRes - 2);
  assert(3 * index < prediffuseRows); 

  EIGEN::getRowsMemory(3 * index, 3, prediffuseFile, cellData, prediffuseRows, prediffuseCols); 
  
}

//////////////////////////////////////////////////////////////////////
// advect a single cell
//////////////////////////////////////////////////////////////////////

// TODO: integrate decoder here!
VectorXd SUBSPACE_FLUID_3D_EIGEN::advectCellStamPeeled(const MatrixXd& U, const Real& dt, const VectorXd& qDot, const int index)
{
  TIMER functionTimer(__FUNCTION__);
  // peeled coordinates were passed in -- need to promote to full grid
  const int decompose = index;
  const int z = decompose / _slabPeeled + 1;
  const int y = (decompose % _slabPeeled) / _xPeeled + 1;
  const int x = (decompose % _slabPeeled) % _xPeeled + 1;

  // try using Eigen's block functionality
  const int index3 = 3 * index;
  const int totalColumns = U.cols();
  VectorXd v = U.block(index3, 0, 3, totalColumns) * qDot;

  // backtrace
  const VEC3F velocity(v[0], v[1], v[2]);
  Real xTrace = x - dt * velocity[0];
  Real yTrace = y - dt * velocity[1];
  Real zTrace = z - dt * velocity[2];

  // clamp backtrace to grid boundaries
  xTrace = (xTrace < 1.5) ? 1.5 : xTrace;
  xTrace = (xTrace > _xRes - 2.5) ? _xRes - 2.5 : xTrace;
  yTrace = (yTrace < 1.5) ? 1.5 : yTrace;
  yTrace = (yTrace > _yRes - 2.5) ? _yRes - 2.5 : yTrace;
  zTrace = (zTrace < 1.5) ? 1.5 : zTrace;
  zTrace = (zTrace > _zRes - 2.5) ? _zRes - 2.5 : zTrace;

  // locate neighbors to interpolate --
  // since we're in peeled coordinates, the lookup needs to be modified slightly
  const int x0 = (int)xTrace - 1;
  const int x1 = (x0 != _xPeeled - 1) ? x0 + 1 : x0;
  const int y0 = (int)yTrace - 1;
  const int y1 = (y0 != _yPeeled - 1) ? y0 + 1 : y0;
  const int z0 = (int)zTrace - 1;
  const int z1 = (z0 != _zPeeled - 1) ? z0 + 1 : z0;

  // get interpolation weights
  const Real s1 = (xTrace - 1) - x0;
  const Real s0 = 1.0f - s1;
  const Real t1 = (yTrace - 1) - y0;
  const Real t0 = 1.0f - t1;
  const Real u1 = (zTrace - 1) - z0;
  const Real u0 = 1.0f - u1;

  const int z0Scaled = z0 * _slabPeeled;
  const int z1Scaled = z1 * _slabPeeled;
  const int y0Scaled = y0 * _xPeeled;
  const int y1Scaled = y1 * _xPeeled;

  const int i000 = 3 * (x0 + y0Scaled + z0Scaled);
  const int i010 = 3 * (x0 + y1Scaled + z0Scaled);
  const int i100 = 3 * (x1 + y0Scaled + z0Scaled);
  const int i110 = 3 * (x1 + y1Scaled + z0Scaled);
  const int i001 = 3 * (x0 + y0Scaled + z1Scaled);
  const int i011 = 3 * (x0 + y1Scaled + z1Scaled);
  const int i101 = 3 * (x1 + y0Scaled + z1Scaled);
  const int i111 = 3 * (x1 + y1Scaled + z1Scaled);

  // NOTE: it spends most of its time (+50%) here,
  const VectorXd v000 = U.block(i000, 0, 3, totalColumns) * qDot;
  // for example,

  // const VectorXd v000 = GetSubmatrix(i000, i000+3) * qDot;
  const VectorXd v010 = U.block(i010, 0, 3, totalColumns) * qDot;
  const VectorXd v100 = U.block(i100, 0, 3, totalColumns) * qDot;
  const VectorXd v110 = U.block(i110, 0, 3, totalColumns) * qDot;
  const VectorXd v001 = U.block(i001, 0, 3, totalColumns) * qDot;
  const VectorXd v011 = U.block(i011, 0, 3, totalColumns) * qDot;
  const VectorXd v101 = U.block(i101, 0, 3, totalColumns) * qDot;
  const VectorXd v111 = U.block(i111, 0, 3, totalColumns) * qDot;

  const Real w000 = u0 * s0 * t0;
  const Real w010 = u0 * s0 * t1;
  const Real w100 = u0 * s1 * t0;
  const Real w110 = u0 * s1 * t1;
  const Real w001 = u1 * s0 * t0;
  const Real w011 = u1 * s0 * t1;
  const Real w101 = u1 * s1 * t0;
  const Real w111 = u1 * s1 * t1;
  
  // interpolate
  // (indices could be computed once)
  //
  // NOTE: it's deceptive to think this cuts down on
  // multiplies, since they will all occur on a VECTOR,
  // not just a scalar
  //
  //return u0 * (s0 * (t0 * v000 + t1 * v010) +
  //             s1 * (t0 * v100 + t1 * v110)) +
  //       u1 * (s0 * (t0 * v001 + t1 * v011) +
  //             s1 * (t0 * v101 + t1 * v111));
  return w000 * v000 + w010 * v010 + w100 * v100 + w110 * v110 +
         w001 * v001 + w011 * v011 + w101 * v101 + w111 * v111;
}
//////////////////////////////////
// advect a single cell (original)
//////////////////////////////////
/*
VectorXd SUBSPACE_FLUID_3D_EIGEN::advectCellStamPeeled(const MatrixXd& U, const Real& dt, const VectorXd& qDot, const int index)
{
  TIMER functionTimer(__FUNCTION__);
  // peeled coordinates were passed in -- need to promote to full grid
  const int decompose = index;
  const int z = decompose / _slabPeeled + 1;
  const int y = (decompose % _slabPeeled) / _xPeeled + 1;
  const int x = (decompose % _slabPeeled) % _xPeeled + 1;

  // try using Eigen's block functionality
  const int index3 = 3 * index;
  const int totalColumns = U.cols();
  VectorXd v = U.block(index3, 0, 3, totalColumns) * qDot;

  // backtrace
  const VEC3F velocity(v[0], v[1], v[2]);
  Real xTrace = x - dt * velocity[0];
  Real yTrace = y - dt * velocity[1];
  Real zTrace = z - dt * velocity[2];

  // clamp backtrace to grid boundaries
  xTrace = (xTrace < 1.5) ? 1.5 : xTrace;
  xTrace = (xTrace > _xRes - 2.5) ? _xRes - 2.5 : xTrace;
  yTrace = (yTrace < 1.5) ? 1.5 : yTrace;
  yTrace = (yTrace > _yRes - 2.5) ? _yRes - 2.5 : yTrace;
  zTrace = (zTrace < 1.5) ? 1.5 : zTrace;
  zTrace = (zTrace > _zRes - 2.5) ? _zRes - 2.5 : zTrace;

  // locate neighbors to interpolate --
  // since we're in peeled coordinates, the lookup needs to be modified slightly
  const int x0 = (int)xTrace - 1;
  const int x1 = (x0 != _xPeeled - 1) ? x0 + 1 : x0;
  const int y0 = (int)yTrace - 1;
  const int y1 = (y0 != _yPeeled - 1) ? y0 + 1 : y0;
  const int z0 = (int)zTrace - 1;
  const int z1 = (z0 != _zPeeled - 1) ? z0 + 1 : z0;

  // get interpolation weights
  const Real s1 = (xTrace - 1) - x0;
  const Real s0 = 1.0f - s1;
  const Real t1 = (yTrace - 1) - y0;
  const Real t0 = 1.0f - t1;
  const Real u1 = (zTrace - 1) - z0;
  const Real u0 = 1.0f - u1;

  const int z0Scaled = z0 * _slabPeeled;
  const int z1Scaled = z1 * _slabPeeled;
  const int y0Scaled = y0 * _xPeeled;
  const int y1Scaled = y1 * _xPeeled;

  const int i000 = 3 * (x0 + y0Scaled + z0Scaled);
  const int i010 = 3 * (x0 + y1Scaled + z0Scaled);
  const int i100 = 3 * (x1 + y0Scaled + z0Scaled);
  const int i110 = 3 * (x1 + y1Scaled + z0Scaled);
  const int i001 = 3 * (x0 + y0Scaled + z1Scaled);
  const int i011 = 3 * (x0 + y1Scaled + z1Scaled);
  const int i101 = 3 * (x1 + y0Scaled + z1Scaled);
  const int i111 = 3 * (x1 + y1Scaled + z1Scaled);

  // NOTE: it spends most of its time (+50%) here,
  const VectorXd v000 = U.block(i000, 0, 3, totalColumns) * qDot;
  const VectorXd v010 = U.block(i010, 0, 3, totalColumns) * qDot;
  const VectorXd v100 = U.block(i100, 0, 3, totalColumns) * qDot;
  const VectorXd v110 = U.block(i110, 0, 3, totalColumns) * qDot;
  const VectorXd v001 = U.block(i001, 0, 3, totalColumns) * qDot;
  const VectorXd v011 = U.block(i011, 0, 3, totalColumns) * qDot;
  const VectorXd v101 = U.block(i101, 0, 3, totalColumns) * qDot;
  const VectorXd v111 = U.block(i111, 0, 3, totalColumns) * qDot;

  const Real w000 = u0 * s0 * t0;
  const Real w010 = u0 * s0 * t1;
  const Real w100 = u0 * s1 * t0;
  const Real w110 = u0 * s1 * t1;
  const Real w001 = u1 * s0 * t0;
  const Real w011 = u1 * s0 * t1;
  const Real w101 = u1 * s1 * t0;
  const Real w111 = u1 * s1 * t1;
  
  // interpolate
  // (indices could be computed once)
  //
  // NOTE: it's deceptive to think this cuts down on
  // multiplies, since they will all occur on a VECTOR,
  // not just a scalar
  //
  //return u0 * (s0 * (t0 * v000 + t1 * v010) +
  //             s1 * (t0 * v100 + t1 * v110)) +
  //       u1 * (s0 * (t0 * v001 + t1 * v011) +
  //             s1 * (t0 * v101 + t1 * v111));
  return w000 * v000 + w010 * v010 + w100 * v100 + w110 * v110 +
         w001 * v001 + w011 * v011 + w101 * v101 + w111 * v111;
}
*/

//////////////////////////////////////////////////////////////////////
// advect a single cell
//////////////////////////////////////////////////////////////////////
VectorXd SUBSPACE_FLUID_3D_EIGEN::advectCellStamPeeled(const MatrixXd& U, const MatrixXd& cellU, const Real& dt, const VectorXd& qDot, const int index)
{
  TIMER functionTimer(__FUNCTION__);
  // peeled coordinates were passed in -- need to promote to full grid
  const int decompose = index;
  const int z = decompose / _slabPeeled + 1;
  const int y = (decompose % _slabPeeled) / _xPeeled + 1;
  const int x = (decompose % _slabPeeled) % _xPeeled + 1;

  // get the velocity to backtrace
  VectorXd v = cellU * qDot;

  // backtrace
  const VEC3F velocity(v[0], v[1], v[2]);
  Real xTrace = x - dt * velocity[0];
  Real yTrace = y - dt * velocity[1];
  Real zTrace = z - dt * velocity[2];

  // clamp backtrace to grid boundaries
  
  // keeping this comment block here for reference
  xTrace = (xTrace < 1.5) ? 1.5 : xTrace;
  xTrace = (xTrace > _xRes - 2.5) ? _xRes - 2.5 : xTrace;
  yTrace = (yTrace < 1.5) ? 1.5 : yTrace;
  yTrace = (yTrace > _yRes - 2.5) ? _yRes - 2.5 : yTrace;
  zTrace = (zTrace < 1.5) ? 1.5 : zTrace;
  zTrace = (zTrace > _zRes - 2.5) ? _zRes - 2.5 : zTrace;

  // locate neighbors to interpolate --
  // since we're in peeled coordinates, the lookup needs to be modified slightly
  
  // keeping this comment block here for reference
  const int x0 = (int)xTrace - 1;
  const int x1 = x0 + 1;
  const int y0 = (int)yTrace - 1;
  const int y1 = y0 + 1;
  const int z0 = (int)zTrace - 1;
  const int z1 = z0 + 1;

  // get interpolation weights
  const Real s1 = (xTrace - 1) - x0;
  const Real s0 = 1.0f - s1;
  const Real t1 = (yTrace - 1) - y0;
  const Real t0 = 1.0f - t1;
  const Real u1 = (zTrace - 1) - z0;
  const Real u0 = 1.0f - u1;

  const int z0Scaled = z0 * _slabPeeled;
  const int z1Scaled = z1 * _slabPeeled;
  const int y0Scaled = y0 * _xPeeled;
  const int y1Scaled = y1 * _xPeeled;

  const int i000 = 3 * (x0 + y0Scaled + z0Scaled);
  const int i010 = 3 * (x0 + y1Scaled + z0Scaled);
  const int i100 = 3 * (x1 + y0Scaled + z0Scaled);
  const int i110 = 3 * (x1 + y1Scaled + z0Scaled);
  const int i001 = 3 * (x0 + y0Scaled + z1Scaled);
  const int i011 = 3 * (x0 + y1Scaled + z1Scaled);
  const int i101 = 3 * (x1 + y0Scaled + z1Scaled);
  const int i111 = 3 * (x1 + y1Scaled + z1Scaled);

  // NOTE: it spends most of its time (+50%) here,
  // unprojecting
  const int totalColumns = U.cols();
  const VectorXd v000 = U.block(i000, 0, 3, totalColumns) * qDot;
  const VectorXd v010 = U.block(i010, 0, 3, totalColumns) * qDot;
  const VectorXd v100 = U.block(i100, 0, 3, totalColumns) * qDot;
  const VectorXd v110 = U.block(i110, 0, 3, totalColumns) * qDot;
  const VectorXd v001 = U.block(i001, 0, 3, totalColumns) * qDot;
  const VectorXd v011 = U.block(i011, 0, 3, totalColumns) * qDot;
  const VectorXd v101 = U.block(i101, 0, 3, totalColumns) * qDot;
  const VectorXd v111 = U.block(i111, 0, 3, totalColumns) * qDot;

  const Real w000 = u0 * s0 * t0;
  const Real w010 = u0 * s0 * t1;
  const Real w100 = u0 * s1 * t0;
  const Real w110 = u0 * s1 * t1;
  const Real w001 = u1 * s0 * t0;
  const Real w011 = u1 * s0 * t1;
  const Real w101 = u1 * s1 * t0;
  const Real w111 = u1 * s1 * t1;
  
  // interpolate
  // (indices could be computed once)
  //
  // NOTE: it's deceptive to think this cuts down on
  // multiplies, since they will all occur on a VECTOR,
  // not just a scalar
  //
  //return u0 * (s0 * (t0 * v000 + t1 * v010) +
  //             s1 * (t0 * v100 + t1 * v110)) +
  //       u1 * (s0 * (t0 * v001 + t1 * v011) +
  //             s1 * (t0 * v101 + t1 * v111));
  return w000 * v000 + w010 * v010 + w100 * v100 + w110 * v110 +
         w001 * v001 + w011 * v011 + w101 * v101 + w111 * v111;
}

//////////////////////////////////////////////////////////////////////
// read in a cubature scheme
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_EIGEN::readAdvectionCubature()
{
  string filename = _reducedPath + string("cubature");
  cout << " Trying to read in cubature file  " << filename.c_str() << " ... "; flush(cout);

  // read in the cubature file
  FILE* file;
  file = fopen(filename.c_str(), "rb");

  if (file == NULL)
  {
    cout << " No cubature file " << filename.c_str() << " found! " << endl;
    cout << " Maybe you're still training? " << endl;
    return;
  }

  // read dimensions
  int size;
  fread((void*)&size, sizeof(int), 1, file);

  // read out the key tet indices
  int cellIndex;
  for (int x = 0; x < size; x++)
  {
    fread((void*)&cellIndex, sizeof(int), 1, file);
    _keyAdvectionCells.push_back(cellIndex);
  }

  // read out the key tet weights
  double weight;
  for (int x = 0; x < size; x++)
  {
    fread((void*)&weight, sizeof(double), 1, file);
    _keyAdvectionWeights.push_back(weight);
  }

  VECTOR weightVector(_keyAdvectionWeights);
  VECTOR::printVertical = false;
  cout << " Read in " << _keyAdvectionCells.size() << " cubature points " << endl;
  fclose(file);

  cout << " done." << endl;

  // look for precached advection matrices
  filename = _reducedPath + string("advection.cubature.cache");
  file = fopen(filename.c_str(), "rb");

  int totalPoints = _keyAdvectionCells.size();

  // if there's no cache, create one
  if (file == NULL)
  {
    //fclose(file);
    cout << " No advection cubature cache found. Precaching all the advection cubature matrices. " << endl;

    // make the before cache
    _advectionCubatureBefore.clear();
    // TODO: insert the decoder for U.preadvect
    string preadvectFile = _reducedPath + string("U.preadvect.matrix");
    EIGEN::read(preadvectFile, _preadvectU);
    for (int x = 0; x < totalPoints; x++)
    {
      const int index = _keyAdvectionCells[x];
      _advectionCubatureBefore.push_back(cellBasisPeeled(_preadvectU, index));
    }
    _preadvectU.resize(0,0);

    // make the after cache
    _advectionCubatureAfter.clear();
    string prediffuseFile = _reducedPath + string("U.prediffuse.matrix");
    EIGEN::read(prediffuseFile, _prediffuseU);
    for (int x = 0; x < totalPoints; x++)
    {
      const int index = _keyAdvectionCells[x];
      const Real weight = _keyAdvectionWeights[x];
      
      const MatrixXd cellU = cellBasisPeeled(_prediffuseU, index);
      _advectionCubatureAfter.push_back(weight * cellU.transpose());
    }
    _prediffuseU.resize(0,0);

    // write the results out to a file
    filename = _reducedPath + string("advection.cubature.cache");
    file = fopen(filename.c_str(), "wb");

    if (file == NULL)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " Couldn't open advection cubature cache file " << filename.c_str() << "!!! " << endl;
      exit(0);
    }

    for (int x = 0; x < totalPoints; x++)
      EIGEN::write(file, _advectionCubatureBefore[x]);
    for (int x = 0; x < totalPoints; x++)
      EIGEN::write(file, _advectionCubatureAfter[x]);
    fclose(file);
  }
  // if there's already a cache, read it in
  else
  {
    cout << " Found advection cubature cache! Reading ... " << flush;
    _advectionCubatureBefore.resize(_keyAdvectionCells.size());
    _advectionCubatureAfter.resize(_keyAdvectionCells.size());

    for (int x = 0; x < totalPoints; x++)
      EIGEN::read(file, _advectionCubatureBefore[x]);
    for (int x = 0; x < totalPoints; x++)
      EIGEN::read(file, _advectionCubatureAfter[x]);
    fclose(file);

    cout << " done. " << endl;
  }
}

//////////////////////////////////////////////////////////////////////
// do Stam-style advection using cubautre
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_EIGEN::reducedAdvectStagedStamFast()
{
  TIMER functionTimer(__FUNCTION__);
  VectorXd final(_advectionCubatureAfter[0].rows());
  final.setZero();
	const Real dt0 = _dt / _dx;

  vector<VectorXd> finals(_keyAdvectionCells.size());
  int totalPoints = _keyAdvectionCells.size();

  assert(_advectionCubatureAfter.size() > 0);
  assert(_advectionCubatureBefore.size() > 0);

  // loop through the cubature points
#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < totalPoints; x++)
  {
    const int index = _keyAdvectionCells[x];
    finals[x] = _advectionCubatureAfter[x] * advectCellStamPeeled(_preadvectU, _advectionCubatureBefore[x], dt0, _qDot, index);
  }

  for (int x = 0; x < totalPoints; x++)
    final += finals[x];
  
  _qDot = final;
}

//////////////////////////////////////////////////////////////////////
// build the staged bases and the projected matrices
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_EIGEN::buildOutOfCoreMatrices()
{
  TIMER functionTimer(__FUNCTION__);

  cout << "=====================================================" << endl;
  cout << " Building huge projected matrices " << endl;
  cout << "=====================================================" << endl;

  // in preparation for writing results to files
  string filename;

  TIMER pressureTimer("Pressure projection");
  // read in the needed matrices
  filename = _reducedPath + string("U.pressure.matrix");
  EIGEN::read(filename, _pressureU);

  int totalCells = (_xRes - 2) * (_yRes - 2) * (_zRes - 2);
  SPARSE_MATRIX_ARRAY sparseA(totalCells, totalCells);
  buildFlatA(sparseA, _obstacles);
  cout << " Projecting the pressure matrix ... " << flush;
  _reducedA = sparseA.projectVerySparse(_pressureU, _pressureU);
  sparseA.clear();
  cout << " done." << endl;

  filename = _reducedPath + string("projected.A.matrix");
  EIGEN::write(filename, _reducedA);
  _Asparse.clear();
  purge();
  pressureTimer.stop();
  TIMER::printTimings();

  TIMER vtodTimer("Velocity to div projection");
  // read in the needed matrices
  filename = _reducedPath + string("U.preproject.matrix");
  EIGEN::read(filename, _preprojectU);

  // build reduced velocity to divergence
  // Need: _pressureU and _preprojectU
  computeVelocityToDivergence();
  cout << " Projecting velocity to divergence ... " << flush;
  _reducedVelocityToDivergence = _velocityToDivergence.project(_pressureU, _preprojectU);
  // for IOP, use this!!
  // //////////
  // _reducedVelocityToDivergence = _velocityToDivergence.project(_pressureU, _IOP basis);
  /////////////
  cout << " done." << endl;

  filename = _reducedPath + string("projected.vtod.matrix");
  EIGEN::write(filename, _reducedVelocityToDivergence);
  _velocityToDivergence.clear();
   
  // stomp matrices to make room in memory
  _pressureU.resize(0,0);
  purge();
  vtodTimer.stop();
  TIMER::printTimings();

  TIMER dampingTimer("Damping matrix projection");
  // read in the needed matrices
  filename = _reducedPath + string("U.prediffuse.matrix");
  EIGEN::read(filename, _prediffuseU);

  // NEW
  cout << " Projecting damping matrix ... " << flush;
  int totalCellsD = 3 * (_xRes - 2) * (_yRes - 2) * (_zRes - 2);
  SPARSE_MATRIX_ARRAY sparseD(totalCellsD, totalCellsD);
  buildPeeledDampingMatrixFlat(sparseD);
  _dampingMatrixReduced = sparseD.projectVerySparse(_preprojectU, _prediffuseU);
  filename = _reducedPath + string("damping.peeled.matrix");
  EIGEN::write(filename, _dampingMatrixReduced);
  sparseD.clear();
  cout << "done. " << endl;

  // stomp matrices to make room in memory
  _prediffuseU.resize(0,0);
  purge();
  dampingTimer.stop();
  TIMER::printTimings();

  TIMER preprojectTimer("Preprojection projection");
  // read in the needed matrices

  // TODO: do we need U.final in full?
  filename = _reducedPath + string("U.final.matrix");
  // EIGEN::readBig(filename, _U);
  EIGEN::read(filename, _U);

  // need: preproject and U
  EIGEN::transposeProduct(_U, _preprojectU, _preprojectToFinal);
  // for iop, need this!
  // ///////////
  // EIGEN::transposeProduct(_U, _iopU, _preprojectToFinal);
  // //////////
  filename = _reducedPath + string("projected.ptof.matrix");
  EIGEN::write(filename, _preprojectToFinal);

  // stomp matrices to make room in memory
  _preprojectU.resize(0,0);
  purge();
  preprojectTimer.stop();
  TIMER::printTimings();

  TIMER ptovTimer("Pressure to velocity projection");
  // read in the needed matrices
  filename = _reducedPath + string("U.pressure.matrix");
  EIGEN::read(filename, _pressureU);

  // build reduced pressure to velocity
  // Need: pressureU and U
  computePressureToVelocity();
  _reducedPressureToVelocity = _pressureToVelocity.project(_U, _pressureU);
  filename = _reducedPath + string("projected.ptov.matrix");
  EIGEN::write(filename, _reducedPressureToVelocity);

  // stomp pressure, just in case
  _pressureU.resize(0,0);
  purge();
  ptovTimer.stop();
  TIMER::printTimings();
  

  // Needs everything prior to be built already
  MatrixXd inverse = _reducedA.inverse();
  _inverseProduct = _reducedPressureToVelocity * inverse * _reducedVelocityToDivergence;
  filename = _reducedPath + string("inverseProduct.matrix");
  EIGEN::write(filename, _inverseProduct);
  TIMER::printTimings();

  // clear out all memory, just to be sure
  stompAllBases();
  purge();
  
  cout << " Done building matrices " << endl;
  TIMER::printTimings();
}

//////////////////////////////////////////////////////////////////////
// build the staged bases and the projected matrices
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_EIGEN::buildOutOfCoreMatricesIOP()
{
  TIMER functionTimer(__FUNCTION__);

  cout << "=====================================================" << endl;
  cout << " Building huge projected matrices " << endl;
  cout << "=====================================================" << endl;

  // in preparation for writing results to files
  string filename;

  TIMER pressureTimer("Pressure projection");
  // read in the needed matrices
  filename = _reducedPath + string("U.pressure.matrix");
  EIGEN::read(filename, _pressureU);

  int totalCells = (_xRes - 2) * (_yRes - 2) * (_zRes - 2);
  SPARSE_MATRIX_ARRAY sparseA(totalCells, totalCells);
  buildFlatA(sparseA, _obstacles);

  // (builds Poisson)
  cout << " Projecting the pressure matrix ... " << flush;
  _reducedA = sparseA.projectVerySparse(_pressureU, _pressureU);
  sparseA.clear();
  cout << " done." << endl;

  filename = _reducedPath + string("projected.A.matrix");
  EIGEN::write(filename, _reducedA);
  _Asparse.clear();
  purge();
  pressureTimer.stop();
  TIMER::printTimings();


  TIMER vtodTimer("Velocity to div projection");
  // read in the needed matrices
  filename = _reducedPath + string("U.iop.matrix");
  EIGEN::read(filename, _iopU);

  // build reduced velocity to divergence
  // Need: _pressureU and _iopU
  // we still have _pressureU from before
  computeVelocityToDivergence();
  cout << " Projecting velocity to divergence ... " << flush;
  // for IOP, use this!!
  _reducedVelocityToDivergence = _velocityToDivergence.project(_pressureU, _iopU);
  cout << " done." << endl;

  filename = _reducedPath + string("projected.vtod.matrix");
  EIGEN::write(filename, _reducedVelocityToDivergence);
  _velocityToDivergence.clear();
   
  // stomp matrices to make room in memory
  // keep iop around for the next one
  _pressureU.resize(0,0);
  purge();
  vtodTimer.stop();
  TIMER::printTimings();

  TIMER preprojectTimer("Preprojection projection");
  // read in the needed matrices

  filename = _reducedPath + string("U.final.matrix");
  EIGEN::read(filename, _U);

  EIGEN::transposeProduct(_U, _iopU, _preprojectToFinal);
  filename = _reducedPath + string("projected.ptof.matrix");
  EIGEN::write(filename, _preprojectToFinal);
  
  // stomp matrices to make room in memory
  _iopU.resize(0,0);
  _U.resize(0,0);
  purge();
  preprojectTimer.stop();
  TIMER::printTimings();

  TIMER dampingTimer("Damping matrix projection");
  // read in the needed matrices
  filename = _reducedPath + string("U.prediffuse.matrix");
  EIGEN::read(filename, _prediffuseU);

  filename = _reducedPath + string("U.preproject.matrix");
  EIGEN::read(filename, _preprojectU);

  cout << " Projecting damping matrix ... " << flush;
  
  int totalCellsD = 3 * (_xRes - 2) * (_yRes - 2) * (_zRes - 2);
  SPARSE_MATRIX_ARRAY sparseD(totalCellsD, totalCellsD);
  buildPeeledDampingMatrixFlat(sparseD);
  cout << " sparseD.rows: " << sparseD.rows() << endl;
  cout << " left -> _preprojectU.rows: " << _preprojectU.rows() << endl;
  cout << " right -> _prediffuseU.rows: " << _prediffuseU.rows() << endl;
  // QUESTION: is this correct for _dampingMatrixReduced?
  _dampingMatrixReduced = sparseD.projectVerySparse(_preprojectU, _prediffuseU);
  filename = _reducedPath + string("damping.peeled.matrix");
  EIGEN::write(filename, _dampingMatrixReduced);
  sparseD.clear();
  cout << "done. " << endl;

  // stomp matrices to make room in memory
  _prediffuseU.resize(0,0);
  _preprojectU.resize(0,0);
  purge();
  dampingTimer.stop();
  TIMER::printTimings();

  // build reduced pressure to velocity
  // Need: pressureU and Ufinal 
  
  TIMER ptovTimer("pressure to velocity");
  filename = _reducedPath + string("U.pressure.matrix");
  EIGEN::read(filename, _pressureU);

  filename = _reducedPath + string("U.final.matrix");
  EIGEN::read(filename, _U);

  computePressureToVelocity();
  _reducedPressureToVelocity = _pressureToVelocity.project(_U, _pressureU);
  filename = _reducedPath + string("projected.ptov.matrix");
  EIGEN::write(filename, _reducedPressureToVelocity);
  // stomp pressure and U, just in case
  _pressureU.resize(0,0);
  _U.resize(0,0);
  purge();
  ptovTimer.stop();
  TIMER::printTimings();
  
  TIMER iopTimer("IOP precomputation"); 
  // read in the projection matrices
  // filename = _reducedPath + string("U.iop.matrix");
  // EIGEN::read(filename, _projectionIOP); 
  // filename = _reducedPath + string("U.preproject.matrix");
  // EIGEN::read(filename, _preprojectU);

  // projection into the subspace

  /*
  cout << " _projectionIOP.rows: " << _projectionIOP.rows() << endl;
  cout << " _projectionIOP.cols: " << _projectionIOP.cols() << endl;
  cout << " _preprojectU.rows: " << _preprojectU.rows() << endl;
  cout << " _preprojectU.cols: " << _preprojectU.cols() << endl;
  _projectionIOP_T_preprojectU = _projectionIOP.transpose() * _preprojectU; 
  filename = _reducedPath + string("U'U.iop.subspace.matrix");
  EIGEN::write(filename, _projectionIOP_T_preprojectU);
  */

  // buildObstacleMatrices();

  // stomp IOP after writing it
  // ADJ: we can't really stomp these in the moving obstacle case, can we?
  // Need to precompute one for each time step
  // _projectionIOP.resize(0, 0);
  // _preprojectU.resize(0, 0);
  // purge();

  iopTimer.stop();
  TIMER::printTimings();

  // Needs everything prior to be built already
  MatrixXd inverse = _reducedA.inverse();
  _inverseProduct = _reducedPressureToVelocity * inverse * _reducedVelocityToDivergence;
  filename = _reducedPath + string("inverseProduct.matrix");
  EIGEN::write(filename, _inverseProduct);
  TIMER::printTimings();

  // clear out all memory, just to be sure
  stompAllBases();
  purge();
  
  cout << " Done building matrices " << endl;
  TIMER::printTimings();
}


//////////////////////////////////////////////////////////////////////////////
// Write out an obstacle matrix for each time step
//////////////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_EIGEN::buildObstacleMatrices()
{
  bool filesBuilt = fileExists(_reducedPath + string("reducedIOP.matrix.0")) &&
                    fileExists(_reducedPath + string("reducedNeumann.vector.0"));

  if (filesBuilt) { 
    puts("Already built first frame of reduced matrices and vectors!");
    puts("Blithely assuming you also have built the rest of them and continuing on.");
    return; 
  }

  // read in the projection matrices
  string filename = _reducedPath + string("U.iop.matrix");
  EIGEN::read(filename, _projectionIOP); 
  filename = _reducedPath + string("U.preproject.matrix");
  EIGEN::read(filename, _preprojectU);

  puts(" Building obstacle matrices...");
  // build U'U matrix first
  _projectionIOP_T_preprojectU = _projectionIOP.transpose() * _preprojectU;

  // allocate buffers for writing to files
  char matrixBuffer[256];
  char vectorBuffer[256];

  printf(" Total reduced steps: %i\n", _totalReducedSteps);

  // loop through each step and write a unique matrix and vector per frame
  for (int step = 0; step < _totalReducedSteps; step++) {
    setPeeledSparseMovingIOPComplement(&_box); 
    _reducedIOP = _projectionIOP_T_preprojectU - 
                  _neumannIOPcomplement.project(_projectionIOP, _preprojectU);

    // generate step-based filename and write 
    sprintf(matrixBuffer, "%sreducedIOP.matrix.%i", _reducedPath.c_str(), step);
    EIGEN::write(matrixBuffer, _reducedIOP);

    // project the homgeneous coordinate
    _reducedNeumannVector = _projectionIOP.transpose() * _neumannVector;

    // generate step-based filename and write
    sprintf(vectorBuffer, "%sreducedNeumann.vector.%i", _reducedPath.c_str(), step);
    EIGEN::write(vectorBuffer, _reducedNeumannVector);
  
    // update the box position before the next step.
    // note that this is the internal box member variable,
    // so it will not mutate the declared box inside reducedMoving.cpp
    _box.translate_center();
    _box.spin();
    _box.update_time();
  }
  // stomp the huge matrices after precaching the subspace matrices is complete
  _projectionIOP.resize(0, 0);
  _preprojectU.resize(0, 0);
  purge();
}


//////////////////////////////////////////////////////////////////////////////
// Write out an obstacle matrix for each time step starting from 'start'
//////////////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_EIGEN::buildRemainingObstacleMatrices(int start)
{

  // read in the projection matrices
  string filename = _reducedPath + string("U.iop.matrix");
  EIGEN::read(filename, _projectionIOP); 
  filename = _reducedPath + string("U.preproject.matrix");
  EIGEN::read(filename, _preprojectU);

  puts(" Building obstacle matrices...");
  // build U'U matrix first
  _projectionIOP_T_preprojectU = _projectionIOP.transpose() * _preprojectU;

  // allocate buffers for writing to files
  char matrixBuffer[256];
  char vectorBuffer[256];

  printf(" Total reduced steps: %i\n", _totalReducedSteps);
  printf(" Starting at step: %i\n", start);

  for (int i = 0; i < start; i++) {
     // update the box position before the next step.
    // note that this is the internal box member variable,
    // so it will not mutate the declared box inside reducedMoving.cpp
    printf(" Updating box position at step: %i\n", i);
    _box.translate_center();
    _box.spin();
    _box.update_time();
  }

  // loop through each step and write a unique matrix and vector per frame
  for (int step = start; step < _totalReducedSteps; step++) {
    printf(" Computing matrices for step: %i...\n", step);
    setPeeledSparseMovingIOPComplement(&_box); 
    _reducedIOP = _projectionIOP_T_preprojectU - 
                  _neumannIOPcomplement.project(_projectionIOP, _preprojectU);

    // generate step-based filename and write 
    sprintf(matrixBuffer, "%sreducedIOP.matrix.%i", _reducedPath.c_str(), step);
    EIGEN::write(matrixBuffer, _reducedIOP);

    // project the homgeneous coordinate
    _reducedNeumannVector = _projectionIOP.transpose() * _neumannVector;

    // generate step-based filename and write
    sprintf(vectorBuffer, "%sreducedNeumann.vector.%i", _reducedPath.c_str(), step);
    EIGEN::write(vectorBuffer, _reducedNeumannVector);
  
    // update the box position before the next step.
    // note that this is the internal box member variable,
    // so it will not mutate the declared box inside reducedMoving.cpp
    _box.translate_center();
    _box.spin();
    _box.update_time();
  }
  // stomp the huge matrices after precaching the subspace matrices is complete
  _projectionIOP.resize(0, 0);
  _preprojectU.resize(0, 0);
  purge();
}

//////////////////////////////////////////////////////////////////////
// Write out the dims of the collection of subspace vectors into a matrix
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_EIGEN::writeSubspaceErrorMatrixDims(int simulationSnapshots)
{
  FILE* matrixFile;
  string filename = _reducedPath + string("qDot.matrix");

  matrixFile = fopen(filename.c_str(), "wb");
  if (matrixFile==NULL) {
    perror("Error opening matrixFile");
    exit(EXIT_FAILURE);
  }

  int rows = simulationSnapshots;
  int cols = simulationSnapshots;
  // DEBUG
  printf("Writing subspace error matrix dims: (%i, %i)\n", rows, cols);

  fwrite((void*)(&rows), sizeof(int), 1, matrixFile);
  fwrite((void*)(&cols), sizeof(int), 1, matrixFile);
  fclose(matrixFile);

}
//////////////////////////////////////////////////////////////////////
// Write the current subspace vector to the subspace matrix file
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_EIGEN::appendSubspaceVectors()
{
  FILE* matrixFile;
  string filename = _reducedPath + string("qDot.matrix");

  matrixFile = fopen(filename.c_str(), "ab");
  if (matrixFile==NULL) {
    perror("Error opening matrixFile");
    exit(EXIT_FAILURE);
  }

  auto count = _qDot.size();
  // DEBUG
  printf("qDot has size: %lu\n", count);
  fwrite((void*)(_qDot.data()), sizeof(double), count, matrixFile);

  fclose(matrixFile);

}


//////////////////////////////////////////////////////////////////////
// check of a file exists
//////////////////////////////////////////////////////////////////////
bool SUBSPACE_FLUID_3D_EIGEN::fileExists(const string& filename)
{
  FILE* file;
  file = fopen(filename.c_str(), "rb");
  
  if (file == NULL)
    return false;

  fclose(file);
  return true;
}

//////////////////////////////////////////////////////////////////////
// stomp all loaded bases
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_EIGEN::stompAllBases()
{
  _pressureU.resize(0,0);
  _preprojectU.resize(0,0);
  _prediffuseU.resize(0,0);
  _preadvectU.resize(0,0);
  _U.resize(0,0);
}

void SUBSPACE_FLUID_3D_EIGEN::loadPreadvect()
{
  string filename = _reducedPath + string("U.preadvect.matrix");
  EIGEN::read(filename, _preadvectU);
}

void SUBSPACE_FLUID_3D_EIGEN::wipePreadvect()
{
  _preadvectU.resize(0,0);
  purge();
}

void SUBSPACE_FLUID_3D_EIGEN::loadPrediffuse()
{
  string filename = _reducedPath + string("U.prediffuse.matrix");
  EIGEN::read(filename, _prediffuseU);
}

void SUBSPACE_FLUID_3D_EIGEN::wipePrediffuse()
{
  _prediffuseU.resize(0,0);
  purge();
}

//////////////////////////////////////////////////////////////////////
// stomp the other matrices and load the ones needed for cubature 
// training
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_EIGEN::loadCubatureTrainingBases()
{
  string filename;
  filename = _reducedPath + string("U.preadvect.matrix");
  EIGEN::read(filename, _preadvectU);

  filename = _reducedPath + string("U.prediffuse.matrix");
  EIGEN::read(filename, _prediffuseU);
}

//////////////////////////////////////////////////////////////////////
// load the bases needed for cubature runtime
//
// the path variable is there in case we want to load off the SSD
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_EIGEN::loadReducedRuntimeBases(string path)
{
  TIMER functionTimer(__FUNCTION__);
  if (path.length() == 0)
    path = _reducedPath;

  string filename;
  
  filename = path + string("U.preadvect.matrix");
  EIGEN::read(filename, _preadvectU);
 
  if (_preadvectU.rows() > 1000000)
    purge();

  filename = path + string("U.final.matrix");
  EIGEN::read(filename, _U);
  
  TIMER::printTimings();
  if (_U.rows() > 1000000)
    purge();
}
//////////////////////////////////////////////////////////////////////
// load ALL the bases needed for cubature runtime debugging
//
// the path variable is there in case we want to load off the SSD
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_EIGEN::loadReducedRuntimeBasesAll(string path)
{
  TIMER functionTimer(__FUNCTION__);
  if (path.length() == 0)
    path = _reducedPath;

  string filename;
  
  filename = path + string("U.preadvect.matrix");
  EIGEN::read(filename, _preadvectU);
  if (_preadvectU.rows() > 1000000)
    purge();
 
  filename = path + string("U.prediffuse.matrix");
  EIGEN::read(filename, _prediffuseU);
  if (_prediffuseU.rows() > 1000000)
    purge();

  filename = path + string("U.preproject.matrix");
  EIGEN::read(filename, _preprojectU);
  if (_preprojectU.rows() > 1000000)
    purge();

  filename = path + string("U.pressure.matrix");
  EIGEN::read(filename, _pressureU);
  if (_pressureU.rows() > 1000000)
    purge();

  filename = path + string("U.final.matrix");
  EIGEN::read(filename, _U);
  if (_U.rows() > 1000000)
    purge();
 
  filename = path + string("U.iop.matrix");
  EIGEN::read(filename, _projectionIOP);
  
  TIMER::printTimings();
  if (_projectionIOP.rows() > 1000000)
    purge(); 

  // this one is in the subspace, so it is relatively small
  /*
  filename = path + string("U'U.iop.subspace.matrix");
  EIGEN::read(filename, _projectionIOP_T_preprojectU);
  */

  TIMER::printTimings();
}

//////////////////////////////////////////////////////////////////////
// compute the velocity-to-divergence matrix
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_EIGEN::computeVelocityToDivergence()
{
  int xPeeled = _xRes - 2;
  int yPeeled = _yRes - 2;
  int zPeeled = _zRes - 2;
  int slabPeeled = xPeeled * yPeeled;
  int peeledDims = xPeeled * yPeeled * zPeeled;
  int totalCells = 3 * peeledDims;

  Real coeff = -_dx * 0.5;
  _velocityToDivergence.resize(peeledDims, totalCells);
  _velocityToDivergence.clearAndStompSparsity();

  for (int z = 0; z < zPeeled; z++)
    for (int y = 0; y < yPeeled; y++)
      for (int x = 0; x < xPeeled; x++)
      {
        int index = x + y * xPeeled + z * slabPeeled;
        int right = index + 1;
        int left  = index - 1;
        int up   = index + xPeeled;
        int down = index - xPeeled;
        int far  = index + slabPeeled;
        int near = index - slabPeeled;

        // this needs to be redone if Neumann is being used
        //
        // the zero will get constant folded away -- it's just there
        // as a reminder.
        if (x != xPeeled - 1)
          _velocityToDivergence(index, 3 * right + 0) += coeff;
        // ADJ: switching this from _domainBcRight, since it is the case when x == xPeeled - 1
        else if (_domainBcRight != 0)
          _velocityToDivergence(index, 3 * index + 0) += -coeff;
        else if (_domainBcRight == 0)
          _velocityToDivergence(index, 3 * left + 0) += coeff;
        else
        {
          cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
          cout << " FELL THROUGH! " << endl;
          exit(0);
        }

        if (x != 0)
          _velocityToDivergence(index, 3 * left + 0)  += -coeff;
        else if (_domainBcLeft != 0)
          _velocityToDivergence(index, 3 * index + 0) += coeff;
        else if (_domainBcLeft == 0)
          _velocityToDivergence(index, 3 * right + 0) += -coeff;
        else
        {
          cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
          cout << " FELL THROUGH! " << endl;
          exit(0);
        }

        if (y != yPeeled - 1)
          _velocityToDivergence(index, 3 * up + 1)   += coeff;
        // ADJ: added the case if _domainBcTop != 0 so as not to fall through
        else if (_domainBcTop != 0)
          _velocityToDivergence(index, 3 * index + 1) += -coeff;
        else if (_domainBcTop == 0)
          _velocityToDivergence(index, 3 * down + 1) += coeff;
        else
        {
          cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
          cout << " FELL THROUGH! " << endl;
          exit(0);
        }

        if (y != 0)
          _velocityToDivergence(index, 3 * down + 1) += -coeff;
        // ADJ: added the case if _domainBc != 0 here as well
        // also switched from top to bottom, since this will be when y == 0
        else if (_domainBcBottom != 0)
          _velocityToDivergence(index, 3 * index + 1) += coeff;
        else if (_domainBcBottom == 0)
          _velocityToDivergence(index, 3 * up + 1)   += -coeff;
        else
        {
          cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
          cout << " FELL THROUGH! " << endl;
          exit(0);
        }
        
        if (z != zPeeled - 1)
          _velocityToDivergence(index, 3 * far + 2)  += coeff;
        else if (_domainBcBack != 0)
          _velocityToDivergence(index, 3 * index + 2) += -coeff;
        // ADJ: added the == case and switched from front to back
        else if (_domainBcBack == 0)
          _velocityToDivergence(index, 3 * near + 2) += coeff;
        else
        {
          cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
          cout << " FELL THROUGH! " << endl;
          exit(0);
        }

        if (z != 0)
          _velocityToDivergence(index, 3 * near + 2) += -coeff;
        else if (_domainBcFront != 0)
          _velocityToDivergence(index, 3 * index + 2)  += coeff;
        // ADJ: added the == case
        else if (_domainBcFront == 0)
          _velocityToDivergence(index, 3 * far + 2) += -coeff;
        else
        {
          cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
          cout << " FELL THROUGH! " << endl;
          exit(0);
        }
      }
}

//////////////////////////////////////////////////////////////////////
// build the peeled damping matrix
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_EIGEN::buildPeeledDampingMatrixFlat(SPARSE_MATRIX_ARRAY& peeledDampingMatrix)
{
  TIMER functionTimer(__FUNCTION__);
  cout << " Building flat peeled damping matrix ... ";flush(cout);
	const Real w = 0.9;
  int xPeeled = _xRes - 2;
  int yPeeled = _yRes - 2;
  int zPeeled = _zRes - 2;
  int slabPeeled = xPeeled * yPeeled;
  const Real wMinus = 1.0 - w;
  const Real sixth = (1./ 6.)* w;

  // initially everybody is identity
  for (int x = 0; x < peeledDampingMatrix.rows(); x++)
    peeledDampingMatrix(x,x) = 1.0;

  // populate the rest
  for (int z = 0; z < zPeeled; z++) 
    for (int y = 0; y < yPeeled; y++)
      for (int x = 0; x < xPeeled; x++) 
      {
        int index = 3 * (x + y * xPeeled + z * slabPeeled);
        int left  = 3 * ((x + 1) + y * xPeeled + z * slabPeeled);
        int right = 3 * ((x - 1) + y * xPeeled + z * slabPeeled);
        int up   = 3 * (x + (y + 1) * xPeeled + z * slabPeeled);
        int down = 3 * (x + (y - 1) * xPeeled + z * slabPeeled);
        int near = 3 * (x + y * xPeeled + (z + 1) * slabPeeled);
        int far  = 3 * (x + y * xPeeled + (z - 1) * slabPeeled);

        for (int i = 0; i < 3; i++)
        {
          peeledDampingMatrix(index + i, index + i) = wMinus;

          if (x != xPeeled - 1)
            peeledDampingMatrix(index + i, left + i) = sixth;

          if (x != 0)
            peeledDampingMatrix(index + i, right + i) = sixth;

          if (y != yPeeled - 1)
            peeledDampingMatrix(index + i, up + i)    = sixth;

          if (y != 0)
            peeledDampingMatrix(index + i, down + i)  = sixth;

          if (z != zPeeled - 1)
            peeledDampingMatrix(index + i, near + i)  = sixth;

          if (z != 0)
            peeledDampingMatrix(index + i, far + i)   = sixth;
        }
      }

  cout << " done. " << endl;
}

//////////////////////////////////////////////////////////////////////
// build a sparse version of the Poisson matrix
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_EIGEN::buildFlatA(SPARSE_MATRIX_ARRAY& sparseA, unsigned char* skip)
{
  int xPeeled = _xRes - 2;
  int yPeeled = _yRes - 2;
  //int zPeeled = _zRes - 2;
  int slabPeeled = xPeeled * yPeeled;
  //int peeledDims = xPeeled * yPeeled * zPeeled;

	int x, y, z, index;
	Real Aoff = 1.0;
	index = _slabSize + _xRes + 1;
	for (z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
		for (y = 1; y < _yRes - 1; y++, index += 2)
			for (x = 1; x < _xRes - 1; x++, index++)
      {
				// if the cell is a variable
				if (!skip[index])
				{
          int peeledIndex = (x - 1) + (y - 1) * xPeeled + (z - 1) * slabPeeled;

					// set the matrix to the Poisson stencil in order
					if (!skip[index + _xs]) sparseA(peeledIndex, peeledIndex) += Aoff;
					if (!skip[index - _xs]) sparseA(peeledIndex, peeledIndex) += Aoff;
					if (!skip[index + _ys]) sparseA(peeledIndex, peeledIndex) += Aoff;
					if (!skip[index - _ys]) sparseA(peeledIndex, peeledIndex) += Aoff;
					if (!skip[index + _zs]) sparseA(peeledIndex, peeledIndex) += Aoff;
					if (!skip[index - _zs]) sparseA(peeledIndex, peeledIndex) += Aoff;

          if (!skip[index - 1] && x != 1)
            sparseA(peeledIndex, peeledIndex - 1) = -Aoff;
          if (!skip[index + 1] && x != _xRes - 2)
            sparseA(peeledIndex, peeledIndex + 1) = -Aoff;
          if (!skip[index - _xRes] && y != 1)
            sparseA(peeledIndex, peeledIndex - xPeeled) = -Aoff;
          if (!skip[index + _xRes] && y != _yRes - 2)
            sparseA(peeledIndex, peeledIndex + xPeeled) = -Aoff;
          if (!skip[index - _slabSize] && z != 1)
            sparseA(peeledIndex, peeledIndex - slabPeeled) = -Aoff;
          if (!skip[index + _slabSize] && z != _zRes - 2)
            sparseA(peeledIndex, peeledIndex + slabPeeled) = -Aoff;
				}
			}
}

// compute the spatial coordinates from integer indexing
VEC3F SUBSPACE_FLUID_3D_EIGEN::cellCenter(int x, int y, int z)
{
  TIMER functionTimer(__FUNCTION__);
  // assuming lengths are all 1.0
  VEC3F halfLengths(0.5, 0.5, 0.5); 

  // set it to the lower corner
  // **********************************************************************
  // ADJ: commenting this out for now. might break the static obstacle code
  // probably the variable _center was not (0.5, 0.5, 0.5) as expected?
  // VEC3F final = _center - halfLengths;
  // **********************************************************************

  VEC3F final(0.0, 0.0, 0.0);

  double dx = 1.0 / _xRes;
  double dy = 1.0 / _yRes;
  double dz = 1.0 / _zRes;
  // displace to the NNN corner
  final[0] += x * dx;
  final[1] += y * dy;
  final[2] += z * dz;

  // displace it to the cell center
  final[0] += dx * 0.5;
  final[1] += dy * 0.5;
  final[2] += dz * 0.5;

  return final;
}

