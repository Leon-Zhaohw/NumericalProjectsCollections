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
// SUBSPACE_TET_MESH.h: interface for the SUBSPACE_TET_MESH class.
//
//////////////////////////////////////////////////////////////////////

#include <string>
#include "SUBSPACE_TET_MESH.h"
#include <INVERTIBLE.h>
#ifdef USING_OPENMP
#include <omp.h>
#endif
#if USING_OSX
#include <Accelerate/Accelerate.h>
#endif

#include <TIMER.h>
#include <float.h>

#include <memory.h>

//////////////////////////////////////////////////////////////////////
// Constructor for tetrahedra
//////////////////////////////////////////////////////////////////////
SUBSPACE_TET_MESH::SUBSPACE_TET_MESH(const char* filename,
                                     MATERIAL** materials,
                                     int totalMaterials,
                                     bool simulateFullspace,
                                     const char* eigenvectors,
                                     const char* eigenvalues,
                                     Real scale) :
  TET_MESH(filename, materials, totalMaterials, simulateFullspace, scale),
  _filename(filename),
  _keyTets(NULL), 
  _keyWeights(0),
  //_q(_x),
  _totalKeyTets(0),
  _stiffnessDensity(NULL),
  _pFpUrepacked( NULL ),
  _tetSubBases( NULL ),
  _tetSubBasesTransposed( NULL ),
  _tetBmatrices( NULL )
{
  // try to find precomputed eigensystem
  if (eigenvectors != NULL && eigenvalues != NULL)
  {
    _eigenvalues.read(eigenvalues);
    _UBasis.read(eigenvectors);
    //cout << " Precomputed eigenvalues found: " << _eigenvalues;
    cout << " eigenvalues size: " << _eigenvalues.size() << endl;
    cout << " eigenvectors size: " << _UBasis.rows() << " " << _UBasis.cols() << endl;
  }
  else
  {
    // try probing the default filenames
    string valuename = string(filename) + string(".eigenvalues.vector");
    string vectorname = string(filename) + string(".eigenvectors.matrix");
    FILE* valuefile = fopen(valuename.c_str(), "rb");
    if (valuefile != NULL)
    {
      _eigenvalues.read(valuename.c_str());
      cout << " Precomputed eigenvalues found: " << endl << _eigenvalues;
    }
    
    FILE* vectorfile = fopen(vectorname.c_str(), "rb");
    if (vectorfile != NULL)
    {
      _UBasis.read(vectorname.c_str());
      cout << " Precomputed eigenvectors found: " << vectorname.c_str() << endl;
    }

    // as a last resort, recompute the values
    if (valuefile == NULL && vectorfile == NULL)
    {
      cout << " No precomputed eigenvalues found." << endl;
      //generateSparseStiffnessMatrix();
      _basisLoaded = false;
    }
    else
      _basisLoaded = true;
    if (valuefile != NULL) fclose(valuefile);
    if (vectorfile != NULL) fclose(vectorfile);
  }

  // initialize the state vector to the correct size
  //int rank = _eigenvalues.size();
  int rank = _UBasis.cols();
  _q.resizeAndWipe(rank);

  // initialize reduced force vector to the correct size
  //_externalForces.resizeAndWipe(rank);

  // tet address to index lookup hash
  for (unsigned int x = 0; x < _tets.size(); x++)
    _tetID[&_tets[x]] = x;

  // compute the reduced mass matrix 
  // (see section 3.3 of [Barbic and James 2005])
  cout << " Computing reduced mass matrix ... ";
  if (_UBasis.rows() != 0)
  {
    cout <<  "Trying to reactivate the correct mass computation " << endl;
    cout << " basis size: " << _UBasis.dims() << endl;
    MATRIX MU = _masses * _UBasis;
    MATRIX UT = _UBasis.transpose();
    _reducedM = UT * MU;
    /*
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Reduced mass breaks on 64 bit, using brute force" << endl;
    MATRIX MU(_UBasis.cols(), _UBasis.cols());
    MU.setToIdentity();
    MU *= _masses(0,0);
    _reducedM = MU;
    */
  }
  cout << "done" << endl;

  // build the subbasis matrix for just the surface vertices
  cout << " Computing surface U ... ";
  flush(cout);
  //string surfacename = string(filename) + string(".surfaceU.matrix");
  string surfacename("");
  surfacename.assign(filename);
  surfacename.append(".surfaceU.matrix");
  FILE* surfacefile = fopen(surfacename.c_str(), "rb");
  if (surfacefile != NULL)
  {
    cout << " cache found ... "; flush(cout);
    _surfaceU.read(surfacefile);
    fclose(surfacefile);
  }
  else
  {
    cout << " no cache found ... "; flush(cout);
    //fclose(surfacefile);
    if (_UBasis.rows() != 0)
    {
      // build it, then write it out
      buildSurfaceU();
      surfacefile = fopen(surfacename.c_str(), "wb");
      _surfaceU.write(surfacefile);
      fclose(surfacefile);
    }
  }
  cout << "done." << endl;

  // try to read in a cubature
  if (_UBasis.rows() != 0)
    readCubature();

  /*
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " CACHING DEACTIVATED" << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  */

  //if (_UBasis.rows() != 0)
  if (_totalKeyTets != 0)
  {
    cout << " Looking for an inertia cache ..."; flush(cout);
    if (readInertiaCache() == false)
    {
      cout << " not found! " << endl;
      // precompute inertia tensor quantities
      cout << " Caching inertia vars ... "; flush(cout);
      cacheInertiaVars();
      cout << " done." << endl;
      cout << " Caching rotation vars ... "; flush(cout);
      cacheRotationDefoVars();
      cout << " done." << endl;
      cacheMassMatrixVars();

      writeInertiaCache();
    }
    else
    {
      cout << " found! " << endl;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// update the full mesh positions using the subspace vector _q
//////////////////////////////////////////////////////////////////////
void SUBSPACE_TET_MESH::updateFullMesh()
{
  VECTOR update = _UBasis * _q;
  int y = 0;

  for (int x = 0; x < _unconstrainedSize; x++)
  {
    VEC3F& restPose = restVertex(&_vertices[x]);
    _vertices[x][0] = restPose[0] + update(y++);
    _vertices[x][1] = restPose[1] + update(y++);
    _vertices[x][2] = restPose[2] + update(y++);
  }
}

//////////////////////////////////////////////////////////////////////
// speculatively update the mesh with an uncommitted q
//////////////////////////////////////////////////////////////////////
void SUBSPACE_TET_MESH::updateFullMesh(VECTOR& q)
{
  VECTOR update = _UBasis * q;
  int y = 0;

  for (int x = 0; x < _unconstrainedSize; x++)
  {
    VEC3F& restPose = restVertex(&_vertices[x]);
    _vertices[x][0] = restPose[0] + update(y++);
    _vertices[x][1] = restPose[1] + update(y++);
    _vertices[x][2] = restPose[2] + update(y++);
  }
}

//////////////////////////////////////////////////////////////////////
// Extract the subrows of _UBasis corresponding to the surface
// vertices and assemble them in _surfaceU
//
// Assumes that computeSurfaceVertices from the parent class has
// been called.
//////////////////////////////////////////////////////////////////////
void SUBSPACE_TET_MESH::buildSurfaceU()
{
  // total surface vertices
  int totalVertices = _surfaceVertices.size();
  _surfaceU.resizeAndWipe(totalVertices * 3, _q.size());

  for (int x = 0; x < totalVertices; x++)
  {
    MATRIX subU = vertexSubBasis(_surfaceVertices[x]);

    // copy the subbasis to _surfaceU
    Real* surfaceURow = _surfaceU.row(3 * x);
    Real* subURow = subU.data();
    for (int y = 0; y < 3 * _q.size(); y++)
      surfaceURow[y] = subURow[y];
  }
}

//////////////////////////////////////////////////////////////////////
// update the vertices on the mesh surface using the subspace 
// vector _q
//////////////////////////////////////////////////////////////////////
void SUBSPACE_TET_MESH::updateSurfaceMesh()
{
  // total surface vertices
  int totalVertices = _surfaceVertices.size();
  VECTOR update = _surfaceU * _q;
  int y = 0;

  for (int x = 0; x < totalVertices; x++)
  {
    VEC3F& restPose = restVertex(_surfaceVertices[x]);
    VEC3F& deformed = (*_surfaceVertices[x]);
    deformed[0] = restPose[0] + update(y++);
    deformed[1] = restPose[1] + update(y++);
    deformed[2] = restPose[2] + update(y++);
  }
}

//////////////////////////////////////////////////////////////////////
// Compute the projected internal force by computing the full
// internal force and then multiplying it by U^T
//////////////////////////////////////////////////////////////////////
VECTOR SUBSPACE_TET_MESH::projectedInternalForce()
{
  // compute full internal forces
  computeInternalForces();

  // flatten forces out into a single vector
  VECTOR forceVector(_internalForces.size() * 3);
  for (unsigned int x = 0; x < _internalForces.size(); x++)
  {
    forceVector(3 * x)     = _internalForces[x][0];
    forceVector(3 * x + 1) = _internalForces[x][1];
    forceVector(3 * x + 2) = _internalForces[x][2];
  }
  
  // project the full force vector
  // Note the carat (^) means multiply by the transpose of _UBasis
  return _UBasis ^ forceVector;
}

//////////////////////////////////////////////////////////////////////
// get submatrix of _UBasis corresponding to vertices in tet
//////////////////////////////////////////////////////////////////////
MATRIX SUBSPACE_TET_MESH::tetSubBasis(TET* tet)
{
  MATRIX subbasis(12, _q.size());

  // populate the matrix for each vertex
  for (int x = 0; x < 4; x++)
  {
    int index = _vertexID[tet->vertices[x]];

    // if the node is constrained, do nothing since
    // the row was initialized to zero
    if (index >= _unconstrainedSize) continue;

    // multiply by 3 since each vertex has an x,y,z
    Real* fullRow = _UBasis.row(3 * index);
    Real* subbasisRow = subbasis.row(3 * x);

    // copy three rows
    for (int y = 0; y < 3 * _q.size(); y++)
      subbasisRow[y] = fullRow[y];
   
   /* 
    for (int y = 0; y < 3; y++)
    {
      int row = 3 * index + y;
      VECTOR rowData = _UBasis.getRow(row);
      for (int z = 0; z < rowData.size(); z++)
        subbasis(y, z) = rowData[z];
    }
    */
  }

  return subbasis;
}

//////////////////////////////////////////////////////////////////////
// get submatrix of _UBasis corresponding to vertices in tet
//////////////////////////////////////////////////////////////////////
void SUBSPACE_TET_MESH::tetSubBasis(TET* tet, MATRIX &subbasis)
{
  if ( subbasis.rows() != 12 || subbasis.cols() != _q.size() )
    subbasis.resizeAndWipe( 12, _q.size() );

  // populate the matrix for each vertex
  for (int x = 0; x < 4; x++)
  {
    int index = _vertexID[tet->vertices[x]];

    // if the node is constrained, do nothing since
    // the row was initialized to zero
    if (index >= _unconstrainedSize) continue;

    // multiply by 3 since each vertex has an x,y,z
    Real* fullRow = _UBasis.row(3 * index);
    Real* subbasisRow = subbasis.row(3 * x);

    // copy three rows
    for (int y = 0; y < 3 * _q.size(); y++)
      subbasisRow[y] = fullRow[y];
  }
}

//////////////////////////////////////////////////////////////////////
// get submatrix of _UBasis corresponding to a vertex
//////////////////////////////////////////////////////////////////////
MATRIX SUBSPACE_TET_MESH::vertexSubBasis(VEC3F* vertex)
{
  MATRIX subbasis(3, _q.size());

  int index = _vertexID[vertex];
  if (index >= _unconstrainedSize) return subbasis;

  // multiply by 3 since each vertex has an x,y,z
  Real* fullRow = _UBasis.row(3 * index);
  Real* subbasisRow = subbasis.data();

  // copy three rows
  for (int y = 0; y < 3 * _q.size(); y++)
    subbasisRow[y] = fullRow[y];

  return subbasis;
}

//////////////////////////////////////////////////////////////////////
// Same as the above, but use a given matrix
//////////////////////////////////////////////////////////////////////
bool SUBSPACE_TET_MESH::vertexSubBasis(int vertexID, MATRIX &subBasis,
                                       bool transpose) const
{
  if ( vertexID >= _unconstrainedSize ) return false;

  if ( !transpose &&
     ( subBasis.rows() != 3 || subBasis.cols() != _q.size() ) )
  {
    subBasis.resizeAndWipe( 3, _q.size() );
  }
  else if ( transpose && 
    ( subBasis.rows() != _q.size() || subBasis.cols() != 3 ) )
  {
    subBasis.resizeAndWipe( _q.size(), 3 );
  }
  Real *subbasisRow = subBasis.data();

  if ( transpose )
  {
    const Real *fullRow = _UBasis.row(0);

    int index = 0;
    int startRow = 3 * vertexID;
    // Copy along columns instead of rows
    for ( int i = 0; i < _q.size(); i++ )
    {
      for ( int j = 0; j < 3; j++ )
      {
        subbasisRow[ index ] = fullRow[ ( startRow + j ) * _q.size() + i ];
        index++;
      }
    }
  }
  else
  {
    const Real* fullRow = _UBasis.row(3 * vertexID);

    memcpy( (void *)subbasisRow, (void *)fullRow,
             3 * _q.size() * sizeof( Real ) );
  }

  return true;
}

//////////////////////////////////////////////////////////////////////
// read in a scheme generated by CUBATURE_GENERATOR
//////////////////////////////////////////////////////////////////////
void SUBSPACE_TET_MESH::readCubature(const char* filename)
{
  string finalfile("");
  
  // if no filename is provided, generate one based on the mesh name
  if (filename == NULL)
    finalfile = _filename + string(".cubature");
  else
    finalfile = string(filename);

  // read in the cubature file
  FILE* file;
  file = fopen(finalfile.c_str(), "rb");

  // if the file does not exist, just return gracefully
  if (file == NULL) {
    cout << " WARNING: Did not find any cubature file! " << endl;
    return;
  }
  cout << " Reading in cubature file: " << finalfile.c_str() << " ";
  flush(cout);

  // write dimensions
  fread((void*)&_totalKeyTets, sizeof(int), 1, file);
  cout << " total key tets: " << _totalKeyTets << endl;

  cout << "Read in " << _totalKeyTets << " key tets" << endl;

  // allocate the key tet indices
  _keyTets = new int[_totalKeyTets];
  _keyWeights.resizeAndWipe(_totalKeyTets);

  _keyTetStiffnesses.resize( _totalKeyTets );
  for ( int x = 0; x < _totalKeyTets; x++ )
    _keyTetStiffnesses[ x ].resizeAndWipe( _UBasis.cols(), _UBasis.cols() );

  _keyTetDeformations.resize( _totalKeyTets );
  for ( int x = 0; x < _totalKeyTets; x++ )
  {
    _keyTetDeformations[ x ]._F.clear();
    _keyTetDeformations[ x ]._U.clear();
    _keyTetDeformations[ x ]._Fhat.clear();
    _keyTetDeformations[ x ]._V.clear();
  }
  
  // write out the key tet indices
  for (int x = 0; x < _totalKeyTets; x++)
    fread((void*)&(_keyTets[x]), sizeof(int), 1, file);

  // read in the key tet weights
  for (int x = 0; x < _totalKeyTets; x++)
  {
    double weight;
    fread((void*)&(weight), sizeof(double), 1, file);
    _keyWeights(x) = weight;
  }
  
  fclose(file);

  // now that we know how many key tets there are, we can allocate a bunch
  // of Newmark variables
  _Freduced.resizeAndWipe(9 * _totalKeyTets);
  _Rreduced.resizeAndWipe(_totalKeyTets);
  _Kreduced.resizeAndWipe(_totalKeyTets, _totalKeyTets);
  _forceDensity.resizeAndWipe(9 * _totalKeyTets);
  _SH.resizeAndWipe(9 * _totalKeyTets, _q.size());

  // allocate _stiffnessDensity
  _stiffnessDensity = new MATRIX*[_totalKeyTets];
  for (int x = 0; x < _totalKeyTets; x++)
    _stiffnessDensity[x] = new MATRIX(9,9);

  // precompute matrices for Newmark later
  generateE();
  generateH();

  // precache the volumes of the rest state tets
  // (needed by diagonalized generateStiffnessMatrix)
  _keyRestVolumes.resizeAndWipe(_totalKeyTets);
  for (int x = 0; x < _totalKeyTets; x++)
    _keyRestVolumes(x) = 1.0 / _tets[_keyTets[x]].volume();
}

SUBSPACE_TET_MESH::~SUBSPACE_TET_MESH()
{
  for (int x = 0; x < _totalKeyTets; x++)
    delete _stiffnessDensity[x];
  delete[] _stiffnessDensity;

  delete[] _pFpUrepacked;
  delete[] _tetSubBases;
  delete[] _tetSubBasesTransposed;
  delete[] _tetBmatrices;

  if (_keyTets)
    delete[] _keyTets;
}

//////////////////////////////////////////////////////////////////////
// Generate deformation gradients F all at once
//////////////////////////////////////////////////////////////////////
void SUBSPACE_TET_MESH::generateF()
{
  /*
  //_Freduced = _E * _q;
  _E.multiplyInplace(_q, _Freduced);

  // Add identity as described in section 6.1
#pragma omp parallel for schedule(static) default(shared)
  for (int x = 0; x < _totalKeyTets; x++)
  {
    _Freduced(x * 9) += 1.0f;
    _Freduced(x * 9 + 4) += 1.0f;
    _Freduced(x * 9 + 8) += 1.0f;
  }
  */

  /*
  for (int x = 0; x < _totalKeyTets; x++)
  {
    // get the key tet
    TET* tet = &_tets[_keyTets[x]];
    
    // get the rest state
    VEC3F vertices[4];
    for (int y = 0; y < 4; y++)
    {
      int vertexID = _vertexID[tet->vertices[y]];
      vertices[y] = _restPose[vertexID];
    }
    
    // get unconstrained displacements of each vertex
    TET newTet(&vertices[0], &vertices[1], &vertices[2], &vertices[3]);
    for (int y = 0; y < 4; y++)
    {
      if (isConstrained(tet->vertices[y]))
      {
        int vertexID = _vertexID[tet->vertices[y]];
        vertices[y] = _vertices[vertexID];
      }
      else
      {
        MATRIX subU = vertexSubBasis(tet->vertices[y]);
        VECTOR displacements = subU * _q;
   
        vertices[y][0] += displacements[0];
        vertices[y][1] += displacements[1];
        vertices[y][2] += displacements[2];
      }
    }
    
    // write to F reduced
    MATRIX3 F = newTet.F();
    int i = 0;
    for (int y = 0; y < 3; y++)
      for (int z = 0; z < 3; z++, i++)
        _Freduced(x * 9 + i) = F(z,y);
  }
  */

  if (rank() == 0) return;

  if (_x.size() == 0) _x.resizeAndWipe(TET_MESH::rank());

  //this->x() = _UBasis * _q;
  //TET_MESH::updateFullMesh();
  for (int x = 0; x < _totalKeyTets; x++)
  {
    // get the key tet
    TET* tet = &_tets[_keyTets[x]];
   
    // unproject the tet state
    for (int y = 0; y < 4; y++)
    {
      if (!(isConstrained(tet->vertices[y])))
      {
        int vertexID = _vertexID[tet->vertices[y]];
        VEC3F& restPose = _restPose[vertexID];
        VEC3F& deformed = _vertices[vertexID];

        MATRIX subU = vertexSubBasis(tet->vertices[y]);
        VECTOR displacements = subU * _q;
   
        deformed[0] = restPose[0] + displacements[0];
        deformed[1] = restPose[1] + displacements[1];
        deformed[2] = restPose[2] + displacements[2];

        int ID3 = 3 * vertexID;
        _x[ID3] = displacements[0];
        _x[ID3 + 1] = displacements[1];
        _x[ID3 + 2] = displacements[2];
      }
    }
 
    // write to F reduced
    MATRIX3 F = tet->F();
    int i = 0;
    int x9 = x * 9;
    for (int y = 0; y < 3; y++)
      for (int z = 0; z < 3; z++, i++)
        _Freduced(x9 + i) = F(z,y);
  }
}

//////////////////////////////////////////////////////////////////////
// Evaluate reduced internal forces for Newmark using a Krysl-style
// projection
//////////////////////////////////////////////////////////////////////
VECTOR& SUBSPACE_TET_MESH::generateKryslInternalForces()
{
  _x = _UBasis * _q;
  TET_MESH::updateFullMesh();

  VECTOR& Rfull = TET_MESH::generateInternalForces();
  _Rreduced = _UBasis ^ Rfull;

  VECTOR diff = Rfull - _UBasis * _Rreduced;
  return _Rreduced;
}

//////////////////////////////////////////////////////////////////////
// Evaluate reduced internal forces for Newmark
// Assumes generateF() has already been called
//////////////////////////////////////////////////////////////////////
VECTOR& SUBSPACE_TET_MESH::generateInternalForces()
{
//#pragma omp parallel for schedule(static) default(shared)
  for (int x = 0; x < _totalKeyTets; x++)
  {
    TET& tet = _tets[_keyTets[x]];
    int materialIndex = tet.materialIndex();
    _materials[materialIndex]->forceDensity(&_Freduced(x * 9), &_forceDensity(x * 9));
  }

  /*
  // DEBUG
  cout << __FILE__ << " " << __LINE__ << " : " << endl;
  for (int x = 0; x < _totalKeyTets; x++)
  {
    if (_tets[_keyTets[x]].inverted())
    {
      cout << __FILE__ << " " << __LINE__ << " : " << endl;
      cout << " A TET IS INVERTED: " << endl;
      cout << _tets[_keyTets[x]] << endl;
      cout << " F: " << _tets[_keyTets[x]].F() << endl;
      cout << " det(F): " << det(_tets[_keyTets[x]].F()) << endl;
    }
  }
  */

  _Rreduced = _leftH ^ _forceDensity;
  //_RreducedGaussKronrod = _leftHGaussKronrod ^ _forceDensity;
  return _Rreduced;
}

//////////////////////////////////////////////////////////////////////
// Evaluate reduced internal forces for Newmark using a Krysl-style
// projection
//////////////////////////////////////////////////////////////////////
VECTOR& SUBSPACE_TET_MESH::generateKryslInternalForces(vector<MATRIX3>& Us, 
                                                  vector<MATRIX3>& Fhats,
                                                  vector<MATRIX3>& Vs)
{
  VECTOR& fullR = TET_MESH::generateInternalForces(Us, Fhats, Vs);
  _Rreduced = _UBasis ^ fullR;
  return _Rreduced;
}

//////////////////////////////////////////////////////////////////////
// Evaluate unreduced internal forces for Newmark
// and then project them down (for debugging purposes)
// 
// Assumes generateF() has already been called
//////////////////////////////////////////////////////////////////////
VECTOR& SUBSPACE_TET_MESH::generateInternalForcesDebug()
{
  /*
#pragma omp parallel for schedule(static) default(shared)
  for (int x = 0; x < _totalKeyTets; x++)
    _material->forceDensity(&_Freduced(x * 9), &_forceDensity(x * 9));

  _Rreduced = _leftH ^ _forceDensity;
  return _Rreduced;
  */
  _Rreduced.resizeAndWipe(_unconstrainedSize * 3);

  // call unreduced version
  TET_MESH::generateInternalForces();

  VECTOR projected = _UBasis * _Rreduced;

  _Rreduced = projected;
  return _Rreduced;
}

//////////////////////////////////////////////////////////////////////
// Evaluate reduced internal forces for Newmark, but do it in a way
// that calls the invertibility filtering
// Assumes SUBSPACE_INTEGRATOR::cacheDiagonalizations() has already 
// been called, and the results are being passed in
//
// FIXME: Dynamic allocations here.  These are probably quite slow
//////////////////////////////////////////////////////////////////////
VECTOR& SUBSPACE_TET_MESH::generateInternalForces(vector<MATRIX3>& Us, 
                                                  vector<MATRIX3>& Fhats,
                                                  vector<MATRIX3>& Vs)
{
  for (int x = 0; x < _totalKeyTets; x++)
  {
    TET& tet = _tets[_keyTets[x]];
    int materialIndex = tet.materialIndex();
    INVERTIBLE* material = (INVERTIBLE*)_materials[materialIndex];
    
    // get the filtered first PK
    MATRIX3 firstPK = material->firstPiolaKirchhoff(Us[x], Fhats[x], Vs[x]);
    //VECTOR flatPK = TET::flattenF(firstPK);

    VECTOR &flatPK = _flatPKworkspace.get();

    TET::flattenF( firstPK, flatPK );

    // negate 1st PK, as is our convention
    for (int i = 0; i < 9; i++)
      _forceDensity(x * 9 + i) = -flatPK(i);
  }
  //_Rreduced = _leftH ^ _forceDensity;

  if ( _Rreduced.size() != _q.size() )
    _Rreduced.resizeAndWipe( _q.size() );

  _leftH.gemvInplace( 1.0, _forceDensity, _Rreduced, 0.0, true );

  return _Rreduced;
}

//////////////////////////////////////////////////////////////////////
// Evaluate reduced stiffness for Newmark
// Assumes generateF() has already been called
//////////////////////////////////////////////////////////////////////
MATRIX& SUBSPACE_TET_MESH::generateStiffnessMatrix()
{
  // first populate the block diagonal matrix
#if USING_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
  for (int x = 0; x < _totalKeyTets; x++)
  {
    TET& tet = _tets[_keyTets[x]];
    int materialIndex = tet.materialIndex();
    _materials[materialIndex]->stiffnessDensity(
            &( _Freduced(x * 9) ),
            (*_stiffnessDensity[x]).data());
  }

  // multiply by _rightH (the one without the weights)
#if USING_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
  for (int x = 0; x < _totalKeyTets; x++)
  {
    // get the rows corresponding to the xth key tet
    SUBMATRIX H(_rightH, 9 * x, 9);

    // get the temp rows to hold the product of 
    // _stiffnessDensity * _rightH in
    SUBMATRIX SH(_SH, 9 * x, 9);

    // SH = _stiffnessDensity * _rightH
    SH.clearingGemm((*_stiffnessDensity[x]), H);
  }
 
  if (_leftH.rows() != 0)
    _Kreduced = _leftH ^ _SH;

  return _Kreduced;
}

//////////////////////////////////////////////////////////////////////
// Perform a Krysl-style evaluation of the stiffness matrix --
// compute the full one and then project it down
//////////////////////////////////////////////////////////////////////
MATRIX& SUBSPACE_TET_MESH::generateKryslStiffnessMatrix()
{
  SPARSE_MATRIX& K = TET_MESH::generateSparseStiffnessMatrix();

  // project down
  MATRIX RHS = K * _UBasis;
  _Kreduced = _UBasis ^ RHS;
  
  // return projection
  return _Kreduced;
}

//////////////////////////////////////////////////////////////////////
// Perform a Krysl-style evaluation of the stiffness matrix --
// compute the full one and then project it down
//////////////////////////////////////////////////////////////////////
MATRIX& SUBSPACE_TET_MESH::generateKryslStiffnessMatrix(vector<MATRIX3>& Us, 
                                                   vector<MATRIX3>& Vs, 
                                                   vector<MATRIX>& stiffnesses, map<string,double>& timingBreakdown)
{
  TIMER fullTimer;
  SPARSE_MATRIX& K = TET_MESH::generateSparseStiffnessMatrix(Us, Vs, stiffnesses);
  timingBreakdown["Krysl-style full stiffness assembly"] += fullTimer.timing();

  // project down
  TIMER projectionTimer;
  MATRIX RHS = K * _UBasis;
  _Kreduced = _UBasis ^ RHS;
  timingBreakdown["Krysl-style projection"] += projectionTimer.timing();
  
  // return projection
  return _Kreduced;
}

//////////////////////////////////////////////////////////////////////
// Evaluate reduced stiffness for Newmark
// Assumes generateF() has already been called
//////////////////////////////////////////////////////////////////////
MATRIX& SUBSPACE_TET_MESH::generateStiffnessMatrix(vector<MATRIX3>& Us, 
                                                   vector<MATRIX3>& Vs, 
                                                   //vector<MATRIX>& stiffnesses)
                                                   vector<MATRIX>& stiffnesses,
                                                   map<string,double>& timingBreakdown)
{
  _Kreduced.clear();
  
  // first populate the block diagonal matrix
  for (int x = 0; x < _totalKeyTets; x++)
  {
    //TIMER preamble;
    TET& tet = _tets[_keyTets[x]];

    // get cached stiffness matrix
    MATRIX& diagonalStiffness = stiffnesses[x];
    MATRIX3& U = Us[x];
    MATRIX3& V = Vs[x];

    // get Bm
    //const VEC3F* b = tet.b();
    
    // get PFPu from the tet - only depends on rest state,
    // so don't need to update tet state at all
    // FIXME: This can be precomputed (doesn't depend on position)
    MATRIX pFpu(9,12);

    _materials[tet.materialIndex()]->computePFPu(tet, pFpu);
    MATRIX stiffnessSubblock(12,12);

    MATRIX3 Utranspose = U.transpose();
    MATRIX3 Vtranspose = V.transpose();

    //timingBreakdown["S. Assem. Preamble"] += preamble.timing();

    //TIMER probing;
    // compute each column of the matrix
    for (int y = 0; y < 12; y++)
    {
      // extract a column from pFpu (ie pretend 'x' is just a 
      // delta with a single 1 and the rest are zero)
      VECTOR deltaF(9);

      for (int z = 0; z < 9; z++)
        deltaF(z) = pFpu(z, y);

      // rotate deltaF
      // FIXME: Don't need to call U.transpose every tiem
      //MATRIX3 rotated = U.transpose() * TET::repackF(deltaF) * V;
      MATRIX3 rotated = Utranspose * TET::repackF(deltaF) * V;
      deltaF = TET::flattenF(rotated);

      VECTOR contraction = diagonalStiffness * deltaF;

      MATRIX3 deltaP = TET::repackF(contraction);

      // rotate deltaP back
      // FIXME; don't need to transpose V every time
      //deltaP = U * deltaP * V.transpose();
      deltaP = U * deltaP * Vtranspose;
      VEC3F forceVecs[4];
      const VEC3F* b = tet.b();
      forceVecs[0] = deltaP * b[0];
      forceVecs[1] = deltaP * b[1];
      forceVecs[2] = deltaP * b[2];
      forceVecs[3] = deltaP * b[3];
     
      // copy result into stiffness column
      for (int z = 0; z < 4; z++)
        for (int a = 0; a < 3; a++)
          stiffnessSubblock(z * 3 + a, y) = forceVecs[z][a];
    }
    //timingBreakdown["S. Assem. Probing"] += probing.timing();

    // FIXME: Precompute this.
    //TIMER getbasis;
    MATRIX tetbasis= tetSubBasis(&tet);

    //timingBreakdown["S. Assem. get basis"] += getbasis.timing();

    /*
    // do a more detailed breakdown
    TIMER final;
    Real weight = _keyRestVolumes(x) * _keyWeights(x);
    MATRIX subK = stiffnessSubblock * tetbasis;
    subK = tetbasis ^ subK;
    _Kreduced += weight * subK;
    timingBreakdown["S. Assem. Final Mult"] += final.timing();
    */

    //TIMER final;
    Real weight = _keyRestVolumes(x) * _keyWeights(x);
    MATRIX subK = stiffnessSubblock * tetbasis;

    //timingBreakdown["S. Assem. Final Mult"] += final.timing();

    //TIMER carat;
    subK = tetbasis ^ subK;

    //timingBreakdown["S. Assem. Carat"] += carat.timing();
    
    //TIMER add;
    _Kreduced += weight * subK;

    //timingBreakdown["S. Assem. Add"] += add.timing();
  }

  return _Kreduced;
}

//////////////////////////////////////////////////////////////////////
// Fast version for more vectorization
//< assumes generateF() was already called!
//////////////////////////////////////////////////////////////////////
MATRIX& SUBSPACE_TET_MESH::generateStiffnessMatrixFast(
                                     vector<MATRIX3>& Us, 
                                     vector<MATRIX3>& Vs, 
                                     vector<MATRIX>& stiffnesses,
                                     map<string, double>& timingBreakdown)
{
  _Kreduced.clear();

  int pFpUsize = 12 * 9;
  int basisSize = 12 * _q.size();
  int BmatrixSize = 3 * 4;

  if ( !_pFpUrepacked )
  {
    initializeAssemblyWorkspace();
  }

  // Store U, V and their respective transposes
  Real U[9];
  Real V[9];
  Real Utranspose[9];
  Real Vtranspose[9];

  // Workspaces
  Real blockWork1[ 12 * 3 * 3 ];
  Real blockWork2[ 3 * 12 * 3 ];
  Real work9x12_1[ 9 * 12 ];
  Real work9x12_2[ 9 * 12 ];
  Real work12x12_1[ 12 * 12 ];
  Real work12x12_2[ 12 * 12 ];

#if 0
  Real work12xq[ 12 * _q.size() ];
#endif

  Real workSpaceRHS[ 12 * _q.size() * _totalKeyTets ];
  
  // first populate the block diagonal matrix
  for (int x = 0; x < _totalKeyTets; x++)
  {
    //TIMER preamble;
    //TET& tet = _tets[_keyTets[x]];

    Real *work12xq = workSpaceRHS + x * 12 * _q.size();

    // get cached stiffness matrix
    MATRIX& diagonalStiffness = stiffnesses[x];
    MATRIX::copy( U, Us[x] );
    MATRIX::copy( V, Vs[x] );

    // get Bm
    //const VEC3F* b = tet.b();
    
    // get PFPu from the tet - only depends on rest state,
    // so don't need to update tet state at all
    // FIXME: This can be precomputed (doesn't depend on position)
    Real *pFpUrepacked =  _pFpUrepacked + x * pFpUsize;
    Real *B = _tetBmatrices + x * BmatrixSize;

    MATRIX::transpose( Utranspose, U, 3, 3 );
    MATRIX::transpose( Vtranspose, V, 3, 3 );

    TIMER vectorPackTimer;

    // Rotate all columns of pFpU
    MATRIX::gemm( pFpUrepacked, V, blockWork1,
                  12 * 3, 3, 3, 3, false, false );

    // Repack this in to a block row vector
    TET::repackBlockColumnToRow( blockWork1, blockWork2, 12 );

    // Rotate all columns
    MATRIX::gemm( Utranspose, blockWork2, blockWork1,
                  3, 3, 3, 12 * 3, false, false );

    // Flatten in to a 9x12 matrix
    TET::flattenBlockRow( blockWork1, work9x12_1, 12 );

    // Multiply by the stiffness density
    MATRIX::gemm( diagonalStiffness.data(), work9x12_1, work9x12_2,
                  9, 9, 9, 12, false, false );

    // Repack to a block column
    TET::repackBlockColumn( work9x12_2, blockWork1, 12 );

    // Rotate all columns
    MATRIX::gemm( blockWork1, Vtranspose, blockWork2,
                  12 * 3, 3, 3, 3, false, false );

    // Repack this in to a block row vector
    TET::repackBlockColumnToRow( blockWork2, blockWork1, 12 );

    // Rotate all columns
    MATRIX::gemm( U, blockWork1, blockWork2,
                  3, 3, 3, 12 * 3, false, false );

    // Repack in to a block column
    TET::repackBlockRowToColumn( blockWork2, blockWork1, 12 );

    // Multiply on the right by B
    MATRIX::gemm( blockWork1, B, work12x12_1,
                  12 * 3, 3, 3, 4, false, false );

    timingBreakdown["S. Assem. Packing"] += vectorPackTimer.timing();

    TIMER localAssembleTimer;

    // Copy final data in to the local stiffness matrix
    for ( int y = 0; y < 12; y++ )
    {
      Real *forceVecs = work12x12_1 + y * 12;

      for ( int z = 0; z < 4; z++ )
      for ( int a = 0; a < 3; a++ )
      {
        MATRIX::access( work12x12_2, 12, 12, z * 3 + a, y )
          = MATRIX::access( forceVecs, 3, 4, a, z );
      }
    }

    timingBreakdown["S. Assem. Local"] += localAssembleTimer.timing();

    Real *tetbasis = _tetSubBases + x * basisSize;

    TIMER final;

    Real weight = _keyRestVolumes(x) * _keyWeights(x);
    MATRIX::gemm( work12x12_2, tetbasis, work12xq,
                  12, 12, 12, _q.size(), false, false, weight );

    timingBreakdown["S. Assem. Final Mult"] += final.timing();

#if 0
    TIMER carat;
    Real *tetbasisTranspose = _tetSubBasesTransposed + x * basisSize;
    MATRIX::gemm( tetbasisTranspose, work12xq, _Kreduced.data(),
                  _q.size(), 12, 12, _q.size(),
                  false, false, 1.0, 1.0 /* Add to _Kreduced */ );
    timingBreakdown["S. Assem. Carat"] += carat.timing();
#endif
  }

  TIMER carat;
  MATRIX::gemm( _tetSubBasesTransposed, workSpaceRHS, _Kreduced.data(),
                _q.size(), 12 * _totalKeyTets, 12 * _totalKeyTets, _q.size(),
                false, false );
  timingBreakdown["S. Assem. Carat"] += carat.timing();

  return _Kreduced;
}

//////////////////////////////////////////////////////////////////////
// This is similar to the above function, but instead of rebuilding
// the stiffness matrix from scratch, updates the existing stiffness
// matrix according to which key tets are dirty.
// Assumes generateF() was already called!
//////////////////////////////////////////////////////////////////////
MATRIX& SUBSPACE_TET_MESH::updateStiffnessMatrixFast(
                                     vector<MATRIX3>& Us, 
                                     vector<MATRIX3>& Vs, 
                                     vector<MATRIX>& stiffnesses,
                                     vector<bool> &keyTetsDirty,
                                     map<string, double>& timingBreakdown )
{
  // NOTE: We don't clear the stiffness matrix here because we
  // are keeping a running total, rather than reconstructing the
  // whole thing.

  int pFpUsize = 12 * 9;
  int basisSize = 12 * _q.size();
  int BmatrixSize = 3 * 4;

  if ( !_pFpUrepacked )
  {
    initializeAssemblyWorkspace( true /* using partial updates */ );
  }

  // Store U, V and their respective transposes
  Real U[9];
  Real V[9];
  Real Utranspose[9];
  Real Vtranspose[9];

  // Workspaces
  Real blockWork1[ 12 * 3 * 3 ];
  Real blockWork2[ 3 * 12 * 3 ];
  Real work9x12_1[ 9 * 12 ];
  Real work9x12_2[ 9 * 12 ];
  Real work12x12_1[ 12 * 12 ];
  Real work12x12_2[ 12 * 12 ];

  Real workSpaceRHS[ 12 * _q.size() * _totalKeyTets ];
  
  // first populate the block diagonal matrix
  for (int x = 0; x < _totalKeyTets; x++)
  {
    // Only update if we are dirty
    if ( !keyTetsDirty[ x ] )
      continue;

#if 0
    cout << "Processing key tet!" << endl;
#endif

    // Subtract the old contribution for this key tet
    // from the running stiffness matrix total
    _Kreduced -= _keyTetStiffnesses[ x ];

    //TIMER preamble;
    TET& tet = _tets[_keyTets[x]];

    Real *work12xq = workSpaceRHS + x * 12 * _q.size();

    // get cached stiffness matrix
    MATRIX& diagonalStiffness = stiffnesses[x];
    MATRIX::copy( U, Us[x] );
    MATRIX::copy( V, Vs[x] );

    // get Bm
    //const VEC3F* b = tet.b();
    
    // get PFPu from the tet - only depends on rest state,
    // so don't need to update tet state at all
    // FIXME: This can be precomputed (doesn't depend on position)
    Real *pFpUrepacked =  _pFpUrepacked + x * pFpUsize;
    Real *B = _tetBmatrices + x * BmatrixSize;

    MATRIX::transpose( Utranspose, U, 3, 3 );
    MATRIX::transpose( Vtranspose, V, 3, 3 );

    TIMER vectorPackTimer;

    // Rotate all columns of pFpU
    MATRIX::gemm( pFpUrepacked, V, blockWork1,
                  12 * 3, 3, 3, 3, false, false );

    // Repack this in to a block row vector
    TET::repackBlockColumnToRow( blockWork1, blockWork2, 12 );

    // Rotate all columns
    MATRIX::gemm( Utranspose, blockWork2, blockWork1,
                  3, 3, 3, 12 * 3, false, false );

    // Flatten in to a 9x12 matrix
    TET::flattenBlockRow( blockWork1, work9x12_1, 12 );

    // Multiply by the stiffness density
    MATRIX::gemm( diagonalStiffness.data(), work9x12_1, work9x12_2,
                  9, 9, 9, 12, false, false );

    // Repack to a block column
    TET::repackBlockColumn( work9x12_2, blockWork1, 12 );

    // Rotate all columns
    MATRIX::gemm( blockWork1, Vtranspose, blockWork2,
                  12 * 3, 3, 3, 3, false, false );

    // Repack this in to a block row vector
    TET::repackBlockColumnToRow( blockWork2, blockWork1, 12 );

    // Rotate all columns
    MATRIX::gemm( U, blockWork1, blockWork2,
                  3, 3, 3, 12 * 3, false, false );

    // Repack in to a block column
    TET::repackBlockRowToColumn( blockWork2, blockWork1, 12 );

    // Multiply on the right by B
    MATRIX::gemm( blockWork1, B, work12x12_1,
                  12 * 3, 3, 3, 4, false, false );

    timingBreakdown["S. Assem. Packing"] += vectorPackTimer.timing();

    TIMER localAssembleTimer;

    // Copy final data in to the local stiffness matrix
    for ( int y = 0; y < 12; y++ )
    {
      Real *forceVecs = work12x12_1 + y * 12;

      for ( int z = 0; z < 4; z++ )
      for ( int a = 0; a < 3; a++ )
      {
        MATRIX::access( work12x12_2, 12, 12, z * 3 + a, y )
          = MATRIX::access( forceVecs, 3, 4, a, z );
      }
    }

    timingBreakdown["S. Assem. Local"] += localAssembleTimer.timing();

    Real *tetbasis = _tetSubBases + x * basisSize;

    TIMER final;

    Real weight = _keyRestVolumes(x) * _keyWeights(x);
    MATRIX::gemm( work12x12_2, tetbasis, work12xq,
                  12, 12, 12, _q.size(), false, false, weight );

    timingBreakdown["S. Assem. Final Mult"] += final.timing();

    TIMER carat;
    Real *tetbasisTranspose = _tetSubBasesTransposed + x * basisSize;

    // Put the outer produce result in to the storage
    // space for this key tet
    MATRIX::gemm( tetbasisTranspose, work12xq, _keyTetStiffnesses[ x ].data(),
                  _q.size(), 12, 12, _q.size(),
                  false, false, 1.0, 0.0 /* Overwrite */ );

    // Add the result to _Kreduced
    _Kreduced += _keyTetStiffnesses[ x ];

    timingBreakdown["S. Assem. Carat"] += carat.timing();
  }

#if 0
  TIMER carat;
  MATRIX::gemm( _tetSubBasesTransposed, workSpaceRHS, _Kreduced.data(),
                _q.size(), 12 * _totalKeyTets, 12 * _totalKeyTets, _q.size(),
                false, false );
  timingBreakdown["S. Assem. Carat"] += carat.timing();
#endif

  return _Kreduced;
}

//////////////////////////////////////////////////////////////////////
// Compute _E for fast evaluation during Newmark integration
//////////////////////////////////////////////////////////////////////
void SUBSPACE_TET_MESH::generateE()
{
  // this matrix is later used to compute all the deformation
  // gradients, F, from q in one shot. The 3x3 F matrix is flattened
  // into a 9-vector, which is why the size is a factor of 9
  _E.resizeAndWipe(9 * _totalKeyTets, _q.size());
 
  // allocate submatrix once 
  MATRIX sum(3, _q.size());
  for (int x = 0; x < _totalKeyTets; x++)
  {
    TET* tet = &(_tets[_keyTets[x]]);
    const MATRIX3 DmInv = tet->DmInv();

    // get the subbases corresponding to each vertex
    MATRIX subbasis0 = vertexSubBasis(tet->vertices[0]);
    MATRIX subbasis1 = vertexSubBasis(tet->vertices[1]);
    MATRIX subbasis2 = vertexSubBasis(tet->vertices[2]);
    MATRIX subbasis3 = vertexSubBasis(tet->vertices[3]);

    // Compute _E
    // This probably doesn't make a lot of sense --
    // There's no compact matrix notation to express this unfortunately,
    // so you'll probably just have to work the math out on your own if
    // you want to understand it
    subbasis1 -= subbasis0;
    subbasis2 -= subbasis0;
    subbasis3 -= subbasis0;
    for (int y = 0; y < 3; y++)
    {
      sum.clearingAxpy(DmInv(0,y), subbasis1);
      sum.axpy(DmInv(1,y), subbasis2);
      sum.axpy(DmInv(2,y), subbasis3);

      // each tet takes up nine rows, with three rows per
      // subbasis
      _E.setSubmatrix(sum, 9 * x + 3 * y);
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Compute _leftH and _rightH for fast evaluation during
// Newmark integration
//////////////////////////////////////////////////////////////////////
void SUBSPACE_TET_MESH::generateH()
{
  // temp matrices used to compute subblocks of leftH and rightH
  MATRIX leftSubmatrix(9, _q.size());
  MATRIX rightSubmatrix(9, _q.size());

  _leftH.resizeAndWipe(9 * _totalKeyTets, _q.size());
  _leftHGaussKronrod.resizeAndWipe(9 * _totalKeyTets, _q.size());
  _rightH.resizeAndWipe(9 * _totalKeyTets, _q.size());
  for (int x = 0; x < _totalKeyTets; x++)
  {
    TET* tet = &(_tets[_keyTets[x]]);
    // \frac{\partial {\bf F}}{\partial {\bf u}}
    MATRIX PFPu(9,12);
    _materials[tet->materialIndex()]->computePFPu(*tet, PFPu);

    // compute the H for this tet
    MATRIX subbasis = tetSubBasis(tet);
    leftSubmatrix.clearingGemm(_keyWeights(x), PFPu, subbasis);

    // note that the -1 is here because our stiffness matrix is
    // negated, and has to be un-negated at some point, so it might
    // as well be here
    rightSubmatrix.clearingGemm(-1.0f, PFPu, subbasis);

    // copy the block into H
    _leftH.setSubmatrix(leftSubmatrix, 9 * x);
    _rightH.setSubmatrix(rightSubmatrix, 9 * x);

    // for first half set weights for error estimator
    if (x < _totalKeyTets / 2 && _kronrodWeights.size() > 0)
    {
      leftSubmatrix.clearingGemm(_kronrodWeights(x), PFPu, subbasis);
      _leftHGaussKronrod.setSubmatrix(leftSubmatrix, 9 * x);
    }
  }
}

//////////////////////////////////////////////////////////////////////
// check if any of the key tets are inverted
//////////////////////////////////////////////////////////////////////
bool SUBSPACE_TET_MESH::keyInverted()
{
  for (int x = 0; x < _totalKeyTets; x++)
  {
    MATRIX3 check(&_Freduced(x * 9));
    if (det(check) < 0.0f) return true;
  }
  return false;
}

//////////////////////////////////////////////////////////////////////
// combine the vertex positions in all the partitions and 
// write out to Matlab format
//////////////////////////////////////////////////////////////////////
void SUBSPACE_TET_MESH::verticesToMatlab(const char* filename, const char* varName, int column)
{
  // update the mesh just in case
  updateFullMesh();
  
  // open a file or writing
  FILE* file = NULL;
  if (column == 0)
    file = fopen(filename, "w");
  else
    file = fopen(filename, "a");
  fprintf(file, "%s(:,%i) = [\n", varName, column + 1);
  for (unsigned int x = 0; x < _vertices.size(); x++)
    fprintf(file, "%f\n%f\n%f\n", _vertices[x][0], _vertices[x][1], _vertices[x][2]);
  fprintf(file, "];\n");
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// combine the vertex positions in all the partitions and 
// write out to binary
//////////////////////////////////////////////////////////////////////
void SUBSPACE_TET_MESH::verticesToBinary(const char* filename, int column)
{
  // update the mesh just in case
  updateFullMesh();
  
  // open a file or writing
  FILE* file = NULL;
  if (column == 0)
  {
    file = fopen(filename, "wb");
    int totalScalars = _vertices.size() * 3;
    fwrite((void*)&totalScalars, sizeof(int), 1, file);
  }
  else
    file = fopen(filename, "ab");
  for (unsigned int x = 0; x < _vertices.size(); x++)
    for (int y = 0; y < 3; y++)
    {
      Real scalar = _vertices[x][y];
      fwrite((void*)&scalar, sizeof(Real), 1, file);
    }
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// update mesh after the basis has changed
//////////////////////////////////////////////////////////////////////
void SUBSPACE_TET_MESH::updateBasis(MATRIX& newBasis)
{
  // change variables over to new basis
  if (_UBasis.cols() > 0)
  {
    MATRIX changeOfBasis = newBasis ^ _UBasis;

    if (_q.size() > 0)
      _q = changeOfBasis * _q;
    else
      _q.resizeAndWipe(newBasis.cols());
    if (_qOld.size() > 0)
      _qOld = changeOfBasis * _qOld;
    else
      _qOld.resizeAndWipe(newBasis.cols());
  }
  else
  {
    //_q.resizeAndWipe(1);
    //_qOld.resizeAndWipe(1);
    _q.resizeAndWipe(newBasis.cols());
    _qOld.resizeAndWipe(newBasis.cols());
  }
 
  // change basis
  _UBasis = newBasis;
  
  // compute the reduced mass matrix 
  // (see section 3.3 of [Barbic and James 2005])
  cout << " Recomputing reduced mass matrix ... ";
  MATRIX MU(_UBasis);
  for (int y = 0; y < MU.cols(); y++)
    for (int x = 0; x < MU.rows(); x++)
      MU(x,y) *= _masses(x,x);
  _reducedM = _UBasis ^ MU;
  cout << "done" << endl;

  // build the subbasis matrix for just the surface vertices
  cout << " Recomputing surface triangles ... ";
  buildSurfaceU();
  cout << "done." << endl;

  // if K hasn't been resized, do it here
  if (_Kreduced.cols() != newBasis.cols())
    _Kreduced.resizeAndWipe(newBasis.cols(), newBasis.cols());
}

//////////////////////////////////////////////////////////////////////
// update cubature quantities after key tets changed
//////////////////////////////////////////////////////////////////////
void SUBSPACE_TET_MESH::updateCubature()
{
  // stomp the old stiffness density
  if (_stiffnessDensity)
  {
    for (int x = 0; x < _totalKeyTets; x++)
      delete _stiffnessDensity[x];
    delete[] _stiffnessDensity;
  }
  
  // cache the new number of key tets
  _totalKeyTets = _keyWeights.size();
  
  // now that we know how many key tets there are, we can allocate a bunch
  // of Newmark variables
  _Freduced.resizeAndWipe(9 * _totalKeyTets);
  _Rreduced.resizeAndWipe(_totalKeyTets);
  //_Kreduced.resizeAndWipe(_totalKeyTets, _totalKeyTets);
  _Kreduced.resizeAndWipe(_q.size(), _q.size());
  _forceDensity.resizeAndWipe(9 * _totalKeyTets);
  _SH.resizeAndWipe(9 * _totalKeyTets, _q.size());

  // allocate _stiffnessDensity
  _stiffnessDensity = new MATRIX*[_totalKeyTets];
  for (int x = 0; x < _totalKeyTets; x++)
    _stiffnessDensity[x] = new MATRIX(9,9);

  // precompute matrices for Newmark later
  generateE();
  generateH();

  // precache the volumes of the rest state tets
  // (needed by diagonalized generateStiffnessMatrix)
  _keyRestVolumes.resizeAndWipe(_totalKeyTets);
  for (int x = 0; x < _totalKeyTets; x++)
  {
    VEC3F* vertices[4];
    for (int y = 0; y < 4; y++)
    {
      int vertexID = _vertexID[_tets[_keyTets[x]].vertices[y]];
      vertices[y] = &_restPose[vertexID];
    }
    
    // build a new tet at the rest state
    TET newTet(vertices[0], vertices[1], vertices[2], vertices[3]);
    
    _keyRestVolumes(x) = 1.0 / newTet.volume();  
  }
}

//////////////////////////////////////////////////////////////////////
// compute error of the cubature using Gauss-Kronrod
//////////////////////////////////////////////////////////////////////
Real SUBSPACE_TET_MESH::cubatureError()
{
  VECTOR diff = _RreducedGaussKronrod - _Rreduced;
  return diff.norm2() / _Rreduced.norm2();
}

//////////////////////////////////////////////////////////////////////
// Blow away the entire basis and start from scratch
//////////////////////////////////////////////////////////////////////
void SUBSPACE_TET_MESH::resetBasis()
{
  // stomp state
  _q.resizeAndWipe(0);
  _qOld.resizeAndWipe(0);
 
  // stomp basis
  _UBasis.resizeAndWipe(0,0);
 
  // stomp cubature
  delete[] _keyTets;
  _keyTets = NULL, 
  _keyWeights.resizeAndWipe(0);
  _kronrodWeights.resizeAndWipe(0);
  _keyRestVolumes.resizeAndWipe(0);

  // stomp anything derived from cubature
  _reducedM.resizeAndWipe(0,0);
  for (int x = 0; x < _totalKeyTets; x++)
    delete _stiffnessDensity[x];
  delete[] _stiffnessDensity;
  _stiffnessDensity = NULL;
  _Freduced.resizeAndWipe(0);
  _Rreduced.resizeAndWipe(0);
  _RreducedGaussKronrod.resizeAndWipe(0);
  _Kreduced.resizeAndWipe(0, 0);
  _forceDensity.resizeAndWipe(0);
  _SH.resizeAndWipe(0,0); 
  _surfaceU.resizeAndWipe(0,0);
  _E.resizeAndWipe(0,0);
  _leftH.resizeAndWipe(0,0);
  _rightH.resizeAndWipe(0,0);
  _leftHGaussKronrod.resizeAndWipe(0,0);

  _totalKeyTets = 0;
}

//////////////////////////////////////////////////////////////////////
// draw a basis vector
//////////////////////////////////////////////////////////////////////
void SUBSPACE_TET_MESH::drawBasis(int whichBasis, Real amplitude)
{
  // backup state
  VECTOR xCopy = _x;

  for (int x = 0; x < _x.size(); x++)
    _x[x] = _UBasis(x, whichBasis) * amplitude;

  TET_MESH::updateSurfaceMesh();
  drawSurfaceFaces();

  // restore backup
  _x = xCopy;
  updateSurfaceMesh();
}

//////////////////////////////////////////////////////////////////////
// write out the state
//////////////////////////////////////////////////////////////////////
void SUBSPACE_TET_MESH::writeState(string filename)
{
  cout << " Dumping state to file " << filename.c_str() << endl;
  FILE* file = fopen(filename.c_str(), "wb");
  _UBasis.write(file);
  _q.write(file);
  _qOld.write(file);
  fwrite((void*)&_totalKeyTets, sizeof(int), 1, file);

  for (int x = 0; x < _totalKeyTets; x++)
    fwrite((void*)&(_keyTets[x]), sizeof(int), 1, file);
  _keyWeights.write(file);
  _kronrodWeights.write(file);
  _keyRestVolumes.write(file);
  _surfaceU.write(file);
  _reducedM.write(file);
  _forceDensity.write(file);
  _SH.write(file);
  _E.write(file);
  _leftH.write(file);
  _rightH.write(file);
  _leftHGaussKronrod.write(file);
  _Freduced.write(file);
  _Rreduced.write(file);
  _Kreduced.write(file);
  _RreducedGaussKronrod.write(file);
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// read in the state
//////////////////////////////////////////////////////////////////////
void SUBSPACE_TET_MESH::readState(string filename)
{
  cout << " Reading state from file " << filename.c_str() << endl;
  FILE* file = fopen(filename.c_str(), "rb");
  _UBasis.read(file);
  _q.read(file);
  _qOld.read(file);
  fread((void*)&_totalKeyTets, sizeof(int), 1, file);
  delete[] _keyTets;
  _keyTets = new int[_totalKeyTets];
  for (int x = 0; x < _totalKeyTets; x++)
    fread((void*)&(_keyTets[x]), sizeof(int), 1, file);
  _keyWeights.read(file);
  _kronrodWeights.read(file);
  _keyRestVolumes.read(file);
  _surfaceU.read(file);
  _reducedM.read(file);
  _forceDensity.read(file);
  _SH.read(file);
  _E.read(file);
  _leftH.read(file);
  _rightH.read(file);
  _leftHGaussKronrod.read(file);
  _Freduced.read(file);
  _Rreduced.read(file);
  _Kreduced.read(file);
  _RreducedGaussKronrod.read(file);
  fclose(file);

  // stomp the old stiffness density
  if (_stiffnessDensity)
  {
    for (int x = 0; x < _totalKeyTets; x++)
      delete _stiffnessDensity[x];
    delete[] _stiffnessDensity;
  }

  // allocate _stiffnessDensity
  _stiffnessDensity = new MATRIX*[_totalKeyTets];
  for (int x = 0; x < _totalKeyTets; x++)
    _stiffnessDensity[x] = new MATRIX(9,9);
}

//////////////////////////////////////////////////////////////////////
// Output colored test to RenderMan
//////////////////////////////////////////////////////////////////////
void SUBSPACE_TET_MESH::drawKeyTetsToRenderMan()
{
#ifdef USING_RENDERMAN
  map<int, bool> isKeyTet;
  for (int x = 0; x < _totalKeyTets; x++)
    isKeyTet[_keyTets[x]] = true;

  cout << " Total Key tets: " << _totalKeyTets << endl;

  Real sliceZ = 1.0;

  // find the maximum density
  Real maxMagnitude = 0.0;
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    // draw a slice, not all the tets
    bool found = false;
    for (int x = 0; x < 4; x++)
      if ((*_tets[i].vertices[x])[2] > sliceZ)
        found = true; 
    if (found) continue;

    // get the force density
    MATERIAL* material = _materials[_tets[i].materialIndex()];
    MATRIX3 F = _tets[i].F();
    MATRIX3 firstPK = material->firstPiolaKirchhoff(F);
    Real magnitude = 0.0;
    for (int x = 0; x < 3; x++)
      for (int y = 0; y < 3; y++)
        magnitude += firstPK(x,y) * firstPK(x,y);
    magnitude = sqrt(magnitude);

    if (magnitude > maxMagnitude)
      maxMagnitude = magnitude;
  }

  // draw on a per-tet basis
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    // draw a slice, not all the tets
    bool found = false;
    for (int x = 0; x < 4; x++)
      if ((*_tets[i].vertices[x])[2] > sliceZ)
        found = true; 
    if (found) continue;

    // only draw key tets
    if (isKeyTet.find(i) == isKeyTet.end())
      continue;

    // get the force density
    MATERIAL* material = _materials[_tets[i].materialIndex()];
    MATRIX3 F = _tets[i].F();
    MATRIX3 firstPK = material->firstPiolaKirchhoff(F);
    Real magnitude = 0.0;
    for (int x = 0; x < 3; x++)
      for (int y = 0; y < 3; y++)
        magnitude += firstPK(x,y) * firstPK(x,y);
    magnitude = sqrt(magnitude);
    magnitude *= 1.0 / maxMagnitude;
    magnitude += 0.00001;

    // get a black run
    //magnitude = 0.0;

    RtColor densityColor;
    rainbowRamp(magnitude, densityColor);
    densityColor[0] = 10.0;
    densityColor[1] = 10.0;
    densityColor[2] = 0.0;
    RiColor(densityColor);

    // init point locations and normals
    RtPoint* P = new RtPoint[4];

    for (int x = 0; x < 4; x++)
    {
      VEC3F vertex = *(_tets[i].vertices[x]);

      P[x][0] = vertex[0];
      P[x][1] = vertex[1];
      P[x][2] = vertex[2];
    }

    // init all to triangles
    int nfaces = 4;
    RtInt* nvertices = new RtInt[nfaces];
    for (int x = 0; x < nfaces; x++)
      nvertices[x] = 3;
    
    // init faces
    RtInt* rifaces = new RtInt[3 * nfaces];
    rifaces[0] = 0;
    rifaces[1] = 1;
    rifaces[2] = 3;
    rifaces[3] = 1;
    rifaces[4] = 2;
    rifaces[5] = 3;
    rifaces[6] = 0;
    rifaces[7] = 3;
    rifaces[8] = 2;
    rifaces[9] = 1;
    rifaces[10] = 0;
    rifaces[11] = 2;
    
    RiPointsPolygons(nfaces, nvertices, rifaces, RI_P, P,  RI_NULL);

    delete[] nvertices;
    delete[] P;
    delete[] rifaces;
  }

  /*
  // draw outlines on a per-tet basis
  RtColor black;
  black[0] = 0.0;
  black[1] = 0.0;
  black[2] = 0.0;
  RiColor(black);
  for (int i = 0; i < _tets.size(); i++)
  {
    // draw a slice, not all the tets
    bool found = false;
    for (int x = 0; x < 4; x++)
      if ((*_tets[i].vertices[x])[2] > sliceZ)
        found = true; 
    if (found) continue;

    // init all to triangles
    int nfaces = 4;
    RtInt* nvertices = new RtInt[nfaces];
    for (int x = 0; x < nfaces; x++)
      nvertices[x] = 3;

    // init point locations and normals
    RtPoint* P = new RtPoint[nfaces * 3];
    Real offset = 0.001;
    for (int y = 0; y < 4; y++)
    {
      TRIANGLE triangle = _tets[i].face(y);
      VEC3F v0 = *(triangle.vertex(0));
      VEC3F v1 = *(triangle.vertex(1));
      VEC3F v2 = *(triangle.vertex(2));

      P[0][0] = v0[0];
      P[0][1] = v0[1] + offset;
      P[0][2] = v0[2];
      P[1][0] = v1[0];
      P[1][1] = v1[1] + offset;
      P[1][2] = v1[2];
      P[2][0] = v2[0];
      P[2][1] = v2[1] + offset;
      P[2][2] = v2[2];
    }

    RtFloat width = 0.002;
    RtToken type = "linear";
    RtToken wrap = "periodic";
    RiDeclare("width", "constant float");
    RiCurves(type, nfaces, nvertices, wrap, RI_P, (RtPointer)P, RI_CONSTANTWIDTH, &width ,RI_NULL);
     
    delete[] P; 
    delete[] nvertices;
  }
  */

#endif
}

//////////////////////////////////////////////////////////////////////
// reset the mass matrix to a specific mass
//////////////////////////////////////////////////////////////////////
void SUBSPACE_TET_MESH::resetMasses(Real mass)
{
  for (int x = 0; x < _masses.rows(); x++)
    _masses(x,x) = mass;

  cout << " Reduced mass breaks on 64 bit, using brute force" << endl;
  MATRIX MU(_UBasis.cols(), _UBasis.cols());
  MU.setToIdentity();
  assert(_masses.rows() > 0 && _masses.cols() > 0);
  MU *= _masses(0,0);
  _reducedM = MU;

  _totalMass = _masses.sum() / 3;
}

//////////////////////////////////////////////////////////////////////
// reset state to zero
//////////////////////////////////////////////////////////////////////
void SUBSPACE_TET_MESH::reset()
{
  _q *= 0.0;
  _qOld *= 0.0;

  updateFullMesh();
}

//////////////////////////////////////////////////////////////////////
// Initializes workspaces for fast reduced stiffness assembly
//////////////////////////////////////////////////////////////////////
void SUBSPACE_TET_MESH::initializeAssemblyWorkspace( bool partialUpdate )
{
  _pFpUrepacked = new Real[ 12 * 9 * _totalKeyTets ];
  _tetSubBases = new Real[ 12 * _q.size() * _totalKeyTets ];
  _tetSubBasesTransposed = new Real[ 12 * _q.size() * _totalKeyTets ];
  _tetBmatrices = new Real[ 3 * 4 * _totalKeyTets ];

  int pFpUsize = 12 * 9;
  int basisSize = 12 * _q.size();
  int BmatrixSize = 3 * 4;

  MATRIX pFpu(9,12);
  MATRIX bMatrix(3,4);

  for ( int i = 0; i < _totalKeyTets; i++ )
  {
    TET& tet = _tets[_keyTets[i]];

    _materials[tet.materialIndex()]->computePFPu( tet, pFpu );

    // Put pFpU in a stacked block vector format, in which each
    // column of the original matrix is repacked in to a 3x3 matrix.
    for ( int j = 0; j < 12; j++ )
    {
      Real *pFpUdata = _pFpUrepacked + i * pFpUsize + j * 9;

      MATRIX::access( pFpUdata, 3, 3, 0, 0 ) = pFpu( 0, j );
      MATRIX::access( pFpUdata, 3, 3, 1, 0 ) = pFpu( 1, j );
      MATRIX::access( pFpUdata, 3, 3, 2, 0 ) = pFpu( 2, j );
      MATRIX::access( pFpUdata, 3, 3, 0, 1 ) = pFpu( 3, j );
      MATRIX::access( pFpUdata, 3, 3, 1, 1 ) = pFpu( 4, j );
      MATRIX::access( pFpUdata, 3, 3, 2, 1 ) = pFpu( 5, j );
      MATRIX::access( pFpUdata, 3, 3, 0, 2 ) = pFpu( 6, j );
      MATRIX::access( pFpUdata, 3, 3, 1, 2 ) = pFpu( 7, j );
      MATRIX::access( pFpUdata, 3, 3, 2, 2 ) = pFpu( 8, j );
    }

    MATRIX subBasis = tetSubBasis( &tet );

    memcpy( (void *)( _tetSubBases + i * basisSize ), subBasis.data(),
             basisSize * sizeof( Real ) );

    subBasis = subBasis.transpose();

    if ( partialUpdate )
    {
      // If we are doing inexact Newton solves we need to be able to
      // easily get at independent transposed sub-bases, so we just
      // stack them on top of one another
      memcpy( (void *)( _tetSubBasesTransposed + i * basisSize ), subBasis.data(),
               basisSize * sizeof( Real ) );
    }
    else
    {
      // Put the bases side by side - eg. in a row vector
      for ( int x = 0; x < subBasis.rows(); x++ )
      for ( int y = 0; y < subBasis.cols(); y++ )
      {
        int row_idx = x;
        int col_idx = i * subBasis.cols() + y;

        MATRIX::access( _tetSubBasesTransposed,
                        subBasis.rows(), _totalKeyTets * subBasis.cols(),
                        row_idx, col_idx ) = subBasis( x, y );
      }
    }

    const VEC3F *b = tet.b();

    for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < 3; i++ )
    {
      bMatrix( i, j ) = b[j][i];
    }

    memcpy( (void *)( _tetBmatrices + i * BmatrixSize ), bMatrix.data(),
            BmatrixSize * sizeof( Real ) );
  }
}

//////////////////////////////////////////////////////////////////////
// precompute all the above inertia quantities
//////////////////////////////////////////////////////////////////////
void SUBSPACE_TET_MESH::cacheInertiaVars()
{
  if (_UBasis.rows() == 0) return;

  vector<VEC3F>& vertices = _vertices;

  _restXrestX = 0.0;
  _restYrestY = 0.0;
  _restZrestZ = 0.0;

  _restXrestY = 0.0;
  _restXrestZ = 0.0;
  _restYrestZ = 0.0;

  _restXUX.resizeAndWipe(_q.size());
  _restYUY.resizeAndWipe(_q.size());
  _restZUZ.resizeAndWipe(_q.size());
  
  _restXUY.resizeAndWipe(_q.size());
  _restXUZ.resizeAndWipe(_q.size());
  _restYUZ.resizeAndWipe(_q.size());

  _restYUX.resizeAndWipe(_q.size());
  _restZUX.resizeAndWipe(_q.size());
  _restZUY.resizeAndWipe(_q.size());

  VECTOR restYUX(_q.size());
  VECTOR restZUX(_q.size());
  VECTOR restZUY(_q.size());

  MATRIX UXs(vertices.size(), _q.size());
  MATRIX UYs(vertices.size(), _q.size());
  MATRIX UZs(vertices.size(), _q.size());

  for (int x = 0; x < _unconstrainedSize; x++)
  {
    Real mass = this->mass(x);
    Real xRest = _restPose[x][0];
    Real yRest = _restPose[x][1];
    Real zRest = _restPose[x][2];

    VECTOR xRow = _UBasis.getRow(3 * x);
    VECTOR yRow = _UBasis.getRow(3 * x + 1);
    VECTOR zRow = _UBasis.getRow(3 * x + 2);
    SUBMATRIX xSubmatrix(_UBasis, 3 * x, 1);
    SUBMATRIX ySubmatrix(_UBasis, 3 * x + 1, 1);
    SUBMATRIX zSubmatrix(_UBasis, 3 * x + 2, 1);

    _restXrestX += mass * _restPose[x][0] * _restPose[x][0];
    _restYrestY += mass * _restPose[x][1] * _restPose[x][1];
    _restZrestZ += mass * _restPose[x][2] * _restPose[x][2];

    _restXrestY += mass * _restPose[x][0] * _restPose[x][1];
    _restXrestZ += mass * _restPose[x][0] * _restPose[x][2];
    _restYrestZ += mass * _restPose[x][1] * _restPose[x][2];

    _restXUX += mass * xRest * xRow;
    _restYUY += mass * yRest * yRow;
    _restZUZ += mass * zRest * zRow;
    
    _restXUY += mass * xRest * yRow;
    _restYUX += mass * yRest * xRow;

    _restXUZ += mass * xRest * zRow;
    _restZUX += mass * zRest * xRow;

    _restYUZ += mass * yRest * zRow;
    _restZUY += mass * zRest * yRow;

    MATRIX scaledX = sqrt(mass) * xSubmatrix;
    MATRIX scaledY = sqrt(mass) * ySubmatrix;
    MATRIX scaledZ = sqrt(mass) * zSubmatrix;

    UXs.setSubmatrix(scaledX, x);
    UYs.setSubmatrix(scaledY, x);
    UZs.setSubmatrix(scaledZ, x);
  }

  _UXTUX = UXs ^ UXs;
  _UYTUY = UYs ^ UYs;
  _UZTUZ = UZs ^ UZs;

  _UXTUY = UXs ^ UYs;
  _UXTUZ = UXs ^ UZs;
  _UYTUZ = UYs ^ UZs;

  _restYUZrestZUY = (_restYUZ + _restZUY);
  _restXUZrestZUX = (_restXUZ + _restZUX);
  _restXUYrestYUX = (_restXUY + _restYUX);

  _workspaceL.resizeAndWipe(_q.size());
  /*
  // pre-cache the 2 in front of the squared terms
  _restXUX *= 2.0;
  _restYUY *= 2.0;
  _restZUZ *= 2.0;
  */
}

//////////////////////////////////////////////////////////////////////
// refresh the inertia tensor
//////////////////////////////////////////////////////////////////////
const MATRIX& SUBSPACE_TET_MESH::refreshInertiaTensor()
{
  // compute tensor using reduced coordinates
  _UXTUX.gemvInplace(1.0, _q, _workspaceL, 0.0);
  Real xSq = _restXrestX + 2.0 * (_restXUX * _q) + _q * _workspaceL;
  _UYTUY.gemvInplace(1.0, _q, _workspaceL, 0.0);
  Real ySq = _restYrestY + 2.0 * (_restYUY * _q) + _q * _workspaceL;
  _UZTUZ.gemvInplace(1.0, _q, _workspaceL, 0.0);
  Real zSq = _restZrestZ + 2.0 * (_restZUZ * _q) + _q * _workspaceL;
  //Real xSq = _restXrestX + 2.0 * (_restXUX * _q) + _q * (_UXTUX * _q);
  //Real ySq = _restYrestY + 2.0 * (_restYUY * _q) + _q * (_UYTUY * _q);
  //Real zSq = _restZrestZ + 2.0 * (_restZUZ * _q) + _q * (_UZTUZ * _q);

  _UYTUZ.gemvInplace(1.0, _q, _workspaceL, 0.0);
  Real yz = _restYrestZ + _restYUZrestZUY * _q + _q * _workspaceL;
  _UXTUZ.gemvInplace(1.0, _q, _workspaceL, 0.0);
  Real xz = _restXrestZ + _restXUZrestZUX * _q + _q * _workspaceL;
  _UXTUY.gemvInplace(1.0, _q, _workspaceL, 0.0);
  Real xy = _restXrestY + _restXUYrestYUX * _q + _q * _workspaceL;
  //Real yz = _restYrestZ + _restYUZrestZUY * _q + _q * (_UYTUZ * _q);
  //Real xz = _restXrestZ + _restXUZrestZUX * _q + _q * (_UXTUZ * _q);
  //Real xy = _restXrestY + _restXUYrestYUX * _q + _q * (_UXTUY * _q);

  _inertiaTensor(0,0) = ySq + zSq;
  _inertiaTensor(1,1) = xSq + zSq;
  _inertiaTensor(2,2) = xSq + ySq;

  _inertiaTensor(0,1) = _inertiaTensor(1,0) = -xy;
  _inertiaTensor(0,2) = _inertiaTensor(2,0) = -xz;
  _inertiaTensor(1,2) = _inertiaTensor(2,1) = -yz;

  return _inertiaTensor;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
const MATRIX& SUBSPACE_TET_MESH::refreshInertiaTensorDt(VECTOR& velocity)
{
  assert(velocity.size() == _q.size());
  _inertiaTensorDt.clear();

  _UXTUX.gemvInplace(1.0, velocity, _workspaceL, 0.0);
  Real XdotX = _restXUX * velocity + _q * _workspaceL;
  _UYTUY.gemvInplace(1.0, velocity, _workspaceL, 0.0);
  Real YdotY = _restYUY * velocity + _q * _workspaceL;
  _UZTUZ.gemvInplace(1.0, velocity, _workspaceL, 0.0);
  Real ZdotZ = _restZUZ * velocity + _q * _workspaceL;
  //Real XdotX = _restXUX * velocity + _q * (_UXTUX) * velocity;
  //Real YdotY = _restYUY * velocity + _q * (_UYTUY) * velocity;
  //Real ZdotZ = _restZUZ * velocity + _q * (_UZTUZ) * velocity;

  _inertiaTensorDt(0,0) = 2 * (YdotY + ZdotZ);
  _inertiaTensorDt(1,1) = 2 * (XdotX + ZdotZ);
  _inertiaTensorDt(2,2) = 2 * (XdotX + YdotY);

  _UXTUY.gemvInplace(1.0, velocity, _workspaceL, 0.0, true);
  Real YdotX = _restYUX * velocity + _q * _workspaceL;
  _UXTUY.gemvInplace(1.0, velocity, _workspaceL, 0.0);
  Real XdotY = _restXUY * velocity + _q * _workspaceL;
  //Real YdotX = _restYUX * velocity + _q * (_UXTUY ^ velocity);
  //Real XdotY = _restXUY * velocity + _q * (_UXTUY * velocity);
  _inertiaTensorDt(0,1) -= YdotX + XdotY;

  _UXTUZ.gemvInplace(1.0, velocity, _workspaceL, 0.0, true);
  Real ZdotX = _restZUX * velocity + _q * _workspaceL;
  _UXTUZ.gemvInplace(1.0, velocity, _workspaceL, 0.0);
  Real XdotZ = _restXUZ * velocity + _q * _workspaceL;
  //Real ZdotX = _restZUX * velocity + _q * (_UXTUZ ^ velocity);
  //Real XdotZ = _restXUZ * velocity + _q * (_UXTUZ * velocity);
  _inertiaTensorDt(0,2) -= ZdotX + XdotZ;

  _UYTUZ.gemvInplace(1.0, velocity, _workspaceL, 0.0, true);
  Real ZdotY = _restZUY * velocity + _q * _workspaceL;
  _UYTUZ.gemvInplace(1.0, velocity, _workspaceL, 0.0);
  Real YdotZ = _restYUZ * velocity + _q * _workspaceL;
  //Real ZdotY = _restZUY * velocity + _q * (_UYTUZ ^ velocity);
  //Real YdotZ = _restYUZ * velocity + _q * (_UYTUZ * velocity);
  _inertiaTensorDt(1,2) -= ZdotY + YdotZ;

  _inertiaTensorDt(1,0) = _inertiaTensorDt(0,1);
  _inertiaTensorDt(2,0) = _inertiaTensorDt(0,2);
  _inertiaTensorDt(2,1) = _inertiaTensorDt(1,2);

  return _inertiaTensorDt;
}

//////////////////////////////////////////////////////////////////////
// precompute all the rotation-defo tensor quantities -- must be public so 
// PARTITIONED_SUBSPACE_TET_MESH can call it after altering all the masses
//////////////////////////////////////////////////////////////////////
void SUBSPACE_TET_MESH::cacheRotationDefoVars()
{
  //if (_UBasis.rows() == 0) return;
  if (_totalKeyTets == 0) return;

  // reduced products
  _restIthetaf.resizeAndWipe(3, _UBasis.cols());
  _UijTildeUij.resizeAndWipe(3, _UBasis.cols(), _UBasis.cols());

  for (int x = 0; x < _UBasis.rows() / 3; x++)
  {
    Real mass = this->mass(x);
    MATRIX restBarTilde = MATRIX::cross(_restPose[x]);
    SUBMATRIX submatrix(_UBasis, 3 * x, 3);

    // rest pose term
    _restIthetaf += mass * restBarTilde * submatrix;

    TENSOR3 UijTilde = TENSOR3::cross(submatrix);
    UijTilde *= mass;
    _UijTildeUij += UijTilde * submatrix;
  }
}

//////////////////////////////////////////////////////////////////////
// compute the rotation-defo tensor in reduced coords
//////////////////////////////////////////////////////////////////////
MATRIX SUBSPACE_TET_MESH::rotationDefoTensor()
{
  MATRIX Ithetaf = _restIthetaf;
  Ithetaf += _UijTildeUij.modeThreeProduct(_q);

  return Ithetaf;
}

//////////////////////////////////////////////////////////////////////
// cache all the values needed for the variable mass matrix
//////////////////////////////////////////////////////////////////////
void SUBSPACE_TET_MESH::cacheMassMatrixVars()
{
  cacheInertiaVars();
  cacheRotationDefoVars();

  // compute center of mass terms
  _SiBar.resizeAndWipe(3, _UBasis.cols()); 
  for (int x = 0; x < _UBasis.rows() / 3; x++)
  {
    Real mass = this->mass(x);
    SUBMATRIX submatrix(_UBasis, 3 * x, 3);
    _SiBar += mass * submatrix;
  }

  _massSummedRestPoses *= 0;
  for (unsigned int x = 0; x < _unconstrainedSize; x++)
  {
    Real mass = this->mass(x);
    _massSummedRestPoses += _restPose[x] * mass;
  }
}

//////////////////////////////////////////////////////////////////////
// compute the center of mass in reduced coords
//////////////////////////////////////////////////////////////////////
VEC3F SUBSPACE_TET_MESH::refreshSitBar()
{
  VEC3F displace(_SiBar * _q);

  _SitBar = _massSummedRestPoses + displace;

  return _SitBar;
}

//////////////////////////////////////////////////////////////////////
// write out an inertia cache
//////////////////////////////////////////////////////////////////////
void SUBSPACE_TET_MESH::writeInertiaCache(string filename)
{
  // if nothing was passed in, try the default
  if (filename.length() == 0)
    filename = _filename + string(".inertia");

  FILE* file;
  file = fopen(filename.c_str(), "wb");

  if (file == NULL)
  {
    cout << " Failed to write to cache " << filename.c_str() << "!" << endl;
    return;
  }

  // cacheInertiaVars
  fwrite((void*)&_restXrestX, sizeof(Real), 1, file);
  fwrite((void*)&_restYrestY, sizeof(Real), 1, file);
  fwrite((void*)&_restZrestZ, sizeof(Real), 1, file);

  fwrite((void*)&_restXrestY, sizeof(Real), 1, file);
  fwrite((void*)&_restXrestZ, sizeof(Real), 1, file);
  fwrite((void*)&_restYrestZ, sizeof(Real), 1, file);

  _restXUX.write(file);
  _restYUY.write(file);
  _restZUZ.write(file);
  
  _restXUY.write(file);
  _restXUZ.write(file);
  _restYUZ.write(file);

  _restYUX.write(file);
  _restZUX.write(file);
  _restZUY.write(file);

  _UXTUX.write(file);
  _UYTUY.write(file);
  _UZTUZ.write(file);

  _UXTUY.write(file);
  _UXTUZ.write(file);
  _UYTUZ.write(file);

  _restYUZrestZUY.write(file);
  _restXUZrestZUX.write(file);
  _restXUYrestYUX.write(file);

  _workspaceL.write(file);

  // cacheRotationDefoVars
  _restIthetaf.write(file);
  _UijTildeUij.write(file);

  // cacheMassMatrixVars
  _SiBar.write(file);
 
  VECTOR massSummed = _massSummedRestPoses.toVector(); 
  massSummed.write(file);

  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// read in an inertia cache
//////////////////////////////////////////////////////////////////////
bool SUBSPACE_TET_MESH::readInertiaCache(string filename)
{
  // if nothing was passed in, try the default
  if (filename.length() == 0)
    filename = _filename + string(".inertia");

  FILE* file;
  file = fopen(filename.c_str(), "rb");

  if (file == NULL)
  {
    //fclose(file);
    return false;
  }

  // cacheInertiaVars
  fread((void*)&_restXrestX, sizeof(Real), 1, file);
  fread((void*)&_restYrestY, sizeof(Real), 1, file);
  fread((void*)&_restZrestZ, sizeof(Real), 1, file);

  fread((void*)&_restXrestY, sizeof(Real), 1, file);
  fread((void*)&_restXrestZ, sizeof(Real), 1, file);
  fread((void*)&_restYrestZ, sizeof(Real), 1, file);

  _restXUX.read(file);
  _restYUY.read(file);
  _restZUZ.read(file);
  
  _restXUY.read(file);
  _restXUZ.read(file);
  _restYUZ.read(file);

  _restYUX.read(file);
  _restZUX.read(file);
  _restZUY.read(file);

  _UXTUX.read(file);
  _UYTUY.read(file);
  _UZTUZ.read(file);

  _UXTUY.read(file);
  _UXTUZ.read(file);
  _UYTUZ.read(file);

  _restYUZrestZUY.read(file);
  _restXUZrestZUX.read(file);
  _restXUYrestYUX.read(file);

  _workspaceL.read(file);

  // cacheRotationDefoVars
  _restIthetaf.read(file);
  _UijTildeUij.read(file);

  // cacheMassMatrixVars
  _SiBar.read(file);
 
  VECTOR massSummed(3);
  massSummed.read(file);
  _massSummedRestPoses = VEC3F(massSummed);

  fclose(file);

  return true;
}

//////////////////////////////////////////////////////////////////////
// recompute the reduced mass matrix
//////////////////////////////////////////////////////////////////////
void SUBSPACE_TET_MESH::recomputeReducedMass()
{
  MATRIX MU = _masses * _UBasis;
  MATRIX UT = _UBasis.transpose();
  _reducedM = UT * MU;
}
