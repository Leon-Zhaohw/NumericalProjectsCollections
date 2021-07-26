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

#ifndef SUBSPACE_TET_MESH_H
#define SUBSPACE_TET_MESH_H

#include <TET_MESH.h>
#include <SUBMATRIX.h>
#include <SETTINGS.h>

#include <ThreadSpecificData.h>

//////////////////////////////////////////////////////////////////////
// Tetrahedron mesh subclass that includes subspace-related
// variables and functions
//////////////////////////////////////////////////////////////////////
class SUBSPACE_TET_MESH : public TET_MESH {

public:
  SUBSPACE_TET_MESH(const char* filename, MATERIAL** materials,
                    int totalMaterials, bool simulateFullspace = false,
                    const char* eigenvectors = NULL, const char* eigenvalues = NULL,
                    Real scale = 1.0);
  virtual ~SUBSPACE_TET_MESH();

public:
  // Helper structure for representing a diagonalized deformation gradient
  struct DeformationDiagonalization {
    DeformationDiagonalization() {}
    DeformationDiagonalization( MATRIX3 &F,
                                MATRIX3 &U, MATRIX3 &Fhat, MATRIX3 &V )
      : _F( F ),
        _U( U ),
        _Fhat( Fhat ),
        _V( V )
    {
    }

    MATRIX3             _F;
    MATRIX3             _U;
    MATRIX3             _Fhat;
    MATRIX3             _V;
  };

  // accessors
  MATRIX& U()           { return _UBasis; };
  VECTOR& eigenvalues() { return _eigenvalues; };
  VECTOR& q()           { return _q; };
  VECTOR& qOld()        { return _qOld; };
  int rank()            { return _q.size(); };
  int totalKeyTets()    { return _totalKeyTets; };
  string filename()     { return _filename; };
  int*& keyTets()       { return _keyTets; };
  VECTOR& keyWeights()  { return _keyWeights; };
  MATRIX& stiffness()   { return _Kreduced; };
  virtual VECTOR& F()   { return _Freduced; };
  VECTOR& reducedInternalForce()    { return _Rreduced; };
  VECTOR& kronrodWeights()  { return _kronrodWeights; };
  MATRIX& leftH()       { return _leftH; };
  VECTOR& forceDensity() { return _forceDensity; };
  MATRIX& surfaceU()	{ return _surfaceU; };
  MATRIX& SiBar()     { return _SiBar; };
  VEC3F& SitBar()     { return _SitBar; };
  VEC3F& massSummedRestPoses() { return _massSummedRestPoses; };

  // Returns the cached deformation state for the given key tet
  DeformationDiagonalization &cachedDeformation( int i )
  {
    return _keyTetDeformations[ i ];
  }

  void updateFullMesh();
  void updateSurfaceMesh();
  VECTOR projectedInternalForce();

  // speculatively update the mesh with an uncommitted q
  void updateFullMesh(VECTOR& qTemporary);

  // get submatrix of _UBasis corresponding to vertices in tet
  MATRIX tetSubBasis(TET* tet);
  void tetSubBasis(TET* tet, MATRIX &subBasis);

  // get submatrix of _UBasis corresponding to vertex
  MATRIX vertexSubBasis(VEC3F* vert);
  MATRIX vertexSubBasis(int vertexID) { return vertexSubBasis(&(_vertices[vertexID])); };

  // Same as the above, but use a given matrix
  bool vertexSubBasis(int vertexID, MATRIX &subBasis,
                      bool transpose = false) const;

  // get the ID corresponding to the tet, needed to index the cubature
  // scheme
  int tetID(TET* tet) { return _tetID[tet]; };

  // read in a cubature scheme
  void readCubature(const char* filename = NULL);

  // see if a basis has already been loaded
  bool basisLoaded() { return _basisLoaded; };

  // compute the external forces
  //void computeExternalForces();

  // dump all vertices to Matlab
  void verticesToMatlab(const char* filename, const char* varName, int column = 0);
  void verticesToBinary(const char* filename, int column);
  
  // Newmark Integrator support functions
  void generateF();

  //< assumes generateF() was already called!
  VECTOR& generateInternalForces();

  //< assumes generateF() was already called!
  VECTOR& generateInternalForcesDebug();

  //< assumes generateF() was already called!
  MATRIX& generateStiffnessMatrix();

  //< assumes generateF() was already called!
  MATRIX& generateStiffnessMatrix(vector<MATRIX3>& Us, 
                                  vector<MATRIX3>& Vs, 
                                  vector<MATRIX>& stiffnesses, map<string,
                                  double>& timingBreakdown);
  
  //< assumes generateF() was already called!
  MATRIX& generateStiffnessMatrixFast( vector<MATRIX3>& Us, 
                                       vector<MATRIX3>& Vs, 
                                       vector<MATRIX>& stiffnesses,
                                       map<string, double>& timingBreakdown);

  // This is similar to the above function, but instead of rebuilding
  // the stiffness matrix from scratch, updates the existing stiffness
  // matrix according to which key tets are dirty.
  // Assumes generateF() was already called!
  MATRIX& updateStiffnessMatrixFast( vector<MATRIX3> &Us,
                                     vector<MATRIX3> &Vs,
                                     vector<MATRIX> &stiffnesses,
                                     vector<bool> &keyTetsDirty,
                                     map<string, double>& timingBreakdown);
  
  //< assumes cacheDiagonalizations() was already called!
  VECTOR& generateInternalForces(vector<MATRIX3>& Us, 
                                 vector<MATRIX3>& Fhats,
                                 vector<MATRIX3>& Vs);
  
  MATRIX& generateKryslStiffnessMatrix();
  
  //< assumes generateF() was already called!
  MATRIX& generateKryslStiffnessMatrix(vector<MATRIX3>& Us, 
                                       vector<MATRIX3>& Vs, 
                                       vector<MATRIX>& stiffnesses,
                                       map<string, double>& timingBreakdown);
  VECTOR& generateKryslInternalForces();
  
  //< assumes cacheDiagonalizations() was already called!
  VECTOR& generateKryslInternalForces(vector<MATRIX3>& Us, 
                                 vector<MATRIX3>& Fhats,
                                 vector<MATRIX3>& Vs);
  MATRIX& reducedMass() { return _reducedM; };

  // compute error of the cubature using Gauss-Kronrod
  Real cubatureError();

  // check if any of the key tets are inverted
  bool keyInverted();

  // update mesh after the basis has changed
  void updateBasis(MATRIX& newBasis);

  // update cubature quantities after key tets changed
  void updateCubature();

  // build the subset of _UBasis that corresponds to the surface
  void buildSurfaceU();

  // blow away the entire basis
  void resetBasis();
  
  // draw a basis vector
  void drawBasis(int whichBasis, Real amplitude);

  // read and write state
  void readState(string filename);
  void writeState(string filename);

  // draw key tets to RenderMan
  void drawKeyTetsToRenderMan();

  // reset the mass matrix to a specific mass
  void resetMasses(Real mass);

  // reset state to zero
  virtual void reset();

  // Initializes workspaces for fast reduced stiffness assembly
  void initializeAssemblyWorkspace( bool partialUpdate = false );

  // reduced order inertia tensor
  virtual const MATRIX& refreshInertiaTensor();
  virtual const MATRIX& refreshInertiaTensorDt(VECTOR& velocity);

  // cache all the values needed for the variable mass matrix
  void cacheMassMatrixVars();

  // compute the rotation-defo tensor in reduced coords
  MATRIX rotationDefoTensor();

  // refesh the center of mass in reduced coords
  VEC3F refreshSitBar();

  // read/write out the inertia vars
  void writeInertiaCache(string filename = string(""));
  bool readInertiaCache(string filename = string(""));

  // recompute the reduced mass matrix
  void recomputeReducedMass();

protected:
  // filename everything is loaded from
  string _filename;
  
  // eigensystem from linear modal analysis
  // I would have named the matrix just "_U" but gcc doesn't like that
  // because it uses _U in STL
  VECTOR _eigenvalues;
  MATRIX _UBasis;

  /*
  // current state
  // This is called 'q' only for notational clarity. It points
  // to the _x vector from the superclass so as not to waste memory.
  // VECTOR _q;
  VECTOR& _q;
  */
  // the online solve needs _q and _x to be separate, so _q is going to
  // be tracked separately from now on
  VECTOR _q;
  VECTOR _qOld;

  // maps tet address to its index in _tets
  map<TET*, int> _tetID;

  // the key tets - this is kept in a flat array instead of a vector
  // because we need direct access for performance reasons
  int* _keyTets;
  VECTOR _keyWeights;
  VECTOR _kronrodWeights;
  int _totalKeyTets;
  VECTOR _keyRestVolumes; // < Note that these are (1.0 / volume) for performance

  // subset of _UBasis rows that correspond only to the surface vertices
  MATRIX _surfaceU;

  // reduced mass matrix;
  MATRIX _reducedM;

  // force densites from MATERIAL
  VECTOR _forceDensity;

  // block diagonal stiffness densities from MATERIAL
  // Don't use a vector here because it's not thread safe for OpenMP --
  // just use an array
  MATRIX** _stiffnessDensity;

  // intermediate value for stiffness density, ie "Stiffness * _rightH"
  MATRIX _SH;

  // Precomputed matrices for Newmark (see paper)
  MATRIX _E;
  MATRIX _leftH;
  MATRIX _rightH;
  MATRIX _leftHGaussKronrod;

  // just checking th eigenvalue size has failed before, so track an
  // explicit bool
  bool _basisLoaded;

  // populate the above matrices -- automatically called when
  // a cubature scheme is loaded
  void generateE();
  void generateH();

  // track reduced F, R and K separately because the online
  // solve needs the unreduced ones around too
  VECTOR _Freduced;
  VECTOR _Rreduced;
  MATRIX _Kreduced;
  VECTOR _RreducedGaussKronrod;

  // Data for stiffness matrix assembly
  // 
  // All pFpU matrices, with each column repacked in to a 3x3 matrix,
  // and then stacked in to a block vector
  Real *_pFpUrepacked;

  // Sub-bases for all key tets
  Real *_tetSubBases;

  // Transposes of all key tet sub-bases
  Real *_tetSubBasesTransposed;

  // B matrix ( eg [ b_0 ... b_3 ] ) stacked for all key tets
  Real *_tetBmatrices;

  // Workspaces for internal force computation
  ThreadSpecificData<VECTOR> _flatPKworkspace;

  //////////////////////////////////////////////////////////////////////
  // Data used to when performing partial updates to the
  // reduced stiffness matrix.
  //////////////////////////////////////////////////////////////////////

  // The current r x r contribution to the stiffness matrix sum for
  // each key tet.  We need this so that we can subtract it from the
  // running stiffness matrix total without having to recompute it
  // from the previous state (which may have occurred at some distant
  // point in time).
  std::vector<MATRIX> _keyTetStiffnesses;

  // Also maintain a list of deformation gradients for each key
  // tet, at the last state for which they were updated.
  std::vector<DeformationDiagonalization> _keyTetDeformations;



  // cached variables for reduced inertia tensor computation
  Real _restXrestX;
  Real _restYrestY;
  Real _restZrestZ;

  Real _restXrestY;
  Real _restXrestZ;
  Real _restYrestZ;

  VECTOR _restXUX;
  VECTOR _restYUY;
  VECTOR _restZUZ;
  
  VECTOR _restXUY;
  VECTOR _restXUZ;
  VECTOR _restYUZ;

  VECTOR _restYUX;
  VECTOR _restZUX;
  VECTOR _restZUY;

  VECTOR _restYUZrestZUY;
  VECTOR _restXUZrestZUX;
  VECTOR _restXUYrestYUX;

  MATRIX _UXTUX;
  MATRIX _UYTUY;
  MATRIX _UZTUZ;

  MATRIX _UXTUY;
  MATRIX _UXTUZ;
  MATRIX _UYTUZ;

  // cached variables for reduced rotation-defo computation
  MATRIX _restIthetaf;
  TENSOR3 _UijTildeUij;

  // cached center of mass variables
  MATRIX _SiBar;
  VEC3F _massSummedRestPoses;
  VEC3F _SitBar;

  // cached workspace
  VECTOR _workspaceL;

  // precompute all the mass matrix quantities -- must be public so 
  // PARTITIONED_SUBSPACE_TET_MESH can call it after altering all the masses
  void cacheInertiaVars();

  // precompute all the rotation-defo tensor quantities -- must be public so 
  // PARTITIONED_SUBSPACE_TET_MESH can call it after altering all the masses
  void cacheRotationDefoVars();
};

#endif
