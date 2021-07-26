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
//------------------------------------------------------------------------------
// GL interface elements are from:
//------------------------------------------------------------------------------
// GLVU : Copyright 1997 - 2002 
//        The University of North Carolina at Chapel Hill
//------------------------------------------------------------------------------
// Permission to use, copy, modify, distribute and sell this software and its 
// documentation for any purpose is hereby granted without fee, provided that 
// the above copyright notice appear in all copies and that both that copyright 
// notice and this permission notice appear in supporting documentation. 
// Binaries may be compiled with this software without any royalties or 
// restrictions. 
//
// The University of North Carolina at Chapel Hill makes no representations 
// about the suitability of this software for any purpose. It is provided 
// "as is" without express or implied warranty.

#include <iostream>
#include <TET_MESH.h>
#include <PARTITIONED_SKINNED_SUBSPACE_TET_MESH.h>
#include <PARTITIONED_TET_MESH.h>
#include <PARTITIONED_CUBATURE_GENERATOR.h>
#include <SPARSE_MATRIX.h>
#include <STVK.h>
#include <MOONEY_RIVLIN.h>
#include <ARRUDA_BOYCE.h>
#include <NEO_HOOKEAN.h>
#include <INVERTIBLE.h>
#include <SIMPLE_PARSER.h>
#include <FULLSPACE_INTEGRATOR.h>
#include <GNUPLOT.h>
#include <SKELETON.h>

using namespace std;

PARTITIONED_SKINNED_SUBSPACE_TET_MESH* partitionedMesh = NULL;
string triangleMeshPath;
string triangleMeshName;
string outputPath;
string posePath("");
string tetMeshName;
Real discardGramSchmidt = 1e-8;
Real discardPCA = 1e-7;

// training parameters
int maxKeyTets;
int candidatesPerTry;
double errorTolerance;
int randSeed;
int rank;
int subrank;
int bccRes = 22;
int which = 0;

// how to sparsely sample the data?
//int snapshotStride = 2;
//int snapshotStride = 4;
int snapshotStride = 1;

FULLSPACE_INTEGRATOR* integrator = NULL;
SKELETON* skeleton = NULL;

PARTITIONED_CUBATURE_GENERATOR* generator = NULL;
int snapshots;
int startingSnapshot;
bool useInvertible = false;
vector<MATRIX3>* snapshotRotations = NULL;

vector<string> dataPaths;
vector<string> mocapNames;

////////////////////////////////////////////////////////////////
// Print integer to a zero-padded string
//////////////////////////////////////////////////////////////////
static std::string itoPaddedString(int frame)
{
  char buffer[256];
  sprintf(buffer, "%i", frame);

  std::string number = std::string(buffer);
  if (frame < 10) number = std::string("0") + number;
  if (frame < 100) number = std::string("0") + number;
  if (frame < 1000) number = std::string("0") + number;

  return number;
}

//////////////////////////////////////////////////////////////////////////////
// Read a material file and return a pointer to it
//////////////////////////////////////////////////////////////////////////////
MATERIAL* readMaterial(SIMPLE_PARSER& parser)
{
  // set material
  MATERIAL* material = NULL;
  string materialType;
  materialType = parser.getString("material type", materialType);

  if (materialType.compare("stvk") == 0)
  {
    double lambda = 10.0;
    double mu = 50.0;
    lambda = parser.getFloat("stvk lambda", lambda);
    mu = parser.getFloat("stvk mu", mu);
    material = new STVK(lambda, mu);
    cout << "==================================================" << endl;
    cout << " Material is St-VK, lambda = " << lambda << " mu = " << mu << endl;
  }
  else if (materialType.compare("mooney-rivlin") == 0)
  {
    double mu01 = 100.0;
    double mu10 = 500.0;
    double k = 100000.0;
    mu01 = parser.getFloat("mooney-rivlin mu01", mu01);
    mu10 = parser.getFloat("mooney-rivlin mu10", mu10);
    k = parser.getFloat("mooney-rivlin k", k);
    material = new MOONEY_RIVLIN(mu01, mu10, k);
    cout << "==================================================" << endl;
    cout << " Material is Mooney-Rivlin, mu01 = " << mu01 << " mu10 = " << mu10 << " k = " << k << endl;
  }
  else if (materialType.compare("arruda-boyce") == 0)
  {
    double nkTheta = 5000.0;
    double N  = 5.0;
    double k = 1000.0;
    nkTheta = parser.getFloat("arruda-boyce nktheta", nkTheta);
    N = parser.getFloat("arruda-boyce n", N);
    k = parser.getFloat("arruda-boyce k", k);
    material = new ARRUDA_BOYCE(nkTheta, N, k);
    cout << "==================================================" << endl;
    cout << " Material is Arruda-Boyce, nkTheta = " << nkTheta << " N = " << N << " k = " << k << endl;
  }
  else if (materialType.compare("neo-hookean") == 0)
  {
    double mu = 50.0;
    double lambda = 10.0;
    mu = parser.getFloat("neo-hookean mu", mu);
    lambda = parser.getFloat("neo-hookean lambda", lambda);
    material = new NEO_HOOKEAN(mu, lambda);
    cout << "==================================================" << endl;
    cout << " Material is Neo-Hookean, mu = " << mu << " lambda = " << lambda << endl;
  }
  else
  {
    cout << " *** Material type undefined! *** " << endl;
    exit(1);
  }

  return material;
}

//////////////////////////////////////////////////////////////////////////////
// undo the rotation in a snapshot for a given partition
//////////////////////////////////////////////////////////////////////////////
void subtractRotation(int partition, VECTOR& snapshot)
{
  UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)partitionedMesh->mesh(partition);
  MATRIX3 rotation = mesh->rotationQuaternion().toExplicitMatrix3x3();
  MATRIX3 transpose = rotation.transpose();

  for (int x = 0; x < snapshot.size() / 3; x++)
  {
    // copy into VEC3s
    VEC3F position3;
    position3[0] = snapshot[3 * x];
    position3[1] = snapshot[3 * x + 1];
    position3[2] = snapshot[3 * x + 2];

    // perform the undo
    position3 = transpose * position3;

    // copy back into the bigger vector
    snapshot[3 * x] = position3[0];
    snapshot[3 * x + 1] = position3[1];
    snapshot[3 * x + 2] = position3[2];
  }
}

//////////////////////////////////////////////////////////////////////////////
// extract the portion of this vector that is in the submesh
//////////////////////////////////////////////////////////////////////////////
VECTOR getSubvector(const int partition, const VECTOR& fullVector)
{
  SUBSPACE_TET_MESH* subMesh= (SUBSPACE_TET_MESH*)partitionedMesh->mesh(partition);
  TET_MESH* fullMesh = partitionedMesh->originalMesh();
 
  assert(fullVector.size() == fullMesh->x().size());

  // set a matrix to the full vector so we can use SUBMATRIX
  MATRIX snapshot(fullMesh->x().size(), 1);
  snapshot.setColumn(fullVector, 0);

  int partitionNodes = subMesh->unconstrainedNodes();
  MATRIX subSnapshotMatrix(partitionNodes * 3, 1);

  // for each vertex in the partition
  for (int y = 0; y < partitionNodes; y++)
  {
    // get the vertex index in the original mesh
    int originalID = partitionedMesh->originalID(partition, y);

    // get the submatrix corresponding to that vertex
    SUBMATRIX vertexBasis(snapshot, originalID * 3, 3);

    // copy the vertex basis into the subbasis
    vertexBasis.copiesInto(subSnapshotMatrix, 3 * y);
  }

  return subSnapshotMatrix.getColumn(0);
}

//////////////////////////////////////////////////////////////////////////////
// compute the mesh force density
//////////////////////////////////////////////////////////////////////////////
VECTOR computeForceDensity(TET_MESH* tetMesh)
{
  vector<TET>& tets = tetMesh->tets();
  VECTOR densities(tets.size() * 12);

  for (unsigned int x = 0; x < tets.size(); x++)
  {
    TET& tet = tets[x];
    MATERIAL* material = tetMesh->materials()[tet.materialIndex()];
    MATRIX3 F3 = tet.F();
    VECTOR F = TET::flattenF(F3);

    VECTOR forceDensity(9);
    Real* data = forceDensity.data();
    material->forceDensity(F.data(), data);
    MATRIX pFpu(9,12);

    // convert density to a force
    material->computePFPu(tet, pFpu);
    VECTOR densityVector = forceDensity * pFpu;

    // store the force
    data = densities.data();
    for (int y = 0; y < 12; y++)
      data[x * 12 + y] = densityVector[y];
  }

  return densities;
}

//////////////////////////////////////////////////////////////////////////////
// Generate training samples the LMA way, but format them like a
// skinned PCA sample -- store constrained displacements and force densities
// as well
//////////////////////////////////////////////////////////////////////////////
void generateLMATrainingSamples(int partition, int totalSamples, Real magnitude, int sampleModes)
{
  SUBSPACE_TET_MESH* mesh = (SUBSPACE_TET_MESH*)(partitionedMesh->mesh(partition));
  VECTOR& eigenvalues = mesh->eigenvalues();
  VECTOR q(eigenvalues.size());
  static MERSENNETWISTER trainingTwister(123456);
  vector<VECTOR>& trainingForces = generator->trainingForcesPartitioned(partition);
  vector<VECTOR>& trainingQs= generator->trainingQsPartitioned(partition);
  vector<VECTOR>& trainingColumns = generator->trainingColumnsPartitioned(partition);
  vector<VECTOR>& trainingConstraints = generator->trainingConstraintsPartitioned(partition);

  // the user wants the first mode to have the given magnitude
  magnitude = magnitude * sqrt(eigenvalues(0));

  int attempts = 0;
  int successes = 0;

  // track if a Nan has been seen before so that the error message
  // doesn't get output too many times
  bool firstNan = true;

  cout << " Generating LMA training set for partition " << partition << " ...";
  for(int i = 0; i < totalSamples; i++)
  {
    int uninvertedTries = 0;
    bool tetInverted = true;

    // try really hard to generate an uninverted sample
    while (tetInverted && uninvertedTries < 100)
    {
      // create random sample
      //for (int x = 0; x < q.size(); x++)
      for (int x = 0; x < sampleModes; x++)
      {
        // scale by inverse freq squared (ie. divide by eig val)
        // 4.0 because.. gauss is usually between -4 and 4
        q(x) = magnitude / sqrt(eigenvalues(x)) / 4.0 * trainingTwister.randNorm();

        // if the eigenvalue is small enough to produce a nan
        // (this has happened before) try to do something graceful
        if (isnan(q(x))) 
        {
          q(x) = 0.0;
          if (firstNan)
          {
            firstNan = false;
            cout << __FILE__ << " " << __LINE__ << " : Very small eigenvalue detected!" << endl
                 << " You should reduce your rank if it is very large "<< endl
                 << " or increase the resolution of your mesh if it is very coarse." << endl;
          }
        }
      }

      mesh->q() = q;
      mesh->updateFullMesh();

      // check for inversion
      tetInverted = mesh->inverted();

      // DEBUG
      tetInverted = false;

      uninvertedTries++;
      attempts++;
    }

    // throw a flag if the tet is still inverted, but forge ahead
    if (tetInverted)
    {
      cout << __FILE__ << " " << __LINE__ << " : Inverted Tet Mesh! "
           << " You should try again with a smaller magnitude" << endl;
    }
    else
      successes++;

    // store q and the force it produces
    mesh->TET_MESH::F().resizeAndWipe(9 * mesh->totalTets());
    mesh->TET_MESH::generateF();  
    VECTOR subInternalForces = mesh->TET_MESH::generateInternalForces();
    VECTOR forceSample = mesh->U() ^ subInternalForces;

    //VECTOR forceSample = mesh->projectedInternalForce();
    trainingForces.push_back(forceSample);
    trainingQs.push_back(q);

    // store constraints and force densities too
    VECTOR forceDensity = computeForceDensity(mesh);
    trainingColumns.push_back(forceDensity);
    VECTOR constraints = mesh->constraintDisplacements();
    trainingConstraints.push_back(constraints);

    // output status every 10%
    if (i % (int)(totalSamples / 10) == 0)
    {
      cout << 100 * ((Real)i / totalSamples) << "% ";
      flush(cout);
    }
  }
  cout << " done. " << endl;
  cout << " " << successes << " of " << attempts << " attempts successful (" 
        << 100.0f * (Real)successes / attempts << "%)" << endl;
}

//////////////////////////////////////////////////////////////////////////////
// Use the new PCA basis to build the training samples in a streaming way
//////////////////////////////////////////////////////////////////////////////
void generatePCATrainingSamples(int partition, map<int, bool>& usedPositions, map<int, bool>& usedVelocities)
{
  SUBSPACE_TET_MESH* subMesh= (SUBSPACE_TET_MESH*)partitionedMesh->mesh(partition);
  vector<VECTOR>& reducedTrainingForces = generator->trainingForcesPartitioned(partition);
  vector<VECTOR>& reducedTrainingPoses = generator->trainingQsPartitioned(partition);
  vector<VECTOR>& reducedTrainingColumns = generator->trainingColumnsPartitioned(partition);
  vector<VECTOR>& reducedTrainingConstraints = generator->trainingConstraintsPartitioned(partition);

  // read in cached, if possible
  string meshFilename = subMesh->filename();
  string trainingFilename = meshFilename + string(".training");
  //FILE* file = fopen(trainingFilename.c_str(), "rb");
  //fclose(file);

  VECTOR restVector = subMesh->restVector();

  TET_MESH* fullMesh = partitionedMesh->originalMesh();
  MATRIX& U = subMesh->U();

  int constrainedNodes = subMesh->totalNodes() - subMesh->unconstrainedNodes();
  bool unconstrained = (constrainedNodes == 0);

  reducedTrainingForces.clear();
  reducedTrainingPoses.clear();
  reducedTrainingColumns.clear();
  reducedTrainingConstraints.clear();

  cout << " Computing PCA training samples ... "; flush(cout);
  for (unsigned int i = 0; i < dataPaths.size(); i++)
  {
    string dataPath = dataPaths[i];
    cout << " Loading data from path " << dataPath.c_str() << endl;

    for (int x = startingSnapshot; x < snapshots; x++)
    {
      if (x % (int)(snapshots / 10) == 0)
      {
        cout << 100 * ((Real)x / snapshots) << "% ";
        flush(cout);
      }

      // Only use the snapshots found by Gram-Schmidt?
      //if (usedPositions.find(x) == usedPositions.end() &&
      //    usedVelocities.find(x) == usedVelocities.end()) continue;
      if (usedPositions[x] == false && usedVelocities[x] == false) continue;

      char buffer[256];
      sprintf(buffer, "%04i", x);
      string skeletonFile = posePath + string("ode.motion.") + string(buffer) + string(".skeleton");
      skeleton->loadOdeFrame(skeletonFile.c_str());
      skeleton->updateOdeSkinning(false);

      // FULLSPACE_INTEGRATOR version
      string integratorFile = dataPath + string("full.integrator.");
      integratorFile += itoPaddedString(x);
      integratorFile += string(".state");
      integrator->readState(integratorFile);

      // load simulation results into mesh
      //partitionedMesh->loadOdeFrame(posePath, x, integrator->position());
      partitionedMesh->loadOdeFrame(dataPath, x, integrator->position());

      // retrieve the position
      VECTOR subPosition = subMesh->x();

      // make sure that the force sample is for a pose that the basis can actually
      // resolve
      if (usedPositions.find(x) != usedPositions.end())
      {
        // CAREFUL TO USE THE PROJECTED POSITION HERE, NOT THE FULL POSITION
        // otherwise the force and force density computed will be out of basis
        VECTOR projectedPosition= U ^ subPosition;
        subMesh->x() = U * projectedPosition;
        subMesh->TET_MESH::updateFullMesh();

        subMesh->TET_MESH::F().resizeAndWipe(9 * subMesh->totalTets());
        subMesh->TET_MESH::generateF();  
        VECTOR subInternalForces = subMesh->TET_MESH::generateInternalForces();
        VECTOR projectedForce  = U ^ subInternalForces;
        reducedTrainingForces.push_back(projectedForce);
        reducedTrainingPoses.push_back(projectedPosition);

        VECTOR forceDensity = computeForceDensity(subMesh);
        reducedTrainingColumns.push_back(forceDensity);

        // get the contraints
        VECTOR constraints = subMesh->constraintDisplacements();
        reducedTrainingConstraints.push_back(constraints);
      }
    }
  }
  cout << "done." << endl;
  cout << "Using " << reducedTrainingForces.size() << " training forces " << endl;

  /*
  // write to a cache
  file = fopen(trainingFilename.c_str(), "wb");
  int size = reducedTrainingForces.size();
  fwrite((void*)&size, sizeof(int), 1, file); 

  for (int x = 0; x < size; x++)
  {
    reducedTrainingForces[x].write(file);
    reducedTrainingPoses[x].write(file);
    reducedTrainingColumns[x].write(file);
    reducedTrainingConstraints[x].write(file);
  }

  fclose(file);
  */
}

//////////////////////////////////////////////////////////////////////////////
// Try to add a snapshot to the basis
//////////////////////////////////////////////////////////////////////////////
bool addSnapshot(vector<VECTOR>& Qs, vector<Real>& Rs, VECTOR& snapshot)
{
  Qs.push_back(snapshot);
  Rs.push_back(1.0);
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// Perform gram schmidt on the snapshots
//////////////////////////////////////////////////////////////////////////////
void computeSkinnedGramSchmidt(int partition,
                               Real discardGramSchmidt,
                               MATRIX& gramSchmidt,
                               VECTOR& Rs,
                               map<int,bool>& usedPositions,
                               map<int,bool>& usedVelocities)
{
  SUBSPACE_TET_MESH* subMesh= (SUBSPACE_TET_MESH*)partitionedMesh->mesh(partition);
  int constrainedNodes = subMesh->totalNodes() - subMesh->unconstrainedNodes();

  // get a vector of all the rest positions
  VECTOR restVector = subMesh->restVector();

  bool unconstrained = (constrainedNodes == 0);
  cout << " Partition " << partition << " is";
  if (unconstrained)
    cout << " unconstrained " << endl;
  else
    cout << " constrained " << endl;
  cout << " Constrained nodes: " << constrainedNodes << endl;

  // stream in the data, discard the ones that are too small according
  // to Gram-Schmidt
  vector<VECTOR> Qs;
  vector<Real> rawRs;
  cout << " Computing Gram-Schmidt ... "; flush(cout);
  usedPositions.clear();
  usedVelocities.clear();

  int velocityRank = 0;
  int positionRank = 0;

  // read in multiple data paths
  for (unsigned int i = 0; i < dataPaths.size(); i++)
  {
    string dataPath = dataPaths[i];
    cout << " Reading in snapshots from " << dataPath << endl;

    for (int x = startingSnapshot; x < snapshots; x++)
    {
      if (x % (int)(snapshots / 20) == 0)
      {
        cout << 100 * ((Real)x / snapshots) << "% ";
        cout << "(Rank: " << Qs.size() << " position: " << positionRank << " velocity: " << velocityRank << ") ";
        flush(cout);
      }

      // sparsely sample the input
      if (x % snapshotStride != 0)
      {
        usedPositions[x] = false;
        continue;
      }

      char buffer[256];
      sprintf(buffer, "%04i", x);
      string skeletonFile = posePath + string("ode.motion.") + string(buffer) + string(".skeleton");
      skeleton->loadOdeFrame(skeletonFile.c_str());
      skeleton->updateOdeSkinning(false);

      VECTOR velocitySnapshot;
      VECTOR positionSnapshot;

      // FULLSPACE_INTEGRATOR velocity-level version
      string integratorFile = dataPath + string("full.integrator.");
      integratorFile += itoPaddedString(x);
      integratorFile += string(".state");
      integrator->readState(integratorFile);

      // load simulation results into mesh
      //partitionedMesh->loadOdeFrame(posePath, x, integrator->position());
      partitionedMesh->loadOdeFrame(dataPath, x, integrator->position());

      // retrieve the position
      TET_MESH* mesh = partitionedMesh->mesh(partition);
      mesh->recoverX();
      positionSnapshot = mesh->x();

      // try adding position -- if it's used, use it as a training sample too
      bool usedPosition = addSnapshot(Qs, rawRs, positionSnapshot);
      if (usedPosition) 
      {
        positionRank++;
        usedPositions[x] = true;
      }
    }
  }
  cout << "done." << endl;
  cout << " Kept " << Qs.size() << " of " << 2 * snapshots << " snapshots" << endl;
  
  float sizeU = Qs.size() * Qs[0].size() * 8.0 / pow(2.0, 20.0);
  cout << " Snapshots consume " << sizeU << " MB " << endl;

  // store the Rs so that we can reweight PCA later as well
  cout << " Building final Gram-Schmidt matrix ... "; flush(cout);
  gramSchmidt = MATRIX(Qs);
  Rs = VECTOR(rawRs);
  cout << " done." << endl;
}

//////////////////////////////////////////////////////////////////////////////
// Build the Gram-Schmidt result
//////////////////////////////////////////////////////////////////////////////
void buildSkinnedGramSchmidtBasis(int partition, int maxRank, map<int, bool>& usedPositions, map<int, bool>& usedVelocities)
{
  SUBSPACE_TET_MESH* subMesh = (SUBSPACE_TET_MESH*)partitionedMesh->mesh(partition);

  cout << endl;
  cout << " Factoring snapshots ... ";

  // make sure TET_MESH put aside room for F and internal forces
  subMesh->TET_MESH::internalForce().resizeAndWipe(subMesh->TET_MESH::rank());
  subMesh->TET_MESH::F().resizeAndWipe(subMesh->TET_MESH::rank() * 9);

  // allocate in advance
  MATRIX gramSchmidt;

  string meshFilename = subMesh->filename();
  string gramSchmidtName = meshFilename + string(".gramschmidt.matrix");
  string rName        = meshFilename + string(".gramschmidt.vector");
  string trainingName = meshFilename + string(".reduced.trainingdata");

  // always do a stomp?
  string rmGram = string("rm ") + gramSchmidtName;
  string rmRs = string("rm ") + gramSchmidtName;
  string rmTraining = string("rm ") + trainingName;
  cout << " Stomping previous Gram-Schmidt ... "; flush(cout);
  system(rmGram.c_str());
  system(rmRs.c_str());
  system(rmTraining.c_str());

  // compute gram schmidt from scratch
  VECTOR Rs;
  computeSkinnedGramSchmidt(partition, discardGramSchmidt, gramSchmidt, Rs, usedPositions, usedVelocities);
  //gramSchmidt.write(gramSchmidtName.c_str());
  //Rs.write(rName.c_str());
  
  subMesh->updateBasis(gramSchmidt);
  MATRIX& U = subMesh->U();

  // cache the final principal values
  VECTOR& eigenvalues = subMesh->eigenvalues();
  eigenvalues.resizeAndWipe(Rs.size());
  for (int x = 0; x < Rs.size(); x++)
    eigenvalues[x] = Rs[x];

  //string meshFilename = subMesh->filename();
  string valuename = meshFilename + string(".eigenvalues.vector");
  string vectorname = meshFilename + string(".eigenvectors.matrix");
  cout << " Writing Gram-Schmidt Qs to file " << valuename.c_str() << endl;
  cout << " Writing Gram-Schmidt Rs to file " << vectorname.c_str() << endl;
  eigenvalues.write(valuename.c_str());
  U.write(vectorname.c_str());

  // stomp any previous surfaceU
  string surfacename = meshFilename + string(".surfaceU.matrix");
  cout << " Stomping any old surface U: " << surfacename.c_str() << endl;
  string rm = string("rm ") + surfacename;
  system(rm.c_str());

  //cout << " Gram-Schmidt Rs: " << endl;
  //cout << Rs << endl;
}

//////////////////////////////////////////////////////////////////////////////
// project and then unproject force snapshots to see what the error looks like
//////////////////////////////////////////////////////////////////////////////
void projectSnapshotForces(PARTITIONED_SKINNED_SUBSPACE_TET_MESH* partitionedMesh)
{
  int totalPartitions = partitionedMesh->partitions();
  for (unsigned int i = 0; i < dataPaths.size(); i++)
  {
    string dataPath = dataPaths[i];

    string mocapFile = mocapNames[i];
    cout << " Loading mocap file " << mocapFile.c_str() << endl;
    skeleton->readMocap(mocapFile);

    for (int x = startingSnapshot; x < snapshots; x++)
    {
      cout << " force snapshot " << x << endl;
      // load up the snapshot
      string integratorFile = dataPath + string("full.integrator.");
      integratorFile += itoPaddedString(x);
      integratorFile += string(".state");
      integrator->readState(integratorFile);

      // load simulation results into mesh
      partitionedMesh->loadSimulationFrame(dataPath, x, integrator->position());

      // look at the projection
      for (int y = 0; y < totalPartitions; y++)
      {
        SUBSPACE_TET_MESH* mesh = (SUBSPACE_TET_MESH*)partitionedMesh->mesh(y);
        mesh->recoverX();
        VECTOR subPosition = mesh->x();

        // resize F for computation
        mesh->TET_MESH::F().resizeAndWipe(9 * mesh->totalTets());

        // get the full internal force
        mesh->TET_MESH::generateF();
        VECTOR fullInternals = mesh->TET_MESH::generateInternalForces();

        // resize reduced F for computation
        mesh->SUBSPACE_TET_MESH::F().resizeAndWipe(9 * mesh->totalKeyTets());

        // get the reduced internal force
        mesh->SUBSPACE_TET_MESH::generateF();
        VECTOR reducedInternals = mesh->SUBSPACE_TET_MESH::generateInternalForces();

        // unproject it
        VECTOR unprojected = mesh->U() * reducedInternals;

        // look at the diff
        VECTOR diff = fullInternals - unprojected;
        Real relative = diff.norm2() / fullInternals.norm2();
        cout << "abs diff: " << diff.norm2() << " relative diff: " << relative << endl;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
// project and then unproject snapshots to see what the error looks like
//////////////////////////////////////////////////////////////////////////////
void projectSnapshotPositions(PARTITIONED_SKINNED_SUBSPACE_TET_MESH* partitionedMesh)
{
  int totalPartitions = partitionedMesh->partitions();
  for (unsigned int i = 0; i < dataPaths.size(); i++)
  {
    string dataPath = dataPaths[i];

    string mocapFile = mocapNames[i];
    cout << " Loading mocap file " << mocapFile.c_str() << endl;
    skeleton->readMocap(mocapFile);

    for (int x = startingSnapshot; x < snapshots; x++)
    {
      //cout << " Snapshot " << x << endl;
      // load up the snapshot
      string integratorFile = dataPath + string("full.integrator.");
      integratorFile += itoPaddedString(x);
      integratorFile += string(".state");
      integrator->readState(integratorFile);

      // load simulation results into mesh
      partitionedMesh->loadSimulationFrame(dataPath, x, integrator->position());

      // look at the projection
      //for (int y = 0; y < totalPartitions; y++)
      for (int y = 0; y < 1; y++)
      {
        SUBSPACE_TET_MESH* mesh = (SUBSPACE_TET_MESH*)partitionedMesh->mesh(y);
        mesh->recoverX();
        VECTOR subPosition = mesh->x();

        // project / unproject
        VECTOR projected = mesh->U() ^ subPosition;
        VECTOR unprojected = mesh->U() * projected;
       
        //Real relative = (subPosition.norm2() - unprojected.norm2()) / subPosition.norm2();
        //cout << " original: " << subPosition.norm2() <<  " projected: " << unprojected.norm2() << " relative: " << relative << endl;
        VECTOR diff = subPosition - unprojected;

        VECTOR::printVertical = false;
        cout << " Snapshot " << x << ": " << endl;
        cout << " q: " << projected << endl;
        cout << " position diff norm: " << diff.norm2() << endl;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
// project and then unproject gravity to see what the error looks like
//////////////////////////////////////////////////////////////////////////////
void projectGravity(PARTITIONED_SUBSPACE_TET_MESH* partitionedMesh)
{
  for (int x = 0; x < partitionedMesh->partitions(); x++)
  {
    cout << " Checking partition " << x << endl;
    SUBSPACE_TET_MESH* subMesh = (SUBSPACE_TET_MESH*)partitionedMesh->mesh(x);
    MATRIX& U = subMesh->U();

    VECTOR gravity(U.rows());
    for (int y = 0; y < U.rows() / 3; y++)
      gravity[3 * y + 1] = 1.0;

    VECTOR projected = U ^ gravity;
    VECTOR unprojected = U * projected;

    //cout << " gravity: " << gravity << endl;
    //cout << " after projection: " << unprojected << endl;
    cout << " gravity before: " << gravity.norm2() << endl;
    cout << " gravity after: " << unprojected.norm2() << endl;
    
  }
}

//////////////////////////////////////////////////////////////////////////////
// load up the lma and modal derivatives matrices
//////////////////////////////////////////////////////////////////////////////
void loadLMA(MATRIX& lmaMatrix, VECTOR& lmaVector, MATRIX& modalMatrix, VECTOR& modalVector)
{
  // load in the LMA modes
  string lmaMatrixName = tetMeshName + string(".eigenvectors.matrix");
  string lmaVectorName = tetMeshName + string(".eigenvalues.vector");

  if (lmaMatrix.read(lmaMatrixName.c_str()) == false)
  {
    cout << " Didn't find any LMA matrix! Run ArpackLMA! " << endl;
    exit(0);
  }
  if (lmaVector.read(lmaVectorName.c_str()) == false)
  {
    cout << " Didn't find any LMA vector! Run ArpackLMA! " << endl;
    exit(0);
  }

  // load in the modal derivatives
  string modalMatrixName = tetMeshName + string(".modal.derivatives.matrix");
  string modalVectorName = tetMeshName + string(".modal.eigenvalues.vector");

  if (modalMatrix.read(modalMatrixName.c_str()) == false)
  {
    cout << " Didn't find any modal derivatives matrix! Run ModalDerivatives! " << endl;
    exit(0);
  }
  if (modalVector.read(modalVectorName.c_str()) == false)
  {
    cout << " Didn't find any modal derivatives vector! Run ModalDerivatives! " << endl;
    exit(0);
  }

  cout << " Loaded LMA matrix of dims " << lmaMatrix.dims() << endl;
  cout << " Vector has size: " << lmaVector.size() << endl;
  cout << " Loaded modal derivative matrix of dims " << modalMatrix.dims() << endl;
  cout << " Vector has size: " << modalVector.size() << endl;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void generateSkinnedPCACubature(MATERIAL**& materials, int totalMaterials)
{
  // read in the partitioned version of the mesh
  cout << "==================================================" << endl;
  cout << " Reading partitions " << endl;
  cout << "==================================================" << endl;
  char buffer[256];

  // load in the LMA and derivative modes
  MATRIX lmaMatrix;
  VECTOR lmaVector;
  MATRIX modalMatrix;
  VECTOR modalVector;
  loadLMA(lmaMatrix, lmaVector, modalMatrix, modalVector);

  // load in the skeleton
  string skeletonFile = posePath + string("ode.motion.0000.skeleton");
  skeleton = new SKELETON(skeletonFile, true);
  int totalPartitions = skeleton->totalBones();

  // load the mesh and integrator
  partitionedMesh = new PARTITIONED_SKINNED_SUBSPACE_TET_MESH(skeleton, tetMeshName.c_str(), materials, totalMaterials, totalPartitions, 0, false, true, tetMeshName);
  integrator = new FULLSPACE_INTEGRATOR(partitionedMesh->originalMesh());

  // associate the mesh with the integrator
  skeleton->tetMesh() = partitionedMesh->originalMesh();
  skeleton->buildOdeSkinning(partitionedMesh->originalMesh());

  // if there is no limit on max key tets, set it to the tet mesh size
  if (maxKeyTets < 0)
    maxKeyTets = partitionedMesh->originalMesh()->totalTets();

  if (generator == NULL)
    generator = new PARTITIONED_CUBATURE_GENERATOR(partitionedMesh);
  else
    generator->reinitialize(partitionedMesh);

  generator->maxKeyTets() = maxKeyTets;
  generator->candidatesPerTry() = candidatesPerTry;
  generator->errorTolerance() = errorTolerance;
  generator->maxIterationsNNLS() = 2;
  generator->randSeed(randSeed);

  // determine the max subrank size
  subrank = (subrank == -1) ? rank: subrank;

  // scatter the basis to the partitions
  cout << " Distributing subbases with max ranks: " << endl;
  vector<int> maxRanks(totalPartitions);
  for (int x = 0; x < totalPartitions; x++)
  {
    maxRanks[x] = subrank;
    cout << maxRanks[x] << endl;
  }
  cout << endl;

  // run Gram-Schmidt on all the sub-basis snapshots
  cout << "==================================================" << endl;
  cout << " Filtering with Gram-Schmidt" << endl;
  cout << "==================================================" << endl;
  map<int, bool>* usedPositions = new map<int,bool>[totalPartitions];
  map<int, bool>* usedVelocities = new map<int,bool>[totalPartitions];

  // do PCA one at a time so full snapshots are only in core one partition at a time
  for (int x = 0; x < totalPartitions; x++)
  {
    // TODO: Make it see multiple paths
    buildSkinnedGramSchmidtBasis(x, maxRanks[x], usedPositions[x], usedVelocities[x]);

    // weight all the subbases by their R values
    SUBSPACE_TET_MESH* subMesh = (SUBSPACE_TET_MESH*)partitionedMesh->mesh(x);
    MATRIX& U = subMesh->U();
    VECTOR& eigenvalues = subMesh->eigenvalues();
  
    for (int j = 0; j < U.cols(); j++)
      for (int i = 0; i < U.rows(); i++)
        U(i,j) *= eigenvalues[j];

    cout << " Calling PCA on partition " << x << endl;
    // run PCA on this subbasis
    generator->pcaSubbasis(maxRanks, x);
  }

  // make a backup of the Us -- we need to swap in the LMA for a moment to PCA
  // those as well
  vector<MATRIX> pcaUs;
  vector<VECTOR> pcaEigs;
  for (int x = 0; x < totalPartitions; x++)
  {
    SUBSPACE_TET_MESH* subMesh = (SUBSPACE_TET_MESH*)partitionedMesh->mesh(x);
    VECTOR& eigenvalues = subMesh->eigenvalues();
    pcaUs.push_back(subMesh->U());
    pcaEigs.push_back(subMesh->eigenvalues());
  }

  // swap in the LMA modes
  SUBSPACE_TET_MESH* originalMesh = (SUBSPACE_TET_MESH*)partitionedMesh->originalMesh();
  originalMesh->U() = lmaMatrix;
  originalMesh->eigenvalues() = lmaVector;

  // PCA the LMA + modal basis
  cout << "===================================================" << endl;
  cout << " Running PCA on LMA and modal derivatives " << endl;
  cout << "===================================================" << endl;
  generator->distributeSubbases(modalMatrix, modalVector, maxRanks);

  cout << "==================================================" << endl;
  cout << " Adding LMA to the basis " << endl;
  cout << "==================================================" << endl;
  // enrich the basis with the LMA and modal derivatives
  for (int x = 0; x < totalPartitions; x++)
  {
    SUBSPACE_TET_MESH* subMesh = (SUBSPACE_TET_MESH*)partitionedMesh->mesh(x);
    MATRIX lmaU = subMesh->U();
    VECTOR lmaEig = subMesh->eigenvalues();
    MATRIX& pcaU = pcaUs[x];
    VECTOR& pcaEig = pcaEigs[x];

    if ((pcaU.cols() % 2) == 1)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " Must select a rank that is even! " << endl;
      cout << " rank: " << pcaU.cols() << endl;
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      exit(0);
    }

    // create the enriched U
    // load the LMA up first to facilitate generating their force samples
    MATRIX enrichedU(lmaU);
    VECTOR enrichedEigs(lmaEig);
    int half = lmaU.cols() / 2;
    for (int y = 0; y < half; y++)
    {
      // write in the LMA column
      VECTOR pcaColumn = pcaU.getColumn(y);
      enrichedU.setColumn(pcaColumn, half + y);

      // clamp the PCA eigs to whatever the last LMA eig was
      enrichedEigs[half + y] = enrichedEigs[half - 1];
      //enrichedEigs[half + y] = pcaEig[y];
    }

    // orthogonalize it
    enrichedU.orthogonalize();

    // save it out the mesh
    subMesh->U() = enrichedU;
    subMesh->eigenvalues() = enrichedEigs;
  }

  cout << "==================================================" << endl;
  cout << " Generating training samples" << endl;
  cout << "==================================================" << endl;

  // stomp any previous samples
  for (int x = 0; x < totalPartitions; x++)
  {
    generator->trainingForcesPartitioned(x).clear();
    generator->trainingQsPartitioned(x).clear();
    generator->trainingColumnsPartitioned(x).clear();
    generator->trainingConstraintsPartitioned(x).clear();
  }

  // generate the LMA-based training samples
  for (int x = 0; x < totalPartitions; x++)
  {
    SUBSPACE_TET_MESH* subMesh = (SUBSPACE_TET_MESH*)partitionedMesh->mesh(x);
    MATRIX lmaU = subMesh->U();
    int all = lmaU.cols();
    int lmaSamples = usedPositions[x].size();
    generateLMATrainingSamples(x, 100, 1.0, all);
  }

  // generate the traning set for each partition
  for (int x = 0; x < totalPartitions; x++)
  {
    cout << " Generating sample for partition " << x << endl;
    generatePCATrainingSamples(x, usedPositions[x], usedVelocities[x]);
  }

  // save the scattered bases to disk
  partitionedMesh->writeSubbases();

  cout << "==================================================" << endl;
  cout << " Training cubature" << endl;
  cout << "==================================================" << endl;
  generator->generatePCACubature();
  generator->writePartitionedCubatures();

  delete partitionedMesh;
  partitionedMesh = NULL;
  delete integrator;
  delete generator;
  generator = NULL;
  delete[] usedPositions;
  delete[] usedVelocities;
  delete skeleton;
}

//////////////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  if (argc <= 1)
  {
    cout << " USAGE: " << argv[0] << " *.cfg " << endl;
    return 0;
  }
  
  // read in different parameters
  string configName(argv[1]);
  SIMPLE_PARSER configFile(configName); 
  triangleMeshPath = configFile.getString("triangle mesh path", triangleMeshPath);
  triangleMeshName = configFile.getString("triangle mesh name", triangleMeshName);
  outputPath       = configFile.getString("output path", outputPath);
  snapshots        = configFile.getInt("snapshots", -1);
  startingSnapshot = configFile.getInt("starting snapshot", -1);
  bccRes           = configFile.getInt("bccRes", bccRes);
  posePath = configFile.getString("pose path", posePath);
  string dataPath;
  dataPath      = configFile.getString("data path", dataPath);
  startingSnapshot++;

  if (snapshots < 0)
  {
    cout << " Config file needs to contain number of snapshots to factor!" << endl;
    exit(0);
  }

  if (configFile.defined("tet mesh name"))
    tetMeshName = configFile.getString("tet mesh name", tetMeshName);
  else
    tetMeshName = outputPath + triangleMeshName + string(".tetmesh");

  tetMeshName = tetMeshName + string(".constrained");

  // read in how many materials there are
  int totalMaterials = 0;
  totalMaterials = configFile.getInt("total materials", totalMaterials);
  if (totalMaterials == 0)
  {
    cout << " NO MATERIALS SPECIFIED!!!!" << endl;
    exit(1);
  }
  
  // read in the actual materials
  MATERIAL** materials = new MATERIAL*[totalMaterials];
  bool invertible = false;
  invertible = configFile.getBool("invertible", invertible);
  for (int x = 0; x < totalMaterials; x++)
  {
    // read in the config file name for the material
    char buffer[256];
    sprintf(buffer, "material %i", x);
    string materialString(buffer);
    string materialFilename;
    materialFilename = configFile.getString(materialString.c_str(), materialFilename);

    // open the config file
    SIMPLE_PARSER materialFile(materialFilename);
    materials[x] = readMaterial(materialFile);
    
    // set the invertible wrapper if necessary
    if (invertible)
    {
      materials[x] = new INVERTIBLE(materials[x]);
      cout << " Setting material to invertible" << endl;
      float inversionEpsilon = configFile.getFloat("inversion epsilon", ((INVERTIBLE*)materials[x])->epsilon());
      cout << " Inversion epsilon: " << inversionEpsilon << endl;
    }
  }

  // set training parameters
  rank                 = configFile.getInt("rank", 20);
  subrank              = configFile.getInt("subrank", -1);
  int trainingSamples      = configFile.getInt("training samples", 100);
  double trainingMagnitude = configFile.getFloat("training magnitude", 1.0);
  maxKeyTets           = configFile.getInt("max key tets", 1000);
  candidatesPerTry     = configFile.getInt("candidates per try", 100);
  errorTolerance    = configFile.getFloat("error tolerance", 0.01);
  randSeed             = configFile.getInt("randomization seed", 123456);

  Real meshScale           = configFile.getFloat("mesh scale", 1.0);

  string originalMeshFile = triangleMeshPath + triangleMeshName;

  // armadillo.ragdoll.16.cfg
  dataPaths.push_back(dataPath);

  cout << "======================================================================" << endl;
  cout << " Generating skeleton cubature " << endl;
  cout << "======================================================================" << endl;
  generateSkinnedPCACubature(materials, totalMaterials);

  for (int x = 0; x < totalMaterials; x++)
    delete materials[x];
  delete[] materials;
}
