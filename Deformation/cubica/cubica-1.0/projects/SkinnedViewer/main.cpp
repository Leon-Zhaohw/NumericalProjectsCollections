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
#include <SUBSPACE_TET_MESH.h>
#include <SPARSE_MATRIX.h>
#include <STVK.h>
#include <MOONEY_RIVLIN.h>
#include <ARRUDA_BOYCE.h>
#include <NEO_HOOKEAN.h>
#include <INVERTIBLE.h>
#include <SIMPLE_PARSER.h>
#include <FULLSPACE_INTEGRATOR.h>
#include <ONLINE_SUBSPACE_INTEGRATOR.h>
#include <SKELETON.h>
#include <TENSOR3.h>

#define KEYFRAMED_MOTION 0

using namespace std;
SUBSPACE_TET_MESH* tetMesh;
ONLINE_SUBSPACE_INTEGRATOR* integrator;
SKELETON* skeleton = NULL;
bool invertible = true;
int basisType = -1;
Real onlineErrorTolerance;
Real consecutiveAmp = 1.0;
string oracleVideo;
string outputPath;

bool regressionTesting = false;

// create a dummy FULLSPACE_INTEGRATOR for rendering the ground truth solution
FULLSPACE_INTEGRATOR* groundTruth;

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

    // do a version check to make sure we're not using the old inversion epsilon
    Real inversionEpsilon = parser.getFloat("inversion epsilon", -1);
    if (inversionEpsilon > 0.0)
    {
      cout << " USING STALE VERSION OF INVERSION " << endl;
      exit(1);
    }
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
//////////////////////////////////////////////////////////////////////////////
Real calcDistanceToOrtho()
{
  MATRIX& U = tetMesh->U();
  MATRIX check = U ^ U;
  // get the max off diagonal entry
  for (int y = 0; y < check.cols(); y++)
    check(y,y) = 0.0;
  return check.maxAbsEntry();
}

//////////////////////////////////////////////////////////////////////////////
// head simulation with Hessian error estimator
//////////////////////////////////////////////////////////////////////////////
void headFoldedErrorSimulation(string renderPath, string dataPath, string videoname, Real discardThreshold, bool invertible)
{
  integrator->resetBasis();
  skeleton->reset();

  // set the discard threshold
  integrator->snapshotDiscardThreshold() = discardThreshold;
  cout << " Setting discard threshold to " << discardThreshold << endl;

  bool fullActive = true;

  string filename;
  char buffer[256];

#if KEYFRAMED_MOTION  
  // new motion
  //int totalSteps = 600;
  int totalSteps = 600 * 0.01 / integrator->dt();
#else
  // old motion
  int totalSteps = 1000.0 * 0.01 / integrator->dt();
#endif
  cout << " Total steps: " << totalSteps << endl;

  float time = 0.0;
  int skippedSteps = 0;
  
  int consecutiveSubsteps = 0;
  int maxConsecutiveSubsteps = 0;

  VECTOR skipped(totalSteps);
  VECTOR addedToBasis(totalSteps);
  VECTOR earlyExits(totalSteps);
  FULLSPACE_INTEGRATOR* fullIntegrator = integrator->fullIntegrator();
  fullIntegrator->reset();
  SUBSPACE_INTEGRATOR* subspaceIntegrator = integrator->subspaceIntegrator();
  Real errorThreshold = integrator->trueErrorThreshold();
  SPARSE_PETSC_MATRIX& petscSolver = fullIntegrator->solver();

  vector<Real> distToOrtho;
  vector<Real> errorNorm;
  bool errorColumnPresent = false;
  bool firstSkip = false;
  for (int x = 0; x < totalSteps; x++)
  {
    // head motion skinning
    Real amp = 1.0;
    bool forceFull = false;
    Real radianTime = 2.0 * M_PI * time;
    BONE* selected = skeleton->selected();
    MATRIX3 rotation = selected->rotation();
#if KEYFRAMED_MOTION
    // new motion
    skeleton->headShakeKeyframed(time);
#else
    // old motion
    rotation = skeleton->headShake(time, 1.0, forceFull);
#endif
    skeleton->updateSkeleton();
    time += integrator->dt();

    int currentRank = tetMesh->rank();

    // step the integrator
    integrator->stepQuasistatic(invertible, false, consecutiveAmp);

    /*
    // dump the state
    string fullfile = dataPath;
    fullfile += string("online.integrator.");
    sprintf(buffer, "%i", x);
    fullfile += string(buffer);
    fullfile += string(".state");
    fullIntegrator->writeState(fullfile);
    */

    // dump the tet mesh state
    VECTOR backup = tetMesh->TET_MESH::x();
    tetMesh->x() = fullIntegrator->position();
    skeleton->addSkinningDisplacements();
    //tetMesh->writeDeformedMesh(x, dataPath);

    skeleton->undoSkinningDisplacements();
    // restore the tet mesh state to whatever it was before
    tetMesh->x() = backup;
    tetMesh->TET_MESH::updateFullMesh();

    if (x == 30 && regressionTesting)
    {
      VECTOR skips(integrator->skippedSteps());

      string skipDump = outputPath + "skips.vector";

      cout << " Dumping skips to file " << skipDump.c_str() << endl;
      skips.write(skipDump.c_str());

      exit(0);
    }

    /*
    // do the render
    filename = renderPath + string("render.position.");
    filename += itoPaddedString(x);
    filename += string(".tif");
    renderMan.drawHeadBegin(filename);
      integrator->forceRenderManDrawFull();
    renderMan.drawEnd(filename);

    if (tetMesh->rank() != currentRank || x % 10 == 0)
      integrator->printTimingBreakdown();
      */
    //if (x == 150) break;
    //if (x == 50) break;
  }
  //cout << " DONE. Discard: " << discardThreshold << " Error: " << errorThreshold << endl;
  //integrator->printTimingBreakdown();

  /*
  string encoder = string("encode ");
  encoder = encoder + renderPath + string(" ");
  encoder = encoder + oracleVideo;
  system(encoder.c_str());
  */
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
  // input parameters
  string triangleMeshPath;
  string triangleMeshName;
  string tetMeshName;
  string renderPath = string("./");
  string dataPath  = string("./");
  string logPath  = string("./");
  bool groundTruthRun = false;

  // read in different parameters
  string configName(argv[1]);
  SIMPLE_PARSER configFile(configName); 
  triangleMeshPath = configFile.getString("triangle mesh path", triangleMeshPath);
  triangleMeshName = configFile.getString("triangle mesh name", triangleMeshName);
  outputPath       = configFile.getString("output path", outputPath);
  renderPath       = configFile.getString("render path", renderPath);
  dataPath         = configFile.getString("data path", dataPath);
  logPath         = configFile.getString("log path", logPath);
  groundTruthRun   = configFile.getBool("ground truth", groundTruthRun);

  if (configFile.defined("tet mesh name"))
    tetMeshName = configFile.getString("tet mesh name", tetMeshName);
  else
    tetMeshName = outputPath + triangleMeshName + string(".tetmesh");

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
    
    // get the material
    MATERIAL* material = readMaterial(materialFile);
    
    // set the invertible wrapper if necessary
    if (invertible)
    {
      cout << " Setting material to invertible" << endl;
      material = new INVERTIBLE(material);
    }

    materials[x] = material;
  }

  // output precision
  cout << " Precision: " << sizeof(Real) * 8 << " bit" << endl;

  // get gravity constant
  Real gravity = 9.8;
  if (configFile.defined("gravity"))
    gravity = configFile.getFloat("gravity", gravity);
  cout << " Gravity: " << gravity << endl;
  
  // see which simulate we are doing
  int whichSimulation = configFile.getInt("which simulation", -1);

  // get rayleigh damping consts
  Real rayleighAlpha;
  Real rayleighBeta;
  if (whichSimulation == 0)
    rayleighAlpha = rayleighBeta = 0.1;
  else if (whichSimulation == 1)
    rayleighAlpha = rayleighBeta = 0.001;
  else if (whichSimulation == 2)
    rayleighAlpha = rayleighBeta = 0.0005;
  rayleighAlpha = configFile.getFloat("rayleigh alpha", rayleighAlpha);
  rayleighBeta  = configFile.getFloat("rayleigh beta", rayleighBeta);

  cout << "==================================================" << endl;
  cout << " Reading mesh file " << tetMeshName.c_str() << endl;
  cout << "==================================================" << endl;

  Real timestep = 1.0f / 60.0f;
  timestep = configFile.getFloat("timestep", timestep);
  cout << " Using timestep: " << timestep << endl;
  tetMesh = new SUBSPACE_TET_MESH(tetMeshName.c_str(), materials, totalMaterials, true);
  integrator = new ONLINE_SUBSPACE_INTEGRATOR(tetMesh, timestep, gravity, rayleighAlpha, rayleighBeta);

  // set solver precision
  int newtonDigits = 2;
  newtonDigits = configFile.getInt("Newton digits of precision", newtonDigits);
  Real precision = pow(10.0, (Real)(-newtonDigits));
  cout << " Setting Newton solver precision to : " << precision << " (" << newtonDigits << " digits)" << endl;
  integrator->setSolverEps(precision);
  int digitsPCG = 8;
  digitsPCG = configFile.getInt("PCG digits of precision", digitsPCG);
  cout << " Setting PCG solver precision to " << digitsPCG << " digits" << endl;
  integrator->setPCGDigits(digitsPCG);

  // do Krysl-style stiffness matrices?
  bool useKrysl = false;
  useKrysl = configFile.getBool("use krysl stiffness", useKrysl);
  if (useKrysl)
    cout << " Using Krysl-style stiffness matrices" << endl;
  else
    cout << " Not using Krysl-style stiffness matrices" << endl;
  integrator->useKryslStiffness() = useKrysl;

  // amount of per-step error tolerated for the subspace
  onlineErrorTolerance = configFile.getFloat("online error tolerance", 0.01);
  integrator->setTrueErrorThreshold(onlineErrorTolerance);
  cout << " Setting online error threshold to " << onlineErrorTolerance << endl;

  // set the max number of Newton iterations allowed
  int maxNewtonSteps = 3;
  maxNewtonSteps = configFile.getInt("max newton steps", maxNewtonSteps);
  integrator->setMaxNewtonSteps(maxNewtonSteps);
  cout << " Max Newton steps allowed: " << maxNewtonSteps << endl;

  // what rendering mode?
  bool renderBoth = true;
  renderBoth = configFile.getBool("render full and reduced", renderBoth);
  integrator->renderBoth() = renderBoth;

  // set training parameters
  int bccRes               = configFile.getInt("bccRes", 32);
  int rank                 = configFile.getInt("rank", 20);
  int trainingSamples      = configFile.getInt("training samples", 100);
  double trainingMagnitude = configFile.getFloat("training magnitude", 1.0);
  int maxKeyTets           = configFile.getInt("max key tets", 1000);
  int candidatesPerTry     = configFile.getInt("candidates per try", 100);
  double errorTolerance    = configFile.getFloat("error tolerance", 0.01);
  int randSeed             = configFile.getInt("randomization seed", 123456);
  ONLINE_CUBATURE_GENERATOR* generator = integrator->cubatureGenerator();
  generator->maxKeyTets() = maxKeyTets;
  generator->candidatesPerTry() = candidatesPerTry;
  //generator->errorTolerance() = errorTolerance;
  generator->errorTolerance() = onlineErrorTolerance;
  generator->randSeed(randSeed);  

  // set rank
  integrator->maxRank() = rank;
  integrator->downdateSize() = integrator->maxRank() / 2;
  
  // factor to increase the mouse-input force by
  double forceMultiplier = 0.1;
  forceMultiplier = configFile.getFloat("force multiplier", forceMultiplier);
  integrator->setForceMultiplier(forceMultiplier);

  // tell integrator the BCC res for fishnet forces
  integrator->bccRes() = bccRes;
  cout << "  Original BCC res: " << bccRes << endl;

  // see if a direct discard threshold was specified
  double discardThreshold = -1.0;
  discardThreshold = configFile.getFloat("snapshot discard threshold", discardThreshold);
  if (discardThreshold > 0.0)
  {
    integrator->snapshotDiscardThreshold() = discardThreshold;
    cout << " Snapshot discard threshold set directly to " << integrator->snapshotDiscardThreshold() << endl;
  }

  // doing a regression test?
  regressionTesting = configFile.getBool("regression testing", regressionTesting);

  // create the data and render directories
  string mkdirRender = string("mkdir ") + renderPath;
  string mkdirData = string("mkdir ") + dataPath;
  system(mkdirRender.c_str());
  system(mkdirData.c_str());
  integrator->logFileDirectory() = logPath;

  // always render to some temp directory
  renderPath = string("./renders/");
  renderPath = configFile.getString("render path", renderPath);

  basisType = configFile.getInt("basis type", basisType);
  if (basisType == -1)
  {
    cout << " MUST SELECT A BASIS TYPE" << endl;
    exit(1);
  }
  if (basisType == 0)
    cout << " Using position basis" << endl;
  if (basisType == 1)
    cout << " Using velocity basis" << endl;
  if (basisType == 2)
    cout << " Using acceleration basis" << endl;

  double pushStrength = configFile.getFloat("push strength", 2.0);
  cout << " Pendulum push strength: " << pushStrength << endl;

  Real gravityStrength = 2.0;
  gravityStrength = configFile.getFloat("gravity strength", gravityStrength);
  cout << " Pendulum gravity strength: " << gravityStrength << endl;

  consecutiveAmp = configFile.getFloat("consecutive amp", consecutiveAmp);
  cout << " Consecutive amplification: " << consecutiveAmp << endl;

  if (whichSimulation == 10)
  {
    // hack the vertex positions
    vector<VEC3F>& vertices = tetMesh->vertices();
    Real maxY = vertices[0][1];
    Real minY = vertices[0][1];
    for (int x = 0; x < vertices.size(); x++)
    {
      vertices[x][0] -= 0.5;
      vertices[x][2] -= 0.5;
      if (vertices[x][1] > maxY) maxY = vertices[x][1];
      if (vertices[x][1] < minY) minY = vertices[x][1];
    }
    vector<VEC3F>& restVertices = tetMesh->restPose();
    for (int x = 0; x < restVertices.size(); x++)
    {
      restVertices[x][0] -= 0.5;
      restVertices[x][2] -= 0.5;
      if (restVertices[x][1] > maxY) maxY = restVertices[x][1];
      if (restVertices[x][1] < minY) minY = restVertices[x][1];
    }

    // scale by the y size
    for (int x = 0; x < restVertices.size(); x++)
    {
      restVertices[x][1] -= minY;
      restVertices[x][1] *= 1.0 / (maxY - minY);
      vertices[x][1] -= minY;
      vertices[x][1] *= 1.0 / (maxY - minY);
    }

    // reinit tets with these positions
    vector<TET>& tets = tetMesh->tets();
    for (int x = 0; x < tets.size(); x++)
      tets[x].init();
  }

  // construct the skeleton
  skeleton = new SKELETON(tetMesh, bccRes);
  skeleton->meshOrigin() = VEC3F(0.0, 0.0, 0.0);
  if (whichSimulation == 10)
    skeleton->buildSkeleton();
  else
    skeleton->buildHeadSkeleton();
  skeleton->readSkinning();
  skeleton->updateSkeleton();

  // select the second joint
  skeleton->cycleSelection();
  integrator->skeleton() = skeleton;
  generator->skeleton() = skeleton;

  ////////////////////////////////////////////////////////////////////
  // do online simulation
  vector<VECTOR> errors;
  vector<int> finalRanks;
  string videoprefix("dummy.avi");
  videoprefix = configFile.getString("video prefix", videoprefix);

  // create a dummy FULLSPACE_INTEGRATOR for rendering the ground truth solution
  groundTruth = new FULLSPACE_INTEGRATOR(tetMesh, integrator->dt(), integrator->rayleighAlpha(), integrator->rayleighBeta(), 0.0);

  oracleVideo = configFile.getString("oracle video", oracleVideo);

  //Real discard = 0.0001;
  headFoldedErrorSimulation(renderPath, dataPath, videoprefix, discardThreshold, invertible);

  delete groundTruth;
  ////////////////////////////////////////////////////////////////////

  /*
  ////////////////////////////////////////////////////////////////////
  // do a full simulation
  //pendulumFullSimulation(renderPath, dataPath, groundTruthRun, bccRes, invertible, pushStrength, gravityStrength);
  //headFullSimulation(renderPath, dataPath, groundTruthRun, bccRes, invertible, pushStrength);
  headFullQuasistaticSimulation(renderPath, dataPath, groundTruthRun, bccRes, invertible, pushStrength);
  //fullSimulation(renderPath, dataPath, groundTruthRun, bccRes, invertible, pushStrength);
  //hessianTestSimulation(renderPath, dataPath, groundTruthRun, bccRes, invertible, pushStrength);
  string oracleVideo("dummy.avi");
  oracleVideo = configFile.getString("oracle video", oracleVideo);

  //string encoder = string("encodeoracle ");
  string encoder = string("encode ");
  //string encoder = string("encodetandem ");
  //string encoder = string("encodereplay ");
  encoder = encoder + renderPath + string(" ");
  encoder = encoder + oracleVideo;
  system(encoder.c_str());
  ////////////////////////////////////////////////////////////////////
  */
}
