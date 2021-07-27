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
#include "EIGEN.h"

#include <iostream>
#include <fstream>

#include "FLUID_3D_MIC.h"
#include "MATRIX.h"
#include "BIG_MATRIX.h"
#include "SIMPLE_PARSER.h"

//////////////////////////////////////////////////////////////////////////////
// Fix up the slashes in a string to play nice with Perl
//////////////////////////////////////////////////////////////////////////////
string fixSlashes(string& input)
{
  string final;

  for (unsigned int x = 0; x < input.size(); x++)
  {
    if (input[x] == '/')
    {
      final.push_back('\\');
    }
    final.push_back(input[x]);
  }
  return final;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void verifyResults(string filenamePrefix, string finalU)
{
  MatrixXd U;
  EIGEN::read(finalU, U);

  string filename = filenamePrefix + string(".transpose");
  FILE* file = fopen(filename.c_str(), "rb");

  for (int x = 0; x < 10; x++)
  {
    VectorXd column(U.rows());
    EIGEN::readRaw(file, U.rows(), column);
    VectorXd diff = column - U * (U.transpose() * column);

    cout << " Diff " << x << ": " << diff.norm() << endl;
  }

  fclose(file);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  cout << " Size of Real: " << sizeof(Real) << endl;

  TIMER functionTimer(__FUNCTION__);
  cout << "=======================================================================================" << endl;
  cout << argv[0] << " " << argv[1] << endl;
  cout << "=======================================================================================" << endl;

  // read in the cfg file
  if (argc != 2)
  {
    cout << " Usage: " << argv[0] << " *.cfg" << endl;
    return 0;
  }
  SIMPLE_PARSER parser(argv[1]);
  string snapshotPath = parser.getString("snapshot path", "./data/dummy/");
  string reducedPath = parser.getString("reduced path", "./data/reduced.dummy/");
  string scratchPath = parser.getString("scratch path", "./data/scratch/");
  int simulationSnapshots = parser.getInt("simulation snapshots", 20);
  int reducedSnapshots = parser.getInt("reduced snapshots", 20);
  double discardThreshold = parser.getFloat("discard threshold", 1e-9);

  bool usingIOP = parser.getBool("iop", false);

  if (usingIOP)
    cout << " Using IOP " << endl;
  else
    cout << " NOT using IOP " << endl;

  cout << " Using discard threshold: " << discardThreshold << endl;

  // set scratch path
  cout << " Using scratch path: " << scratchPath.c_str() << endl;

  if (reducedSnapshots > simulationSnapshots)
  {
    cout << " You asked to use more snapshots than were simulated! " << endl;
    exit(0);
  }

  // make sure the reduced directory in fact exists
  string mkdir("mkdir ");
  mkdir = mkdir + reducedPath;
  system(mkdir.c_str());

  vector<string> filenamePrefixes;
  vector<string> finalFilenames;
  filenamePrefixes.push_back(snapshotPath + string("velocity.final.matrix"));
  filenamePrefixes.push_back(snapshotPath + string("velocity.preproject.matrix"));
  filenamePrefixes.push_back(snapshotPath + string("velocity.prediffuse.matrix"));
  filenamePrefixes.push_back(snapshotPath + string("velocity.preadvect.matrix"));
  filenamePrefixes.push_back(snapshotPath + string("pressure.matrix"));
  if (usingIOP)
    filenamePrefixes.push_back(snapshotPath + string("velocity.iop.matrix"));
  
  finalFilenames.push_back(reducedPath + string("U.final.matrix"));
  finalFilenames.push_back(reducedPath + string("U.preproject.matrix"));
  finalFilenames.push_back(reducedPath + string("U.prediffuse.matrix"));
  finalFilenames.push_back(reducedPath + string("U.preadvect.matrix"));
  finalFilenames.push_back(reducedPath + string("U.pressure.matrix"));
  if (usingIOP)
    finalFilenames.push_back(reducedPath + string("U.iop.matrix"));

  vector<string> pcaFilenames;
  pcaFilenames.push_back(reducedPath + string("pca.final.vector"));
  pcaFilenames.push_back(reducedPath + string("pca.preproject.vector"));
  pcaFilenames.push_back(reducedPath + string("pca.prediffuse.vector"));
  pcaFilenames.push_back(reducedPath + string("pca.preadvect.vector"));
  pcaFilenames.push_back(reducedPath + string("pca.pressure.vector"));
  if (usingIOP)
    pcaFilenames.push_back(reducedPath + string("pca.iop.vector"));

  for (unsigned int x = 0; x < filenamePrefixes.size(); x++)
  {
    BIG_MATRIX::outOfCoreSVD(filenamePrefixes[x], reducedPath, discardThreshold);
  
    // back up the final U matrix SVD values so we can build reduced accuracy versions later
    string cp("cp ");
    cp = cp + BIG_MATRIX::scratchPath() + string("svd.vector ") + pcaFilenames[x];
    system(cp.c_str());

    BIG_MATRIX::writeFinalU(finalFilenames[x]);
    TIMER::printTimings();

    verifyResults(filenamePrefixes[x], finalFilenames[x]);
  }

  return EXIT_SUCCESS;
}
