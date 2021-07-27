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

/*
 * Multi-dimensional DCT testing
 * Aaron Demby Jones
 * Fall 2014
 */

#include <iostream>
#include <fftw3.h>
#include "EIGEN.h"
#include "SUBSPACE_FLUID_3D_EIGEN.h"
#include "FLUID_3D_MIC.h"
#include "CUBATURE_GENERATOR_EIGEN.h"
#include "MATRIX.h"
#include "BIG_MATRIX.h"
#include "SIMPLE_PARSER.h"
#include "COMPRESSION.h"
#include <string>
#include <cmath>
#include <cfenv>
#include <climits>
#include "VECTOR3_FIELD_3D.h"

using std::vector;
using std::string;


////////////////////////////////////////////////////////
// Function Declarations
////////////////////////////////////////////////////////

// set the damping matrix and compute the number of blocks
void PreprocessEncoder(COMPRESSION_DATA* data0, COMPRESSION_DATA* data1, COMPRESSION_DATA* data2, int maxIterations, const char* filename);

// rescale the singular values to use for damping
void PreprocessSingularValues(const char* filename, double threshold); 

// check if a file exists
bool fileExists(const string& filename);

// get the length of a file in bytes
unsigned long long FileSize(const string& filename);

// compute and print the compression ratios
double GetCompressionRatios(const string& preadvectFilename, const string& finalFilename);

// update the cfg file to point to the correct compression path
void UpdateCfgFile(int roundedOverallCompression);

string preadvectPath;
string finalPath;
string cfgFilename;

////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {
  TIMER functionTimer(__FUNCTION__);
  
  // read in the cfg file
  if (argc != 2) {
    cout << " Usage: " << argv[0] << " *.cfg" << endl;
    return 0;
  }
  
  SIMPLE_PARSER parser(argv[1]);
  cfgFilename = argv[1];
  string reducedPath = parser.getString("reduced path", "./data/reduced.dummy/"); 
  int xRes = parser.getInt("xRes", 48);
  int yRes = parser.getInt("yRes", 64);
  int zRes = parser.getInt("zRes", 48);
  int numCols = parser.getInt("reduced snapshots", 50);
  bool usingIOP = parser.getBool("iop", 0);
  cout << " Using IOP: " << usingIOP << endl;
  cout << " Block size: " << BLOCK_SIZE << endl; 
  // we want the peeled resolutions for the matrices
  xRes -= 2;
  yRes -= 2;
  zRes -= 2;

  VEC3I dims(xRes, yRes, zRes);
  
  // times 3 since it is a VELOCITY3_FIELD_3D flattened out
  int numRows = 3 * xRes * yRes * zRes;
  cout << "numRows: " << numRows << endl;
  cout << "numCols: "<< numCols << endl;

  MatrixXd U_preadvect(numRows, numCols);
  MatrixXd U_final(numRows, numCols);

  int nBits = parser.getInt("nBits", 24); 
  cout << " nBits: " << nBits << endl;
  double percent = parser.getFloat("percent", 0.99);
  cout << " percent: " << percent << endl;
  int maxIterations = parser.getInt("maxIterations", 32);
  cout << " maxIterations: " << maxIterations << endl;
  bool debug = parser.getBool("debug", false);
  cout << "debug: " << debug << endl;

  bool usingFastPow = parser.getBool("fastPow", false);
  cout << " fastPow: " << usingFastPow << endl;
  FIELD_3D::usingFastPow() = usingFastPow;

  preadvectPath = reducedPath + string("U.preadvect.matrix");
  finalPath = reducedPath + string("U.final.matrix");

  EIGEN::read(preadvectPath, U_preadvect);
  EIGEN::read(finalPath, U_final);

  // set the parameters in compression data
  COMPRESSION_DATA preadvect_compression_data0(dims, numCols, nBits, percent);
  COMPRESSION_DATA preadvect_compression_data1(dims, numCols, nBits, percent);
  COMPRESSION_DATA preadvect_compression_data2(dims, numCols, nBits, percent);
  COMPRESSION_DATA final_compression_data0(dims, numCols, nBits, percent);
  COMPRESSION_DATA final_compression_data1(dims, numCols, nBits, percent);
  COMPRESSION_DATA final_compression_data2(dims, numCols, nBits, percent);

  /*
  // compute some additional parameters for compression data
  const char* preadvectSingularFilename = "singularValues_preadvect.vector";
  PreprocessEncoder(&preadvect_compression_data0, &preadvect_compression_data1, &preadvect_compression_data2, 
      maxIterations, preadvectSingularFilename);
  const char* finalSingularFilename = "singularValues_final.vector";
  PreprocessEncoder(&final_compression_data0, &final_compression_data1, &final_compression_data2,
      maxIterations, finalSingularFilename);
  */

  // ADJ: change this threshold to modulate singular value damping
  
  
  const double threshold = 1.0; 
  // ADJ: change the scratch path to SSD for big runs!
  string scratchPath = "./scratch/";
  
  string preadvectSingularFilename = scratchPath + string("velocity.preadvect.matrix.singularValues.vector");

  /*
  string preadvectProcessed = preadvectSingularFilename + string(".processed");
  if (!fileExists(preadvectProcessed)) {
    puts("Preadvect singular values are unprocessed; processing now...");
    printf("Threshold equals: %f\n", threshold);
    PreprocessSingularValues(preadvectSingularFilename.c_str(), threshold);
    puts("Done.");
  }
  
  else { 
    puts("Rewriting over previously pre-processed preadvect singular values...");
    printf("Threshold equals: %f\n", threshold);
    PreprocessSingularValues(preadvectSingularFilename.c_str(), threshold);
  }
  */

  PreprocessEncoder(&preadvect_compression_data0, &preadvect_compression_data1, &preadvect_compression_data2, 
      maxIterations, preadvectSingularFilename.c_str());
  

  string finalSingularFilename = scratchPath + string("velocity.final.matrix.singularValues.vector");
  /*
  string finalProcessed = finalSingularFilename + string(".processed");
  if (!fileExists(finalProcessed)) {
    puts("FInal singular values are unprocessed; processing now...");
    printf("Threshold equals: %f\n", threshold);
    PreprocessSingularValues(finalSingularFilename.c_str(), threshold);
  }

  else {
    puts("Rewriting over previously pre-processed final singular values...");
    printf("Threshold equals: %f\n", threshold);
    PreprocessSingularValues(finalSingularFilename.c_str(), threshold);
  }
  */

  PreprocessEncoder(&final_compression_data0, &final_compression_data1, &final_compression_data2,
      maxIterations, finalSingularFilename.c_str());

  // write a binary file for each scalar field component

  string command = string("mkdir ") + reducedPath + string("tmp");
  system(command.c_str());
  string preadvectFilename = reducedPath + string("tmp/U.preadvect.component");
  string finalFilename = reducedPath + string("tmp/U.final.component");


  // write out the compressed matrix files

  if (debug) {
    CompressAndWriteMatrixComponentsDebug(preadvectFilename.c_str(), U_preadvect, &preadvect_compression_data0,
      &preadvect_compression_data1, &preadvect_compression_data2);

    CompressAndWriteMatrixComponentsDebug(finalFilename.c_str(), U_final, &final_compression_data0, 
      &final_compression_data1, &final_compression_data2);
  }

    else {
      CompressAndWriteMatrixComponents(preadvectFilename.c_str(), U_preadvect, &preadvect_compression_data0,
        &preadvect_compression_data1, &preadvect_compression_data2);
    
      CompressAndWriteMatrixComponents(finalFilename.c_str(), U_final, &final_compression_data0, 
        &final_compression_data1, &final_compression_data2);
      

      // ADJ: this is old experimental code testing out the different DCT/DST hybrids
      /*
      for (int planType = 0; planType < 8; planType++) {

        CompressAndWriteMatrixComponentsDST(preadvectFilename.c_str(), planType, U_preadvect, &preadvect_compression_data0,
          &preadvect_compression_data1, &preadvect_compression_data2);

        CompressAndWriteMatrixComponentsDST(finalFilename.c_str(), planType, U_final, &final_compression_data0,
          &final_compression_data1, &final_compression_data2);
      } 
      */

    }

  double ratio = GetCompressionRatios(preadvectFilename, finalFilename);
  int roundedRatio = rint(ratio);
  string newName = reducedPath + to_string(roundedRatio) + string("to1");
  string rename = string("mv ") + reducedPath + string("tmp ") + newName;
  system(rename.c_str());
  string mkdir = string("mkdir ") + newName + string("/pbrt");
  system(mkdir.c_str()); 

  UpdateCfgFile(roundedRatio);

  TIMER::printTimings();
  
  return 0;
}

void PreprocessEncoder(COMPRESSION_DATA* data0, COMPRESSION_DATA* data1, COMPRESSION_DATA* data2, int maxIterations, const char* filename)
{
  // set integer rounding 'to nearest' 
  fesetround(FE_TONEAREST);
  
  // precompute and set the damping  and zigzag arrays
  data0->set_dampingArray();
  data1->set_dampingArray();
  data2->set_dampingArray();

  data0->set_zigzagArray();
  data1->set_zigzagArray();
  data2->set_zigzagArray();
   
  // fill in the appropriate paddings
  const VEC3I& dims = data0->get_dims();
  VEC3I paddings;
  GetPaddings(dims, &paddings);
  paddings += dims;

  data0->set_paddedDims(paddings);
  data1->set_paddedDims(paddings);
  data2->set_paddedDims(paddings);

  // calculates number of blocks, assuming BLOCK_SIZE 
  int numBlocks = paddings[0] * paddings[1] * paddings[2] / (BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE);

  data0->set_numBlocks(numBlocks);
  data1->set_numBlocks(numBlocks);
  data2->set_numBlocks(numBlocks);

  // set the maxIterations
  data0->set_maxIterations(maxIterations);
  data1->set_maxIterations(maxIterations);
  data2->set_maxIterations(maxIterations);
  
  // set the singular values
  // DEBUG: let's leave this out for now.
  
  /*
  data0->set_singularValues(filename);
  data1->set_singularValues(filename);
  data2->set_singularValues(filename);
  */
}

void PreprocessSingularValues(const char* filename, double threshold)
{
  VECTOR singularValues;
  singularValues.read(filename);
  for (int i = 0; i < singularValues.size(); i++) {
    singularValues[i] = log(singularValues[i]);
  }
  double min = singularValues.min();
  for (int i = 0; i < singularValues.size(); i++) {
    singularValues[i] -= min;
  }
  double s0inv = 1.0 / singularValues[0];
  singularValues *= s0inv;

  for (int i = 0; i < singularValues.size(); i++) {
    singularValues[i] *= (1 - threshold);
    singularValues[i] += threshold;
  }

  string output(filename);
  output += string(".processed");
  singularValues.write(output.c_str());
  printf("Wrote out rescaled singular values!\n");
}
   

//////////////////////////////////////////////////////////////////////
// check if a file exists
//////////////////////////////////////////////////////////////////////
bool fileExists(const string& filename)
{
  FILE* file;
  file = fopen(filename.c_str(), "rb");
  
  if (file == NULL)
    return false;

  fclose(file);
  return true;
}


//////////////////////////////////////////////////////////////////////
// get the size of a file in bytes
//////////////////////////////////////////////////////////////////////
unsigned long long FileSize(const string& filename)
{
  TIMER functionTimer(__FUNCTION__);
  FILE* file;
  file = fopen(filename.c_str(), "rb");
  unsigned long long size;

  if (file == NULL) perror ("Error opening file");
  else {
  fseek(file, 0, SEEK_END);   // non-portable
  size = (unsigned long long) ftell(file);
  fclose (file);
  }
  return size;
}

//////////////////////////////////////////////////////////////////////
// compute and print the compression ratios
//////////////////////////////////////////////////////////////////////
double GetCompressionRatios(const string& preadvectFilename, const string& finalFilename)
{
  TIMER functionTimer(__FUNCTION__);
  puts("Computing compression ratios...");

  auto preadvectSize = FileSize(preadvectPath);
  auto finalSize     = FileSize(finalPath);

  string compressedPreadvect0 = preadvectFilename + '0';
  string compressedPreadvect1 = preadvectFilename + '1';
  string compressedPreadvect2 = preadvectFilename + '2';
  auto preadvectSize0 = FileSize(compressedPreadvect0);
  auto preadvectSize1 = FileSize(compressedPreadvect1);
  auto preadvectSize2 = FileSize(compressedPreadvect2);

  string compressedFinal0 = finalFilename + '0';
  string compressedFinal1 = finalFilename + '1';
  string compressedFinal2 = finalFilename + '2';
  auto finalSize0 = FileSize(compressedFinal0);
  auto finalSize1 = FileSize(compressedFinal1);
  auto finalSize2 = FileSize(compressedFinal2);

  double preadvectCompression = preadvectSize / (double) (preadvectSize0 + preadvectSize1 + preadvectSize2);
  double finalCompression = finalSize / (double) (finalSize0 + finalSize1 + finalSize2);
  double overallCompression = 0.5 * (preadvectCompression + finalCompression);

  printf("U.preadvect compression ratio is %f : 1\n", preadvectCompression);
  printf("U.final compression ratio is %f : 1\n", finalCompression);
  printf("Overall compression ratio is %f: 1\n", overallCompression);
  return overallCompression;
}

// automatically update the compression and movie paths based on compression ratio. calls a python script inside ./cfg
void UpdateCfgFile(int roundedOverallCompression)
{
  string cmd = string("python ./cfg/findReplace.py ") + cfgFilename + string(" ") + to_string(roundedOverallCompression);
  system(cmd.c_str());
}
