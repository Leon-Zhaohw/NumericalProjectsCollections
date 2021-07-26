/*

  PAT dipole sources evaluation
  This is supplemental material to our SIGGRAPH 2006 paper on real-time sound simulation.
  The code evaluates the given dipole sources at the specified listening locations, summing the
  contributions of all the dipoles. The output are the absolute values of the sound pressure, at
  the specified listening locations, per each frequency.
  
  The dipole data is provided for three models: bell, dragon, chair. Also, the .obj mesh for
  each model is provided.

  The main() routine is a driver that takes in the microphone locations and a text file listing
  the dipole datafiles per each frequency, and produces the pressure file. 
  See 'example' file for an example command-line.

  Jernej Barbic
  Carnegie Mellon University
  2006

  See the paper:

  Doug James, Jernej Barbic, Dinesh Pai:
  Precomputed Acoustic Transfer:
  Output-sensitive, accurate sound generation for
  geometrically complex vibration sources
  ACM Transactions on Graphics (SIGGRAPH 2006),
  Boston, MA, August 2006

  Note: this code can be accelerated by using vector performance libraries,
  such as Intel MKL. Efficient hardware implementations are most likely possible too.

*/

/*
  Note on matrices used in this code:

  All matrices loaded/written to disk are stored in the following **binary** format:
  <number of rows> (signed integer, 4 bytes)
  <number of columns> (signed integer, 4 bytes)
  data in double precision (8 bytes), in COLUMN-major order (first column first, then second column, and so on)

  So, a m x n matrix will occupy 2*sizeof(int) + m*n*sizeof(double) bytes on disk. 

  The routines in matrixIO.h and all the dipole datafiles conform to this format.
  Likewise, the input microphone position file should be given in this format.
  You can use the matrixIO routines to generate input microphone position files in this format. 
  Also, you can use matrixIO to read the .sources and .k files (as also demonstrated in this code in the main() routine).
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include "matrixIO.h"
#include "performanceCounter.h"

int numFrequencies; // total number of frequencies
int T; // total number of distinct locations where PAT pressures need to be evaluated
       // T is determined by reading the header of the input file specifying microphone positions
double * pressureAbs; // output pressure buffer
double * microphonePos; // the 3xT 2D array holding 3D locations where PAT pressures need to be evaluated
std::vector<double*> patObjects; // vector of pointers to dipole data, one pointer per frequency
std::vector<double> kTable; // vector of wave numbers, k=2*pi/lambda, one per frequency
std::vector<int> numSources; // num sources per each frequency
int maxNumSources; // max num sources over all frequencies
double * buffer; // workspace
PerformanceCounter performanceCounter; // to time PAT evaluation
// note: these global variables could be wrapped in a class

// evaluates dipoles for T microphone positions from the array microphonePos,
// using data from patObjects, kTable, numSources, numFrequencies
// stores computed pressures into a numFrequencies x T matrix 'pressureAbs'
void EvaluateDipoles()
{
  // each .sources file contains the data for all the dipoles for that frequency
  // the .sources file is an 11 x (num sources) matrix
  // each dipole is described by 11 coefficients, double precision each
  // each dipole corresponds to one column of the .sources matrix
  // there will be one .sources and .k file for every frequency

  // structure of data for each individual dipole (11 entries):
  // first 3 coefficients: world-coordinate dipole position (x,y,z)
  // then, monopole term (complex number, first Re, then Im)
  // then m=0 dipole term (complex number, first Re, then Im)
  // then m=-1 dipole term (complex number, first Re, then Im)
  // then m=+1 dipole term (complex number, first Re, then Im)

  // note: the 0,-1,+1 order might seem counter-intuitive, but it
  // naturally generalizes to higher-order sources: 0,-1,+1,-2,+2,-3,+3, etc.

  const int ssize = 11;

  // over all frequencies:
  for(int frequencyIndex = 0; frequencyIndex < numFrequencies; frequencyIndex++)
  {
    printf("%d ",frequencyIndex+1);fflush(stdout);
    double k = kTable[frequencyIndex]; // get the wave number for this frequency
    double pressureRe, pressureIm;
    double * sourceData = patObjects[frequencyIndex]; // point to the sources for this frequency

    // over all microphone positions (i.e. over all time samples):
    for(int t=0; t<T; t++)
    {
      double * pos = &microphonePos[3*t]; // get microphone position at this time-step
      pressureRe = 0;
      pressureIm = 0;

      #define COMPLEX_MULTIPLY(aRe,aIm,bRe,bIm,oRe,oIm)\
        oRe = (aRe) * (bRe) - (aIm) * (bIm);\
        oIm = (aRe) * (bIm) + (aIm) * (bRe);

      #define COMPLEX_MULTIPLY_ADD(aRe,aIm,bRe,bIm,oRe,oIm)\
        oRe += (aRe) * (bRe) - (aIm) * (bIm);\
        oIm += (aRe) * (bIm) + (aIm) * (bRe);

      double * coef = sourceData;
      // over all dipole sources at this frequency:
      for(int source=0; source < numSources[frequencyIndex]; source++)
      {
        // relative position of microphone with respect to the source
        double x = pos[0] - coef[0];
        double y = pos[1] - coef[1];
        double z = pos[2] - coef[2];

        // auxiliary quantities
        double planarR2 = x*x + y*y;
        double planarR = sqrt(planarR2);
        double r = sqrt(planarR2 + z*z);
        double cosTheta = z / r;
        double sinTheta = planarR / r;
        double cosPhi = x / planarR;
        double sinPhi = y / planarR;

        double kr = k * r;

        double invKr = 1.0 / kr;
        double sinKr = sin(kr);
        double cosKr = cos(kr);

        // monopole term:

        double bufferRe = sinKr * invKr;
        double bufferIm = cosKr * invKr;

        COMPLEX_MULTIPLY_ADD(bufferRe,bufferIm,coef[3],coef[4],pressureRe,pressureIm);

        // dipole terms:

        double radialRe = invKr * (-cosKr + invKr * sinKr);
        double radialIm = invKr * ( sinKr + invKr * cosKr);

        // m = 0
        bufferRe = radialRe * cosTheta;
        bufferIm = radialIm * cosTheta;
        COMPLEX_MULTIPLY_ADD(bufferRe,bufferIm,coef[5],coef[6],pressureRe,pressureIm);

        double cosPhiSinTheta = cosPhi * sinTheta;
        double sinPhiSinTheta = sinPhi * sinTheta;

        // m = -1
        COMPLEX_MULTIPLY(radialRe,radialIm,cosPhiSinTheta,-sinPhiSinTheta,bufferRe,bufferIm);
        COMPLEX_MULTIPLY_ADD(bufferRe,bufferIm,coef[7],coef[8],pressureRe,pressureIm);

        // m = 1
        COMPLEX_MULTIPLY(radialRe,radialIm,-cosPhiSinTheta,-sinPhiSinTheta,bufferRe,bufferIm);
        COMPLEX_MULTIPLY_ADD(bufferRe,bufferIm,coef[9],coef[10],pressureRe,pressureIm);

        coef += ssize; // increment coef to point to the next source
      }

      // store |p| :
      pressureAbs[numFrequencies * t + frequencyIndex] = sqrt(pressureRe*pressureRe + pressureIm*pressureIm);

    }
  }
  printf("\n");
}

// the driver routine

int main(int argc, char** argv)
{
  const int numFixedArg = 4;
  if ( argc < numFixedArg )
  {
    printf("Evaluates PAT functions at the given locations.\n");
    printf("Usage: %s [3xT matrix of microphone positions] [text file listing PAT dipole files] [output pressure file]\n",argv[0]);
    return 1;
  }

  char * posFile = argv[1];  // microphone positions file
  char * listFile = argv[2]; // a text file listing the dipole files, one per each frequency
    // the dipole files are listed without the .sources or .k prefix, one frequency per line
  char * outputFile = argv[3]; // output file to store evaluated pressures

  // load microphone positions
  int m;
  ReadMatrixFromDisk_(posFile,&m,&T,&microphonePos);
  Assert_(m,3,0);

  // read and load PAT dipole files, one file per frequency
  FILE * fin;
  OpenFile_(listFile,&fin,"ra");
  char patEntry[4096];
  char s[4096];
  int totalNumSources = 0;
  maxNumSources = 0;
  while (fgets(patEntry,4096,fin) != NULL)
  {
    patEntry[strlen(patEntry)-1] = 0; // kill newline at the end

    // patEntry now contains the name of the file containing
    // dipole information for the next frequency
    // load data from that file into main memory
    double * patObject;
    int m1, n1;
    sprintf(s,"%s.sources",patEntry);
    printf("%s\n",s);
    ReadMatrixFromDisk_(s,&m1,&n1,&patObject);

    Assert_(m1,11,1); // matrix must have 11 rows, 1 is just to identify assertion location
    patObjects.push_back(patObject);
    numSources.push_back(n1);

    // on-the-fly computation of the max number of dipoles over all frequencies
    // (necessary to allocate the appropriate size workspace buffer)
    if (n1 > maxNumSources)
      maxNumSources = n1;

    totalNumSources += n1;

    // read the wave numbers k = 2pi/lambda
    sprintf(s,"%s.k",patEntry);
    double * kk;
    ReadMatrixFromDisk_(s,&m1,&n1,&kk);
    Assert_(m1,1,2); // the .k matrix must have exactly 1 row
    Assert_(n1,1,3); // the .k matrix must have exactly 1 column; note: 2,3 are error indices, just to be able to locate the assertion location
    kTable.push_back(*kk);
    free(kk);
  }
  fclose(fin);

  numFrequencies = patObjects.size();

  // allocate space for pressure output (numFrequencies x T matrix)
  pressureAbs = (double*) calloc (numFrequencies * T ,sizeof(double));

  // allocate workspace buffer
  buffer = (double*) calloc (4 * maxNumSources ,sizeof(double));

  printf("Total num sources: %d\n", totalNumSources);
  printf("Evaluating %d PAT models at %d positions...\n",numFrequencies,T);

  // evaluate the dipoles
  performanceCounter.StartCounter();
  EvaluateDipoles();
  performanceCounter.StopCounter();

  double evaluationTime = performanceCounter.GetElapsedTime();
  printf("Total evaluation time: %G .\n", evaluationTime);
  printf("Evaluation time per source per position: %G msec .\n", 1000 * evaluationTime / totalNumSources / T);

  // write pressures to output file
  WriteMatrixToDisk_(outputFile,numFrequencies,T,pressureAbs);

  // cleanup
  for(int i=0; i<numFrequencies; i++)
    free(patObjects[i]);

  free(pressureAbs);
  free(microphonePos);

  free(buffer);

  return 0;
}
