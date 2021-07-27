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
#include "BIG_MATRIX.h"
#include <iomanip>

// IMPORTANT! Must be set based on the system's memory capacity.
// There will be *two* blocks in memory at a time, so set accordingly
int BIG_MATRIX::_blockSize = 115;

// IMPORTANT! Set this path to the SSD drive
string BIG_MATRIX::_scratchPath = string("./scratch/");

//////////////////////////////////////////////////////////////////////
// Constructor for the full matrix
//////////////////////////////////////////////////////////////////////
BIG_MATRIX::BIG_MATRIX() :
  _rows(-1), _cols(-1)
{
}

BIG_MATRIX::BIG_MATRIX(int rows, int cols) :
  _rows(rows), _cols(cols)
{
  // allocate placeholders for the columns, but don't allocate them
  _columns.resize(_cols);
}

BIG_MATRIX::BIG_MATRIX(const BIG_MATRIX& m) 
{
  _cols = m._cols;
  _rows = m._rows;

  for (int x = 0; x < _cols; x++)
    _columns.push_back(m._columns[x]);
}

BIG_MATRIX::BIG_MATRIX(const MATRIX& m)
{
  _cols = m.cols();
  _rows = m.rows();

  for (int x = 0; x < _cols; x++)
    _columns.push_back(VECTOR(_rows));

  for (int x = 0; x < _cols; x++)
    for (int y = 0; y < _rows; y++)
      (*this)(y,x) = m(y,x);
}

BIG_MATRIX::BIG_MATRIX(const string& filename)
{
  read(filename);
}

//////////////////////////////////////////////////////////////////////
// Print matrix to stream, assuming the matrix is fully in core
//////////////////////////////////////////////////////////////////////
ostream& operator<<(ostream &out, const BIG_MATRIX& matrix)
{
  int oldPrecision = out.precision();

  out << "[" << endl;
  for (int row = 0; row < matrix.rows(); row++)
  {
    for (int col = 0; col < matrix.cols(); col++)
    {
      //out << matrix(row, col) << " ";
      out.width(MATRIX::coutWidth());
      out << setprecision(MATRIX::coutPrecision()) << matrix(row, col) << " ";
    }
    out << endl;
  }
  out << "]" << endl;
  out << resetiosflags(ios::floatfield);
  out.precision(oldPrecision);
  return out;
}

//////////////////////////////////////////////////////////////////////
// do an out-of-core Gram-Schmidt QR factorization
//////////////////////////////////////////////////////////////////////
MATRIX BIG_MATRIX::outOfCoreQR(const string& filenamePrefix, int& qRows, int& qCols)
{
  TIMER functionTimer(__FUNCTION__);
  string filename = filenamePrefix + string(".transpose");

  // make sure it's not compressed
  int size = filename.size();
  if (filename[size - 1] == 'z' && filename[size - 2] == 'g')
  {
    cout << " Don't send the out-of-core factorer a compressed file! " << endl;
    exit(0);
  }

  // read dimensions file
  FILE* file;
  string dimsFilename = filenamePrefix + string(".dims");
  file = fopen(dimsFilename.c_str(), "rb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __LINE__ << " : File " << dimsFilename.c_str() << " not found! " << endl;
    exit(0);
  }
  int rows, cols;
  fread((void*)&rows, sizeof(int), 1, file);
  fread((void*)&cols, sizeof(int), 1, file);
  fclose(file);

  // open the snapshots file
  file = fopen(filename.c_str(), "rb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __LINE__ << " : File " << filename.c_str() << " not found! " << endl;
    exit(0);
  }
  cout << "======================================================================" << endl;
  cout << " Out-of-Core QR factoring file " << filename.c_str() << endl;
  cout << "======================================================================" << endl;

  // swap them, since we're reading in a transpose file
  int temp = rows;
  rows = cols;
  cols = temp;

  // store them to be returned
  qRows = rows;
  qCols = cols;

  cout << " Big matrix dims: " << rows << " " << cols << endl;
  _blockSize = computeBlockSize(rows, cols);

  MATRIX R(cols, cols);
  int totalBlocks = cols / _blockSize;
  if (cols % _blockSize != 0)
    totalBlocks++;

  cout << " Snapshot matrix dimensions are (" << rows << ", " << cols << ")" << endl;
  cout << " There are " << totalBlocks << " blocks to process " << endl;

  // write the dimensions of the Q file to the scratch space
  writeDimensions("Q.dims", qRows, qCols); 

  for (int b = 0; b < totalBlocks; b++)
  {
    // figure out which columns we're looking at
    int currentBegin = b * _blockSize;
    int currentEnd = (b + 1) * _blockSize;
    currentEnd = (currentEnd > cols) ? cols : currentEnd;

    // if we're on the last block, scale back the number of blocks to
    // read if needed
    int totalReadColumns = currentEnd - currentBegin;

    // read in the current columns
    BIG_MATRIX currentColumns;
    system("purge");

    // indent the reads for output readability
    TIMER currentRead("QR current read");
    cout << " Reading block " << b << " ... " << flush;
    currentColumns.readColumns(rows, totalReadColumns, file);
    cout << " done. " << endl;
    currentRead.stop();
    system("purge");

    // subtract off columns from previous blocks
    for (int p = 0; p < b; p++)
    {
      // read in the block
      TIMER previousRead("QR previous read");
      char buffer[256];
      sprintf(buffer, "%i", p);
      string previousFilename = _scratchPath + string(buffer) + string(".block");
      // indent the reads for output readability
      cout << "  ";
      BIG_MATRIX previousColumns(previousFilename);
      previousRead.stop();

      TIMER subtractTimer("QR previous subtract");
      int previousBegin = p * _blockSize;

      // subtract off the upstream column
      cout << "     Subtracting off previous block " << p << " ... " << flush;
      for (int i = 0; i < previousColumns.cols(); i++)
      {
        const VECTOR& previous = previousColumns[i];

        // for each column in the current block
        // NOTE: this can be done in parallel
#pragma omp parallel
#pragma omp for schedule(static)
        for (int j = 0; j < currentColumns.cols(); j++)
        {
          Real dot = previous * currentColumns[j];
          R(previousBegin + i, currentBegin + j) = dot;

          //currentColumns[j].axpySSE(-dot, previous);
          currentColumns[j].axpy(-dot, previous);
        }
      }
      cout << " done. " << endl;
      subtractTimer.stop();
      system("purge");
    }

    // perform QR on the current block
    TIMER factorTimer("QR current factor");
    cout << " Factoring block " << b << " ... " << flush;
    int totalColumns = currentColumns.cols();
    for (int i = 0; i < totalColumns; i++)
    {
      cout << " column " << i << " of " << totalColumns << flush;

      // get the 2-norm of the column
      Real norm = currentColumns[i] * currentColumns[i];
      norm = sqrt(norm);
      
      // populate the R diagonal
      R(currentBegin + i, currentBegin + i) = norm;

      // normalize the column
      Real invNorm = 1.0 / norm;
      currentColumns[i] *= invNorm;

      // subtract from everybody downstream
      // note we can do this in parallel
      const VECTOR& column = currentColumns[i];
#pragma omp parallel
#pragma omp for schedule(static)
      for (int j = i + 1; j < currentColumns.cols(); j++)
      {
        Real dot = column * currentColumns[j];
        R(currentBegin + i, currentBegin + j) = dot;
        currentColumns[j].axpy(-dot, column);
      }
    }
    cout << " done." << endl;
    factorTimer.stop();

    // write out the block
    TIMER writeTimer("QR current write");
    char buffer[256];
    sprintf(buffer, "%i", b);
    string currentBlockFilename = _scratchPath + string(buffer) + string(".block");
    currentColumns.write(currentBlockFilename);
    writeTimer.stop();

    TIMER::printTimings();
  }
  system("purge");

  return R;
}

//////////////////////////////////////////////////////////////////////
// do an out-of-core SVD using the ooc QR
//////////////////////////////////////////////////////////////////////
void BIG_MATRIX::outOfCoreSVD(const string& filenamePrefix, const string& reducedPath, const Real& discardThreshold)
{
  TIMER functionTimer(__FUNCTION__);

  // do the out of core QR, assume the block files are in the
  // scratch space
  int qRows, qCols;
  MATRIX R = outOfCoreQR(filenamePrefix, qRows, qCols);

  // do the SVD on R (can replace this with Matlab later)
  TIMER svdTimer("SVD call");
#if __APPLE__
  MATRIX U;
  MATRIX VT;
  VECTOR S;

  R.SVD(U, S, VT);
#else
  MATRIX U;
  MATRIX VT;
  VECTOR S;
  MatrixXd eigenR(R.rows(), R.cols());
  for (int y = 0; y < R.rows(); y++)
    for (int x = 0; x < R.cols(); x++)
      eigenR(x,y) = R(x,y);

  JacobiSVD<MatrixXd> svd(eigenR, ComputeThinU | ComputeThinV);

  S.resizeAndWipe(svd.singularValues().size());
  for (int x = 0; x < svd.singularValues().size(); x++)
    S[x] = svd.singularValues()[x];

  U.resizeAndWipe(svd.matrixU().rows(), svd.matrixU().cols());
  for (int y = 0; y < svd.matrixU().cols(); y++)
    for (int x = 0; x < svd.matrixU().rows(); x++)
      U(x,y) = svd.matrixU()(x,y);
  
  VT.resizeAndWipe(svd.matrixV().rows(), svd.matrixV().cols());
  for (int y = 0; y < svd.matrixV().cols(); y++)
    for (int x = 0; x < svd.matrixV().rows(); x++)
      VT(x,y) = svd.matrixV()(x,y);
#endif

  // actually, give back PCA values
  for (int x =0; x < S.size(); x++)
    S[x] = S[x] * S[x];

  svdTimer.stop();

  // figure out how many columns to keep
  int keepingColumns = 0;
  for (int x = 0; x < S.size(); x++)
  {
    if (S[x] / S[0] < discardThreshold)
      break;
    keepingColumns++;
  }
  cout << " Keeping " << keepingColumns << " of " << S.size() << " columns " << endl;
  cout << "=======================================" << endl;
  cout << " Building U from SVD of R" << endl;
  cout << "=======================================" << endl;

  // how many blocks are there in Q?
  int blocksQ = qCols / _blockSize;
  if (qCols % _blockSize != 0)
    blocksQ++;

  // do the out-of-core multiply -- build a block at a time
  int blocksU = keepingColumns / _blockSize;
  int lastBlockSize = _blockSize;
  if (keepingColumns % _blockSize != 0)
  {
    blocksU++;
    lastBlockSize = keepingColumns % _blockSize;
  }
 
  // write out the dims of the U so we can recombine it later 
  writeDimensions("U.dims", qRows, keepingColumns); 

  // get a version of the SVD's U (not the full version) with easily accessible columns
  BIG_MATRIX bigU(U);

  // go through each block in the U
  for (int u = 0; u < blocksU; u++)
  {
    // get the columns of U to multiply this time
    int uBegin = u * _blockSize;
    int uEnd   = (u + 1) * _blockSize;
    uEnd = (uEnd < keepingColumns) ? uEnd : keepingColumns;
    int uTotalCols = uEnd - uBegin;

    // build the current matrix block to multiply into
    BIG_MATRIX blockFinalU(qRows, uTotalCols);

    // force U to allocate all the memory
    for (int x = 0; x < uTotalCols; x++)
      blockFinalU[x].resizeAndWipe(qRows);

    // pull in each block of Q from disk to form the final U
    for (int q = 0; q < blocksQ; q++)
    {
      // read in the Q block
      TIMER qRead("SVD reading in Q");
      char buffer[256];
      sprintf(buffer, "%i", q);
      string qBlockFilename = _scratchPath + string(buffer) + string(".block");
      BIG_MATRIX qColumns;
      // indent the reads for output readability
      cout << "  ";
      qColumns.read(qBlockFilename);
      qRead.stop();

      // which columns of Q are we looking at?
      int qBegin = q * _blockSize;
      int qEnd = (q + 1) * _blockSize;
      qEnd = (qEnd < qCols) ? qEnd : qCols;
      int qTotalCols = qEnd - qBegin;

      // for each column in the just pulled in Q block
      TIMER uMultiply("SVD Q times U");
      cout << "    Multiplying block " << q << " of Q to form U, " << qTotalCols << " columns total " << endl;
      cout << "      Computing column ... " << flush;
      for (int qCol = 0; qCol < qTotalCols; qCol++)
      {
        cout << qCol << " " << flush;

        // use it to add the appropriate U column
        // should be able to parallelize this?
#pragma omp parallel
#pragma omp for  schedule(static)
        for (int uCol = 0; uCol < uTotalCols; uCol++)
        {
          // extract the scalar from u
          Real uEntry = bigU[uBegin + uCol][qBegin + qCol];
          
          // scale the column of Q and add it to the final U
          blockFinalU[uCol].axpy(uEntry, qColumns[qCol]);
        }
      }
      cout << " done. " << endl;

      uMultiply.stop();
    }

    // write the new U to disk
    TIMER uWrite("SVD U write");
    char buffer[256];
    sprintf(buffer, "%i", u);
    string uBlockFilename = _scratchPath + string(buffer) + string(".U.block");
    blockFinalU.write(uBlockFilename);
    uWrite.stop();
    
    TIMER::printTimings();
  }
}

//////////////////////////////////////////////////////////////////////
// read in columns from a stream
//////////////////////////////////////////////////////////////////////
bool BIG_MATRIX::readColumns(const int rows, const int totalColumns, FILE* file)
{
  // set dims
  _rows = rows;
  _cols = totalColumns;
  _columns.clear();
  _columns.resize(_cols);

  for (int x = 0; x < totalColumns; x++)
  {
    _columns[x].resizeAndWipe(rows);
    fread((void*)(_columns[x].data()), sizeof(double), rows, file);
  }

  return true;
}

//////////////////////////////////////////////////////////////////////
// write the matrix to a binary file
//////////////////////////////////////////////////////////////////////
void BIG_MATRIX::write(const string& filename)
{
  // open the file
  FILE* file;
  file = fopen(filename.c_str(), "wb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : Couldn't open file " << filename.c_str() << "! " << endl;
    exit(0);
  }
  cout << " Writing file " << filename.c_str() << " ... ";flush(cout);

  // write dimensions
  fwrite((void*)&_rows, sizeof(int), 1, file);
  fwrite((void*)&_cols, sizeof(int), 1, file);

  // dump out the columns 
  for (int x = 0; x < _cols; x++)
    _columns[x].write(file);

  fclose(file);

  cout << " done." << endl;
}

//////////////////////////////////////////////////////////////////////
// read the matrix from a binary file
//////////////////////////////////////////////////////////////////////
void BIG_MATRIX::read(const string& filename)
{
  // open the file
  FILE* file;
  file = fopen(filename.c_str(), "rb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __LINE__ << " : File " << filename.c_str() << " not found! " << endl;
    exit(0);
  }
  cout << " Reading file " << filename.c_str() << " ... ";flush(cout);

  // read dimensions
  fread((void*)&_rows, sizeof(int), 1, file);
  fread((void*)&_cols, sizeof(int), 1, file);
  
  _columns.resize(_cols);

  // read in the columns 
  for (int x = 0; x < _cols; x++)
    _columns[x].read(file);

  fclose(file);

  cout << " done." << endl;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
BIG_MATRIX& BIG_MATRIX::operator=(const BIG_MATRIX m)
{
  _rows = m._rows;
  _cols = m._cols;

  _columns.clear();
  for (int x = 0; x < _cols; x++)
    _columns.push_back(m._columns[x]);
  return *this;
}

//////////////////////////////////////////////////////////////////////
// write a dimensions file to the scratch space
//////////////////////////////////////////////////////////////////////
void BIG_MATRIX::writeDimensions(const string filename, const int& rows, const int& cols)
{
  string dimsFilename = _scratchPath + filename;
  FILE* file = fopen(dimsFilename.c_str(), "wb");
  
  if (file == NULL)
  {
    cout << " Couldn't open dimensions file " << dimsFilename.c_str() << " for writing" << endl;
    exit(0);
  }  

  // write dimensions
  fwrite((void*)&rows, sizeof(int), 1, file);
  fwrite((void*)&cols, sizeof(int), 1, file);
  fwrite((void*)&_blockSize, sizeof(int), 1, file);

  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// read write a dimensions file to the scratch space
//////////////////////////////////////////////////////////////////////
void BIG_MATRIX::readDimensions(const string filename, int& rows, int& cols)
{
  string dimsFilename = _scratchPath + filename;
  FILE* file = fopen(dimsFilename.c_str(), "rb");
  
  if (file == NULL)
  {
    cout << " Couldn't open dimensions file " << dimsFilename.c_str() << " for reading" << endl;
    exit(0);
  }

  // read dimensions
  fread((void*)&rows, sizeof(int), 1, file);
  fread((void*)&cols, sizeof(int), 1, file);
  fread((void*)&_blockSize, sizeof(int), 1, file);

  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// compute a feasbile block size
//////////////////////////////////////////////////////////////////////
int BIG_MATRIX::computeBlockSize(int rows, int cols)
{
  double columnBytes = rows * 8;
  
  // this is hard coded to kraken/mako
  double systemMemory =  92.0 * pow(10.0, 9.0);

  double totalColumns = systemMemory / columnBytes;
  cout << " Total columns that fit in memory: " << (int)totalColumns << endl;

  int half = totalColumns / 2.0;

  // round it off to the nearest 10
  int rounded = (half - half % 10);

  cout << " Block size, rounded to nearest 10: " << rounded << endl;

  return rounded;
}

//////////////////////////////////////////////////////////////////////
// write out the final big matrix
//////////////////////////////////////////////////////////////////////
void BIG_MATRIX::writeFinalU(const string& filename)
{
  TIMER functionTimer(__FUNCTION__);
  cout << " ========================================== " << endl;
  cout << "  Writing out final U from " << _scratchPath.c_str() << endl;
  cout << " ========================================== " << endl;

  int uRows, uCols;
  readDimensions("U.dims", uRows, uCols);

  // how many blocks are there in Q?
  int blocksU = uCols / _blockSize;
  if (uCols % _blockSize != 0)
    blocksU++;

  vector<BIG_MATRIX> matrices(blocksU);
  for (int u = 0; u < blocksU; u++)
  {
    // read in the Q block
    char buffer[256];
    sprintf(buffer, "%i", u);
    string uBlockFilename = _scratchPath + string(buffer) + string(".U.block");
    matrices[u].read(uBlockFilename);
  }

  cout << " Final U of size " << uRows << " " << uCols << endl;

  cout << " Writing out final file " << filename.c_str() << " ... " << flush;
  FILE* file;
  file = fopen(filename.c_str(), "wb");

  fwrite((void*)&uRows, sizeof(int), 1, file);
  fwrite((void*)&uCols, sizeof(int), 1, file);

  for (int y = 0; y < uRows; y++)
  {
    double row[uCols];
    int i = 0;
    for (unsigned int u = 0; u < matrices.size(); u++)
      for (int x = 0; x < matrices[u].cols(); x++, i++)
        row[i] = matrices[u](y, x);
    fwrite((void*)row, sizeof(double), uCols, file);
  }
  cout << endl;

  fclose(file);
}
