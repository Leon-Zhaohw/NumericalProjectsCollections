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

//////////////////////////////////////////////////////////////////////
// read matrix from a file
//////////////////////////////////////////////////////////////////////
bool EIGEN::read(const string& filename, MatrixXd& input)
{
  string copy = filename;

  // make sure it's names gz
  int size = copy.size();
  if (copy[size - 1] == 'z' || copy[size - 2] == 'g')
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Don't send EIGEN::read a gz file! " << endl;
    exit(0);
  }

  FILE* file;
  file = fopen(copy.c_str(), "rb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __LINE__ << " : File " << copy << " not found! " << endl;
    return false;
  }
  cout << " Reading file: " << copy.c_str() << " ... ";flush(cout);

  // read dimensions
  int rows, cols;
  fread((void*)&rows, sizeof(int), 1, file);
  fread((void*)&cols, sizeof(int), 1, file);
  input.resize(rows, cols);
  input.setZero();

  cout << " rows: " << rows << " cols: " << cols << " " << flush;

  int tenth = rows / 10;

  // fread appears to fail after 2 GB?
  // need to read entry-by-entry for this reason
  double dummy;
  int index = 0;
  for (int x = 0; x < rows; x++)
  {
    for (int y = 0; y < cols; y++, index++)
    {
      fread((void*)&dummy, sizeof(double), 1, file);
      input(x,y) = dummy;
    }

    if (rows > 10000000 && x > 0 && x % tenth == 0)
    {
      cout << " Purging at " << x << " of " << rows << " rows... " << flush;
      system("purge");
      cout << " done." << flush;
    }
  }

  fclose(file);
  cout << "done." << endl;

  return true;
}

//////////////////////////////////////////////////////////////////////
// read bigmatrix from a file
//////////////////////////////////////////////////////////////////////
bool EIGEN::readBig(const string& filename, MatrixXd& input)
{
  // open the file
  FILE* file;
  file = fopen(filename.c_str(), "rb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __LINE__ << " : File " << filename.c_str() << " not found! " << endl;
    return false;
  }
  cout << " Reading file " << filename.c_str() << " ... ";flush(cout);

  // read dimensions
  int rows, cols;
  fread((void*)&rows, sizeof(int), 1, file);
  fread((void*)&cols, sizeof(int), 1, file); 
  input.resize(rows, cols);
  input.setZero();

  cout << "Num rows: " << rows << endl;
  cout << "Num cols: " << cols << endl;

  int tenth = rows / 10;
 
  // read in entry-by-entry, but every (numRow + 1)th entry is numRows, so skip those
  double dummy;
  int index = 0;
  for (int x = 0; x < rows; x++)
  {
    for (int y = 0; y < cols; y++, index++)
    {
      if (index % (rows + 1) == 0) {
        fseek(file, sizeof(int), SEEK_CUR);  
      }

      fread((void*)&dummy, sizeof(double), 1, file);
      input(x,y) = dummy;
    }

    if (rows > 10000000 && x > 0 && x % tenth == 0)
    {
      cout << " Purging at " << x << " of " << rows << " rows... " << flush;
      system("purge");
      cout << " done." << flush;
    }
  }
  
  fclose(file);
  cout << " done." << endl;
  return true;
}
//////////////////////////////////////////////////////////////////////
// read matrix from a file
//////////////////////////////////////////////////////////////////////
void EIGEN::readRaw(FILE* file, const int size, VectorXd& input)
{
  input.resize(size);

  // fread appears to fail after 2 GB?
  // need to read entry-by-entry for this reason
  double dummy;
  int index = 0;
  for (int x = 0; x < size; x++)
  {
    fread((void*)&dummy, sizeof(double), 1, file);
    input[x] = dummy;
  }
}

//////////////////////////////////////////////////////////////////////
// read matrix from a file
//////////////////////////////////////////////////////////////////////
void EIGEN::read(FILE* file, MatrixXd& input)
{
  // read dimensions
  int rows, cols;
  fread((void*)&rows, sizeof(int), 1, file);
  fread((void*)&cols, sizeof(int), 1, file);
  input.resize(rows, cols);

  // fread appears to fail after 2 GB?
  // need to read entry-by-entry for this reason
  double dummy;
  int index = 0;
  for (int x = 0; x < rows; x++)
    for (int y = 0; y < cols; y++, index++)
    {
      fread((void*)&dummy, sizeof(double), 1, file);
      input(x,y) = dummy;
    }
}

//////////////////////////////////////////////////////////////////////
// read vector from passed in c-string filename
//////////////////////////////////////////////////////////////////////
void EIGEN::read(const char* filename, VectorXd& input)
{
  FILE* file;
  file = fopen(filename, "rb");
  if (file == NULL) 
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " VECTOR read failed! " << endl;
    cout << " Could not open file " << filename << endl;
    exit(EXIT_FAILURE);
  }
  // read dimension
  int size;
  fread((void*)&size, sizeof(int), 1, file);
  input.resize(size);
  
  double dummy;
  for (int x = 0; x < size; x++) {
    fread((void*)&dummy, sizeof(double), 1, file);
    input[x] = dummy;
  }
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// write the matrix to a file
//////////////////////////////////////////////////////////////////////
bool EIGEN::write(const string& filename, const MatrixXd& input)
{
  string copy = filename;

  // make sure it's not named gz
  int size = copy.size();
  if (copy[size - 1] == 'z' || copy[size - 2] == 'g')
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Don't try to write gz file from EIGEN::write!!!" << endl;
    exit(0);
  }

  FILE* file;
  file = fopen(copy.c_str(), "wb"); 
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " MATRIX write failed! " << endl;
    cout << " Could not open file " << copy.c_str() << endl;
    return false;
  }
  cout << " Writing file: " << copy.c_str() << " ... "; flush(cout);

  // write dimensions
  int rows = input.rows();
  int cols = input.cols();
  fwrite((void*)&rows, sizeof(int), 1, file);
  fwrite((void*)&cols, sizeof(int), 1, file);

  // just do the dumb but memory stingy way
  for (int x = 0; x < rows; x++)
    for (int y = 0; y < cols; y++)
    {
      double dummy = input(x,y);
      fwrite((void*)&dummy, sizeof(double), 1, file);
    }
  fclose(file);

  cout << " done. " << endl;

  return true;
}

//////////////////////////////////////////////////////////////////////
// write the gzipped matrix to a file
//////////////////////////////////////////////////////////////////////
void EIGEN::write(FILE* file, const MatrixXd& input)
{
  TIMER functionTimer(__FUNCTION__);
  // write dimensions
  int rows = input.rows();
  int cols = input.cols();
  fwrite((void*)&rows, sizeof(int), 1, file);
  fwrite((void*)&cols, sizeof(int), 1, file);

  // just do the dumb but memory stingy way
  for (int x = 0; x < rows; x++)
    for (int y = 0; y < cols; y++)
    {
      double dummy = input(x,y);
      fwrite((void*)&dummy, sizeof(double), 1, file);
    }
}

//////////////////////////////////////////////////////////////////////
// write the gzipped matrix to a file
//////////////////////////////////////////////////////////////////////
void EIGEN::write(FILE* file, const VectorXd& input)
{
  // write dimensions
  int size = input.size();
  fwrite((void*)&size, sizeof(int), 1, file);

  // just do the dumb but memory stingy way
  for (int x = 0; x < size; x++)
  {
    double dummy = input[x];
    fwrite((void*)&dummy, sizeof(double), 1, file);
  }
}

//////////////////////////////////////////////////////////////////////
// write matrix to a file using passed in C string
//////////////////////////////////////////////////////////////////////
void EIGEN::write(const char* filename, const MatrixXd& input)
{
  TIMER functionTimer(__FUNCTION__);
  FILE* file;
  file = fopen(filename, "wb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " MATRIX write failed! " << endl;
    cout << " Could not open file " << filename << endl; 
    exit(EXIT_FAILURE); 
  }
  cout << " Writing file: " << filename << " ... "; flush(cout);
  // write dimensions
  int rows = input.rows();
  int cols = input.cols();
  fwrite((void*)&rows, sizeof(int), 1, file);
  fwrite((void*)&cols, sizeof(int), 1, file);

  // just do the dumb but memory stingy way
  for (int x = 0; x < rows; x++) {
    for (int y = 0; y < cols; y++) {
      double dummy = input(x,y);
      fwrite((void*)&dummy, sizeof(double), 1, file);
    }
  }
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// write vector to a file using passed in C string
//////////////////////////////////////////////////////////////////////
void EIGEN::write(const char* filename, const VectorXd& input)
{
  TIMER functionTimer(__FUNCTION__);
  FILE* file;
  file = fopen(filename, "wb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " VECTOR write failed! " << endl;
    cout << " Could not open file " << filename << endl; 
    exit(EXIT_FAILURE); 
  }
  cout << " Writing file: " << filename << " ... "; flush(cout);
  // write dimensions
  int size = input.size();
  fwrite((void*)&size, sizeof(int), 1, file);

  // just do the dumb but memory stingy way
  for (int x = 0; x < size; x++) {
    double dummy = input[x];
    fwrite((void*)&dummy, sizeof(double), 1, file);
  }
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// Build a matrix from an array of columns
//////////////////////////////////////////////////////////////////////
MatrixXd EIGEN::buildFromColumns(const vector<VectorXd>& columns)
{
  assert(columns.size() > 0);
  int rows = columns[0].size();
  int cols = columns.size();

  MatrixXd final(columns[0].size(), columns.size());

  for (int y = 0; y < cols; y++)
    for (int x = 0; x < rows; x++)
      final(x,y) = columns[y][x];

  return final;
}

//////////////////////////////////////////////////////////////////////
// Get a submatrix of the current matrix
//////////////////////////////////////////////////////////////////////
MatrixXd EIGEN::getRows(const int rowBegin, const int totalRows, const MatrixXd& input)
{
  const int rows = input.rows();
  const int cols = input.cols();

  assert(rowBegin < rows);
  assert(rowBegin >= 0);
  assert(totalRows > 0);
  assert(rowBegin + totalRows <= rows);

  MatrixXd toReturn(totalRows, cols);

  for (int x = 0; x < totalRows; x++)
    for (int y = 0; y < cols; y++)
      toReturn(x,y) = input(x + rowBegin, y);

  return toReturn;
}


//////////////////////////////////////////////////////////////////////
// Get a submatrix from reading a large file
//////////////////////////////////////////////////////////////////////
void EIGEN::getRowsMemory(const int rowBegin, const int totalRows, FILE* matrixInfo, double* matrixData, int rows, int cols)
{
  // check bounds

  assert(rowBegin < rows);
  assert(rowBegin >= 0);
  assert(totalRows > 0);
  assert(rowBegin + totalRows <= rows);

  // ADJ: we're only ever using it for 3 at a time
  assert(totalRows == 3);

  // get the offset from the rows, cols at BOF
  auto beginningOffset = 2 * sizeof(int);

  // DEBUG
  auto offset = beginningOffset + rowBegin*cols*sizeof(double);


  // seek from BOF (SEEK_SET)
  rewind(matrixInfo);
  auto checkError = fseek(matrixInfo, beginningOffset + rowBegin*cols*sizeof(double), SEEK_SET);
  if (checkError!=0) {perror("fseek failed!"); exit(EXIT_FAILURE);};
  int numElements = totalRows * cols;
  auto numRead = fread(matrixData, sizeof(double), numElements, matrixInfo);
  if (numRead != numElements) {fputs ("Reading error",stderr); exit (3);}
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void EIGEN::transposeProduct(const MatrixXd& left, const MatrixXd& right, MatrixXd& output)
{
  output.resize(left.cols(), right.cols());

#pragma omp parallel
#pragma omp for  schedule(static)
  for (int r = 0; r < right.cols(); r++)
    for (int l = 0; l < left.cols(); l++)
      output(l,r) = left.col(l).dot(right.col(r));
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
VECTOR EIGEN::convert(const VectorXd& input)
{
  VECTOR final(input.size());

  for (int x = 0; x < input.size(); x++)
    final[x] = input[x];

  return final;
}
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

VECTOR EIGEN::convertInt(const VectorXi& input)
{
  VECTOR final(input.size());

  for (int x = 0; x < input.size(); x++)
    final[x] = input[x];

  return final;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
VectorXd EIGEN::convert(const VECTOR& input)
{
  VectorXd final(input.size());

  for (int x = 0; x < input.size(); x++)
    final[x] = input[x];

  return final;
}
