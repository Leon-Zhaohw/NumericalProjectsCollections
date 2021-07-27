#include "INTEGER_FIELD_3D.h"
#include "VECTOR.h"
#include "FIELD_3D.h"
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;
using std::vector;

INTEGER_FIELD_3D::INTEGER_FIELD_3D() : 
  _xRes(-1), _yRes(-1), _zRes(-1), _data(NULL)
{
}


INTEGER_FIELD_3D::INTEGER_FIELD_3D(const int& xRes, const int& yRes, const int& zRes) :
  _xRes(xRes), _yRes(yRes), _zRes(zRes)
{
  _totalCells = _xRes * _yRes * _zRes;
 
  try {
    _data = new int[_totalCells];
  }
  catch(std::bad_alloc& exc)
  {
    cout << " Failed to allocate " << _xRes << " " << _yRes << " " << _zRes << " INTEGER_FIELD_3D!" << endl;
    int bytes = _totalCells * sizeof(int);
    cout << (double) bytes / pow(2.0, 20.0) << " MB needed" << endl;
    exit(0);
  }

  for (int x = 0; x < _totalCells; x++)
    _data[x] = 0;
}

INTEGER_FIELD_3D::INTEGER_FIELD_3D(const int* data, const int& xRes, const int& yRes, const int& zRes) :
  _xRes(xRes), _yRes(yRes), _zRes(zRes)
{
  _totalCells = _xRes * _yRes * _zRes;
  try {
    _data = new int[_totalCells];
  }
  catch(std::bad_alloc& exc)
  {
    cout << " Failed to allocate " << _xRes << " " << _yRes << " " << _zRes << " INTEGER_FIELD_3D!" << endl;
    int bytes = _totalCells * sizeof(int);
    cout << (double) bytes / pow(2.0, 20.0) << " MB needed" << endl;
    exit(0);
  }

  for (int x = 0; x < _totalCells; x++) {
    _data[x] = data[x];
  }
}
INTEGER_FIELD_3D::INTEGER_FIELD_3D(const INTEGER_FIELD_3D& m) :
  _xRes(m.xRes()), _yRes(m.yRes()), _zRes(m.zRes()) 
{
  _totalCells = _xRes * _yRes * _zRes;
  
  try {
    _data = new int[_totalCells];
  }
  catch(std::bad_alloc& exc)
  {
    cout << " Failed to allocate " << _xRes << " " << _yRes << " " << _zRes << " INTEGER_FIELD_3D!" << endl;
    double bytes = _xRes * _yRes * _zRes * sizeof(int);
    cout <<  bytes / pow(2.0,20.0) << " MB needed" << endl;
    exit(0);
  }


  for (int x = 0; x < _totalCells; x++)
    _data[x] = m[x];

}
 
INTEGER_FIELD_3D::~INTEGER_FIELD_3D() {
  if (_data) {
    delete[] _data;
    _data = NULL;
  }
}

// Return a flattened VECTOR
VECTOR INTEGER_FIELD_3D::flattened() const {
  VECTOR final(_totalCells);
  int index = 0;
  for (int z = 0; z < _zRes; z++) {
    for (int y = 0; y < _yRes; y++) {
      for (int x = 0; x < _xRes; x++, index++) {
        final[index] = (*this)(x, y, z);
      }
    }
  }
  return final;
}

// Return a flattened VECTOR in row order
VECTOR INTEGER_FIELD_3D::flattenedRow() const {
  VECTOR final(_totalCells);
  int index = 0;
  for (int x = 0; x < _xRes; x++) {
    for (int y = 0; y < _yRes; y++) {
      for (int z = 0; z < _zRes; z++, index++) {
          final[index] = (*this)(x, y, z);
      }
    }
  }
  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


void INTEGER_FIELD_3D::clear()
{
  //for (int x = 0; x < _totalCells; x++)
  //  _data[x] = 0;

  memset(_data, 0, sizeof(int) * _totalCells);
}

void INTEGER_FIELD_3D::clear(const vector<int>& nonZeros)
{
  for (unsigned int x = 0; x < nonZeros.size(); x++)
    _data[nonZeros[x]] = 0;
}

void INTEGER_FIELD_3D::clear(const vector<int>& nonZeros, const int size)
{
  //TIMER functionTimer("INTEGER_FIELD_3D clear");
  for (unsigned int x = 0; x < size; x++)
    _data[nonZeros[x]] = 0;
}

void INTEGER_FIELD_3D::resizeAndWipe(int xRes, int yRes, int zRes)
{
  if (_xRes == xRes && _yRes == yRes && _zRes == zRes)
  {
    clear();
    return;
  }

  if (_data) delete[] _data;

  _xRes = xRes;
  _yRes = yRes;
  _zRes = zRes;
  _totalCells = _xRes * _yRes * _zRes;

  try {
    _data = new int[_totalCells];
  }
  catch(std::bad_alloc& exc)
  {
    cout << " Failed to allocate " << _xRes << " " << _yRes << " " << _zRes << " INTEGER_FIELD_3D!" << endl;
    double bytes = _xRes * _yRes * _zRes * sizeof(int);
    cout <<  bytes / pow(2.0,20.0) << " MB needed" << endl;
    exit(0);
  }

  for (int x = 0; x < _totalCells; x++)
    _data[x] = 0;
}

int INTEGER_FIELD_3D::energySum() {
  int sum = 0;
  for (int i = 0; i < _totalCells; i++) {
    sum += (_data[i] * _data[i]);
  }
  return sum;
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
// Overloaded operators
INTEGER_FIELD_3D& INTEGER_FIELD_3D::operator*=(const double& alpha)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] *= alpha;

  return *this;
}

// Get a subfield
INTEGER_FIELD_3D INTEGER_FIELD_3D::subfield(const int xBegin, const int xEnd, 
                            const int yBegin, const int yEnd,
                            const int zBegin, const int zEnd) const {
  assert(xBegin >= 0);
  assert(yBegin >= 0);
  assert(zBegin >= 0);
  assert(xEnd <= _xRes);
  assert(yEnd <= _yRes);
  assert(zEnd <= _zRes);
  assert(xBegin < xEnd);
  assert(yBegin < yEnd);
  assert(zBegin < zEnd);

  int xInterval = xEnd - xBegin;
  int yInterval = yEnd - yBegin;
  int zInterval = zEnd - zBegin;

  INTEGER_FIELD_3D final(xInterval, yInterval, zInterval);

  /* cout << "size: " << final.xRes() << ", " << final.yRes() << ", " << 
     final.zRes() << endl; */
  for (int z = 0; z < zInterval; z++) {
    for (int y = 0; y < yInterval; y++) {
      for (int x = 0; x < xInterval; x++) {
        final(x, y, z) = (*this)(xBegin + x, yBegin + y, zBegin + z);
      }
    }
  }
  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
INTEGER_FIELD_3D& INTEGER_FIELD_3D::operator=(const INTEGER_FIELD_3D& A)
{
  resizeAndWipe(A.xRes(), A.yRes(), A.zRes());

  for (int x = 0; x < _totalCells; x++)
    _data[x] = A[x];

  return *this;
}

// Swap contents with annother 3d integer field
void INTEGER_FIELD_3D::swapPointers(INTEGER_FIELD_3D& field)
{
  assert(field.xRes() == _xRes);
  assert(field.yRes() == _yRes);
  assert(field.zRes() == _zRes);
  
  int* temp = _data;
  _data = field._data;
  field._data = temp;
}

