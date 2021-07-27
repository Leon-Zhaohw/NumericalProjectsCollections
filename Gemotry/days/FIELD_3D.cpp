/*
QUIJIBO: Source code for the paper Symposium on Geometry Processing
         2015 paper "Quaternion Julia Set Shape Optimization"
Copyright (C) 2015  Theodore Kim

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "FIELD_3D.h"
#include <zlib.h>
#include "TIMER.h"

int FIELD_3D::_quinticClamps = 0;

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D::FIELD_3D(const int& xRes, const int& yRes, const int& zRes,
    const VEC3F& center, const VEC3F& lengths) :
  _xRes(xRes), _yRes(yRes), _zRes(zRes), _center(center), _lengths(lengths), _rotation(1,0,0,0)
{
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
 
  try {
    _data = new Real[_totalCells];
  }
  catch(std::bad_alloc& exc)
  {
    cout << " Failed to allocate " << _xRes << " " << _yRes << " " << _zRes << " FIELD_3D!" << endl;
    double bytes = _xRes * _yRes * _zRes * sizeof(Real);
    cout <<  bytes / pow(2.0,20.0) << " MB needed" << endl;
    exit(0);
  }

  for (int x = 0; x < _totalCells; x++)
    _data[x] = 0;

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
  _invDz = 1.0 / _dz;

  _outside = maxRes() * maxRes();
}

FIELD_3D::FIELD_3D(const Real* data, const int& xRes, const int& yRes, const int& zRes,
    const VEC3F& center, const VEC3F& lengths) :
  _xRes(xRes), _yRes(yRes), _zRes(zRes), _center(center), _lengths(lengths), _rotation(1,0,0,0)
{
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  try {
    _data = new Real[_totalCells];
  }
  catch(std::bad_alloc& exc)
  {
    cout << " Failed to allocate " << _xRes << " " << _yRes << " " << _zRes << " FIELD_3D!" << endl;
    double bytes = _xRes * _yRes * _zRes * sizeof(Real);
    cout <<  bytes / pow(2.0,20.0) << " MB needed" << endl;
    exit(0);
  }

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
  _invDz = 1.0 / _dz;

  _outside = maxRes() * maxRes();

  for (int x = 0; x < _totalCells; x++)
    _data[x] = data[x];
}

FIELD_3D::FIELD_3D(const bool* data, const int& xRes, const int& yRes, const int& zRes) :
  _xRes(xRes), _yRes(yRes), _zRes(zRes), _center(0,0,0), _lengths(1,1,1), _rotation(1,0,0,0)
{
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  try {
    _data = new Real[_totalCells];
  }
  catch(std::bad_alloc& exc)
  {
    cout << " Failed to allocate " << _xRes << " " << _yRes << " " << _zRes << " FIELD_3D!" << endl;
    double bytes = _xRes * _yRes * _zRes * sizeof(Real);
    cout <<  bytes / pow(2.0,20.0) << " MB needed" << endl;
    exit(0);
  }

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
  _invDz = 1.0 / _dz;

  _outside = maxRes() * maxRes();

  for (int x = 0; x < _totalCells; x++)
    _data[x] = data[x];
}

FIELD_3D::FIELD_3D(const unsigned char* data, const int& xRes, const int& yRes, const int& zRes) :
  _xRes(xRes), _yRes(yRes), _zRes(zRes), _center(0,0,0), _lengths(1,1,1), _rotation(1,0,0,0)
{
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  try {
    _data = new Real[_totalCells];
  }
  catch(std::bad_alloc& exc)
  {
    cout << " Failed to allocate " << _xRes << " " << _yRes << " " << _zRes << " FIELD_3D!" << endl;
    double bytes = _xRes * _yRes * _zRes * sizeof(Real);
    cout <<  bytes / pow(2.0,20.0) << " MB needed" << endl;
    exit(0);
  }

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
  _invDz = 1.0 / _dz;

  _outside = maxRes() * maxRes();

  for (int x = 0; x < _totalCells; x++)
    _data[x] = data[x];
}

FIELD_3D::FIELD_3D(const FIELD_3D& m) :
  _xRes(m.xRes()), _yRes(m.yRes()), _zRes(m.zRes()),
  _center(m.center()), _lengths(m.lengths()), _rotation(m.rotation())
{
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  try {
    _data = new Real[_totalCells];
  }
  catch(std::bad_alloc& exc)
  {
    cout << " Failed to allocate " << _xRes << " " << _yRes << " " << _zRes << " FIELD_3D!" << endl;
    double bytes = _xRes * _yRes * _zRes * sizeof(Real);
    cout <<  bytes / pow(2.0,20.0) << " MB needed" << endl;
    exit(0);
  }

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
  _invDz = 1.0 / _dz;

  for (int x = 0; x < _totalCells; x++)
    _data[x] = m[x];

  _outside = maxRes() * maxRes();
}

FIELD_3D::FIELD_3D() :
  _xRes(-1), _yRes(-1), _zRes(-1), _totalCells(-1), _data(NULL), _rotation(1,0,0,0)
{
}

FIELD_3D::FIELD_3D(const char* filename) :
  _xRes(-1), _yRes(-1), _zRes(-1), _totalCells(-1), _data(NULL)
{
  cout << " Reading file " << filename << endl;
  int size = string(filename).size();
  if (filename[size - 1] == 'z' && filename[size - 2] == 'g')
    readGz(filename);
  else
    read(filename);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D::~FIELD_3D()
{
  if (_data)
    delete[] _data;
}
  
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_3D::clear()
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] = 0.0;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_3D::write(string filename) const
{
  FILE* file;
  file = fopen(filename.c_str(), "wb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " FIELD_3D write failed! " << endl;
    cout << " Could not open file " << filename.c_str() << endl;
    exit(0);
  }

  cout << " Writing file " << filename.c_str() << " ... "; flush(cout);

  // write to the stream
  write(file);

  // close the stream
  fclose(file);

  cout << " done. " << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_3D::writeGz(string filename) const
{
  TIMER functionTimer(__FUNCTION__);
  // make sure it's named gz
  int size = filename.size();
  if (filename[size - 1] != 'z' || filename[size - 2] != 'g')
    filename = filename + string(".gz");

  gzFile file;
  file = gzopen(filename.c_str(), "wb1"); 
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " FIELD_3D write failed! " << endl;
    cout << " Could not open file " << filename.c_str() << endl;
    exit(0);
  }
  cout << " Writing file " << filename.c_str() << " ... "; flush(cout);

  writeGz(file);

  gzclose(file);

  cout << " done. " << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_3D::read(string filename)
{
  int size = filename.size();
  if (filename[size - 1] == 'z' && filename[size - 2] == 'g')
  {
    readGz(filename);
    return;
  }

  FILE* file;
  file = fopen(filename.c_str(), "rb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " FIELD_3D read failed! " << endl;
    cout << " Could not open file " << filename.c_str() << endl;
    exit(0);
  }

  // read from the stream
  read(file);

  // close the file
  fclose(file);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_3D::readGz(string filename)
{
  // make sure it's names gz
  int size = filename.size();
  if (filename[size - 1] != 'z' || filename[size - 2] != 'g')
    filename = filename + string(".gz");

  gzFile file;
  file = gzopen(filename.c_str(), "rb1");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " FIELD_3D read failed! " << endl;
    cout << " Could not open file " << filename.c_str() << endl;
    exit(0);
  }

  cout << " Reading file " << filename.c_str() << " ... ";flush(cout);

  readGz(file);
  gzclose(file);

  cout << " done." << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_3D::resizeAndWipe(int xRes, int yRes, int zRes, const VEC3F& center, const VEC3F& lengths)
{
  if (_xRes == xRes && _yRes == yRes && _zRes == zRes)
  {
    _center = center;
    _lengths = lengths;
    clear();

    _dx = _lengths[0] / _xRes;
    _dy = _lengths[1] / _yRes;
    _dz = _lengths[2] / _zRes;
    _invDx = 1.0 / _dx;
    _invDy = 1.0 / _dy;
    _invDz = 1.0 / _dz;
    return;
  }

  if (_data) delete[] _data;

  _xRes = xRes;
  _yRes = yRes;
  _zRes = zRes;
  _slabSize = _xRes * _yRes;
  _totalCells = _xRes * _yRes * _zRes;
  _outside = maxRes() * maxRes();
  _center = center;
  _lengths = lengths;

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
  _invDz = 1.0 / _dz;

  try {
    _data = new Real[_totalCells];
  }
  catch(std::bad_alloc& exc)
  {
    cout << " Failed to allocate " << _xRes << " " << _yRes << " " << _zRes << " FIELD_3D!" << endl;
    double bytes = _xRes * _yRes * _zRes * sizeof(Real);
    cout <<  bytes / pow(2.0,20.0) << " MB needed" << endl;
    exit(0);
  }

  for (int x = 0; x < _totalCells; x++)
    _data[x] = 0;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator=(const Real& alpha)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] = alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator*=(const Real& alpha)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] *= alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator/=(const Real& alpha)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] /= alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator+=(const Real& alpha)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] += alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator-=(const Real& alpha)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] -= alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator-=(const FIELD_3D& input)
{
  assert(input.xRes() == _xRes);
  assert(input.yRes() == _yRes);
  assert(input.zRes() == _zRes);
  for (int x = 0; x < _totalCells; x++)
    _data[x] -= input[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator+=(const FIELD_3D& input)
{
  assert(input.xRes() == _xRes);
  assert(input.yRes() == _yRes);
  assert(input.zRes() == _zRes);
  for (int x = 0; x < _totalCells; x++)
    _data[x] += input[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator*=(const FIELD_3D& input)
{
  assert(input.xRes() == _xRes);
  assert(input.yRes() == _yRes);
  assert(input.zRes() == _zRes);

  for (int x = 0; x < _totalCells; x++)
    _data[x] *= input[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator^(const FIELD_3D& A, const Real alpha)
{
  FIELD_3D final(A);

  for (int x = 0; x < final.totalCells(); x++)
    final[x] = pow(final[x], alpha);

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator*(const FIELD_3D& A, const Real alpha)
{
  FIELD_3D final(A);
  final *= alpha;
  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator/(const FIELD_3D& A, const FIELD_3D& B)
{
  FIELD_3D final(A);

  for (int x = 0; x < A.totalCells(); x++)
    final[x] *= 1.0 / B[x];

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator/(const FIELD_3D& A, const Real alpha)
{
  FIELD_3D final(A);
  final /= alpha;
  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator-(const FIELD_3D& A, const FIELD_3D& B)
{
  assert(A.xRes() == B.xRes());
  assert(A.yRes() == B.yRes());
  assert(A.zRes() == B.zRes());

  FIELD_3D final(A);
  final -= B;
  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator+(const FIELD_3D& A, const FIELD_3D& B)
{
  FIELD_3D final(A);
  final += B;
  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator*(const FIELD_3D& A, const FIELD_3D& B)
{
  FIELD_3D final(A);
  final *= B;
  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator+(const FIELD_3D& A, const Real alpha)
{
  FIELD_3D final(A);
  final += alpha;
  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator*(const Real alpha, const FIELD_3D& A)
{
  return A * alpha;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator+(const Real alpha, const FIELD_3D& A)
{
  return A + alpha;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator-(const Real alpha, const FIELD_3D& A)
{
  FIELD_3D final(A);

  for (int x = 0; x < final.totalCells(); x++)
    final[x] = alpha - final[x];

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator=(const FIELD_3D& A)
{
  if (A.totalCells() < 0)
  {
    _xRes = -1;
    _yRes = -1;
    _zRes = -1;
    _totalCells = -1;

    if (_data) delete[] _data;
    _data = NULL;
  }

  resizeAndWipe(A.xRes(), A.yRes(), A.zRes(), A.center(), A.lengths());

  for (int x = 0; x < _totalCells; x++)
    _data[x] = A[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
// real-valued cell center coordinates
///////////////////////////////////////////////////////////////////////
VEC3F FIELD_3D::cellCenter(int index) const
{
  int z = index / _slabSize;
  index -= z * _slabSize;
  int y = index / _xRes;
  int x = index - y * _xRes;

  return cellCenter(x,y,z);
}

///////////////////////////////////////////////////////////////////////
// real-valued cell center coordinates
///////////////////////////////////////////////////////////////////////
VEC3F FIELD_3D::cellCenter(int x, int y, int z) const
{
  VEC3F halfLengths = (Real)0.5 * _lengths;

  // set it to the lower corner
  VEC3F final = _center - halfLengths;

  // displace to the NNN corner
  final[0] += x * _dx;
  final[1] += y * _dy;
  final[2] += z * _dz;

  // displace it to the cell center
  final[0] += _dx * 0.5;
  final[1] += _dy * 0.5;
  final[2] += _dz * 0.5;

  // apply the rotation
  final -= _center;
  final = _rotation.toRotationMatrix() * final;
  final += _center;

  return final;
}

///////////////////////////////////////////////////////////////////////
// do a union with 'field', assuming both this and field are signed i
// distance fields
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::signedDistanceUnion(const FIELD_3D& field)
{
  assert(field.xRes() == _xRes);
  assert(field.yRes() == _yRes);
  assert(field.zRes() == _zRes);

  FIELD_3D final = field;

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        if ((*this)(x,y,z) > 0.0 && field(x,y,z) > 0.0)
          final(x,y,z) = (*this)(x,y,z) < field(x,y,z) ? (*this)(x,y,z) : field(x,y,z);
        else if ((*this)(x,y,z) < 0.0 && field(x,y,z) < 0.0)
          final(x,y,z) = (*this)(x,y,z) < field(x,y,z) ? (*this)(x,y,z) : field(x,y,z);
        else if ((*this)(x,y,z) < 0.0 && field(x,y,z) > 0.0)
          final(x,y,z) = (*this)(x,y,z);
        else 
          final(x,y,z) = field(x,y,z);
      }

  return final;
}

///////////////////////////////////////////////////////////////////////
// what's the maximum resolution in any direction?
///////////////////////////////////////////////////////////////////////
int FIELD_3D::maxRes()
{
  int final = _xRes;
  if (_yRes > final) final = _yRes;
  if (_zRes > final) final = _zRes;
  return final;
}


///////////////////////////////////////////////////////////////////////
// return the indices of the grid points insides a world-space bounding box
///////////////////////////////////////////////////////////////////////
void FIELD_3D::boundingBoxIndices(const VEC3F& mins, const VEC3F& maxs, VEC3I& iMins, VEC3I& iMaxs)
{
  //VEC3F halfLengths = (Real)0.5 * _lengths;

  VEC3F ds(_dx, _dy, _dz);

  VEC3F fracMins = (mins - (ds * (Real)0.5) - _center + (_lengths * (Real)0.5));
  VEC3F fracMaxs = (maxs - (ds * (Real)0.5) - _center + (_lengths * (Real)0.5));

  // fatten things to account for the integer cast
  iMins[0] = fracMins[0] / _dx - 2;
  iMins[1] = fracMins[1] / _dy - 2;
  iMins[2] = fracMins[2] / _dz - 2;

  iMaxs[0] = fracMaxs[0] / _dx + 2;
  iMaxs[1] = fracMaxs[1] / _dy + 2;
  iMaxs[2] = fracMaxs[2] / _dz + 2;

  // clamp them
  for (int x = 0; x < 3; x++)
  {
    if (iMins[x] < 0) iMins[x] = 0;
    if (iMaxs[x] < 0) iMaxs[x] = 0;
  }
  if (iMins[0] > _xRes - 1) iMins[0] = _xRes - 1;
  if (iMins[1] > _yRes - 1) iMins[1] = _yRes - 1;
  if (iMins[2] > _zRes - 1) iMins[2] = _zRes - 1;

  if (iMaxs[0] > _xRes - 1) iMaxs[0] = _xRes - 1;
  if (iMaxs[1] > _yRes - 1) iMaxs[1] = _yRes - 1;
  if (iMaxs[2] > _zRes - 1) iMaxs[2] = _zRes - 1;
}

///////////////////////////////////////////////////////////////////////
// copy out the boundary
///////////////////////////////////////////////////////////////////////
void FIELD_3D::copyBorderAll()
{
  int index;
  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
    {
      // front slab
      index = x + y * _xRes;
      _data[index] = _data[index + _slabSize];

      // back slab
      index += _totalCells - _slabSize;
      _data[index] = _data[index - _slabSize];
    }

  for (int z = 0; z < _zRes; z++)
    for (int x = 0; x < _xRes; x++)
    {
      // bottom slab
      index = x + z * _slabSize;
      _data[index] = _data[index + _xRes];

      // top slab
      index += _slabSize - _xRes;
      _data[index] = _data[index - _xRes];
    }

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
    {
      // left slab
      index = y * _xRes + z * _slabSize;
      _data[index] = _data[index + 1];

      // right slab
      index += _xRes - 1;
      _data[index] = _data[index - 1];
    }
}

///////////////////////////////////////////////////////////////////////
// lookup value at some real-valued spatial position
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::operator()(const VEC3F& position) const
{
  VEC3F positionCopy = position;

  // get the lower corner position
  VEC3F corner = _center - (Real)0.5 * _lengths;
  VEC3F dxs(_dx, _dy, _dz);
  corner += (Real)0.5 * dxs;

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  int x0 = (int)positionCopy[0];
  int x1    = x0 + 1;
  int y0 = (int)positionCopy[1];
  int y1    = y0 + 1;
  int z0 = (int)positionCopy[2];
  int z1    = z0 + 1;

  // clamp everything
  x0 = (x0 < 0) ? 0 : x0;
  y0 = (y0 < 0) ? 0 : y0;
  z0 = (z0 < 0) ? 0 : z0;
  
  x1 = (x1 < 0) ? 0 : x1;
  y1 = (y1 < 0) ? 0 : y1;
  z1 = (z1 < 0) ? 0 : z1;

  x0 = (x0 > _xRes - 1) ? _xRes - 1 : x0;
  y0 = (y0 > _yRes - 1) ? _yRes - 1 : y0;
  z0 = (z0 > _zRes - 1) ? _zRes - 1 : z0;

  x1 = (x1 > _xRes - 1) ? _xRes - 1 : x1;
  y1 = (y1 > _yRes - 1) ? _yRes - 1 : y1;
  z1 = (z1 > _zRes - 1) ? _zRes - 1 : z1;

  // get interpolation weights
  const float s1 = positionCopy[0]- x0;
  const float s0 = 1.0f - s1;
  const float t1 = positionCopy[1]- y0;
  const float t0 = 1.0f - t1;
  const float u1 = positionCopy[2]- z0;
  const float u0 = 1.0f - u1;

  const int i000 = x0 + y0 * _xRes + z0 * _slabSize;
  const int i010 = x0 + y1 * _xRes + z0 * _slabSize;
  const int i100 = x1 + y0 * _xRes + z0 * _slabSize;
  const int i110 = x1 + y1 * _xRes + z0 * _slabSize;
  const int i001 = x0 + y0 * _xRes + z1 * _slabSize;
  const int i011 = x0 + y1 * _xRes + z1 * _slabSize;
  const int i101 = x1 + y0 * _xRes + z1 * _slabSize;
  const int i111 = x1 + y1 * _xRes + z1 * _slabSize;

  // interpolate
  // (indices could be computed once)
  return u0 * (s0 * (t0 * _data[i000] + t1 * _data[i010]) +
               s1 * (t0 * _data[i100] + t1 * _data[i110])) +
         u1 * (s0 * (t0 * _data[i001] + t1 * _data[i011]) +
               s1 * (t0 * _data[i101] + t1 * _data[i111]));
}

///////////////////////////////////////////////////////////////////////
// summed squared entries
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::sumSq()
{
  Real final = 0;
  for (int i = 0; i < _totalCells; i++)
    final += _data[i] * _data[i];

  return final;
}

///////////////////////////////////////////////////////////////////////
// maximum entry
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::max()
{
  Real final = _data[0];

  for (int i = 0; i < _totalCells; i++)
    if (_data[i] > final)
      final = _data[i];

  return final;
}

///////////////////////////////////////////////////////////////////////
// maximum entry
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::absMax()
{
  Real final = fabs(_data[0]);

  for (int i = 0; i < _totalCells; i++)
    if (fabs(_data[i]) > final)
      final = fabs(_data[i]);

  return final;
}

///////////////////////////////////////////////////////////////////////
// build a const field with the given dims
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::constField(const FIELD_3D& dims, Real value)
{
  FIELD_3D final(dims);
  final = value;
  return final;
}

///////////////////////////////////////////////////////////////////////
// check if any entry is a nan
///////////////////////////////////////////////////////////////////////
bool FIELD_3D::isNan()
{
  int i = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, i++)
#ifdef QUAD_PRECISION
        if (std::isnanq(_data[i]))
#else
        if (std::isnan(_data[i]))
#endif
        {
          cout << " Nan found at: " << x << ", " << y << ", " << z << endl;
          return true;
        }

  return false;
}

///////////////////////////////////////////////////////////////////////
// compute the inverse
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::inverse()
{
  FIELD_3D final(*this);

  for (int x = 0; x < _totalCells; x++)
    final[x] = 1.0 / _data[x];

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_3D::printNeighborhood(int index) const
{
  int z = index / _slabSize;
  int y = (index - z * _slabSize) / _xRes;
  int x = (index - z * _slabSize - y * _xRes);

  printNeighborhood(x,y,z);
}

///////////////////////////////////////////////////////////////////////
// print the neighborhood of a cell for debugging
///////////////////////////////////////////////////////////////////////
void FIELD_3D::printNeighborhood(int x, int y, int z) const
{
  cout << " Neighborhood of (" << x << ", " << y << ", " << z << ")" << endl;
  for (int k = z - 1; k <= z + 1; k++)
  {
    cout.precision(8);
    cout << "[ " << endl;
    for (int j = y - 1; j <= y + 1; j++)
    {
      for (int i = x - 1; i <= x + 1; i++)
      {
        int index = i + j * _xRes + k * _slabSize;
        cout << _data[index] << " ";
      }
      cout << ";" << endl;
    }
    cout << "]" << endl;
  }
}

///////////////////////////////////////////////////////////////////////
// clamp nans to some specified value
///////////////////////////////////////////////////////////////////////
void FIELD_3D::clampNans(Real value)
{
  for (int x = 0; x < _totalCells; x++)
#ifdef QUAD_PRECISION
    if (std::isnanq(_data[x]))
#else
    if (std::isnan(_data[x]))
#endif
      _data[x] = value;
}

///////////////////////////////////////////////////////////////////////
// clamp nans to some specified value
///////////////////////////////////////////////////////////////////////
void FIELD_3D::clampInfs(Real value)
{
  for (int x = 0; x < _totalCells; x++)
#ifdef QUAD_PRECISION
    if (std::isinfq(_data[x]))
#else
    if (std::isinf(_data[x]))
#endif
      _data[x] = value;
}

///////////////////////////////////////////////////////////////////////
// clamp infinities to values in this field
///////////////////////////////////////////////////////////////////////
void FIELD_3D::clampInfs(FIELD_3D& clampField)
{
  for (int x = 0; x < _totalCells; x++)
#ifdef QUAD_PRECISION
    if (std::isinfq(_data[x]))
#else
    if (std::isinf(_data[x]))
#endif
      _data[x] = clampField[x];
}

///////////////////////////////////////////////////////////////////////
// clamp nans to values in this field
///////////////////////////////////////////////////////////////////////
void FIELD_3D::clampNans(FIELD_3D& clampField)
{
  for (int x = 0; x < _totalCells; x++)
#ifdef QUAD_PRECISION
    if (std::isnanq(_data[x]))
#else
    if (std::isnan(_data[x]))
#endif
      _data[x] = clampField[x];
}

///////////////////////////////////////////////////////////////////////
// cell index of a real-valued position
///////////////////////////////////////////////////////////////////////
void FIELD_3D::cellIndex(const VEC3F& position, VEC3I& indices) const
{
  VEC3F positionCopy = position;

  // get the lower corner position
  VEC3F corner = _center - (Real)0.5 * _lengths;
  VEC3F dxs(_dx, _dy, _dz);
  corner += (Real)0.5 * dxs;

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  int x0 = (int)positionCopy[0];
  int x1    = x0 + 1;
  int y0 = (int)positionCopy[1];
  int y1    = y0 + 1;
  int z0 = (int)positionCopy[2];
  int z1    = z0 + 1;

  // clamp everything
  x0 = (x0 < 0) ? 0 : x0;
  y0 = (y0 < 0) ? 0 : y0;
  z0 = (z0 < 0) ? 0 : z0;
  
  x1 = (x1 < 0) ? 0 : x1;
  y1 = (y1 < 0) ? 0 : y1;
  z1 = (z1 < 0) ? 0 : z1;

  x0 = (x0 > _xRes - 1) ? _xRes - 1 : x0;
  y0 = (y0 > _yRes - 1) ? _yRes - 1 : y0;
  z0 = (z0 > _zRes - 1) ? _zRes - 1 : z0;

  x1 = (x1 > _xRes - 1) ? _xRes - 1 : x1;
  y1 = (y1 > _yRes - 1) ? _yRes - 1 : y1;
  z1 = (z1 > _zRes - 1) ? _zRes - 1 : z1;

  // get interpolation weights
  const float s1 = positionCopy[0]- x0;
  const float s0 = 1.0f - s1;
  const float t1 = positionCopy[1]- y0;
  const float t0 = 1.0f - t1;
  const float u1 = positionCopy[2]- z0;
  const float u0 = 1.0f - u1;

  indices[0] = (s0 > s1) ? x0 : x1;
  indices[1] = (t0 > t1) ? y0 : y1;
  indices[2] = (u0 > u1) ? z0 : z1;
}

///////////////////////////////////////////////////////////////////////
// cell index of a real-valued position
///////////////////////////////////////////////////////////////////////
void FIELD_3D::cellNeighborhood(const VEC3F& position, vector<VEC3I>& indices, vector<Real>& weights) const
{
  indices.clear();
  weights.clear();

  VEC3F positionCopy = position;

  // get the lower corner position
  VEC3F corner = _center - (Real)0.5 * _lengths;
  VEC3F dxs(_dx, _dy, _dz);
  corner += (Real)0.5 * dxs;

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  int x0 = (int)positionCopy[0];
  int x1    = x0 + 1;
  int y0 = (int)positionCopy[1];
  int y1    = y0 + 1;
  int z0 = (int)positionCopy[2];
  int z1    = z0 + 1;

  // clamp everything
  x0 = (x0 < 0) ? 0 : x0;
  y0 = (y0 < 0) ? 0 : y0;
  z0 = (z0 < 0) ? 0 : z0;
  
  x1 = (x1 < 0) ? 0 : x1;
  y1 = (y1 < 0) ? 0 : y1;
  z1 = (z1 < 0) ? 0 : z1;

  x0 = (x0 > _xRes - 1) ? _xRes - 1 : x0;
  y0 = (y0 > _yRes - 1) ? _yRes - 1 : y0;
  z0 = (z0 > _zRes - 1) ? _zRes - 1 : z0;

  x1 = (x1 > _xRes - 1) ? _xRes - 1 : x1;
  y1 = (y1 > _yRes - 1) ? _yRes - 1 : y1;
  z1 = (z1 > _zRes - 1) ? _zRes - 1 : z1;

  // get interpolation weights
  const float s1 = positionCopy[0]- x0;
  const float s0 = 1.0f - s1;
  const float t1 = positionCopy[1]- y0;
  const float t0 = 1.0f - t1;
  const float u1 = positionCopy[2]- z0;
  const float u0 = 1.0f - u1;

  indices.push_back(VEC3I(x0, y0, z0));
  indices.push_back(VEC3I(x0, y0, z1));
  indices.push_back(VEC3I(x0, y1, z0));
  indices.push_back(VEC3I(x0, y1, z1));
  indices.push_back(VEC3I(x1, y0, z0));
  indices.push_back(VEC3I(x1, y0, z1));
  indices.push_back(VEC3I(x1, y1, z0));
  indices.push_back(VEC3I(x1, y1, z1));

  weights.push_back(s0 * t0 * u0);
  weights.push_back(s0 * t0 * u1);
  weights.push_back(s0 * t1 * u0);
  weights.push_back(s0 * t1 * u1);
  weights.push_back(s1 * t0 * u0);
  weights.push_back(s1 * t0 * u1);
  weights.push_back(s1 * t1 * u0);
  weights.push_back(s1 * t1 * u1);
}

///////////////////////////////////////////////////////////////////////
// cell index of a real-valued position
///////////////////////////////////////////////////////////////////////
int FIELD_3D::cellIndex(VEC3F& position) const
{
  VEC3F positionCopy = position;

  // get the lower corner position
  VEC3F corner = _center - (Real)0.5 * _lengths;
  VEC3F dxs(_dx, _dy, _dz);
  corner += (Real)0.5 * dxs;

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  int x0 = (int)positionCopy[0];
  int x1    = x0 + 1;
  int y0 = (int)positionCopy[1];
  int y1    = y0 + 1;
  int z0 = (int)positionCopy[2];
  int z1    = z0 + 1;

  // clamp everything
  x0 = (x0 < 0) ? 0 : x0;
  y0 = (y0 < 0) ? 0 : y0;
  z0 = (z0 < 0) ? 0 : z0;
  
  x1 = (x1 < 0) ? 0 : x1;
  y1 = (y1 < 0) ? 0 : y1;
  z1 = (z1 < 0) ? 0 : z1;

  x0 = (x0 > _xRes - 1) ? _xRes - 1 : x0;
  y0 = (y0 > _yRes - 1) ? _yRes - 1 : y0;
  z0 = (z0 > _zRes - 1) ? _zRes - 1 : z0;

  x1 = (x1 > _xRes - 1) ? _xRes - 1 : x1;
  y1 = (y1 > _yRes - 1) ? _yRes - 1 : y1;
  z1 = (z1 > _zRes - 1) ? _zRes - 1 : z1;

  // get interpolation weights
  const float s1 = positionCopy[0]- x0;
  const float s0 = 1.0f - s1;
  const float t1 = positionCopy[1]- y0;
  const float t0 = 1.0f - t1;
  const float u1 = positionCopy[2]- z0;
  const float u0 = 1.0f - u1;

  int xFinal = (s0 > s1) ? x0 : x1;
  int yFinal = (t0 > t1) ? y0 : y1;
  int zFinal = (u0 > u1) ? z0 : z1;

  return xFinal + yFinal * _xRes + zFinal * _slabSize;
}

///////////////////////////////////////////////////////////////////////
// pass a field to fieldViewer3D
///////////////////////////////////////////////////////////////////////
void FIELD_3D::fieldViewer(const FIELD_3D& field, string name)
{
  field.write("./temp/temp3d.field");
  string execute("./bin/fieldViewer3D ./temp/temp3d.field \"");
  execute = execute + name + string("\" &");
  cout << " Executing " << execute.c_str() << endl;
  system(execute.c_str());
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_3D::volumeRender(const FIELD_3D& field, const string name)
{
  field.write("./temp/temp3d.field");
  string execute("./bin/fieldVolumeRender3D ./temp/temp3d.field \"");
  execute = execute + name + string("\" &");
  cout << " Executing " << execute.c_str() << endl;
  system(execute.c_str());
}

///////////////////////////////////////////////////////////////////////
// pass a field to fieldViewer3D
///////////////////////////////////////////////////////////////////////
void FIELD_3D::fieldViewerYZ(const FIELD_3D& field, string name)
{
  field.write("./temp/temp3d.field");
  string execute("./bin/fieldViewer3D ./temp/temp3d.field \"");
  execute = execute + name + string("\" yz &");
  cout << " Executing " << execute.c_str() << endl;
  system(execute.c_str());
}

///////////////////////////////////////////////////////////////////////
// pass a field with a distance field (inside/outside) overlay
// to a overlayFieldViewer3D
///////////////////////////////////////////////////////////////////////
void FIELD_3D::overlayFieldViewer(const FIELD_3D& field, const FIELD_3D& distance, string name)
{
  field.write("./temp/temp3d.field");
  distance.write("./temp/temp3d.distance.field");
  string execute("./bin/overlayFieldViewer3D ./temp/temp3d.field ./temp/temp3d.distance.field \"");
  execute = execute + name + string("\" &");
  cout << " Executing " << execute.c_str() << endl;
  system(execute.c_str());
}

///////////////////////////////////////////////////////////////////////
// pass a field with a distance field (inside/outside) overlay
// to a overlayFieldViewer3D
///////////////////////////////////////////////////////////////////////
void FIELD_3D::overlayFieldViewerYZ(const FIELD_3D& field, const FIELD_3D& distance, string name)
{
  field.write("./temp/temp3d.field");
  distance.write("./temp/temp3d.distance.field");
  string execute("./bin/overlayFieldViewer3D ./temp/temp3d.field ./temp/temp3d.distance.field \"");
  execute = execute + name + string("\" yz &");
  cout << " Executing " << execute.c_str() << endl;
  system(execute.c_str());
}

///////////////////////////////////////////////////////////////////////
// get the integer indices of a spatial position
///////////////////////////////////////////////////////////////////////
void FIELD_3D::indices(const VEC3F& position, int* x)
{
  VEC3F positionCopy = position;

  // get the lower corner position
  VEC3F corner = _center - (Real)0.5 * _lengths;
  VEC3F dxs(_dx, _dy, _dz);
  corner += (Real)0.5 * dxs;

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  x[0] = (int)positionCopy[0];
  x[1] = (int)positionCopy[1];
  x[2] = (int)positionCopy[2];
}

///////////////////////////////////////////////////////////////////////
// load a PhysBAM level set
///////////////////////////////////////////////////////////////////////
void FIELD_3D::readPhysBAM(const char* filename)
{
  if (_data) delete[] _data;

  // level set contains
  //
  // counts (TV_INT) vector<int, 3>
  // domain (RANGE) probably 6 floats
  // mac_offset (T) single floating point, 0 or 0.5
  //
  // then is reads in a scalar array, which contains
  //
  // length2 - an int, which seems to always equal 1
  // domain (RANGE<TV>) not TV_INT, float ranges, again?
  // the entries - scalars, 
  // everything appears to be single precision

  FILE* file = fopen(filename, "rb");

  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Could not open " << filename << "!" << endl;
    exit(0);
  }

  int xRes, yRes, zRes;

  // counts - grid resolutions without padding
  fread((void*)&xRes, sizeof(int), 1, file);
  fread((void*)&yRes, sizeof(int), 1, file);
  fread((void*)&zRes, sizeof(int), 1, file);

  // domain
  float xMin, yMin, zMin;
  float xMax, yMax, zMax;
  fread((void*)&xMin, sizeof(float), 1, file);
  fread((void*)&yMin, sizeof(float), 1, file);
  fread((void*)&zMin, sizeof(float), 1, file);
  fread((void*)&xMax, sizeof(float), 1, file);
  fread((void*)&yMax, sizeof(float), 1, file);
  fread((void*)&zMax, sizeof(float), 1, file);

  // MAC offset
  float macOffset;
  fread((void*)&macOffset, sizeof(float), 1, file);

  // length2
  int length2;
  fread((void*)&length2, sizeof(int), 1, file);

  // domain (grid resolutions with padding)
  int xPaddedMin, xPaddedMax;
  int yPaddedMin, yPaddedMax;
  int zPaddedMin, zPaddedMax;
  fread((void*)&xPaddedMin, sizeof(int), 1, file);
  fread((void*)&xPaddedMax, sizeof(int), 1, file);
  fread((void*)&yPaddedMin, sizeof(int), 1, file);
  fread((void*)&yPaddedMax, sizeof(int), 1, file);
  fread((void*)&zPaddedMin, sizeof(int), 1, file);
  fread((void*)&zPaddedMax, sizeof(int), 1, file);

  _xRes = (xPaddedMax - xPaddedMin) + 1;
  _yRes = (yPaddedMax - yPaddedMin) + 1;
  _zRes = (zPaddedMax - zPaddedMin) + 1;

  int temp = _xRes;
  _xRes = _zRes;
  _zRes = temp;

  // set longest dimension to 1
  int biggest = (_xRes > _yRes) ? _xRes : _yRes;
  biggest = (_zRes > biggest) ? _zRes : biggest;
  _lengths[0] = (Real)_xRes / biggest;
  _lengths[1] = (Real)_yRes / biggest;
  _lengths[2] = (Real)_zRes / biggest;

  cout << " lengths: " << _lengths << endl;

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
  _invDz = 1.0 / _dz;

  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  _data = new Real[_totalCells];

  if (sizeof(Real) == sizeof(float))
    fread((void*)_data, sizeof(float), _totalCells, file);
  else
    for (int x = 0; x < _totalCells; x++)
    {
      float single;
      fread((void*)&single, sizeof(float), 1, file);
      _data[x] = single;
    }
  fclose(file);

  // find the biggest length
  Real maxLength = (_lengths[1] > _lengths[0]) ? _lengths[1] : _lengths[0];
  maxLength = (_lengths[2] > maxLength) ? _lengths[2] : maxLength;

  (*this) *= 1.0 / maxLength;
}

///////////////////////////////////////////////////////////////////////
// load a PhysBAM level set
///////////////////////////////////////////////////////////////////////
void FIELD_3D::readPhysBAMGz(const char* filename)
{
  // level set contains
  //
  // counts (TV_INT) vector<int, 3>
  // domain (RANGE) probably 6 floats
  // mac_offset (T) single floating point, 0 or 0.5
  //
  // then is reads in a scalar array, which contains
  //
  // length2 - an int, which seems to always equal 1
  // domain (RANGE<TV>) not TV_INT, float ranges, again?
  // the entries - scalars, 
  // everything appears to be single precision

  gzFile file;
  file = gzopen(filename, "rb1");

  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Could not open " << filename << "!" << endl;
    exit(0);
  }
  cout << " Reading file " << filename << " ... "; flush(cout);

  // counts - grid resolutions without padding
  gzread(file, (void*)&_xRes, sizeof(int));
  gzread(file, (void*)&_yRes, sizeof(int));
  gzread(file, (void*)&_zRes, sizeof(int));

  cout << " Resolution: " << _xRes << " " << _yRes << " " << _zRes << " "; flush(cout);

  // domain
  float xMin, yMin, zMin;
  float xMax, yMax, zMax;
  gzread(file, (void*)&xMin, sizeof(float));
  gzread(file, (void*)&yMin, sizeof(float));
  gzread(file, (void*)&zMin, sizeof(float));
  gzread(file, (void*)&xMax, sizeof(float));
  gzread(file, (void*)&yMax, sizeof(float));
  gzread(file, (void*)&zMax, sizeof(float));

  // MAC offset
  float macOffset;
  gzread(file, (void*)&macOffset, sizeof(float));

  // length2
  int length2;
  gzread(file, (void*)&length2, sizeof(int));

  // domain (grid resolutions with padding)
  int xPaddedMin, xPaddedMax;
  int yPaddedMin, yPaddedMax;
  int zPaddedMin, zPaddedMax;
  gzread(file, (void*)&xPaddedMin, sizeof(int));
  gzread(file, (void*)&xPaddedMax, sizeof(int));
  gzread(file, (void*)&yPaddedMin, sizeof(int));
  gzread(file, (void*)&yPaddedMax, sizeof(int));
  gzread(file, (void*)&zPaddedMin, sizeof(int));
  gzread(file, (void*)&zPaddedMax, sizeof(int));

  _xRes = (xPaddedMax - xPaddedMin) + 1;
  _yRes = (yPaddedMax - yPaddedMin) + 1;
  _zRes = (zPaddedMax - zPaddedMin) + 1;

  int temp = _xRes;
  _xRes = _zRes;
  _zRes = temp;

  _lengths[0] = 1;
  _lengths[1] = 1;
  _lengths[2] = 1;

  // set longest dimension to 1
  int biggest = (_xRes > _yRes) ? _xRes : _yRes;
  biggest = (_zRes > biggest) ? _zRes : biggest;
  _lengths[0] = (Real)_xRes / biggest;
  _lengths[1] = (Real)_yRes / biggest;
  _lengths[2] = (Real)_zRes / biggest;

  cout << " lengths: " << _lengths << endl;

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
  _invDz = 1.0 / _dz;

  cout << " dx: " << _dx << " dy: " << _dy << " dz: " << _dz << endl;

  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;

  if (_data) delete[] _data;
  _data = new Real[_totalCells];
  
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        float single;
        gzread(file, (void*)&single, sizeof(float));
        (*this)(x,y,z) = single;
      }

  cout << " done." << endl;
  gzclose(file);
}

///////////////////////////////////////////////////////////////////////
// set to a checkerboard solid texture
///////////////////////////////////////////////////////////////////////
void FIELD_3D::setToSolidCheckboard(int xChecks, int yChecks, int zChecks)
{
  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, index++)
      {
        int xMod = (x / (_xRes / xChecks)) % 2;
        int yMod = (y / (_yRes / yChecks)) % 2;
        int zMod = (z / (_zRes / zChecks)) % 2;

        if (((xMod && yMod) || (!xMod && !yMod)) && zMod)
          _data[index] = 1;
        
        if (!((xMod && yMod) || (!xMod && !yMod)) && !zMod)
          _data[index] = 1;
      }
}

///////////////////////////////////////////////////////////////////////
// set to a checkerboard solid texture
///////////////////////////////////////////////////////////////////////
void FIELD_3D::setToGrayCheckerboard(int xChecks, int yChecks, int zChecks)
{
  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, index++)
      {
        int xMod = (x / (_xRes / xChecks)) % 2;
        int yMod = (y / (_yRes / yChecks)) % 2;
        int zMod = (z / (_zRes / zChecks)) % 2;

        if (((xMod && yMod) || (!xMod && !yMod)) && zMod)
          //_data[index] = 0.75;
          _data[index] = 0.25;
        else if (!((xMod && yMod) || (!xMod && !yMod)) && !zMod)
          //_data[index] = 0.75;
          _data[index] = 0.25;
        else
          _data[index] = -0.25;
      }
}

///////////////////////////////////////////////////////////////////////
// triqunitic interpolation lookup
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::quinticLookup(const VEC3F& position) const
{
  VEC3F positionCopy = position;

  // get the lower corner position
  VEC3F corner = _center - (Real)0.5 * _lengths;
  VEC3F dxs(_dx, _dy, _dz);
  corner += (Real)0.5 * dxs;

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  const int x1 = (int)positionCopy[0];
  const int x2    = x1 + 1;
  const int x3    = x1 + 2;
  const int x0    = x1 - 1;

  const int y1 = (int)positionCopy[1];
  const int y2    = y1 + 1;
  const int y3    = y1 + 2;
  const int y0    = y1 - 1;
  
  const int z1 = (int)positionCopy[2];
  const int z2    = z1 + 1;
  const int z3    = z1 + 2;
  const int z0    = z1 - 1;

  if (x0 < 0 || y0 < 0 || z0 < 0 ||
      x3 >= _xRes || y3 >= _yRes || z3 >= _zRes)
    return (*this)(position);

  const float xInterp = positionCopy[0]- x1;
  const float yInterp = positionCopy[1]- y1;
  const float zInterp = positionCopy[2]- z1;

  const int z0Slab = z0 * _slabSize;
  const int z1Slab = z1 * _slabSize;
  const int z2Slab = z2 * _slabSize;
  const int z3Slab = z3 * _slabSize;
  
  const int y0x = y0 * _xRes;
  const int y1x = y1 * _xRes;
  const int y2x = y2 * _xRes;
  const int y3x = y3 * _xRes;

  const int y0z0 = y0x + z0Slab;
  const int y1z0 = y1x + z0Slab;
  const int y2z0 = y2x + z0Slab;
  const int y3z0 = y3x + z0Slab;

  const int y0z1 = y0x + z1Slab;
  const int y1z1 = y1x + z1Slab;
  const int y2z1 = y2x + z1Slab;
  const int y3z1 = y3x + z1Slab;

  const int y0z2 = y0x + z2Slab;
  const int y1z2 = y1x + z2Slab;
  const int y2z2 = y2x + z2Slab;
  const int y3z2 = y3x + z2Slab;

  const int y0z3 = y0x + z3Slab;
  const int y1z3 = y1x + z3Slab;
  const int y2z3 = y2x + z3Slab;
  const int y3z3 = y3x + z3Slab;

  // do the z0 slice
  const Real p0[] = {_data[x0 + y0z0], _data[x1 + y0z0], _data[x2 + y0z0], _data[x3 + y0z0]};
  const Real p1[] = {_data[x0 + y1z0], _data[x1 + y1z0], _data[x2 + y1z0], _data[x3 + y1z0]};
  const Real p2[] = {_data[x0 + y2z0], _data[x1 + y2z0], _data[x2 + y2z0], _data[x3 + y2z0]};
  const Real p3[] = {_data[x0 + y3z0], _data[x1 + y3z0], _data[x2 + y3z0], _data[x3 + y3z0]};

  // do the z1 slice
  const Real p4[] = {_data[x0 + y0z1], _data[x1 + y0z1], _data[x2 + y0z1], _data[x3 + y0z1]};
  const Real p5[]= {_data[x0 + y1z1], _data[x1 + y1z1], _data[x2 + y1z1], _data[x3 + y1z1]};
  const Real p6[] = {_data[x0 + y2z1], _data[x1 + y2z1], _data[x2 + y2z1], _data[x3 + y2z1]};
  const Real p7[] = {_data[x0 + y3z1], _data[x1 + y3z1], _data[x2 + y3z1], _data[x3 + y3z1]};

  // do the z2 slice
  const Real p8[] = {_data[x0 + y0z2], _data[x1 + y0z2], _data[x2 + y0z2], _data[x3 + y0z2]};
  const Real p9[] = {_data[x0 + y1z2], _data[x1 + y1z2], _data[x2 + y1z2], _data[x3 + y1z2]};
  const Real p10[] = {_data[x0 + y2z2], _data[x1 + y2z2], _data[x2 + y2z2], _data[x3 + y2z2]};
  const Real p11[] = {_data[x0 + y3z2], _data[x1 + y3z2], _data[x2 + y3z2], _data[x3 + y3z2]};

  // do the z3 slice
  const Real p12[] = {_data[x0 + y0z3], _data[x1 + y0z3], _data[x2 + y0z3], _data[x3 + y0z3]};
  const Real p13[] = {_data[x0 + y1z3], _data[x1 + y1z3], _data[x2 + y1z3], _data[x3 + y1z3]};
  const Real p14[] = {_data[x0 + y2z3], _data[x1 + y2z3], _data[x2 + y2z3], _data[x3 + y2z3]};
  const Real p15[] = {_data[x0 + y3z3], _data[x1 + y3z3], _data[x2 + y3z3], _data[x3 + y3z3]};
  
  const Real z0Points[] = {quinticInterp(xInterp, p0), quinticInterp(xInterp, p1), quinticInterp(xInterp, p2), quinticInterp(xInterp, p3)};
  const Real z1Points[] = {quinticInterp(xInterp, p4), quinticInterp(xInterp, p5), quinticInterp(xInterp, p6), quinticInterp(xInterp, p7)};
  const Real z2Points[] = {quinticInterp(xInterp, p8), quinticInterp(xInterp, p9), quinticInterp(xInterp, p10), quinticInterp(xInterp, p11)};
  const Real z3Points[] = {quinticInterp(xInterp, p12), quinticInterp(xInterp, p13), quinticInterp(xInterp, p14), quinticInterp(xInterp, p15)};

  const Real finalPoints[] = {quinticInterp(yInterp, z0Points), quinticInterp(yInterp, z1Points), quinticInterp(yInterp, z2Points), quinticInterp(yInterp, z3Points)};

  return quinticInterp(zInterp, finalPoints);
}

///////////////////////////////////////////////////////////////////////
// tricubic interpolation lookup
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::cubicLookup(const VEC3F& position) const
{
  VEC3F positionCopy = position;

  // get the lower corner position
  const VEC3F corner = _center - (Real)0.5 * _lengths + (Real)0.5 * dxs();

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  const int x1 = (int)positionCopy[0];
  const int x2    = x1 + 1;
  const int x3    = x1 + 2;
  const int x0    = x1 - 1;

  const int y1 = (int)positionCopy[1];
  const int y2    = y1 + 1;
  const int y3    = y1 + 2;
  const int y0    = y1 - 1;
  
  const int z1 = (int)positionCopy[2];
  const int z2    = z1 + 1;
  const int z3    = z1 + 2;
  const int z0    = z1 - 1;

  if (x0 < 0 || y0 < 0 || z0 < 0 ||
      x3 >= _xRes || y3 >= _yRes || z3 >= _zRes)
    return (*this)(position);

  const float xInterp = positionCopy[0] - x1;
  const float yInterp = positionCopy[1] - y1;
  const float zInterp = positionCopy[2] - z1;

  const int z0Slab = z0 * _slabSize;
  const int z1Slab = z1 * _slabSize;
  const int z2Slab = z2 * _slabSize;
  const int z3Slab = z3 * _slabSize;
  
  const int y0x = y0 * _xRes;
  const int y1x = y1 * _xRes;
  const int y2x = y2 * _xRes;
  const int y3x = y3 * _xRes;

  const int y0z0 = y0x + z0Slab;
  const int y1z0 = y1x + z0Slab;
  const int y2z0 = y2x + z0Slab;
  const int y3z0 = y3x + z0Slab;

  const int y0z1 = y0x + z1Slab;
  const int y1z1 = y1x + z1Slab;
  const int y2z1 = y2x + z1Slab;
  const int y3z1 = y3x + z1Slab;

  const int y0z2 = y0x + z2Slab;
  const int y1z2 = y1x + z2Slab;
  const int y2z2 = y2x + z2Slab;
  const int y3z2 = y3x + z2Slab;

  const int y0z3 = y0x + z3Slab;
  const int y1z3 = y1x + z3Slab;
  const int y2z3 = y2x + z3Slab;
  const int y3z3 = y3x + z3Slab;

  // do the z0 slice
  const Real p0[] = {_data[x0 + y0z0], _data[x1 + y0z0], _data[x2 + y0z0], _data[x3 + y0z0]};
  const Real p1[] = {_data[x0 + y1z0], _data[x1 + y1z0], _data[x2 + y1z0], _data[x3 + y1z0]};
  const Real p2[] = {_data[x0 + y2z0], _data[x1 + y2z0], _data[x2 + y2z0], _data[x3 + y2z0]};
  const Real p3[] = {_data[x0 + y3z0], _data[x1 + y3z0], _data[x2 + y3z0], _data[x3 + y3z0]};

  // do the z1 slice
  const Real p4[] = {_data[x0 + y0z1], _data[x1 + y0z1], _data[x2 + y0z1], _data[x3 + y0z1]};
  const Real p5[]= {_data[x0 + y1z1], _data[x1 + y1z1], _data[x2 + y1z1], _data[x3 + y1z1]};
  const Real p6[] = {_data[x0 + y2z1], _data[x1 + y2z1], _data[x2 + y2z1], _data[x3 + y2z1]};
  const Real p7[] = {_data[x0 + y3z1], _data[x1 + y3z1], _data[x2 + y3z1], _data[x3 + y3z1]};

  // do the z2 slice
  const Real p8[] = {_data[x0 + y0z2], _data[x1 + y0z2], _data[x2 + y0z2], _data[x3 + y0z2]};
  const Real p9[] = {_data[x0 + y1z2], _data[x1 + y1z2], _data[x2 + y1z2], _data[x3 + y1z2]};
  const Real p10[] = {_data[x0 + y2z2], _data[x1 + y2z2], _data[x2 + y2z2], _data[x3 + y2z2]};
  const Real p11[] = {_data[x0 + y3z2], _data[x1 + y3z2], _data[x2 + y3z2], _data[x3 + y3z2]};

  // do the z3 slice
  const Real p12[] = {_data[x0 + y0z3], _data[x1 + y0z3], _data[x2 + y0z3], _data[x3 + y0z3]};
  const Real p13[] = {_data[x0 + y1z3], _data[x1 + y1z3], _data[x2 + y1z3], _data[x3 + y1z3]};
  const Real p14[] = {_data[x0 + y2z3], _data[x1 + y2z3], _data[x2 + y2z3], _data[x3 + y2z3]};
  const Real p15[] = {_data[x0 + y3z3], _data[x1 + y3z3], _data[x2 + y3z3], _data[x3 + y3z3]};
  
  const Real z0Points[] = {cubicInterp(xInterp, p0), cubicInterp(xInterp, p1), cubicInterp(xInterp, p2), cubicInterp(xInterp, p3)};
  const Real z1Points[] = {cubicInterp(xInterp, p4), cubicInterp(xInterp, p5), cubicInterp(xInterp, p6), cubicInterp(xInterp, p7)};
  const Real z2Points[] = {cubicInterp(xInterp, p8), cubicInterp(xInterp, p9), cubicInterp(xInterp, p10), cubicInterp(xInterp, p11)};
  const Real z3Points[] = {cubicInterp(xInterp, p12), cubicInterp(xInterp, p13), cubicInterp(xInterp, p14), cubicInterp(xInterp, p15)};

  const Real finalPoints[] = {cubicInterp(yInterp, z0Points), cubicInterp(yInterp, z1Points), cubicInterp(yInterp, z2Points), cubicInterp(yInterp, z3Points)};

  return cubicInterp(zInterp, finalPoints);
}

///////////////////////////////////////////////////////////////////////
// triquartic interpolation lookup
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::quarticLookup(const VEC3F& position) const
{
  VEC3F positionCopy = position;

  // get the lower corner position
  const VEC3F corner = _center - (Real)0.5 * _lengths + (Real)0.5 * dxs();

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  const int x1 = (int)positionCopy[0];
  const int x2    = x1 + 1;
  const int x3    = x1 + 2;
  const int x0    = x1 - 1;

  const int y1 = (int)positionCopy[1];
  const int y2    = y1 + 1;
  const int y3    = y1 + 2;
  const int y0    = y1 - 1;
  
  const int z1 = (int)positionCopy[2];
  const int z2    = z1 + 1;
  const int z3    = z1 + 2;
  const int z0    = z1 - 1;

  if (x0 < 0 || y0 < 0 || z0 < 0 ||
      x3 >= _xRes || y3 >= _yRes || z3 >= _zRes)
    return (*this)(position);

  const float xInterp = positionCopy[0] - x1;
  const float yInterp = positionCopy[1] - y1;
  const float zInterp = positionCopy[2] - z1;

  const int z0Slab = z0 * _slabSize;
  const int z1Slab = z1 * _slabSize;
  const int z2Slab = z2 * _slabSize;
  const int z3Slab = z3 * _slabSize;
  
  const int y0x = y0 * _xRes;
  const int y1x = y1 * _xRes;
  const int y2x = y2 * _xRes;
  const int y3x = y3 * _xRes;

  const int y0z0 = y0x + z0Slab;
  const int y1z0 = y1x + z0Slab;
  const int y2z0 = y2x + z0Slab;
  const int y3z0 = y3x + z0Slab;

  const int y0z1 = y0x + z1Slab;
  const int y1z1 = y1x + z1Slab;
  const int y2z1 = y2x + z1Slab;
  const int y3z1 = y3x + z1Slab;

  const int y0z2 = y0x + z2Slab;
  const int y1z2 = y1x + z2Slab;
  const int y2z2 = y2x + z2Slab;
  const int y3z2 = y3x + z2Slab;

  const int y0z3 = y0x + z3Slab;
  const int y1z3 = y1x + z3Slab;
  const int y2z3 = y2x + z3Slab;
  const int y3z3 = y3x + z3Slab;

  // do the z0 slice
  const Real p0[] = {_data[x0 + y0z0], _data[x1 + y0z0], _data[x2 + y0z0], _data[x3 + y0z0]};
  const Real p1[] = {_data[x0 + y1z0], _data[x1 + y1z0], _data[x2 + y1z0], _data[x3 + y1z0]};
  const Real p2[] = {_data[x0 + y2z0], _data[x1 + y2z0], _data[x2 + y2z0], _data[x3 + y2z0]};
  const Real p3[] = {_data[x0 + y3z0], _data[x1 + y3z0], _data[x2 + y3z0], _data[x3 + y3z0]};

  // do the z1 slice
  const Real p4[] = {_data[x0 + y0z1], _data[x1 + y0z1], _data[x2 + y0z1], _data[x3 + y0z1]};
  const Real p5[]= {_data[x0 + y1z1], _data[x1 + y1z1], _data[x2 + y1z1], _data[x3 + y1z1]};
  const Real p6[] = {_data[x0 + y2z1], _data[x1 + y2z1], _data[x2 + y2z1], _data[x3 + y2z1]};
  const Real p7[] = {_data[x0 + y3z1], _data[x1 + y3z1], _data[x2 + y3z1], _data[x3 + y3z1]};

  // do the z2 slice
  const Real p8[] = {_data[x0 + y0z2], _data[x1 + y0z2], _data[x2 + y0z2], _data[x3 + y0z2]};
  const Real p9[] = {_data[x0 + y1z2], _data[x1 + y1z2], _data[x2 + y1z2], _data[x3 + y1z2]};
  const Real p10[] = {_data[x0 + y2z2], _data[x1 + y2z2], _data[x2 + y2z2], _data[x3 + y2z2]};
  const Real p11[] = {_data[x0 + y3z2], _data[x1 + y3z2], _data[x2 + y3z2], _data[x3 + y3z2]};

  // do the z3 slice
  const Real p12[] = {_data[x0 + y0z3], _data[x1 + y0z3], _data[x2 + y0z3], _data[x3 + y0z3]};
  const Real p13[] = {_data[x0 + y1z3], _data[x1 + y1z3], _data[x2 + y1z3], _data[x3 + y1z3]};
  const Real p14[] = {_data[x0 + y2z3], _data[x1 + y2z3], _data[x2 + y2z3], _data[x3 + y2z3]};
  const Real p15[] = {_data[x0 + y3z3], _data[x1 + y3z3], _data[x2 + y3z3], _data[x3 + y3z3]};
  
  const Real z0Points[] = {quarticInterp(xInterp, p0), quarticInterp(xInterp, p1), quarticInterp(xInterp, p2), quarticInterp(xInterp, p3)};
  const Real z1Points[] = {quarticInterp(xInterp, p4), quarticInterp(xInterp, p5), quarticInterp(xInterp, p6), quarticInterp(xInterp, p7)};
  const Real z2Points[] = {quarticInterp(xInterp, p8), quarticInterp(xInterp, p9), quarticInterp(xInterp, p10), quarticInterp(xInterp, p11)};
  const Real z3Points[] = {quarticInterp(xInterp, p12), quarticInterp(xInterp, p13), quarticInterp(xInterp, p14), quarticInterp(xInterp, p15)};

  const Real finalPoints[] = {quarticInterp(yInterp, z0Points), quarticInterp(yInterp, z1Points), quarticInterp(yInterp, z2Points), quarticInterp(yInterp, z3Points)};

  return quarticInterp(zInterp, finalPoints);
}

///////////////////////////////////////////////////////////////////////
// tricubic interpolation lookup
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::cubicLookupUnclamped(const VEC3F& position) const
{
  VEC3F positionCopy = position;

  // get the lower corner position
  const VEC3F corner = _center - (Real)0.5 * _lengths + (Real)0.5 * dxs();

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  const int x1 = (int)positionCopy[0];
  const int x2    = x1 + 1;
  const int x3    = x1 + 2;
  const int x0    = x1 - 1;

  const int y1 = (int)positionCopy[1];
  const int y2    = y1 + 1;
  const int y3    = y1 + 2;
  const int y0    = y1 - 1;
  
  const int z1 = (int)positionCopy[2];
  const int z2    = z1 + 1;
  const int z3    = z1 + 2;
  const int z0    = z1 - 1;

  if (x0 < 0 || y0 < 0 || z0 < 0 ||
      x3 >= _xRes || y3 >= _yRes || z3 >= _zRes)
    return (*this)(position);

  const float xInterp = positionCopy[0] - x1;
  const float yInterp = positionCopy[1] - y1;
  const float zInterp = positionCopy[2] - z1;

  const int z0Slab = z0 * _slabSize;
  const int z1Slab = z1 * _slabSize;
  const int z2Slab = z2 * _slabSize;
  const int z3Slab = z3 * _slabSize;
  
  const int y0x = y0 * _xRes;
  const int y1x = y1 * _xRes;
  const int y2x = y2 * _xRes;
  const int y3x = y3 * _xRes;

  const int y0z0 = y0x + z0Slab;
  const int y1z0 = y1x + z0Slab;
  const int y2z0 = y2x + z0Slab;
  const int y3z0 = y3x + z0Slab;

  const int y0z1 = y0x + z1Slab;
  const int y1z1 = y1x + z1Slab;
  const int y2z1 = y2x + z1Slab;
  const int y3z1 = y3x + z1Slab;

  const int y0z2 = y0x + z2Slab;
  const int y1z2 = y1x + z2Slab;
  const int y2z2 = y2x + z2Slab;
  const int y3z2 = y3x + z2Slab;

  const int y0z3 = y0x + z3Slab;
  const int y1z3 = y1x + z3Slab;
  const int y2z3 = y2x + z3Slab;
  const int y3z3 = y3x + z3Slab;

  // do the z0 slice
  const Real p0[] = {_data[x0 + y0z0], _data[x1 + y0z0], _data[x2 + y0z0], _data[x3 + y0z0]};
  const Real p1[] = {_data[x0 + y1z0], _data[x1 + y1z0], _data[x2 + y1z0], _data[x3 + y1z0]};
  const Real p2[] = {_data[x0 + y2z0], _data[x1 + y2z0], _data[x2 + y2z0], _data[x3 + y2z0]};
  const Real p3[] = {_data[x0 + y3z0], _data[x1 + y3z0], _data[x2 + y3z0], _data[x3 + y3z0]};

  // do the z1 slice
  const Real p4[] = {_data[x0 + y0z1], _data[x1 + y0z1], _data[x2 + y0z1], _data[x3 + y0z1]};
  const Real p5[]= {_data[x0 + y1z1], _data[x1 + y1z1], _data[x2 + y1z1], _data[x3 + y1z1]};
  const Real p6[] = {_data[x0 + y2z1], _data[x1 + y2z1], _data[x2 + y2z1], _data[x3 + y2z1]};
  const Real p7[] = {_data[x0 + y3z1], _data[x1 + y3z1], _data[x2 + y3z1], _data[x3 + y3z1]};

  // do the z2 slice
  const Real p8[] = {_data[x0 + y0z2], _data[x1 + y0z2], _data[x2 + y0z2], _data[x3 + y0z2]};
  const Real p9[] = {_data[x0 + y1z2], _data[x1 + y1z2], _data[x2 + y1z2], _data[x3 + y1z2]};
  const Real p10[] = {_data[x0 + y2z2], _data[x1 + y2z2], _data[x2 + y2z2], _data[x3 + y2z2]};
  const Real p11[] = {_data[x0 + y3z2], _data[x1 + y3z2], _data[x2 + y3z2], _data[x3 + y3z2]};

  // do the z3 slice
  const Real p12[] = {_data[x0 + y0z3], _data[x1 + y0z3], _data[x2 + y0z3], _data[x3 + y0z3]};
  const Real p13[] = {_data[x0 + y1z3], _data[x1 + y1z3], _data[x2 + y1z3], _data[x3 + y1z3]};
  const Real p14[] = {_data[x0 + y2z3], _data[x1 + y2z3], _data[x2 + y2z3], _data[x3 + y2z3]};
  const Real p15[] = {_data[x0 + y3z3], _data[x1 + y3z3], _data[x2 + y3z3], _data[x3 + y3z3]};
  
  const Real z0Points[] = {cubicInterpUnclamped(xInterp, p0), cubicInterp(xInterp, p1), cubicInterp(xInterp, p2), cubicInterp(xInterp, p3)};
  const Real z1Points[] = {cubicInterpUnclamped(xInterp, p4), cubicInterp(xInterp, p5), cubicInterp(xInterp, p6), cubicInterp(xInterp, p7)};
  const Real z2Points[] = {cubicInterpUnclamped(xInterp, p8), cubicInterp(xInterp, p9), cubicInterp(xInterp, p10), cubicInterp(xInterp, p11)};
  const Real z3Points[] = {cubicInterpUnclamped(xInterp, p12), cubicInterp(xInterp, p13), cubicInterp(xInterp, p14), cubicInterp(xInterp, p15)};

  const Real finalPoints[] = {cubicInterpUnclamped(yInterp, z0Points), cubicInterp(yInterp, z1Points), cubicInterp(yInterp, z2Points), cubicInterp(yInterp, z3Points)};

  return cubicInterpUnclamped(zInterp, finalPoints);
}

///////////////////////////////////////////////////////////////////////
// tricubic interpolation lookup
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::cubicNewtonLookup(const VEC3F& position) const
{
  VEC3F positionCopy = position;

  // get the lower corner position
  VEC3F corner = _center - (Real)0.5 * _lengths;
  VEC3F dxs(_dx, _dy, _dz);
  corner += (Real)0.5 * dxs;

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  int x1 = (int)positionCopy[0];
  int x2    = x1 + 1;
  int x3    = x1 + 2;
  int x0    = x1 - 1;

  int y1 = (int)positionCopy[1];
  int y2    = y1 + 1;
  int y3    = y1 + 2;
  int y0    = y1 - 1;
  
  int z1 = (int)positionCopy[2];
  int z2    = z1 + 1;
  int z3    = z1 + 2;
  int z0    = z1 - 1;

  if (x0 < 0 || y0 < 0 || z0 < 0 ||
      x3 >= _xRes || y3 >= _yRes || z3 >= _zRes)
    return (*this)(position);

  Real points[4];
  const float xInterp = positionCopy[0]- x1;
  const float yInterp = positionCopy[1]- y1;
  const float zInterp = positionCopy[2]- z1;

  Real finalPoints[4];

  // do the z0 slice
  Real z0Points[4];
  points[0] = _data[x0 + y0 * _xRes + z0 * _slabSize];
  points[1] = _data[x1 + y0 * _xRes + z0 * _slabSize];
  points[2] = _data[x2 + y0 * _xRes + z0 * _slabSize];
  points[3] = _data[x3 + y0 * _xRes + z0 * _slabSize];
  z0Points[0] = cubicNewtonInterp(xInterp, points);

  points[0] = _data[x0 + y1 * _xRes + z0 * _slabSize];
  points[1] = _data[x1 + y1 * _xRes + z0 * _slabSize];
  points[2] = _data[x2 + y1 * _xRes + z0 * _slabSize];
  points[3] = _data[x3 + y1 * _xRes + z0 * _slabSize];
  z0Points[1] = cubicNewtonInterp(xInterp, points);
  
  points[0] = _data[x0 + y2 * _xRes + z0 * _slabSize];
  points[1] = _data[x1 + y2 * _xRes + z0 * _slabSize];
  points[2] = _data[x2 + y2 * _xRes + z0 * _slabSize];
  points[3] = _data[x3 + y2 * _xRes + z0 * _slabSize];
  z0Points[2] = cubicNewtonInterp(xInterp, points);

  points[0] = _data[x0 + y3 * _xRes + z0 * _slabSize];
  points[1] = _data[x1 + y3 * _xRes + z0 * _slabSize];
  points[2] = _data[x2 + y3 * _xRes + z0 * _slabSize];
  points[3] = _data[x3 + y3 * _xRes + z0 * _slabSize];
  z0Points[3] = cubicNewtonInterp(xInterp, points);
  finalPoints[0] = cubicNewtonInterp(yInterp, z0Points);

  // do the z1 slice
  Real z1Points[4];
  points[0] = _data[x0 + y0 * _xRes + z1 * _slabSize];
  points[1] = _data[x1 + y0 * _xRes + z1 * _slabSize];
  points[2] = _data[x2 + y0 * _xRes + z1 * _slabSize];
  points[3] = _data[x3 + y0 * _xRes + z1 * _slabSize];
  z1Points[0] = cubicNewtonInterp(xInterp, points);

  points[0] = _data[x0 + y1 * _xRes + z1 * _slabSize];
  points[1] = _data[x1 + y1 * _xRes + z1 * _slabSize];
  points[2] = _data[x2 + y1 * _xRes + z1 * _slabSize];
  points[3] = _data[x3 + y1 * _xRes + z1 * _slabSize];
  z1Points[1] = cubicNewtonInterp(xInterp, points);
  
  points[0] = _data[x0 + y2 * _xRes + z1 * _slabSize];
  points[1] = _data[x1 + y2 * _xRes + z1 * _slabSize];
  points[2] = _data[x2 + y2 * _xRes + z1 * _slabSize];
  points[3] = _data[x3 + y2 * _xRes + z1 * _slabSize];
  z1Points[2] = cubicNewtonInterp(xInterp, points);

  points[0] = _data[x0 + y3 * _xRes + z1 * _slabSize];
  points[1] = _data[x1 + y3 * _xRes + z1 * _slabSize];
  points[2] = _data[x2 + y3 * _xRes + z1 * _slabSize];
  points[3] = _data[x3 + y3 * _xRes + z1 * _slabSize];
  z1Points[3] = cubicNewtonInterp(xInterp, points);
  finalPoints[1] = cubicNewtonInterp(yInterp, z1Points);

  // do the z2 slice
  Real z2Points[4];
  points[0] = _data[x0 + y0 * _xRes + z2 * _slabSize];
  points[1] = _data[x1 + y0 * _xRes + z2 * _slabSize];
  points[2] = _data[x2 + y0 * _xRes + z2 * _slabSize];
  points[3] = _data[x3 + y0 * _xRes + z2 * _slabSize];
  z2Points[0] = cubicNewtonInterp(xInterp, points);

  points[0] = _data[x0 + y1 * _xRes + z2 * _slabSize];
  points[1] = _data[x1 + y1 * _xRes + z2 * _slabSize];
  points[2] = _data[x2 + y1 * _xRes + z2 * _slabSize];
  points[3] = _data[x3 + y1 * _xRes + z2 * _slabSize];
  z2Points[1] = cubicNewtonInterp(xInterp, points);
  
  points[0] = _data[x0 + y2 * _xRes + z2 * _slabSize];
  points[1] = _data[x1 + y2 * _xRes + z2 * _slabSize];
  points[2] = _data[x2 + y2 * _xRes + z2 * _slabSize];
  points[3] = _data[x3 + y2 * _xRes + z2 * _slabSize];
  z2Points[2] = cubicNewtonInterp(xInterp, points);

  points[0] = _data[x0 + y3 * _xRes + z2 * _slabSize];
  points[1] = _data[x1 + y3 * _xRes + z2 * _slabSize];
  points[2] = _data[x2 + y3 * _xRes + z2 * _slabSize];
  points[3] = _data[x3 + y3 * _xRes + z2 * _slabSize];
  z2Points[3] = cubicNewtonInterp(xInterp, points);
  finalPoints[2] = cubicNewtonInterp(yInterp, z2Points);

  // do the z3 slice
  Real z3Points[4];
  points[0] = _data[x0 + y0 * _xRes + z3 * _slabSize];
  points[1] = _data[x1 + y0 * _xRes + z3 * _slabSize];
  points[2] = _data[x2 + y0 * _xRes + z3 * _slabSize];
  points[3] = _data[x3 + y0 * _xRes + z3 * _slabSize];
  z3Points[0] = cubicNewtonInterp(xInterp, points);

  points[0] = _data[x0 + y1 * _xRes + z3 * _slabSize];
  points[1] = _data[x1 + y1 * _xRes + z3 * _slabSize];
  points[2] = _data[x2 + y1 * _xRes + z3 * _slabSize];
  points[3] = _data[x3 + y1 * _xRes + z3 * _slabSize];
  z3Points[1] = cubicNewtonInterp(xInterp, points);
  
  points[0] = _data[x0 + y2 * _xRes + z3 * _slabSize];
  points[1] = _data[x1 + y2 * _xRes + z3 * _slabSize];
  points[2] = _data[x2 + y2 * _xRes + z3 * _slabSize];
  points[3] = _data[x3 + y2 * _xRes + z3 * _slabSize];
  z3Points[2] = cubicNewtonInterp(xInterp, points);

  points[0] = _data[x0 + y3 * _xRes + z3 * _slabSize];
  points[1] = _data[x1 + y3 * _xRes + z3 * _slabSize];
  points[2] = _data[x2 + y3 * _xRes + z3 * _slabSize];
  points[3] = _data[x3 + y3 * _xRes + z3 * _slabSize];
  z3Points[3] = cubicNewtonInterp(xInterp, points);
  finalPoints[3] = cubicNewtonInterp(yInterp, z3Points);

  return cubicNewtonInterp(zInterp, finalPoints);
}

///////////////////////////////////////////////////////////////////////
// do a cubic Hermite that clamps to the immediate neighborhood
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::cubicInterpClamped(const Real interp, const Real* points)
{
  Real d0 = (points[2] - points[0]) * 0.5;
  Real d1 = (points[3] - points[1]) * 0.5;

  const Real deltak = (points[2] - points[1]);

  // do monotonic interpolation
  if (deltak * d0 < 0.0)
    d0 = 0;
  if (deltak * d1 < 0.0)
    d1 = 0;

  const Real a0 = points[1];
  const Real a1 = d0;
  const Real a2 = 3.0 * deltak - 2.0 * d0 - d1;
  const Real a3 = -2.0 * deltak + d0 + d1;

  const Real squared = interp * interp;
  const Real cubed = squared * interp;
  Real trial = a3 * cubed + a2 * squared + a1 * interp + a0;

  const Real intervalMax = (points[1] > points[2]) ? points[1] : points[2];
  const Real intervalMin = (points[1] > points[2]) ? points[2] : points[1];

  trial = (trial > intervalMax) ? intervalMax : trial;
  trial = (trial < intervalMin) ? intervalMin : trial;

  return trial;
}

///////////////////////////////////////////////////////////////////////
// do a cubic Hermite interpolation
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::cubicInterp(const Real interp, const Real* points)
{
  Real d0 = (points[2] - points[0]) * 0.5;
  Real d1 = (points[3] - points[1]) * 0.5;

  Real deltak = (points[2] - points[1]);

  // do monotonic interpolation
  if (deltak * d0 < 0.0)
    d0 = 0;
  if (deltak * d1 < 0.0)
    d1 = 0;

  Real a0 = points[1];
  Real a1 = d0;
  Real a2 = 3.0 * deltak - 2.0 * d0 - d1;
  Real a3 = -2.0 * deltak + d0 + d1;

  Real squared = interp * interp;
  Real cubed = squared * interp;
  return a3 * cubed + a2 * squared + a1 * interp + a0;
}

///////////////////////////////////////////////////////////////////////
// do a cubic Hermite interpolation
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::cubicInterpUnclamped(const Real interp, const Real* points)
{
  Real d0 = (points[2] - points[0]) * 0.5;
  Real d1 = (points[3] - points[1]) * 0.5;

  Real deltak = (points[2] - points[1]);

  Real a0 = points[1];
  Real a1 = d0;
  Real a2 = 3.0 * deltak - 2.0 * d0 - d1;
  Real a3 = -2.0 * deltak + d0 + d1;

  Real squared = interp * interp;
  Real cubed = squared * interp;
  return a3 * cubed + a2 * squared + a1 * interp + a0;
}

///////////////////////////////////////////////////////////////////////
// do a quintic Hermite interpolation
//
// The constraint matrix is:
// A = [ 0 0 0 0 0 1;0 0 0 0 1 0;0 0 0 2 0 0 ; 1 1 1 1 1 1 ;5 4 3 2 1 0; 20 12 6 2 0 0 ]
//
// where b = [y0 dy0 ddy0 d1 dy1 ddy1]';
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::quinticInterp(Real interp, const Real* points)
{
  // do monotone increasing
  if (points[1] > points[2])
  {
    const Real newPts[] = {points[3], points[2], points[1], points[0]};
    return quinticInterp(1 - interp, newPts); 
  }

  const Real dy0 = (points[2] - points[0]) * 0.5;
  const Real dy1= (points[3] - points[1]) * 0.5;
  const Real ddy0 = (points[2] - 2.0 * points[1] + points[0]);
  const Real ddy1 = (points[3] - 2.0 * points[2] + points[1]);
  const Real y0 = points[1];
  const Real y1 = points[2];

  // make it monotonic -- the pow calls cause a 10x slowdown!
  const Real A = 2 * dy0;
  const Real B = 2 * dy1;
  const Real C = ddy0;
  const Real D = ddy1;

  const Real alpha = 4.f * (B - D) / (pow(A, (Real)0.25) * pow(B, (Real)0.75));
  const Real beta  = 6.f * (D - C - 2.f * B - 2.f * A + 5.f) / (sqrt(A) * sqrt(B));
  const Real gamma = 4.f * (A + C) / (pow(A, (Real)0.75) * pow(B, (Real)0.25));

  const Real test = (beta > 6.f) ? -2.f * sqrt(beta - 2.f) : -(beta + 2.f) * 0.5f;
  if (!(alpha > test && gamma > test))
  {
    _quinticClamps++;
    return cubicInterp(interp, points);
  }

  // the quintic interpolation matrix
  const Real a =  -6. * y0 - 3. * dy0 - 0.5 * ddy0 + 6.  * y1 - 3. * dy1 + 0.5 * ddy1;
  const Real b =  15. * y0 + 8. * dy0 + 1.5 * ddy0 - 15. * y1 + 7. * dy1 - ddy1;
  const Real c = -10. * y0 - 6. * dy0 - 1.5 * ddy0 + 10. * y1 - 4. * dy1 + 0.5 * ddy1;
  const Real d =  0.5 *ddy0;
  const Real e = dy0; 
  const Real f = y0;

  const Real x2 = interp * interp;
  const Real x3 = x2 * interp;
  const Real x4 = x2 * x2;
  const Real x5 = x3 * x2;

  return a * x5 + b * x4 + c * x3 + d * x2 + e * interp + f;
}

///////////////////////////////////////////////////////////////////////
// do a quartic WENO interpolation
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::quarticInterp(const Real interp, const Real* points)
{
  const Real fim1 = points[0];
  const Real fi   = points[1];
  const Real fip1 = points[2];
  const Real fip2 = points[3];
  const Real x = interp;

  const Real p1 = fi + ((fip1 - fim1) + (fip1 - 2.0 * fi + fim1) * x) * 0.5 * x;
  const Real p2 = fi + ((-fip2 + 4.0 * fip1 - 3 * fi) + (fip2 - 2.0 * fip1 + fi) * x) * 0.5 * x;

  const Real C1 = (2 - x) / 3.0;
  const Real C2 = (x + 1) / 3.0;

  const Real middle = -76 * fip1 * fi;
  const Real fip1Sq = fip1 * fip1;
  const Real fiSq = fi * fi;

  const Real eps = 1e-6;
  const Real IS1 = (26 * fip1 * fim1 - 52 * fi * fim1 + middle + 25 * fip1Sq + 64 * fiSq + 13 * fim1 * fim1) / 12.0 + eps;
  const Real IS2 = (26 * fip2 * fi - 52 * fip2 * fip1 + middle + 25 * fiSq + 64 * fip1Sq + 13 * fip2 * fip2) / 12.0 + eps;

  const Real alpha1 = C1 / (IS1 * IS1);
  const Real alpha2 = C2 / (IS2 * IS2);

  const Real sum = alpha1 + alpha2;
  const Real w1 = alpha1 / sum;
  const Real w2 = alpha2 / sum;

  const Real final = w1 * p1 + w2 * p2;

  const Real intervalMax = (points[1] > points[2]) ? points[1] : points[2];
  const Real intervalMin = (points[1] > points[2]) ? points[2] : points[1];

  const Real linear = x * points[2] + (1.0 - x) * points[1];
  return (final < intervalMax) ? ((final > intervalMin) ? final : linear) : linear;
}

///////////////////////////////////////////////////////////////////////
// do a cubic Newton interpolation
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::cubicNewtonInterp(Real interp, Real* points)
{
  Real x = interp;
  const Real x1 = -1;
  const Real x2 = 0;
  const Real x3 = 1;
  const Real x4 = 2;

  Real y1 = points[0];
  Real y2 = points[1];
  Real y3 = points[2];
  Real y4 = points[3];

  Real c1 = y1;
  Real c2 = (y2 - c1) / (x2 - x1);
  Real c3 = (y3 - (c1 + c2 * (x3 - x1))) / ((x3 - x1) * (x3 - x2));
  Real c4 = (y4 - ((c1 + c2 * (x4 - x1) + c3 * (x4 - x1) * (x4 - x2)))) / ((x4 - x1) * (x4 - x2) * (x4 - x3));

  return c1 + c2 * (x - x1) + c3 * (x - x1) * (x - x2) + c4 * (x - x1) * (x - x2) * (x - x3);
}

///////////////////////////////////////////////////////////////////////
// Compute the elements of the vertical derivative convolution kernel
///////////////////////////////////////////////////////////////////////
void FIELD_3D::setToVerticalDerivativeKernel(double kMax, double dk, double sigma, double L)
{
  assert(_xRes % 2);
  assert(_xRes == _yRes);
  assert(_xRes == _zRes);

  double norm = 0;

  for (double k = 0; k < kMax; k += dk)
    norm += k * k * exp(-sigma * k * k);

  int halfWidth = _xRes / 2;

  for (int h = -halfWidth; h <= halfWidth; h++)
    for (int i = -halfWidth; i <= halfWidth; i++)
      for (int j = -halfWidth; j <= halfWidth; j++)
      {
        double r = sqrt((float)(i * i + j * j + h * h));
        double kern = 0;
        for (double k = 0; k < kMax; k += dk)
          kern += k * k * (sqrt(1.0 + k * k * L * L)) * exp(-sigma * k * k) * j0(r * k);
          //kern += k * k * (sqrt(k * k * L * L)) * exp(-sigma * k * k) * j0(r * k);

        (*this)(i + halfWidth, j + halfWidth, h + halfWidth) = kern / norm;
      }
}

///////////////////////////////////////////////////////////////////////
// set whole field to a Gaussian
///////////////////////////////////////////////////////////////////////
void FIELD_3D::setToGaussian(Real amplitude, VEC3F sigmas)
{
  float dx = 1.0f / (_xRes - 1);
  float dy = 1.0f / (_yRes - 1);
  float dz = 1.0f / (_zRes - 1);

  float x0 = 0.5;
  float y0 = 0.5;
  float z0 = 0.5;
  float xReal = 0;
  float yReal = 0;
  float zReal = 0;

  for (int z = 0; z < _zRes; z++)
  {
    for (int y = 0; y < _yRes; y++)
    {
      for (int x = 0; x < _xRes; x++)
      {
        int index = x + y * _xRes + z * _slabSize;

        float xValue = xReal;
        float yValue = yReal;
        float zValue = zReal;

        xValue = xValue - x0;
        xValue *= xValue;
        xValue *= 1.0 / (2.0 * sigmas[0] * sigmas[0]);

        yValue = yValue - y0;
        yValue *= yValue;
        yValue *= 1.0 / (2.0 * sigmas[1] * sigmas[1]);
        
        zValue = zValue - z0;
        zValue *= zValue;
        zValue *= 1.0 / (2.0 * sigmas[2] * sigmas[2]);

        _data[index] = exp(-(xValue + yValue + zValue));

        xReal += dx;
      }
      xReal = 0;
      yReal += dy;
    }
    yReal = 0;
    zReal += dz;
  }

  // normalize
  (*this) *= amplitude / sum();
}

///////////////////////////////////////////////////////////////////////
// convolve this field with a smaller field
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::convolve(const FIELD_3D& filter)
{
  TIMER functionTimer(__FUNCTION__);
  FIELD_3D final(*this);

  assert(filter.xRes() < _xRes);
  assert(filter.yRes() < _yRes);
  assert(filter.zRes() < _zRes);

  assert(filter.xRes() % 2);
  assert(filter.yRes() % 2);
  assert(filter.zRes() % 2);

  int xHalf = filter.xRes() / 2;
  int yHalf = filter.yRes() / 2;
  int zHalf = filter.zRes() / 2;

#pragma omp parallel
#pragma omp for  schedule(static)
  for (int z = zHalf; z < _zRes - zHalf; z++)
    for (int y = yHalf; y < _yRes - yHalf; y++)
      for (int x = xHalf; x < _xRes - xHalf; x++)
      {
        final(x,y,z) = 0;
        for (int fz = 0; fz < filter.zRes(); fz++)
          for (int fy = 0; fy < filter.yRes(); fy++)
            for (int fx = 0; fx < filter.xRes(); fx++)
              final(x,y,z) += filter(fx,fy,fz) * (*this)(x + (fx - xHalf), 
                                                         y + (fy - yHalf), 
                                                         z + (fz - zHalf)); 
      }

  return final;
}

///////////////////////////////////////////////////////////////////////
// convolve this field with a smaller field
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::convolveToroidal(const FIELD_3D& filter)
{
  FIELD_3D final(*this);

  assert(filter.xRes() < _xRes);
  assert(filter.yRes() < _yRes);
  assert(filter.zRes() < _zRes);

  assert(filter.xRes() % 2);
  assert(filter.yRes() % 2);
  assert(filter.zRes() % 2);

  int xHalf = filter.xRes() / 2;
  int yHalf = filter.yRes() / 2;
  int zHalf = filter.zRes() / 2;

#pragma omp parallel
  {
#pragma omp for  schedule(static)
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        final(x,y,z) = 0;
        for (int fz = 0; fz < filter.zRes(); fz++)
          for (int fy = 0; fy < filter.yRes(); fy++)
            for (int fx = 0; fx < filter.xRes(); fx++)
            {
              int ffx = x + fx - xHalf;
              int ffy = y + fy - yHalf;
              int ffz = z + fz - zHalf;

              if (ffx < 0) ffx += _xRes;
              if (ffy < 0) ffy += _yRes;
              if (ffz < 0) ffz += _zRes;

              ffx = ffx % _xRes;
              ffy = ffy % _yRes;
              ffz = ffz % _zRes;

              final(x,y,z) += filter(fx,fy,fz) * (*this)(ffx, ffy, ffz);
            }
      }
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
// convolve this field with a smaller field
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::convolveNarrowBand(const FIELD_3D& filter, const FIELD_3D& distance, int maxCells)
{
  FIELD_3D final(*this);
  final = 0;

  assert(filter.xRes() < _xRes);
  assert(filter.yRes() < _yRes);
  assert(filter.zRes() < _zRes);

  assert(filter.xRes() % 2);
  assert(filter.yRes() % 2);
  assert(filter.zRes() % 2);

  // for the time being, assume the fields are the same red
  assert(_xRes == distance.xRes());
  assert(_yRes == distance.yRes());
  assert(_zRes == distance.zRes());

  int xHalf = filter.xRes() / 2;
  int yHalf = filter.yRes() / 2;
  int zHalf = filter.zRes() / 2;

  Real invDx = 1.0 / distance.dx();

  // find which cells are inside the band
  vector<int> xs;
  vector<int> ys;
  vector<int> zs;
  for (int z = zHalf; z < _zRes - zHalf; z++)
    for (int y = yHalf; y < _yRes - yHalf; y++)
      for (int x = xHalf; x < _xRes - xHalf; x++)
      {
        Real currentDistance = fabs(distance(x,y,z) * invDx);
        if (currentDistance < maxCells)
        {
          xs.push_back(x);
          ys.push_back(y);
          zs.push_back(z);
        }
      }

  int size = xs.size();
#pragma omp parallel
#pragma omp for  schedule(static)
  for (int i = 0; i < size; i++)
  {
    int x = xs[i];
    int y = ys[i];
    int z = zs[i];

    for (int fz = 0; fz < filter.zRes(); fz++)
      for (int fy = 0; fy < filter.yRes(); fy++)
        for (int fx = 0; fx < filter.xRes(); fx++)
          final(x,y,z) += filter(fx,fy,fz) * (*this)(x + (fx - xHalf), 
                                                     y + (fy - yHalf), 
                                                     z + (fz - zHalf)); 
  }
  cout << endl;

  return final;
}

///////////////////////////////////////////////////////////////////////
// stomp all the value outside a narrow band to zero in order to 
// boost compression
///////////////////////////////////////////////////////////////////////
void FIELD_3D::stompOutsideNarrowBand(const FIELD_3D& distance, int maxCells)
{
  Real invDx = 1.0 / distance.dx();

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        Real currentDistance = fabs(distance(x,y,z) * invDx);

        if (currentDistance > maxCells)
          (*this)(x,y,z) = 0;
      }
}

///////////////////////////////////////////////////////////////////////
// get the sum of the field
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::sum() const
{
  Real final = 0;
  for (int x = 0; x < _totalCells; x++)
    final += _data[x];

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the sum of the field
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::absSum() const
{
  Real final = 0;
  for (int x = 0; x < _totalCells; x++)
    final += fabs(_data[x]);

  return final;
}

///////////////////////////////////////////////////////////////////////
// insert a Gaussian at a specific point
///////////////////////////////////////////////////////////////////////
void FIELD_3D::insertGaussian(const VEC3F& center, const Real amplitude, const VEC3F sigmas)
{
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        VEC3F cell = cellCenter(x,y,z);

        float xValue = cell[0] - center[0];
        float yValue = cell[1] - center[1];
        float zValue = cell[2] - center[2];

        xValue *= xValue;
        xValue *= 1.0 / (2.0 * sigmas[0] * sigmas[0]);

        yValue *= yValue;
        yValue *= 1.0 / (2.0 * sigmas[1] * sigmas[1]);
        
        zValue *= zValue;
        zValue *= 1.0 / (2.0 * sigmas[2] * sigmas[2]);

        int index = x + y * _xRes + z * _slabSize;
        _data[index] += exp(-(xValue + yValue + zValue));
      }

  (*this) *= amplitude;
}

///////////////////////////////////////////////////////////////////////
// reset dimensions
///////////////////////////////////////////////////////////////////////
void FIELD_3D::setLengths(const VEC3F& lengths)
{
  _lengths = lengths;

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
  _invDz = 1.0 / _dz;
}

///////////////////////////////////////////////////////////////////////
// determine how many non-zero entries are in the filter
///////////////////////////////////////////////////////////////////////
int FIELD_3D::nonZeroEntries()
{
  int nonZero = 0;

  for (int x = 0; x < _totalCells; x++)
    if (fabs(_data[x]) > 0)
      nonZero++;

  return nonZero;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::fieldMax() const
{
  assert(_totalCells > 0);

  Real final = _data[0];

  for (int x = 0; x < _totalCells; x++)
    if (_data[x] > final)
      final = _data[x];

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VEC3I FIELD_3D::maxIndex() const
{
  Real maxFound = _data[0];

  VEC3I maxFoundIndex;
  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, index++)
        if (_data[index] > maxFound)
        {
          maxFound = _data[index];

          maxFoundIndex[0] = x;
          maxFoundIndex[1] = y;
          maxFoundIndex[2] = z;
        }

  return maxFoundIndex;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VEC3F FIELD_3D::maxCellCenter() const
{
  Real maxFound = _data[0];

  VEC3I maxFoundIndex;
  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, index++)
        if (_data[index] > maxFound)
        {
          maxFound = _data[index];

          maxFoundIndex[0] = x;
          maxFoundIndex[1] = y;
          maxFoundIndex[2] = z;
        }

  return cellCenter(maxFoundIndex[0], maxFoundIndex[1], maxFoundIndex[2]);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::fieldMin() const
{
  assert(_totalCells > 0);

  Real final = _data[0];

  for (int x = 0; x < _totalCells; x++)
    if (_data[x] < final)
      final = _data[x];

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VEC3I FIELD_3D::minIndex() const
{
  Real minFound = _data[0];

  VEC3I minFoundIndex;
  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, index++)
        if (_data[index] < minFound)
        {
          minFound = _data[index];

          minFoundIndex[0] = x;
          minFoundIndex[1] = y;
          minFoundIndex[2] = z;
        }

  return minFoundIndex;
}

///////////////////////////////////////////////////////////////////////
// normalize the data
///////////////////////////////////////////////////////////////////////
void FIELD_3D::normalize()
{
  Real totalMin = fieldMin();
  Real totalMax = fieldMax();

  for (int x = 0; x < _totalCells; x++)
    _data[x] = (_data[x] - totalMin) / (totalMax - totalMin);
}

///////////////////////////////////////////////////////////////////////
// flip the z and y coordinates
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::flipZY() const
{
  FIELD_3D final(_xRes, _zRes, _yRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        final(x,z,y) = (*this)(x,y,z);

  return final;
}

///////////////////////////////////////////////////////////////////////
// flip the x and y coordinates
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::flipXY() const
{
  FIELD_3D final(_yRes, _xRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        final(y,x,z) = (*this)(x,y,z);

  return final;
}

///////////////////////////////////////////////////////////////////////
// flip the x and z coordinates
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::flipXZ() const
{
  FIELD_3D final(_zRes, _yRes, _xRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        final(z,y,x) = (*this)(x,y,z);

  return final;
}

///////////////////////////////////////////////////////////////////////
// create a mirror image along Z
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::mirrorZ() const
{
  FIELD_3D final(_xRes, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        final(x,y,z) = (*this)(x,y,(_zRes - 1) - z);

  return final;
}

///////////////////////////////////////////////////////////////////////
// create a mirror image along Y
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::mirrorY() const
{
  FIELD_3D final(_xRes, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        final(x,y,z) = (*this)(x,(_yRes - 1) - y,z);

  return final;
}

///////////////////////////////////////////////////////////////////////
// create a mirror image along Z
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::mirrorX() const
{
  FIELD_3D final(_xRes, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        final(x,y,z) = (*this)((_xRes - 1) - x,y,z);

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_3D::writeMatlab(string filename, string variableName) const
{
  string dummy = filename;
  dummy = variableName;
}

///////////////////////////////////////////////////////////////////////
// unit tests
///////////////////////////////////////////////////////////////////////
void FIELD_3D::cubicUnitTest(Real& ground, Real& computed, const Real x)
{
  //-0.352667 2.08522 -0.542286 0.128062
  Real yn1 = -0.3526673018932343;
  Real y0 =  2.085222005844116;
  Real y1 =  -0.542286217212677;
  Real y2 =  0.1280621886253357;

  Real dy0 = (y1 - yn1) * 0.5;
  Real dy1 = (y2 - y0) * 0.5;

  // compute the coefficients
  Real a = 2.0 * (y0 - y1) + dy0 - dy1;
  Real b = dy1 - 2.0 * dy0 + 3.0 * (y1 - y0);
  Real c = dy0;
  Real d = y0;

  // compute the ground truth
  ground = a * x * x * x + b * x * x + c * x + d;
  cout << " ground: " << ground << endl;

  // compute grid points
  Real points[] = {yn1, y0, y1, y2};

  computed = cubicInterp(x, points);

  cout << " computed: " << computed << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_3D::readHoudini(string filename)
{
  FILE* file;
  file = fopen(filename.c_str(), "r");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " FIELD_3D read failed! " << endl;
    cout << " Could not open file " << filename.c_str() << endl;
    exit(0);
  }

  bool volumeFound = false;
  char buffer[256];
  string search("Volume");
  while (!volumeFound)
  {
    fscanf(file, "%s", buffer);
    if (search.compare(buffer) == 0)
      volumeFound = true;
  }

  // read off the collision, vel.x, vel.y, and vel.z volumes
  readHoudiniField(file, true);
  readHoudiniField(file, false);
  readHoudiniField(file, false);
  readHoudiniField(file, false);
  
  // get the distance field
  readHoudiniField(file, false);

  fclose(file);
}
	  
///////////////////////////////////////////////////////////////////////
// read in a triplet of fields for velocity
///////////////////////////////////////////////////////////////////////
void FIELD_3D::readHoudiniVel(const string filename, FIELD_3D& xVelocity, FIELD_3D& yVelocity, FIELD_3D& zVelocity)
{
  FILE* file;
  file = fopen(filename.c_str(), "r");
  if (file == NULL)
  {
	  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
	  cout << " FIELD_3D read failed! " << endl;
	  cout << " Could not open file " << filename.c_str() << endl;
	  exit(0);
  }
  
  bool volumeFound = false;
  char buffer[256];
  string search("Volume");
  while (!volumeFound)
  {
	  fscanf(file, "%s", buffer);
	  if (search.compare(buffer) == 0)
		  volumeFound = true;
  }

  xVelocity = FIELD_3D::readHoudiniField(file);
  yVelocity = FIELD_3D::readHoudiniField(file);
  zVelocity = FIELD_3D::readHoudiniField(file);

  fclose(file);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_3D::readHoudiniSurf(string filename)
{
  int size = filename.size();
  if (filename[size - 1] == 'z' && filename[size - 2] == 'g')
  {
    readHoudiniSurfGz(filename);
    return;
  }

  FILE* file;
  file = fopen(filename.c_str(), "r");
  if (file == NULL)
  {
	  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
	  cout << " FIELD_3D read failed! " << endl;
	  cout << " Could not open file " << filename.c_str() << endl;
	  exit(0);
  }
  
  bool volumeFound = false;
  char buffer[256];
  string search("Volume");
  while (!volumeFound)
  {
	  fscanf(file, "%s", buffer);
	  if (search.compare(buffer) == 0)
		  volumeFound = true;
  }
  
  // get the distance field
  readHoudiniField(file, true);
  fclose(file);
}
	  
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_3D::readHoudiniSurfGz(string filename)
{
  cout << " Reading Houdini file ... "; flush(cout);

  // unzip the file
  string gunzip = string("gunzip ") + filename + string(".gz");
  cout << " Unzipping ... "; flush(cout);
  system(gunzip.c_str());

  readHoudiniSurf(filename);

  string gzip = string("gzip ") + filename;
  cout << " Re-zipping ... "; flush(cout);
  system(gzip.c_str());
  cout << " done. " << endl;
}

///////////////////////////////////////////////////////////////////////
// read a Houdini field off from a file stream -- it is assumed that the
// file is already advanced to the beginning of the field
///////////////////////////////////////////////////////////////////////
void FIELD_3D::readHoudiniField(FILE* file, bool storeValues)
{
  char buffer[256];
  //int index;
  VEC3F v0, v1, v2;
  int xRes, yRes, zRes;
  //Real negativeTwo;
#ifdef SINGLE_PRECISION
  fscanf(file, "%i %f %f %f %f %f %f %f %f %f %f %i %i %i", &index, 
      &v0[0], &v0[1], &v0[2],
      &v1[0], &v1[1], &v1[2],
      &v2[0], &v2[1], &v2[2],
      &negativeTwo,
      &xRes, &yRes, &zRes);
#else
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " Need to implement this for double precision! " << endl;
  exit(0);
#endif
	_lengths[0] = v0[0];
	_lengths[1] = v1[1];
	_lengths[2] = v2[2];	
  // initialize things if we aren't going to throw the field away
  if (storeValues)
  {
    if (_data) delete[] _data;
    _xRes = xRes; 
    _yRes = yRes; 
    _zRes = zRes; 

    _totalCells = _xRes * _yRes * _zRes;
    _slabSize = _xRes * _yRes;
    _data = new Real[_totalCells];

    _dx = _lengths[0] / _xRes;
    _dy = _lengths[1] / _yRes;
    _dz = _lengths[2] / _zRes;
    _invDx = 1.0 / _dx;
    _invDy = 1.0 / _dy;
    _invDz = 1.0 / _dz;

    _outside = maxRes() * maxRes();
  }

  // read in "streak"
  fscanf(file, "%s", buffer);

  // read in two zeros
  int zero1, zero2;
  fscanf(file, "%i %i", &zero1, &zero2);

  // read in "iso" or "invisible"
  fscanf(file, "%s", buffer);

  // read in zero and ten
  fscanf(file, "%i %i", &zero1, &zero2);

  if (storeValues)
  {
	  for (int y = 0; y < xRes * yRes * zRes; y++){
#ifdef SINGLE_PRECISION
      fscanf(file, "%f", &_data[y]);
#else
      {
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        cout << " Need to implement this for double precision! " << endl;
        exit(0);
      }
#endif
	  }
  }
  else
  {
    //Real data;
	  for (int y = 0; y < xRes * yRes * zRes; y++){
#ifdef SINGLE_PRECISION
      fscanf(file, "%f", &data);
#else
      {
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        cout << " Need to implement this for double precision! " << endl;
        exit(0);
      }
#endif
	  }
  }

  // chomp [0 4]
  fscanf(file, "%s", buffer);
  fscanf(file, "%s", buffer);
}

///////////////////////////////////////////////////////////////////////
// read a Houdini field off from a file stream -- it is assumed that the
// file is already advanced to the beginning of the field
//
// a static version for the velocity field
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::readHoudiniField(FILE* file)
{
  char buffer[256];
  //int index;
  VEC3F v0, v1, v2;
  int xRes, yRes, zRes;
  //Real negativeTwo;
#ifdef SINGLE_PRECISION
  fscanf(file, "%i %f %f %f %f %f %f %f %f %f %f %i %i %i", &index, 
      &v0[0], &v0[1], &v0[2],
      &v1[0], &v1[1], &v1[2],
      &v2[0], &v2[1], &v2[2],
      &negativeTwo,
      &xRes, &yRes, &zRes);
#else
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " Need to implement this for double precision! " << endl;
  exit(0);
#endif

  VEC3F lengths(v0[0], v1[1], v2[2]);
  VEC3F center;

  // initialize things if we aren't going to throw the field away
  FIELD_3D final(xRes, yRes, zRes, center, lengths);

  // read in "streak"
  fscanf(file, "%s", buffer);

  // read in two zeros
  int zero1, zero2;
  fscanf(file, "%i %i", &zero1, &zero2);

  // read in "iso" or "invisible"
  fscanf(file, "%s", buffer);

  // read in zero and ten
  fscanf(file, "%i %i", &zero1, &zero2);

  //Real* data = final.data();
	for (int y = 0; y < xRes * yRes * zRes; y++)
  {
#ifdef SINGLE_PRECISION
    fscanf(file, "%f", &data[y]);
#else
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Need to implement this for double precision! " << endl;
    exit(0);
#endif
  }

  // chomp [0 4]
  fscanf(file, "%s", buffer);
  fscanf(file, "%s", buffer);

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::Dx(const VEC3F& point) const
{
  vector<VEC3I> indices;
  vector<Real> weights;

  cellNeighborhood(point, indices, weights);

  assert(indices.size() == weights.size());

  Real final = 0;
  for (unsigned int x = 0; x < weights.size(); x++)
    final += weights[x] * Dx(indices[x][0], indices[x][1], indices[x][2]);

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::Dy(const VEC3F& point) const
{
  vector<VEC3I> indices;
  vector<Real> weights;

  cellNeighborhood(point, indices, weights);

  assert(indices.size() == weights.size());

  Real final = 0;
  for (unsigned int x = 0; x < weights.size(); x++)
    final += weights[x] * Dy(indices[x][0], indices[x][1], indices[x][2]);

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::Dz(const VEC3F& point) const
{
  vector<VEC3I> indices;
  vector<Real> weights;

  cellNeighborhood(point, indices, weights);

  assert(indices.size() == weights.size());

  Real final = 0;
  for (unsigned int x = 0; x < weights.size(); x++)
    final += weights[x] * Dz(indices[x][0], indices[x][1], indices[x][2]);

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::Dx(int x, int y, int z) const
{
  assert(x >= 0);
  assert(x < _xRes);
  assert(y >= 0);
  assert(y < _yRes);
  assert(z >= 0);
  assert(z < _zRes);

  int index = x + y * _xRes + z * _slabSize;

  const Real right = (x < _xRes - 1) ? _data[index + 1] : _data[index];
  const Real left  = (x > 0)         ? _data[index - 1] : _data[index];
  const Real denom = (x > 0 && x < _xRes -1) ? 1.0 / (2.0 * _dx) : 1.0 / _dx;
  return (right - left) * denom;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::DDx(int x, int y, int z) const
{
  assert(x >= 0);
  assert(x < _xRes);
  assert(y >= 0);
  assert(y < _yRes);
  assert(z >= 0);
  assert(z < _zRes);

  int index = x + y * _xRes + z * _slabSize;

  const Real right = (x < _xRes - 1) ? _data[index + 1] : _data[index];
  const Real left  = (x > 0)         ? _data[index - 1] : _data[index];
  const Real denom = 1.0 / (_dx * _dx);
  return (right - 2.0 * _data[index] + left) * denom;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::Dy(int x, int y, int z) const
{
  assert(x >= 0);
  assert(x < _xRes);
  assert(y >= 0);
  assert(y < _yRes);
  assert(z >= 0);
  assert(z < _zRes);

  int index = x + y * _xRes + z * _slabSize;

  const Real up   = (y < _yRes - 1) ? _data[index + _xRes] : _data[index];
  const Real down = (y > 0)         ? _data[index - _xRes] : _data[index];
  const Real denom = (y > 0 && y < _yRes -1) ? 1.0 / (2.0 * _dy) : 1.0 / _dy;
  return (up - down) * denom;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::DDy(int x, int y, int z) const
{
  assert(x >= 0);
  assert(x < _xRes);
  assert(y >= 0);
  assert(y < _yRes);
  assert(z >= 0);
  assert(z < _zRes);

  int index = x + y * _xRes + z * _slabSize;

  const Real up   = (y < _yRes - 1) ? _data[index + _xRes] : _data[index];
  const Real down = (y > 0)         ? _data[index - _xRes] : _data[index];
  const Real denom = 1.0 / (_dy * _dy);
  return (up - 2.0 * _data[index] + down) * denom;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::Dz(int x, int y, int z) const
{
  assert(x >= 0);
  assert(x < _xRes);
  assert(y >= 0);
  assert(y < _yRes);
  assert(z >= 0);
  assert(z < _zRes);

  int index = x + y * _xRes + z * _slabSize;

  const Real in  = (z < _zRes - 1) ? _data[index + _slabSize] : _data[index];
  const Real out = (z > 0)         ? _data[index - _slabSize] : _data[index];
  const Real denom = (z > 0 && z < _zRes -1) ? 1.0 / (2.0 * _dz) : 1.0 / _dz;
  return (in - out) * denom;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::DDz(int x, int y, int z) const
{
  assert(x >= 0);
  assert(x < _xRes);
  assert(y >= 0);
  assert(y < _yRes);
  assert(z >= 0);
  assert(z < _zRes);

  int index = x + y * _xRes + z * _slabSize;

  const Real in  = (z < _zRes - 1) ? _data[index + _slabSize] : _data[index];
  const Real out = (z > 0)         ? _data[index - _slabSize] : _data[index];
  const Real denom = 1.0 / (_dz * _dz);
  return (in - out) * denom;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::DDxy(int x, int y, int z) const
{
  assert(x >= 0);
  assert(x < _xRes);
  assert(y >= 0);
  assert(y < _yRes);
  assert(z >= 0);
  assert(z < _zRes);

  int index = x + y * _xRes + z * _slabSize;

  // if it's on the border, just use the center cell
  const int xPlus  = (x < _xRes - 1) ? 1 : 0;
  const int xMinus = (x > 0) ? -1 : 0;
  const int yPlus  = (y < _yRes - 1) ? _xRes : 0;
  const int yMinus = (y > 0) ? -_xRes : 0;

  const Real plusPlus   = _data[index + xPlus + yPlus];
  const Real plusMinus  = _data[index + xPlus + yMinus];
  const Real minusPlus  = _data[index + xMinus + yPlus];
  const Real minusMinus = _data[index + xMinus + yMinus];

  return (plusPlus - plusMinus - minusPlus + minusMinus) / (_dx * _dy);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::DDxz(int x, int y, int z) const
{
  assert(x >= 0);
  assert(x < _xRes);
  assert(y >= 0);
  assert(y < _yRes);
  assert(z >= 0);
  assert(z < _zRes);

  int index = x + y * _xRes + z * _slabSize;

  // if it's on the border, just use the center cell
  const int xPlus  = (x < _xRes - 1) ? 1 : 0;
  const int xMinus = (x > 0) ? -1 : 0;
  const int zPlus  = (z < _zRes - 1) ? _slabSize : 0;
  const int zMinus = (z > 0) ? -_slabSize : 0;

  const Real plusPlus   = _data[index + xPlus + zPlus];
  const Real plusMinus  = _data[index + xPlus + zMinus];
  const Real minusPlus  = _data[index + xMinus + zPlus];
  const Real minusMinus = _data[index + xMinus + zMinus];

  return (plusPlus - plusMinus - minusPlus + minusMinus) / (_dx * _dz);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::DDyz(int x, int y, int z) const
{
  assert(x >= 0);
  assert(x < _xRes);
  assert(y >= 0);
  assert(y < _yRes);
  assert(z >= 0);
  assert(z < _zRes);

  int index = x + y * _xRes + z * _slabSize;

  // if it's on the border, just use the center cell
  const int yPlus  = (y < _yRes - 1) ? _xRes : 0;
  const int yMinus = (y > 0) ? -_xRes : 0;
  const int zPlus  = (z < _zRes - 1) ? _slabSize : 0;
  const int zMinus = (z > 0) ? -_slabSize : 0;

  const Real plusPlus   = _data[index + yPlus + zPlus];
  const Real plusMinus  = _data[index + yPlus + zMinus];
  const Real minusPlus  = _data[index + yMinus + zPlus];
  const Real minusMinus = _data[index + yMinus + zMinus];

  return (plusPlus - plusMinus - minusPlus + minusMinus) / (_dy * _dz);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::Dx() const
{
  FIELD_3D final(_xRes, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        final(x,y,z) = Dx(x,y,z);

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::Dy() const
{
  FIELD_3D final(_xRes, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        final(x,y,z) = Dy(x,y,z);

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::Dz() const
{
  FIELD_3D final(_xRes, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        final(x,y,z) = Dz(x,y,z);

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::DDx() const
{
  FIELD_3D final(_xRes, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        final(x,y,z) = DDx(x,y,z);

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::DDy() const
{
  FIELD_3D final(_xRes, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        final(x,y,z) = DDy(x,y,z);

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::DDz() const
{
  FIELD_3D final(_xRes, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        final(x,y,z) = DDz(x,y,z);

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::DDxy() const
{
  FIELD_3D final(_xRes, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        final(x,y,z) = DDxy(x,y,z);

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::DDxz() const
{
  FIELD_3D final(_xRes, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        final(x,y,z) = DDxz(x,y,z);

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::DDyz() const
{
  FIELD_3D final(_xRes, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        final(x,y,z) = DDyz(x,y,z);

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the mean curvature
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::meanCurvature() const
{
  FIELD_3D final(_xRes, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        Real px = Dx(x,y,z);
        Real py = Dy(x,y,z);
        Real pz = Dz(x,y,z);
        
        Real pxx = DDx(x,y,z);
        Real pyy = DDy(x,y,z);
        Real pzz = DDz(x,y,z);
        
        Real pxy = DDxy(x,y,z);
        Real pxz = DDxz(x,y,z);
        Real pyz = DDyz(x,y,z);

        final(x,y,z) = (pyy + pzz) * px * px +
                       (pxx + pzz) * py * py +
                       (pxx + pyy) * pz * pz - 
                       2.0 * (px * py * pxy + px * pz * pxz + py * pz * pyz);

        Real denom = pow(px + py + pz, (Real)1.5);
        denom = (fabs(denom) > 1e-6) ? 1.0 / denom : 0;

        final(x,y,z) *= denom;
      }

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the Gaussian curvature
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::gaussianCurvature() const
{
  FIELD_3D final(_xRes, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        Real px = Dx(x,y,z);
        Real py = Dy(x,y,z);
        Real pz = Dz(x,y,z);
        
        Real pxx = DDx(x,y,z);
        Real pyy = DDy(x,y,z);
        Real pzz = DDz(x,y,z);
        
        Real pxy = DDxy(x,y,z);
        Real pxz = DDxz(x,y,z);
        Real pyz = DDyz(x,y,z);

        final(x,y,z) = (pyy + pzz - pyz * pyz) * px * px +
                       (pxx + pzz - pxz * pxz) * py * py +
                       (pxx + pyy - pxy * pxy) * pz * pz + 
                       2.0 * (px * py * (pxz * pyz - pxy * pzz) + 
                              py * pz * (pxy * pxz - pyz * pxx) +
                              px * pz * (pxy * pyz - pxz * pyy));

        Real denom = pow(px + py + pz, (Real)2.0);
        denom = (fabs(denom) > 1e-6) ? 1.0 / denom : 0;

        final(x,y,z) *= denom;
      }

  return final;
}

///////////////////////////////////////////////////////////////////////
// mask out any values past a certain distance
///////////////////////////////////////////////////////////////////////
void FIELD_3D::maskByDistance(const FIELD_3D& distanceField, const Real distance)
{
  assert(_xRes == distanceField.xRes());
  assert(_yRes == distanceField.yRes());
  assert(_zRes == distanceField.zRes());

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        if (fabs(distanceField(x,y,z)) > distance)
          (*this)(x,y,z) = 0;
      }
}

///////////////////////////////////////////////////////////////////////
// clamp the field to a min and max
///////////////////////////////////////////////////////////////////////
void FIELD_3D::clamp(const Real minValue, const Real maxValue)
{
  assert(minValue < maxValue);

  for (int x = 0; x < _totalCells; x++)
  {
    if (_data[x] < minValue)
      _data[x] = minValue;
    if (_data[x] > maxValue)
      _data[x] = maxValue;
  }
}

///////////////////////////////////////////////////////////////////////
// get a resampled version
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::resampleCubic(int xRes, int yRes, int zRes) const
{
  TIMER functionTimer(__FUNCTION__);
  FIELD_3D final(xRes, yRes, zRes, _center, _lengths);

  const int slabSize = xRes * yRes;
#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        int index = x + y * xRes + z * slabSize;
        VEC3F center = final.cellCenter(x,y,z);
        final[index] = this->cubicLookup(center);
      }

  return final;
}

///////////////////////////////////////////////////////////////////////
// get a resampled version
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::resampleQuintic(int xRes, int yRes, int zRes) const
{
  TIMER functionTimer(__FUNCTION__);
  FIELD_3D final(xRes, yRes, zRes, _center, _lengths);

  int index = 0;
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++, index++)
      {
        VEC3F center = final.cellCenter(x,y,z);
        final[index] = this->quinticLookup(center);
      }

  return final;
}

///////////////////////////////////////////////////////////////////////
// take the square root of the field
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::squareRoot()
{
  TIMER functionTimer(__FUNCTION__);
  FIELD_3D final(*this);

  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, index++)
        final[index] = sqrt(_data[index]);

  return final;
}

///////////////////////////////////////////////////////////////////////
// from "Level Set Surface Editing Operators", Museth et al. 2002
// "Geometric Surface Processing via Normal Maps", Tasdizen 2003
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::principalCurvature() const
{
  FIELD_3D final(*this);

  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, index++)
      {
        MATRIX3 N = Dnormal(x,y,z);
        VEC3F n = normal(x,y,z);
        MATRIX3 outer = MATRIX3::outer_product(n, n);

        MATRIX3 B = N * (MATRIX3::I() - outer);

        Real D = sqrt(B.squaredSum());
        Real H = trace(B) * 0.5;

        Real discrim = D * D * 0.5 - H * H;

        if (discrim < 0.0)
        {
          final[index] = 0;
          continue;
        }

        Real root = sqrt(discrim);
        Real k1 = H + root;
        Real k2 = H - root;

        final[index] = (fabs(k1) > fabs(k2)) ? k1 : k2;
      }

  return final;
}

///////////////////////////////////////////////////////////////////////
// from "Level Set Surface Editing Operators", Museth et al. 2002
// "Geometric Surface Processing via Normal Maps", Tasdizen 2003
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::principalCurvature(const VEC3F& point) const
{
  vector<VEC3I> indices;
  vector<Real> weights;

  cellNeighborhood(point, indices, weights);

  Real final = 0;
  //int index = 0;
  for (unsigned int i = 0; i < indices.size(); i++)
  {
    VEC3I index = indices[i];

    MATRIX3 N = Dnormal(index[0], index[1], index[2]);
    VEC3F n = normal(index[0], index[1], index[2]);
    MATRIX3 outer = MATRIX3::outer_product(n, n);

    MATRIX3 B = N * (MATRIX3::I() - outer);

    Real D = sqrt(B.squaredSum());
    Real H = trace(B) * 0.5;

    Real discrim = D * D * 0.5 - H * H;

    if (discrim < 0.0)
      continue;

    Real root = sqrt(discrim);
    Real k1 = H + root;
    Real k2 = H - root;

    Real curvature = (fabs(k1) > fabs(k2)) ? k1 : k2;
    final += weights[i] * curvature;
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the normal at a point
///////////////////////////////////////////////////////////////////////
VEC3F FIELD_3D::normal(int x, int y, int z) const
{
  assert (x >= 0 && x < _xRes);
  assert (y >= 0 && y < _yRes);
  assert (z >= 0 && z < _zRes);

  VEC3F final;
  final[0] = Dx(x,y,z);
  final[1] = Dy(x,y,z);
  final[2] = Dz(x,y,z);

  final.normalize();

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the normal at a point
///////////////////////////////////////////////////////////////////////
VEC3F FIELD_3D::normal(const VEC3F& point) const
{
  VEC3F final;
  final[0] = Dx(point);
  final[1] = Dy(point);
  final[2] = Dz(point);

  final.normalize();

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the normal at a point
///////////////////////////////////////////////////////////////////////
MATRIX3 FIELD_3D::Dnormal(int x, int y, int z) const
{
  const VEC3F left  = (x == 0)         ? normal(x,y,z) : normal(x-1,y,z);
  const VEC3F right = (x == _xRes - 1) ? normal(x,y,z) : normal(x+1,y,z);
  const Real dx     = (x == 0 || x == _xRes - 1) ? 1.0 / _dx : 0.5 / _dx;

  const VEC3F down = (y == 0)         ? normal(x,y,z) : normal(x,y-1,z);
  const VEC3F up   = (y == _yRes - 1) ? normal(x,y,z) : normal(x,y+1,z);
  const Real dy    = (y == 0 || y == _yRes - 1) ? 1.0 / _dy : 0.5 / _dy;

  const VEC3F out = (z == 0)         ? normal(x,y,z) : normal(x,y,z-1);
  const VEC3F in  = (z == _zRes - 1) ? normal(x,y,z) : normal(x,y,z+1);
  const Real dz   = (z == 0 || z == _zRes - 1) ? 1.0 / _dz : 0.5 / _dz;

  MATRIX3 final((right - left) * dx, (up - down) * dy, (in - out) * dz);

  return final;
}

///////////////////////////////////////////////////////////////////////
// BLAS-like interface, output += alpha * input
///////////////////////////////////////////////////////////////////////
void FIELD_3D::axpy(const Real& alpha, const FIELD_3D& input, FIELD_3D& output)
{
  assert(input.xRes() == output.xRes());
  assert(input.yRes() == output.yRes());
  assert(input.zRes() == output.zRes());

  int totalCells = input.totalCells();
  for (int x = 0; x < totalCells; x++)
    output[x] += alpha * input[x];
}

///////////////////////////////////////////////////////////////////////
// BLAS-like interface, output += alpha * input
///////////////////////////////////////////////////////////////////////
void FIELD_3D::axpy(const Real& alpha, const FIELD_3D& input)
{
  assert(input.xRes() == _xRes);
  assert(input.yRes() == _yRes);
  assert(input.zRes() == _zRes);

  int totalCells = input.totalCells();
  for (int x = 0; x < totalCells; x++)
    _data[x] += alpha * input[x];
}

///////////////////////////////////////////////////////////////////////
// clamp the field to a min and max
///////////////////////////////////////////////////////////////////////
void FIELD_3D::bandPass(const Real minValue, const Real maxValue)
{
  assert(minValue < maxValue);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        if ((*this)(x,y,z) < minValue)
          (*this)(x,y,z) = 0;
        if ((*this)(x,y,z) > maxValue)
          (*this)(x,y,z) = 0;
      }
}

///////////////////////////////////////////////////////////////////////
// isolate values near the current value, within a certain width
///////////////////////////////////////////////////////////////////////
void FIELD_3D::isolateBand(const Real target, const Real width)
{
  int unstomped = 0;

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        Real value = fabs((*this)(x,y,z));

        if (value > target + width || value < target - width)
          (*this)(x,y,z) = 0;
        else
          unstomped++;
      }
}

///////////////////////////////////////////////////////////////////////
// single explicit diffusion step
///////////////////////////////////////////////////////////////////////
void FIELD_3D::blur(const Real dt)
{
  FIELD_3D temp(*this);

  for (int z = 1; z < _zRes - 1; z++)
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++)
      {
        int index = x + y * _xRes + z * _slabSize;
        Real mean = -6.0 * _data[index] + 
                    _data[index + 1] + _data[index - 1] +
                    _data[index + _xRes] + _data[index - _xRes] +
                    _data[index + _slabSize] + _data[index - _slabSize];

        temp(x,y,z) = _data[index] + dt * mean;
      }

  *this = temp;
}

///////////////////////////////////////////////////////////////////////
// get band-limited curvature, with some diffusion
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::bandLimitedCurvature(const Real target, const Real width)
{
  FIELD_3D final(*this);

  final = this->principalCurvature();

  final.maskByDistance(*this, _dx * 5);

  final.isolateBand(target, width);

  return final;
}

///////////////////////////////////////////////////////////////////////
// set to the absolute value
///////////////////////////////////////////////////////////////////////
void FIELD_3D::absoluteValue()
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] = fabs(_data[x]);
}

///////////////////////////////////////////////////////////////////////
// do a soft bandpass where there's a gradual cubic falloff
///////////////////////////////////////////////////////////////////////
void FIELD_3D::softBandPass(const Real band, const Real falloff)
{
  for (int x = 0; x < _totalCells; x++)
    if (fabs(fabs(_data[x]) - band) > falloff)
      _data[x] = 0;
    else
    {
      Real interp = fabs(fabs(_data[x]) - band) / falloff;
      Real squared = interp * interp;
      Real scaling = 2 * squared * interp - 3 * squared + 1;
      _data[x] *= scaling;
    }
}

///////////////////////////////////////////////////////////////////////
// write out a field to a file stream
///////////////////////////////////////////////////////////////////////
void FIELD_3D::write(FILE* file) const
{
  // write dimensions
  fwrite((void*)&_xRes, sizeof(int), 1, file);
  fwrite((void*)&_yRes, sizeof(int), 1, file);
  fwrite((void*)&_zRes, sizeof(int), 1, file);
  _center.write(file);
  _lengths.write(file);

  if (_totalCells <= 0)
    return;

  // always write out as a double
  if (sizeof(Real) != sizeof(double))
  {
    double* dataDouble = new double[_totalCells];
    for (int x = 0; x < _totalCells; x++)
      dataDouble[x] = _data[x];

    fwrite((void*)dataDouble, sizeof(double), _totalCells, file);
    delete[] dataDouble;
  }
  else
    fwrite((void*)_data, sizeof(Real), _totalCells, file);
}

///////////////////////////////////////////////////////////////////////
// read in a field from a file stream
///////////////////////////////////////////////////////////////////////
void FIELD_3D::read(FILE* file)
{
  // read dimensions
  fread((void*)&_xRes, sizeof(int), 1, file);
  fread((void*)&_yRes, sizeof(int), 1, file);
  fread((void*)&_zRes, sizeof(int), 1, file);
  _center.read(file);
  _lengths.read(file);
  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
  _invDz = 1.0 / _dz;
  _outside = maxRes() * maxRes();

  _slabSize = _xRes * _yRes;
  _totalCells = _xRes * _yRes * _zRes;

  if (_xRes < 0 || _yRes < 0 || _zRes < 0)
  {
    _data = NULL;
    return;
  }

  if (_data) delete[] _data;
  try {
    _data = new Real[_totalCells];
  }
  catch(std::bad_alloc& exc)
  {
    cout << " Failed to allocate " << _xRes << " " << _yRes << " " << _zRes << " FIELD_3D!" << endl;
    double bytes = _xRes * _yRes * _zRes * sizeof(Real);
    cout <<  bytes / pow(2.0,20.0) << " MB needed" << endl;
    exit(0);
  }

  // always read in as a double
  if (sizeof(Real) != sizeof(double))
  {
    double* dataDouble = new double[_totalCells];
    fread((void*)dataDouble, sizeof(double), _totalCells, file);

    for (int x = 0; x < _totalCells; x++)
      _data[x] = dataDouble[x];

    delete[] dataDouble;
  }
  else
    fread((void*)_data, sizeof(Real), _totalCells, file);
}

///////////////////////////////////////////////////////////////////////
// see if the field has any data in it yet
///////////////////////////////////////////////////////////////////////
bool FIELD_3D::initialized() const
{
  if (_xRes < 0 || _yRes < 0 || _zRes < 0 || _totalCells < 0 || _data == NULL)
    return false;

  return true;
}

///////////////////////////////////////////////////////////////////////
// write out a field to a file stream
///////////////////////////////////////////////////////////////////////
void FIELD_3D::writeGz(gzFile& file) const
{
  // write dimensions
  gzwrite(file, (void*)&_xRes, sizeof(int));
  gzwrite(file, (void*)&_yRes, sizeof(int));
  gzwrite(file, (void*)&_zRes, sizeof(int));
  _center.writeGz(file);
  _lengths.writeGz(file);

  // always write out as a double
  if (sizeof(Real) != sizeof(double))
  {
    double* dataDouble = new double[_totalCells];
    for (int x = 0; x < _totalCells; x++)
      dataDouble[x] = _data[x];

    gzwrite(file, (void*)dataDouble, _totalCells * sizeof(double));
    delete[] dataDouble;
  }
  else
    gzwrite(file, (void*)_data, _totalCells * sizeof(Real));
}

///////////////////////////////////////////////////////////////////////
// read in a field from a file stream
///////////////////////////////////////////////////////////////////////
void FIELD_3D::readGz(gzFile& file)
{
  // read dimensions
  gzread(file, (void*)&_xRes, sizeof(int));
  gzread(file, (void*)&_yRes, sizeof(int));
  gzread(file, (void*)&_zRes, sizeof(int));
  _center.readGz(file);
  _lengths.readGz(file);
  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
  _invDz = 1.0 / _dz;
  _outside = maxRes() * maxRes();

  _slabSize = _xRes * _yRes;
  _totalCells = _xRes * _yRes * _zRes;
  if (_data) delete[] _data;
  _data = new Real[_totalCells];

  // always read in as a double
  if (sizeof(Real) != sizeof(double))
  {
    double* dataDouble = new double[_totalCells];
    gzread(file, (void*)dataDouble, _totalCells * sizeof(double));

    for (int x = 0; x < _totalCells; x++)
      _data[x] = dataDouble[x];

    delete[] dataDouble;
  }
  else
    gzread(file, (void*)_data, _totalCells * sizeof(Real));
}

///////////////////////////////////////////////////////////////////////
// copy values out into the border, assuming that "borderSize" is the 
// width of the grid padding
///////////////////////////////////////////////////////////////////////
void FIELD_3D::copyIntoBorder(int borderSize)
{
  TIMER functionTimer(__FUNCTION__);
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        if (x == borderSize)
        {
          Real value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x - i,y,z) = value;
        }
        if (x == _xRes - 1 - borderSize)
        {
          Real value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x + i,y,z) = value;
        }              
        if (y == borderSize)
        {
          Real value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x,y - i,z) = value;
        }
        if (y == _yRes - 1 - borderSize)
        {
          Real value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x,y+i,z) = value;
        }              
        if (z == borderSize)
        {
          Real value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x,y,z - i) = value;
        }
        if (z == _zRes - 1 - borderSize)
        {
          Real value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x,y,z+i) = value;
        }

        if (x == borderSize && y == borderSize)
        {
          Real value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x - i,y - j,z) = value;
        }
        if (x == _xRes - 1 - borderSize && y == _yRes - 1 - borderSize)
        {
          Real value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x + i,y + j,z) = value;
        }

        // handle the corners
        if (x == borderSize && z == borderSize)
        {
          Real value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x - i,y,z-j) = value;
        }
        if (x == _xRes - 1 - borderSize && z == _zRes - 1 - borderSize)
        {
          Real value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x + i,y,z+j) = value;
        }

        if (z == borderSize && y == borderSize)
        {
          Real value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x,y - j,z -i) = value;
        }
        if (z == _xRes - 1 - borderSize && y == _yRes - 1 - borderSize)
        {
          Real value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x,y + j,z+i) = value;
        }

      }
}

///////////////////////////////////////////////////////////////////////
// pass back a field with a new padding of size "paddingSize"
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::withAddedPadding(int paddingSize) const
{
  // new length, with padding
  VEC3F newLengths = _lengths;
  newLengths[0] += paddingSize * 2 * _dx;
  newLengths[1] += paddingSize * 2 * _dx;
  newLengths[2] += paddingSize * 2 * _dx;

  FIELD_3D final(_xRes + 2 * paddingSize, 
                 _yRes + 2 * paddingSize, 
                 _zRes + 2 * paddingSize, _center, newLengths);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        final(x + paddingSize,
              y + paddingSize,
              z + paddingSize) = (*this)(x,y,z);

  final.copyIntoBorder(paddingSize);

  return final;
}

///////////////////////////////////////////////////////////////////////
// stomp the border to zero
///////////////////////////////////////////////////////////////////////
void FIELD_3D::stompBorder(int borderSize)
{
  TIMER functionTimer(__FUNCTION__);
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        if (x < borderSize || y < borderSize || z < borderSize ||
            x >= _xRes - borderSize || y >= _yRes - borderSize || z >= _zRes - borderSize)
          (*this)(x,y,z) = 0;
      }
}

///////////////////////////////////////////////////////////////////////
// set a border of size 1 to zero
///////////////////////////////////////////////////////////////////////
void FIELD_3D::setZeroBorder()
{
	int index;
	for (int z = 0; z < _zRes; z++)
		for (int y = 0; y < _yRes; y++)
		{
			// left slab
			index = y * _xRes + z * _slabSize;
			_data[index] = 0.0f;

			// right slab
			index += _xRes - 1;
			_data[index] = 0.0f;
		}
	for (int z = 0; z < _zRes; z++)
		for (int x = 0; x < _xRes; x++)
		{
			// bottom slab
			index = x + z * _slabSize;
			_data[index] = 0.0f;

			// top slab
			index += _slabSize - _xRes;
			_data[index] = 0.0f;
		}
	for (int y = 0; y < _yRes; y++)
		for (int x = 0; x < _xRes; x++)
		{
			// front slab
			index = x + y * _xRes;
			_data[index] = 0.0f;

			// back slab
			index += _totalCells - _slabSize;
			_data[index] = 0.0f;
		}
}

//////////////////////////////////////////////////////////////////////
// swap the contents with another object
//////////////////////////////////////////////////////////////////////
void FIELD_3D::swapPointers(FIELD_3D& field)
{
  assert(field.xRes() == _xRes);
  assert(field.yRes() == _yRes);
  assert(field.zRes() == _zRes);
  
  Real* temp = _data;
  _data = field._data;
  field._data = temp;
}

//////////////////////////////////////////////////////////////////////
// take the dot product with respect to another field
//////////////////////////////////////////////////////////////////////
Real FIELD_3D::dot(const FIELD_3D& input) const
{
  assert(input.xRes() == _xRes);
  assert(input.yRes() == _yRes);
  assert(input.zRes() == _zRes);
  Real final = 0;

  for (int x = 0; x < _totalCells; x++)
    final += input[x] * _data[x];

  return final;
}

//////////////////////////////////////////////////////////////////////
// get a subfield of this field
//////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::subfield(const int xBegin, const int xEnd, 
                            const int yBegin, const int yEnd, 
                            const int zBegin, const int zEnd) const
{
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

  FIELD_3D final(xInterval, yInterval, zInterval);

  cout << " size: " << final.dims() << endl;

  for (int z = 0; z < zInterval; z++)
    for (int y = 0; y < yInterval; y++)
      for (int x = 0; x < xInterval; x++)
        final(x,y,z) = (*this)(xBegin + x, yBegin + y, zBegin + z);

  return final;
}

//////////////////////////////////////////////////////////////////////
// peel off the outer boundary of grid cells
//////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::peelBoundary() const
{
  FIELD_3D final(_xRes - 2, _yRes - 2, _zRes - 2);

  for (int z = 1; z < _zRes - 1; z++)
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++)
        final(x - 1, y - 1, z - 1) = (*this)(x,y,z);

  return final;
}

//////////////////////////////////////////////////////////////////////
// is this point outside the bounding box?
//////////////////////////////////////////////////////////////////////
bool FIELD_3D::inside(const VEC3F& point) const
{
  const VEC3F halfLengths = (Real)0.5 * _lengths;

  if (point[0] < _center[0] - halfLengths[0])
    return false;
  if (point[0] > _center[0] + halfLengths[0])
    return false;
  if (point[1] < _center[1] - halfLengths[1])
    return false;
  if (point[1] > _center[1] + halfLengths[1])
    return false;
  if (point[2] < _center[2] - halfLengths[2])
    return false;
  if (point[2] > _center[2] + halfLengths[2])
    return false;

  return true;
}

//////////////////////////////////////////////////////////////////////
// splat a density particle to this grid
//////////////////////////////////////////////////////////////////////
void FIELD_3D::splat(const Real weight, const VEC3F& point, const int cellRadius)
{
  VEC3F positionCopy = point;

  // get the lower corner position
  VEC3F corner = _center - (Real)0.5 * _lengths;
  VEC3F dxs(_dx, _dy, _dz);
  corner += (Real)0.5 * dxs;

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  int x0 = (int)positionCopy[0];
  int x1    = x0 + 1;
  int y0 = (int)positionCopy[1];
  int y1    = y0 + 1;
  int z0 = (int)positionCopy[2];
  int z1    = z0 + 1;

  x0 -= cellRadius;
  x1 += cellRadius;
  y0 -= cellRadius;
  y1 += cellRadius;
  z0 -= cellRadius;
  z1 += cellRadius;

  // clamp everything
  x0 = (x0 < 0) ? 0 : x0;
  y0 = (y0 < 0) ? 0 : y0;
  z0 = (z0 < 0) ? 0 : z0;
  
  x1 = (x1 < 0) ? 0 : x1;
  y1 = (y1 < 0) ? 0 : y1;
  z1 = (z1 < 0) ? 0 : z1;

  x0 = (x0 > _xRes - 1) ? _xRes - 1 : x0;
  y0 = (y0 > _yRes - 1) ? _yRes - 1 : y0;
  z0 = (z0 > _zRes - 1) ? _zRes - 1 : z0;

  x1 = (x1 > _xRes - 1) ? _xRes - 1 : x1;
  y1 = (y1 > _yRes - 1) ? _yRes - 1 : y1;
  z1 = (z1 > _zRes - 1) ? _zRes - 1 : z1;

  for (int z = z0; z <= z1; z++)
  {
    const int zStride = z * _slabSize;
    for (int y = y0; y <= y1; y++)
    {
      const int yStride = y * _xRes + zStride;
      for (int x = x0; x <= x1; x++)
      {
        VEC3F center = cellCenter(x,y,z);
        Real distance = fabs(norm(point - center) / ((Real)(1 + cellRadius) * _dx));
        int index = x + yStride;

        //_data[index] += weight * (3.0 * distance * distance +
        //                          2.0 * distance * distance * distance);
        //_data[index] += distance;
        distance = 1.0 - distance;
        distance = (distance < 0.0) ? 0.0 : distance;
        distance = (distance > 1.0) ? 1.0 : distance;

        Real final = weight * (3.0 * distance * distance +
                               2.0 * distance * distance * distance);

#pragma omp atomic
        _data[index] += final;
      }
    }
  }
}

static Real smooth_step(Real r)
{
   if(r<0) return 0;
   else if(r>1) return 1;
   return r*r*r*(10+r*(-15+r*6));
}
static Real ramp(Real r)
{
  return smooth_step((r+1)/2)*2-1;
}

//////////////////////////////////////////////////////////////////////
// compute a ramp function in the y direction
//////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::yRampField(const FIELD_3D& example, const Real plumeBase)
{
  FIELD_3D final(example);
  final.clear();

  int xRes = example.xRes();
  int yRes = example.yRes();
  int zRes = example.zRes();

  int biggest = (xRes > yRes) ? xRes : yRes;
  biggest = (biggest > zRes) ? biggest : zRes;

  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        //float xFrac = (float)x / biggest;
        float yFrac = (float)y / biggest;
        //float zFrac = (float)z / biggest;
        
        float plumeHeight = (1);
        float heightFactor= ramp((yFrac - plumeBase) / plumeHeight);

        if (heightFactor< 0) heightFactor = 0;
        
        final(x,y,z) = heightFactor;
      }
  
  return final;
}

//////////////////////////////////////////////////////////////////////
// get the sum of gradients
//////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::divergence() const
{
  FIELD_3D x = this->Dx();
  x.absoluteValue();
  FIELD_3D y = this->Dy();
  y.absoluteValue();
  FIELD_3D z = this->Dz();
  z.absoluteValue();

  return x + y + z;
}
  
//////////////////////////////////////////////////////////////////////
// get the Laplacian
//////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::laplacian() const
{
  FIELD_3D x = DDx();
  FIELD_3D y = DDy();
  FIELD_3D z = DDz();

  return x + y + z;
}

//////////////////////////////////////////////////////////////////////
// linearly downsampled -- halve the size of the field
//////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::linearHalfSample()
{
  assert(_xRes % 2 == 0);
  assert(_yRes % 2 == 0);
  assert(_zRes % 2 == 0);
  VEC3F lengths = _lengths;
  lengths[0] -= (lengths[0] / _xRes);
  lengths[1] -= (lengths[1] / _yRes);
  lengths[2] -= (lengths[2] / _zRes);

  int xHalf = _xRes / 2;
  int yHalf = _yRes / 2;
  int zHalf = _zRes / 2;

  FIELD_3D half(xHalf, yHalf, zHalf, _center, lengths);

  for (int z = 0; z < zHalf; z++)
    for (int y = 0; y < yHalf; y++)
      for (int x = 0; x < xHalf; x++)
      {
        Real sum = 0;
        sum += (*this)(2 * x    , 2 * y    , 2 * z);
        sum += (*this)(2 * x    , 2 * y    , 2 * z + 1);
        sum += (*this)(2 * x    , 2 * y + 1, 2 * z);
        sum += (*this)(2 * x    , 2 * y + 1, 2 * z + 1);
        sum += (*this)(2 * x + 1, 2 * y    , 2 * z);
        sum += (*this)(2 * x + 1, 2 * y    , 2 * z + 1);
        sum += (*this)(2 * x + 1, 2 * y + 1, 2 * z);
        sum += (*this)(2 * x + 1, 2 * y + 1, 2 * z + 1);
        sum *= 1.0 / 8.0;
        half(x,y,z) = sum;
      }
  return half;
}

//////////////////////////////////////////////////////////////////////
// downsample multigrid-style so that half the sample still remain 
// on the grid
//////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::multigridDownsample()
{
  assert(_xRes % 2 == 1);
  assert(_yRes % 2 == 1);
  assert(_zRes % 2 == 1);

  int xResNew = (int)(_xRes / 2) + 1;
  int yResNew = (int)(_yRes / 2) + 1;
  int zResNew = (int)(_zRes / 2) + 1;

  FIELD_3D final(xResNew, yResNew, zResNew, _center, _lengths);

  for (int z = 0; z < zResNew; z++)
    for (int y = 0; y < yResNew; y++)
      for (int x = 0; x < xResNew; x++)
        final(x,y,z) = (*this)(2 * x, 2 * y, 2 * z);

  return final;
}

//////////////////////////////////////////////////////////////////////
// take the power of a field
//////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::power(const Real& exponent) const
{
  FIELD_3D final = *this;

  for (int x = 0; x < _totalCells; x++)
    final[x] = pow(final[x], exponent);

  return final;
}

//////////////////////////////////////////////////////////////////////
// assuming that this field is a signed distance field, get the shell
// of cells that are immediately on the surface
//////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::outerShell() const
{
  FIELD_3D final(*this);
  final = 0;

#pragma omp parallel
#pragma omp for  schedule(static)
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        Real center = (*this)(x,y,z);
        if (center < 0) continue;

        bool flip = false;
        for (int i = -1; i <= 1; i++)
          for (int j = -1; j <= 1; j++)
            for (int k = -1; k <= 1; k++)
            {
              if (i == 0 && j == 0 && k == 0) continue;
              int xNew = x + i;
              int yNew = y + j;
              int zNew = z + k;

              if (xNew < 0 || xNew > _xRes - 1) continue;
              if (yNew < 0 || yNew > _yRes - 1) continue;
              if (zNew < 0 || zNew > _zRes - 1) continue;

              if ((*this)(xNew, yNew, zNew) * center < 0)
                flip = true;
            }
         final(x,y,z) = flip ? 1 : 0;
      }

  return final;
}

//////////////////////////////////////////////////////////////////////
// assuming that this field is a signed distance field, get the shell
// of cells that are immediately on the surface
//////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::innerShell() const
{
  FIELD_3D final(*this);
  final = 0;

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        Real center = (*this)(x,y,z);
        if (center > 0) continue;

        bool flip = false;
        for (int i = -1; i <= 1; i++)
          for (int j = -1; j <= 1; j++)
            for (int k = -1; k <= 1; k++)
            {
              if (i == 0 && j == 0 && k == 0) continue;
              int xNew = x + i;
              int yNew = y + j;
              int zNew = z + k;

              if (xNew < 0 || xNew > _xRes - 1) continue;
              if (yNew < 0 || yNew > _yRes - 1) continue;
              if (zNew < 0 || zNew > _zRes - 1) continue;

              if ((*this)(xNew, yNew, zNew) * center < 0)
                flip = true;
            }
         final(x,y,z) = flip ? 1 : 0;
      }

  return final;
}

//////////////////////////////////////////////////////////////////////
// distance from a point to the surface of the bounding box
//////////////////////////////////////////////////////////////////////
Real FIELD_3D::distanceToBoundingBox(const VEC3F& v)
{
  Real xMinus = _center[0] - _lengths[0] * 0.5;
  Real xPlus = _center[0] + _lengths[0] * 0.5;
  Real yMinus = _center[1] - _lengths[1] * 0.5;
  Real yPlus = _center[1] + _lengths[1] * 0.5;
  Real zMinus = _center[2] - _lengths[2] * 0.5;
  Real zPlus = _center[2] + _lengths[2] * 0.5;

  Real diff[] = {0.0f, 0.0f, 0.0f};

  // handle edges and points
  if (v[2] > zPlus) 
    diff[2] = v[2] - zPlus;
  else if (v[2] < zMinus)
    diff[2] = zMinus - v[2];
  
  if (v[1] > yPlus)
    diff[1] = v[1] - yPlus;
  else if (v[1] < yMinus)
  diff[1] = yMinus - v[1];
  
  if (v[0] > xPlus)
    diff[0] = v[0] - xPlus;
  else if (v[0] < xMinus)
    diff[0] = xMinus - v[0];
      
  return sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
}

//////////////////////////////////////////////////////////////////////
// clamp infinities to mean of the field
//////////////////////////////////////////////////////////////////////
void FIELD_3D::clampInfsToMean()
{
  Real sum = 0;
  int nonInfsSeen = 0;

  for (int x = 0; x < _totalCells; x++)
#ifdef QUAD_PRECISION
    if (!std::isinfq(_data[x]))
#else
    if (!std::isinf(_data[x]))
#endif
    {
      sum += _data[x];
      nonInfsSeen++;
    }

  Real mean = sum / nonInfsSeen;
  for (int x = 0; x < _totalCells; x++)
#ifdef QUAD_PRECISION
    if (std::isinfq(_data[x]))
#else
    if (std::isinf(_data[x]))
#endif
      _data[x] = mean;
}

//////////////////////////////////////////////////////////////////////
// clamp infinities to mean of the field
//////////////////////////////////////////////////////////////////////
void FIELD_3D::clampNansToMean()
{
  Real sum = 0;
  int nonInfsSeen = 0;

  for (int x = 0; x < _totalCells; x++)
#ifdef QUAD_PRECISION
    if (!std::isnanq(_data[x]))
#else
    if (!std::isnan(_data[x]))
#endif
    {
      sum += _data[x];
      nonInfsSeen++;
    }

  Real mean = sum / nonInfsSeen;
  for (int x = 0; x < _totalCells; x++)
#ifdef QUAD_PRECISION
    if (std::isnanq(_data[x]))
#else
    if (std::isnan(_data[x]))
#endif
      _data[x] = mean;
}

//////////////////////////////////////////////////////////////////////
// invert a distance field so that there is 1 at the border, and 
// tapers off afterwards
//////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::invertedDistance() const
{
  TIMER functionTimer(__FUNCTION__);
  FIELD_3D final = this->outerShell() + this->innerShell();

  Real eps = 1e-5;

  bool stillZeros = true;
  while (stillZeros)
  {
    stillZeros = false;
    int remainingZeros = 0;
    FIELD_3D current = final;
    for (int z = 0; z < _zRes; z++)
      for (int y = 0; y < _yRes; y++)
        for (int x = 0; x < _xRes; x++)
        {
          if (final(x,y,z) > eps) continue;
          stillZeros = true;

          Real neighbors[6];
          neighbors[0] = (x != 0)         ? final(x-1,y,z) : 0;
          neighbors[1] = (x != _xRes - 1) ? final(x+1,y,z) : 0;
          neighbors[2] = (y != 0)         ? final(x,y-1,z) : 0;
          neighbors[3] = (y != _yRes - 1) ? final(x,y+1,z) : 0;
          neighbors[4] = (z != 0)         ? final(x,y,z-1) : 0;
          neighbors[5] = (z != _zRes - 1) ? final(x,y,z+1) : 0;

          Real biggestNonzero = _totalCells + 1;

          for (int i = 0; i < 6; i++)
            if (neighbors[i] > 1e-5 && neighbors[i] < biggestNonzero)
              biggestNonzero = neighbors[i];
         
          if (biggestNonzero < _totalCells)
            current(x,y,z) = biggestNonzero + 1;
          else
            remainingZeros++;
        }
    final = current;
  }

  // take the inverse
  final = final.inverse();

  // copy over the signs
  for (int x = 0; x < final.totalCells(); x++)
    if ((*this)[x] < 0)
      final[x] *= -1;
    // penalize the outside less
    else
      final[x] *= 0.25;

  return final;
}

//////////////////////////////////////////////////////////////////////
// scale the length of the field, update the _dx and _invDx as well
//////////////////////////////////////////////////////////////////////
void FIELD_3D::scaleLengths(const Real& scale)
{
  _lengths *= scale;
  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
  _invDz = 1.0 / _dz;
}

//////////////////////////////////////////////////////////////////////
// invert a distance field so that there is 1 at the border, and 
// tapers off afterwards, making sure inside and outside both i
// sum to 0.5
//////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::invertedDistanceNormalized() const
{
  FIELD_3D final = invertedDistance();

  // normalize the whole field
  FIELD_3D temp = final;
  temp.absoluteValue();
  final *= 1.0 / temp.sum();

  Real insideSum = 0;
  Real outsideSum = 0;

  for (int x = 0; x < final.totalCells(); x++)
  {
    if ((*this)[x] < 0)
      insideSum += final[x];
    else
      outsideSum += final[x];
  }

  cout << " inside sum: " << insideSum << endl;
  cout << " outside sum: " << outsideSum << endl;

  Real insideNormalize = 0.5 / fabs(insideSum);
  Real outsideNormalize = 0.5 / fabs(outsideSum);

  for (int x = 0; x < final.totalCells(); x++)
  {
    if ((*this)[x] < 0)
      final[x] *= insideNormalize;
    else
      final[x] *= outsideNormalize;
  }

  //FIELDVIEW3D(final);
  return final;
}

//////////////////////////////////////////////////////////////////////
// set all the borders to a constant
//////////////////////////////////////////////////////////////////////
void FIELD_3D::setBorders(const Real& value)
{
  for (int x = 0; x < _xRes; x++)
    for (int y = 0; y < _yRes; y++)
    {
      (*this)(x,y,0) = value;
      (*this)(x,y,_zRes - 1) = value;
    }
  
  for (int x = 0; x < _xRes; x++)
    for (int z = 0; z < _zRes; z++)
    {
      (*this)(x,0,z) = value;
      (*this)(x,_yRes - 1,z) = value;
    }

  for (int y = 0; y < _yRes; y++)
    for (int z = 0; z < _zRes; z++)
    {
      (*this)(0,y,z) = value;
      (*this)(_xRes - 1,y,z) = value;
    }
}

//////////////////////////////////////////////////////////////////////
// count the number of infs in the field
//////////////////////////////////////////////////////////////////////
int FIELD_3D::totalInfs() const
{
  int final = 0;
  for (int x = 0; x < _totalCells; x++)
    if (std::isinf(_data[x])) final++;

  return final;
}

//////////////////////////////////////////////////////////////////////
// count the number of NaNs in the field
//////////////////////////////////////////////////////////////////////
int FIELD_3D::totalNans() const
{
  int final = 0;
  for (int x = 0; x < _totalCells; x++)
    if (std::isnan(_data[x])) final++;

  return final;
}

void FIELD_3D::rotateX(const Real& angle)
{
  MATRIX3 R = MATRIX3::rotation(VEC3F(1,0,0), angle);
  _rotation = _rotation * QUATERNION(R);
  _rotation.normalize();
}

void FIELD_3D::rotateY(const Real& angle)
{
  MATRIX3 R = MATRIX3::rotation(VEC3F(0,1,0), angle);
  _rotation = _rotation * QUATERNION(R);
  _rotation.normalize();
}

void FIELD_3D::rotateZ(const Real& angle)
{
  MATRIX3 R = MATRIX3::rotation(VEC3F(0,0,1), angle);
  _rotation = _rotation * QUATERNION(R);
  _rotation.normalize();
}
