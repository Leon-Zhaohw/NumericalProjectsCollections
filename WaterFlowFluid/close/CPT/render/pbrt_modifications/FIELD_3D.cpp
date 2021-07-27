#include "FIELD_3D.h"
#include <zlib.h>

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D::FIELD_3D(const int& xRes, const int& yRes, const int& zRes,
    const VEC3F& center, const VEC3F& lengths) :
  _xRes(xRes), _yRes(yRes), _zRes(zRes), _center(center), _lengths(lengths)
{
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  _data = new Real[_totalCells];

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;

  _outside = maxRes() * maxRes();
}

FIELD_3D::FIELD_3D(const double* data, const int& xRes, const int& yRes, const int& zRes,
    const VEC3F& center, const VEC3F& lengths) :
  _xRes(xRes), _yRes(yRes), _zRes(zRes), _center(center), _lengths(lengths)
{
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  _data = new Real[_totalCells];

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;

  _outside = maxRes() * maxRes();

  for (int x = 0; x < _totalCells; x++)
    _data[x] = data[x];
}

FIELD_3D::FIELD_3D(const FIELD_3D& m) :
  _xRes(m.xRes()), _yRes(m.yRes()), _zRes(m.zRes()),
  _center(m.center()), _lengths(m.lengths())
{
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  _data = new Real[_totalCells];

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;

  for (int x = 0; x < _totalCells; x++)
    _data[x] = m[x];

  _outside = maxRes() * maxRes();
}

FIELD_3D::FIELD_3D() :
  _xRes(-1), _yRes(-1), _zRes(-1), _totalCells(-1), _data(NULL)
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

  /*
  // write dimensions
  fwrite((void*)&_xRes, sizeof(int), 1, file);
  fwrite((void*)&_yRes, sizeof(int), 1, file);
  fwrite((void*)&_zRes, sizeof(int), 1, file);
  _center.write(file);
  _lengths.write(file);

  // always write out as a double
  if (sizeof(Real) != sizeof(double))
  {
    double* dataDouble = new double[_totalCells];
    for (int x = 0; x < _totalCells; x++)
      dataDouble[x] = _data[x];

    fwrite((void*)dataDouble, sizeof(double), _totalCells, file);
    delete[] dataDouble;
    fclose(file);
  }
  else
    fwrite((void*)_data, sizeof(Real), _totalCells, file);
  fclose(file);
  */

  cout << " done. " << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_3D::writeGz(string filename) const
{
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

  /*
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
    */
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

  /*
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
  if (_data) delete[] _data;
  _data = new Real[_totalCells];

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
  */

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

  /*
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
  */

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
    return;
  }

  if (_data) delete[] _data;

  _xRes = xRes;
  _yRes = yRes;
  _zRes = zRes;
  _slabSize = _xRes * _yRes;
  _totalCells = _xRes * _yRes * _zRes;
  _outside = maxRes() * maxRes();

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;

  _data = new Real[_totalCells];
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
FIELD_3D& FIELD_3D::operator=(const FIELD_3D& A)
{
  resizeAndWipe(A.xRes(), A.yRes(), A.zRes(), A.center(), A.lengths());

  for (int x = 0; x < _totalCells; x++)
    _data[x] = A[x];

  return *this;
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
const Real FIELD_3D::operator()(const VEC3F& position) const
{
  VEC3F positionCopy = position;

  // get the lower corner position
  VEC3F corner = _center - (Real)0.5 * _lengths;
  VEC3F dxs(_dx, _dy, _dz);
  corner += (Real)0.5 * dxs;

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= 1.0 / _dx;
  positionCopy[1] *= 1.0 / _dy;
  positionCopy[2] *= 1.0 / _dz;

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

  positionCopy[0] *= 1.0 / _dx;
  positionCopy[1] *= 1.0 / _dy;
  positionCopy[2] *= 1.0 / _dz;

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
// reset dimensions
///////////////////////////////////////////////////////////////////////
void FIELD_3D::setLengths(const VEC3F& lengths)
{
  _lengths = lengths;

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
}

///////////////////////////////////////////////////////////////////////
// see if the field has any data in it yet
///////////////////////////////////////////////////////////////////////
const bool FIELD_3D::initialized() const
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

  assert(_xRes > 0);
  assert(_yRes > 0);
  assert(_zRes > 0);

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
// write out a field to a file stream
///////////////////////////////////////////////////////////////////////
void FIELD_3D::write(FILE* file) const
{
  assert(_xRes > 0);
  assert(_yRes > 0);
  assert(_zRes > 0);

  // write dimensions
  fwrite((void*)&_xRes, sizeof(int), 1, file);
  fwrite((void*)&_yRes, sizeof(int), 1, file);
  fwrite((void*)&_zRes, sizeof(int), 1, file);
  _center.write(file);
  _lengths.write(file);

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
// do a quartic WENO interpolation
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::quarticInterp(const Real interp, const Real* points)
{
  const Real& fim1 = points[0];
  const Real& fi   = points[1];
  const Real& fip1 = points[2];
  const Real& fip2 = points[3];
  const Real& x = interp;
  const Real xHalf = interp * 0.5;

  const Real p1 = fi + ((fip1 - fim1) + (fip1 - 2.0f * fi + fim1) * x) * xHalf;
  const Real p2 = fi + ((-fip2 + 4.0f * fip1 - 3.0f * fi) + (fip2 - 2.0f * fip1 + fi) * x) * xHalf;

  const Real third = 1.0 / 3.0;
  const Real C1 = (2.0f - x) * third;
  const Real C2 = (x + 1.0f) * third;

  const Real middle = -76.0f * fip1 * fi;
  const Real fip1Sq = fip1 * fip1;
  const Real fiSq = fi * fi;

  const Real twelfth = 1.0 / 12.0;
  const Real IS1 = (26.0f * fip1 * fim1 - 52.0f * fi * fim1 + middle + 25.0f * fip1Sq + 64.0f * fiSq + 13.0f * fim1 * fim1) * twelfth + 1e-6f;
  const Real IS2 = (26.0f * fip2 * fi - 52.0f * fip2 * fip1 + middle + 25.0f * fiSq + 64.0f * fip1Sq + 13.0f * fip2 * fip2) * twelfth + 1e-6f;

  const Real alpha1 = C1 / (IS1 * IS1);
  const Real alpha2 = C2 / (IS2 * IS2);

  return (alpha1 * p1 + alpha2 * p2) / (alpha1 + alpha2);
}
