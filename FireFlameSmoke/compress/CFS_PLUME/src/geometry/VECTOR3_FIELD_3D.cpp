#include "VECTOR3_FIELD_3D.h"
//#include <omp.h>
#include "WAVELET_NOISE.h"

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/gl.h> // OpenGL itself.
//#include <GL/glu.h> // GLU support library.
#include <GL/glut.h> // GLUT support library.
#endif


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D::VECTOR3_FIELD_3D(const int& xRes, const int& yRes, const int& zRes,
    const VEC3F& center, const VEC3F& lengths) :
  _xRes(xRes), _yRes(yRes), _zRes(zRes), _center(center), _lengths(lengths), _initialized(true)
{
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  _data = new VEC3F[_totalCells];

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
}

VECTOR3_FIELD_3D::VECTOR3_FIELD_3D(double* data, const int& xRes, const int& yRes, const int& zRes,
    const VEC3F& center, const VEC3F& lengths) :
  _xRes(xRes), _yRes(yRes), _zRes(zRes), _center(center), _lengths(lengths), _initialized(true)
{
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  _data = new VEC3F[_totalCells];

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;

  for (int x = 0; x < _totalCells; x++)
  {
    _data[x][0] = data[3 * x];
    _data[x][1] = data[3 * x + 1];
    _data[x][2] = data[3 * x + 2];
  }
}

VECTOR3_FIELD_3D::VECTOR3_FIELD_3D(const VECTOR& data, const int& xRes, const int& yRes, const int& zRes,
    const VEC3F& center, const VEC3F& lengths) :
  _xRes(xRes), _yRes(yRes), _zRes(zRes), _center(center), _lengths(lengths), _initialized(true)
{
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  _data = new VEC3F[_totalCells];

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;

  for (int x = 0; x < _totalCells; x++)
  {
    _data[x][0] = data[3 * x];
    _data[x][1] = data[3 * x + 1];
    _data[x][2] = data[3 * x + 2];
  }
}

VECTOR3_FIELD_3D::VECTOR3_FIELD_3D(const VectorXd& data, const int& xRes, const int& yRes, const int& zRes,
    const VEC3F& center, const VEC3F& lengths) :
  _xRes(xRes), _yRes(yRes), _zRes(zRes), _center(center), _lengths(lengths), _initialized(true)
{
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  _data = new VEC3F[_totalCells];

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;

  for (int x = 0; x < _totalCells; x++)
  {
    _data[x][0] = data[3 * x];
    _data[x][1] = data[3 * x + 1];
    _data[x][2] = data[3 * x + 2];
  }
}

VECTOR3_FIELD_3D::VECTOR3_FIELD_3D(Real* xData, Real* yData, Real* zData, const int& xRes, const int& yRes, const int& zRes,
    const VEC3F& center, const VEC3F& lengths) :
  _xRes(xRes), _yRes(yRes), _zRes(zRes), _center(center), _lengths(lengths), _initialized(true)
{
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  _data = new VEC3F[_totalCells];

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;

  for (int x = 0; x < _totalCells; x++)
  {
    _data[x][0] = xData[x];
    _data[x][1] = yData[x];
    _data[x][2] = zData[x];
  }
}

VECTOR3_FIELD_3D::VECTOR3_FIELD_3D(const VECTOR3_FIELD_3D& m) :
  _xRes(m.xRes()), _yRes(m.yRes()), _zRes(m.zRes()),
  _center(m.center()), _lengths(m.lengths()), _initialized(true)
{
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  _data = new VEC3F[_totalCells];

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;

  for (int x = 0; x < _totalCells; x++)
    _data[x] = m[x];
}

VECTOR3_FIELD_3D::VECTOR3_FIELD_3D(const FIELD_3D& m) :
  _xRes(m.xRes()), _yRes(m.yRes()), _zRes(m.zRes()),
  _center(m.center()), _lengths(m.lengths()), _initialized(true)
{
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  _data = new VEC3F[_totalCells];

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;

  for (int x = 0; x < _totalCells; x++)
    _data[x] = m[x];
}

VECTOR3_FIELD_3D::VECTOR3_FIELD_3D() :
  //_xRes(-1), _yRes(-1), _zRes(-1), _totalCells(-1), _data(NULL), _initialized(false)
  _xRes(0), _yRes(0), _zRes(0), _totalCells(0), _data(NULL), _initialized(false)
{
}

///////////////////////////////////////////////////////////////////////
// strip the dimensions from m, but populate using data
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D::VECTOR3_FIELD_3D(const VECTOR& data, const VECTOR3_FIELD_3D& m) :
  _xRes(m.xRes()), _yRes(m.yRes()), _zRes(m.zRes()),
  _center(m.center()), _lengths(m.lengths()), _initialized(true)
{
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  _data = new VEC3F[_totalCells];

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;

  assert(data.size() == _totalCells * 3);

  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        for (int i = 0; i < 3; i++, index++)
          (*this)(x,y,z)[i] = data[index];
}

///////////////////////////////////////////////////////////////////////
// strip the dimensions from m, but populate using data
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D::VECTOR3_FIELD_3D(const VectorXd& data, const VECTOR3_FIELD_3D& m) :
  _xRes(m.xRes()), _yRes(m.yRes()), _zRes(m.zRes()),
  _center(m.center()), _lengths(m.lengths()), _initialized(true)
{
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  _data = new VEC3F[_totalCells];

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;

  assert(data.size() == _totalCells * 3);

  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        for (int i = 0; i < 3; i++, index++)
          (*this)(x,y,z)[i] = data[index];
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D::~VECTOR3_FIELD_3D()
{
  delete[] _data;
}
  
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::clear()
{
  TIMER functionTimer(__FUNCTION__);
  for (int x = 0; x < _totalCells; x++)
    _data[x] = 0.0;
}

///////////////////////////////////////////////////////////////////////
// reset the lengths to something else, and recompute all the 
// dimesions as well
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::setLengths(const VEC3F& lengths)
{
  _lengths = lengths;

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
}

///////////////////////////////////////////////////////////////////////
// create a field of the grid positions of the passed in grid
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::cellCenters(const FIELD_3D& input)
{
  int xRes = input.xRes();
  int yRes = input.yRes();
  int zRes = input.zRes();
  const VEC3F& center = input.center();
  const VEC3F& lengths = input.lengths();
  VECTOR3_FIELD_3D final(xRes, yRes, zRes, center, lengths);

  int index = 0;
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++, index++)
        final[index] = input.cellCenter(x,y,z);

  return final;
}

///////////////////////////////////////////////////////////////////////
// take the gradient of a scalar field
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::gradient(const FIELD_3D& input)
{
  int xRes = input.xRes();
  int yRes = input.yRes();
  int zRes = input.zRes();
  int slabSize = xRes * yRes;
  int totalCells = xRes * yRes * zRes;
  const VEC3F& center = input.center();
  const VEC3F& lengths = input.lengths();
  VECTOR3_FIELD_3D final(xRes, yRes, zRes, center, lengths);
 
  Real dxHalfInv = 0.5 / input.dx();
  Real dyHalfInv = 0.5 / input.dy();
  Real dzHalfInv = 0.5 / input.dz();

  // do the x middles
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 1; x < xRes - 1; x++)
      {
        int index = x + y * xRes + z * slabSize;
        final[index][0] = (input[index + 1] - input[index - 1]) * dxHalfInv;
      }

  // do the y middles
  for (int z = 0; z < zRes; z++)
    for (int y = 1; y < yRes - 1; y++)
      for (int x = 0; x < xRes; x++)
      {
        int index = x + y * xRes + z * slabSize;
        final[index][1] = (input[index + xRes] - input[index - xRes]) * dyHalfInv;
      }

  // do the z middles
  for (int z = 1; z < zRes - 1; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        int index = x + y * xRes + z * slabSize;
        final[index][2] = (input[index + slabSize] - input[index - slabSize]) * dzHalfInv;
      }

  // reset dx's to a single cell
  Real dxInv = 1.0 / input.dx();
  Real dyInv = 1.0 / input.dy();
  Real dzInv = 1.0 / input.dz();

  int index;
  for (int y = 0; y < yRes; y++)
    for (int x = 0; x < xRes; x++)
    {
      // front slab
      index = x + y * xRes;
      final[index][2] = (input[index + slabSize] - input[index]) * dzInv;

      // back slab
      index += totalCells - slabSize;
      final[index][2] = (input[index] - input[index - slabSize]) * dzInv;
    }

  for (int z = 0; z < zRes; z++)
    for (int x = 0; x < xRes; x++)
    {
      // bottom slab
      index = x + z * slabSize;
      final[index][1] = (input[index + xRes] - input[index]) * dyInv;

      // top slab
      index += slabSize - xRes;
      final[index][1] = (input[index] - input[index - xRes]) * dyInv;
    }

  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
    {
      // left slab
      index = y * xRes + z * slabSize;
      final[index][0] = (input[index + 1] - input[index]) * dxInv;

      // right slab
      index += xRes - 1;
      final[index][0] = (input[index] - input[index - 1]) * dxInv;
    }

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D VECTOR3_FIELD_3D::scalarField(int component) const
{
  assert(component >= 0 && component < 3);

  FIELD_3D final(_xRes, _yRes, _zRes, _center, _lengths);

  for (int x = 0; x < _totalCells; x++)
    final[x] = _data[x][component];

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D VECTOR3_FIELD_3D::magnitudeField() const
{
  FIELD_3D final(_xRes, _yRes, _zRes, _center, _lengths);

  for (int x = 0; x < _totalCells; x++)
    final[x] = norm(_data[x]);

  return final;
}

///////////////////////////////////////////////////////////////////////
// take the field dot product
///////////////////////////////////////////////////////////////////////
FIELD_3D operator*(const VECTOR3_FIELD_3D&u, const VECTOR3_FIELD_3D& v)
{
  assert(u.xRes() == v.xRes());
  assert(u.yRes() == v.yRes());
  assert(u.zRes() == v.zRes());

  FIELD_3D final(u.xRes(), u.yRes(), u.zRes(), u.center(), u.lengths());

  int totalCells = u.totalCells();
  for (int x = 0; x < totalCells; x++)
    final[x] = u[x] * v[x];

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D operator*(const VECTOR3_FIELD_3D& v, const Real& a)
{
  return a * v;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D operator*(const Real& a, const VECTOR3_FIELD_3D& v)
{
  VECTOR3_FIELD_3D final(v);
  
  int totalCells = v.totalCells();
  for (int x = 0; x < totalCells; x++)
    final[x] = a * v[x];

  return final;
}

/*
///////////////////////////////////////////////////////////////////////
// take the field dot product
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D operator*(const VECTOR3_FIELD_3D& u, const FIELD_3D& v)
{
  assert(u.xRes() == v.xRes());
  assert(u.yRes() == v.yRes());
  assert(u.zRes() == v.zRes());

  VECTOR3_FIELD_3D final(v);

  int totalCells = u.totalCells();
  for (int x = 0; x < totalCells; x++)
    final[x] = u[x] * v[x];

  return final;
}
*/

///////////////////////////////////////////////////////////////////////
// take the field dot product
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D operator*(const FIELD_3D& u, const VECTOR3_FIELD_3D& v)
{
  assert(u.xRes() == v.xRes());
  assert(u.yRes() == v.yRes());
  assert(u.zRes() == v.zRes());

  VECTOR3_FIELD_3D final(v);

  int totalCells = u.totalCells();
  for (int x = 0; x < totalCells; x++)
    final[x] = u[x] * v[x];

  return final;
}

///////////////////////////////////////////////////////////////////////
// sum two vector fields
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D operator+(const VECTOR3_FIELD_3D&u, const VECTOR3_FIELD_3D& v)
{
  assert(u.xRes() == v.xRes());
  assert(u.yRes() == v.yRes());
  assert(u.zRes() == v.zRes());

  VECTOR3_FIELD_3D final(u.xRes(), u.yRes(), u.zRes(), u.center(), u.lengths());
  int totalCells = u.totalCells();
  for (int x = 0; x < totalCells; x++)
    final[x] = u[x] + v[x];

  return final;
}

///////////////////////////////////////////////////////////////////////
// diff two vector fields
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D operator-(const VECTOR3_FIELD_3D&u, const VECTOR3_FIELD_3D& v)
{
  assert(u.xRes() == v.xRes());
  assert(u.yRes() == v.yRes());
  assert(u.zRes() == v.zRes());

  VECTOR3_FIELD_3D final(u.xRes(), u.yRes(), u.zRes(), u.center(), u.lengths());
  int totalCells = u.totalCells();
  for (int x = 0; x < totalCells; x++)
    final[x] = u[x] - v[x];

  return final;
}

///////////////////////////////////////////////////////////////////////
// return a grid of values at the given spatial positions
///////////////////////////////////////////////////////////////////////
FIELD_3D VECTOR3_FIELD_3D::compose(const FIELD_3D& values, const VECTOR3_FIELD_3D& positions)
{
  assert(positions.xRes() == values.xRes());
  assert(positions.yRes() == values.yRes());
  assert(positions.zRes() == values.zRes());

  // intialize to the same dims as the input
  FIELD_3D final(values);

  int i = 0;
  for (int z = 0; z < values.zRes(); z++)
    for (int y = 0; y < values.yRes(); y++)
      for (int x = 0; x < values.xRes(); x++, i++)
        final[i] = values(positions[i]);

  return final;
}

///////////////////////////////////////////////////////////////////////
// return a grid of values at the given spatial positions
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::compose(const VECTOR3_FIELD_3D& values, const VECTOR3_FIELD_3D& positions)
{
  assert(positions.xRes() == values.xRes());
  assert(positions.yRes() == values.yRes());
  assert(positions.zRes() == values.zRes());

  // intialize to the same dims as the input
  VECTOR3_FIELD_3D final(values);

  int i = 0;
  for (int z = 0; z < values.zRes(); z++)
    for (int y = 0; y < values.yRes(); y++)
      for (int x = 0; x < values.xRes(); x++, i++)
        final[i] = values(positions[i]);

  return final;
}

///////////////////////////////////////////////////////////////////////
// lookup value at some real-valued spatial position
///////////////////////////////////////////////////////////////////////
const VEC3F VECTOR3_FIELD_3D::operator()(const VEC3F& position) const
{
  /*
  VEC3F positionCopy = position;

  // get the lower corner position
  VEC3F corner = _center - (Real)0.5 * _lengths;
  VEC3F dxs(_dx, _dy, _dz);
  corner += (Real)0.5 * dxs;

  // recenter position
  positionCopy -= corner;
  */
  VEC3F positionCopy = position - _center + (Real)0.5 * _lengths - (Real)0.5 * dxs();

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
  const Real s1 = positionCopy[0]- x0;
  const Real s0 = 1.0f - s1;
  const Real t1 = positionCopy[1]- y0;
  const Real t0 = 1.0f - t1;
  const Real u1 = positionCopy[2]- z0;
  const Real u0 = 1.0f - u1;

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
// lookup value at some real-valued spatial position
///////////////////////////////////////////////////////////////////////
VEC3F VECTOR3_FIELD_3D::debugPositionOperator(const VEC3F& position) const
{
  /*
  VEC3F positionCopy = position;

  // get the lower corner position
  VEC3F corner = _center - (Real)0.5 * _lengths;
  VEC3F dxs(_dx, _dy, _dz);
  corner += (Real)0.5 * dxs;

  // recenter position
  positionCopy -= corner;
  */
  VEC3F positionCopy = position - _center + (Real)0.5 * _lengths - (Real)0.5 * dxs();

  positionCopy[0] *= 1.0 / _dx;
  positionCopy[1] *= 1.0 / _dy;
  positionCopy[2] *= 1.0 / _dz;

  int x0 = (int)positionCopy[0];
  int x1    = x0 + 1;
  int y0 = (int)positionCopy[1];
  int y1    = y0 + 1;
  int z0 = (int)positionCopy[2];
  int z1    = z0 + 1;

  cout << " recomputed fine: " << x0 << " " << y0 << " " << z0 << endl;
  cout << " recomputed dx: " << _lengths[0] / _xRes << endl;
  cout << " length: " << _lengths << endl;
  cout << " center: " << _center << endl;

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
  const Real s1 = positionCopy[0]- x0;
  const Real s0 = 1.0f - s1;
  const Real t1 = positionCopy[1]- y0;
  const Real t0 = 1.0f - t1;
  const Real u1 = positionCopy[2]- z0;
  const Real u0 = 1.0f - u1;

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
// norms
///////////////////////////////////////////////////////////////////////
Real VECTOR3_FIELD_3D::sumMagnitudes()
{
  Real final = 0;
  for (int i = 0; i < _totalCells; i++)
    final += norm(_data[i]);

  return final;
}

///////////////////////////////////////////////////////////////////////
// norms
///////////////////////////////////////////////////////////////////////
Real VECTOR3_FIELD_3D::maxMagnitude()
{
  Real final = norm(_data[0]);
  for (int i = 0; i < _totalCells; i++)
    if (norm(_data[i]) > final)
      final = norm(_data[i]);

  return final;
}

///////////////////////////////////////////////////////////////////////
// norms
///////////////////////////////////////////////////////////////////////
Real VECTOR3_FIELD_3D::maxMagnitude(int& index)
{
  Real final = norm(_data[0]);
  index = 0;
  for (int i = 0; i < _totalCells; i++)
    if (norm(_data[i]) > final)
    {
      final = norm(_data[i]);
      index = i;
    }

  return final;
}

///////////////////////////////////////////////////////////////////////
// check if any entry is a nan
///////////////////////////////////////////////////////////////////////
bool VECTOR3_FIELD_3D::isNan()
{
  for (int x = 0; x < _totalCells; x++)
    for (int y = 0; y < 3; y++)
      if (isnan(_data[x][y]))
        return true;

  return false;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D& VECTOR3_FIELD_3D::operator-=(const VECTOR3_FIELD_3D& input)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] -= input[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D& VECTOR3_FIELD_3D::operator+=(const VECTOR3_FIELD_3D& input)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] += input[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D& VECTOR3_FIELD_3D::operator*=(const FIELD_3D& input)
{
  assert(_xRes == input.xRes());
  assert(_yRes == input.yRes());
  assert(_zRes == input.zRes());

  for (int x = 0; x < _totalCells; x++)
  {
    _data[x][0] *= input[x];
    _data[x][1] *= input[x];
    _data[x][2] *= input[x];
  }

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D& VECTOR3_FIELD_3D::operator*=(const VEC3F& alpha)
{
  for (int x = 0; x < _totalCells; x++)
  {
    _data[x][0] *= alpha[0];
    _data[x][1] *= alpha[1];
    _data[x][2] *= alpha[2];
  }

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D& VECTOR3_FIELD_3D::operator+=(const Real& alpha)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] += alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D& VECTOR3_FIELD_3D::operator*=(const Real& alpha)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] *= alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D& VECTOR3_FIELD_3D::operator=(const Real& value)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] = value;

  return *this;
}

///////////////////////////////////////////////////////////////////////
// extend some vector quantity off of a front, given a signed distance function
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::fastExtension(const FIELD_3D& signedDistance)
{
  cout << " Extending vectors ... "; flush(cout);
  // assume that initializeExtensionScalars has been called elsewhere

  // insert the front for marching onto the heap 
  MIN_HEAP minHeap;

  // clear the retired nodes
  _retired.clear();

  FIELD_3D distanceCopy(signedDistance);
  insertFront(true, distanceCopy, minHeap);

  // march forward
  extendOneway(true, distanceCopy, minHeap);

  insertFront(false, distanceCopy, minHeap);

  extendOneway(false, distanceCopy, minHeap);

  cout << " done. " << endl;
}

///////////////////////////////////////////////////////////////////////
// insert the front in preparation for reinitialization or extension
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::insertFront(const bool forward, FIELD_3D& distance, MIN_HEAP& minHeap)
{
  Real outside = distance.outside();

  // insert the ones for the forward marching
  minHeap.clear();
  for (int i = 0; i < _totalCells; i++)
  {
    bool compare = (forward) ? distance[i] >= 0.0f && distance[i] < outside
                             : distance[i] <= 0.0f;

    Real sum = 0.0f;
    if (compare)
    {
      int zIndex = i / _slabSize;
      int yIndex = (i - zIndex * _slabSize) / _xRes;
      int xIndex = (i - zIndex * _slabSize - yIndex * _xRes);
      Real center = distance[i];

      // x plus and minus
      Real xPlus  = (xIndex < _xRes - 1) ? distance[i + 1] : outside;
      Real xMinus = (xIndex > 0)         ? distance[i - 1] : outside;

      // y plus and minus
      Real yPlus  = (yIndex < _yRes - 1) ? distance[i + _xRes] : outside;
      Real yMinus = (yIndex > 0)         ? distance[i - _xRes] : outside;

      // z plus and minus
      Real zPlus  = (zIndex < _zRes - 1) ? distance[i + _slabSize] : outside;
      Real zMinus = (zIndex > 0)         ? distance[i - _slabSize] : outside;

      Real interpolate;
      compare = (forward) ? true : yPlus < outside;
      if (yPlus * center <= 0.0f && compare)
      {
        interpolate = center / (center - yPlus);
        sum += 1.0f / (interpolate * interpolate);
      }
      compare = (forward) ? true : yMinus < outside;
      if (yMinus * center <= 0.0f && compare)
      {
        interpolate = center / (center - yMinus);
        sum += 1.0f / (interpolate * interpolate);
      }
      compare = (forward) ? true : xMinus < outside;
      if (xMinus * center <= 0.0f && compare)
      {
        interpolate = center / (center - xMinus);
        sum += 1.0f / (interpolate * interpolate);
      }
      compare = (forward) ? true : xPlus < outside;
      if (xPlus * center <= 0.0f && compare)
      {
        interpolate = center / (center - xPlus);
        sum += 1.0f / (interpolate * interpolate);
      }
      compare = (forward) ? true : zMinus < outside;
      if (zMinus * center <= 0.0f && compare)
      {
        interpolate = center / (center - zMinus);
        sum += 1.0f / (interpolate * interpolate);
      }
      compare = (forward) ? true : zPlus < outside;
      if (zPlus * center <= 0.0f && compare)
      {
        interpolate = center / (center - zPlus);
        sum += 1.0f / (interpolate * interpolate);
      }

      Real sign = (forward) ? 1.0f : -1.0f;

      if (sum > 0.0f)
      {
        Real finalDistance = sign / sqrtf(sum);
       
        HEAP_ENTRY entry;
        entry.distance = finalDistance;
        entry.index = i;
        distance[i] = finalDistance;
        minHeap.insert(entry);
      }
      else
        distance[i] = sign * outside;
    }
    if (!forward && distance[i] >= outside)
      distance[i] = -outside;
  }
}

//////////////////////////////////////////////////////////////////////
// do extension in one direction
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::extendOneway(bool forward, FIELD_3D& distance, MIN_HEAP& minHeap)
{
  int pops = 0;
  Real outside = distance.outside();
  
  // do the forward marching
  while (!minHeap.empty())
  {
    // pop off the top and retire it
    HEAP_ENTRY popped = minHeap.popMin();
    _retired[popped.index] = true;

    pops++;

    // popped was just retired, calculate a new distance value for all its neighbors
    int candidate = popped.index;
    int zIndex = candidate / _slabSize;
    int yIndex = (candidate - zIndex * _slabSize) / _xRes;
    int xIndex = (candidate - zIndex * _slabSize - yIndex * _xRes);
  
    for (int y = 0; y < 6; y++)
    {
      candidate = popped.index;

      if (y == 0)
        if (xIndex < _xRes - 1) candidate++;
        else continue;

      if (y == 1) 
        if (xIndex > 0) candidate--;
        else continue;

      if (y == 2)
        if (yIndex < _yRes - 1) candidate += _xRes;
        else continue;

      if (y == 3) 
        if (yIndex > 0) candidate -= _xRes;
        else continue;

      if (y == 4) 
        if (zIndex < _zRes - 1) candidate += _slabSize;
        else continue;

      if (y == 5)
        if (zIndex > 0) candidate -= _slabSize;
        else continue;

      // account for marching direction
      bool distanceTest = distance[candidate] >= 0.0 ? true : false;
      if (!forward) distanceTest = !distanceTest;

      if (distanceTest && (_retired.find(candidate) == _retired.end()))
      {
        // get the distance values
        Real distances[3];
        VEC3F extensions[3];
        Real neighbors[6];
        VEC3F extensionNeighbors[6];

        int zIndex = candidate / _slabSize;
        int yIndex = (candidate - zIndex * _slabSize) / _xRes;
        int xIndex = (candidate - zIndex * _slabSize - yIndex * _xRes);

        // x plus and minus
        extensionNeighbors[0] = (xIndex < _xRes - 1) ? _data[candidate + 1] : _data[candidate];
        extensionNeighbors[1] = (xIndex > 0)         ? _data[candidate - 1] : _data[candidate];
        neighbors[0] = (xIndex < _xRes - 1) ? distance[candidate + 1] : distance[candidate];
        neighbors[1] = (xIndex > 0)         ? distance[candidate - 1] : distance[candidate];

        // y plus and minus
        extensionNeighbors[2] = (yIndex < _yRes - 1) ? _data[candidate + _xRes] : _data[candidate];
        extensionNeighbors[3] = (yIndex > 0)         ? _data[candidate - _xRes] : _data[candidate];
        neighbors[2] = (yIndex < _yRes - 1) ? distance[candidate + _xRes] : distance[candidate];
        neighbors[3] = (yIndex > 0)         ? distance[candidate - _xRes] : distance[candidate];

        // z plus and minus
        extensionNeighbors[4] = (zIndex < _zRes - 1) ? _data[candidate + _slabSize] : _data[candidate];
        extensionNeighbors[5] = (zIndex > 0)         ? _data[candidate - _slabSize] : _data[candidate];
        neighbors[4] = (zIndex < _zRes - 1) ? distance[candidate + _slabSize] : distance[candidate];
        neighbors[5] = (zIndex > 0)         ? distance[candidate - _slabSize] : distance[candidate];
      
        for (int x = 0; x < 3; x++)
        {
          // store the positive direction
          distances[x] = neighbors[2 * x];
          extensions[x] = extensionNeighbors[2 * x];
          
          // clamp out the wrong side
          if (forward && distances[x] < 0.0f)  distances[x] = outside;
          if (!forward && distances[x] > 0.0f) distances[x] = -outside;
          
          // see if the negative direction is more upwind
          if (forward)
          {
            if (neighbors[2 * x + 1] < distances[x] && neighbors[2 * x + 1] > 0.0f)
            {
              distances[x] = neighbors[2 * x + 1];
              extensions[x] = extensionNeighbors[2 * x + 1];
            }
          }
          else if (neighbors[2 * x + 1] > distances[x] && neighbors[2 * x + 1] < 0.0f)
          {
            distances[x] = neighbors[2 * x + 1];
            extensions[x] = extensionNeighbors[2 * x + 1];
          }
        }

        // clamp out values from the wrong side
        for (int x = 0; x < 3; x++)
          if (forward)
            distances[x] = (distances[x] >= 0.0f) ? distances[x] : outside;
          else
            distances[x] = (distances[x] <= 0.0f) ? -distances[x] : outside;
      
        // do a 3-element bubble sort 
        Real temp;
        VEC3F tempPoint;
        if (distances[0] > distances[1])
        {
          temp = distances[0];
          distances[0] = distances[1];
          distances[1] = temp;
          tempPoint = extensions[0];
          extensions[0] = extensions[1];
          extensions[1] = tempPoint;
        }
        if (distances[1] > distances[2])
        {
          temp = distances[1];
          distances[1] = distances[2];
          distances[2] = temp;
          tempPoint = extensions[1];
          extensions[1] = extensions[2];
          extensions[2] = tempPoint;
        }
        if (distances[0] > distances[1])
        {
          temp = distances[0];
          distances[0] = distances[1];
          distances[1] = temp;
          tempPoint = extensions[0];
          extensions[0] = extensions[1];
          extensions[1] = tempPoint;
        }
        
        // set up the quadratics David's way
        Real b[3];
        Real c[3];
        b[0] = distances[0];
        c[0] = distances[0] * distances[0] - 1.0f;
        for (int x = 1; x < 3; x++)
        {
          b[x] = distances[x] + b[x-1];
          c[x] = distances[x] * distances[x] + c[x-1];
        }

        // solve for the right one
        int i = 2;
        Real discrim = b[i] * b[i] - (Real)(i + 1) * c[i];
        Real newDist = (b[i] + sqrtf(discrim)) / (Real)(i + 1);
        while ((discrim < 0.0f || newDist < distances[i]) && i != -1)
        {
          i--;
          discrim = b[i] * b[i] - (Real)(i + 1) * c[i];
          
          if (discrim > 0.0f)
            newDist = (b[i] + sqrtf(discrim)) / (Real)(i + 1);
        }
        if (i == -1)
          cout << __FILE__ << " " << __LINE__ 
               << " Couldn't solve the distance field quadratic!" << endl;
       
        // get the extension
        if (newDist < fabs(distance[candidate]))
        {
          if (i == 2)
          {
            Real interpolate = 1.0f / (3.0f * newDist - distances[2] - distances[1] - distances[0]);
            Real diffs[] = {(newDist - distances[0]), (newDist - distances[1]), (newDist - distances[2])};
            _data[candidate] = (diffs[0] * extensions[0] +
                                diffs[1] * extensions[1] +
                                diffs[2] * extensions[2]) * interpolate;
          }
          else if (i == 1)
          {
            Real interpolate = 1.0f / (2.0f * newDist - distances[1] - distances[0]);
            Real diffs[] = {(newDist - distances[0]), (newDist - distances[1])};
            _data[candidate] = (diffs[0] * extensions[0] +
                                diffs[1] * extensions[1]) * interpolate;
          }
          else
            _data[candidate] = extensions[0];
        }

        // try and insert the new distance value
        if (newDist < outside)
        {
          Real minus = forward ? 1.0f : -1.0f;

          // if it has never been on the heap
          if (fabs(distance[candidate]) >= outside)
          {
            distance[candidate] = minus * newDist;
            HEAP_ENTRY entry;
            entry.distance = distance[candidate];
            entry.index = candidate;
            minHeap.insert(entry);
          }
          else if (minus * distance[candidate] > newDist)
          {
            distance[candidate] = minus * newDist;
            minHeap.decreaseKey(candidate, distance[candidate]);
          }
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D& VECTOR3_FIELD_3D::operator=(const VECTOR3_FIELD_3D& input)
{
  if (input.xRes() != _xRes ||
      input.yRes() != _yRes ||
      input.zRes() != _zRes)
  {
    delete[] _data;

    _xRes = input.xRes();
    _yRes = input.yRes();
    _zRes = input.zRes();

    _totalCells = _xRes * _yRes * _zRes;
    _slabSize = _xRes * _yRes;
    _data = new VEC3F[_totalCells];
  }

  _center = input.center();
  _lengths = input.lengths();

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  
  for (int x = 0; x < _totalCells; x++)
    _data[x] = input[x];

  _initialized = input._initialized;

  return *this;
}


//////////////////////////////////////////////////////////////////////
// set the values in the field to uniform random doubles from 0 to 1
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::setToRandom()
{
  for (int index = 0; index < _totalCells; index++) {
   (_data[index])[0] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
    (_data[index])[1] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
    (_data[index])[2] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
  }
}
//////////////////////////////////////////////////////////////////////
// set the values in the field to the values at the closest points
//////////////////////////////////////////////////////////////////////
FIELD_3D VECTOR3_FIELD_3D::setToClosestPointValues(const FIELD_3D& input, const VECTOR3_FIELD_3D& closestPoints)
{
  TIMER functionTimer(__FUNCTION__);
  FIELD_3D final(input);

  //int before = input.quinticClamps();
  cout << " Using quinitic "; flush(cout);
  //cout << " Using cubic"; flush(cout);
  //cout << " Using linear"; flush(cout);

#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < closestPoints.zRes(); z++)
    for (int y = 0; y < closestPoints.yRes(); y++)
      for (int x = 0; x < closestPoints.xRes(); x++)
        //final(x,y,z) = input(closestPoints(x,y,z));
        //final(x,y,z) = input.cubicLookup(closestPoints(x,y,z));
        final(x,y,z) = input.quinticLookup(closestPoints(x,y,z));
        //final(x,y,z) = input.cubicNewtonLookup(closestPoints(x,y,z));

  /*
  int totalCalls = 21 * closestPoints.totalCells();
  int after = input.quinticClamps();
  int clamps = after - before;

  cout << " Quintic clamps: " << clamps << " of " << totalCalls << " " << 100.0 * clamps / (Real)totalCalls << "%" << endl;
  */

  return final;
}

//////////////////////////////////////////////////////////////////////
// set the values in the field to the values at the closest points
//////////////////////////////////////////////////////////////////////
FIELD_3D VECTOR3_FIELD_3D::setToClosestPointValuesNarrowBand(const FIELD_3D& input, const VECTOR3_FIELD_3D& closestPoints, const FIELD_3D& distance, int maxCells)
{
  FIELD_3D final(input);
  Real invDx = 1.0 / distance.dx();

  // find which cells are inside the band
  vector<int> xs;
  vector<int> ys;
  vector<int> zs;
  for (int z = 0; z < closestPoints.zRes(); z++)
    for (int y = 0; y < closestPoints.yRes(); y++)
      for (int x = 0; x < closestPoints.xRes(); x++)
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
#pragma omp for  schedule(dynamic)
  for (int i = 0; i < size; i++)
  {
    int x = xs[i];
    int y = ys[i];
    int z = zs[i];
        //final(x,y,z) = input(closestPoints(x,y,z));
        final(x,y,z) = input.cubicLookup(closestPoints(x,y,z));
        //final(x,y,z) = input.quinticLookup(closestPoints(x,y,z));
        //final(x,y,z) = input.cubicNewtonLookup(closestPoints(x,y,z));
  }

  return final;
}

//////////////////////////////////////////////////////////////////////
// file IO
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::write(const string filename) const
{
  FILE* file;
  file = fopen(filename.c_str(), "wb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " VECTOR3_FIELD_3D write failed! " << endl;
    cout << " Could not open file " << filename.c_str() << endl;
    exit(0);
  }

  /*
  // write dimensions
  fwrite((void*)&_xRes, sizeof(int), 1, file);
  fwrite((void*)&_yRes, sizeof(int), 1, file);
  fwrite((void*)&_zRes, sizeof(int), 1, file);
  _center.write(file);
  _lengths.write(file);

  double* dataDouble = new double[3 * _totalCells];
  for (int x = 0; x < _totalCells; x++)
  {
    dataDouble[3 * x] = _data[x][0];
    dataDouble[3 * x + 1] = _data[x][1];
    dataDouble[3 * x + 2] = _data[x][2];
  }

  fwrite((void*)dataDouble, sizeof(double), 3 * _totalCells, file);
  delete[] dataDouble;
  */
  write(file);
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::write(FILE* file) const
{
  // write dimensions
  fwrite((void*)&_xRes, sizeof(int), 1, file);
  fwrite((void*)&_yRes, sizeof(int), 1, file);
  fwrite((void*)&_zRes, sizeof(int), 1, file);
  _center.write(file);
  _lengths.write(file);

  double* dataDouble = new double[3 * _totalCells];
  for (int x = 0; x < _totalCells; x++)
  {
    dataDouble[3 * x] = _data[x][0];
    dataDouble[3 * x + 1] = _data[x][1];
    dataDouble[3 * x + 2] = _data[x][2];
  }

  fwrite((void*)dataDouble, sizeof(double), 3 * _totalCells, file);
  delete[] dataDouble;
}

///////////////////////////////////////////////////////////////////////
// write out a field to a file stream
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::writeGz(gzFile& file) const
{
  // write dimensions
  gzwrite(file, (void*)&_xRes, sizeof(int));
  gzwrite(file, (void*)&_yRes, sizeof(int));
  gzwrite(file, (void*)&_zRes, sizeof(int));
  _center.writeGz(file);
  _lengths.writeGz(file);

  double* dataDouble = new double[3 * _totalCells];
  for (int x = 0; x < _totalCells; x++)
  {
    dataDouble[3 * x] = _data[x][0];
    dataDouble[3 * x + 1] = _data[x][1];
    dataDouble[3 * x + 2] = _data[x][2];
  }

  gzwrite(file, (void*)dataDouble, sizeof(double) * 3 * _totalCells);
  delete[] dataDouble;
}

///////////////////////////////////////////////////////////////////////
// read in a field from a file stream
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::readGz(gzFile& file)
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

  _slabSize = _xRes * _yRes;
  _totalCells = _xRes * _yRes * _zRes;
  if (_data) delete[] _data;
  _data = new VEC3F[_totalCells];

  double* dataDouble = new double[3 * _totalCells];

  int totalBytes = sizeof(double) * 3 * _totalCells;
  int readBytes = gzread(file, (void*)dataDouble, totalBytes);

  assert(readBytes == totalBytes);

  for (int x = 0; x < _totalCells; x++)
  {
    _data[x][0]= dataDouble[3 * x];
    _data[x][1]= dataDouble[3 * x + 1];
    _data[x][2]= dataDouble[3 * x + 2];
  }

  delete[] dataDouble;
}

//////////////////////////////////////////////////////////////////////
// file IO
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::read(const string filename)
{
  FILE* file;
  file = fopen(filename.c_str(), "rb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " VECTOR3_FIELD_3D read failed! " << endl;
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

  _slabSize = _xRes * _yRes;
  _totalCells = _xRes * _yRes * _zRes;
  if (_data) delete[] _data;
  _data = new VEC3F[_totalCells];

  double* dataDouble = new double[3 * _totalCells];
  fread((void*)dataDouble, sizeof(double), 3 * _totalCells, file);
  for (int x = 0; x < _totalCells; x++)
  {
    _data[x][0]= dataDouble[3 * x];
    _data[x][1]= dataDouble[3 * x + 1];
    _data[x][2]= dataDouble[3 * x + 2];
  }
  delete[] dataDouble;
  */
  read(file);
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::read(FILE* file)
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

  _slabSize = _xRes * _yRes;
  _totalCells = _xRes * _yRes * _zRes;
  if (_data) delete[] _data;
  _data = new VEC3F[_totalCells];

  double* dataDouble = new double[3 * _totalCells];
  fread((void*)dataDouble, sizeof(double), 3 * _totalCells, file);
  for (int x = 0; x < _totalCells; x++)
  {
    _data[x][0]= dataDouble[3 * x];
    _data[x][1]= dataDouble[3 * x + 1];
    _data[x][2]= dataDouble[3 * x + 2];
  }
  delete[] dataDouble;
}
//////////////////////////////////////////////////////////////////////
// draw to GL
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::draw(const VECTOR3_FIELD_3D& origins) const
{
  int stride = 6;

  glPointSize(1);

  /*
  glBegin(GL_LINES);
  glColor4f(0,0,0,1);
  for (int z = 0; z < _zRes; z+=stride)
    for (int y = 0; y < _yRes; y+=stride)
      for (int x = 0; x < _xRes; x+=stride)
      {
        const VEC3F& origin = origins(x,y,z);
        const VEC3F& endpoint = (*this)(x,y,z);
        glVertex3f(origin[0], origin[1], origin[2]);
        glVertex3f(endpoint[0], endpoint[1], endpoint[2]);
      }
  glEnd();
  */

  glBegin(GL_POINTS);
  glColor4f(1,0,0,1);
  for (int z = 0; z < _zRes; z+=stride)
    for (int y = 0; y < _yRes; y+=stride)
      for (int x = 0; x < _xRes; x+=stride)
      {
        const VEC3F& origin = origins(x,y,z);
        glVertex3f(origin[0], origin[1], origin[2]);
      }
  glEnd();

  glBegin(GL_POINTS);
  glColor4f(0,0,1,1);
  for (int z = 0; z < _zRes; z+=stride)
    for (int y = 0; y < _yRes; y+=stride)
      for (int x = 0; x < _xRes; x+=stride)
      {
        const VEC3F& endpoint = (*this)(x,y,z);
        glVertex3f(endpoint[0], endpoint[1], endpoint[2]);
      }
  glEnd();
}

//////////////////////////////////////////////////////////////////////
// draw to GL
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::drawZSlice(const VECTOR3_FIELD_3D& origins, const int zSlice, const Real scale, const int stride) const
{
  float radius = 0.001;
  
  glBegin(GL_LINES);
  for (int z = 0; z < _zRes; z++)
  {
    if (z != zSlice) continue;

    for (int y = 0; y < _yRes; y+=stride)
      for (int x = 0; x < _xRes; x+=stride)
      {
        const VEC3F& origin = origins(x,y,z);
        VEC3F endpoint = origin + scale * (*this)(x,y,z);

        glColor4f(1,1,1,1);
        glVertex3f(origin[0], origin[1], origin[2]);
        glColor4f(0,0,0,0);
        glVertex3f(endpoint[0], endpoint[1], endpoint[2]);
      }
  }
  glEnd();

  glColor4f(1,0,0,1);
  for (int z = 0; z < _zRes; z++)
  {
    if (z != zSlice) continue;

    for (int y = 0; y < _yRes; y+=stride)
      for (int x = 0; x < _xRes; x+=stride)
      {
        const VEC3F& origin = origins(x,y,z);
        glPushMatrix();
        glTranslatef(origin[0], origin[1], origin[2]);
        glutSolidSphere(radius, 10,10);
        glPopMatrix();
      }
  }

  /*
  glColor4f(0,0,1,1);
  for (int z = 0; z < _zRes; z++)
  {
    if (z != zSlice) continue;

    for (int y = 0; y < _yRes; y+=stride)
      for (int x = 0; x < _xRes; x+=stride)
      {
        //const VEC3F& endpoint = (*this)(x,y,z);
        const VEC3F& origin = origins(x,y,z);
        VEC3F endpoint = origin + scale * (*this)(x,y,z);
        glPushMatrix();
        glTranslatef(endpoint[0], endpoint[1], endpoint[2]);
        glutSolidSphere(radius, 10,10);
        glPopMatrix();
      }
  }
  */

  /*
  glBegin(GL_POINTS);
  glColor4f(1,0,0,1);
  for (int z = 0; z < _zRes; z++)
  {
    if (z != zSlice) continue;

    for (int y = 0; y < _yRes; y+=stride)
      for (int x = 0; x < _xRes; x+=stride)
      {
        const VEC3F& origin = origins(x,y,z);
        glVertex3f(origin[0], origin[1], origin[2]);
      }
  }
  glEnd();

  glBegin(GL_POINTS);
  glColor4f(0,0,1,1);
  for (int z = 0; z < _zRes; z++)
  {
    if (z != zSlice) continue;

    for (int y = 0; y < _yRes; y+=stride)
      for (int x = 0; x < _xRes; x+=stride)
      {
        const VEC3F& endpoint = (*this)(x,y,z);
        glVertex3f(endpoint[0], endpoint[1], endpoint[2]);
      }
  }
  glEnd();
  */
}

//////////////////////////////////////////////////////////////////////
// draw to GL, but assume the field represents closest points,
// so don't add the origin to the endpoint
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::drawClosestPointZSlice(const VECTOR3_FIELD_3D& origins, const int zSlice, const Real scale, const int stride) const
{
  float radius = 0.00125;
  
  glBegin(GL_LINES);
  for (int z = 0; z < _zRes; z++)
  {
    if (z != zSlice) continue;

    for (int y = 0; y < _yRes; y+=stride)
      for (int x = 0; x < _xRes; x+=stride)
      {
        const VEC3F& origin = origins(x,y,z);
        VEC3F endpoint = (*this)(x,y,z);

        if (norm2(endpoint) < 1e-5)
          continue;

        //glColor4f(1,1,1,1);
        glColor4f(0.25,0.25,0.25,1);
        glVertex3f(origin[0], origin[1], origin[2]);
        glVertex3f(endpoint[0], endpoint[1], 0);
      }
  }
  glEnd();

  for (int z = 0; z < _zRes; z++)
  {
    if (z != zSlice) continue;

    for (int y = 0; y < _yRes; y+=stride)
      for (int x = 0; x < _xRes; x+=stride)
      {
        const VEC3F& origin = origins(x,y,z);
        const VEC3F endpoint = (*this)(x,y,z);

        if (norm2(endpoint) < 1e-5)
          continue;

        glColor4f(1,0,0,1);
        glPushMatrix();
        glTranslatef(origin[0], origin[1], origin[2]);
        glutSolidSphere(radius, 10,10);
        glPopMatrix();
        
        glColor4f(0,0,1,1);
        glPushMatrix();
        glTranslatef(endpoint[0], endpoint[1], endpoint[2]);
        glutSolidSphere(radius, 10,10);
        glPopMatrix();

      }
  }
}

//////////////////////////////////////////////////////////////////////
// dump to a viewer
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::fieldViewer(const VECTOR3_FIELD_3D& field, const FIELD_3D& distanceField, string name)
{
  field.write("temp3d.vector.field");
  distanceField.write("temp3d.field");
  string execute("./bin/vectorFieldViewer3D temp3d.vector.field temp3d.field \"");
  execute = execute + name + string("\" &");
  cout << " Executing " << execute.c_str() << endl;
  system(execute.c_str());
}

//////////////////////////////////////////////////////////////////////
// dump to a viewer
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::gradientViewer(const VECTOR3_FIELD_3D& field, const VECTOR3_FIELD_3D& gradientField, string name)
{
  field.write("temp3d.vector.field");
  gradientField.write("temp3d.gradient.field");
  string execute("./bin/gradientFieldViewer3D temp3d.vector.field temp3d.gradient.field \"");
  execute = execute + name + string("\" &");
  cout << " Executing " << execute.c_str() << endl;
  system(execute.c_str());
}

//////////////////////////////////////////////////////////////////////
// dump to a viewer
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::overlayFieldViewer(const VECTOR3_FIELD_3D& field, const FIELD_3D& distanceField, string name)
{
  field.write("temp3d.vector.field");
  distanceField.write("temp3d.field");
  string execute("./bin/overlayVectorFieldViewer3D temp3d.vector.field temp3d.field \"");
  execute = execute + name + string("\" &");
  cout << " Executing " << execute.c_str() << endl;
  system(execute.c_str());
}

//////////////////////////////////////////////////////////////////////
// dump to a viewer
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::closestPointViewer(const VECTOR3_FIELD_3D& field, const FIELD_3D& distanceField, string name)
{
  field.write("temp3d.vector.field");
  distanceField.write("temp3d.field");
  string execute("./bin/closestPointViewer3D temp3d.vector.field temp3d.field \"");
  execute = execute + name + string("\" &");
  cout << " Executing " << execute.c_str() << endl;
  system(execute.c_str());
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::readPhysBAMGz(const FIELD_3D& distanceField, const string filename)
{
  VECTOR3_FIELD_3D unpadded;
  unpadded.readPhysBAMGz(filename);

  *this = VECTOR3_FIELD_3D(distanceField);

  // TODO: is the padding consistent?
  //
  // appears to be, at least for these boundary conditions
  for (int z = 0; z < unpadded.zRes(); z++)
    for (int y = 0; y < unpadded.yRes(); y++)
      for (int x = 0; x < unpadded.xRes(); x++)
        (*this)(x + 2, y + 2, z + 2) = unpadded(x,y,z);
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::readPhysBAMGz(const string filename)
{
  gzFile file;
  file = gzopen(filename.c_str(), "rb1");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " VECTOR3_FIELD_3D read failed! " << endl;
    cout << " Could not open file " << filename.c_str() << endl;
    exit(0);
  }
  cout << " Reading PhysBAM velocity " << filename << " ... "; flush(cout);

  // read in a RANGE<TV_INT>, i.e. index bounds
  int uBegin;
  int uEnd;
  int vBegin;
  int vEnd;
  int wBegin;
  int wEnd;
    
  gzread(file, (void*)&uBegin, sizeof(int));
  gzread(file, (void*)&uEnd, sizeof(int));
  gzread(file, (void*)&vBegin, sizeof(int));
  gzread(file, (void*)&vEnd, sizeof(int));
  gzread(file, (void*)&wBegin, sizeof(int));
  gzread(file, (void*)&wEnd, sizeof(int));

  //cout << " u: " << uBegin << " " << uEnd << endl;
  //cout << " v: " << vBegin << " " << vEnd << endl;
  //cout << " w: " << wBegin << " " << wEnd << endl;

  // read in a buffer size
  //int guess = 3 * ((uEnd + 1) * vEnd * wEnd);
  //int guess = ((uEnd + 1) * vEnd * wEnd) + ((vEnd + 1) * uEnd * wEnd) + ((wEnd + 1) * vEnd * uEnd);

  int temp = uEnd;
  uEnd = wEnd;
  wEnd = temp;

  int scanTotalX = (uEnd + 1) * vEnd * wEnd;
  int scanTotalY = uEnd * (vEnd + 1) * wEnd;
  int scanTotalZ = uEnd * vEnd * (wEnd + 1);
  int guess = scanTotalX+scanTotalY+scanTotalZ;
  //cout << " total guessed cells: " << guess << endl;
  int totalCells;
  gzread(file, (void*)&totalCells, sizeof(int));
  //cout << " total read cells: " << totalCells << endl;

  // am I computing the total cells right?
  assert(guess == totalCells);

  float* data = new float[totalCells];

  gzread(file, data, sizeof(float) * totalCells);
 
  //cout << " first: " << data[0] << endl; 

  FIELD_3D zCompMAC(uEnd, vEnd, wEnd + 1);
  FIELD_3D yCompMAC(uEnd, vEnd + 1, wEnd);
  FIELD_3D xCompMAC(uEnd + 1, vEnd, wEnd);

  cout << " uRes: " << uEnd << " vEnd: " << vEnd << " wEnd: " << wEnd << endl;

  // Z first!?
  //
  // TK: this is correct -- PhysBAM does a zyx (reverse) ordering for some reason

  for (int x = 0; x < scanTotalZ; x++)
    zCompMAC[x] = data[x];

  for (int x = 0; x < scanTotalY; x++)
    yCompMAC[x] = data[scanTotalZ + x];
  
  for (int x = 0; x < scanTotalX; x++)
    xCompMAC[x] = data[scanTotalZ + scanTotalY + x];

  FIELD_3D xComp(uEnd, vEnd, wEnd);
  FIELD_3D yComp(uEnd, vEnd, wEnd);
  FIELD_3D zComp(uEnd, vEnd, wEnd);

  for (int z = 0; z < wEnd; z++)
    for (int y = 0; y < vEnd; y++)
      for (int x = 0; x < uEnd; x++)
        xComp(x,y,z) = (xCompMAC(x,y,z) + xCompMAC(x+1,y,z)) * 0.5;
  
  for (int z = 0; z < wEnd; z++)
    for (int y = 0; y < vEnd; y++)
      for (int x = 0; x < uEnd; x++)
        yComp(x,y,z) = (yCompMAC(x,y,z) + yCompMAC(x,y+1,z)) * 0.5;
  
  for (int z = 0; z < wEnd; z++)
    for (int y = 0; y < vEnd; y++)
      for (int x = 0; x < uEnd; x++)
        zComp(x,y,z) = (zCompMAC(x,y,z) + zCompMAC(x,y,z+1)) * 0.5;

  _xRes = uEnd;
  _yRes = vEnd;
  _zRes = wEnd;

  // set longest dimension to 1
  int biggest = (_xRes > _yRes) ? _xRes : _yRes;
  biggest = (_zRes > biggest) ? _zRes : biggest;
  _lengths[0] = (Real)_xRes / biggest;
  _lengths[1] = (Real)_yRes / biggest;
  _lengths[2] = (Real)_zRes / biggest;

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;

  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;

  if (_data) delete[] _data;
  _data = new VEC3F[_totalCells];

  //int index = 0;
  int index = 0;
  for (int z = 0; z < wEnd; z++)
    for (int y = 0; y < vEnd; y++)
      for (int x = 0; x < uEnd; x++, index++)
      {
        _data[index][0] = xComp(x,y,z);
        _data[index][1] = yComp(x,y,z);
        _data[index][2] = zComp(x,y,z);
      }

  delete[] data;
  gzclose(file);

  cout << " done. " << endl;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::readHoudiniGz(const string filename)
{
  cout << " Reading Houdini velocity file ... "; flush(cout);

  // unzip the file
  string gunzip = string("gunzip ") + filename + string(".gz");
  cout << " Unzipping ... "; flush(cout);
  system(gunzip.c_str());

  readHoudini(filename);

  string gzip = string("gzip ") + filename;
  cout << " Re-zipping ... "; flush(cout);
  system(gzip.c_str());
  cout << " done. " << endl;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::readHoudini(const string filename)
{
  FIELD_3D xCompMAC;
  FIELD_3D yCompMAC;
  FIELD_3D zCompMAC;

  FIELD_3D::readHoudiniVel(filename, xCompMAC, yCompMAC, zCompMAC);

  int uEnd = zCompMAC.xRes();
  int vEnd = zCompMAC.yRes();
  int wEnd = yCompMAC.zRes();

  FIELD_3D xComp(uEnd, vEnd, wEnd);
  FIELD_3D yComp(uEnd, vEnd, wEnd);
  FIELD_3D zComp(uEnd, vEnd, wEnd);

  for (int z = 0; z < wEnd; z++)
    for (int y = 0; y < vEnd; y++)
      for (int x = 0; x < uEnd; x++)
        xComp(x,y,z) = (xCompMAC(x,y,z) + xCompMAC(x+1,y,z)) * 0.5;
  
  for (int z = 0; z < wEnd; z++)
    for (int y = 0; y < vEnd; y++)
      for (int x = 0; x < uEnd; x++)
        yComp(x,y,z) = (yCompMAC(x,y,z) + yCompMAC(x,y+1,z)) * 0.5;
  
  for (int z = 0; z < wEnd; z++)
    for (int y = 0; y < vEnd; y++)
      for (int x = 0; x < uEnd; x++)
        zComp(x,y,z) = (zCompMAC(x,y,z) + zCompMAC(x,y,z+1)) * 0.5;

  _xRes = uEnd;
  _yRes = vEnd;
  _zRes = wEnd;

  _lengths[0] = 1;
  _lengths[1] = 1;
  _lengths[2] = 1;
  
  _dx = 1.0 / _xRes;
  _dy = 1.0 / _yRes;
  _dz = 1.0 / _zRes;

  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;

  if (_data) delete[] _data;
  _data = new VEC3F[_totalCells];

  int index = 0;
  for (int z = 0; z < wEnd; z++)
    for (int y = 0; y < vEnd; y++)
      for (int x = 0; x < uEnd; x++, index++)
      {
        _data[index][0] = xComp(x,y,z);
        _data[index][1] = yComp(x,y,z);
        _data[index][2] = zComp(x,y,z);
      }
}

///////////////////////////////////////////////////////////////////////
// flip the x and y coordinates
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::flipXY() const
{
  VECTOR3_FIELD_3D final(_yRes, _xRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        final(y,x,z) = (*this)(x,y,z);

  return final;
}

///////////////////////////////////////////////////////////////////////
// flip the x and y coordinates
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::flipXZ() const
{
  VECTOR3_FIELD_3D final(_zRes, _yRes, _xRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        final(z,y,x) = (*this)(x,y,z);

        Real temp = final(z,y,x)[0];
        final(z,y,x)[0] = final(z,y,x)[2];
        final(z,y,x)[2] = temp;
      }

  return final;
}

///////////////////////////////////////////////////////////////////////
// flip the z and y coordinates
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::flipZY() const
{
  VECTOR3_FIELD_3D final(_xRes, _zRes, _yRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        final(x,z,y) = (*this)(x,y,z);

  return final;
}

///////////////////////////////////////////////////////////////////////
// advect using first order semi-Lagrangian
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::advect(const Real dt, const VECTOR3_FIELD_3D& velocityGrid, const VECTOR3_FIELD_3D& oldField, VECTOR3_FIELD_3D& newField)
{
  const int xRes = velocityGrid.xRes();
  const int yRes = velocityGrid.yRes();
  const int zRes = velocityGrid.zRes();
  const int slabSize = velocityGrid.slabSize();

  // scale dt up to grid resolution
#pragma omp parallel
#pragma omp for schedule(static)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        const int index = x + y * xRes + z * slabSize;
        
        // backtrace
        const VEC3F velocity = velocityGrid[index];
        Real xTrace = x - dt * velocity[0];
        Real yTrace = y - dt * velocity[1];
        Real zTrace = z - dt * velocity[2];

        // clamp backtrace to grid boundaries
        xTrace = (xTrace < 1.5) ? 1.5 : xTrace;
        xTrace = (xTrace > xRes - 2.5) ? xRes - 2.5 : xTrace;
        yTrace = (yTrace < 1.5) ? 1.5 : yTrace;
        yTrace = (yTrace > yRes - 2.5) ? yRes - 2.5 : yTrace;
        zTrace = (zTrace < 1.5) ? 1.5 : zTrace;
        zTrace = (zTrace > zRes - 2.5) ? zRes - 2.5 : zTrace;

        /*
        // clamp backtrace to grid boundaries
        xTrace = (xTrace < 0.5) ? 0.5 : xTrace;
        xTrace = (xTrace > xRes - 1.5) ? xRes - 1.5 : xTrace;
        yTrace = (yTrace < 0.5) ? 0.5 : yTrace;
        yTrace = (yTrace > yRes - 1.5) ? yRes - 1.5 : yTrace;
        zTrace = (zTrace < 0.5) ? 0.5 : zTrace;
        zTrace = (zTrace > zRes - 1.5) ? zRes - 1.5 : zTrace;
        */

        // locate neighbors to interpolate
        const int x0 = (int)xTrace;
        const int x1 = x0 + 1;
        const int y0 = (int)yTrace;
        const int y1 = y0 + 1;
        const int z0 = (int)zTrace;
        const int z1 = z0 + 1;

        // get interpolation weights
        const Real s1 = xTrace - x0;
        const Real s0 = 1.0f - s1;
        const Real t1 = yTrace - y0;
        const Real t0 = 1.0f - t1;
        const Real u1 = zTrace - z0;
        const Real u0 = 1.0f - u1;

        const int i000 = x0 + y0 * xRes + z0 * slabSize;
        const int i010 = x0 + y1 * xRes + z0 * slabSize;
        const int i100 = x1 + y0 * xRes + z0 * slabSize;
        const int i110 = x1 + y1 * xRes + z0 * slabSize;
        const int i001 = x0 + y0 * xRes + z1 * slabSize;
        const int i011 = x0 + y1 * xRes + z1 * slabSize;
        const int i101 = x1 + y0 * xRes + z1 * slabSize;
        const int i111 = x1 + y1 * xRes + z1 * slabSize;

        // interpolate
        // (indices could be computed once)
        newField[index] = u0 * (s0 * (t0 * oldField[i000] +
                                      t1 * oldField[i010]) +
                                s1 * (t0 * oldField[i100] +
                                      t1 * oldField[i110])) +
                          u1 * (s0 * (t0 * oldField[i001] +
                                      t1 * oldField[i011]) +
                                s1 * (t0 * oldField[i101] +
                                      t1 * oldField[i111]));
      }
}

///////////////////////////////////////////////////////////////////////
// advect using first order semi-Lagrangian
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::advect(const Real dt, const VECTOR3_FIELD_3D& velocityGrid, const FIELD_3D& oldField, FIELD_3D& newField)
{
  const int xRes = velocityGrid.xRes();
  const int yRes = velocityGrid.yRes();
  const int zRes = velocityGrid.zRes();
  const int slabSize = velocityGrid.slabSize();

  // scale dt up to grid resolution
#pragma omp parallel
#pragma omp for schedule(static)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        const int index = x + y * xRes + z * slabSize;
        
        // backtrace
        const VEC3F velocity = velocityGrid[index];
        Real xTrace = x - dt * velocity[0];
        Real yTrace = y - dt * velocity[1];
        Real zTrace = z - dt * velocity[2];

        // clamp backtrace to grid boundaries
        xTrace = (xTrace < 1.5) ? 1.5 : xTrace;
        xTrace = (xTrace > xRes - 2.5) ? xRes - 2.5 : xTrace;
        yTrace = (yTrace < 1.5) ? 1.5 : yTrace;
        yTrace = (yTrace > yRes - 2.5) ? yRes - 2.5 : yTrace;
        zTrace = (zTrace < 1.5) ? 1.5 : zTrace;
        zTrace = (zTrace > zRes - 2.5) ? zRes - 2.5 : zTrace;

        /*
        // clamp backtrace to grid boundaries
        xTrace = (xTrace < 0.5) ? 0.5 : xTrace;
        xTrace = (xTrace > xRes - 1.5) ? xRes - 1.5 : xTrace;
        yTrace = (yTrace < 0.5) ? 0.5 : yTrace;
        yTrace = (yTrace > yRes - 1.5) ? yRes - 1.5 : yTrace;
        zTrace = (zTrace < 0.5) ? 0.5 : zTrace;
        zTrace = (zTrace > zRes - 1.5) ? zRes - 1.5 : zTrace;
        */

        // locate neighbors to interpolate
        const int x0 = (int)xTrace;
        const int x1 = x0 + 1;
        const int y0 = (int)yTrace;
        const int y1 = y0 + 1;
        const int z0 = (int)zTrace;
        const int z1 = z0 + 1;

        // get interpolation weights
        const Real s1 = xTrace - x0;
        const Real s0 = 1.0f - s1;
        const Real t1 = yTrace - y0;
        const Real t0 = 1.0f - t1;
        const Real u1 = zTrace - z0;
        const Real u0 = 1.0f - u1;

        const int i000 = x0 + y0 * xRes + z0 * slabSize;
        const int i010 = x0 + y1 * xRes + z0 * slabSize;
        const int i100 = x1 + y0 * xRes + z0 * slabSize;
        const int i110 = x1 + y1 * xRes + z0 * slabSize;
        const int i001 = x0 + y0 * xRes + z1 * slabSize;
        const int i011 = x0 + y1 * xRes + z1 * slabSize;
        const int i101 = x1 + y0 * xRes + z1 * slabSize;
        const int i111 = x1 + y1 * xRes + z1 * slabSize;

        // interpolate
        // (indices could be computed once)
        newField[index] = u0 * (s0 * (t0 * oldField[i000] +
                                      t1 * oldField[i010]) +
                                s1 * (t0 * oldField[i100] +
                                      t1 * oldField[i110])) +
                          u1 * (s0 * (t0 * oldField[i001] +
                                      t1 * oldField[i011]) +
                                s1 * (t0 * oldField[i101] +
                                      t1 * oldField[i111]));
      }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::advectMacCormack(const Real dt, const VECTOR3_FIELD_3D& velocityGrid, FIELD_3D& oldField, FIELD_3D& newField, FIELD_3D& temp1, FIELD_3D& temp2)
{
	FIELD_3D& phiHatN  = temp1;
	FIELD_3D& phiHatN1 = temp2;

	const int sx = oldField.xRes();
	const int sy = oldField.yRes();
	const int sz = oldField.zRes();

	for (int x = 0; x < sx * sy * sz; x++)
		phiHatN[x] = phiHatN1[x] = oldField[x];

	FIELD_3D& phiN    = oldField;
	FIELD_3D& phiN1   = newField;

	// phiHatN1 = A(phiN)
	advect(dt, velocityGrid, phiN, phiHatN1);

	// phiHatN = A^R(phiHatN1)
	advect(-1.0 * dt, velocityGrid, phiHatN1, phiHatN);

	// phiN1 = phiHatN1 + (phiN - phiHatN) / 2
	const int border = 0; 
	for (int z = border; z < sz - border; z++)
		for (int y = border; y < sy - border; y++)
			for (int x = border; x < sx - border; x++) {
				int index = x + y * sx + z * sx*sy;
				phiN1[index] = phiHatN1[index] + (phiN[index] - phiHatN[index]) * 0.50f;
			}

  phiN1.copyBorderAll();

	// clamp any newly created extrema
	clampExtrema(dt, velocityGrid, oldField, newField);

	// if the error estimate was bad, revert to first order
	clampOutsideRays(dt, velocityGrid, oldField, phiHatN1, newField);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::advectMacCormack(const Real dt, const VECTOR3_FIELD_3D& velocityGrid, VECTOR3_FIELD_3D& oldField, VECTOR3_FIELD_3D& newField, VECTOR3_FIELD_3D& temp1, VECTOR3_FIELD_3D& temp2)
{
	VECTOR3_FIELD_3D& phiHatN  = temp1;
	VECTOR3_FIELD_3D& phiHatN1 = temp2;

	const int sx = oldField.xRes();
	const int sy = oldField.yRes();
	const int sz = oldField.zRes();

	for (int x = 0; x < sx * sy * sz; x++)
		phiHatN[x] = phiHatN1[x] = oldField[x];

	VECTOR3_FIELD_3D& phiN    = oldField;
	VECTOR3_FIELD_3D& phiN1   = newField;

	// phiHatN1 = A(phiN)
	//advect(dt, velocity field, old field, new field)
	advect(dt, velocityGrid, phiN, phiHatN1);

	// phiHatN = A^R(phiHatN1)
	advect(-1.0 * dt, velocityGrid, phiHatN1, phiHatN);

	// phiN1 = phiHatN1 + (phiN - phiHatN) / 2
	const int border = 0; 
	for (int z = border; z < sz - border; z++)
		for (int y = border; y < sy - border; y++)
			for (int x = border; x < sx - border; x++) {
				int index = x + y * sx + z * sx*sy;
				phiN1[index] = phiHatN1[index] + (phiN[index] - phiHatN[index]) * 0.5;
			}

  phiN1.copyBorderAll();

	// clamp any newly created extrema
	//clampExtrema(dt, velocityGrid, oldField, newField);

	// if the error estimate was bad, revert to first order
	//clampOutsideRays(dt, velocityGrid, oldField, phiHatN1, newField);
}

///////////////////////////////////////////////////////////////////////
// Clamp the extrema generated by the BFECC error correction
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::clampExtrema(const Real dt, const VECTOR3_FIELD_3D& velocityField, const VECTOR3_FIELD_3D& oldField, VECTOR3_FIELD_3D& newField)
{
	const int xRes = velocityField.xRes();
	const int yRes = velocityField.yRes();
	const int zRes = velocityField.zRes();
	const int slabSize = velocityField.slabSize();

	for (int z = 1; z < zRes-1; z++)
		for (int y = 1; y < yRes-1; y++)
			for (int x = 1; x < xRes-1; x++)
			{
				const int index = x + y * xRes + z * slabSize;
				// backtrace
        const VEC3F velocity = velocityField[index];
        Real xTrace = x - dt * velocity[0];
        Real yTrace = y - dt * velocity[1];
        Real zTrace = z - dt * velocity[2];

				// clamp backtrace to grid boundaries
        xTrace = (xTrace < 0.5) ? 0.5 : xTrace;
        xTrace = (xTrace > xRes - 1.5) ? xRes - 1.5 : xTrace;
        yTrace = (yTrace < 0.5) ? 0.5 : yTrace;
        yTrace = (yTrace > yRes - 1.5) ? yRes - 1.5 : yTrace;
        zTrace = (zTrace < 0.5) ? 0.5 : zTrace;
        zTrace = (zTrace > zRes - 1.5) ? zRes - 1.5 : zTrace;

				// locate neighbors to interpolate
				const int x0 = (int)xTrace;
				const int x1 = x0 + 1;
				const int y0 = (int)yTrace;
				const int y1 = y0 + 1;
				const int z0 = (int)zTrace;
				const int z1 = z0 + 1;

				const int i000 = x0 + y0 * xRes + z0 * slabSize;
				const int i010 = x0 + y1 * xRes + z0 * slabSize;
				const int i100 = x1 + y0 * xRes + z0 * slabSize;
				const int i110 = x1 + y1 * xRes + z0 * slabSize;
				const int i001 = x0 + y0 * xRes + z1 * slabSize;
				const int i011 = x0 + y1 * xRes + z1 * slabSize;
				const int i101 = x1 + y0 * xRes + z1 * slabSize;
				const int i111 = x1 + y1 * xRes + z1 * slabSize;

        /*
        //if (x == 15 && y == 1 && z == 24)
        if (x == 28 && y == 1 && z == 24)
        {
          cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
          cout << " Before: " << newField[index] << endl;
          cout << " olds: " << endl;
          cout << oldField[i000] << endl;
          cout << oldField[i001] << endl;
          cout << oldField[i010] << endl;
          cout << oldField[i011] << endl;
          cout << oldField[i100] << endl;
          cout << oldField[i101] << endl;
          cout << oldField[i110] << endl;
          cout << oldField[i111] << endl;
          cout << " indices: " << endl;
          cout << x0 << " " << y0 << " " << z0 << endl;
          cout << x0 << " " << y0 << " " << z1 << endl;
          cout << x0 << " " << y1 << " " << z0 << endl;
          cout << x0 << " " << y1 << " " << z1 << endl;
          cout << x1 << " " << y0 << " " << z0 << endl;
          cout << x1 << " " << y0 << " " << z1 << endl;
          cout << x1 << " " << y1 << " " << z0 << endl;
          cout << x1 << " " << y1 << " " << z1 << endl;
        }
        */

        for (int i = 0; i < 3; i++)
        {
          Real minField = oldField[i000][i];
          Real maxField = oldField[i000][i];

          minField = (oldField[i010][i] < minField) ? oldField[i010][i] : minField;
          maxField = (oldField[i010][i] > maxField) ? oldField[i010][i] : maxField;

          minField = (oldField[i100][i] < minField) ? oldField[i100][i] : minField;
          maxField = (oldField[i100][i] > maxField) ? oldField[i100][i] : maxField;

          minField = (oldField[i110][i] < minField) ? oldField[i110][i] : minField;
          maxField = (oldField[i110][i] > maxField) ? oldField[i110][i] : maxField;

          minField = (oldField[i001][i] < minField) ? oldField[i001][i] : minField;
          maxField = (oldField[i001][i] > maxField) ? oldField[i001][i] : maxField;

          minField = (oldField[i011][i] < minField) ? oldField[i011][i] : minField;
          maxField = (oldField[i011][i] > maxField) ? oldField[i011][i] : maxField;

          minField = (oldField[i101][i] < minField) ? oldField[i101][i] : minField;
          maxField = (oldField[i101][i] > maxField) ? oldField[i101][i] : maxField;

          minField = (oldField[i111][i] < minField) ? oldField[i111][i] : minField;
          maxField = (oldField[i111][i] > maxField) ? oldField[i111][i] : maxField;

          newField[index][i] = (newField[index][i] > maxField) ? maxField : newField[index][i];
          newField[index][i] = (newField[index][i] < minField) ? minField : newField[index][i];
        }
        
        /*
        //if (x == 15 && y == 1 && z == 24)
        if (x == 28 && y == 1 && z == 24)
        {
          cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
          cout << " After: " << newField[index] << endl;
        }
        */
			}
}

///////////////////////////////////////////////////////////////////////
// Clamp the extrema generated by the BFECC error correction
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::clampExtrema(const Real dt, const VECTOR3_FIELD_3D& velocityField, const FIELD_3D& oldField, FIELD_3D& newField)
{
	const int xRes = velocityField.xRes();
	const int yRes = velocityField.yRes();
	const int zRes = velocityField.zRes();
	const int slabSize = velocityField.slabSize();

	for (int z = 1; z < zRes-1; z++)
		for (int y = 1; y < yRes-1; y++)
			for (int x = 1; x < xRes-1; x++)
			{
				const int index = x + y * xRes + z * slabSize;
				// backtrace
        const VEC3F velocity = velocityField[index];
        Real xTrace = x - dt * velocity[0];
        Real yTrace = y - dt * velocity[1];
        Real zTrace = z - dt * velocity[2];

				// clamp backtrace to grid boundaries
        xTrace = (xTrace < 0.5) ? 0.5 : xTrace;
        xTrace = (xTrace > xRes - 1.5) ? xRes - 1.5 : xTrace;
        yTrace = (yTrace < 0.5) ? 0.5 : yTrace;
        yTrace = (yTrace > yRes - 1.5) ? yRes - 1.5 : yTrace;
        zTrace = (zTrace < 0.5) ? 0.5 : zTrace;
        zTrace = (zTrace > zRes - 1.5) ? zRes - 1.5 : zTrace;

				// locate neighbors to interpolate
				const int x0 = (int)xTrace;
				const int x1 = x0 + 1;
				const int y0 = (int)yTrace;
				const int y1 = y0 + 1;
				const int z0 = (int)zTrace;
				const int z1 = z0 + 1;

				const int i000 = x0 + y0 * xRes + z0 * slabSize;
				const int i010 = x0 + y1 * xRes + z0 * slabSize;
				const int i100 = x1 + y0 * xRes + z0 * slabSize;
				const int i110 = x1 + y1 * xRes + z0 * slabSize;
				const int i001 = x0 + y0 * xRes + z1 * slabSize;
				const int i011 = x0 + y1 * xRes + z1 * slabSize;
				const int i101 = x1 + y0 * xRes + z1 * slabSize;
				const int i111 = x1 + y1 * xRes + z1 * slabSize;

				Real minField = oldField[i000];
				Real maxField = oldField[i000];

				minField = (oldField[i010] < minField) ? oldField[i010] : minField;
				maxField = (oldField[i010] > maxField) ? oldField[i010] : maxField;

				minField = (oldField[i100] < minField) ? oldField[i100] : minField;
				maxField = (oldField[i100] > maxField) ? oldField[i100] : maxField;

				minField = (oldField[i110] < minField) ? oldField[i110] : minField;
				maxField = (oldField[i110] > maxField) ? oldField[i110] : maxField;

				minField = (oldField[i001] < minField) ? oldField[i001] : minField;
				maxField = (oldField[i001] > maxField) ? oldField[i001] : maxField;

				minField = (oldField[i011] < minField) ? oldField[i011] : minField;
				maxField = (oldField[i011] > maxField) ? oldField[i011] : maxField;

				minField = (oldField[i101] < minField) ? oldField[i101] : minField;
				maxField = (oldField[i101] > maxField) ? oldField[i101] : maxField;

				minField = (oldField[i111] < minField) ? oldField[i111] : minField;
				maxField = (oldField[i111] > maxField) ? oldField[i111] : maxField;

				newField[index] = (newField[index] > maxField) ? maxField : newField[index];
				newField[index] = (newField[index] < minField) ? minField : newField[index];
			}
}

//////////////////////////////////////////////////////////////////////
// Reverts any backtraces that go into boundaries back to first 
// order -- in this case the error correction term was totally
// incorrect
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::clampOutsideRays(const Real dt, const VECTOR3_FIELD_3D& velocityField, const VECTOR3_FIELD_3D& oldField, const VECTOR3_FIELD_3D& oldAdvection, VECTOR3_FIELD_3D& newField)
{
	const int sx = velocityField.xRes();
	const int sy = velocityField.yRes();
	const int sz = velocityField.zRes();
	const int slabSize = velocityField.slabSize();

	for (int z = 1; z < sz - 1; z++)
		for (int y = 1; y < sy - 1; y++)
			for (int x = 1; x < sx - 1; x++)
			{
				const int index = x + y * sx + z * slabSize;
				// backtrace
        VEC3F velocity = velocityField[index];
        velocity *= dt;
				float xBackward = x + velocity[0];
				float yBackward = y + velocity[1];
				float zBackward = z + velocity[2];

				float xTrace    = x - velocity[0];
				float yTrace    = y - velocity[1];
				float zTrace    = z - velocity[2];

				// see if it goes outside the boundaries
				bool hasObstacle = 
					(zTrace < 1.0f)    || (zTrace > sz - 2.0f) ||
					(yTrace < 1.0f)    || (yTrace > sy - 2.0f) ||
					(xTrace < 1.0f)    || (xTrace > sx - 2.0f) ||
					(zBackward < 1.0f) || (zBackward > sz - 2.0f) ||
					(yBackward < 1.0f) || (yBackward > sy - 2.0f) ||
					(xBackward < 1.0f) || (xBackward > sx - 2.0f);

        if (x == 46 && y == 60 && z == 24)
        {
          cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
          cout << " Before:   " << newField[index] << endl;
          cout << " trace:    " << xTrace << " " << yTrace << " " << zTrace << endl;
          cout << " backward: " << xBackward << " " << yBackward << " " << zBackward << endl;
          cout << " obstacle: " << hasObstacle << endl;
          cout << (zTrace < 1.0f)   << endl; 
          cout << (zTrace > sz - 2.0f) << endl; 
          cout << (yTrace < 1.0f)    << endl; 
          cout << (yTrace > sy - 2.0f) << endl; 
          cout << (xTrace < 1.0f)    << endl; 
          cout << (xTrace > sx - 2.0f) << endl; 
          cout << (zBackward < 1.0f) << endl; 
          cout << (zBackward > sz - 2.0f) << endl; 
          cout << (yBackward < 1.0f) << endl; 
          cout << (yBackward > sy - 2.0f) << endl; 
          cout << (xBackward < 1.0f) << endl; 
          //cout << (xBackward > sx - 2.0f) << endl; 
          cout << (xBackward - (sx - 2.0f)) << endl;
        }

				// reuse old advection instead of doing another one...
				if(hasObstacle) 
          newField[index] = oldAdvection[index];

        if (x == 46 && y == 60 && z == 24)
        {
          cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
          cout << " After: " << newField[index] << endl;
        }
			}
}

//////////////////////////////////////////////////////////////////////
// Reverts any backtraces that go into boundaries back to first 
// order -- in this case the error correction term was totally
// incorrect
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::clampOutsideRays(const Real dt, const VECTOR3_FIELD_3D& velocityField, const FIELD_3D& oldField, const FIELD_3D& oldAdvection, FIELD_3D& newField)
{
	const int sx = velocityField.xRes();
	const int sy = velocityField.yRes();
	const int sz = velocityField.zRes();
	const int slabSize = velocityField.slabSize();

	for (int z = 1; z < sz - 1; z++)
		for (int y = 1; y < sy - 1; y++)
			for (int x = 1; x < sx - 1; x++)
			{
				const int index = x + y * sx + z * slabSize;
				// backtrace
        VEC3F velocity = velocityField[index];
        velocity *= dt;
				float xBackward = x + velocity[0];
				float yBackward = y + velocity[1];
				float zBackward = z + velocity[2];

				float xTrace    = x - velocity[0];
				float yTrace    = x - velocity[1];
				float zTrace    = x - velocity[2];

				// see if it goes outside the boundaries
				bool hasObstacle = 
					(zTrace < 1.0f)    || (zTrace > sz - 2.0f) ||
					(yTrace < 1.0f)    || (yTrace > sy - 2.0f) ||
					(xTrace < 1.0f)    || (xTrace > sx - 2.0f) ||
					(zBackward < 1.0f) || (zBackward > sz - 2.0f) ||
					(yBackward < 1.0f) || (yBackward > sy - 2.0f) ||
					(xBackward < 1.0f) || (xBackward > sx - 2.0f);

				// reuse old advection instead of doing another one...
				if(hasObstacle) { newField[index] = oldAdvection[index]; }
			}
}

///////////////////////////////////////////////////////////////////////
// real-valued cell center coordinates
///////////////////////////////////////////////////////////////////////
VEC3F VECTOR3_FIELD_3D::cellCenter(int x, int y, int z) const
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
// draw bounding box to OpenGL
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::drawBoundingBox()
{
  glPushMatrix();
    VEC3F v000 = cellCenter(0,0,0);
    VEC3F v001 = cellCenter(0,0,_zRes);
    VEC3F v010 = cellCenter(0,_yRes,0);
    VEC3F v011 = cellCenter(0,_yRes,_zRes);
    VEC3F v100 = cellCenter(_xRes,0,0);
    VEC3F v101 = cellCenter(_xRes,0,_zRes);
    VEC3F v110 = cellCenter(_xRes,_yRes,0);
    VEC3F v111 = cellCenter(_xRes, _yRes, _zRes);

    glColor4f(1,0,0,1);

    glBegin(GL_LINES);
#ifdef SINGLE_PRECISION
      glVertex3fv(v000);
      glVertex3fv(v001);

      glVertex3fv(v001);
      glVertex3fv(v011);
      
      glVertex3fv(v011);
      glVertex3fv(v010);

      glVertex3fv(v010);
      glVertex3fv(v000);

      glVertex3fv(v110);
      glVertex3fv(v111);
      
      glVertex3fv(v111);
      glVertex3fv(v101);

      glVertex3fv(v101);
      glVertex3fv(v100);
      
      glVertex3fv(v100);
      glVertex3fv(v110);

      glVertex3fv(v010);
      glVertex3fv(v110);

      glVertex3fv(v000);
      glVertex3fv(v100);
      
      glVertex3fv(v111);
      glVertex3fv(v011);
      
      glVertex3fv(v101);
      glVertex3fv(v001);
#else
      glVertex3dv(v000);
      glVertex3dv(v001);

      glVertex3dv(v001);
      glVertex3dv(v011);
      
      glVertex3dv(v011);
      glVertex3dv(v010);

      glVertex3dv(v010);
      glVertex3dv(v000);

      glVertex3dv(v110);
      glVertex3dv(v111);
      
      glVertex3dv(v111);
      glVertex3dv(v101);

      glVertex3dv(v101);
      glVertex3dv(v100);
      
      glVertex3dv(v100);
      glVertex3dv(v110);

      glVertex3dv(v010);
      glVertex3dv(v110);

      glVertex3dv(v000);
      glVertex3dv(v100);
      
      glVertex3dv(v111);
      glVertex3dv(v011);
      
      glVertex3dv(v101);
      glVertex3dv(v001);
#endif
    glEnd();
  glPopMatrix();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::resizeAndWipe(int xRes, int yRes, int zRes, const VEC3F& center, const VEC3F& lengths)
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
  _lengths = lengths;
  _center = center;

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;

  try {
    _data = new VEC3F[_totalCells];
  }
  catch(std::bad_alloc& exc)
  {
    cout << " Failed to allocate " << _xRes << " " << _yRes << " " << _zRes << " VECTOR3_FIELD_3D!" << endl;
    double bytes = _xRes * _yRes * _zRes * sizeof(VEC3F);
    cout <<  bytes / pow(2.0,20.0) << " MB needed" << endl;
    exit(0);
  }

  _initialized = true;
}

///////////////////////////////////////////////////////////////////////
// normalize all the vectors in the field
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::normalize()
{
  for (int x = 0; x < _totalCells; x++)
    _data[x].normalize();
}

///////////////////////////////////////////////////////////////////////
// copy values out into the border, assuming that "borderSize" is the 
// width of the grid padding
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::copyIntoBorder(int borderSize)
{
  TIMER functionTimer(__FUNCTION__);
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        if (x == borderSize)
        {
          VEC3F value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x - i,y,z) = value;
        }
        if (x == _xRes - 1 - borderSize)
        {
          VEC3F value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x + i,y,z) = value;
        }              
        if (y == borderSize)
        {
          VEC3F value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x,y - i,z) = value;
        }
        if (y == _yRes - 1 - borderSize)
        {
          VEC3F value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x,y+i,z) = value;
        }              
        if (z == borderSize)
        {
          VEC3F value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x,y,z - i) = value;
        }
        if (z == _zRes - 1 - borderSize)
        {
          VEC3F value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x,y,z+i) = value;
        }

        // handle the corners
        if (x == borderSize && z == borderSize)
        {
          VEC3F value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x - i,y,z-j) = value;
        }
        if (x == _xRes - 1 - borderSize && z == _zRes - 1 - borderSize)
        {
          VEC3F value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x + i,y,z+j) = value;
        }

        if (z == borderSize && y == borderSize)
        {
          VEC3F value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x,y - j,z -i) = value;
        }
        if (z == _xRes - 1 - borderSize && y == _yRes - 1 - borderSize)
        {
          VEC3F value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x,y + j,z+i) = value;
        }

      }
}

///////////////////////////////////////////////////////////////////////
// pass back a field with a new padding of size "paddingSize"
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::withAddedPadding(int paddingSize) const
{
  // new length, with padding
  VEC3F newLengths = _lengths;
  newLengths[0] += paddingSize * 2 * _dx;
  newLengths[1] += paddingSize * 2 * _dx;
  newLengths[2] += paddingSize * 2 * _dx;

  VECTOR3_FIELD_3D final(_xRes + 2 * paddingSize, 
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

//////////////////////////////////////////////////////////////////////
// set everything on the border to zero
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::setZeroBorder()
{
  TIMER functionTimer(__FUNCTION__);
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
// set x direction to zero
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::setZeroX()
{
	int index;
	for (int z = 0; z < _zRes; z++)
		for (int y = 0; y < _yRes; y++)
		{
			// left slab
			index = y * _xRes + z * _slabSize;
			_data[index][0] = 0.0f;

			// right slab
			index += _xRes - 1;
			_data[index][0] = 0.0f;
		}
}

//////////////////////////////////////////////////////////////////////
// set y direction to zero
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::setZeroY()
{
	int index;
	for (int z = 0; z < _zRes; z++)
		for (int x = 0; x < _xRes; x++)
		{
			// bottom slab
			index = x + z * _slabSize;
			_data[index][1] = 0.0f;

			// top slab
			index += _slabSize - _xRes;
			_data[index][1] = 0.0f;
		}
}

//////////////////////////////////////////////////////////////////////
// set z direction to zero
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::setZeroZ()
{
	int index;
	for (int y = 0; y < _yRes; y++)
		for (int x = 0; x < _xRes; x++)
		{
			// front slab
			index = x + y * _xRes;
			_data[index][2] = 0.0f;

			// back slab
			index += _totalCells - _slabSize;
			_data[index][2] = 0.0f;
		}
}

//////////////////////////////////////////////////////////////////////
// set x direction to Neumann boundary conditions
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::setNeumannX()
{
	int index;
	for (int z = 0; z < _zRes; z++)
		for (int y = 0; y < _yRes; y++)
		{
			// left slab
			index = y * _xRes + z * _slabSize;
			_data[index][0] = _data[index + 2][0];

			// right slab
			index += _xRes - 1;
			_data[index][0] = _data[index - 2][0];
		}
}


//////////////////////////////////////////////////////////////////////
// set y direction to Neumann boundary conditions
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::setNeumannY()
{
	int index;
	for (int z = 0; z < _zRes; z++)
		for (int x = 0; x < _xRes; x++)
		{
			// bottom slab
			index = x + z * _slabSize;
			_data[index][1] = _data[index + 2 * _xRes][1];

			// top slab
			index += _slabSize - _xRes;
			_data[index][1] = _data[index - 2 * _xRes][1];
		}

  // TK: Putting this here is a terrible idea!
  //
  // fix, force top slab to only allow outwards flux
  /*
	for (int z = 0; z < _zRes; z++)
		for (int x = 0; x < _xRes; x++)
		{
			// top slab
			int index = x + z * _slabSize;

			index += _slabSize - _xRes;
			if(_data[index][1] < 0.) _data[index][1] = 0.;

			index -= _xRes;
			if(_data[index][1] < 0.) _data[index][1] = 0.;
		}
    */
}

//////////////////////////////////////////////////////////////////////
// set z direction to Neumann boundary conditions
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::setNeumannZ()
{
	int index;
	for (int y = 0; y < _yRes; y++)
		for (int x = 0; x < _xRes; x++)
		{
			// front slab
			index = x + y * _xRes;
			_data[index][2] = _data[index + 2 * _slabSize][2];

			// back slab
			index += _totalCells - _slabSize;
			_data[index][2] = _data[index - 2 * _slabSize][2];
		}
}

//////////////////////////////////////////////////////////////////////
// copy grid boundary
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::copyBorderX()
{
	int index;
	for (int z = 0; z < _zRes; z++)
		for (int y = 0; y < _yRes; y++)
		{
			// left slab
			index = y * _xRes + z * _slabSize;
			_data[index][0] = _data[index + 1][0];

			// right slab
			index += _xRes - 1;
			_data[index][0] = _data[index - 1][0];
		}
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::copyBorderY()
{
	int index;
	for (int z = 0; z < _zRes; z++)
		for (int x = 0; x < _xRes; x++)
		{
			// bottom slab
			index = x + z * _slabSize;
			_data[index][1] = _data[index + _xRes][1]; 
			// top slab
			index += _slabSize - _xRes;
			_data[index][1] = _data[index - _xRes][1];
		}
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::copyBorderZ()
{
	int index;
	for (int y = 0; y < _yRes; y++)
		for (int x = 0; x < _xRes; x++)
		{
			// front slab
			index = x + y * _xRes;
			_data[index][2] = _data[index + _slabSize][2]; 
			// back slab
			index += _totalCells - _slabSize;
			_data[index][2] = _data[index - _slabSize][2];
		}
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::copyBorderAll()
{
	int index;
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
}

void VECTOR3_FIELD_3D::setZeroSphere(const VEC3I& center, double radius) 
{
  assert ( radius < _lengths.maxElement() ); 
  VEC3F centerCoords = cellCenter(center[0], center[1], center[2]);
  int index = 0;
  for (int z = 0; z < _zRes; z++) {
    for (int y = 0; y < _yRes; y++) {
      for (int x = 0; x < _xRes; x++, index++) {
        if ( norm2(cellCenter(x, y, z) - centerCoords) < radius*radius ) {
          // this sets the entire VEC3F to (0, 0, 0)
          _data[index] = 0.0f;
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// BLAS-like interface, output += alpha * input
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::axpy(const Real& alpha, const VECTOR3_FIELD_3D& field)
{
  assert(field.xRes() == _xRes);
  assert(field.yRes() == _yRes);
  assert(field.zRes() == _zRes);
  for (int x = 0; x < _totalCells; x++)
    _data[x] += alpha * field[x];
}

//////////////////////////////////////////////////////////////////////
// swap the contents with another object
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::swapPointers(VECTOR3_FIELD_3D& field)
{
  assert(field.xRes() == _xRes);
  assert(field.yRes() == _yRes);
  assert(field.zRes() == _zRes);
  
  VEC3F* temp = _data;
  _data = field._data;
  field._data = temp;
}

//////////////////////////////////////////////////////////////////////
// return a flattened array of all the field contents
//////////////////////////////////////////////////////////////////////
VECTOR VECTOR3_FIELD_3D::flattened() const
{
  VECTOR final(_totalCells * 3);

  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        for (int i = 0; i < 3; i++, index++)
          final[index] = (*this)(x,y,z)[i];

  return final;
}

//////////////////////////////////////////////////////////////////////
// return a flattened array of all the field contents
//////////////////////////////////////////////////////////////////////
VectorXd VECTOR3_FIELD_3D::flattenedEigen() const
{
  VectorXd final(_totalCells * 3);

  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        for (int i = 0; i < 3; i++, index++)
          final(index) = (*this)(x,y,z)[i];

  return final;
}

//////////////////////////////////////////////////////////////////////
// return a flattened array of all the field contents
//////////////////////////////////////////////////////////////////////
VectorXd VECTOR3_FIELD_3D::flattenedEigenHomogeneous() const
{
  VectorXd final(_totalCells * 3 + 1);

  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        for (int i = 0; i < 3; i++, index++)
          final(index) = (*this)(x,y,z)[i];

  final[final.size() - 1] = 1.0;

  return final;
}

//////////////////////////////////////////////////////////////////////
// get gradient field
//////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::gradient()
{
  VECTOR3_FIELD_3D final(*this);
  final.clear();

  float halfInv = 1.0 / (2.0 * _dx);

  VECTOR3_FIELD_3D& input = *this;

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 1; x < _xRes - 1; x++)
        final(x,y,z)[0] = (input(x + 1,y,z)[0] - input(x - 1,y,z)[0]) * halfInv;

  for (int z = 0; z < _zRes; z++)
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 0; x < _xRes; x++)
        final(x,y,z)[1] = (input(x,y + 1,z)[1] - input(x,y - 1,z)[1]) * halfInv;

  for (int z = 1; z < _zRes - 1; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        final(x,y,z)[2] = (input(x,y,z + 1)[2] - input(x,y,z - 1)[2]) * halfInv;

  return final;
}

//////////////////////////////////////////////////////////////////////
// peel off the outer boundary of grid cells in preparation for PCA
//////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::peelBoundary() const
{
  VECTOR3_FIELD_3D final(_xRes - 2, _yRes - 2, _zRes - 2);

  for (int z = 1; z < _zRes - 1; z++)
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++)
        final(x - 1, y - 1, z - 1) = (*this)(x,y,z);

  return final;
}

//////////////////////////////////////////////////////////////////////
// peel off and return a single boundary
//////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::leftBoundary() const
{
  VECTOR3_FIELD_3D final(1, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      final(0,y,z) = (*this)(0,y,z);

  return final;
}

//////////////////////////////////////////////////////////////////////
// peel off and return a single boundary
//////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::rightBoundary() const
{
  VECTOR3_FIELD_3D final(1, _yRes, _zRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      final(0,y,z) = (*this)(_xRes - 1,y,z);

  return final;
}

//////////////////////////////////////////////////////////////////////
// peel off and return a single boundary
//////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::topBoundary() const
{
  VECTOR3_FIELD_3D final(_xRes - 2, 1, _zRes);

  for (int x = 1; x < _xRes - 1; x++)
    for (int z = 0; z < _zRes; z++)
      final(x - 1,0,z) = (*this)(x,_yRes - 1,z);

  return final;
}

//////////////////////////////////////////////////////////////////////
// peel off and return a single boundary
//////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::bottomBoundary() const
{
  VECTOR3_FIELD_3D final(_xRes - 2, 1, _zRes);

  for (int x = 1; x < _xRes - 1; x++)
    for (int z = 0; z < _zRes; z++)
      final(x - 1,0,z) = (*this)(x,0,z);

  return final;
}

//////////////////////////////////////////////////////////////////////
// peel off and return a single boundary
//////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::nearBoundary() const
{
  VECTOR3_FIELD_3D final(_xRes - 2, _yRes - 2, 1);

  for (int x = 1; x < _xRes - 1; x++)
    for (int y = 1; y < _yRes - 1; y++)
      final(x - 1,y - 1,0) = (*this)(x,y,0);

  return final;
}

//////////////////////////////////////////////////////////////////////
// peel off and return a single boundary
//////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::farBoundary() const
{
  VECTOR3_FIELD_3D final(_xRes - 2, _yRes - 2, 1);

  for (int x = 1; x < _xRes - 1; x++)
    for (int y = 1; y < _yRes - 1; y++)
      final(x - 1,y - 1,0) = (*this)(x,y,_zRes - 1);

  return final;
}

//////////////////////////////////////////////////////////////////////
// retrieve the boundaries in an array
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::getPeeledBoundaries(vector<VECTOR3_FIELD_3D>& boundaries) const
{
  boundaries.clear();

  boundaries.push_back(leftBoundary()); 
  boundaries.push_back(rightBoundary());
  boundaries.push_back(bottomBoundary());
  boundaries.push_back(topBoundary());
  boundaries.push_back(nearBoundary());
  boundaries.push_back(farBoundary());
}

//////////////////////////////////////////////////////////////////////
// reassemble a velocity field from a bunch of peeled boundaries 
// and a center
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::reassemblePeels(const VECTOR3_FIELD_3D& middle, const vector<VECTOR3_FIELD_3D>& boundaries)
{
  // do the middle
  assert(_xRes == middle.xRes() + 2);
  assert(_yRes == middle.yRes() + 2);
  assert(_zRes == middle.zRes() + 2);

  for (int z = 1; z < _zRes - 1; z++)
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++)
        (*this)(x,y,z) = middle(x - 1, y - 1, z - 1);

  // left boundary
  assert(boundaries[0].yRes() == _yRes);
  assert(boundaries[0].zRes() == _zRes);
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      (*this)(0,y,z) = boundaries[0](0,y,z);

  // right boundary
  assert(boundaries[1].yRes() == _yRes);
  assert(boundaries[1].zRes() == _zRes);
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      (*this)(_xRes - 1,y,z) = boundaries[1](0,y,z);

  // bottom boundary
  assert(boundaries[2].xRes() + 2 == _xRes);
  assert(boundaries[2].zRes() == _zRes);
  for (int x = 1; x < _xRes - 1; x++)
    for (int z = 0; z < _zRes; z++)
      (*this)(x,0,z) = boundaries[2](x - 1,0,z);

  // top boundary
  assert(boundaries[3].xRes() + 2 == _xRes);
  assert(boundaries[3].zRes() == _zRes);
  for (int x = 1; x < _xRes - 1; x++)
    for (int z = 0; z < _zRes; z++)
      (*this)(x,_yRes - 1,z) = boundaries[3](x - 1,0,z);

  // near boundary
  assert(boundaries[4].xRes() + 2 == _xRes);
  assert(boundaries[4].yRes() + 2 == _yRes);
  for (int x = 1; x < _xRes - 1; x++)
    for (int y = 1; y < _yRes - 1; y++)
      (*this)(x,y,0) = boundaries[4](x - 1, y - 1, 0);

  // far boundary
  assert(boundaries[5].xRes() + 2 == _xRes);
  assert(boundaries[5].yRes() + 2 == _yRes);
  for (int x = 1; x < _xRes - 1; x++)
    for (int y = 1; y < _yRes - 1; y++)
      (*this)(x,y,_zRes - 1) = boundaries[5](x - 1, y - 1, 0);
}

//////////////////////////////////////////////////////////////////////
// set the field innards to a peeled version
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::setWithPeeled(const VECTOR& data)
{
  TIMER functionTimer(__FUNCTION__);
  assert(data.size() == 3 * (_xRes - 2) * (_yRes - 2) * (_zRes - 2) ||
         data.size() == 3 * (_xRes - 2) * (_yRes - 2) * (_zRes - 2) + 1);

  int index = 0;
  for (int z = 1; z < _zRes - 1; z++)
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++, index++)
      {
        (*this)(x,y,z)[0] = data[3 * index];
        (*this)(x,y,z)[1] = data[3 * index + 1];
        (*this)(x,y,z)[2] = data[3 * index + 2];
      }
}

//////////////////////////////////////////////////////////////////////
// set the field innards to a peeled version
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::setWithPeeled(const VectorXd& data)
{
  TIMER functionTimer(__FUNCTION__);
  assert(data.size() == 3 * (_xRes - 2) * (_yRes - 2) * (_zRes - 2));

  int index = 0;
  for (int z = 1; z < _zRes - 1; z++)
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++, index++)
      {
        (*this)(x,y,z)[0] = data[3 * index];
        (*this)(x,y,z)[1] = data[3 * index + 1];
        (*this)(x,y,z)[2] = data[3 * index + 2];
      }
}

//////////////////////////////////////////////////////////////////////
// do a projection of the peeled field, using the passed in basis
//////////////////////////////////////////////////////////////////////
VectorXd VECTOR3_FIELD_3D::peeledProject(const MatrixXd& U)
{
  TIMER functionTimer(__FUNCTION__);

  /*
  int index = 0;
  for (int z = 1; z < _zRes - 1; z++)
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++, index++)
      {
        (*this)(x,y,z)[0] = data[3 * index];
        (*this)(x,y,z)[1] = data[3 * index + 1];
        (*this)(x,y,z)[2] = data[3 * index + 2];
      }
      */
  

  const int xPeeled = _xRes - 2;
  const int slabPeeled = (_xRes - 2) * (_yRes - 2);
  //const int totalThreads = omp_get_max_threads();
  const int totalThreads = 1; 

  vector<VectorXd> finals(totalThreads);
  for (unsigned int x = 0; x < finals.size(); x++)
  {
    finals[x].resize(U.cols());
    finals[x].setZero();
  }

  const int totalColumns = U.cols();
#pragma omp parallel
#pragma omp for  schedule(static)
  for (int z = 0; z < _zRes - 2; z++)
  {
    const int zCached = 3 * z * slabPeeled;
    //const int threadID = omp_get_thread_num();
    const int threadID = 0; 
    VectorXd& currentFinal = finals[threadID];

    int stride = 4;
    VectorXd unroll0(finals[threadID].size());
    VectorXd unroll1(finals[threadID].size());
    VectorXd unroll2(finals[threadID].size());
    VectorXd unroll3(finals[threadID].size());

    for (int y = 0; y < _yRes - 2; y++)
    {
      const int yzCached = 3 * y * xPeeled + zCached;

      int xEnd = (_xRes - 2) / stride * stride;

      unroll0.setZero();
      unroll1.setZero();
      unroll2.setZero();
      unroll3.setZero();
      for (int x = 0; x < xEnd; x += stride)
      {

        const int index0 = 3 * x + yzCached;
        const int index1 = 3 * (x + 1) + yzCached;
        const int index2 = 3 * (x + 2) + yzCached;
        const int index3 = 3 * (x + 3) + yzCached;

        const VEC3F v30 = (*this)(x + 1, y + 1, z + 1);
        const VEC3F v31 = (*this)(x + 2, y + 1, z + 1);
        const VEC3F v32 = (*this)(x + 3, y + 1, z + 1);
        const VEC3F v33 = (*this)(x + 4, y + 1, z + 1);

        const Vector3d v0(v30[0], v30[1], v30[2]);
        const Vector3d v1(v31[0], v31[1], v31[2]);
        const Vector3d v2(v32[0], v32[1], v32[2]);
        const Vector3d v3(v33[0], v33[1], v33[2]);

        unroll0.noalias() += U.block(index0, 0, 3, totalColumns).transpose() * v0;
        unroll1.noalias() += U.block(index1, 0, 3, totalColumns).transpose() * v1;
        unroll2.noalias() += U.block(index2, 0, 3, totalColumns).transpose() * v2;
        unroll3.noalias() += U.block(index3, 0, 3, totalColumns).transpose() * v3;

      }


      currentFinal.noalias() += unroll0;
      currentFinal.noalias() += unroll1;
      currentFinal.noalias() += unroll2;
      currentFinal.noalias() += unroll3;
      
      // cout << "current Final prior to unroll leftovers: " << endl << currentFinal << endl;

      // unrolling leftovers
      for (int x = xEnd; x < _xRes - 2; x++)
      {
        const int index = 3 * x + yzCached;
        const VEC3F v3 = (*this)(x + 1, y + 1, z + 1);
        const Vector3d v(v3[0], v3[1], v3[2]);
        currentFinal.noalias() += U.block(index, 0, 3, totalColumns).transpose() * v;
      }

      /*
      for (int x = 0; x < _xRes - 2; x++)
      {
        const int index = 3 * x + yzCached;
        const VEC3F v3 = (*this)(x + 1, y + 1, z + 1);
        const Vector3d v(v3[0], v3[1], v3[2]);
        currentFinal.noalias() += U.block(index, 0, 3, totalColumns).transpose() * v;
      }
      */
    }
  }

  TIMER finalReduce("peeledProject final reduce");

  // do a sum of the thread results
  VectorXd final(U.cols());
  final.setZero();
  for (unsigned int x = 0; x < finals.size(); x++)
    final += finals[x];
  finalReduce.stop();

  return final;
}

//////////////////////////////////////////////////////////////////////
// unproject the reduced coordinate into the peeled cells in this field
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::peeledUnproject(const MatrixXd& U, const VectorXd& q)
{
  TIMER functionTimer(__FUNCTION__);
  const int xPeeled = _xRes - 2;
  const int slabPeeled = (_xRes - 2) * (_yRes - 2);
  const int totalColumns = U.cols();

#pragma omp parallel
#pragma omp for  schedule(static)
  for (int z = 0; z < _zRes - 2; z++)
  {
    const int zCached = 3 * z * slabPeeled;
    for (int y = 0; y < _yRes - 2; y++)
    {
      const int yzCached = 3 * y * xPeeled + zCached;

      int stride = 2;
      int xEnd = (_xRes - 2) / stride * stride;

      // do the unrolling
      for (int x = 0; x < xEnd; x += stride)
      {
        const int index0 = 3 * x + yzCached;
        const int index1 = 3 * (x + 1) + yzCached;

        Real finalX0 = U.row(index0).dot(q);
        Real finalY0 = U.row(index0 + 1).dot(q);
        Real finalZ0 = U.row(index0 + 2).dot(q);
        Real finalX1 = U.row(index1).dot(q);
        Real finalY1 = U.row(index1 + 1).dot(q);
        Real finalZ1 = U.row(index1 + 2).dot(q);

        (*this)(x + 1, y + 1, z + 1) = VEC3F(finalX0, finalY0, finalZ0);
        (*this)(x + 2, y + 1, z + 1) = VEC3F(finalX1, finalY1, finalZ1);
      }

      // catch the leftovers
      for (int x = xEnd; x < _xRes - 2; x++)
      {
        const int index = 3 * x + yzCached;
        Real finalX = U.row(index).dot(q);
        Real finalY = U.row(index + 1).dot(q);
        Real finalZ = U.row(index + 2).dot(q);
        (*this)(x + 1, y + 1, z + 1) = VEC3F(finalX, finalY, finalZ);
      }

    }
  }
}
#if 0
{
  const int xPeeled = _xRes - 2;
  const int slabPeeled = (_xRes - 2) * (_yRes - 2);
  const int totalColumns = U.cols();

#pragma omp parallel
#pragma omp for  schedule(static)
  for (int z = 0; z < _zRes - 2; z++)
  {
    const int zCached = 3 * z * slabPeeled;
    for (int y = 0; y < _yRes - 2; y++)
    {
      const int yzCached = 3 * y * xPeeled + zCached;
      for (int x = 0; x < _xRes - 2; x++)
      {
        const int index = 3 * x + yzCached;

        /*
        VectorXd v = U.block(index, 0, 3, totalColumns) * q;
        (*this)(x + 1, y + 1, z + 1) = VEC3F(v[0], v[1], v[2]);
        */
        Real finalX = U.row(index).dot(q);
        Real finalY = U.row(index + 1).dot(q);
        Real finalZ = U.row(index + 2).dot(q);
        (*this)(x + 1, y + 1, z + 1) = VEC3F(finalX, finalY, finalZ);
      }
    }
  }
}
#endif

///////////////////////////////////////////////////////////////////////
// do a cubic Hermite interpolation
///////////////////////////////////////////////////////////////////////
Real VECTOR3_FIELD_3D::cubicInterp(const Real interp, const Real* points)
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
// tricubic interpolation lookup
///////////////////////////////////////////////////////////////////////
VEC3F VECTOR3_FIELD_3D::cubicLookup(const VEC3F& position) const
{
  VEC3F positionCopy = position;

  // get the lower corner position
  const VEC3F corner = _center - (Real)0.5 * _lengths + (Real)0.5 * dxs();

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= 1.0 / _dx;
  positionCopy[1] *= 1.0 / _dy;
  positionCopy[2] *= 1.0 / _dz;

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

  // do the interp for all three components
  VEC3F final;
  for (int x = 0; x < 3; x++)
  {
    // do the z0 slice
    const Real p0[] = {_data[x0 + y0z0][x], _data[x1 + y0z0][x], _data[x2 + y0z0][x], _data[x3 + y0z0][x]};
    const Real p1[] = {_data[x0 + y1z0][x], _data[x1 + y1z0][x], _data[x2 + y1z0][x], _data[x3 + y1z0][x]};
    const Real p2[] = {_data[x0 + y2z0][x], _data[x1 + y2z0][x], _data[x2 + y2z0][x], _data[x3 + y2z0][x]};
    const Real p3[] = {_data[x0 + y3z0][x], _data[x1 + y3z0][x], _data[x2 + y3z0][x], _data[x3 + y3z0][x]};

    // do the z1 slice
    const Real p4[] = {_data[x0 + y0z1][x], _data[x1 + y0z1][x], _data[x2 + y0z1][x], _data[x3 + y0z1][x]};
    const Real p5[] = {_data[x0 + y1z1][x], _data[x1 + y1z1][x], _data[x2 + y1z1][x], _data[x3 + y1z1][x]};
    const Real p6[] = {_data[x0 + y2z1][x], _data[x1 + y2z1][x], _data[x2 + y2z1][x], _data[x3 + y2z1][x]};
    const Real p7[] = {_data[x0 + y3z1][x], _data[x1 + y3z1][x], _data[x2 + y3z1][x], _data[x3 + y3z1][x]};

    // do the z2 slice
    const Real p8[] = {_data[x0 + y0z2][x], _data[x1 + y0z2][x], _data[x2 + y0z2][x], _data[x3 + y0z2][x]};
    const Real p9[] = {_data[x0 + y1z2][x], _data[x1 + y1z2][x], _data[x2 + y1z2][x], _data[x3 + y1z2][x]};
    const Real p10[] = {_data[x0 + y2z2][x], _data[x1 + y2z2][x], _data[x2 + y2z2][x], _data[x3 + y2z2][x]};
    const Real p11[] = {_data[x0 + y3z2][x], _data[x1 + y3z2][x], _data[x2 + y3z2][x], _data[x3 + y3z2][x]};

    // do the z3 slice
    const Real p12[] = {_data[x0 + y0z3][x], _data[x1 + y0z3][x], _data[x2 + y0z3][x], _data[x3 + y0z3][x]};
    const Real p13[] = {_data[x0 + y1z3][x], _data[x1 + y1z3][x], _data[x2 + y1z3][x], _data[x3 + y1z3][x]};
    const Real p14[] = {_data[x0 + y2z3][x], _data[x1 + y2z3][x], _data[x2 + y2z3][x], _data[x3 + y2z3][x]};
    const Real p15[] = {_data[x0 + y3z3][x], _data[x1 + y3z3][x], _data[x2 + y3z3][x], _data[x3 + y3z3][x]};
    
    const Real z0Points[] = {cubicInterp(xInterp, p0), cubicInterp(xInterp, p1), cubicInterp(xInterp, p2), cubicInterp(xInterp, p3)};
    const Real z1Points[] = {cubicInterp(xInterp, p4), cubicInterp(xInterp, p5), cubicInterp(xInterp, p6), cubicInterp(xInterp, p7)};
    const Real z2Points[] = {cubicInterp(xInterp, p8), cubicInterp(xInterp, p9), cubicInterp(xInterp, p10), cubicInterp(xInterp, p11)};
    const Real z3Points[] = {cubicInterp(xInterp, p12), cubicInterp(xInterp, p13), cubicInterp(xInterp, p14), cubicInterp(xInterp, p15)};

    const Real finalPoints[] = {cubicInterp(yInterp, z0Points), cubicInterp(yInterp, z1Points), cubicInterp(yInterp, z2Points), cubicInterp(yInterp, z3Points)};
    final[x] = cubicInterp(zInterp, finalPoints);
  }
  return final;
}

///////////////////////////////////////////////////////////////////////
// triquartic interpolation lookup
///////////////////////////////////////////////////////////////////////
VEC3F VECTOR3_FIELD_3D::quarticLookup(const VEC3F& position) const
{
  VEC3F positionCopy = position;

  // get the lower corner position
  const VEC3F corner = _center - (Real)0.5 * _lengths + (Real)0.5 * dxs();

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= 1.0 / _dx;
  positionCopy[1] *= 1.0 / _dy;
  positionCopy[2] *= 1.0 / _dz;

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

  // do the interp for all three components
  VEC3F final;
  for (int x = 0; x < 3; x++)
  {
    // do the z0 slice
    const Real p0[] = {_data[x0 + y0z0][x], _data[x1 + y0z0][x], _data[x2 + y0z0][x], _data[x3 + y0z0][x]};
    const Real p1[] = {_data[x0 + y1z0][x], _data[x1 + y1z0][x], _data[x2 + y1z0][x], _data[x3 + y1z0][x]};
    const Real p2[] = {_data[x0 + y2z0][x], _data[x1 + y2z0][x], _data[x2 + y2z0][x], _data[x3 + y2z0][x]};
    const Real p3[] = {_data[x0 + y3z0][x], _data[x1 + y3z0][x], _data[x2 + y3z0][x], _data[x3 + y3z0][x]};

    // do the z1 slice
    const Real p4[] = {_data[x0 + y0z1][x], _data[x1 + y0z1][x], _data[x2 + y0z1][x], _data[x3 + y0z1][x]};
    const Real p5[] = {_data[x0 + y1z1][x], _data[x1 + y1z1][x], _data[x2 + y1z1][x], _data[x3 + y1z1][x]};
    const Real p6[] = {_data[x0 + y2z1][x], _data[x1 + y2z1][x], _data[x2 + y2z1][x], _data[x3 + y2z1][x]};
    const Real p7[] = {_data[x0 + y3z1][x], _data[x1 + y3z1][x], _data[x2 + y3z1][x], _data[x3 + y3z1][x]};

    // do the z2 slice
    const Real p8[] = {_data[x0 + y0z2][x], _data[x1 + y0z2][x], _data[x2 + y0z2][x], _data[x3 + y0z2][x]};
    const Real p9[] = {_data[x0 + y1z2][x], _data[x1 + y1z2][x], _data[x2 + y1z2][x], _data[x3 + y1z2][x]};
    const Real p10[] = {_data[x0 + y2z2][x], _data[x1 + y2z2][x], _data[x2 + y2z2][x], _data[x3 + y2z2][x]};
    const Real p11[] = {_data[x0 + y3z2][x], _data[x1 + y3z2][x], _data[x2 + y3z2][x], _data[x3 + y3z2][x]};

    // do the z3 slice
    const Real p12[] = {_data[x0 + y0z3][x], _data[x1 + y0z3][x], _data[x2 + y0z3][x], _data[x3 + y0z3][x]};
    const Real p13[] = {_data[x0 + y1z3][x], _data[x1 + y1z3][x], _data[x2 + y1z3][x], _data[x3 + y1z3][x]};
    const Real p14[] = {_data[x0 + y2z3][x], _data[x1 + y2z3][x], _data[x2 + y2z3][x], _data[x3 + y2z3][x]};
    const Real p15[] = {_data[x0 + y3z3][x], _data[x1 + y3z3][x], _data[x2 + y3z3][x], _data[x3 + y3z3][x]};
    
    const Real z0Points[] = {quarticInterp(xInterp, p0), quarticInterp(xInterp, p1), quarticInterp(xInterp, p2), quarticInterp(xInterp, p3)};
    const Real z1Points[] = {quarticInterp(xInterp, p4), quarticInterp(xInterp, p5), quarticInterp(xInterp, p6), quarticInterp(xInterp, p7)};
    const Real z2Points[] = {quarticInterp(xInterp, p8), quarticInterp(xInterp, p9), quarticInterp(xInterp, p10), quarticInterp(xInterp, p11)};
    const Real z3Points[] = {quarticInterp(xInterp, p12), quarticInterp(xInterp, p13), quarticInterp(xInterp, p14), quarticInterp(xInterp, p15)};

    const Real finalPoints[] = {quarticInterp(yInterp, z0Points), quarticInterp(yInterp, z1Points), quarticInterp(yInterp, z2Points), quarticInterp(yInterp, z3Points)};
    final[x] = quarticInterp(zInterp, finalPoints);
  }
  return final;
  /*
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
  */
}

///////////////////////////////////////////////////////////////////////
// do a quartic WENO interpolation
///////////////////////////////////////////////////////////////////////
Real VECTOR3_FIELD_3D::quarticInterp(const Real interp, const Real* points)
{
  /*
  Real fim1 = points[0];
  Real fi   = points[1];
  Real fip1 = points[2];
  Real fip2 = points[3];
  Real x = interp;
  Real xSq = x * x;

  Real p1 = fi + (fip1 - fim1) * 0.5 * x +
                 (fip1 - 2 * fi - fim1) * 0.5 * xSq;

  Real p2 = fi + (-fip2 + 4.0 * fip1 - 3 * fi) * 0.5 * x +
                 (fip2 - 2.0 * fip1 + fi) * 0.5 * xSq;

  Real C1 = (2 - x) / 3.0;
  Real C2 = (x + 1) / 3.0;

  Real IS1 = (26 * fip1 * fim1 - 52 * fi * fim1 - 76 * fip1 * fi + 25 * fip1 * fip1 + 64 * fi * fi + 13 * fim1 * fim1) / 12.0;
  Real IS2 = (26 * fip2 * fi - 52 * fip2 * fip1 - 76 * fip1 * fi + 25 * fi * fi + 64 * fip1 * fip1 + 13 * fip2 * fip2) / 12.0;

  Real eps = 1e-6;
  Real alpha1 = C1 / ((eps + IS1) * (eps + IS1));
  Real alpha2 = C2 / ((eps + IS2) * (eps + IS2));

  Real w1 = alpha1 / (alpha1 + alpha2);
  Real w2 = alpha2 / (alpha1 + alpha2);

  return w1 * p1 + w2 * p2;
  */

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

  //return final;

  const Real intervalMax = (points[1] > points[2]) ? points[1] : points[2];
  const Real intervalMin = (points[1] > points[2]) ? points[2] : points[1];

  //return (final < intervalMax) ? ((final > intervalMin) ? final : intervalMin) : intervalMax;
  
  const Real linear = x * points[2] + (1.0 - x) * points[1];
  return (final < intervalMax) ? ((final > intervalMin) ? final : linear) : linear;
}

///////////////////////////////////////////////////////////////////////
// take the curl to get a divergence free velocity field
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::curl(const FIELD_3D& scalar1, const FIELD_3D& scalar2, const FIELD_3D& scalar3)
{
  VECTOR3_FIELD_3D final(scalar1);
  int zRes = scalar1.zRes();
  int yRes = scalar1.yRes();
  int xRes = scalar1.xRes();

  for (int z = 0; z < zRes; z ++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        Real s3Dy = scalar3.Dy(x,y,z);
        Real s2Dz = scalar2.Dz(x,y,z);

        Real s1Dz = scalar1.Dz(x,y,z);
        Real s3Dx = scalar3.Dx(x,y,z);
        
        Real s2Dx = scalar2.Dx(x,y,z);
        Real s1Dy = scalar1.Dy(x,y,z);

        final(x,y,z) = VEC3F(s3Dy - s2Dz,
                             s1Dz - s3Dx,
                             s2Dx - s1Dy);
      }

  return final;
}

///////////////////////////////////////////////////////////////////////
// take the curl to get a divergence free velocity field
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::curl(const VECTOR3_FIELD_3D& potential)
{
  VECTOR3_FIELD_3D final(potential);
  final.clear();
  int zRes = potential.zRes();
  int yRes = potential.yRes();
  int xRes = potential.xRes();

  FIELD_3D scalar1 = potential.scalarField(0);
  FIELD_3D scalar2 = potential.scalarField(1);
  FIELD_3D scalar3 = potential.scalarField(2);

  for (int z = 0; z < zRes; z ++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        Real s3Dy = scalar3.Dy(x,y,z);
        Real s2Dz = scalar2.Dz(x,y,z);

        Real s1Dz = scalar1.Dz(x,y,z);
        Real s3Dx = scalar3.Dx(x,y,z);
        
        Real s2Dx = scalar2.Dx(x,y,z);
        Real s1Dy = scalar1.Dy(x,y,z);

        final(x,y,z) = VEC3F(s3Dy - s2Dz,
                             s1Dz - s3Dx,
                             s2Dx - s1Dy);
      }

  return final;
}

//////////////////////////////////////////////////////////////////////
// return a turbulence field
//////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::noisePotential(const FIELD_3D& example)
{
  int xRes = example.xRes(); 
  int yRes = example.yRes(); 
  int zRes = example.zRes(); 

  // load up the noise tile, or if none exists, create it
  int gridSize = noiseTileSize * noiseTileSize * noiseTileSize;

  Real* noiseTile[3];
  for (int seed = 0; seed < 3; seed++)
  {
    noiseTile[seed] = new Real[gridSize];
    string filename("./data/wavelet.noise.seed.");
    char buffer[256];
    sprintf(buffer, "%i", seed);
    filename = filename + string(buffer);
    bool success =  loadTile(noiseTile[seed], filename);
    if (!success)
      generateTile(noiseTile[seed], filename, seed);
  }

  VECTOR3_FIELD_3D final(example);
  final.clear();
  int octaves = 1;

  int biggest = (xRes > yRes) ? xRes : yRes;
  biggest = (biggest > zRes) ? biggest : zRes;

  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        float xFrac = (float)x / biggest;
        float yFrac = (float)y / biggest;
        float zFrac = (float)z / biggest;
        VEC3F point(xFrac, yFrac, zFrac);

        for (int seed = 0; seed < 3; seed++)
        {
          float accum = 0.0;
          for (int i = 0; i < octaves; i++)
          {
            //float scale = 100.0 * pow(2.0, i);
            float scale = 10.0 * pow(2.0, i);

            //if (scale > biggest) continue;
            float coefficient = pow(0.5f, i + 1);

            VEC3F point(scale * xFrac, scale * yFrac, scale * zFrac);
            accum += coefficient * WNoise(point, noiseTile[seed]);
          } 
          final(x,y,z)[seed] += accum;
        }
      }
  //final.normalize();

  for (int x = 0; x < 3; x++)
    delete[] noiseTile[x];

  return final;
}

//////////////////////////////////////////////////////////////////////
// take the dot product of the current field with another vector field
// and return the scalar field
//////////////////////////////////////////////////////////////////////
FIELD_3D VECTOR3_FIELD_3D::dot(const VECTOR3_FIELD_3D& rhs)
{
  FIELD_3D final(_xRes, _yRes, _zRes);

  assert(rhs.xRes() == _xRes);
  assert(rhs.yRes() == _yRes);
  assert(rhs.zRes() == _zRes);

  //for (int x = 0; x < _totalCells; x++)
  //  final[x] = rhs[x] * (*this)[x];
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        final(x,y,z) = rhs(x,y,z) * (*this)(x,y,z);
        if (x == 0 && y == 0 && z == 0)
        {
          cout << "RHS: " << rhs(x,y,z) << endl;
          cout << "this: " << (*this)(x,y,z) << endl;
          cout << "final: " << final(x,y,z) << endl;
        }

      }

  return final;
}

//////////////////////////////////////////////////////////////////////
// get a z slice
//////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_2D VECTOR3_FIELD_3D::zSlice(const int z) const
{
  VECTOR3_FIELD_2D final(_xRes, _yRes);

  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
      final(x,y) = (*this)(x,y,z);

  return final;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::normalizeToLargest()
{
  Real maxMagnitude = 0.0;
  for (int x = 0; x < _totalCells; x++)
  {
    Real candidate = norm(_data[x]);
    if (candidate > maxMagnitude)
      maxMagnitude = candidate;
  }

  Real inverse = 1.0 / maxMagnitude;
  for (int x = 0; x < _totalCells; x++)
    _data[x] *= inverse;
}
