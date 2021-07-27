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
#include "FIELD_2D.h"
#ifndef WIN32
#include <jpeglib.h>
#endif
#include <png.h>
//#include <WAVELET_NOISE.h>
#include "TIMER.h"

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D::FIELD_2D(const int& rows, const int& cols, const VEC3F& center, const VEC3F& lengths) :
  _xRes(rows), _yRes(cols), _center(center), _lengths(lengths)
{
  _totalCells = _xRes * _yRes;
  _data = new Real[_totalCells];
  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;

  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;

  for (int x = 0; x < _totalCells; x++)
    _data[x] = 0.0;
}

FIELD_2D::FIELD_2D(const FIELD_2D& m) :
  _xRes(m.xRes()), _yRes(m.yRes()), _center(m.center()), _lengths(m.lengths()),
  _dx(m.dx()), _dy(m.dy()), _invDx(m.invDx()), _invDy(m.invDy())
{
  _totalCells = _xRes * _yRes;
  _data = new Real[_totalCells];

  for (int x = 0; x < _totalCells; x++)
    _data[x] = m[x];
}

FIELD_2D::FIELD_2D() :
  _xRes(0), _yRes(0), _totalCells(0), _data(NULL), _dx(0), _dy(0),
  _invDx(0), _invDy(0)
{
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D::~FIELD_2D()
{
  delete[] _data;
}
  
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_2D::clear()
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] = 0.0;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_2D::write(string filename) const
{
  FILE* file;
  file = fopen(filename.c_str(), "wb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " FIELD_2D write failed! " << endl;
    cout << " Could not open file " << filename.c_str() << endl;
    exit(0);
  }

  // write dimensions
  fwrite((void*)&_xRes, sizeof(int), 1, file);
  fwrite((void*)&_yRes, sizeof(int), 1, file);

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
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_2D::write(FILE* file) const
{
  // write dimensions
  fwrite((void*)&_xRes, sizeof(int), 1, file);
  fwrite((void*)&_yRes, sizeof(int), 1, file);

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
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_2D::read(string filename)
{
  FILE* file;
  file = fopen(filename.c_str(), "rb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " FIELD_2D read failed! " << endl;
    cout << " Could not open file " << filename.c_str() << endl;
    exit(0);
  }

  // read dimensions
  fread((void*)&_xRes, sizeof(int), 1, file);
  fread((void*)&_yRes, sizeof(int), 1, file);
  _totalCells = _xRes * _yRes;
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

  fclose(file);

  _lengths = VEC3F(1,1,0);

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_2D::read(FILE* file)
{
  // read dimensions
  fread((void*)&_xRes, sizeof(int), 1, file);
  fread((void*)&_yRes, sizeof(int), 1, file);
  _totalCells = _xRes * _yRes;
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

  _lengths = VEC3F(1,1,0);

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
}

///////////////////////////////////////////////////////////////////////
// jpeglib code based on:
// http://andrewewhite.net/wordpress/2008/09/02/very-simple-jpeg-writer-in-c-c
///////////////////////////////////////////////////////////////////////
void FIELD_2D::writeJPG(string filename)
{
#ifndef WIN32
  FILE* file = fopen(filename.c_str(), "wb");
 
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Could not open file: " << filename.c_str() << endl;
    exit(0);
  }

  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr       jerr;
   
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo, file);
   
  cinfo.image_width      = _xRes;
  cinfo.image_height     = _yRes;
  cinfo.input_components = 3;
  cinfo.in_color_space   = JCS_RGB;

  jpeg_set_defaults(&cinfo);
  /*set the quality [0..100]  */
  jpeg_set_quality (&cinfo, 75, true);
  jpeg_start_compress(&cinfo, true);

  // copy data to a char buffer
  unsigned char* buffer = new unsigned char[3 * _totalCells];
  for (int x = 0; x < _totalCells; x++)
  {
    Real entry = _data[x];
    entry = (entry < 0.0) ? 0.0 : entry;
    entry = (entry > 1.0) ? 1.0 : entry;

    buffer[3 * x] = (unsigned char) (255 * entry);
    buffer[3 * x + 1] = (unsigned char) (255 * entry);
    buffer[3 * x + 2] = (unsigned char) (255 * entry);
  }

  JSAMPROW row_pointer;
 
  while (cinfo.next_scanline < cinfo.image_height) {
    int index = cinfo.next_scanline * 3 * _xRes;
    row_pointer = (JSAMPROW)&buffer[index];
    jpeg_write_scanlines(&cinfo, &row_pointer, 1);
  }
  jpeg_finish_compress(&cinfo);

  delete[] buffer;

  fclose(file);
#endif
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_2D::writePPM(string filename)
{
  FILE *fp;
  unsigned char* pixels = new unsigned char[3 * _totalCells];

  for (int x = 0; x < _totalCells; x++)
  {
    pixels[3 * x] = 255 * _data[x];
    pixels[3 * x + 1] = 255 * _data[x];
    pixels[3 * x + 2] = 255 * _data[x];
  }

  fp = fopen(filename.c_str(), "wb");
  fprintf(fp, "P6\n%d %d\n255\n", _xRes, _yRes);
  fwrite(pixels, 1, _totalCells * 3, fp);
  fclose(fp);
  delete[] pixels;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_2D::writeMatlab(string filename, string variableName) const
{
  FILE* file;
  file = fopen(filename.c_str(), "w");
  fprintf(file, "%s = [", variableName.c_str());
  for (int y = 0; y < _yRes; y++)
  {
    for (int x = 0; x < _xRes; x++)
      fprintf(file, "%f ", (double)(*this)(x,y));
    fprintf(file, "; ");
  }
  fprintf(file, "];\n");

  fclose(file);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_2D::shiftFFT()
{
  assert(_xRes % 2 == 1);
  assert(_yRes % 2 == 1);
  Real* scratch = new Real[_xRes * _yRes];

  int xHalf = _xRes / 2;
  int yHalf = _yRes / 2;

  int xMod = _xRes % 2;
  int yMod = _yRes % 2;

  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < xHalf + xMod; x++)
    {
      int index = x + y * _xRes;
      scratch[index] = _data[index + xHalf + xMod];
    }
  for (int y = 0; y < _yRes; y++)
    for (int x = xHalf; x < _xRes; x++)
    {
      int index = x + y * _xRes;
      scratch[index] = _data[index - xHalf];
    }

  for (int y = 0; y < yHalf + yMod; y++)
    for (int x = 0; x < _xRes; x++)
    {
      int original = x + y * _xRes;
      int copy = x + (y + yHalf) * _xRes;
      _data[copy] = scratch[original];
    }

  for (int y = yHalf; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
    {
      int original = x + y * _xRes;
      int copy = x + (y - yHalf - yMod) * _xRes;
      _data[copy] = scratch[original];
    }

  delete[] scratch;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_2D::normalize()
{
  Real maxFound = 0.0;
  Real minFound = _data[0];
  for (int x = 0; x < _totalCells; x++)
  {
    maxFound = (_data[x] > maxFound) ? _data[x] : maxFound;
    minFound = (_data[x] < minFound) ? _data[x] : minFound;
  }

  float range = 1.0 / (maxFound - minFound);
  for (int x = 0; x < _totalCells; x++)
    _data[x] = (_data[x] - minFound) * range;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D& FIELD_2D::abs()
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] = fabs(_data[x]);

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_2D::resizeAndWipe(int xRes, int yRes, const VEC3F& center, const VEC3F& lengths)
{
  if (_xRes == xRes && _yRes == yRes)
  {
    clear();
    return;
  }

  if (_data)
    delete[] _data;

  _xRes = xRes;
  _yRes = yRes;
  _totalCells = _xRes * _yRes;

  _data = new Real[_xRes * _yRes];
 
  _center = center;
  _lengths = lengths; 
  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
}

///////////////////////////////////////////////////////////////////////
// set this field to the result of convolving filter and input
///////////////////////////////////////////////////////////////////////
void FIELD_2D::convolve(const FIELD_2D& filter, const FIELD_2D& input)
{
  TIMER functionTimer(__FUNCTION__);

  int filterWidth = filter.xRes();

  assert(filterWidth == filter.yRes());
  assert(filterWidth % 2 == 1);
  assert(input.xRes() == _xRes);
  assert(input.yRes() == _yRes);

  int iwidth = input.xRes();
  int iheight = input.yRes();

  // periodic boundaries
  int halfWidth = filterWidth / 2;
  for (int ix = 0; ix < iwidth; ix++)
    for (int iy = 0; iy < iheight; iy++)
    {
      int index = ix + iwidth * iy;
      float vd = 0;
      for (int iix = -halfWidth; iix <= halfWidth; iix++)
        for (int iiy = -halfWidth; iiy <= halfWidth; iiy++)
        {
          int xIndex = ix + iix;
          int yIndex = iy + iiy;

          if (xIndex < 0)
            xIndex += iwidth;
          if (yIndex < 0)
            yIndex += iheight;

          if (xIndex >= iwidth)
            xIndex -= iwidth;
          if (yIndex >= iheight)
            yIndex -= iheight;
          
          vd += filter(iix + halfWidth,iiy + halfWidth) * input(xIndex, yIndex);
        }
      _data[index] = vd;
    }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D& FIELD_2D::operator=(const Real& alpha)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] = alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D& FIELD_2D::operator*=(const Real& alpha)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] *= alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D& FIELD_2D::operator/=(const Real& alpha)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] /= alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D& FIELD_2D::operator+=(const Real& alpha)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] += alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D& FIELD_2D::operator-=(const Real& alpha)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] -= alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D& FIELD_2D::operator-=(const FIELD_2D& input)
{
  assert(input.xRes() == _xRes);
  assert(input.yRes() == _yRes);
  for (int x = 0; x < _totalCells; x++)
    _data[x] -= input[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D& FIELD_2D::operator+=(const FIELD_2D& input)
{
  assert(input.xRes() == _xRes);
  assert(input.yRes() == _yRes);
  for (int x = 0; x < _totalCells; x++)
    _data[x] += input[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D& FIELD_2D::operator*=(const FIELD_2D& input)
{
  assert(input.xRes() == _xRes);
  assert(input.yRes() == _yRes);

  for (int x = 0; x < _totalCells; x++)
    _data[x] *= input[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D& FIELD_2D::operator/=(const FIELD_2D& input)
{
  assert(input.xRes() == _xRes);
  assert(input.yRes() == _yRes);

  for (int x = 0; x < _totalCells; x++)
    if (fabs(input[x]) > 1e-6)
      _data[x] /= input[x];
    else
      _data[x] = 0;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D operator*(const FIELD_2D& A, const Real alpha)
{
  FIELD_2D final(A);
  final *= alpha;
  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D operator/(const FIELD_2D& A, const Real alpha)
{
  FIELD_2D final(A);
  final /= alpha;
  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D operator+(const FIELD_2D& A, const FIELD_2D& B)
{
  FIELD_2D final(A);
  final += B;
  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D operator-(const FIELD_2D& A, const FIELD_2D& B)
{
  FIELD_2D final(A);
  final -= B;
  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D operator*(const FIELD_2D& A, const FIELD_2D& B)
{
  FIELD_2D final(A);
  final *= B;
  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D operator+(const FIELD_2D& A, const Real alpha)
{
  FIELD_2D final(A);
  final += alpha;
  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D operator*(const Real alpha, const FIELD_2D& A)
{
  return A * alpha;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D operator+(const Real alpha, const FIELD_2D& A)
{
  return A + alpha;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D& FIELD_2D::operator=(const FIELD_2D& A)
{
  resizeAndWipe(A.xRes(), A.yRes(), A.center(), A.lengths());

  for (int x = 0; x < _totalCells; x++)
    _data[x] = A[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
// sum of all entries
///////////////////////////////////////////////////////////////////////
Real FIELD_2D::sum()
{
  Real total = 0;
  for (int x = 0; x < _totalCells; x++)
    total += _data[x];

  return total;
}

///////////////////////////////////////////////////////////////////////
// take the log
///////////////////////////////////////////////////////////////////////
void FIELD_2D::log(Real base)
{
#if QUAD_PRECISION
  Real scale = 1.0 / logq(base);
  for (int x = 0; x < _totalCells; x++)
    _data[x] = logq(_data[x]) * scale;
#else
  Real scale = 1.0 / std::log(base);
  for (int x = 0; x < _totalCells; x++)
    _data[x] = std::log(_data[x]) * scale;
#endif
}

///////////////////////////////////////////////////////////////////////
// Compute the elements of the vertical derivative convolution kernel
///////////////////////////////////////////////////////////////////////
void FIELD_2D::verticalDerivativeKernel(double kMax, double dk, double sigma, double L)
{
  assert(_xRes % 2);
  assert(_xRes == _yRes);

  double norm = 0;

  for (double k = 0; k < kMax; k += dk)
    norm += k * k * exp(-sigma * k * k);

  int halfWidth = _xRes / 2;

  for (int i = -halfWidth; i <= halfWidth; i++)
    for (int j = -halfWidth; j <= halfWidth; j++)
    {
      double r = sqrt((float)(i * i + j * j));
      double kern = 0;
      for (double k = 0; k < kMax; k += dk)
        kern += k * k * (sqrt(k * k * L * L)) * exp(-sigma * k * k) * j0(r * k);

      (*this)(i + halfWidth, j + halfWidth) = kern / norm;
    }
}

///////////////////////////////////////////////////////////////////////
// Compute a radial Bessel function
///////////////////////////////////////////////////////////////////////
void FIELD_2D::radialBessel()
{
  assert(_xRes % 2);
  assert(_xRes == _yRes);

  int halfWidth = _xRes / 2;
  //int kMax = 2;
  int kMin = 0;
  int kMax = 5;
  double dk = 0.01;
  //double dk = 1;

  for (int i = -halfWidth; i <= halfWidth; i++)
    for (int j = -halfWidth; j <= halfWidth; j++)
    {
      double r = sqrt((float)(i * i + j * j)) / (float)_xRes * 20;
      double kern = 0.0;
      for (double k = kMin; k < kMax; k += dk)
        kern += j0(r * k);

      (*this)(i + halfWidth, j + halfWidth) = kern;
    }
}

///////////////////////////////////////////////////////////////////////
// set to a bessel function
///////////////////////////////////////////////////////////////////////
void FIELD_2D::setToBessel(float k)
{
  assert(_xRes % 2);
  assert(_xRes == _yRes);

  int halfWidth = _xRes / 2;

  for (int i = -halfWidth; i <= halfWidth; i++)
    for (int j = -halfWidth; j <= halfWidth; j++)
    {
      double r = sqrt((float)(i * i + j * j)) / (float)_xRes * 20;
      double kern = 0.0;
      kern += j0(r * k);

      (*this)(i + halfWidth, j + halfWidth) = kern;
    }
}

///////////////////////////////////////////////////////////////////////
// upsample the texture by a certain factor
///////////////////////////////////////////////////////////////////////
FIELD_2D FIELD_2D::nearestNeighborUpsample(int factor)
{
  FIELD_2D final(_xRes * factor, _yRes * factor);

  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
    {
      double value = (*this)(x,y);
      for (int i = 0; i < factor; i++)
        for (int j = 0; j < factor; j++)
          final(factor * x + i, factor * y + j) = value;
    }
  return final;
}

///////////////////////////////////////////////////////////////////////
// get the min of the field
///////////////////////////////////////////////////////////////////////
Real FIELD_2D::min()
{
  assert(_xRes > 0);
  assert(_yRes > 0);
  Real final = _data[0];

  for (int i = 0; i < _xRes * _yRes; i++)
    final = (_data[i] < final) ? _data[i] : final;

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the max of the field
///////////////////////////////////////////////////////////////////////
Real FIELD_2D::max()
{
  assert(_xRes > 0);
  assert(_yRes > 0);
  Real final = _data[0];

  for (int i = 0; i < _xRes * _yRes; i++)
    final = (_data[i] > final) ? _data[i] : final;

  return final;
}

///////////////////////////////////////////////////////////////////////
// set to a checkboard for debugging
///////////////////////////////////////////////////////////////////////
void FIELD_2D::setToCheckerboard(int xChecks, int yChecks)
{
  int index = 0;
  for (int x = 0; x < _xRes; x++)
    for (int y = 0; y < _yRes; y++, index++)
    {
      int xMod = (x / (_xRes / xChecks)) % 2;
      int yMod = (y / (_yRes / yChecks)) % 2;

      if ((xMod && yMod) || (!xMod && !yMod))
        _data[index] = 1;
    }
}

///////////////////////////////////////////////////////////////////////
// set to a checkboard for debugging
///////////////////////////////////////////////////////////////////////
void FIELD_2D::setToRampedCheckerboard(int xChecks, int yChecks)
{
  int index = 0;
  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++, index++)
    {
      int xMod = (x / (_xRes / xChecks)) % 2;
      int yMod = (y / (_yRes / yChecks)) % 2;

      if ((xMod && yMod) || (!xMod && !yMod))
        _data[index] = (float)y / _yRes;
    }
}

///////////////////////////////////////////////////////////////////////
// pass a field to fieldViewer2D
///////////////////////////////////////////////////////////////////////
void FIELD_2D::fieldViewer(const FIELD_2D& field, string name)
{
  field.write("temp.field");
  string execute("./bin/fieldViewer temp.field \"");
  execute = execute + name + string("\" &");
  cout << " Executing " << execute.c_str() << endl;
  system(execute.c_str());
}

///////////////////////////////////////////////////////////////////////
// pass a field to fieldViewer2D
///////////////////////////////////////////////////////////////////////
void FIELD_2D::overlayFieldViewer(const FIELD_2D& distanceField, const FIELD_2D& field, string name)
{
  distanceField.write("temp.distance.field");
  field.write("temp.field");
  string execute("./bin/overlayFieldViewer2D temp.distance.field temp.field \"");
  execute = execute + name + string("\" &");
  cout << " Executing " << execute.c_str() << endl;
  system(execute.c_str());
}

///////////////////////////////////////////////////////////////////////
// set to a ramp for debugging
///////////////////////////////////////////////////////////////////////
void FIELD_2D::setToRampX()
{
  int index = 0;
  for (int x = 0; x < _xRes; x++)
    for (int y = 0; y < _yRes; y++, index++)
      _data[index] = (float)x / _xRes;
}

///////////////////////////////////////////////////////////////////////
// set to a ramp for debugging
///////////////////////////////////////////////////////////////////////
void FIELD_2D::setToRampY()
{
  int index = 0;
  for (int x = 0; x < _xRes; x++)
    for (int y = 0; y < _yRes; y++, index++)
      _data[index] = (float)y / _yRes;
}

///////////////////////////////////////////////////////////////////////
// return a field for the Laplacian of this field
///////////////////////////////////////////////////////////////////////
FIELD_2D FIELD_2D::laplacian()
{
  FIELD_2D final(_xRes, _yRes);

  // assume a unit field
  Real dx = 1.0 / _xRes;
  Real dy = 1.0 / _yRes;
  Real dx2 = dx * dx;
  Real dy2 = dy * dy;

  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
    {
      Real yPlus  = safe(x, y+1); 
      Real yMinus = safe(x, y-1); 
      Real xPlus  = safe(x+1, y); 
      Real xMinus = safe(x-1, y);

      final(x,y) = (-2.0 * (*this)(x,y) + xPlus + xMinus) / dx2 + (-2.0 * (*this)(x,y) + yPlus + yMinus) / dy2; 
    }
  return final;
}

///////////////////////////////////////////////////////////////////////
// return a field for the Laplacian of this field
///////////////////////////////////////////////////////////////////////
FIELD_2D FIELD_2D::laplacian4th()
{
  FIELD_2D final(_xRes, _yRes);

  // assume a unit field
  Real dx = 1.0 / _xRes;
  Real dy = 1.0 / _yRes;
  Real dx2 = dx * dx;
  Real dy2 = dy * dy;

  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
    {
      Real yPlus  = safe(x, y+1); 
      Real yPlusPlus  = safe(x, y+2); 
      Real yMinus = safe(x, y-1); 
      Real yMinusMinus = safe(x, y-2); 
      Real xPlus  = safe(x+1, y); 
      Real xPlusPlus  = safe(x+2, y); 
      Real xMinus = safe(x-1, y);
      Real xMinusMinus = safe(x-2, y);

      final(x,y) = ((-1.0 / 12.0) * xMinusMinus + (4.0 / 3.0) * xMinus + (-5.0 / 2.0) * (*this)(x,y) + (4.0 / 3.0) * xPlus + (-1.0 / 12.0) * xPlusPlus) / dx2;
      final(x,y) += ((-1.0 / 12.0) * yMinusMinus + (4.0 / 3.0) * yMinus + (-5.0 / 2.0) * (*this)(x,y) + (4.0 / 3.0) * yPlus + (-1.0 / 12.0) * yPlusPlus) / dy2;
    }
  return final;
}

///////////////////////////////////////////////////////////////////////
// a safe, toroidal data accessor -- does all bounds checking for you
///////////////////////////////////////////////////////////////////////
Real FIELD_2D::safe(int x, int y)
{
  while (x < 0)
    x += _xRes;
  while (y < 0)
    y += _yRes;

  while (x > _xRes - 1)
    x -= _xRes;
  while (y > _yRes - 1)
    y -= _yRes;

  return (*this)(x,y);
}

///////////////////////////////////////////////////////////////////////
// return a field for the gradient of this field
///////////////////////////////////////////////////////////////////////
FIELD_2D FIELD_2D::gradient()
{
  FIELD_2D final(_xRes, _yRes);

  // assume a unit field
  Real dx = 1.0 / _xRes;
  Real dy = 1.0 / _yRes;

  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
    {
      Real yPlus  = safe(x, y+1); 
      Real yMinus = safe(x, y-1); 
      Real xPlus  = safe(x+1, y); 
      Real xMinus = safe(x-1, y);

      final(x,y) = (xPlus - xMinus) / dx + (yPlus - yMinus) / dy;
    }
  return final;
}

///////////////////////////////////////////////////////////////////////
// return the transpose (flip x and y)
///////////////////////////////////////////////////////////////////////
FIELD_2D FIELD_2D::transpose() const
{
  FIELD_2D final(_yRes, _xRes);

  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
      final(y,x) = (*this)(x,y);

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VEC3F FIELD_2D::maxIndex()
{
  Real maxFound = _data[0];

  VEC3F maxFoundIndex;
  int index = 0;
  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++, index++)
      if (_data[index] > maxFound)
      {
        maxFound = _data[index];

        maxFoundIndex[0] = x;
        maxFoundIndex[1] = y;
      }

  return maxFoundIndex;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VEC3F FIELD_2D::minIndex()
{
  Real minFound = _data[0];

  VEC3F minFoundIndex;
  int index = 0;
  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++, index++)
      if (_data[index] < minFound)
      {
        minFound = _data[index];

        minFoundIndex[0] = x;
        minFoundIndex[1] = y;
      }

  return minFoundIndex;
}

///////////////////////////////////////////////////////////////////////
// get the real-space location of a cell
///////////////////////////////////////////////////////////////////////
VEC3F FIELD_2D::cellCenter(int x, int y) const
{
  VEC3F halfLengths = (Real)0.5 * _lengths;

  // set it to the lower corner
  VEC3F final = _center - halfLengths;

  // displace to the NNN corner
  final[0] += x * _dx;
  final[1] += y * _dy;

  // displace it to the cell center
  final[0] += _dx * 0.5;
  final[1] += _dy * 0.5;

  final[2] = 0;

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the real-space location of a cell
///////////////////////////////////////////////////////////////////////
VEC3F FIELD_2D::cellCenter(int index) const
{
  int x = index % _xRes;
  int y = index / _xRes;
  return cellCenter(x,y);
}

///////////////////////////////////////////////////////////////////////
// get the index of the nearest cell
///////////////////////////////////////////////////////////////////////
int FIELD_2D::cellIndex(const VEC3F& position) const
{
  VEC3F positionCopy = position;

  // get the lower corner position
  VEC3F corner = _center - (Real)0.5 * _lengths;
  VEC3F dxs(_dx, _dy, 0);
  corner += (Real)0.5 * dxs;

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  
  int x0 = (int)positionCopy[0];
  int x1    = x0 + 1;
  int y0 = (int)positionCopy[1];
  int y1    = y0 + 1;

  // clamp everything
  x0 = (x0 < 0) ? 0 : x0;
  y0 = (y0 < 0) ? 0 : y0;
  
  x1 = (x1 < 0) ? 0 : x1;
  y1 = (y1 < 0) ? 0 : y1;

  x0 = (x0 > _xRes - 1) ? _xRes - 1 : x0;
  y0 = (y0 > _yRes - 1) ? _yRes - 1 : y0;

  x1 = (x1 > _xRes - 1) ? _xRes - 1 : x1;
  y1 = (y1 > _yRes - 1) ? _yRes - 1 : y1;

  // get interpolation weights
  const float s1 = positionCopy[0]- x0;
  const float s0 = 1.0f - s1;
  const float t1 = positionCopy[1]- y0;
  const float t0 = 1.0f - t1;

  int xFinal = (s0 > s1) ? x0 : x1;
  int yFinal = (t0 > t1) ? y0 : y1;

  return xFinal + yFinal * _xRes;
}

///////////////////////////////////////////////////////////////////////
// get the index of the nearest cell
///////////////////////////////////////////////////////////////////////
void FIELD_2D::cellIndex(const VEC3F& position, int& x, int& y) const
{
  int index = cellIndex(position);

  x = index % _xRes;
  y = index / _xRes;
}

///////////////////////////////////////////////////////////////////////
// get the sum of the absolute value of all the entries
///////////////////////////////////////////////////////////////////////
Real FIELD_2D::fabsSum()
{
  Real final = 0;

  for (int x = 0; x < _totalCells; x++)
    final += fabs(_data[x]);

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the mean curvature
///////////////////////////////////////////////////////////////////////
FIELD_2D FIELD_2D::meanCurvature() const
{
  FIELD_2D final(_xRes, _yRes);

  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
    {
      Real px = Dx(x,y);
      Real py = Dy(x,y);
      
      Real pxx = DDx(x,y);
      Real pyy = DDy(x,y);
      
      Real pxy = DDxy(x,y);

      final(x,y) = (pyy) * px * px + (pxx) * py * py -
                   2.0 * (px * py * pxy);
    }

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real FIELD_2D::Dx(int x, int y) const
{
  assert(x >= 0);
  assert(x < _xRes);
  assert(y >= 0);
  assert(y < _yRes);

  int index = x + y * _xRes;

  const Real right = (x < _xRes - 1) ? _data[index + 1] : _data[index];
  const Real left  = (x > 0)         ? _data[index - 1] : _data[index];
  const Real denom = (x > 0 && x < _xRes -1) ? 1.0 / (2.0 * _dx) : 1.0 / _dx;
  return (right - left) * denom;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real FIELD_2D::DDx(int x, int y) const
{
  assert(x >= 0);
  assert(x < _xRes);
  assert(y >= 0);
  assert(y < _yRes);

  int index = x + y * _xRes;

  const Real right = (x < _xRes - 1) ? _data[index + 1] : _data[index];
  const Real left  = (x > 0)         ? _data[index - 1] : _data[index];
  const Real denom = 1.0 / (_dx * _dx);
  return (right - 2.0 * _data[index] + left) * denom;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real FIELD_2D::Dy(int x, int y) const
{
  assert(x >= 0);
  assert(x < _xRes);
  assert(y >= 0);
  assert(y < _yRes);

  int index = x + y * _xRes;

  const Real up   = (y < _yRes - 1) ? _data[index + _xRes] : _data[index];
  const Real down = (y > 0)         ? _data[index - _xRes] : _data[index];
  const Real denom = (y > 0 && y < _yRes -1) ? 1.0 / (2.0 * _dy) : 1.0 / _dy;
  return (up - down) * denom;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real FIELD_2D::DDy(int x, int y) const
{
  assert(x >= 0);
  assert(x < _xRes);
  assert(y >= 0);
  assert(y < _yRes);

  int index = x + y * _xRes;

  const Real up   = (y < _yRes - 1) ? _data[index + _xRes] : _data[index];
  const Real down = (y > 0)         ? _data[index - _xRes] : _data[index];
  const Real denom = 1.0 / (_dy * _dy);
  return (up - 2.0 * _data[index] + down) * denom;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real FIELD_2D::DDxy(int x, int y) const
{
  assert(x >= 0);
  assert(x < _xRes);
  assert(y >= 0);
  assert(y < _yRes);

  int index = x + y * _xRes;

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
// get the Gaussian curvature
///////////////////////////////////////////////////////////////////////
FIELD_2D FIELD_2D::gaussianCurvature() const
{
  FIELD_2D final(_xRes, _yRes);

  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
    {
      Real px = Dx(x,y);
      Real py = Dy(x,y);
      
      Real pxx = DDx(x,y);
      Real pyy = DDy(x,y);
      
      //Real pxy = DDxy(x,y);

      final(x,y) = (pyy) * px * px +
                   (pxx) * py * py; 
    }

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the normal at a point
///////////////////////////////////////////////////////////////////////
VEC3F FIELD_2D::normal(int x, int y) const
{
  assert (x >= 0 && x < _xRes);
  assert (y >= 0 && y < _yRes);

  VEC3F final;
  final[0] = Dx(x,y);
  final[1] = Dy(x,y);
  final[2] = 0;

  final.normalize();

  return final;
}

//////////////////////////////////////////////////////////////////////
// count the number of infs in the field
//////////////////////////////////////////////////////////////////////
int FIELD_2D::totalInfs() const
{
  int final = 0;
  for (int x = 0; x < _totalCells; x++)
    if (std::isinf(_data[x])) final++;

  return final;
}

//////////////////////////////////////////////////////////////////////
// count the number of NaNs in the field
//////////////////////////////////////////////////////////////////////
int FIELD_2D::totalNans() const
{
  int final = 0;
  for (int x = 0; x < _totalCells; x++)
    if (std::isnan(_data[x])) final++;

  return final;
}
