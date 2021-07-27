#include "VECTOR3_FIELD_2D.h"
//#include <omp.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/gl.h> // OpenGL itself.
//#include <GL/glu.h> // GLU support library.
#include <GL/glut.h> // GLUT support library.
#endif



///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_2D::VECTOR3_FIELD_2D(const int& xRes, const int& yRes, const VEC3F& center, const VEC3F& lengths) :
  _xRes(xRes), _yRes(yRes), _center(center), _lengths(lengths), _initialized(true)
{
  _totalCells = _xRes * _yRes;
  _data = new VEC3F[_totalCells];

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
}

VECTOR3_FIELD_2D::VECTOR3_FIELD_2D(double* data, const int& xRes, const int& yRes, const VEC3F& center, const VEC3F& lengths) :
  _xRes(xRes), _yRes(yRes), _center(center), _lengths(lengths), _initialized(true)
{
  _totalCells = _xRes * _yRes;
  _data = new VEC3F[_totalCells];

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;

  for (int x = 0; x < _totalCells; x++)
  {
    _data[x][0] = data[3 * x];
    _data[x][1] = data[3 * x + 1];
  }
}

VECTOR3_FIELD_2D::VECTOR3_FIELD_2D(float* xData, float* yData, const int& xRes, const int& yRes, 
    const VEC3F& center, const VEC3F& lengths) :
  _xRes(xRes), _yRes(yRes), _center(center), _lengths(lengths), _initialized(true)
{
  _totalCells = _xRes * _yRes;
  _data = new VEC3F[_totalCells];

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;

  for (int x = 0; x < _totalCells; x++)
  {
    _data[x][0] = xData[x];
    _data[x][1] = yData[x];
  }
}

VECTOR3_FIELD_2D::VECTOR3_FIELD_2D(const VECTOR3_FIELD_2D& m) :
  _xRes(m.xRes()), _yRes(m.yRes()), _center(m.center()), _lengths(m.lengths()), _initialized(true)
{
  _totalCells = _xRes * _yRes;
  _data = new VEC3F[_totalCells];

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;

  for (int x = 0; x < _totalCells; x++)
    _data[x] = m[x];
}

VECTOR3_FIELD_2D::VECTOR3_FIELD_2D(const FIELD_2D& m) :
  _xRes(m.xRes()), _yRes(m.yRes()), _initialized(true)
{
  _totalCells = _xRes * _yRes;
  _data = new VEC3F[_totalCells];

  _lengths[0] = 1.0;
  _lengths[1] = 1.0;
  _lengths[2] = 1.0;

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;

  for (int x = 0; x < _totalCells; x++)
    _data[x] = 0;

  cout << " cells: " << _totalCells << endl;
  cout << " zero entry: " << _data[0] << endl;

  // do the middle x gradient
  for (int y = 0; y < _yRes; y++)
    for (int x = 1; x < _xRes - 1; x++)
      (*this)(x,y)[0] = (m(x+1, y) - m(x-1,y)) * 0.5;

  // do the left x gradient
  for (int y = 0; y < _yRes; y++)
    (*this)(0,y)[0] = m(1, y) - m(0,y);

  // do the right x gradient
  for (int y = 0; y < _yRes; y++)
    (*this)(_xRes - 1,y)[0] = m(_xRes - 1, y) - m(_xRes - 2,y);

  // do the middle y gradient
  for (int y = 1; y < _yRes - 1; y++)
    for (int x = 0; x < _xRes; x++)
      (*this)(x,y)[1] = (m(x, y+1) - m(x,y-1)) * 0.5;

  // do the bottom y gradient
  for (int x = 0; x < _xRes; x++)
    (*this)(x,0)[1] = m(x,1) - m(x,0);
  
  // do the top y gradient
  for (int x = 0; x < _xRes; x++)
    (*this)(x,_yRes - 1)[1] = m(x,_yRes - 1) - m(x,_yRes - 2);

  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
    {
      (*this)(x,y) *= 1.0 / _xRes;

      float xComp = (*this)(x,y)[0];

      (*this)(x,y)[0] = (*this)(x,y)[1];
      (*this)(x,y)[1] = -xComp;
    }
}


VECTOR3_FIELD_2D::VECTOR3_FIELD_2D() :
  _xRes(-1), _yRes(-1), _totalCells(-1), _data(NULL), _initialized(false)
{
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_2D::~VECTOR3_FIELD_2D()
{
  delete[] _data;
}
  
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_2D::clear()
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] = 0.0;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_2D& VECTOR3_FIELD_2D::operator=(const VECTOR3_FIELD_2D& input)
{
  if (input.xRes() != _xRes || input.yRes() != _yRes) 
  {
    delete[] _data;

    _xRes = input.xRes();
    _yRes = input.yRes();

    _totalCells = _xRes * _yRes;
    _data = new VEC3F[_totalCells];
  }

  _center = input.center();
  _lengths = input.lengths();

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  
  for (int x = 0; x < _totalCells; x++)
    _data[x] = input[x];

  _initialized = input._initialized;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D VECTOR3_FIELD_2D::scalarField(int component) const
{
  assert(component >= 0 && component < 2);

  FIELD_2D final(_xRes, _yRes);

  for (int x = 0; x < _totalCells; x++)
    final[x] = _data[x][component];

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_2D VECTOR3_FIELD_2D::magnitudeField() const
{
  FIELD_2D final(_xRes, _yRes);

  for (int x = 0; x < _totalCells; x++)
    final[x] = norm(_data[x]);

  return final;
}

///////////////////////////////////////////////////////////////////////
// take the field dot product
///////////////////////////////////////////////////////////////////////
FIELD_2D operator*(const VECTOR3_FIELD_2D&u, const VECTOR3_FIELD_2D& v)
{
  assert(u.xRes() == v.xRes());
  assert(u.yRes() == v.yRes());

  FIELD_2D final(u.xRes(), u.yRes());

  int totalCells = u.totalCells();
  for (int x = 0; x < totalCells; x++)
    final[x] = u[x] * v[x];

  return final;
}

///////////////////////////////////////////////////////////////////////
// take the field dot product
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_2D operator*(const FIELD_2D&u, const VECTOR3_FIELD_2D& v)
{
  assert(u.xRes() == v.xRes());
  assert(u.yRes() == v.yRes());

  VECTOR3_FIELD_2D final(v);

  int totalCells = u.totalCells();
  for (int x = 0; x < totalCells; x++)
  {
    Real uReal = u[x];
    VEC3F vVec = v[x];
    VEC3F product = uReal * vVec;
    final[x] = product;
  }

  return final;
}

///////////////////////////////////////////////////////////////////////
// lookup value at some real-valued spatial position
///////////////////////////////////////////////////////////////////////
const VEC3F VECTOR3_FIELD_2D::operator()(const VEC3F& position) const
{
  int x0 = (int)position[0];
  int x1    = x0 + 1;
  int y0 = (int)position[1];
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
  const Real s1 = position[0]- x0;
  const Real s0 = 1.0f - s1;
  const Real t1 = position[1]- y0;
  const Real t0 = 1.0f - t1;

  const int i00 = x0 + y0 * _xRes;
  const int i01 = x0 + y1 * _xRes;
  const int i10 = x1 + y0 * _xRes;
  const int i11 = x1 + y1 * _xRes;

  // interpolate
  // (indices could be computed once)
  return s0 * (t0 * _data[i00] + t1 * _data[i01]) +
         s1 * (t0 * _data[i10] + t1 * _data[i11]);
}

///////////////////////////////////////////////////////////////////////
// check if any entry is a nan
///////////////////////////////////////////////////////////////////////
bool VECTOR3_FIELD_2D::isNan()
{
  for (int x = 0; x < _totalCells; x++)
    for (int y = 0; y < 3; y++)
      if (isnan(_data[x][y]))
        return true;

  return false;
}

///////////////////////////////////////////////////////////////////////
// real-valued cell center coordinates
///////////////////////////////////////////////////////////////////////
VEC3F VECTOR3_FIELD_2D::cellCenter(int x, int y) const
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

  return final;
}

///////////////////////////////////////////////////////////////////////
// normalize all the vectors in the field
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_2D::normalize()
{
  for (int x = 0; x < _totalCells; x++)
    _data[x].normalize();
}

//////////////////////////////////////////////////////////////////////
// draw to GL
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_2D::draw()
{
  //float radius = 0.001;
  int stride = 1;
  Real scale = _dx;
  
  glBegin(GL_LINES);
    for (int y = 0; y < _yRes; y+=stride)
      for (int x = 0; x < _xRes; x+=stride)
      {
        const VEC3F& origin = cellCenter(x,y);
        //VEC3F endpoint = origin + scale * (*this)(x,y);
        VEC3F endpoint = origin + scale * (*this)[x + y * _xRes];

        glColor4f(1,1,1,1);
        glVertex3f(origin[0], origin[1], origin[2]);
        glColor4f(0,0,0,0);
        glVertex3f(endpoint[0], endpoint[1], endpoint[2]);
      }
  glEnd();
}

//////////////////////////////////////////////////////////////////////
// write a LIC image
// http://www.zhanpingliu.org/research/flowvis/LIC/MiniLIC.c
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_2D::writeLIC(int scaleUp, const char* filename)
{
  // create the filters
  int filterSize = 64;
  float* forwardFilter = new float[filterSize];
  float* backwardFilter = new float[filterSize];
  for(int i = 0; i < filterSize; i++)  
    forwardFilter[i] = backwardFilter[i] = i;

  float kernelLength = 10;
  VECTOR3_FIELD_2D& vectorField = *this;
  int xRes = vectorField.xRes() * scaleUp;
  int yRes = vectorField.yRes() * scaleUp;

  cout << " Writing LIC resolution: " << xRes << " " << yRes << endl;

  // generate a noise field
  FIELD_2D noiseField(xRes, yRes);
  for(int j = 0; j < yRes; j++)
    for(int i = 0; i < xRes; i++)
    { 
      int  r = rand();
      r = (  (r & 0xff) + ( (r & 0xff00) >> 8 )  ) & 0xff;
      noiseField(i,j) = (unsigned char) r;
    }

  // create the final image field
  FIELD_2D finalImage(xRes, yRes);

  int   maxAdvects = kernelLength * 3;
  float len2ID = (filterSize - 1) / kernelLength; ///map a curve LENgth TO an ID in the LUT

  ///for each pixel in the 2D output LIC image///
  for (int j = 0; j < yRes; j++)
    for (int i = 0; i < xRes; i++)
    { 
      ///init the composite texture accumulators and the weight accumulators///
      float textureAccum[] = {0,0};
      float weightAccum[]  = {0,0};
      float textureValue = 0;
    
      ///for either advection direction///
      for(int advectionDirection = 0; advectionDirection < 2; advectionDirection++)
      { 
        ///init the step counter, curve-length measurer, and streamline seed///
        int advects = 0;
        float currentLength = 0.0f;
        float clippedX0 = i + 0.5f;
        float clippedY0 = j + 0.5f;

        ///access the target filter LUT///
        float* weightLUT = (advectionDirection == 0) ? forwardFilter : backwardFilter;

        /// until the streamline is advected long enough or a tightly  spiralling center / focus is encountered///
        while (currentLength < kernelLength && advects < maxAdvects) 
        {
          ///access the vector at the sample
          VEC3F position((float)i / scaleUp, (float)j / scaleUp, 0);
          float vectorX = vectorField(position)[0];
          float vectorY = vectorField(position)[1];

          /// negate the vector for the backward-advection case///
          vectorX = (advectionDirection == 0) ? vectorX : -vectorX;
          vectorY = (advectionDirection == 0) ? vectorY : -vectorY;

          ///clip the segment against the pixel boundaries --- find the shorter from the two clipped segments///
          ///replace  all  if-statements  whenever  possible  as  they  might  affect the computational speed///
          const float lineSquare = 100000;
          const float vectorMin = 0.05;
          float segmentLength = lineSquare;
          segmentLength = (vectorX < -vectorMin) ? ( int(     clippedX0         ) - clippedX0 ) / vectorX : segmentLength;
          segmentLength = (vectorX >  vectorMin) ? ( int( int(clippedX0) + 1.5f ) - clippedX0 ) / vectorX : segmentLength;

          if (vectorY < -vectorMin)
          {
            float tmpLength = (int(clippedY0) - clippedY0) / vectorY;
            
            if (tmpLength < segmentLength) 
              segmentLength = tmpLength;
          }

          if (vectorY > vectorMin)
          {
            float tmpLength = (int(int(clippedY0) + 1.5f) - clippedY0) / vectorY;
            if (tmpLength <  segmentLength)
              segmentLength = tmpLength;
          }
          
          ///update the curve-length measurers///
          float previousLength  = currentLength;
          currentLength += segmentLength;
          segmentLength += 0.0004f;
         
          ///check if the filter has reached either end///
          segmentLength = (currentLength > kernelLength) ? ( (currentLength = kernelLength) - previousLength ) : segmentLength;

          ///obtain the next clip point///
          float clippedX1 = clippedX0 + vectorX * segmentLength;
          float clippedY1 = clippedY0 + vectorY * segmentLength;

          ///obtain the middle point of the segment as the texture-contributing sample///
          float sampleX = (clippedX0 + clippedX1) * 0.5f;
          float sampleY = (clippedY0 + clippedY1) * 0.5f;

          ///obtain the texture value of the sample///
          textureValue = noiseField(sampleX, sampleY);

          ///update the accumulated weight and the accumulated composite texture (texture x weight)
          float currentWeightAccum = weightLUT[ int(currentLength * len2ID) ];
          float sampleWeight = currentWeightAccum - weightAccum[advectionDirection];     
          weightAccum[advectionDirection] = currentWeightAccum;               
          textureAccum[advectionDirection] += textureValue * sampleWeight;
        
          ///update the step counter and the "current" clip point
          advects++;
          clippedX0 = clippedX1;
          clippedY0 = clippedY1;

          ///check if the streamline has gone beyond the flow field
          if (clippedX0 < 0.0f || clippedX0 >= xRes ||
              clippedY0 < 0.0f || clippedY0 >= yRes)  break;
        } 
      }

      ///normalize the accumulated composite texture
      textureValue = (textureAccum[0] + textureAccum[1]) / (weightAccum[0] + weightAccum[1]);

      ///clamp the texture value against the displayable intensity range [0, 255]
      textureValue = (textureValue <   0.0f) ?   0.0f : textureValue;
      textureValue = (textureValue > 255.0f) ? 255.0f : textureValue; 
      finalImage(i,j) = textureValue / 255.0;
    }

  finalImage.writeJPG(filename);
  delete[] forwardFilter;
  delete[] backwardFilter;
}

//////////////////////////////////////////////////////////////////////
// write a LIC image
// http://www.zhanpingliu.org/research/flowvis/LIC/MiniLIC.c
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_2D::writeLIC(const char* filename)
{
  // create the filters
  int filterSize = 64;
  float* forwardFilter = new float[filterSize];
  float* backwardFilter = new float[filterSize];
  for(int i = 0; i < filterSize; i++)  
    forwardFilter[i] = backwardFilter[i] = i;

  float kernelLength = 10;
  VECTOR3_FIELD_2D& vectorField = *this;
  int xRes = vectorField.xRes();
  int yRes = vectorField.yRes();

  // generate a noise field
  FIELD_2D noiseField(xRes, yRes);
  for(int j = 0; j < yRes; j++)
    for(int i = 0; i < xRes; i++)
    { 
      int  r = rand();
      r = (  (r & 0xff) + ( (r & 0xff00) >> 8 )  ) & 0xff;
      noiseField(i,j) = (unsigned char) r;
    }

  // create the final image field
  FIELD_2D finalImage(xRes, yRes);

  int   maxAdvects = kernelLength * 3;
  float len2ID = (filterSize - 1) / kernelLength; ///map a curve LENgth TO an ID in the LUT

  ///for each pixel in the 2D output LIC image///
  for (int j = 0; j < yRes; j++)
    for (int i = 0; i < xRes; i++)
    { 
      ///init the composite texture accumulators and the weight accumulators///
      float textureAccum[] = {0,0};
      float weightAccum[]  = {0,0};
      float textureValue = 0;
    
      ///for either advection direction///
      for(int advectionDirection = 0; advectionDirection < 2; advectionDirection++)
      { 
        ///init the step counter, curve-length measurer, and streamline seed///
        int advects = 0;
        float currentLength = 0.0f;
        float clippedX0 = i + 0.5f;
        float clippedY0 = j + 0.5f;

        ///access the target filter LUT///
        float* weightLUT = (advectionDirection == 0) ? forwardFilter : backwardFilter;

        /// until the streamline is advected long enough or a tightly  spiralling center / focus is encountered///
        while (currentLength < kernelLength && advects < maxAdvects) 
        {
          ///access the vector at the sample
          float vectorX = vectorField(i,j)[0];
          float vectorY = vectorField(i,j)[1];

          /// negate the vector for the backward-advection case///
          vectorX = (advectionDirection == 0) ? vectorX : -vectorX;
          vectorY = (advectionDirection == 0) ? vectorY : -vectorY;

          ///clip the segment against the pixel boundaries --- find the shorter from the two clipped segments///
          ///replace  all  if-statements  whenever  possible  as  they  might  affect the computational speed///
          const float lineSquare = 100000;
          const float vectorMin = 0.05;
          float segmentLength = lineSquare;
          segmentLength = (vectorX < -vectorMin) ? ( int(     clippedX0         ) - clippedX0 ) / vectorX : segmentLength;
          segmentLength = (vectorX >  vectorMin) ? ( int( int(clippedX0) + 1.5f ) - clippedX0 ) / vectorX : segmentLength;

          if (vectorY < -vectorMin)
          {
            float tmpLength = (int(clippedY0) - clippedY0) / vectorY;
            
            if (tmpLength < segmentLength) 
              segmentLength = tmpLength;
          }

          if (vectorY > vectorMin)
          {
            float tmpLength = (int(int(clippedY0) + 1.5f) - clippedY0) / vectorY;
            if (tmpLength <  segmentLength)
              segmentLength = tmpLength;
          }
          
          ///update the curve-length measurers///
          float previousLength  = currentLength;
          currentLength += segmentLength;
          segmentLength += 0.0004f;
         
          ///check if the filter has reached either end///
          segmentLength = (currentLength > kernelLength) ? ( (currentLength = kernelLength) - previousLength ) : segmentLength;

          ///obtain the next clip point///
          float clippedX1 = clippedX0 + vectorX * segmentLength;
          float clippedY1 = clippedY0 + vectorY * segmentLength;

          ///obtain the middle point of the segment as the texture-contributing sample///
          float sampleX = (clippedX0 + clippedX1) * 0.5f;
          float sampleY = (clippedY0 + clippedY1) * 0.5f;

          ///obtain the texture value of the sample///
          textureValue = noiseField(sampleX, sampleY);

          ///update the accumulated weight and the accumulated composite texture (texture x weight)
          float currentWeightAccum = weightLUT[ int(currentLength * len2ID) ];
          float sampleWeight = currentWeightAccum - weightAccum[advectionDirection];     
          weightAccum[advectionDirection] = currentWeightAccum;               
          textureAccum[advectionDirection] += textureValue * sampleWeight;
        
          ///update the step counter and the "current" clip point
          advects++;
          clippedX0 = clippedX1;
          clippedY0 = clippedY1;

          ///check if the streamline has gone beyond the flow field
          if (clippedX0 < 0.0f || clippedX0 >= xRes ||
              clippedY0 < 0.0f || clippedY0 >= yRes)  break;
        } 
      }

      ///normalize the accumulated composite texture
      textureValue = (textureAccum[0] + textureAccum[1]) / (weightAccum[0] + weightAccum[1]);

      ///clamp the texture value against the displayable intensity range [0, 255]
      textureValue = (textureValue <   0.0f) ?   0.0f : textureValue;
      textureValue = (textureValue > 255.0f) ? 255.0f : textureValue; 
      finalImage(i,j) = textureValue / 255.0;
    }

  finalImage.writeJPG(filename);
  delete[] forwardFilter;
  delete[] backwardFilter;
}
