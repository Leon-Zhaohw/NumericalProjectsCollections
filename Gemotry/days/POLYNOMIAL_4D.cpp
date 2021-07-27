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
#include "POLYNOMIAL_4D.h"
#include <assert.h>
#include "TIMER.h"

POLYNOMIAL_4D::POLYNOMIAL_4D(const vector<QUATERNION>& roots)
{
  _powerScalar = 1.0;

  _roots = roots;
  _totalRoots = _roots.size();
  _coeffs.resize(_totalRoots + 1);

  if (roots.size() == 0)
    return;

  _rootPowers.resize(_totalRoots);
  for (int x = 0; x < _totalRoots; x++)
    _rootPowers[x] = 1.0;

  computeCoeffsFast();
  computeDerivativeCoeffs();
}

POLYNOMIAL_4D::POLYNOMIAL_4D(const vector<QUATERNION>& roots, const vector<Real> powers)
{
  assert(roots.size() == powers.size());
  _powerScalar = 1.0;

  _roots = roots;
  _totalRoots = _roots.size();
  _coeffs.resize(_totalRoots + 1);

  if (roots.size() == 0)
    return;

  _rootPowers.resize(_totalRoots);
  for (int x = 0; x < _totalRoots; x++)
    _rootPowers[x] = powers[x];

  computeCoeffsFast();
  computeDerivativeCoeffs();
}

///////////////////////////////////////////////////////////////////////
// set the coefficients directly - powers increase with index, i.e.
// coeffs[0] -> x^0
// coeffs[1] -> x^1
///////////////////////////////////////////////////////////////////////
POLYNOMIAL_4D::POLYNOMIAL_4D(const vector<float>& coeffs)
{
  _coeffs.resize(coeffs.size());
  _powerScalar = 1.0;

  for (unsigned int x = 0; x < coeffs.size(); x++)
    _coeffs[x] = QUATERNION(coeffs[x], 0);
  _totalRoots = _coeffs.size() - 1;

  _rootPowers.resize(_totalRoots);
  for (int x = 0; x < _totalRoots; x++)
    _rootPowers[x] = 1.0;

  computeDerivativeCoeffs();
}

POLYNOMIAL_4D::POLYNOMIAL_4D()
{
  _totalRoots = -1;
  _powerScalar = 1.0;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
POLYNOMIAL_4D::~POLYNOMIAL_4D()
{
}

///////////////////////////////////////////////////////////////////////
// m choose n
///////////////////////////////////////////////////////////////////////
vector<QUATERNION> POLYNOMIAL_4D::choose(const vector<QUATERNION> m, const int n)
{
  if (n == 1)
    return m;

  vector<QUATERNION> final;
  for (unsigned int x = 0; x < m.size(); x++)
  {
    vector<QUATERNION> subset;
    for (unsigned int y = x + 1; y < m.size(); y++)
      subset.push_back(m[y]);

    vector<QUATERNION> chooseOneLess = choose(subset, n - 1);

    for (unsigned int y = 0; y < chooseOneLess.size(); y++)
      final.push_back(chooseOneLess[y] * m[x]);
  }
  return final;
}

//////////////////////////////////////////////////////////////////////
// resize the polynomial to a given number of roots
//////////////////////////////////////////////////////////////////////
void POLYNOMIAL_4D::resizeAndWipe(int totalRoots)
{
  _roots.clear();
  for (int x = 0; x < totalRoots; x++)
    _roots.push_back(QUATERNION(0,0,0,0));

  _totalRoots = totalRoots;
  _rootPowers.resize(_totalRoots);
  for (int x = 0; x < totalRoots; x++)
    _rootPowers[x] = 1.0;

  computeCoeffsFast();
  computeDerivativeCoeffs();
}

//////////////////////////////////////////////////////////////////////
// generically compute new root coefficients
//////////////////////////////////////////////////////////////////////
void POLYNOMIAL_4D::computeCoeffs()
{
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " SHOULDN'T BE CALLING THIS " << endl;
  assert(_totalRoots > 0);

  vector<QUATERNION> chosen;

  for (int x = 0; x < _totalRoots; x++)
    _coeffs[x] = 0;

  for (int x = 0; x < _totalRoots; x++)
  {
    chosen = choose(_roots, _totalRoots - x);
    for (unsigned int y = 0; y < chosen.size(); y++)
      _coeffs[x] += chosen[y];

    if (x % 2 != _totalRoots % 2)
      _coeffs[x] *= -1;
  }
  //_coeffs[_totalRoots] = QUATERNION(1.0,0);
  _coeffs[_totalRoots] = QUATERNION(0,0,0,1.0);
}

//////////////////////////////////////////////////////////////////////
// compute the polynomial coefficients
//////////////////////////////////////////////////////////////////////
void POLYNOMIAL_4D::computeCoeffsFast()
{
  //TIMER functionTimer(__FUNCTION__);
  vector<QUATERNION> coeffs;
  vector<QUATERNION> coeffsOld;

  coeffs.push_back(QUATERNION(1.0, 0.0));
  coeffs.push_back(-1.0 * _roots[0]);

  for (int x = 1; x < _totalRoots; x++)
  {
    coeffsOld = coeffs;
    coeffs.clear();

    QUATERNION alpha = _roots[x];
    coeffs.push_back(coeffsOld[0]);

    for (unsigned int y = 1; y < coeffsOld.size(); y++)
      coeffs.push_back(coeffsOld[y] - alpha * coeffsOld[y - 1]);

    coeffs.push_back(-1.0 * alpha * coeffsOld.back());
  }

  _coeffs.clear();
  for (unsigned int x = 0; x < coeffs.size(); x++)
    _coeffs.push_back(coeffs[coeffs.size() - 1 - x]);
}

//////////////////////////////////////////////////////////////////////
// compute the derivative coefficients
//////////////////////////////////////////////////////////////////////
void POLYNOMIAL_4D::computeDerivativeCoeffs()
{
  _derivs.resize(_totalRoots);
  _secondDerivs.resize(_totalRoots - 1);

  for (int x = 0; x < _totalRoots; x++)
    _derivs[x] = (x + 1.0) * _coeffs[x + 1];
  for (int x = 0; x < _totalRoots - 1; x++)
    _secondDerivs[x] = (x + 1.0) * _derivs[x + 1];
}

//////////////////////////////////////////////////////////////////////
// evaluate the polynomial at a point
//////////////////////////////////////////////////////////////////////
QUATERNION POLYNOMIAL_4D::evaluate(const QUATERNION& point) const
{
  assert(_totalRoots > 0);

  int roots = totalRoots();
  QUATERNION final(_coeffs[roots]);
  for (int x = roots - 1; x >= 0; x--)
    //g = g * point + topCoeffs[x];
    final.multiplyAdd(point, _coeffs[x]);

  return final;
}

//////////////////////////////////////////////////////////////////////
// evaluate the derivative at a point
//////////////////////////////////////////////////////////////////////
QUATERNION POLYNOMIAL_4D::evaluateDerivative(const QUATERNION& point) const
{
  assert(_totalRoots > 0);

  int roots = totalRoots();
  QUATERNION final(_derivs[roots - 1]);
  for (int x = roots - 2; x >= 0; x--)
    //gPrime = gPrime * point + topDerivs[x];
    final.multiplyAdd(point, _derivs[x]);

  return final;
}

//////////////////////////////////////////////////////////////////////
// evaluate the second derivative at a point
//////////////////////////////////////////////////////////////////////
QUATERNION POLYNOMIAL_4D::evaluateSecondDerivative(const QUATERNION& point) const
{
  assert(_totalRoots > 0);

  // compute all the needed powers
  vector<QUATERNION> powers(_totalRoots);
  powers[0] = point;
  for (int x = 1; x < _totalRoots; x++)
    powers[x] = powers[x - 1] * point;

  // compute the second derivative
  QUATERNION final = _secondDerivs[0];
  for (int x = 1; x < _totalRoots - 1; x++)
    final += _secondDerivs[x] * powers[x-1];

  return final;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void POLYNOMIAL_4D::evaluateMultiple(const QUATERNION& point, QUATERNION& poly, QUATERNION& deriv) const
{
  assert(_totalRoots > 0);

  // compute all the needed powers
  vector<QUATERNION> powers(_totalRoots);
  powers[0] = point;
  for (int x = 1; x < _totalRoots; x++)
    powers[x] = powers[x - 1] * point;

  // compute the polynomial
  poly = _coeffs[0];
  for (int x = 1; x < _totalRoots + 1; x++)
    poly += _coeffs[x] * powers[x - 1];
  
  // compute the derivative
  deriv = _derivs[0];
  for (int x = 1; x < _totalRoots; x++)
    deriv += _derivs[x] * powers[x-1];
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void POLYNOMIAL_4D::evaluateMultiple(const QUATERNION& point, QUATERNION& poly, QUATERNION& deriv, QUATERNION& secondDeriv) const
{
  assert(_totalRoots > 0);

  // compute all the needed powers
  vector<QUATERNION> powers(_totalRoots);
  powers[0] = point;
  for (int x = 1; x < _totalRoots; x++)
    powers[x] = powers[x - 1] * point;

  // compute the polynomial
  poly = _coeffs[0];
  for (int x = 1; x < _totalRoots + 1; x++)
    poly += _coeffs[x] * powers[x - 1];
  
  // compute the derivative
  deriv = _derivs[0];
  for (int x = 1; x < _totalRoots; x++)
    deriv += _derivs[x] * powers[x-1];
  
  // compute the second derivative
  secondDeriv = _secondDerivs[0];
  for (int x = 1; x < _totalRoots - 1; x++)
    secondDeriv += _secondDerivs[x] * powers[x-1];
}

//////////////////////////////////////////////////////////////////////
// evaluate rational function
//////////////////////////////////////////////////////////////////////
void POLYNOMIAL_4D::evaluateRational(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, const QUATERNION& point, QUATERNION& p, QUATERNION& pPrime, QUATERNION& pPrimePrime)
{
  QUATERNION g, gPrime, gPrimePrime;
  QUATERNION h, hPrime, hPrimePrime;

  const int topRoots = top.totalRoots();
  const int bottomRoots = bottom.totalRoots();
  const int maxRoots = (topRoots > bottomRoots) ? topRoots : bottomRoots;

  // compute all the needed powers
  vector<QUATERNION> powers(maxRoots);
  powers[0] = point;
  for (int x = 1; x < maxRoots; x++)
    powers[x] = powers[x - 1] * point;
  
  // compute the g polynomial
  const vector<QUATERNION>& topCoeffs = top.coeffs();
  const vector<QUATERNION>& topDerivs = top.derivs();
  const vector<QUATERNION>& topSecond = top.secondDerivs();
  g = topCoeffs[0];
  for (int x = 1; x < topRoots + 1; x++)
    g += topCoeffs[x] * powers[x - 1];
  gPrime = topDerivs[0];
  for (int x = 1; x < topRoots; x++)
    gPrime += topDerivs[x] * powers[x-1];
  gPrimePrime = topSecond[0];
  for (int x = 1; x < topRoots- 1; x++)
    gPrimePrime += topSecond[x] * powers[x-1];

  // compute the h polynomial
  const vector<QUATERNION>& bottomCoeffs = bottom.coeffs();
  const vector<QUATERNION>& bottomDerivs = bottom.derivs();
  const vector<QUATERNION>& bottomSecond = bottom.secondDerivs();
  h = bottomCoeffs[0];
  for (int x = 1; x < bottomRoots + 1; x++)
    h += bottomCoeffs[x] * powers[x - 1];
  hPrime = bottomDerivs[0];
  for (int x = 1; x < bottomRoots; x++)
    hPrime += bottomDerivs[x] * powers[x-1];
  hPrimePrime = bottomSecond[0];
  for (int x = 1; x < bottomRoots- 1; x++)
    hPrimePrime += bottomSecond[x] * powers[x-1];
 
  // compute the rational derivatives 
  const QUATERNION numerator = gPrime * h - g * hPrime;
  const QUATERNION hSq = h * h;

  const QUATERNION inverse = h.inverse();
  const QUATERNION inverseSq = inverse * inverse;
  const QUATERNION inverseCubed = inverseSq * inverse;

  p = g * inverse;
  pPrime = numerator * inverseSq;

  //const QUATERNION hCubed = hSq * h;
  const QUATERNION first = gPrimePrime * hSq;
  const QUATERNION second = -2.0 * gPrime * hPrime * h;
  const QUATERNION third = g * (2.0 * hPrime * hPrime - h * hPrimePrime);

  pPrimePrime = (first + second + third) * inverseCubed;
}

//////////////////////////////////////////////////////////////////////
// evaluate rational function
//////////////////////////////////////////////////////////////////////
void POLYNOMIAL_4D::evaluateRational(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, const QUATERNION& point, QUATERNION& p, QUATERNION& pPrime)
{
  QUATERNION g = top.evaluate(point);
  QUATERNION gPrime = top.evaluateDerivative(point);
  QUATERNION h = bottom.evaluate(point);

  QUATERNION hPrime = bottom.evaluateDerivative(point);
 
  // compute the rational derivatives 
  const QUATERNION numerator = gPrime * h - g * hPrime;
  const QUATERNION inverse = h.inverse();
  const QUATERNION inverseSq = inverse * inverse;
  p = g * inverse;
  pPrime = numerator * inverseSq;
}

//////////////////////////////////////////////////////////////////////
// run a unit test for a known rational function
//////////////////////////////////////////////////////////////////////
void POLYNOMIAL_4D::rationalTest()
{
  vector<float> coeffs;
  coeffs.push_back(-1.0);
  coeffs.push_back(0.0);
  coeffs.push_back(1.0);
  POLYNOMIAL_4D top = POLYNOMIAL_4D(coeffs);

  coeffs.clear();
  coeffs.push_back(0);
  coeffs.push_back(0);
  coeffs.push_back(1);
  POLYNOMIAL_4D bottom = POLYNOMIAL_4D(coeffs);

  QUATERNION point(1.23,4.56);

  QUATERNION g, gPrime, gPrimePrime;
  top.evaluateMultiple(point, g, gPrime, gPrimePrime);
  QUATERNION h, hPrime, hPrimePrime;
  bottom.evaluateMultiple(point, h, hPrime, hPrimePrime);

  QUATERNION rational = g / h;

  // compute (1 - 1 / z^2) directly
  QUATERNION ground = QUATERNION(1,0) - QUATERNION(1,0) / (point * point);

  QUATERNION diff = rational - ground;

  cout << " ground: " << ground << endl;
  cout << " computed: " << rational << endl;
  cout << " diff: " << diff.magnitude() << endl;

  QUATERNION numerator = gPrime * h - g * hPrime;
  QUATERNION hSq = h * h;
  QUATERNION pPrime = numerator / hSq;

  // compute (2 / z^3) directly
  ground = QUATERNION(2,0) / (point * point * point);
  diff = pPrime - ground;
  
  cout << " ground derivative: " << ground << endl;
  cout << " computed derivative: " << pPrime << endl;
  cout << " diff: " << diff.magnitude() << endl;

  QUATERNION hCubed = hSq * h;
  QUATERNION first = gPrimePrime * hSq;
  QUATERNION second = 2.0 * gPrime * hPrime * h;
  QUATERNION third = g * (2.0 * hPrime * hPrime - h * hPrimePrime);
  QUATERNION pPrimePrime = (first - second + third) / hCubed;

  // compute (-6 / z^4) directly
  ground = QUATERNION(-6,0) / (point * point * point * point);
  diff = pPrimePrime - ground;
  cout << " ground second derivative: " << ground << endl;
  cout << " computed second derivative: " << pPrimePrime << endl;
  cout << " diff: " << diff.magnitude() << endl;

}

//////////////////////////////////////////////////////////////////////
// add a new root
//////////////////////////////////////////////////////////////////////
void POLYNOMIAL_4D::addRoot(const QUATERNION& newRoot)
{
  _roots.push_back(newRoot);
  _totalRoots = _roots.size();
  _coeffs.resize(_totalRoots + 1);

  _rootPowers.push_back(1.0);

  computeCoeffsFast();
  computeDerivativeCoeffs();
}

//////////////////////////////////////////////////////////////////////
// add a new root
//////////////////////////////////////////////////////////////////////
void POLYNOMIAL_4D::addFrontRoot(const QUATERNION& newRoot)
{
  vector<QUATERNION> roots;
  roots.push_back(newRoot);
  for (unsigned int x = 0; x < _roots.size(); x++)
    roots.push_back(_roots[x]);
  _roots = roots;

  _totalRoots = _roots.size();
  _coeffs.resize(_totalRoots + 1);

  vector<Real> rootPowers;
  rootPowers.push_back(1.0);
  for (unsigned int x = 0; x < _rootPowers.size(); x++)
    rootPowers.push_back(_rootPowers[x]);
  _rootPowers = rootPowers;

  computeCoeffsFast();
  computeDerivativeCoeffs();
}

//////////////////////////////////////////////////////////////////////
// add a new root
//////////////////////////////////////////////////////////////////////
void POLYNOMIAL_4D::addRoot(const QUATERNION& newRoot, const Real& power)
{
  _roots.push_back(newRoot);
  _totalRoots = _roots.size();
  _coeffs.resize(_totalRoots + 1);

  _rootPowers.push_back(power);

  computeCoeffsFast();
  computeDerivativeCoeffs();
}

//////////////////////////////////////////////////////////////////////
// modify an existing root
//////////////////////////////////////////////////////////////////////
void POLYNOMIAL_4D::modifyRoot(const int whichRoot, const QUATERNION& newRoot)
{
  assert(whichRoot < _totalRoots);

  if (whichRoot >= _totalRoots) return;

  _roots[whichRoot] = newRoot;
  computeCoeffsFast();
  computeDerivativeCoeffs();
}

//////////////////////////////////////////////////////////////////////
// file IO
//////////////////////////////////////////////////////////////////////
void POLYNOMIAL_4D::writeTextfile(FILE* file) const
{
  fprintf(file, "total roots = %i\n", _totalRoots);

  int size = _roots.size();
  fprintf(file, "roots size = %i\n", size);
  for (int x = 0; x < size; x++)
    _roots[x].writeTextfile(file);

  size = _rootPowers.size();
  for (int x = 0; x < size; x++)
  {
    double power = _rootPowers[x];
    fprintf(file, "%.20e\n", power);
  }

  double powerScalar = _powerScalar;
  fprintf(file, "power scale = %.20e\n", powerScalar);
}

//////////////////////////////////////////////////////////////////////
// file IO
//////////////////////////////////////////////////////////////////////
void POLYNOMIAL_4D::readTextfile(FILE* file)
{
  fscanf(file, "total roots = %i\n", &_totalRoots);

  int size;
  fscanf(file, "roots size = %i\n", &size);
  _roots.resize(size);
  _rootPowers.resize(size);

  for (int x = 0; x < size; x++)
    _roots[x].readTextfile(file);

  for (int x = 0; x < size; x++)
    fscanf(file, "%lf\n", &(_rootPowers[x]));

  fscanf(file, "power scale = %lf\n", &_powerScalar);
}

//////////////////////////////////////////////////////////////////////
// file IO
//////////////////////////////////////////////////////////////////////
void POLYNOMIAL_4D::write(FILE* file) const
{
  fwrite((void*)&_totalRoots, sizeof(int), 1, file);

  int size = _coeffs.size();
  fwrite((void*)&size, sizeof(int), 1, file);
  for (int x = 0; x < size; x++)
    _coeffs[x].write(file);

  size = _derivs.size();
  fwrite((void*)&size, sizeof(int), 1, file);
  for (int x = 0; x < size; x++)
    _derivs[x].write(file);

  size = _secondDerivs.size();
  fwrite((void*)&size, sizeof(int), 1, file);
  for (int x = 0; x < size; x++)
    _secondDerivs[x].write(file);
  
  size = _roots.size();
  fwrite((void*)&size, sizeof(int), 1, file);
  for (int x = 0; x < size; x++)
    _roots[x].write(file);

  size = _rootPowers.size();
  fwrite((void*)&size, sizeof(int), 1, file);
  for (int x = 0; x < size; x++)
  {
    double power = _rootPowers[x];
    fwrite((void*)&power, sizeof(double), 1, file);
  }

  double powerScalar = _powerScalar;
  fwrite((void*)&powerScalar, sizeof(double), 1, file);
}

//////////////////////////////////////////////////////////////////////
// file IO
//////////////////////////////////////////////////////////////////////
void POLYNOMIAL_4D::read(FILE* file)
{
  fread((void*)&_totalRoots, sizeof(int), 1, file);
  cout << " Reading in 4D polynomial with " << _totalRoots << " roots " << endl;

  int size;
  fread((void*)&size, sizeof(int), 1, file);
  _coeffs.resize(size);
  for (int x = 0; x < size; x++)
    _coeffs[x].read(file);

  fread((void*)&size, sizeof(int), 1, file);
  _derivs.resize(size);
  for (int x = 0; x < size; x++)
    _derivs[x].read(file);

  fread((void*)&size, sizeof(int), 1, file);
  _secondDerivs.resize(size);
  for (int x = 0; x < size; x++)
    _secondDerivs[x].read(file);

  fread((void*)&size, sizeof(int), 1, file);
  _roots.resize(size);
  for (int x = 0; x < size; x++)
    _roots[x].read(file);

  fread((void*)&size, sizeof(int), 1, file);
  _rootPowers.resize(size);
  for (int x = 0; x < _totalRoots; x++)
  {
    double power = 0;
    fread((void*)&power, sizeof(double), 1, file);
    _rootPowers[x] =power;
  }

  double powerScalar;
  fread((void*)&powerScalar, sizeof(double), 1, file);
  _powerScalar = powerScalar;
}

//////////////////////////////////////////////////////////////////////
// Print polynomial to stream
//////////////////////////////////////////////////////////////////////
ostream& operator<<(ostream &out, const POLYNOMIAL_4D& poly)
{
  for (unsigned int x = 0; x < poly.roots().size(); x++)
    out << " Root " << x << ": " << poly.roots()[x] << endl;

  return out;
}

//////////////////////////////////////////////////////////////////////
// get the condition number of the polynomial
//////////////////////////////////////////////////////////////////////
Real POLYNOMIAL_4D::conditionNumber()
{
  vector<Real> derivs;

  for (int x = 0; x < _totalRoots; x++)
  {
    QUATERNION derivative = evaluateDerivative(_roots[x]);
    derivs.push_back(derivative.magnitude());
  }

  vector<QUATERNION> rootPowers = _roots;
  Real maxFound = 0;
  for (int i = 0; i < _totalRoots; i++)
  {
    for (int j = 0; j < _totalRoots; j++)
    {
      Real condition = (_coeffs[i] * rootPowers[j]).magnitude() / derivs[j];
      maxFound = (condition < maxFound) ? maxFound : condition;
    }

    // add another power to the roots
    for (int j = 0; j < _totalRoots; j++)
      rootPowers[j] = rootPowers[j] * _roots[j];
  }

  return maxFound;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void POLYNOMIAL_4D::computeNestedCoeffs()
{
  _ws.clear();

  // the first w is the same as the root
  _ws.push_back(_roots[0]);

  QUATERNION w = _roots[1];
  _ws.push_back(w);
  w = _roots[2];
  _ws.push_back(w);
}

//////////////////////////////////////////////////////////////////////
// use the brute force nested formulation
//////////////////////////////////////////////////////////////////////
QUATERNION POLYNOMIAL_4D::evaluateFactoredPositive(const QUATERNION& point) const
{
  QUATERNION result = point + _roots[0];

  for (int x = 1; x < _totalRoots; x++)
  {
    result = result * (point + _roots[x]);
  }

  return result;
}

//////////////////////////////////////////////////////////////////////
// use the brute force nested formulation
//////////////////////////////////////////////////////////////////////
QUATERNION POLYNOMIAL_4D::evaluateFactoredDouble(const QUATERNION& point) const
{
  QUATERNION result = point - _roots[0];

  int power = 2;
  for (int x = 0; x < power; x++)
    result = result * result;

  for (int x = 1; x < _totalRoots; x++)
  {
    QUATERNION currentRoot = (point - _roots[x]);
    currentRoot = currentRoot * currentRoot;
    result = currentRoot * result;
  }

  return result;
}

//////////////////////////////////////////////////////////////////////
// use the brute force nested formulation
//
// this is where the code spends most of its time.
// Worthwhile to unroll this? - Doesn't seem to make a difference
//////////////////////////////////////////////////////////////////////
QUATERNION POLYNOMIAL_4D::evaluateScaledPowerFactored(const QUATERNION& point) const
{
  QUATERNION result = point - _roots[0];
  result = result.pow(_powerScalar * _rootPowers[0]);

  for (int x = 1; x < _totalRoots; x++)
  {
    QUATERNION term = (point - _roots[x]);
    result *= term.pow(_powerScalar * _rootPowers[x]);
  }

  return result;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
QUATERNION POLYNOMIAL_4D::evaluateRoot(const QUATERNION& point, const int index, const bool debug) const
{
  QUATERNION result = point - _roots[index];
  result = result.pow(_powerScalar * _rootPowers[index]);
  if (debug)
  {
    cout << "\t diff:  " << point - _roots[index]<< endl;
    cout << "\t power: " << _powerScalar * _rootPowers[index] << endl;
    cout << "\t final: " << result << endl;
  }

  return result;
}

//////////////////////////////////////////////////////////////////////
// use the brute force nested formulation
//////////////////////////////////////////////////////////////////////
QUATERNION POLYNOMIAL_4D::evaluateScaledPowerFactoredVerbose(const QUATERNION& point) const
{
  assert(_roots.size() == _rootPowers.size());
  QUATERNION result = point - _roots[0];
  result = result.pow(_powerScalar * _rootPowers[0]);

  for (int x = 1; x < _totalRoots; x++)
  {
    cout << " result: " << result << endl;
    QUATERNION term = (point - _roots[x]);
    cout << " diff: " << term << endl;
    result *= term.pow(_powerScalar * _rootPowers[x]);
  }
  cout << " final: " << result << endl;

  return result;
}

//////////////////////////////////////////////////////////////////////
// use the brute force nested formulation
//////////////////////////////////////////////////////////////////////
QUATERNION POLYNOMIAL_4D::evaluatePowerFactored(const QUATERNION& point) const
{
  assert(_roots.size() == _rootPowers.size());
  QUATERNION result = point - _roots[0];
  result = result.pow(_rootPowers[0]);

  for (int x = 1; x < _totalRoots; x++)
  {
    QUATERNION term = (point - _roots[x]);
    result *= term.pow(_rootPowers[x]);
  }

  return result;
}

//////////////////////////////////////////////////////////////////////
// use the brute force nested formulation, but cache the
// multiplies
//////////////////////////////////////////////////////////////////////
QUATERNION POLYNOMIAL_4D::evaluatePowerFactored(const QUATERNION& point, vector<QUATERNION>& forward, vector<QUATERNION>& backward) const
{
  assert(_roots.size() == _rootPowers.size());
  vector<QUATERNION> powers;
  forward.clear();
  backward.clear();

  QUATERNION result = point - _roots[0];
  result = result.pow(_rootPowers[0]);
  powers.push_back(result);
  forward.push_back(result);

  for (int x = 1; x < _totalRoots; x++)
  {
    QUATERNION term = (point - _roots[x]).pow(_rootPowers[x]);
    powers.push_back(term);
    result *= term;
    forward.push_back(result);
  }

  // cache out the multiplies
  QUATERNION current = powers[_totalRoots - 1];
  backward.resize(_totalRoots);
  backward[_totalRoots - 1] = current;
  for (int x = _totalRoots - 2; x >= 0; x--)
  {
    current = powers[x] * current;
    backward[x] = current;
  }

  return result;
}

//////////////////////////////////////////////////////////////////////
// use the brute force nested formulation, but cache the
// multiplies
//////////////////////////////////////////////////////////////////////
QUATERNION POLYNOMIAL_4D::evaluateScaledPowerFactored(const QUATERNION& point, vector<QUATERNION>& forward, vector<QUATERNION>& backward) const
{
  assert(_roots.size() == _rootPowers.size());
  vector<QUATERNION> powers;
  forward.clear();
  backward.clear();

  QUATERNION result = point - _roots[0];
  result = result.pow(_powerScalar * _rootPowers[0]);
  powers.push_back(result);
  forward.push_back(result);

  for (int x = 1; x < _totalRoots; x++)
  {
    QUATERNION term = (point - _roots[x]).pow(_powerScalar * _rootPowers[x]);
    powers.push_back(term);
    result *= term;
    forward.push_back(result);
  }

  // cache out the multiplies
  QUATERNION current = powers[_totalRoots - 1];
  backward.resize(_totalRoots);
  backward[_totalRoots - 1] = current;
  for (int x = _totalRoots - 2; x >= 0; x--)
  {
    current = powers[x] * current;
    backward[x] = current;
  }

  return result;
}

//////////////////////////////////////////////////////////////////////
// use the brute force nested formulation
//////////////////////////////////////////////////////////////////////
QUATERNION POLYNOMIAL_4D::evaluateFactored(const QUATERNION& point) const
{
  QUATERNION result = point - _roots[0];

  for (int x = 1; x < _totalRoots; x++)
    result = result * (point - _roots[x]);

  return result;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void POLYNOMIAL_4D::evaluateFactoredRational(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, const QUATERNION& point, QUATERNION& p, QUATERNION& pPrime)
{
  const vector<QUATERNION> topRoots = top.roots();
  assert(topRoots.size() > 0);

  // assume ordering is (x - _roots[0]) * (x - _roots[1]) ...
  QUATERNION g = (point - topRoots[0]);
  for (unsigned int x = 1; x < topRoots.size(); x++)
    g = g * (point - topRoots[x]);

  const vector<QUATERNION> bottomRoots = bottom.roots();
  assert(bottomRoots.size() > 0);

  // assume ordering is (x - _roots[0]) * (x - _roots[1]) ...
  QUATERNION h = (point - bottomRoots[0]);
  for (unsigned int x = 1; x < bottomRoots.size(); x++)
    h = h * (point - bottomRoots[x]);
  
  QUATERNION gPrime = top.evaluateFactoredDerivative(point);
  QUATERNION hPrime = bottom.evaluateFactoredDerivative(point);
  
  const QUATERNION numerator = gPrime * h - g * hPrime;
  const QUATERNION inverse = h.inverse();
  const QUATERNION inverseSq = inverse * inverse;
  p = g * inverse;
  pPrime = numerator * inverseSq;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void POLYNOMIAL_4D::evaluateFactoredQuadratic(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, const QUATERNION& point, QUATERNION& p, QUATERNION& pPrime)
{
  // assume ordering is (x - _roots[0]) * (x - _roots[1]) ...
  const vector<QUATERNION> topRoots = top.roots();
  assert(topRoots.size() > 0);
  QUATERNION g = (point - topRoots[0]) * (point - topRoots[1]);

  const vector<QUATERNION> bottomRoots = bottom.roots();
  assert(bottomRoots.size() > 0);
  QUATERNION h = (point - bottomRoots[0]) * (point - bottomRoots[1]);
  
  QUATERNION gPrime = (point - topRoots[0]) + (point - topRoots[1]);
  QUATERNION hPrime = (point - bottomRoots[0]) + (point -  bottomRoots[1]);
  
  const QUATERNION numerator = gPrime * h - g * hPrime;
  const QUATERNION inverse = h.inverse();
  const QUATERNION inverseSq = inverse * inverse;
  p = g * inverse;
  pPrime = numerator * inverseSq;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
QUATERNION POLYNOMIAL_4D::evaluateFactoredDerivative(const QUATERNION& point) const
{
  QUATERNION final;

  //  which one to knock out?
  for (int x = 0; x < _totalRoots; x++)
  {
    QUATERNION knockedOut(1,0,0,0);
    for (int y = 0; y < _totalRoots; y++)
    {
      if (y == x) continue;
      knockedOut = knockedOut * (point - _roots[y]);
    }

    final += knockedOut;
  }

  return final;
}

//////////////////////////////////////////////////////////////////////
// sum of all the root coefficients
//////////////////////////////////////////////////////////////////////
Real POLYNOMIAL_4D::rootSum() const
{
  Real final = 0;
  for (int x = 0; x < _totalRoots; x++)
  {
    final += _roots[x].x();
    final += _roots[x].y();
    final += _roots[x].z();
    final += _roots[x].w();
  }
  return final; 
}

//////////////////////////////////////////////////////////////////////
// sum of all the root coefficients
//////////////////////////////////////////////////////////////////////
Real POLYNOMIAL_4D::powerSum() const
{
  Real result = 0;
  for (int x = 0; x < _totalRoots; x++)
    result += _rootPowers[x];
  return result; 
}

//////////////////////////////////////////////////////////////////////
// sum of all the root coefficients
//////////////////////////////////////////////////////////////////////
Real POLYNOMIAL_4D::scaledPowerSum() const
{
  Real result = 0;
  for (int x = 0; x < _totalRoots; x++)
    result += _powerScalar * _rootPowers[x];
  return result; 
}

//////////////////////////////////////////////////////////////////////
// overloaded operators
//////////////////////////////////////////////////////////////////////
POLYNOMIAL_4D& POLYNOMIAL_4D::operator*=(const Real& alpha)
{
  for (int x = 0; x < _totalRoots; x++)
    _roots[x] *= alpha;

  return *this;
}

//////////////////////////////////////////////////////////////////////
// overloaded operators
//////////////////////////////////////////////////////////////////////
POLYNOMIAL_4D& POLYNOMIAL_4D::operator-=(const VEC3F& v)
{
  QUATERNION q(v[0], v[1], v[2], 0);

  for (int x = 0; x < _totalRoots; x++)
    _roots[x] -= q;

  return *this;
}

//////////////////////////////////////////////////////////////////////
// overloaded operators
//////////////////////////////////////////////////////////////////////
POLYNOMIAL_4D& POLYNOMIAL_4D::operator+=(const VEC3F& v)
{
  QUATERNION q(v[0], v[1], v[2], 0);

  for (int x = 0; x < _totalRoots; x++)
    _roots[x] += q;

  return *this;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void POLYNOMIAL_4D::translateExceptFirst(const VEC3F& v)
{
  QUATERNION q(v[0], v[1], v[2], 0);

  for (int x = 1; x < _totalRoots; x++)
    _roots[x] += q;
}

//////////////////////////////////////////////////////////////////////
// test out taking the derivative of a power
//////////////////////////////////////////////////////////////////////
void POLYNOMIAL_4D::testSingleDerivative()
{
  //int res = 100;
  vector<QUATERNION> roots;
  roots.push_back(QUATERNION(1.0, 2.0, 3.0, 0));
  POLYNOMIAL_4D polynomial(roots);

  vector<Real> xs;
  vector<Real> ys;
  vector<Real> squares;
  vector<Real> gradients;
  vector<Real> temps;

  QUATERNION root(2, 2,3,4);
  //for (Real x = -10; x < 10; x += 0.25)
  //for (Real x = -10; x < 10; x += 0.1)
  for (Real x = -10; x < 10; x += 0.1)
  {
    QUATERNION value = root.pow(x);
    QUATERNION derivative = value * root.log();

    Real gradient = (1.0 / sqrt(value.dot(value))) * value.dot(derivative);
    assert(!(value.anyNans()) && !(derivative.anyNans()));

    xs.push_back(x);
    ys.push_back(sqrt(value.dot(value)));
    gradients.push_back(gradient);
  }

  // try to do some convergence
  Real dx = 0.1;
  Real fixedPoint = 2.0;

  QUATERNION value = root.pow(fixedPoint);
  QUATERNION derivative = value * root.log();

  Real gradient = pow(value.dot(value), (Real)-0.5) * value.dot(derivative);
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " Computed analytic: " << gradient << endl;

  vector<Real> dxs;
  vector<Real> diffs;

  for (int x = 0; x < 10; x++)
  {
    Real R = sqrt(value.dot(value));

    QUATERNION displaced = root.pow(fixedPoint + dx);
    Real Rdx = sqrt(displaced.dot(displaced));
    
    Real numerical = (Rdx - R) / dx;
    dx *= 0.1;

    dxs.push_back(dx);
    diffs.push_back((numerical - gradient) / gradient);

    cout << " diff: " << (numerical - gradient) / gradient << "\t dx: " << dx << "\t numerical: " << numerical << endl;
  }

  return;
}

//////////////////////////////////////////////////////////////////////
// take the derivative with respect to a root
//////////////////////////////////////////////////////////////////////
QUATERNION POLYNOMIAL_4D::powerDerivative(const QUATERNION& point, const int whichRoot) const
{
  TIMER functionTimer(__FUNCTION__);
  assert(whichRoot >= 0);
  assert(whichRoot < _totalRoots);

  QUATERNION left(1,0,0,0);
  QUATERNION right(1,0,0,0);
  
  for (int x = 0; x < whichRoot; x++)
    left *= (point - _roots[x]).pow(_rootPowers[x]);
  
  for (int x = whichRoot + 1; x < _totalRoots; x++)
    right *= (point - _roots[x]).pow(_rootPowers[x]);

  // get the derivative of the root in question
  QUATERNION base = (point - _roots[whichRoot]);
  QUATERNION power = base.pow(_rootPowers[whichRoot]);
  QUATERNION middle = power * base.log();

  return left * middle * right;
}

//////////////////////////////////////////////////////////////////////
// take the derivative with respect to a root
//////////////////////////////////////////////////////////////////////
QUATERNION POLYNOMIAL_4D::inversePowerDerivative(const QUATERNION& point, const int whichRoot) const
{
  TIMER functionTimer(__FUNCTION__);
  assert(whichRoot >= 0);
  assert(whichRoot < _totalRoots);

  QUATERNION left(1,0,0,0);
  QUATERNION right(1,0,0,0);
  
  for (int x = 0; x < whichRoot; x++)
    left *= (point - _roots[x]).pow(_rootPowers[x]);
  
  for (int x = whichRoot + 1; x < _totalRoots; x++)
    right *= (point - _roots[x]).pow(_rootPowers[x]);

  // get the derivative of the root in question
  QUATERNION base = (point - _roots[whichRoot]);
  QUATERNION power = base.pow(_rootPowers[whichRoot]);
  QUATERNION middle = power * base.log().inverse();

  return QUATERNION(-1,0,0,0) * left * middle * right;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void POLYNOMIAL_4D::testBulkDerivative()
{
  cout << " ========================================================================= " << endl;
  cout << "  TESTING DERIVATIVE OF ALL POWERS AT ONCE " << endl;
  cout << " ========================================================================= " << endl;
  vector<QUATERNION> roots;
  roots.push_back(QUATERNION(1.0, 2.0, 3.0, 4.0));
  roots.push_back(QUATERNION(4.0, 5.0, 6.0, 7.0));
  roots.push_back(QUATERNION(0.8, 0.9, 0.10, 0.11));

  //Real fixedPower = 3;
  vector<Real> powers;
  powers.push_back(1);
  powers.push_back(1);
  powers.push_back(1);

  POLYNOMIAL_4D polynomial(roots, powers);

  //QUATERNION fixedPoint(1,1,1,1);
  QUATERNION fixedPoint(-1,0.3,-0.1234,2);

  QUATERNION value = fixedPoint * polynomial.evaluatePowerFactored(fixedPoint);

  // analytical derivative
  Real totalGradient = 0;
  for (unsigned int x = 0; x < powers.size(); x++)
  {
    QUATERNION derivative = fixedPoint * polynomial.powerDerivative(fixedPoint, x);
    Real gradient = pow(value.dot(value), (Real)-0.5) * value.dot(derivative);

    totalGradient += gradient;
  }

  cout << " Analytical gradient:   " << totalGradient << endl;

  // does the numerical converge?
  Real dx = 0.1;

  vector<Real>& rootPowers = polynomial.powersMutable();

  Real R = sqrt(value.dot(value));
  for (int x = 0; x < 10; x++)
  {
    // do the dx displacement
    vector<Real> originalPowers;
    for (unsigned int y = 0; y < powers.size(); y++)
    {
      originalPowers.push_back(rootPowers[y]);
      rootPowers[y] += dx;
    }
    QUATERNION displaced = fixedPoint * polynomial.evaluatePowerFactored(fixedPoint);

    // wind things back
    for (unsigned int y = 0; y < powers.size(); y++)
      rootPowers[y] = originalPowers[y];

    Real Rdx = sqrt(displaced.dot(displaced));
    Real numerical = (Rdx - R) / dx;
    dx *= 0.1;

    cout << " diff: " << (numerical - totalGradient) / totalGradient << "\t dx: " << dx << "\t numerical: " << numerical << endl;
  }

  return;
}

//////////////////////////////////////////////////////////////////////
// the average root position
//////////////////////////////////////////////////////////////////////
QUATERNION POLYNOMIAL_4D::meanRootPosition() const
{
  QUATERNION result = _roots[0];
  for (unsigned int x = 1; x < _roots.size(); x++)
    result += _roots[x];

  return result * (1.0 / _roots.size());
}

//////////////////////////////////////////////////////////////////////
// stomp everything
//////////////////////////////////////////////////////////////////////
void POLYNOMIAL_4D::clear()
{
  _totalRoots = 0;
  _coeffs.clear();
  _derivs.clear();
  _secondDerivs.clear();
  _roots.clear();
  _rootPowers.clear();
  _ws.clear();
}
