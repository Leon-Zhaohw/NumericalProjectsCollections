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
#ifndef POLYNOMIAL_4D_H
#define POLYNOMIAL_4D_H

#include <cmath>
#include <vector>
#include "VEC3.h"
#include "QUATERNION.h"

using namespace std;

class POLYNOMIAL_4D {
public:
  POLYNOMIAL_4D();
  ~POLYNOMIAL_4D();

  // compute coefficients from the roots
  POLYNOMIAL_4D(const vector<QUATERNION>& roots);
  POLYNOMIAL_4D(const vector<QUATERNION>& roots, const vector<Real> powers);

  // set the coefficients directly - powers increase with index, i.e.
  // coeffs[0] -> x^0
  // coeffs[1] -> x^1
  POLYNOMIAL_4D(const vector<float>& coeffs);

  QUATERNION evaluate(const QUATERNION& point) const;
  QUATERNION evaluateDerivative(const QUATERNION& point) const;
  QUATERNION evaluateFactoredDerivative(const QUATERNION& point) const;
  QUATERNION evaluateSecondDerivative(const QUATERNION& point) const;
  void evaluateMultiple(const QUATERNION& point, QUATERNION& poly, QUATERNION& deriv) const;
  void evaluateMultiple(const QUATERNION& point, QUATERNION& poly, QUATERNION& deriv, QUATERNION& secondDeriv) const;

  // use the brute force nested formulation
  QUATERNION evaluateFactored(const QUATERNION& point) const;
  QUATERNION evaluatePowerFactored(const QUATERNION& point) const;
  QUATERNION evaluateScaledPowerFactored(const QUATERNION& point) const;
  QUATERNION evaluateScaledPowerFactoredVerbose(const QUATERNION& point) const;
  QUATERNION evaluateFactoredDouble(const QUATERNION& point) const;
  QUATERNION evaluateFactoredPositive(const QUATERNION& point) const;
  void computeNestedCoeffs();

  // use the brute force nested formulation, but cache the multiplies
  QUATERNION evaluatePowerFactored(const QUATERNION& point, vector<QUATERNION>& forward, vector<QUATERNION>& backward) const;
  QUATERNION evaluateScaledPowerFactored(const QUATERNION& point, vector<QUATERNION>& forward, vector<QUATERNION>& backward) const;
  QUATERNION evaluateRoot(const QUATERNION& point, const int index, const bool debug = false) const;

  // add a new root
  void addRoot(const QUATERNION& newRoot);
  void addRoot(const QUATERNION& newRoot, const Real& power);
  void addFrontRoot(const QUATERNION& newRoot);
  
  // modify an existing root
  void modifyRoot(const int whichRoot, const QUATERNION& newRoot);

  // accessors
  const int totalRoots() const { return _totalRoots; };
  const vector<QUATERNION>& roots() const { return _roots; };
  const vector<QUATERNION>& coeffs() const { return _coeffs; };
  const vector<QUATERNION>& derivs() const { return _derivs; };
  const vector<QUATERNION>& secondDerivs() const { return _secondDerivs; };
  const vector<Real>& powers() const { return _rootPowers; };
  vector<QUATERNION>& rootsMutable() { return _roots; };
  vector<Real>& powersMutable() { return _rootPowers; };
  Real& powerScalar() { return _powerScalar; };
  const Real powerScalar() const { return _powerScalar; };

  // run a unit test for a known rational function
  static void rationalTest();

  // evaluate rational function
  static void evaluateRational(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, const QUATERNION& point, QUATERNION& p, QUATERNION& pPrime, QUATERNION& pPrimePrime);
  static void evaluateRational(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, const QUATERNION& point, QUATERNION& p, QUATERNION& pPrime);
  
  static void evaluateFactoredRational(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, const QUATERNION& point, QUATERNION& p, QUATERNION& pPrime);

  // for debugging purposes
  static void evaluateFactoredQuadratic(const POLYNOMIAL_4D& top, const POLYNOMIAL_4D& bottom, const QUATERNION& point, QUATERNION& p, QUATERNION& pPrime);

  // resize the polynomial to a given number of roots
  void resizeAndWipe(int totalRoots);

  // file IO
  void write(FILE* file) const;
  void writeTextfile(FILE* file) const;
  void read(FILE* file);
  void readTextfile(FILE* file);

  // get the condition number of the polynomial
  Real conditionNumber();

  // sum of all the root coefficients
  Real rootSum() const;
  
  // sum of all the powers
  Real powerSum() const;
  Real scaledPowerSum() const;

  // overloaded operators
  POLYNOMIAL_4D& operator*=(const Real& alpha);
  POLYNOMIAL_4D& operator-=(const VEC3F& v);
  POLYNOMIAL_4D& operator+=(const VEC3F& v);

  void translateExceptFirst(const VEC3F& v);

  // change the power of a root
  void changePower(const int& whichRoot, const Real& newPower) {
    assert(whichRoot < _totalRoots);
    _rootPowers[whichRoot] = newPower;
  };

  // the average root position
  QUATERNION meanRootPosition() const;

  // take the derivative with respect to a root
  QUATERNION powerDerivative(const QUATERNION& point, const int whichRoot) const;
  QUATERNION inversePowerDerivative(const QUATERNION& point, const int whichRoot) const;

  // test out taking the derivative of a power
  static void testSingleDerivative();
  static void testBulkDerivative();

  // stomp everything
  void clear();

private:
  int _totalRoots;
  vector<QUATERNION> _coeffs;
  vector<QUATERNION> _derivs;
  vector<QUATERNION> _secondDerivs;
  vector<QUATERNION> _roots;

  vector<Real> _rootPowers;

  // the nested linear factorization
  vector<QUATERNION> _ws;

  // a uniform scalar to multiply all the powers by
  Real _powerScalar;

  // compute the polynomial coefficients
  void computeCoeffs();
 
  // compute the derivative coefficients
  void computeDerivativeCoeffs();

  // compute the polynomial coefficients
  void computeCoeffsFast();

  // support function for coefficients
  vector<QUATERNION> choose(const vector<QUATERNION> m, const int n);
};

ostream& operator<<(ostream &out, const POLYNOMIAL_4D& poly);

#endif
