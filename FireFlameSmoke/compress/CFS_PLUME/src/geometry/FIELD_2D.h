#ifndef FIELD_2D_H
#define FIELD_2D_H

#include <cmath>
#include <string>
#include <iostream>
//#include "VEC3F.h"
#include "VEC3.h"

using namespace std;

class FIELD_2D {
public:
  FIELD_2D();
  FIELD_2D(const int& rows, const int& cols);
  FIELD_2D(const FIELD_2D& m);
  ~FIELD_2D();

  // accessors
  inline float& operator()(int x, int y) { return _data[y * _xRes + x]; };
  const float operator()(int x, int y) const { return _data[y * _xRes + x]; };
  inline float& operator[](int x) { return _data[x]; };
  const float operator[](int x) const { return _data[x]; };
  float* data() { return _data; };
  const int xRes() const { return _xRes; };
  const int yRes() const { return _yRes; };
  const int totalCells() const { return _totalCells; };

  // common field operations
  void clear();
  void normalize();
  FIELD_2D& abs();

  float min();
  float max();

  // take the log
  void log(float base = 2.0);
 
  // generic IO functions
  void writeMatlab(string filename, string variableName) const;
  void write(string filename) const;
  void read(string filename);

  // some image file support
  void writePPM(string filename);
  void writeJPG(string filename);
  void writePNG(string filename);
  void readPNG(string filename);

  void resizeAndWipe(int xRes, int yRes);

  // overloaded operators
  FIELD_2D& operator=(const float& alpha);
  FIELD_2D& operator=(const FIELD_2D& A);
  FIELD_2D& operator*=(const float& alpha);
  FIELD_2D& operator/=(const float& alpha);
  FIELD_2D& operator+=(const float& alpha);
  FIELD_2D& operator-=(const FIELD_2D& input);
  FIELD_2D& operator+=(const FIELD_2D& input);
  FIELD_2D& operator*=(const FIELD_2D& input);
  FIELD_2D& operator/=(const FIELD_2D& input);

  // sum of all entries
  float sum();
  
  // set to a checkboard for debugging
  void setToCheckerboard(int xChecks = 10, int yChecks = 10);
  
private:
  int _xRes;
  int _yRes;
  int _totalCells;
  float* _data;
};

FIELD_2D operator*(const FIELD_2D& A, const float alpha);
FIELD_2D operator/(const FIELD_2D& A, const float alpha);
FIELD_2D operator+(const FIELD_2D& A, const float alpha);
FIELD_2D operator*(const float alpha, const FIELD_2D& A);
FIELD_2D operator+(const float alpha, const FIELD_2D& A);
FIELD_2D operator-(const FIELD_2D& A, const FIELD_2D& B);
FIELD_2D operator+(const FIELD_2D& A, const FIELD_2D& B);

#endif
