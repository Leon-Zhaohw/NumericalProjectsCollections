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
#ifndef EIGEN_H
#define EIGEN_H
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

#include <string>
#include <zlib.h>
#include <VECTOR.h>
#include <vector>

using namespace std;
using namespace Eigen;

class EIGEN {
public:
  static bool read(const string& filename, MatrixXd& input);
  static bool readBig(const string& filename, MatrixXd& input);
  static void read(FILE* file, MatrixXd& input);
  static void readRaw(FILE* file, const int size, VectorXd& input);
  static void read(const char* filename, VectorXd& inpute);

  static bool write(const string& filename, const MatrixXd& input);
  static void write(FILE* file, const MatrixXd& input);
  static void write(FILE* file, const VectorXd& input);
  
  static void write(const char* filename, const MatrixXd& input);
  static void write(const char* filename, const VectorXd& input);

  static MatrixXd buildFromColumns(const vector<VectorXd>& columns);
  static MatrixXd getRows(const int rowBegin, const int totalRows, const MatrixXd& input);
  static void getRowsMemory(const int rowBegin, const int totalRows, FILE* matrixInfo, double* matrixData, int rows, int cols);
  static void transposeProduct(const MatrixXd& left, const MatrixXd& right, MatrixXd& output);
  
  static VECTOR convert(const VectorXd& input);
  static VECTOR convertInt(const VectorXi& input);
  static VectorXd convert(const VECTOR& input);
};

#endif
