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
 
 Copyright 2018 Qiaodong Cui (qiaodong@ucsb.edu)
 */

#include "Eigen"
// #include <unsupported/Eigen/SparseExtra>
#include <fstream>
// #include <glog/logging.h>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>

#include "util/util.h"
#include "util/read_write_tensor.h"
#include "util/timer.h"
#include "util/stringprintf.h"

void ReadTensor(std::ifstream& infile, const int expected_dim ,const int expected_type, 
                std::vector<Adv_Tensor_Type>* Adv_tensor) {
  if (!infile.is_open()) {
    std::cout << "read_write_tensor.cpp " << __LINE__ << " FATAL: " <<  "Cannot read in stream: " << std::endl; exit(0);
  }
  int tensor_dim;
  int tensor_type;
  // First read in the tensor dimention.
  infile.read(reinterpret_cast<char*>(&tensor_dim), sizeof(int));
  infile.read(reinterpret_cast<char*>(&tensor_type), sizeof(int));
  if (! tensor_dim == expected_dim) {std::cout << "read_write_tensor.cpp " << __LINE__ << " FATAL: "  << "Read the tensor with wrong dimention: " << tensor_dim << " The expected dimention: " << expected_dim << std::endl; exit(0);}
  if (! tensor_type == expected_type) {std::cout << "read_write_tensor.cpp " << __LINE__ << " FATAL: "  << "Read the tensor with wrong type: " << tensor_type << " The expected dimention: " << expected_type << std::endl; exit(0);}
  if (Adv_tensor != nullptr ){ Adv_tensor->clear();}
  Adv_tensor->reserve(tensor_dim);
  std::cout <<  "start to read..." << std::endl;
  // Read each matrix.
  for (int i = 0; i < tensor_dim; i++) {
    uint64_t nnzs, outS,innS;
    infile.read(reinterpret_cast<char*>(&nnzs), sizeof(uint64_t));
    infile.read(reinterpret_cast<char*>(&outS), sizeof(uint64_t));
    infile.read(reinterpret_cast<char*>(&innS), sizeof(uint64_t));
    Adv_Tensor_Type Ck(tensor_dim, tensor_dim);
    Ck.makeCompressed();
    Ck.resizeNonZeros(nnzs);
    infile.read(reinterpret_cast<char*>(Ck.valuePtr()), sizeof(double)*nnzs);
    infile.read(reinterpret_cast<char*>(Ck.outerIndexPtr()), sizeof(int)*outS);
    infile.read(reinterpret_cast<char*>(Ck.innerIndexPtr()), sizeof(int)*nnzs);
    Ck.finalize();
    Adv_tensor->emplace_back(Ck);
  }
}

void ReadTensor(const std::string file, const int expected_dim ,const int expected_type, 
                   std::vector<Adv_Tensor_Type>* Adv_tensor) {
  std::ifstream infile(file.c_str());
  if (!infile.is_open()) {
    std::cout << "read_write_tensor.cpp " << __LINE__ << " FATAL: " <<  "Cannot open tensor file: " << file << std::endl; exit(0);
  }
  ReadTensor(infile, expected_dim, expected_type, Adv_tensor);
  std::cout <<  "Read completed.." << std::endl;
  infile.close();
}

void WriteTensor(const std::vector<Adv_Tensor_Type>& Adv_tensor,
                 const int tensor_type, std::ofstream& out) {
  if (!out.is_open()) {
    std::cout << "read_write_tensor.cpp " << __LINE__ << " FATAL: " <<  "Cannot write to stream:" << std::endl; exit(0);
  }
  int tensor_dim = Adv_tensor.size();
  int type_ = tensor_type;
  
  // first write out the tensor dimention.
  out.write(reinterpret_cast<const char *>(&tensor_dim), sizeof(int));
  // Then write out the tensor type. encoded as int.
  out.write(reinterpret_cast<const char *>(&type_), sizeof(int));
  // Write each sparse matrix.
  for (int i = 0; i < tensor_dim; i++) {
    uint64_t nnzs = Adv_tensor[i].nonZeros();
    uint64_t outS = Adv_tensor[i].outerSize();
    uint64_t innS = Adv_tensor[i].innerSize();
    out.write(reinterpret_cast<const char *>(&nnzs), sizeof(uint64_t));
    out.write(reinterpret_cast<const char *>(&outS), sizeof(uint64_t));
    out.write(reinterpret_cast<const char *>(&innS), sizeof(uint64_t));
    // Write the matrix values.
    out.write(reinterpret_cast<const char *>(Adv_tensor[i].valuePtr()), sizeof(double)*nnzs);
    // Output the index values, the default index value type for eigen is int.
    out.write(reinterpret_cast<const char *>(Adv_tensor[i].outerIndexPtr()), sizeof(int)*outS);
    out.write(reinterpret_cast<const char *>(Adv_tensor[i].innerIndexPtr()), sizeof(int)*nnzs);
  }
}


void WriteTensor(const std::vector<Adv_Tensor_Type>& Adv_tensor,
                 const int tensor_type, const std::string file) {
  
  std::ofstream out(file);
  if (!out.is_open()) {
    std::cout << "read_write_tensor.cpp " << __LINE__ << " FATAL: " <<  "Cannot write to file: " << file << std::endl; exit(0);
  }
  WriteTensor(Adv_tensor, tensor_type, out);
  out.close();
}

void ReadTensorFromMatlab(const std::string& file_name, const int basis_dim, 
std::vector<Adv_Tensor_Type>* new_tensor) {
  
  new_tensor->reserve(basis_dim);
  for (int k = 0; k < basis_dim; k++) {
    Adv_Tensor_Type Ck(basis_dim, basis_dim);
    Ck.setZero();
    new_tensor->emplace_back(Ck);
  }
  
  const int N = basis_dim*basis_dim*basis_dim;
  double *U = (double*)malloc(N * sizeof(double));
  std::ifstream infile;
  infile.open(file_name.c_str(), std::ios::in | std::ios::binary);
  infile.read((char*)U, N*sizeof(double));
  infile.close();
  for (int k = 0; k < basis_dim; k++) {
    for (int j = 0; j < basis_dim; j++) {
      for (int i = 0; i < basis_dim; i++) {
        int index = k * basis_dim * basis_dim + j * basis_dim + i;
        if (std::abs(U[index]) > 1e-10) {
          SetMatrixEntry(&(*new_tensor)[k], i, j, U[index]);
        }
      }
    }
  }
  delete [] U;
  infile.close();
}

