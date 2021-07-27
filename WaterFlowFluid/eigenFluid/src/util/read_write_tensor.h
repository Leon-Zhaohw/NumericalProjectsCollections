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

#ifndef READ_WRITE_TENSOR_H
#define READ_WRITE_TENSOR_H

#include <fstream>
#include "Eigen"
#include <vector>
#include <string>

#include "setting.h"

void ReadTensorFromMatlab(const std::string& file_name, const int basis_dim, 
                          std::vector<Adv_Tensor_Type>* new_tensor);

void ReadTensor(const std::string file, const int expected_dim ,const int expected_type, 
                std::vector<Adv_Tensor_Type>* Adv_tensor);

void ReadTensor(std::ifstream& infile, const int expected_dim ,const int expected_type, 
                std::vector<Adv_Tensor_Type>* Adv_tensor);

void WriteTensor(const std::vector<Adv_Tensor_Type>& Adv_tensor,
                 const int tensor_type, std::ofstream& out);

void WriteTensor(const std::vector<Adv_Tensor_Type>& Adv_tensor,
                    const int tensor_type, const std::string file);

#endif  // READ_WRITE_TENSOR_H
