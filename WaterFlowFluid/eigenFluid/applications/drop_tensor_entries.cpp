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
#include <fstream>
//#include <glog/logging.h>
#include <omp.h>
#include <string>
#include <sstream>
#include <math.h>
#include <memory>
#include <float.h>

#include "3D/laplacian_basis_set_3d.h"
#include "3D/dirichlet_basis_set_3d.h"
#include "3D/one_neumann_basis_set_3d.h"
#include "3D/two_neumann_x_3d_basis_set.h"
#include "3D/four_neumann_basis_set_3d.h"
#include "3D/six_neumann_basis_set_3d.h"
#include "setting.h"
#include "util/util.h"
#include "util/read_write_tensor.h"
#include "util/timer.h"

namespace {
// Compress the tensor by dropping some small entries.
void DropTensorEntries ( std::vector<Adv_Tensor_Type> * Adv_tensor) {
  const uint basis_dim_ = Adv_tensor->size();
  // First get the maximum and minimum value.
  double absMin = FLT_MAX;
  double absMax = -FLT_MAX;
  int64_t non_zeros_before = 0;
  
  for (int k = 0; k < basis_dim_; k++) {
    const Adv_Tensor_Type& mat = (*Adv_tensor)[k];
    for (int j = 0; j < mat.outerSize(); j++) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(mat,j); it; ++it) {
        const double Cijk = it.value();
        if (Cijk != 0) {
            const double absval = std::abs(Cijk);
          if (absval < absMin )
            absMin = absval;
          if (absval > absMax)
            absMax = absval;
          non_zeros_before ++;
        }
      }
    }
  }

  
  std::cout <<  "Min entry: " << absMin << std::endl;
  std::cout <<  "Max entry: " << absMax << std::endl;
  std::cout <<  "Non-zeros before dropping: " << non_zeros_before << std::endl;
  // Drop the tensor entries. Let's first drop five percent.
  const double dropping_threshold = absMax*0.47;
  
  for (int k = 0; k < basis_dim_; k++) {
    const Adv_Tensor_Type& mat = (*Adv_tensor)[k];
    for (int j = 0; j < mat.outerSize(); j++) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(mat,j); it; ++it) {
        const double Cijk = it.value();
        if (Cijk != 0) {
          const double absval = std::abs(Cijk);
          if (absval < dropping_threshold) {
            it.valueRef() = 0;
          }
        }
      }
    }
  }
  
  int64_t non_zeros_before_after = 0;
  for (int i = 0; i < basis_dim_; i++) {
    (*Adv_tensor)[i].makeCompressed();
    (*Adv_tensor)[i].prune(0, 0);
    non_zeros_before_after += (*Adv_tensor)[i].nonZeros();
  }
  std::cout <<  "Non-zeros after dropping: " << non_zeros_before_after << std::endl;
  std::cout <<  "Entries left: " << static_cast<float>(non_zeros_before_after) / non_zeros_before << std::endl;
}

void GetPerSliceMaxMin(const std::vector<Adv_Tensor_Type> & Adv_tensor) {
  const int basis_dim_ = Adv_tensor.size();

  for (int k = 0; k < basis_dim_; k++) {
    const Adv_Tensor_Type& mat = Adv_tensor[k];
    double absMin = FLT_MAX;
    double absMax = -FLT_MAX;
    for (int j = 0; j < mat.outerSize(); j++) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(mat,j); it; ++it) {
        const double Cijk = it.value();
        const double absval = std::abs(Cijk);
          if (absval < absMin )
            absMin = absval;
          if (absval > absMax)
            absMax = absval;
      }
    }
    
    std::cout <<  "slice: " << k << " Min: " << absMin << " Max: " << absMax << std::endl;
  }
}

void OutputHistogram(const Eigen::VectorXd& sigvals, const std::string& fname) {
  std::ofstream out(fname);
  if (!out.is_open()) {
    std::cout << "drop_tensor_entries.cpp " << __LINE__ << " FATAL: " <<  "Cannot open: " << fname << std::endl; exit(0);
  }
  
  int hist[11];
  for (int i = 0; i < 11; i++) hist[i] = 0;
  const double step_size = sigvals[0] / 10.0;
  int num_zeros = 0;
  for (int i = 0; i < sigvals.size(); i++) {
    const int index = static_cast<int> (sigvals[i] / step_size);
    hist[index] ++;
    if (std::abs(sigvals[i]) < 1e-10) {
      num_zeros ++;
    }
    out << sigvals[i] << "\n";
  }
  out.close();
  
  for (int i = 0; i < 11; i++) {
    std::cout <<  "Bin: " << i << " numberof eigvalues: " << hist[i] << std::endl;
  }
  std::cout <<  "Number of eig values close to zero: " << num_zeros << std::endl;
}

void ComputeTensorSliceSVD(const Adv_Tensor_Type& T, const std::string& fname) {
  
  // Convert sparse to dense.
  Eigen::MatrixXd dMat(T);
  
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(dMat);
  const Eigen::VectorXd sigvals = svd.singularValues();
  
  const double max_eigen = svd.singularValues()(0);
  const double min_eigen = svd.singularValues()(svd.singularValues().size() - 1);
  std::cout <<  "max eigenvalue: " << max_eigen << std::endl;
  std::cout <<  "min eigenvalue: " << min_eigen << std::endl;
  
  OutputHistogram(sigvals, fname);
  
}

}

int main(int argc, char ** argv){
//  google::ParseCommandLineFlags(&argc, &argv, true);
//  google::InitGoogleLogging(argv[0]);
  std::unique_ptr<LaplacianBasisSet3D> basis_set_;
  

  const std::string tensorName = "T3DType2Dim8162X128Y64Z64PX";
  const std::string basis_tensor_file_ = "./Tensor/Test_T/" + tensorName;
  const std::string basis_type_ = "two_neumann_x";

  const std::string constant_init_strategy_ = "principle_x";
  const std::string outTensorFile = "/media/qiaodong/Storage1/EigenFluid-Tensor/" + tensorName + "ZIPXX%";
  
  const int des_basis_dim_ = 2048;
  const int xRes_ = 64, yRes_ = 64, zRes_ = 64;
  
  int basis_typeint_;
  
  if (basis_type_ == "all_dirichlet") {
    basis_set_.reset(new DirichletBasisSet3D(des_basis_dim_, xRes_, yRes_, zRes_, 
                     constant_init_strategy_));
    basis_typeint_ = 0;
  } else if (basis_type_ == "one_neumann") {
    basis_set_.reset(new OneNeumannBasisSet3D(des_basis_dim_, xRes_, yRes_, zRes_, 
                     constant_init_strategy_));
    basis_typeint_ = 1;
  } else if (basis_type_ == "two_neumann_x") {
    basis_set_.reset(new TwoNeumannXBasisSet3D(des_basis_dim_, xRes_, yRes_, zRes_, 
                     constant_init_strategy_));
    basis_typeint_ = 2;
  } else if (basis_type_ == "four_neumann_xz") {
    basis_set_.reset(new FourNeumannBasisSet3D(des_basis_dim_, xRes_, yRes_, zRes_, 
                     constant_init_strategy_));
    basis_typeint_ = 4;
  } else if (basis_type_ == "six_neumann") {
    basis_set_.reset(new SixNeumannBasisSet3D(des_basis_dim_, xRes_, yRes_, zRes_, 
                     constant_init_strategy_));
    basis_typeint_ = 6;
  } else {
    std::cout << "drop_tensor_entries.cpp " << __LINE__ << " FATAL: " <<  "Unknown basis type: " << basis_type_ << std::endl; exit(0);
  }
  int basis_dim_;
  
  std::ifstream infile(basis_tensor_file_.c_str());
  if (!infile.is_open()) {
    std::cout << "drop_tensor_entries.cpp " << __LINE__ << " FATAL: " <<  "Cannot open tensor file: " << basis_tensor_file_ << std::endl; exit(0);
    return false;
  }
  basis_dim_ = basis_set_.get()->ReadBasis(infile);
  std::cout <<  "tensor dimension: " << basis_dim_ << std::endl;
  // Read tensor.
  std::vector<Adv_Tensor_Type> Adv_tensor_;
  ReadTensor(infile, basis_dim_, basis_typeint_, &Adv_tensor_);
  
  // ComputeTensorSliceSVD(Adv_tensor_[0], "./Eigs0.txt");
  // ComputeTensorSliceSVD(Adv_tensor_[100], "./Eigs100.txt");
  // ComputeTensorSliceSVD(Adv_tensor_[1600], "./Eigs1600.txt");
  
  
  DropTensorEntries(&Adv_tensor_);
  
  // Verify the tensor is still symmetric.
  // basis_set_.get()->VerifyAntisymmetric(Adv_tensor_);
  
  // std::stringstream file_name;
  // file_name << basis_tensor_file_ << "ZIPXX%";
  std::ofstream out(outTensorFile);
  if (!out.is_open()) {
    std::cout << "drop_tensor_entries.cpp " << __LINE__ << " FATAL: " <<  "Cannot write to file: " << outTensorFile << std::endl; exit(0);
  }
  
  basis_set_.get()->WriteBasis(out);
  WriteTensor(Adv_tensor_, basis_typeint_, out);
  out.close();
  
  return 0;
}
