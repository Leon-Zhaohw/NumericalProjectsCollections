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
// #include <glog/logging.h>
#include <sstream>
#include <math.h>

#include "2D/laplacian_basis_2D.h"
#include "2D/all_dirichlet_basis_2D.h"
#include "2D/three_dirichlet_one_neumann.h"
// #include "2D/two_neumann_x.h"
#include "2D/two_neumann_x_2D.h"
#include "setting.h"
#include "util/util.h"
#include "util/read_write_tensor.h"
#include "util/timer.h"
#include "util/SIMPLE_PARSER.h"

/*
DEFINE_int32(basis_dim_root, 20, "The root square of the basis.");
DEFINE_string(basis_type, "all_dirichlet", "The basis type. Can be all_dirichlet,"
              "three_dirichlet_one_neumann, two_neumann_x");
DEFINE_string(folder_name,"", "The folder to put the tensor, with slash.");
*/

// Type number for all_dirichlet is 0, three_dirichlet_one_neumann is 1, two_neumann_x is 2.

int main(int argc, char ** argv){
  // google::ParseCommandLineFlags(&argc, &argv, true);
  // google::InitGoogleLogging(argv[0]);
  
  if (argc != 2)
  {
    std::cout << " Usage: " << argv[0] << " *.cfg" << std::endl;
    return 0;
  }
  SIMPLE_PARSER parser(argv[1]);
  
  int FLAGS_basis_dim_root = parser.getInt("basis_dim_root", 20);
  std::string  FLAGS_basis_type = parser.getString("basis_type", "all_dirichlet");
  std::string FLAGS_folder_name = parser.getString("folder_name", "");
  
  const int basis_dim = FLAGS_basis_dim_root*FLAGS_basis_dim_root;
  
  std::cout << "Precompute the tensor with dimention: " << basis_dim << " Type: "
      << FLAGS_basis_type;
      
  int basis_type;
  std::unique_ptr<LaplacianBasis2D> basis_;
  if (FLAGS_basis_type == "all_dirichlet") {
    basis_.reset(new AllDirichletBasis2D(128, FLAGS_basis_dim_root));
    basis_type = 0;
  } else if (FLAGS_basis_type == "three_dirichlet_one_neumann") {
    basis_.reset(new ThreeDirichletOneNeumann(128, FLAGS_basis_dim_root));
    basis_type = 1;
  } else if (FLAGS_basis_type == "two_neumann_x") {
    basis_.reset(new TwoNeumannX2D(128, FLAGS_basis_dim_root));
    basis_type = 2;
  } else {
    std::cout << "precompute_2D_tensor.cpp " << __LINE__ << " FATAL: " <<  "Unknow basis type." << std::endl; exit(0);
  }
  std::stringstream file_name;
  file_name << FLAGS_folder_name << "T2D" << "Type" << basis_type << "Dim" << basis_dim;
  std::cout <<  "Writting tensor to file: " << file_name.str() << std::endl;
  
  std::vector<Adv_Tensor_Type> Adv_tensor;
  basis_.get()->FillVariationalTensor(&Adv_tensor);
  // WriteTensor(Adv_tensor,basis_type ,file_name.str());
  WriteTensor(Adv_tensor, basis_type, file_name.str());
  
  return 0;
}
