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
#include <string>
#include <sstream>
#include <math.h>
#include <memory>

#include "3D/dirichlet_basis_set_3d.h"
#include "3D/laplacian_basis_set_3d.h"
#include "3D/one_neumann_basis_set_3d.h"
#include "3D/two_neumann_x_3d_basis_set.h"
#include "3D/four_neumann_basis_set_3d.h"
#include "3D/six_neumann_basis_set_3d.h"

#include "setting.h"
#include "util/util.h"
#include "util/read_write_tensor.h"
#include "util/timer.h"
#include "util/SIMPLE_PARSER.h"
/*
DEFINE_int32(basis_dim_des, 20, "The number of basis desired.");
DEFINE_string(basis_type, "all_dirichlet", "The basis type. Can be all_dirichlet,"
              "one_neumann, two_neumann_x");
DEFINE_string(folder_name,"", "The folder to put the tensor, with slash.");
DEFINE_int32(xRes, 64, "The x resolution of fluid.");
DEFINE_int32(yRes, 64, "The y resolution of fluid.");
DEFINE_int32(zRes, 64, "The z resolution of fluid.");

DEFINE_string(constant_init_strategy, "principle_x", "Choose which direction as principle"
              " propagation direction. Can be principle_x, principle_y, principle_z, random.");
*/

int main(int argc, char ** argv){
  // google::ParseCommandLineFlags(&argc, &argv, true);
  // google::InitGoogleLogging(argv[0]);
  
  // read in the cfg file
  if (argc != 2)
  {
    cout << " Usage: " << argv[0] << " *.cfg" << endl;
    return 0;
  }
  SIMPLE_PARSER parser(argv[1]);

  int FLAGS_xRes = parser.getInt("xRes", 48);
  int FLAGS_yRes = parser.getInt("yRes", 64);
  int FLAGS_zRes = parser.getInt("zRes", 48);
  std::string FLAGS_folder_name = parser.getString("folder_name", "");
  std::string FLAGS_constant_init_strategy = parser.getString("constant_init_strategy", "principle_x");
  std::string FLAGS_basis_type = parser.getString("basis_type", "all_dirichlet");
  int FLAGS_basis_dim_des = parser.getInt("basis_dim_des", 20);
  
  int basis_type = 0;
  std::cout <<  "xRes: " << FLAGS_xRes << " yRes : " << FLAGS_yRes << " zRes: " << FLAGS_zRes << std::endl;
  
  std::unique_ptr<LaplacianBasisSet3D> basis_set;
  if (FLAGS_basis_type == "all_dirichlet") {
    basis_set.reset(new DirichletBasisSet3D(FLAGS_basis_dim_des,
        FLAGS_xRes,FLAGS_yRes, FLAGS_zRes, FLAGS_constant_init_strategy));

    basis_type = 0;
  } else if (FLAGS_basis_type == "one_neumann") {
    basis_set.reset(new OneNeumannBasisSet3D(FLAGS_basis_dim_des,
        FLAGS_xRes,FLAGS_yRes, FLAGS_zRes, FLAGS_constant_init_strategy));
    
    basis_type = 1;
  } else if (FLAGS_basis_type == "two_neumann_x") {
    basis_set.reset(new TwoNeumannXBasisSet3D(FLAGS_basis_dim_des,
        FLAGS_xRes,FLAGS_yRes, FLAGS_zRes, FLAGS_constant_init_strategy));
    
    basis_type = 2;
  } else if (FLAGS_basis_type == "four_neumann_xz") {
    basis_set.reset(new FourNeumannBasisSet3D(FLAGS_basis_dim_des,
        FLAGS_xRes,FLAGS_yRes, FLAGS_zRes, FLAGS_constant_init_strategy));
    
    basis_type = 4;
  } else if (FLAGS_basis_type == "six_neumann") {
    basis_set.reset(new SixNeumannBasisSet3D(FLAGS_basis_dim_des,
        FLAGS_xRes,FLAGS_yRes, FLAGS_zRes, FLAGS_constant_init_strategy));
    basis_type = 6;
  } else {
    std::cout << "precompute_3D_tensor.cpp " << __LINE__ << " FATAL: " <<  "Unknow basis type: " << FLAGS_basis_type << std::endl; exit(0);
  }
  std::string const_affix;
  if (FLAGS_constant_init_strategy == "principle_x") {
    const_affix = "PX";
  } else if (FLAGS_constant_init_strategy == "principle_y") {
    const_affix = "PY";
  } else if (FLAGS_constant_init_strategy == "principle_z") {
    const_affix = "PZ";
  } else if (FLAGS_constant_init_strategy == "random") {
    const_affix = "Rnd";
  } else if (FLAGS_constant_init_strategy == "uniform") {
    const_affix = "uniform";
  } else {
    std::cout << "precompute_3D_tensor.cpp " << __LINE__ << " FATAL: " <<  "Unknow constant_init_strategy: " << FLAGS_constant_init_strategy << std::endl; exit(0);
  }
  
  int basis_allocated = basis_set.get()->AllocateAllBasis();
  std::cout <<  "basis allocated actually: " << basis_allocated << std::endl;
  std::vector<Adv_Tensor_Type> Adv_tensor;
  std::cout <<  "Filling adv tensor." << std::endl;
  Timer timer;
  timer.Reset();
  basis_set.get()->FillVariationalTensor(&Adv_tensor);
  std::cout <<  "Time to compute the tensor: " << timer.ElapsedTimeInSeconds() << std::endl;
  
  std::stringstream file_name;
  file_name << FLAGS_folder_name << "T3D" << "Type" << basis_type << "Dim" << basis_allocated
   << "X" << FLAGS_xRes << "Y" << FLAGS_yRes << "Z" << FLAGS_zRes << const_affix;
  std::cout <<  "Writting tensor to file: " << file_name.str() << std::endl;
  
  std::ofstream out(file_name.str());
  if (!out.is_open()) {
    std::cout << "precompute_3D_tensor.cpp " << __LINE__ << " FATAL: " <<  "Cannot write to file: " << file_name.str() << std::endl; exit(0);
  }
  
  basis_set.get()->WriteBasis(out);
  WriteTensor(Adv_tensor, basis_type, out);
  out.close();
  
  return 0;
}
