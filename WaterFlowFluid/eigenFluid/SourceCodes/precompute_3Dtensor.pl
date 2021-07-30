use strict;
use warnings;
use Time::HiRes;
my $BIN_PATH = "bin/precompute_3D_tensor";
#
#This will determine how many basis is allocated to each direction. The number
#of basis allocated to that direction is proportional to the resolution
#of that direction. The resolution specificed here does not need to match
#the resolution specificed in the scene cfg file. Only the ratio of xRes yRes and zRes matters here.

my $xRes = 64;
my $yRes = 64;
my $zRes = 64;
#'''
#basis_dim_des: The total number of basis desired. The actually computed basis would be smaller than that,
#so the number of basis in each direction would be exactly an integer. E.g. if xRes = yRes = zRes,
#and basis_dim_des = 30, then 27 basis will be allocated, because 3^3 = 27 < 30
#'''
my $basis_dim_des = 3000;
#'''
#Which type of basis.
#'''
#all_dirichlet, one_neumann, two_neumann_x, four_neumann_xz, six_neumann
my $basis_type = "all_dirichlet";
#principle_x, principle_y. principle_z, random, uniform
my $constant_init_strategy = "principle_x";
my $folder_name = "./Tensor/";

# convert to flags.
$basis_dim_des = "basis_dim_des=".$basis_dim_des;
$basis_type = "basis_type=".$basis_type;
$folder_name = "folder_name=".$folder_name;
$constant_init_strategy = "constant_init_strategy=".$constant_init_strategy;
$xRes = "xRes=".$xRes;
$yRes = "yRes=".$yRes;
$zRes = "zRes=".$zRes;

my $scriptname = 'precompute3DTensor.cfg';

open(my $fh, '>', $scriptname) or die "Could not open file '$scriptname' $!";
print $fh "$basis_dim_des\n$basis_type\n$folder_name\n$constant_init_strategy\n$xRes\n$yRes\n$zRes";
close $fh;

system("$BIN_PATH $scriptname");

system("rm ".$scriptname);