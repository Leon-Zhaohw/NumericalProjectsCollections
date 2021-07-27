use strict;
use warnings;
use Time::HiRes;
my $BIN_PATH = "bin/precompute_2D_tensor";

#
my $basis_dim_root = 20;
#all_dirichlet, three_dirichlet_one_neumann, two_neumann_x
my $basis_type = "all_dirichlet";
my $folder_name = "./Tensor/";

# convert to flags.
$basis_dim_root = "basis_dim_root=".$basis_dim_root;
$basis_type = "basis_type=".$basis_type;
$folder_name = "folder_name=".$folder_name;

my $scriptname = 'precompute2DTensor.cfg';

open(my $fh, '>', $scriptname) or die "Could not open file '$scriptname' $!";
print $fh "$basis_dim_root\n$basis_type\n$folder_name";
close $fh;

system("$BIN_PATH $scriptname");

system("rm ".$scriptname);