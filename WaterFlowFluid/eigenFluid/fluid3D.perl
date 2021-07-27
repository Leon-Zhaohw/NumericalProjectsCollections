#!/usr/bin/perl -w
use strict;
use warnings;
use Time::HiRes;

my $BIN_PATH = "bin/laplacian_fluid_sim3D";

my $script_file = "./fluid_moving_sphere.cfg";
#my $script_file = "./fluid_3D_paddle.cfg";
#my $script_file = "./fluid_two_phase_pulse.cfg";

system("$BIN_PATH $script_file");