#!/usr/bin/perl -w
use strict;
use warnings;
use Time::HiRes;

my $BIN_PATH = "bin/laplacian_fluid_sim2D";
my $script_file = "./fluid2D.cfg";

system("$BIN_PATH $script_file");