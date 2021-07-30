#!/usr/bin/perl -w
use strict;
use warnings;
use Time::HiRes;

my $BIN_PATH = "bin/fluid_DCT_sim3D";
my $script_file = "./ted_cfgs/fluid_DCT.cfg";

system($BIN_PATH." ".$script_file." --logtostderr");