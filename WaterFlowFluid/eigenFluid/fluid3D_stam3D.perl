#!/usr/bin/perl -w
use strict;
use warnings;
use Time::HiRes;

my $BIN_PATH = "bin/fluidStam3D";
my $script_file = "./ted_cfgs/paddle.cfg";

system($BIN_PATH." ".$script_file." --logtostderr");