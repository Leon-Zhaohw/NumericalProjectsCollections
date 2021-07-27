#!/usr/bin/perl -w
use strict;
use warnings;
use Time::HiRes;

my $BIN_PATH = "bin/fluid_stam";
my $script_file = "./Dedalus.flags";

system($BIN_PATH." --flagfile=".$script_file );