#!/usr/bin/perl
use warnings;
use strict;
require file;

if($#ARGV == -1){
    die "Given the prefix of two files (out.all: containing the stats of all groups; out.txt: containing the ILs and PBs corresponding to a group, with the group number appearing before all following #).\nUsage: $0 <prefix> <group> <out_file>\n";
}
file::generate_group_file(@ARGV);
