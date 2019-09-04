#!/usr/bin/perl
use warnings;
use strict;

require m4;

if(@ARGV == 0){
    die "This scripts takes a m4 file, a file containing ILs and PBs, each in one line, and produces a new m4 file that contains a subset of the lines in the original m4 file, so that the IL and PB are in the given file simultaneously. \nUsage: $0 <m4_file> <new_m4_file> <out_txt_file>\n";
}

m4::get_new_m4_from_out(@ARGV);

