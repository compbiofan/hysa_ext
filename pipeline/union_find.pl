#!/usr/bin/perl
use warnings;
use strict;

require union_find;

if($#ARGV == -1){
    die "This uses union find to connect components in a m4 file, output with the prefix the files of read names that are in the same connected component. Starting from 09042015, output to output_prefix.txt, with a stat file in output_prefix.all. Old version is in bitbucket.\nUsage:union_find.pl <m4_file> <output_prefix>\n"; 
}

union_find::union_find(@ARGV);
