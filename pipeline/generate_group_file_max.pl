#!/usr/bin/perl
use warnings;
use strict;
require file;

if($#ARGV == -1){
    die "Given the prefix of two files (out.all: containing the stats of all groups; out.txt: containing the ILs and PBs corresponding to a group, with the group number appearing before all following #), generate a file with all the reads whose group has > max pbs. \nUsage: $0 <prefix> <max> <out_file>\n";
}
my ($prefix, $threshold, $out) = @ARGV;
my $tmp = $out . ".tmp"; 
file::get_max($prefix, $threshold, $tmp);
file::generate_groups_file($prefix, $tmp, $out);
`rm $tmp`;
