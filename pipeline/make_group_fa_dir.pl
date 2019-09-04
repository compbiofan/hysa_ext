#!/usr/bin/perl
use warnings;
use strict;
require union_find;

if($#ARGV == -1){
    die "This takes a file with read names of both short and long reads, and make a fasta file from the original fastas and a hash of long. The long file name is x.y, in which x is a number denoting where the original file comes from, and y is the index in that file. The prefix should contain everything until x, and long_fa should be everything after x. If there are more than one resources of short read, concatenate the prefix by semicolon, and concatenate the string in the suffix of readnames and the prefix by colon. If not suffix in readname, then leave it blank and without colon. \nUsage:$0 <group> <short_fa_prefix> <long_prefix> <long_fa>\n";
}
union_find::readname2fa_dir(@ARGV);
