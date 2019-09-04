#!/usr/bin/perl
use warnings;
use strict;
require fa;

if($#ARGV == -1){
    die "This converts an fq file to fa. Prefix is the whole path + fq filename excluding suffix such as fq or fastq. Print to prefix.fa. \nUsage:$0 <fq_file_prefix> <suffix_of_readname>\n";
}
fa::fq_to_fa(@ARGV);
