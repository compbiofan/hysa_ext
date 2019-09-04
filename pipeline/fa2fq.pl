#!/usr/bin/perl
use warnings;
use strict;

require fa;

if($#ARGV == -1){
    die "This converts an fa file to fq file, with quals all I.\nUsage:fa2fq.pl <fa_file> > <fq_file>\n"; 
}

fa::fa2fq(@ARGV);
