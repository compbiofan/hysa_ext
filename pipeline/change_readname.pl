#!/usr/bin/perl
use warnings;
use strict;
require fa;

if($#ARGV == -1){
    die "Insert string at the pos (head/default:tail) of readname, overwrite the original fasta.\nUsage: $0 <fasta> <str> <pos>\n";
}
fa::change_readname(@ARGV);
