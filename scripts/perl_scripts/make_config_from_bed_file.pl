#!/usr/bin/perl
use warnings;
use strict;
require fa;

if(@ARGV == 0){
    die "Given the bed file, reference, flanking region size and the mode (support del and ins now), reconstruct the alternative allele and write to given fa file. \nUsage: $0 <bed_file> <ref> <out> <flank> <mode>\n";
}

my ($bed, $ref, $out, $flank, $mode) = @ARGV;
fa::config_fa_bed_file($bed, $ref, $out, $flank, $mode);
