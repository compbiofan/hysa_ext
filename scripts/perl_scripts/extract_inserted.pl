#!/usr/bin/perl
use warnings;
use strict;
require fa;

if(scalar(@ARGV) == 0){
    die "This extracted a subsequence from a fa file, with the corresponding name and positions in the txt file, col4-6 (col4 without the *_* in the suffix for read names in fa), 5 and 6 are start and end coordinates on the sequence. If the subset txt file exists, which contains each row as a read name, then extract only those in this file. Output to stdout. \n$0 <fa> <SV_call_txt> <subset_read_txt>\n";
}

my ($fa, $SV_call_txt, ) = @ARGV;
my $flank = 1000;
$flank = $ARGV[2] if(defined $ARGV[2]);
my $subset_read_txt= "NA";
$subset_read_txt = $ARGV[3] if(defined $ARGV[3]);
if(defined $ARGV[2] && $ARGV[2] !~ /^\d+$/){
    $subset_read_txt = $ARGV[2];
    $flank = 0;
}
fa::get_substr_by_call($fa, $SV_call_txt, $flank, $subset_read_txt);
