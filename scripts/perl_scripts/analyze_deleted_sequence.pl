#!/usr/bin/perl
use warnings;
use strict;
require sam;

if(scalar(@ARGV) == 0){
    die "match/ref > per && match/ctg > per, largest gap < indel_size, largest clip < clip_size. This looks at a sam file, with each read aligned to at most 3 places by blasr, classify the reads into well-aligned and poor-aligned from the best alignment. Well-aligned: > per length are aligned (not clipped), the largest clipped length is < min_clip, and there is no large ins or del (size < indel_size). Output is a two-line stdout, with #well_aligned: in the first, and #poor_aligned: in the second. \nUsage: $0 <sam> <per> <min_clip> <indel_size>\n";
}

my ($sam) = $ARGV[0];
my $per = 0.9;
$per = $ARGV[1] if(defined $ARGV[1]);
my $min_clip = 500;
$min_clip = $ARGV[2] if(defined $ARGV[2]);
my $indel_size = 300;
$indel_size = $ARGV[3] if(defined $ARGV[3]);

sam::classify_aln_del($sam, $per, $min_clip, $indel_size);
