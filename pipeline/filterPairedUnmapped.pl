#!/usr/bin/perl
use warnings;
use strict;
require base;
use File::Basename;

# 08252015 Without changing the details, now add the indexing with the starting index to combine one_mapped and both_unmapped together. 
# 07242015 After done with both unmapped test, in the rest of the test, change || to && for the pair. This is based on that if one end has problem due to low qual, the pair probably does not contain breakpoint and is not meaningful to keep it.

if($#ARGV == -1){
    die "This reads a sam file, with pairs in consecutive rows, check the base quality of both, if both pass the filter, output them with the readname changed to the index starting from the startingIndex, and write the index file with the suffix .id in the output file. To avoid the missing of reads containing breakpoints due to low quality small ends, make len_of#_to_filter large, > 50% read length. If not wanting to rewrite the read name by new index, use NA as startingIndex. \nUsage: $0 <sam> <output_prefix> <startingIndex> <base_qual_threshold> <len_of_#_to_filter> <min_clipped_len>\n";
}
my $version = "0.0_10212015";
my $sam = $ARGV[0];

my $sam_dir = dirname($sam);
my $output_prefix = "unmapped.both.filtered";
$output_prefix = $ARGV[1] if(defined $ARGV[1]);
my $startingIndex = 0;
$startingIndex = $ARGV[2] if(defined $ARGV[2]);
my $baseQual = 35;
$baseQual = $ARGV[3] if(defined $ARGV[3]);
my $len = 20;
$len = $ARGV[4] if(defined $ARGV[4]);
my $min_clipped = 20;
$min_clipped = $ARGV[5] if(defined $ARGV[5]);

my $base_pm_version = base::get_version();
print "#Command: $0 $sam $output_prefix $startingIndex $baseQual $len $min_clipped\n#Version: $version\n#base.pm version: $base_pm_version\n";
my $tag = 1;
my $line;
my $id = "";
if($startingIndex ne "NA"){
    $id = $startingIndex;
}
open SAM_out, ">$sam_dir/$output_prefix.sam" or die $!;
open ID_out, ">$sam_dir/index.$output_prefix.txt" or die $!;
open SAM, "<$sam" or die $!;
while(<SAM>){
    next if($_ =~ /^@/);
    if($tag == 1){
        $line = $_;
        $tag = 0;
    }
    else{
        my @a1 = split(/\s+/, $line);
        my @a2 = split(/\s+/, $_);
        if(base::check_basequal($a1[10], $a1[5], $baseQual, $len, $min_clipped) && base::check_basequal($a2[10], $a2[5], $baseQual, $len, $min_clipped)){
            if($startingIndex eq "NA"){
                print SAM_out $line, $_;
            }
            else{
                print SAM_out join("\t", $id, @a1[1 .. $#a1]) . "\n", join("\t", $id, @a2[1 .. $#a2])."\n";
                print ID_out join("\t", $a1[0], $id) . "\n";
                $id ++;
            }
        }
        $tag = 1;
    }
}
close SAM;
close SAM_out;
close ID_out;
