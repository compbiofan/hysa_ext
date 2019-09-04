#!/usr/bin/perl
use warnings;
use strict;
require math;

if(@ARGV == 0){
    die "This takes a call file with the last column the Pacbio gap starting position and size, separated by ; and further ',', summarize each call's gap starting positiona and size's mean standard deviation. The added five columns are 1) # of supporting PB; 2) PB starting position absolute avg; 3) PB starting position std; 4) PB size difference absolute avg; 5) PB size difference std. \nUsage: $0 <file> <col (0-based)>\n";
}

my ($f, $col) = @ARGV;
my $avg_ab_p;
my $std_p;
my $avg_ab_s;
my $std_s;
open fh_, "<$f" or die $!;
while(<fh_>){
    chomp;
    my @aa = split(/\t/, $_);
    if($aa[$col] eq "-1;-1" || $aa[$col] eq "NA;NA"){
        print join("\t", @aa, 0, "NA", "NA", "NA", "NA") . "\n"; 
    }
    else{
        my @bb = split(/\;/, $aa[$col]);
        # distance to the estimated starting position
        my @c1 = split(/\,/, $bb[0]);
        if(scalar(@c1) <= 5){
            print join("\t", @aa, scalar(@c1), "NA", "NA", "NA", "NA") . "\n";
            next;
        }
        else{
            $avg_ab_p = math::avg_ab(\@c1);
            $std_p = math::std(\@c1);
        }
        # distance to the estimated size
        my @c2 = split(/\,/, $bb[1]);
        $avg_ab_s = math::avg_ab(\@c2);
        $std_s = math::std(\@c2);
        print join("\t", @aa, scalar(@c1), $avg_ab_p, $std_p, $avg_ab_s, $std_s) . "\n";
    }
}
close fh_;

