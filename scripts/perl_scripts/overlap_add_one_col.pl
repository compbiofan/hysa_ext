#!/usr/bin/perl
use warnings;
use strict;

if(@ARGV == 0){
    die "Given file A and file B in bed format, overlap the two and if A can find overlap in B, add extra column 1, otherwise 0. \nUsage: <fA> <fB> <out> <mode>\n";
}

my ($fA, $fB, $out, $mode) = @ARGV;

if($mode eq "reciprocal50"){
    `intersectBed -a $fA -b $fB -f 0.5 -r -wa > $fA.tmp`;
}
elsif($mode eq "touch"){
    `intersectBed -a $fA -b $fB -wa > $fA.tmp`;
}

my $h;
open fh_, "<$fA.tmp" or die $!;
while(<fh_>){
    chomp;
    $h->{$_} = 1;
}
close fh_;

open OUT, ">$out" or die $!;
open fh_, "<$fA" or die $!;
while(<fh_>){
    chomp;
    if(defined $h->{$_}){
        print OUT join("\t", $_, 1) . "\n";
    }
    else{
        print OUT join("\t", $_, 0) . "\n";
    }
}
close OUT;
close fh_;
`rm $fA.tmp`;
