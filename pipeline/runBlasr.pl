#!/usr/bin/perl
use warnings;
use strict;


my ($short_readname, $PB_readname, $blasr, $outPrefix, $suffix, $nproc, $option) = @ARGV;

my @ops = split(/\_/, $option);
my $op = "";
foreach (@ops){
    $op .= " " . join(" ", split(/:/, $_));
}
if($suffix eq "m4"){
    my $com = "$blasr $short_readname $PB_readname $op -m 4 -out $outPrefix.m4";
    print $com;
    system($com);
}
elsif($suffix eq "sam"){
    my $com = "$blasr $short_readname $PB_readname $op -sam -out $outPrefix.sam";
    print $com;
    #$DB::single = 1;
    system($com);
}


