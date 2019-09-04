#!/usr/bin/perl
use warnings;
use strict;

if(@ARGV == 0){
    die "Given a txt file with the group number, and the txt file with group number and the coordinates, print the coordinates with the group number appearing in the first txt file. Note: second txt is in bed + txt (for ins, -50, +50 for the two breakpoints). \nUsage: $0 <first_txt> <second_txt>\n";
}

my ($txt1, $txt2) = @ARGV;
my $h;
open fh_, "<$txt2" or die $!;
while(<fh_>){
    my @x = split(/\t/, $_);
    my @y = split(/\//, $x[5]);
    $h->{$y[0]} = join("\t", @x[0 .. 3]);
}
close fh_;
open fh_, "<$txt1" or die $!;
while(<fh_>){
    my @x = split(/:/, $_);
    print $h->{$x[0]} . "\n";
}
close fh_;
