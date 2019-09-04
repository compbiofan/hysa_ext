#!/usr/bin/perl
use warnings;
use strict;

if(@ARGV == 0){
    die "Given a txt file with each line one column with the group number, find the actual txt from the second txt file. Used in insertion analysis.\nUSage: $0 <txt1> <txt2>\n";
}

my ($txt1, $txt2)  = @ARGV;
open fh_, "<$txt2" or die $!;
my $h;
while(<fh_>){
    chomp;
    my @a = split(/\t/, $_);
    my @b = split(/\//, $a[4]);

    $h->{$b[0]} = $_;
}
close fh_;
open fh_, "<$txt1" or die $!;
while(<fh_>){
    chomp;
    my @a = split(/:/, $_);
    if(defined $h->{$a[0]}){
        print $h->{$a[0]} . "\n";
    }
}
close fh_;

