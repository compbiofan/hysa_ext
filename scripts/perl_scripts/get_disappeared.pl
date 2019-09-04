#!/usr/bin/perl
use warnings;
use strict;

if(scalar(@ARGV) == 0){
    die "Given an fa file and a sam file, output the contigs in fa file that does not appear in sam to stdout. \nUsage: $0 <sam> <fa>\n";
}

my ($sam, $fa) = @ARGV;
my $h_sam;
# read sam
open sam_fh, "<$sam" or die $!;
while(<sam_fh>){
    next if($_ =~ /^@/);
    my @a = split(/\t/, $_);
    my @b = split(/\//, $a[0]);
    $h_sam->{join("/", @b[0 .. $#b-1])} = 1;
}
close sam_fh;

# read fa
open fa_fh, "<$fa" or die $!;
while(<fa_fh>){
    if($_ =~ /^>(\S+)/){
        if(!defined $h_sam->{$1}){
            print $1 . "\n";
        }
    }
}
close fa_fh;


