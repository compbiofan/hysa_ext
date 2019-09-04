#!/usr/bin/perl
use warnings;
use strict;

if(@ARGV == 0){
    die "This is to to remove the duplicated calls from one bed file. By default mode is reciprocal50. It could also be touch.\nUsage: $0 <bed> <mode>\n";
}

my $mode = "reciprocal50";
my $bed_file = $ARGV[0];
$mode = $ARGV[1] if(defined $ARGV[1]);
my $bin = "intersectBed";
my $c = 0.5;
if($mode eq "reciprocal50"){
    `$bin -a $bed_file -b $bed_file -wo -f $c -r > tmp.intersect.bed`;
}
elsif($mode eq "touch"){
    `$bin -a $bed_file -b $bed_file -wo > tmp.intersect.bed`;
}

# pointing to the parent
my $h_p;
# record the samples
open tmp_fh, "<tmp.intersect.bed" or die $!;
while(<tmp_fh>){
    chomp;
    my @a = split(/\t/, $_);
    my $num = scalar(@a)/2;
    my $str = "$a[0]:$a[1]-$a[2]";
    my $passenger = "$a[$num]:$a[$num + 1]-$a[$num + 2]";
    my $root;
    if(defined $h_p->{$passenger}){
        my $root_pas = &get_root($h_p, $passenger);
        my $root_this = &get_root($h_p, $str);
        if($root_pas ne $root_this){
            $h_p->{$root_pas} = $root_this;
        }
        $root = $root_this;
    }
    elsif(!defined $h_p->{$passenger}){
        # passenger's root never defined before, add its root the str
        $root = &get_root($h_p, $str);
        $h_p->{$passenger} = $root if($root ne $passenger);
    }
}
close tmp_fh;

open tmp_fh, "<$bed_file" or die $!;
my $root_fh;
while(<tmp_fh>){
    chomp;
    my @a = split(/\t/, $_);
    my $str = "$a[0]:$a[1]-$a[2]";
    my $root = &get_root($h_p, $str);
    if($root eq $str && !defined $root_fh->{$root}){
        print $_ . "\n";
        $root_fh->{$root} = 1;;
    }
}
close tmp_fh;
1;
sub get_root{
    my ($h_p, $str) = @_;
    while(defined $h_p->{$str}){
        $str = $h_p->{$str};
    }
    return $str;
}





