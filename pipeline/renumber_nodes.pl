#!/usr/bin/perl
use warnings;
use strict;

if(@ARGV == 0){
    die "Given a graph with each line the link between two nodes, reorder them and output the corresponding node names. \nUsage: $0 <m4> <prefix>\n";
}

my ($m4, $prefix) = @ARGV;
my $num_il = 1;
my $num_pb = 1; 
my $h_il;
my $h_pb;
my $h_link;
open fh_, "<$m4" or die $!;
while(<fh_>){
    chomp;
    my @a = split(/\s+/, $_);
    my @b = split(/\//, $a[0]);

    if(!defined $h_il->{$b[0]}){
        $h_il->{$b[0]} = $num_il;
        $num_il ++;
    }
    if(!defined $h_pb->{$a[1]}){
        $h_pb->{$a[1]} = $num_pb;
        $num_pb ++;
    }
    my $il = $h_il->{$b[0]};
    my $pb = $h_pb->{$a[1]};
    $h_link->{$il}->{$pb} = 1;
}
close fh_;
open OUT_link, ">$prefix.link" or die $!;
foreach my $il (sort {$a <=> $b} keys %$h_link){
    foreach my $pb (sort {$a <=> $b} keys %{$h_link->{$il}}){
        print OUT_link join("\t", $il, ($pb + $num_il - 1)) . "\n";
    }
}

close OUT_link;
open OUT_node, ">$prefix.node" or die $!;
foreach my $il (keys %$h_il){
    print OUT_node join("\t", $h_il->{$il}, $il) . "\n";
}
foreach my $pb (keys %$h_pb){
    print OUT_node join("\t", $h_pb->{$pb} + $num_il - 1, $pb) . "\n";
}
close OUT_link;


