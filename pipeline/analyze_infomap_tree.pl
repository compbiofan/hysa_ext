#!/usr/bin/perl
use strict;
use warnings;

if(@ARGV == 0){
    die "Given an output file from the infomap, make a new union find all and txt file. Union find all file is to indicate the total cluster number from union find. If it is given, the indexing should begin from this number. \nUsage: $0 <tree_file> <node_corresponding_file> <prefix> <base>\n";
}

my ($tree_file, $node_file, $prefix, $union_find_all) = @ARGV;
my $base = 0;
if(defined $union_find_all){
    my @bases = split(/\s+/, `wc -l $union_find_all`);
    $base = $bases[0];
}

my $node_h = &read_node($node_file);
my ($tree_h, $tree_h_keys) = &read_tree($tree_file, $node_h);

if(-e "$prefix.all"){
    die "$prefix.all already exists.\n";
}
# change from > to >>, to merge infomap and union find together
open fh_all, ">>$prefix.all";
open fh_txt, ">>$prefix.txt";
foreach my $tree_node_id (sort {$a <=> $b} keys %$tree_h){
    my $str = join("\t", ($tree_node_id + $base), scalar(keys %{$tree_h->{$tree_node_id}->{il}}), scalar(keys %{$tree_h->{$tree_node_id}->{pb}}));
    print fh_all "$str\n";
    print fh_txt "#$str\n";
    foreach my $il (keys %{$tree_h->{$tree_node_id}->{il}}){
        print fh_txt $il . "\n";
    }
    foreach my $pb (keys %{$tree_h->{$tree_node_id}->{pb}}){
        print fh_txt $pb . "\n";
    }
}
close fh_all;
close fh_txt;
       
sub read_tree{
    my ($tree, $node_h) = @_;
    my $h;
    my $h_keys;
    my $num = 0;
    open tr_fh, "<$tree" or die $!;
    while(<tr_fh>){
        chomp;
        next if($_ =~ /^#/);
        my @a = split(/\s+/, $_);
        my @b = split(/:/, $a[0]);
        my $str = join(":", @b[0 .. $#b - 1]);
        if(!defined $h_keys->{$str}){
            $h_keys->{$str} = $num;
            $num ++;
        }
        my $cluster = $h_keys->{$str};
        my $read = $node_h->{$a[$#a]};
        if($read =~ /\./){
            $h->{$cluster}->{pb}->{$read} = 1;
        }
        else{
            $h->{$cluster}->{il}->{$read} = 1;
        }
    }
    close tr_fh;
    return ($h, $h_keys);
}

sub read_node{
    my ($node_file) = @_;
    my $h;
    open fh_, "<$node_file" or die $!;
    while(<fh_>){
        chomp;
        my @a = split(/\s+/, $_);
        $h->{$a[0]} = $a[1];
    }
    close fh_;
    return $h;
}
