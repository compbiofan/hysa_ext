#!/usr/bin/perl
use warnings;
use strict;

if(scalar(@ARGV) == 0){
    die "This script takes in an fa and reorder according to the cluster index, to which the Illumina reads in fa belong. The output is a large \"fa\" file, but with separator \"#group IL#\" before the reads of each cluster. \nUsage: $0 group_prefix Illumina_fa_prefix output_file whether_include_ref_by_default_off\n";
}

my ($group_f, $fa_f_prefix, $out_f) = @ARGV;
my $include_ref = 0;
$include_ref = 1 if(defined $ARGV[3] && $ARGV[3] eq "ref");

# Store all reads in a hash according to their cluster
my $h;
# reverse retrieval
my $hh;
my $g;
open g_fh, "<$group_f.txt" or die $!;
while(<g_fh>){
    if($_ =~ /^#(\d+)/){
        $g = $1;
    }
    else{
        if($_ =~ /^(\d+)$/){
            $h->{$g}->{$1} = 1;
            $hh->{$1} = $g;
            if($include_ref == 1){
                my $rn = $1 . "_ref";
                $h->{$g}->{$rn} = 1;
                $hh->{$rn} = $g;
            }
        }
    }
}
close g_fh;

# Read fa file and store the sequence in hash
my $fa = "$fa_f_prefix" . "_1.fa";
my $fa2 = "$fa_f_prefix" . "_2.fa";
my $id;
# used to control if the next line is to be analyzed, in the face of pseudo ref reads
my $tag = 0; 
open fa_fh, "<$fa" or die $!;
open fa2_fh, "<$fa2" or die $!;
while(<fa_fh>){
    my $l = <fa2_fh>;
    if($_ =~ /^>(\S+)$/){
        $id = $1;
        $tag = 1;
    }
    elsif($tag == 1){
        chomp;
        if(defined $hh->{$id}){
            my $g = $hh->{$id};
            chomp $l;
            $h->{$g}->{$id} = join("#", $_, $l);
        }
        $tag = 0;
    }
}
close fa_fh;
close fa2_fh;

if($include_ref == 0){
    open out_fh, ">$out_f" or die $!;
    foreach my $g (sort {$a <=> $b} keys %$h){
        my $h_t = $h->{$g};
        my @k = keys %{$h_t};
        print out_fh "#$g\t" . scalar(@k) . "\n";
        foreach my $id (sort {$a <=> $b} @k){
            my $line = $h_t->{$id};
            if($line =~ /#/){
                my @ls = split(/#/, $line);
                print out_fh join("\n", $id . "_1", $ls[0], $id . "_2", $ls[1]) . "\n";
            }
        }
    }
    close out_fh;
}
else{
    open out_fh, ">$out_f" or die $!;
    foreach my $g (sort {$a <=> $b} keys %$h){
        my $h_t = $h->{$g};
        my @k = keys %{$h_t};
        print out_fh "#$g\t" . scalar(@k) . "\n";
        foreach my $id (sort @k){
            my $line = $h_t->{$id};
            if($line =~ /#/){
                my @ls = split(/#/, $line);
                print out_fh join("\n", $id . "_1", $ls[0], $id . "_2", $ls[1]) . "\n";
            }
        }
    }
    close out_fh;
}
