#!/usr/bin/perl
use warnings;
use strict;

if(scalar(@ARGV) == 0){
    die "This reads input fa file containing pairs of reads, which contains header lines beginning with # followed by group number, print to output_prefix.1.fa and output_prefix.2.fa the corresponding reads in that group. \nUsage: $0 group_num sorted_fa output_prefix\n";
}

my ($g, $fa, $o_pre) = @ARGV;

my $IL_o_f_1 = $o_pre . ".1.fa";
my $IL_o_f_2 = $o_pre . ".2.fa";
&extract_n_print($fa, $g, $IL_o_f_1, $IL_o_f_2);
1;

sub extract_n_print{
    my ($IL_f, $g, $IL_o_f_1, $IL_o_f_2) = @_;
    open IL_fh, "<$IL_f" or die $!;
    open IL_o_fh_1, ">$IL_o_f_1" or die $!;
    open IL_o_fh_2, ">$IL_o_f_2" or die $!;
    my $tag = 0;
    my $tag1 = -1;
    while(<IL_fh>){
        if($_ =~ /^#$g\s+/){
            $tag = 1;
        }
        elsif($tag == 1 && $_ !~ /^#/){
            if($_ =~ /_1$/ || $tag1 == 1){
                if($_ =~ /_1$/){
                    $tag1 = 1;
                    $_ =~ s/_1$//;
                    $_ = ">$_" if($_ !~ /^>/);
                }
                else{
                    $tag1 = 0;
                }
                print IL_o_fh_1 $_;
            }
            elsif($_ =~ /_2$/ || $tag1 == 2){
                if($_ =~ /_2$/){
                    $tag1 = 2;
                    $_ =~ s/_2$//;
                    $_ = ">$_" if($_ !~ /^>/);
                }
                else{
                    $tag1 = 0;
                }
                print IL_o_fh_2 $_;
            }
        }
        elsif($tag == 1 && $_ =~ /^#/){
            last;
        }
    }
    close IL_fh;
    close IL_o_fh_1;
    close IL_o_fh_2;
}
