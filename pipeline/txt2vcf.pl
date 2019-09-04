#!/usr/bin/perl
use warnings;
use strict;


if($#ARGV == -1){
    die "This converts a txt file into a vcf file. \nUsage: txt2vcf.pl <txt_file> <vcf_file>\n";
}

my ($txt, $vcf) = @ARGV;

my $str;
my $chr;
my $pos;
my $sz;
my $end;
my $sp;
my $filter;
my $h;
open out_fh, ">$vcf" or die $!;
open fh_, "<$txt" or die $!;
while(<fh_>){
    chomp;
    my @x = split(/\t/, $_);
    $chr = $x[0];
    $pos = $x[1];
    $sz = $x[2];
    if($x[3] eq "I"){
        $str = "INS";
        $end = $x[1] + 1;
    }
    elsif($x[3] eq "D"){
        $str = "DEL";
        $end = $x[1] + $sz;
        $sz *= -1;
    }
    $sp = $x[$#x];
    $filter="PASS";
    if($sp <= 2){
        $filter = "LOWQUAL";
    }
    $h->{$chr}->{$pos} = join("\t", $chr, $pos, ".", ".", $str, ".", $filter, "CIPOS=-50,50;CIEND=-50,50;SVTYPE=$str;END=$end;SVLEN=$sz", "SPNUM", $sp) . "\n"; 
}
close fh_;
        print out_fh qq(##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion">
##FORMAT=<ID=SPNUM,Number=1,Type=Integer,Description="Number of supporting Illumina reads.">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT  SAMPLE
);

foreach my $c (sort keys %$h){
    foreach my $p (sort {$a <=> $b} keys %{$h->{$c}}){

        print out_fh $h->{$c}->{$p};
    }
}
close out_fh;


