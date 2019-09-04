#!/usr/bin/perl
use warnings;
use strict;
require cigar;

if(@ARGV == 0){
    die "Given a bam file, a bed file, look for the reads that have a gap (del) correspondign to the bed (gap size within a threshold of bed's gap size, gap starting position within a threshold of bed's starting position). For each such read, record its difference from the bed's starting position, and its difference from the bed's size. In the output file, add one column, separated by semicolon are the entries of differences (starting position and size); each further separated by colon. Each read should have good alignment (< clip). \nUsage: $0 <bam> <bed> <out> <pos_t> <size_t> <clip_t> <outer> <inner> <match_perc_t>\n";
}

my ($bam, $bed, $out, $pos_t, $size_t, $clip_t, $outer, $inner, $match_perc_t) = @ARGV;
open fh_, "<$bed" or die $!;
open OUT, ">$out" or die $!;
while(<fh_>){
    chomp;
    my @p_refs;
    my @sizes;
    my $str;
    my @aa = split(/\t/, $_);
    my $s = $aa[1];
    my $size = $aa[2] - $aa[1];
    my $l_s = $s - $outer;
    my $r_s = $s + $inner;
    my $x = `samtools view -c $bam $aa[0]:$l_s-$r_s`;
    if($x > 10000){
        $str = "-1;-1";
        print OUT join("\t", @aa, $str) . "\n";
        next;
    }
    my @r1 = split(/\n/, `samtools view $bam $aa[0]:$l_s-$r_s`);
    foreach my $r (@r1){
        next if($r =~ /^@/);
        my @x = split(/\t/, $r);
        my $cigar = $x[5];
        next if($cigar eq "*");
        my ($matchlen, $readlen) = split(/\t/, cigar::get_matchlen_on_read_of_all($cigar)); 
        if($matchlen/$readlen < $match_perc_t || $readlen - $matchlen > $clip_t){
            next;
        }
        my ($p_ref, $p_read, $tag, $num) = cigar::get_pos_from_cigar_input($r);
        foreach my $i (0 .. scalar(@$p_ref) - 1){
            my $p_ref_ = $p_ref->[$i];
            my $tag_ = $tag->[$i];
            my $num_ = $num->[$i];
            if(abs($p_ref_ - $s) < $pos_t && $tag_ eq "D" && abs($num_ - $size) < $size_t){
                push @p_refs, ($p_ref_ - $s);
                push @sizes, ($num_ - $size);
            }
        }
    }
    if(@p_refs == 0){
        push @p_refs, "NA";
    }
    if(@sizes == 0){
        push @sizes, "NA";
    }
    print OUT join("\t", @aa, join(";", join(",", @p_refs), join(",", @sizes))) . "\n";
}
close fh_;
close OUT;



