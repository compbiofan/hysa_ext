#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;
use File::Basename;

my %opts = (d=>"drp.log", t=>150, m=>200, s=>10, l=>30);
getopts('d:t:m:s:l:', \%opts);
die("
Function: Given a bam file, this script extracts discordant read pairs, clipped reads and reads with gaps, and write the SAM file to the output directory. 
Usage:  $0 <m4> <out_filename> <drp_f> 
Options:
    -d  drp log file, the output [$opts{d}]
    -t  the maximum distance between the read and its reference read's aligned position on the long read, for them to be considered as a pair for comparison. [$opts{t}]
    -m  the minimum difference of insert size between the short read pair's alignment to the reference and the long read, for this long read to be counted as discordant. [$opts{m}]
    -s  If the short read's aligned length on the long read is more than this number of its reference read, the long read is considered as from the alternative allele. [$opts{s}]
    -l  If the reference read aligns to the long read with more than this number of difference between the aligned length on the read itself and that on the long read, the long read is considered as from the alternative allele. [$opts{l}]
") unless (@ARGV);

my ($m4, $out_filename, $drp_f) = @ARGV;
   
my $version = "V0.0_06172016";
print "$0: filter by reference. Version: $version\n";
print "Adapt to sort -h. Previous version assumes 1 and 2 interlace. But with -h they are not. Also a bug losing abs in the previous version. \n";
my $out_drp_filename = $opts{d};
my $align_t = $opts{t};
my $min_len_diff = $opts{m};
my $thresholds = $opts{s};
my $thresholdl = $opts{l};

my $h_drp;
if(-e $drp_f){
    $h_drp = &open_drp($drp_f);
}
my $pre = "NA";
my @bl;
open m4_fh, "<$m4" or die $!;
open OUT, ">$out_filename" or die $!;
open DRP_out, ">$out_drp_filename" or die $!;
while(<m4_fh>){
    chomp;
    my @a = split(/\s+/, $_);
    my @b = split(/\//, $a[0]);
    # set up blocks, each block containing the same read (including _ref)

    if($b[0] !~ /ref/ && $b[0] ne $pre){ 
        # conclude the previous
        if($pre ne "NA"){
            &print_previous();
        }
        # start a new block
        $pre = $b[0];
        @bl = ();
        push @bl, $_;
    }
    elsif($b[0] !~ /ref/ && $b[0] eq $pre){
        push @bl, $_;
    }
    elsif($b[0] =~ /^(\S+)_ref/ && $1 eq $pre){
        push @bl, $_;
    }
    elsif($b[0] =~ /^(\S+)_ref/ && $1 ne $pre){
        # do nothing, ignore this line
    }
}
&print_previous();
close m4_fh;
close OUT;
close DRP_out;
1;


sub print_previous{
    # need to separate pacbio
    my $h;
    my $ret;
    # categorize into read and ref
    foreach my $x (@bl){
        my @xx = split(/\s+/, $x);
        my @xxx = split(/\//, $xx[0]);
        $xxx[0] =~ s/_ref//;

        if($x =~ /ref/){
            if(!defined $h->{$xx[1]}->{$xxx[0]}->{ref}->{$xxx[2]}){
                $h->{$xx[1]}->{$xxx[0]}->{ref}->{$xxx[2]} = $x;
            }
            else{
                $h->{$xx[1]}->{$xxx[0]}->{ref}->{$xxx[2]} .= "\n" . $x;
            }
        }
        else{
            if(!defined $h->{$xx[1]}->{$xxx[0]}->{re}->{$xxx[2]}){
                $h->{$xx[1]}->{$xxx[0]}->{re}->{$xxx[2]} = $x;
            }
            else{
                $h->{$xx[1]}->{$xxx[0]}->{re}->{$xxx[2]} .= "\n" . $x;
            }
        }
    }
    foreach my $pb (keys %$h){
        foreach my $il (keys %{$h->{$pb}}){
            if(defined $h->{$pb}->{$il}->{re} && !defined $h->{$pb}->{$il}->{ref}){
                # check if drp
                if(defined $h_drp->{$il}){
                    my $reads1 = $h->{$pb}->{$il}->{re}->{1};
                    my $reads2 = $h->{$pb}->{$il}->{re}->{2};
                    my @reads1_ = split(/\n/, $reads1);
                    my @reads2_ = split(/\n/, $reads2);
                    my $this_tag = 0;
                    foreach my $r1 (@reads1_){
                        foreach my $r2 (@reads2_){
                            my $insert = &get_insert_size($r1, $r2);
                            if($insert < $h_drp->{$il} - $min_len_diff){
                                print OUT join("\n", $r1, $r2) . "\n";
                                print DRP_out join("\t", $il, $h_drp->{$il}, $insert) . "\n";
                                $this_tag = 1;
                                last;
                            }
                        }
                        last if($this_tag == 1);
                    }
                }
                else{
                    print OUT $h->{$pb}->{$il}->{re}->{1} . "\n";
                    print OUT $h->{$pb}->{$il}->{re}->{2} . "\n";
                }
                next;
            }
            elsif(!defined $h->{$pb}->{$il}->{re} && defined $h->{$pb}->{$il}->{ref}){
                # do nothing
                next;
            }
# for the rest
            $ret = &compare($h->{$pb}->{$il}->{re}, $h->{$pb}->{$il}->{ref});
            print OUT $ret if($ret ne "NA");
        }
    }
    undef %$h;
}

sub compare{
    my ($read, $ref) = @_;
    my $h;
    my @read_1 = split(/\n/, $read->{1});
    my @read_2 = split(/\n/, $read->{2});
    my @ref_1 = split(/\n/, $ref->{1});
    my @ref_2 = split(/\n/, $ref->{2});

    #my ($ret1, $better1) = &compare_read_two_criteria_all_aln(\@read_1, \@ref_1, $threshold);
    #my ($ret2, $better2) = &compare_read_two_criteria_all_aln(\@read_2, \@ref_2, $threshold);
    my ($better1) = &compare_read_two_criteria_all_aln(\@read_1, \@ref_1, $thresholds, $thresholdl);
    my ($better2) = &compare_read_two_criteria_all_aln(\@read_2, \@ref_2, $thresholds, $thresholdl);

    if($better1 > 0 && $better2 >= 0 || $better1 >= 0 && $better2 > 0){
        #return join("\n", $ret1->[0], $ret2->[0]) . "\n";
        return join("\n", $read->{1}, $read->{2}) . "\n";
    }

    return "NA";
}

# The closer the length on the read and pb, the better the alignment -> no indel
#
# The longer the read is mapped to pb, the better the alignment -> short clip 
#
# But if it does not clip, but span, the second criterion does not work.
#
# If always clip, use criterion 1. 
#
# If always span, use criterion 2. 
#
# If clip + span,
#
# diff_re = read alignment length - PB alignment length
# DEL:
#  if pb is ref, 
#  if read span, ref clip (error), use criterion 1 -> small diff_(re_aln - ref_aln), large diff_re  (-), small diff_ref
#  if read clip, ref clip (error), use criterion 2 -> large diff_(re_aln - ref_aln) (-), small diff_re, small diff_ref
#  if pb is var,
#  if read clip (error), ref clip, use criterion 2 -> large diff_(re_aln - ref_aln) (+), small diff_re, small diff_ref
#  if read clip (error), ref insertion span, use criterion 1 -> small diff_(re_aln - ref_aln), small diff_re, large diff_ref (+)
#
#  INS:
#  if pb is ref,
#  if read clip, ref clip (error), use criterion 2 -> large diff_(re_aln - ref_aln) (-), small diff_re, small diff_ref
#  if read insertion span, ref clip (error), use criterion 1 -> small diff_(re_aln - ref_aln), large diff_re (+), small diff_ref
#  if pb is var,
#  if read clip (error), ref clip, use criterion 2 -> large diff_(re_aln - ref_aln) (+), small diff_re, small diff_ref
#  if read clip (error), ref span, use criterion 1 -> small diff_(re_aln - ref_aln), small diff_re, large diff_ref (-)
#
#  in all, 
#  to filter
#  large (-), [small, small]
#  small, large (+, -), small
#
#
#  The first only need large(-). The second criterion needs to look at
#
#  first is small (not satisfying the first criterion)
#  third is small 
#  x = length(re) - length(ref) difference of alignment length on the read, y = length(re_PB) - legnth(ref_PB) difference of alignment length on PB
#  if |x - y| > threshold, PB is from ref allele


sub compare_read_two_criteria_all_aln{
    my ($read, $ref, $thresholds, $thresholdl) = @_;

    my $h_re = &read_read($read);
    my $h_ref = &read_read($ref);

    my $better = 1;
    # as long as there is one var not better, return it (<= 0)
    foreach my $ori (keys %$h_re){
        if(defined $h_ref->{$ori}){
            foreach my $re_range (keys %{$h_re->{$ori}}){
                my ($re_s, $re_e) = split(/\./, $re_range);
                foreach my $ref_range (keys %{$h_ref->{$ori}}){
                    my ($ref_s, $ref_e) = split(/\./, $ref_range);
                    if($ref_s <= $re_s && $re_s <= $ref_e || $re_s <= $ref_s && $ref_s <= $re_e){
                        # overlap
                        $better = &compare_read_two_criteria($h_re->{$ori}->{$re_range}, $h_ref->{$ori}->{$ref_range}, $thresholds, $thresholdl);
                        return $better if($better != 1);
                    }
                }
            }
        }
        else{
            # no ori in ref
        }
    }


    return $better; 


}

sub read_read{
    my ($read) = @_;
    my $h;
    foreach my $r (@$read){
        my @xx = split(/\s+/, $r);
        if($xx[8] == 0){
            $xx[9] -= $xx[5];
            $xx[10] += ($xx[7] - $xx[6]);
        }
        else{
            $xx[9] -= ($xx[7] - $xx[6]);
            $xx[10] += $xx[5];
        }
        #$h->{$xx[8]}->{"$xx[9].$xx[10]"} = $xx[6] - $xx[5];
        $h->{$xx[8]}->{"$xx[9].$xx[10]"} = $r;
    }
    return $h;
}

sub compare_read_two_criteria{
    my ($best_re, $best_ref, $thresholds, $thresholdl) = @_;

    my $better = 1;
    # 1. Look for the best alignment (longest) on read, and compare with ref for each PB.
    # 2. For the aligned length on read (ref - re), set up a threshold. 
    # 3. If > threshold, PB is ref.
    # 4. If <= threshold, check if ref contain indel: abs(length(ref) - length(ref_PB)) > threshold
    # 5. If no, check if read contains indel based on ref's alignment: x=length(ref) - length(re); y=length(ref_PB) - length(re_PB); if |x - y| > threshold, PB is ref.  

    # get the position where the two are closest

    my @xx = split(/\s+/, $best_re);
    my @xx1 = split(/\s+/, $best_ref);

    return $better if($xx[8] != $xx1[8]);
    return $better if(abs($xx[9] - $xx1[9]) >= $align_t);

    my $re_aln = $xx[6] - $xx[5];
    my $ref_aln = $xx1[6] - $xx1[5];

    my $re_PB_aln = abs($xx[9] - $xx[10]);
    my $ref_PB_aln = abs($xx1[9] - $xx1[10]);

    my $re_diff = abs($re_aln - $re_PB_aln);
    my $ref_diff = abs($ref_aln - $ref_PB_aln);

    # 2.
    if($ref_aln - $re_aln > $thresholds){
        $better = -1;
        return $better;
    }
    elsif($ref_aln - $re_aln < -$thresholds){
        return $better;
    }
    else{
        # 4
        if($ref_diff <= $thresholdl){ 
            # 5
            my $x = $ref_aln - $re_aln;
            my $y = $ref_PB_aln - $re_PB_aln;
            if(abs($x - $y) > $thresholds){
                $better = -1;
                return $better;
            }
            else{
                # neutral
                $better = 0;
                return $better;
            }
        }
        else{
            return $better;
        }
    }

    # nothing should be left
    return $better;
}

# Compare to see if any of the read is aligned longer than ref, or not shorter than ref
sub compare_read{
    my ($read, $ref) = @_;
    my $better = -1;
    my @ret1;
    foreach my $r (@$read){
        my @xx = split(/\s+/, $r);
        my $len = $xx[6] - $xx[5];
        my $fail = 0;
        my $tag1 = 0;
        my $ever_enter = 0;
        foreach my $r1 (@$ref){
            my @xx1 = split(/\s+/, $r1);
            next if($xx[8] != $xx1[8]);
            if(abs($xx[9] - $xx1[9]) < $align_t){
                $ever_enter = 1;
                if($len > $xx1[6] - $xx1[5]){
                    $tag1 = 1;
                }
                if($len < $xx1[6] - $xx1[5]){
                    $fail = 1;
                }
            }
        }
        if($ever_enter == 0 || $fail == 0){
            $better = 0;
            push @ret1, $r;
            if($tag1 == 1 || $ever_enter == 0){
                $better = 1;
            }
        }
    }
    return (\@ret1, $better);
}



# return the insert size from a pair of m4 alignments
sub get_insert_size{
    my ($this, $pre) = @_;
    my @a = split(/\s+/, $this);
    my @b = split(/\s+/, $pre);
    my $a_s;
    my $b_s;
    if($a[8] == 0){
        $a_s = $a[9] - $a[5];
    }
    else{
        $a_s = $a[11] - ($a[10] + $a[5]);
        # debug: 06302016
        #$a_s = $a[11] - ($a[9] - $a[5]);
    }
    if($b[8] == 0){
        $b_s = $b[9] - $b[5];
    }
    else{
        #$b_s = $b[11] - ($b[9] - $b[5]);
        $b_s = $b[11] - ($b[10] + $b[5]);
    }
    return abs($a_s - $b_s);
}

sub open_split{
    my ($split_f) = @_;
    my $h_sp;
    open split_fh, "<$split_f" or die $!;
    while(<split_fh>){
        my @a = split(/\t/, $_);
        $h_sp->{$a[0]} = 1;
    }
    close split_fh;

    return $h_sp;
}

sub open_drp{
    my ($drp_f) = @_;
    my $h_drp;
    return $h_drp if(!-e $drp_f);
    open drp_fh, "<$drp_f" or die $!;
    while(<drp_fh>){
        chomp;
        my @a = split(/\t/, $_);
        # added 06302016, drp file starts from 1
        $a[0] -= 1;
        $h_drp->{$a[0]} = $a[1];
    }
    return $h_drp;
}
