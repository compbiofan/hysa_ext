#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;
use File::Basename;
require sam;
require seq;

#if($#ARGV < 1){
#    die "Function: This file runs blasr alignment from the contigs to the reference, infer both SV and INDEL from the alignment and output in a breakdancer format." 
#}

# The option g should be the uniq ctg number counting from 1. If we are processing the 2nd ctg, then the range is 2:2. 
# we use bestn 3 for both blasr alignment (from PB to ref, from IL to PB)
my %opts = (d=>"out_ctg_aln", P=>"out", I=>"NA", R=>"NA", a=>"NA", g=>"NA", m=>500, q=>500, B=>"NA", b=>"blasr", O=>"NA", C=>4, K=>"NA", S=>1200, c=>20, r=>20, i=>10, n=>3, p=>70, M=>100, s=>10, X=>"NA", N=>1, Y=>"NA", x=>100, E=>"./out.all.txt", F=>"./out.all.acc.txt", T=>"NA", e=>10, y=>0.1, Z=>50, J=>"NA", D=>"pairHMM", z=>20, Q=>50);
getopts('d:P:I:AR:a:g:B:b:O:C:S:K:c:r:m:q:i:n:p:M:s:X:N:Y:x:E:F:T:e:y:Z:J:D:z:Q:', \%opts);
die("
    Usage:  $0 <ctg_file> > <output>
    Options:
    -d  temporary directory [$opts{d}]
    -P  file prefix [$opts{P}]
    -I  Illumina read prefix (followed by _1.fa and _2.fa) [$opts{I}]
    -A  skip the rest of the code after Blasr alignment from ctg to reference is done [false]
    -R  reference fasta file [$opts{R}]
    -a  Union find result with Illumina read index [$opts{a}]
    -g  Processing range for parallelization [$opts{g}]
    -B  BLASR option: blasr binary with full path [$opts{B}]
    -b  BLASR option: new blasr binary for aligning short reads to ctg [$opts{b}]
    -O  BLASR option: blasr options if not using default ones. NA if using default. [$opts{O}]
    -C  BLASR option: processor number [$opts{C}]
    -S  BLASR option: insert size in filter for paired-end reads [$opts{S}]
    -K  BLASR option: skip running blasr from PB to reference if this option is given with sam result [$opts{K}]
    -c  SV option: max gap on the contig [$opts{c}]
    -r  SV option: max gap on the reference [$opts{r}]
    -m  SV option: minimum clip to confirm the contig has two alignments coming up with an SV [$opts{m}]
    -q  INDEL option: maximum allowed clip to confirm the contig has a head to nose alignment for INDEL [$opts{q}]
    -i  INDEL option: minimum indel size [$opts{i}]
    -n  INDEL option: minimum num of split reads supporting the breakpoint [$opts{n}]
    -p  INDEL option: short read percentage of identity to PB surrounding breakpoints [$opts{p}]
    -M  INDEL option: size of match surrounding breakpoints until hit a large INDEL [$opts{M}]
    -s  INDEL option: flanking region size surrounding breakpoints to be checked for short read alignment [$opts{s}]
    -X  Illumina option: Ordered Illumina read file. This will save some retrieval time of Illumina reads [$opts{X}]
    -x  Illumina option: search breakpoint indicated by alignment of Illumina to reference, within [$opts{x}]bp of those called by Pacbio. 
    -N  Illumina option: used in the same step as -x, but the minimum number of support Illumina reads on each end. [$opts{N}]
    -Y  Illumina option: used in the same step as -x. Ordered Illumina sam file. [$opts{Y}]
    -E  If no Illumina boundary, this is the output breakpoint bed file. [$opts{E}]
    -F  If Illumina boundary, this is the final output of with more accurate breakpoint. [$opts{F}]
    -T  Trio sample names. For example, YRI trio samples are 38:39:40. [$opts{T}]
    -e  INDEL option: the maximum size of small matching strings in between gaps to qualify the gaps for pairHMM. Larger leads to more pairHMM. [$opts{e}]
    -y  pairHMM option: modification error penalty. The larger, the lower the penalty. [$opts{y}]
    -J  pairHMM option: pairHMM trans matrix file. [$opts{J}] 
    -Z  post-processing option: the length of microinsertion allowed in deletion. [$opts{Z}] 
    -D  the mode of local alignment, including pairHMM and bwa. If complex DEL, choose bwa. [$opts{D}]
    -z  BWA option: Minimum length of M in bwa mem when D is bwa. [$opts{z}] 
    -Q  BWA option: Padding of local sequences for bwa mem when D is bwa. Should be small to avoid multiple small INDELs inside the anchors. [$opts{z}] 
    \n
    ") unless (@ARGV);

my $path = dirname($0);
my ($ctg) = @ARGV;
my $output_dir = $opts{d};
my $prefix = $opts{P};
my $short_prefix = $opts{I};
my $skip_infer_SV = 0;
$skip_infer_SV = 1 if(defined $opts{A});
my $ref = $opts{R};
my $union_find_txt = $opts{a};
my $processing_range = $opts{g};
my $min_clip_SV = $opts{m};
my $min_clip_INDEL = $opts{q};
my $blasr_binary = $opts{B};
my $blasr_binary_new = $opts{b};
my $blasr_option = $opts{O};
my $nproc = $opts{C};
my $given_sam = $opts{K};
my $insert_size = $opts{S};
my $max_gap_on_ctg = $opts{c};
my $max_gap_on_ref = $opts{r};
my $min_INDEL_sz = $opts{i};
my $min_spt_reads = $opts{n};
my $perc = $opts{p};
my $sz_M = $opts{M};
my $s = $opts{s};
my $random_access = 1;
my $short_prefix_sort_by_group = $opts{X};
my $sz_m = $opts{e};
if($short_prefix_sort_by_group ne "NA"){
    $random_access = 0;
}
# if no IL boundary, then this b_bed file is the output
my $b_bed = $opts{E};
my $b_bed_2 = $opts{F};
# the following options are used when random access is false
my $flanking_size = $opts{x};
my $num_sp_reads = $opts{N};
my $IL_sam = $opts{Y};
my $IL_boundary = 1;
my $trio = $opts{T};
my $modification_error = $opts{y};
my $infer_INDEL_PB_wPairHMM_memEffi = $opts{J};
my $microins = $opts{Z};
my $aln_mode = $opts{D};
my $bwa_M_threshold = $opts{z};
my $bwa_padding = $opts{Q};

if($IL_sam eq "NA"){
    print STDERR "#Warnings: Ordered Illumina sam file is not given. Cannot use Illumina read alignment to reference get more accurate breakpoint.\n";
    $IL_boundary = 0;
}

if($blasr_option eq "NA"){
    $blasr_option = "-bestn 3 -maxAnchorsPerPosition 100 -advanceExactMatches 10 -affineAlign -affineOpen 100 -affineExtend 0 -insertion 5 -deletion 5 -extend -maxExtendDropoff 20 -clipping subread -nproc $nproc";
}
# beginning
`mkdir $output_dir` if(!-d $output_dir);
my $output_file = "$output_dir/$prefix.sv";
open output_fh, ">$output_file" or die $!;
# to store putative calls before IL check

# First align the contig fasta to the reference, submit a job using blasr
my $sam = "NA";
if($given_sam eq "NA"){
    &runBLASR_ctg_to_ref($ctg, $ref, $blasr_binary, $blasr_option, $output_dir, $prefix);
    $sam = "$output_dir/$prefix.sam";
}
else{
    $sam = $given_sam;
}

exit 0 if($skip_infer_SV == 1); 
# Second for all alignments to one contig, infer SV/INDEL, extract Illumina reads corresponding to the cluster of this contig, align them to the ctg and confirm breakpoints. Report the SV/INDEL if confirmed.

# extract the alignments of one contig and write to a sam file
# record the size, the type of SV/INDELs in hash
my $h;
#my $sam = "$output_dir/$prefix.sam";
my $prev = "NA";
my $lSAM = "$output_dir/$prefix.tmp.sam";
open SAM, "<$sam" or die $!;
open lSAM, ">$lSAM" or die $!;
close lSAM;
my $start_l = 0;
my $end_l = "NA";
($start_l, $end_l) = split(/:/, $processing_range) if($processing_range ne "NA");
my $count = 0;
my $tag = 0;
#open b_bed, ">$b_bed" or die $!;
while(<SAM>){
    next if($_ =~ /^@/);
    my $line = $_;
    chomp $line;
    my @a = split(/\t/, $line);
    if($a[0] ne $prev && $prev ne "NA"){
        # process previous ones
        #close lSAM;
        $count ++;
        $prev = $a[0];
        last if($end_l ne "NA" && $count > $end_l);
        next if($count < $start_l);
        $tag = 1;
        if($count > $start_l){
            # not the first time coming here
            $h = &call_SV_INDEL($h, $lSAM, $union_find_txt, $short_prefix, $output_dir, $prefix, $output_file);
            `rm $lSAM`;
        }
        # begin a new session
        #open lSAM, ">$lSAM" or die $!;
        #print lSAM $line;
        `echo "$line" >> $lSAM`;
    }
    elsif($prev eq "NA"){
        $count ++;
        $prev = $a[0];
        last if($end_l ne "NA" && $count > $end_l);
        next if($count < $start_l);
        $tag = 1;
        #open lSAM, ">$lSAM" or die $!;
        #print lSAM $line;
        `echo "$line" >> $lSAM`;
    }
    elsif($a[0] eq $prev){
        #print lSAM $line;
        `echo "$line" >> $lSAM` if($tag == 1);
    }
}
close SAM;
$h = &call_SV_INDEL($h, $lSAM, $union_find_txt, $short_prefix, $output_dir, $prefix, $output_file);
#close b_bed;
#open b_bed_2, ">$b_bed_2" or die $!;
#close b_bed_2;
close output_fh;
if($IL_boundary == 1){
    `perl $path/tune_bp_by_IL.pl $b_bed $IL_sam $flanking_size $num_sp_reads > $b_bed_2`;
    #print "perl $path/tune_bp_by_IL.pl $b_bed $IL_sam $flanking_size $num_sp_reads\n";
}

#close lSAM;
`rm $lSAM`;
&print_SV_INDEL_range($h);

print STDERR "Successfully infer SVs. Done. \n"; 

1;
# now print SV/INDEL ranges
sub print_SV_INDEL_range{
    my ($h) = @_;
    foreach my $str (sort keys %$h){
        foreach my $type (sort keys %{$h->{$str}}){
            foreach my $size (sort {$a <=> $b} keys %{$h->{$str}->{$type}}){
                print join("\t", "#" . $str, $type, $size, $h->{$str}->{$type}->{$size}) . "\n";
            }
        }
    }
}

=cut
sub runBLASR_ctg_to_ref{
    my ($ctg, $ref, $blasr_binary, $blasr_option, $output_dir, $output_prefix) = @_;
    print "perl ~/combinePBIL/scripts/runBlasr.pl $ctg $ref $blasr_binary $output_dir/$output_prefix sam $nproc \"$blasr_option\"";
    `perl ~/combinePBIL/scripts/runBlasr.pl $ctg $ref $blasr_binary $output_dir/$output_prefix sam $nproc "$blasr_option"`;
}
=cut

# call on each contig
sub call_SV_INDEL{
    my ($h, $lSAM, $union_find_txt, $short_prefix, $output_dir, $output_prefix) = @_;
    # First extract short reads and run alignment to ctg 
    my $ctg_fa = "$output_dir/$output_prefix.fa";
    #print "$ctg_fa\n";
    my $g = "NA";
    open lSAM, "<$lSAM" or die $!;
    while(<lSAM>){
        next if($_ =~ /^@/);
        #print $_;
        my @a = split(/\t/, $_);
        # write the contig to a fasta file
        if($a[5] !~ /^\d+H/ && $a[5] !~ /\d+H$/){
            # write the contig in the original sequence
            open FA, ">$ctg_fa" or die $!;
            if($a[1] & 0x0010){
                print FA join("\n", ">$a[0]", seq::revcom($a[9])) . "\n";
            }
            else{
                print FA join("\n", ">$a[0]", $a[9]) . "\n";
            }
            close FA;
        }
        if(-e $ctg_fa && $_ =~ /^g(\d+)/){
            $g = $1;
            `mkdir $output_dir/g$g`;
            my $group_file = "$output_dir/g$g/out.g$g.txt";
            my $small_short_prefix = "$group_file.IL";
            #last if(-e "$group_file.IL.1.fa");
            if($random_access == 1){
                # generate the group file
                `perl $path/generate_group_file.pl $union_find_txt $g $group_file`;
                #print "perl ~/u/generate_group_file.pl $union_find_txt $g $group_file\n";
                # generate the reads
                `perl $path/make_group_fa_dir.pl $group_file $short_prefix NA NA IL_only`;
                #print "perl ~/u/make_group_fa_dir.pl $group_file $short_prefix NA NA IL_only\n";
            }
            else{
                `perl $path/get_IL_by_group.pl $1 $short_prefix_sort_by_group $small_short_prefix`;
                #print "perl ~/combinePBIL/scripts/get_IL_by_group.pl $1 $short_prefix_sort_by_group $small_short_prefix\n";
            }
            # make a sam aligning extracted short read to the ctg
            my ($ctg_prefix) = ($ctg_fa =~ /^(.+)\.fa/);
            # results are in out.sam
            #`perl ~/combinePBIL/scripts/runBLASR.bigPipe.lite.pl -b $blasr_binary_new -c 1 -p $nproc -f sam -n 1 -i $insert_size -o $output_dir/g$g $small_short_prefix $ctg_prefix`;
            `perl $path/runBLASR.bigPipe.lite.pl -b $blasr_binary_new -c 1 -p $nproc -f sam -n 1 -i $insert_size -o $output_dir/g$g $small_short_prefix $ctg_prefix`;
            #print "perl ~/combinePBIL/scripts/runBLASR.bigPipe.lite.pl -b $blasr_binary -c 1 -p $nproc -f sam -n 1 -i $insert_size -o $output_dir/g$g $small_short_prefix $ctg_prefix\n";
            #print "perl ~/combinePBIL/scripts/runBLASR.bigPipe.lite.pl -b $blasr_binary_new -c 2 -p $nproc -f sam -n 2 -i $insert_size -o $output_dir/g$g $small_short_prefix $ctg_prefix";

            # clean up
            #`rm $ctg_fa`;
            #`rm $group_file*`;
            #`rm $output_dir/g$g/blasr*`;
            last;
        }
    }
    close lSAM;
    my $cSAM = "$output_dir/g$g/out.sam";
    # Second analyze the alignment of PB to ref (lSAM) and IL to PB (cSAM) and infer/confirm SV/INDEL
    print STDERR "g$g\n";
    if($g eq "NA"){
        `rm -r $output_dir/g$g`;
        `rm $ctg_fa`;
        return $h;
    }
    my $a = sam::infer_SV_PB($lSAM, $min_clip_SV, $max_gap_on_ctg, $max_gap_on_ref);
    $a = sam::check_INDEL_bp($cSAM, $a, $min_spt_reads, $perc, $s);
    print join("\n", @$a) . "\tSV\n" if(scalar(@$a) != 0);
    my $b;
    if($aln_mode eq "pairHMM"){
        $b = sam::infer_INDEL_PB_wPairHMM_memEffi($lSAM, $min_clip_INDEL, $min_INDEL_sz, $sz_M, $sz_m, $ctg, $ref, "$output_dir/g$g", $modification_error, $microins, $infer_INDEL_PB_wPairHMM_memEffi, $infer_INDEL_PB_wPairHMM_memEffi);
    }
    elsif($aln_mode eq "bwa"){
        $DB::single = 1;
        $b = sam::infer_INDEL_PB_wLocalBWA($lSAM, $min_clip_INDEL, $min_INDEL_sz, $sz_M, $sz_m, $ctg, $ref, "$output_dir/g$g", $bwa_M_threshold, $bwa_padding);
    }
    #print join("\n", @$b) . "\n";
    print output_fh join("\n", @$b) . "\n" if(scalar(@$b) != 0);
    #my $length_diff_t = 20;
    #my $max_match = 500;
    #my $b = sam::infer_INDEL_PB_wRealign($lSAM, $min_clip_INDEL, $min_INDEL_sz, $sz_M, $ref, $length_diff_t, $max_match);
    #print b_bed "Putative: " . join("\n", @$b) . "\n" if(scalar(@$b) != 0);
    $b = sam::check_INDEL_bp($cSAM, $b, $min_spt_reads, $perc, $s, $trio);
    print join("\n", @$b) . "\n" if(scalar(@$b) != 0);
    #`rm -r $output_dir/g*`;
    # TODO this print will be moved after tune by IL
    #print b_bed join("\n", @$b) . "\n" if(scalar(@$b) != 0);
    # get stats from a and b
    #$h = &get_stats($a, $b, $h);
    #print join("\n", @$a) . "\n" if(scalar(@$a) != 0);
    `rm -r $output_dir/g$g`;
    `rm $ctg_fa`;
    return $h;
}

# make a stat count of the calls of SV ($a), INDEL ($b) and store them in hash $h
sub get_stats{
    my ($a, $b, $h) = @_;
    foreach (@$b){
        my @x = split(/\t/, $_);
        if(!defined $h->{INDEL}->{$x[3]}->{$x[2]}){
            $h->{INDEL}->{$x[3]}->{$x[2]}  = 1;
        }
        else{
            $h->{INDEL}->{$x[3]}->{$x[2]} ++;
        }
    }
    foreach (@$a){
        my @x = split(/\t/, $_);
        if(!defined $h->{SV}->{$x[6]}){
            $h->{SV}->{$x[6]}->{$x[8]} = 1;
        }
        else{
            $h->{SV}->{$x[6]}->{$x[8]} ++;
        }
    }
    return $h;
}


