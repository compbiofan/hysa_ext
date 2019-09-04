#!/risapps/rhel6/perl/5.10.1/bin/perl
use warnings;
use strict;
use Getopt::Std;
use File::Basename;
require base;

if(scalar(@ARGV) < 2){
    die "This takes chromosome as the fourth input argument, so that memories of jobs are limited. This extracts discordant read pairs (including unmapped ones), index read names, filter by insert size and base qualities, convert to a pair of fasta files. For those extracted due to split reads, write extra lines with read name changed to id_ref, but now the split part is on the reference side. For those from abnormal insert size, write an extra file with each the id and insert size separated by tab. \nUsage: $0 <bam> <output_dir> <reference_of_bam> <round_num, 1: extracting putative read names per chromosome; 2: whole genome examining putative read pairs; 3: whole genome examining inter-chromosomal alignments> <chr> <region_for_stat:same for all chromosomes> <if_get_ref> <base_qual_threshold> <low_quality_end_len> <insert_size_std> <min_SV> <min_clipped>\n";
}

my $debug = 0;
my $bam = $ARGV[0];
my $out_dir = $ARGV[1];
my $ref = $ARGV[2];
die "Reference is required!\n" if(!defined $ref);
my $round = $ARGV[3];
my $chr = $ARGV[4];
my $stat_region = $ARGV[5];

my $op_base = 6;
my $get_ref = 1;
$get_ref = $ARGV[$op_base] if(defined $ARGV[$op_base]);
$op_base ++;
my $ref_block_size = 100000;
$ref_block_size = $ARGV[$op_base] if(defined $ARGV[$op_base]);
$op_base ++;
my $base_t = 35;
$base_t = $ARGV[$op_base] if(defined $ARGV[$op_base]);
$op_base ++;
my $lowqual_num = 20;
$lowqual_num = $ARGV[$op_base] if(defined $ARGV[$op_base]);
$op_base ++;
my $std_times = 4;
$std_times = $ARGV[$op_base] if(defined $ARGV[$op_base]);
$op_base ++;
my $min_SV = 1000;
$min_SV = $ARGV[$op_base] if(defined $ARGV[$op_base]);
$op_base ++;
my $min_clipped = 5;
$min_clipped = $ARGV[$op_base] if(defined $ARGV[$op_base]);


# estimate min insert size as a threshold
# read the same region for all bams
my ($upper, $readlen);
if($round == 1 || $round == 2){
    ($upper, $readlen) = &estimate_upper($bam, $stat_region);
}

my ($ref_string, $ref_start, $ref_end) = ("", -1, -1);
# fetch discordant

my $min_size = $min_SV;
$min_size = $upper if($upper > $min_SV);
# Note: round 1, 2 and 3 should be in strict order.
if($round == 1){
    # extract variant read names
    my $tmp_file = "tmp.$chr.txt";
    &extract_discordant_first_round($bam, $out_dir, $min_size, $readlen, $min_clipped, $chr, $tmp_file);
}
elsif($round == 2){
    # extract intra-variant reads
    my $tmp_file = "tmp.$chr.txt";
    my $tmp_inter_file = "tmp.inter.$chr.txt";
    my $n = &extract_discordant_second_round($bam, $out_dir, $chr, $tmp_file, $min_size, $tmp_inter_file);
    print "There are in all $n discordant reads in $chr.\n";
}
elsif($round == 3){
    # extract inter-variant reads
    my $tmp_inter_file = "tmp.inter.txt";
    my $n = &extract_discordant_third_round($bam, $out_dir, $tmp_inter_file, $min_size);
    print "There are in all $n discordant reads for inter-chromosomes.\n";
}


# fetch both unmapped
#my ($hi, $n2) = &extract_unmapped($bam, $out_dir, $n1, $hi);
#print "There are in all $n2 extracted reads, including both discordant and unmapped reads.\n";
1;

# This go through the bam file at a particular chromosome once, get the read names of those of interest and save to a temporary file.
sub extract_discordant_first_round{
    my ($bam, $out_dir, $min_size, $readlen, $min_clipped, $chr, $tmp_file) = @_;
    open tmp_fh, ">$out_dir/$tmp_file" or die $!;
    open IN, "samtools view $bam $chr | env x=$readlen y=$min_size z=$min_clipped perl -lane 'if(!(\$F[4] eq \"0\" || \$F[1] & 0x0100 || \$F[1] & 0x0800 || \$F[5] =~ /H/) && (((\$F[1] & 0x0004) || (\$F[1] & 0x0008)) && ! ((\$F[1] & 0x0004) && (\$F[1] & 0x0008)) || \$F[5] =~ /(\\d+)S/ && \$1 > \$ENV{z} && abs(\$F[8]) > 2*\$ENV{x} || (\$F[5] =~ /(\\d+)[DI]/ && \$1 > \$ENV{z} && abs(\$F[8]) > 2*\$ENV{x}) || \$F[6] ne \"=\" || abs(\$F[8]) > \$ENV{y} || ((\$F[1] & 0x0010) && (\$F[1] & 0x0020) || !(\$F[1] & 0x0010) && !(\$F[1] & 0x0020)) && abs(\$F[8]) > \$ENV{y})){print \$F[0]}' | " or die $!;
    while(<IN>){
        print tmp_fh $_;
    }
    close IN;
    close tmp_fh;
}

sub extract_discordant_third_round{
    my ($bam, $out_dir, $tmp_inter_file) = @_;
    $tmp_inter_file = $out_dir . "/" . $tmp_inter_file;
    my $h;
    open OUT, ">>$out_dir/discordant.inter.sam";
    open OUT_index, ">>$out_dir/index.inter.txt";
    open IN, "samtools view $bam | env x=$tmp_inter_file perl -lane 'BEGIN{open FH_, \"<\$ENV{x}\"; while(<FH_>){chomp; \$h->{\$_}=1};close FH_}print \$_ if(defined \$h->{\$F[0]})' |" or die $!;
    my $n = 0;
    while(<IN>){
        my @a = split(/\t/, $_);
        next if($a[5] =~ /H/ || $a[1] & 0x0100 || $a[1] & 0x0800);
        push @{$h->{$a[0]}}, $_;
        if(scalar(@{$h->{$a[0]}}) == 2){
            if(&check_qualify($h->{$a[0]}, $base_t, $lowqual_num, $min_clipped)){
                print OUT @{$h->{$a[0]}};
                print OUT_index join("\t", $a[0], $n) . "\n";
                $n ++;
                delete $h->{$a[0]};
            }
        }
    }
    close IN;
    close OUT;
    close OUT_index;
    return $n;
}


sub extract_discordant_second_round{
    my ($bam, $out_dir, $chr, $tmp_file, $min_size, $tmp_inter_file) = @_;
    $tmp_file = $out_dir . "/" . $tmp_file;
    # read ref chr size to avoid overflow
    my ($ref_chr_size) = &read_ref_chr_size($bam);
    # hash for sorted (paired end reads are consecutive) reads
    my $h;
    open OUT, ">$out_dir/discordant.$chr.sam";
    open OUT_index, ">$out_dir/index.$chr.txt";
    open OUT_drp_insert, ">$out_dir/drp_insert.$chr.txt";
    open OUT_sp, ">$out_dir/sp.$chr.txt";
    open IN, "samtools view $bam $chr | env x=$tmp_file perl -lane 'BEGIN{open FH_, \"<\$ENV{x}\"; while(<FH_>){chomp;\$h->{\$_}=1};close FH_}print \$_ if(defined \$h->{\$F[0]})' |" or die $!;
    my $n = 0;
    while(<IN>){
        my @a = split(/\t/, $_);
        next if($a[5] =~ /H/ || $a[1] & 0x0100 || $a[1] & 0x0800);
        push @{$h->{$a[0]}}, $_;
        if(scalar(@{$h->{$a[0]}}) == 2){
            # print the read pair, delete in hash;
            if(&check_qualify($h->{$a[0]}, $base_t, $lowqual_num, $min_clipped)){
                chomp $h->{$a[0]}->[0];
                chomp $h->{$a[0]}->[1];
                my @l1 = split(/\t/, $h->{$a[0]}->[0]);
                my @l2 = split(/\t/, $h->{$a[0]}->[1]);
                print OUT join("\n", join("\t", $n, @l1[1 .. $#l1]), join("\t", $n, @l2[1 .. $#l2])) . "\n"; 
                print OUT_index join("\t", $a[0], $n)."\n";
                $n ++;
                # check if drp, if sr
                # a bug found. Since for INV, the current read is the second one, it should be reverse and its partner forward, which indicate no INV. Fine for the non-chr version. 
                if( !(($a[1] & 0x0004) || ($a[1] & 0x0008) || $a[6] ne "=" || ! (($a[1] & 0x0010) && !($a[1] & 0x0020)))){
                    # possible for drp and sr
                    if(abs($a[8]) > $min_size){
                        print OUT_drp_insert join("\t", $n, abs($a[8])) . "\n";
                    }
                    else{
                        my @t1 = split(/\t/, $h->{$a[0]}->[0]);
                        my @t2 = split(/\t/, $h->{$a[0]}->[1]);
                        if($t1[5] =~ /(\d+)S/ && $1 > $min_clipped || $t2[5] =~ /(\d+)S/ && $1 > $min_clipped){
                            print OUT_sp join("\t", $n, $t1[5], $t2[5]) . "\n";
                            if($get_ref == 1){
                                # deal with ref on the clipped part to tolerate alignment errors in diploid genome due to high error profile of Pacbio
                                if($t1[3] > $ref_end - 1000 || $t2[3] > $ref_end - 1000 || $t1[3] < $ref_start + 1000 || $t2[3] < $ref_start + 1000 || $ref eq ""){
                                    # build ref_file
                                    ($ref_string, $ref_start, $ref_end) = &build_ref($t1[2], $t1[3], $t2[3], $ref, $ref_block_size, $ref_chr_size);
                                }
                                my ($ref_line1, $ref_line2) = &get_reference_allele($h->{$a[0]}->[0], $h->{$a[0]}->[1], $ref_string, $ref_start);
                                my @lr1 = split(/\t/, $ref_line1);
                                my @lr2 = split(/\t/, $ref_line2);
                                print OUT join("\n", join("\t", ($n - 1) . "_ref", @lr1[1 .. $#lr1]), join("\t", ($n - 1) . "_ref", @lr2[1 .. $#lr2])) . "\n";
                            }

                        }
                    }
                }
            }
            delete $h->{$a[0]};
        }
    }
    # save rest of the reads to inter chromosomal text
    open INT, ">$out_dir/$tmp_inter_file" or die $!;
    foreach my $key (keys %$h){
        if(&check_qualify($h->{$key}, $base_t, $lowqual_num, $min_clipped)){
            print INT $key . "\n";
        }
    }
    close INT;
    close IN;
    close OUT;
    close OUT_sp;
    close OUT_drp_insert;
    return $n;
}

sub check_qualify{
    my ($reads, $base_t, $lowqual_num, $min_clipped) = @_;
    my @readss = @$reads;
    if($debug == 1){
        if(!defined $readss[1]){
            print $readss[0] . "\n"; 
        }
    }
    my @e1;
    @e1 = split(/\t/, $readss[0]);
    my @e2; 
    @e2 = split(/\t/, $readss[1]) if(scalar(@readss) >= 2);

    if(scalar(@readss) >= 2 && base::check_basequal($e1[10], $e1[5], $base_t, $lowqual_num, $min_clipped) && base::check_basequal($e2[10], $e2[5], $base_t, $lowqual_num, $min_clipped) || scalar(@readss) == 1 && base::check_basequal($e1[10], $e1[5], $base_t, $lowqual_num, $min_clipped)){
        return 1;
    }
    return 0;
}

sub read_ref_chr_size{
    my ($bam) = @_;
    my $ref_chr_size;
    open IN, "samtools view -H $bam | " or die $!;
    while(<IN>){
        if($_ =~ /^@/){
            my @a = split(/\t/, $_);
            if($a[0] =~ /SQ/ && $a[1] =~ /SN:(\S+)$/){
                my $chr = $1;
                if($a[2] =~ /LN:(\d+)$/){
                    $ref_chr_size->{$chr} = $1;
                }
            }
        }
    }
    return $ref_chr_size;
}


# build reference block file to avoid retrieving the file
sub build_ref{
    my ($chr, $pos1, $pos2, $ref_file, $ref_block_size, $ref_chr_size) = @_;
    my ($start, $end) = (-1, -1);
    if($pos2 - $pos1 > $ref_block_size){
        $start = $pos1 - 2000;
        $end = $pos2 + 2000;
    }
    else{
        $start = $pos1 - 2000;
        $end = $start + $ref_block_size + 4000;
    }
    $start = 1 if($start < 1);
    $end = $ref_chr_size->{$chr} - 1 if($end >= $ref_chr_size->{$chr});
    my @seq = split(/\n/, `samtools faidx $ref_file $chr:$start-$end`);
    return (join("", @seq[1 .. $#seq]), $start, $end);
}


# for a read pair, one of which has split, replace the split part with reference, and return with everything else the same
sub get_reference_allele{
    my ($line1, $line2, $ref_string, $ref_start) = @_;
    my ($ref_line1, $ref_line2);
    my @a = split(/\t/, $line1);
    my @b = split(/\t/, $line2);
    $a[9] = &get_ref_read($a[5], $a[2], $a[3], $a[9], $ref_string, $ref_start);
    $b[9] = &get_ref_read($b[5], $b[2], $b[3], $b[9], $ref_string, $ref_start);
    $ref_line1 = join("\t", $a[0] . "_ref", @a[1 .. $#a]);
    $ref_line2 = join("\t", $b[0] . "_ref", @b[1 .. $#b]);
    return ($ref_line1, $ref_line2);
}


sub get_ref_seq{
    my ($chr, $p, $len, $ref_seq, $ref_start) = @_;
    if(defined $ref_seq){
        my $start = $p - $ref_start;
        return substr($ref_seq, $start, $len);
    }
    else{
        my @seq;
        my $end = $p + $len;
        @seq = split(/\n/, `samtools faidx $ref $chr:$p-$end`);
        return join("", @seq[1 .. $#seq]);
    }
}

sub get_ref_read{
    my ($cigar, $chr, $pos, $seq, $ref_seq, $ref_start) = @_;
    if($cigar =~ /S/){
        my $p = $pos;
        if($cigar =~ /^(\d+)S/){
            $p = $pos - $1;
        }
        $seq = &get_ref_seq($chr, $p, length($seq), $ref_seq, $ref_start);
    }
    return $seq;
}
# get reference seq according to chr, pos and seq (to get length)
#sum up the elements in an array
sub sum{
    my ($a) = @_;
    my $s = 0;
    foreach my $t (@$a){
        $s += $t;
    }
    return $s;
}



# estimate min insert size 
sub estimate_upper{
    my ($bam, $stat_region) = @_;
    open IN, "samtools view $bam $stat_region |" or die $!;

    my @ins = ();
    my @readlens = ();
    while (<IN>) {
        last if (scalar(@ins) == 10000);
        next if /^@/;
        my @e = split;
        next unless ($e[1] & 66);
        next if ($e[5] =~ /S/);
        next if ($e[4] < 20);
        next if ($e[8] <= 100 or $e[8] > 1000);
        push @ins, $e[8];
        push @readlens, length($e[9]);
        #print "Pushed $e[8]\n";
    }

    #my @ins1 = sort { $a <=> $b } @ins;
    my $mean;
    my $sd;
    if(scalar(@ins) != 0){
        $mean = &sum(\@ins)/scalar(@ins);
        my $sqsum = 0;
        foreach (@ins) {
            $sqsum += ($_-$mean)**2;
        }
        $sd = sqrt($sqsum/scalar(@ins));

    }
    close IN;
    my $upper = $mean + $std_times*$sd;
#my $lower = $mean - 3*$sd>0?$mean-3*$sd:0;
    my $readl = &sum(\@readlens)/scalar(@ins);
    print STDERR "mean= ", $mean, " sd=", $sd, " upper=$upper", " readlen=$readl" , "\n"; #, " lower=$lower", "\n";
    return ($upper, $readl);
}


