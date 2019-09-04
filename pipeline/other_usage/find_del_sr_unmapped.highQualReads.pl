use warnings;
use strict;
require fa;
require cigar;

if(@ARGV == 0){
    die "Given a bed file with deletion calls, and a bam file, this script reconstructs the alternative allele, extract soft clipped reads surrounding the breakpoints and align them to the alternative allele. If there are no soft clip, the read is counted as a supporting read. Output a new bed file with the 4th column the number of soft clip supporting read of the call. Now support any read for confirmation, including Illumina and Pacbio. \nUsage: $0 <bed> <bam> <outer_flanking_region> <inner_flanking_region> <reference> <out>\n";
}

my $mode = "short";
my ($bed, $bam, $outer, $inner, $reference, $out_file, ) = @ARGV;
$mode = "long" if(defined $ARGV[6] && $ARGV[6] eq "long");
open OUT, ">$out_file" or die $!;
open fh_, "<$bed" or die $!;
my $h_record;
while(<fh_>){
    my $h;
    my $num = 0;
    chomp;
    my @a = split(/\t/, $_);
    my $l_s = $a[1] - $outer;
    my $r_s = $a[1] + $inner;
    my $l_e = $a[2] - $inner;
    my $r_e = $a[2] + $outer;
    my $x = `samtools view -c $bam $a[0]:$l_s-$r_s`;
    if($x > 50000){
        # ignore this call
        print OUT join("\t", @a, -1) . "\n";
        next;
    }
    $x = `samtools view -c $bam $a[0]:$l_e-$r_e`;
    if($x > 50000){
        print OUT join("\t", @a, -1) . "\n";
        next;
    }
    my @r1 = split(/\n/, `samtools view $bam $a[0]:$l_s-$r_s`);
    my @r2 = split(/\n/, `samtools view $bam $a[0]:$l_e-$r_e`);
    foreach (@r1){
        my @x = split(/\t/, $_);
        next if($x[4] <= 0);
        if($mode eq "short"){
            if($x[5] =~ /(\d+)S$/){
                if($1 < 10){
                    next;
                }
            }
            elsif($x[5] ne "*"){
                next;
            }
            #next if(! ($x[5] =~ /(\d+)S$/ && $1 >= 10));
            # it should not have very large number of mismatches
            next if($_ =~ /NM:i:(\d+)/ && $1 >= 3);
        }

        if($mode eq "long"){
            my $cigar = $x[5];
            # check if matchlen < 0.5 totallen
            my ($matchlen, $readlen) = split(/\t/, cigar::get_matchlen_on_read_of_all($cigar));
            if($matchlen / $readlen < 0.5){
                next;
            }
        }
        $h->{"$x[0].1.$x[5]"} = $x[9];
        $h_record->{$x[0]} = $x[1];
    }
    foreach (@r2){
        my @x = split(/\t/, $_);
        next if($x[4] <= 0);

        if($mode eq "short"){
            #next if($x[5] !~ /^\d+S/);
            if($x[5] =~ /^(\d+)S/){
                if($1 < 10){
                    next;
                }
            }
            elsif($x[5] ne "*"){
                next;
            }
            #next if(!($x[5] =~ /^(\d+)S/ && $1 >= 10));
            next if(defined $h_record->{$x[0]} && ($h_record->{$x[0]} & 0x0040) == ($x[1] & 0x0040));
            # it should not have very large number of mismatches
            next if($_ =~ /NM:i:(\d+)/ && $1 >= 3);
        }

        if($mode eq "long"){
            next if(defined $h_record->{$x[0]});
            my $cigar = $x[5];
            # check if matchlen < 0.5 totallen
            my ($matchlen, $readlen) = split(/\t/, cigar::get_matchlen_on_read_of_all($cigar));
            if($matchlen / $readlen < 0.5){
                next;
            }
        }
        $h->{"$x[0].2.$x[5]"} = $x[9];
    }
    my $bed_str = "$a[0]:$l_s-$a[1];$a[0]:$a[2]-$r_e";
    fa::config_fa_bed($bed_str, $reference, "tmp.alt.ctg.fasta");
    fa::write_fa($h, [keys %$h], "tmp.sr.fasta");
    `~/pkg/bwa/bwa index tmp.alt.ctg.fasta`;
    my @b =  split(/\n/, `~/pkg/bwa/bwa mem tmp.alt.ctg.fasta tmp.sr.fasta`);
    my $n = 0;
    foreach my $b_ (@b){
        next if($b_ =~ /^@/);
        chomp;
        my @x = split(/\t/, $b_);
        #if($x[5] !~ /[SH]/ && $x[5] =~ /M/ && 501 - $x[3] < 102 && $x[3] < 501){
        # instead of checking by fixed readlength 102, now check whether the read overlap with the breakpoint, as well as check whether there is soft clip at the end while the ctg is still overhanging
        # criteria:
        # 1. the breakpoint 500 is within the matched (not in S or H) string
        # 2. FOr short reads, no SH
        # 3. For long reads, check if within 450:550, > 80% percentage of identity w.r.t. both ref and reads. Also check for the entire ctg, whether > 80% percentage of identity w.r.t. ref.  
        if(&check_reads($b_, $mode)){
            $n ++;
        }
    }
    print OUT join("\t", @a, $n) . "\n";
}
close fh_;
close OUT;

1;

sub check_reads{
    my ($line, $mode) = @_;
    my @x = split(/\t/, $line);
    my $cigar = $x[5];
    if(! cigar::check_crossing_bp($line, 500)){
        return 0;
    }

    if($mode eq "short"){
        # short
        if($cigar =~ /[SH]/){
            return 0;
        }
        else{
            return 1;
        }
    }
    else{
        # long
        my $matchlen = cigar::get_matchlen_on_read_of_all($cigar);
        my $c1 = cigar::test_perc_identity("450:550", $x[3], $cigar, 0.8, "I");
        my $c2 = cigar::test_perc_identity("450:550", $x[3], $cigar, 0.8, "D");
        my $c3 = cigar::test_perc_identity("10:990", $x[3], $cigar, 0.8, "D");
        my $c4 = cigar::test_perc_identity("10:990", $x[3], $cigar, 0.8, "I");
        if($c1 == 1 && $c2 == 1 && $c3 == 1 && $c4 == 1){
            return 1;
        }
        else{
            return 0;
        }
    }
    return 0;
}



