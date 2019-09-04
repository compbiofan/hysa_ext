#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;
use File::Basename;
# Instead of submitting jobs, this script runs blasr and filter directly in a sequential order.
#if($#ARGV < 1){
#    die "Function: This file runs the blasr alignment from short read pairs on fasta format (_1.fa,_2.fa) to pacbio file. It examines if any normal alignment and print them. It outputs the m4 format.Ready for further filtering. Usage: $0 <short_read_prefix> <PB_read.fa> <output_dir> <bestn> <suffix_sam_or_m4> <insert size threshold>\n";
#}

my %opts = (b=>"~/blasrn", m=>4, M=>13, j=>40, c=>40, s=>0, p=>12, l=>0, f=>"m4", o=>"blasr", i=>800, P=>70, L=>70, n=>40, F=>3);
getopts('b:m:M:i:c:s:p:f:F:i:P:L:l:o:n:u', \%opts);
die("
Usage:   runBLASR.bigPipe.pl <short_read_prefix> <PB_read_prefix> 
Options:
        -b  STR     BLASR binary

    BLASR options:
        -m  INT     minMatch [$opts{m}]
        -M  INT     maxMatch [$opts{M}]
        -j  INT     minInterval [$opts{j}]
        -c  INT     nCandidates [$opts{c}]
        -s  INT     maxScore (<=0) [$opts{s}]
        -p  INT     nproc (>=1) [$opts{p}]
        -f  STR     output format (m4 or sam) [$opts{f}]
        -n  INT     bestn [$opts{n}]

    Pairing/Filtering options:
        -F  INT     filtering stage. 0 if no pairing, -P and -L applied to each read separately; 1 if pairing but with no orientation and insert size requirement; 2 if pairing and with orientation requirement; 3 if pairing and with both orientation and insert size -i. [$opts{F}]
        -i  INT     insert size threshold (>0) [$opts{i}]
        -P  INT     percentage of identity threshold (<100) [$opts{P}]
        -L  INT     alignment length threshold [$opts{L}]

    Pipeline options
        -l  INT     level to start the pipeline. 1 if from submitting blasr jobs; 2 if from matching pairs; 3 if from clustering; 4 if from assembly; 5 if from aligning contigs; 6 if calling variants [$opts{l}]
        -o  STR     output directory [$opts{o}]
        -u          unpaired; by default off
\n
") unless (@ARGV);

my ($short_prefix, $PB_prefix) = @ARGV;

my $blasr = $opts{b};
my $minMatch = $opts{m};
my $maxMatch = $opts{M};
my $minInterval = $opts{j};
my $nCandidates = $opts{c};
my $maxScore = $opts{s};
my $nproc = $opts{p};
my $level = $opts{l};
my $suffix = $opts{f};
my $output_dir = $opts{o};
my $insert_t = $opts{i};
my $perc_t = $opts{P};
my $alignLen_t = $opts{L};
my $filterL = $opts{F};
my $bestn = $opts{n};
my $unpaired = 0;
$unpaired = 1 if($opts{u});
#my $option = join("_", join(":", "-minInterval", $minInterval), join(":", "-nCandidates", $nCandidates), join(":", "-maxScore", $maxScore), join(":", "-minMatch", $minMatch), join(":", "-maxMatch", $maxMatch), join(":", "-nproc", $nproc), join(":", "-bestn", $bestn));
my $option = join("_", join(":", "-maxScore", $maxScore), join(":", "-minMatch", $minMatch), join(":", "-maxMatch", $maxMatch), join(":", "-nproc", $nproc), join(":", "-bestn", $bestn));
#my $option = "-minInterval $minInterval -nCandidates $nCandidates -maxScore $maxScore -minMatch $minMatch -maxMatch $maxMatch -nproc $nproc -bestn $bestn";
#my $option = "-maxMatch 11  -minMatch 4 -indelRate 0 -sdpTupleSize 5 -extend -maxExtendDropoff 40";
my $path = dirname($0);

# beginning
`mkdir $output_dir` if(!-d $output_dir);

if($unpaired == 1){
    my $blasr_jobIDs = "NA";
    if($level <= 1){
        # run BLASR for a single file alone
        $blasr_jobIDs = &runBLASR_single($short_prefix, $PB_prefix, $output_dir);
    }
    print $blasr_jobIDs . "\n";

    my $match_jobID = "NA";
# run pairing and filtering
    if($level <= 2){
        $blasr_jobIDs = "" if($level == 2);
        $match_jobID = &runFilter_single($output_dir, $blasr_jobIDs, $suffix);
    }

}
else{
# run BLASR
    my $blasr_jobIDs = "NA";
    if($level <= 1){
        # run BLASR for all PBs suppose the PB file is local
        $blasr_jobIDs = &runBLASR($short_prefix, $PB_prefix, $output_dir);
    }

    my $match_jobID = "NA";
# run pairing and filtering
    if($level <= 2){
        $blasr_jobIDs = "" if($level == 2);
        $match_jobID = &runMatch($output_dir, $blasr_jobIDs, $suffix);
    }
}

sub runMatch{
    # now only supports m4
    my ($output_dir, $blasr_jobIDs, $suffix) = @_;
    my $dep = "";
    if($blasr_jobIDs =~ /\d+:\d+/){
        $dep = "-W x=depend=afterok:$blasr_jobIDs";
    }
    `perl $path/runFilter.pl $output_dir $suffix $perc_t $alignLen_t $insert_t $short_prefix`;
    #`gen_jb_sub.pl "perl $path/runFilter.pl $output_dir $suffix $perc_t $alignLen_t $insert_t $short_prefix" filter 1 36:00 40gb`;
    #my $match_jobID = `msub $dep filter.PBS`;
    #return $match_jobID;
    return 1;
}

sub runBLASR{
    my ($short_prefix, $PB_prefix, $output_dir) = @_;
    my $PB_fa = $PB_prefix . ".fa";
    if(! -e $PB_fa && -e $PB_fa . "sta"){
        $PB_fa .= "sta";
    }
    my $short;
    if(-e "${short_prefix}_1.fa"){
        $short = $short_prefix . "_1.fa";
    }
    elsif(-e "${short_prefix}.1.fa"){
        $short = $short_prefix . ".1.fa";
    }
    elsif(-e "${short_prefix}_1.fq"){
        `perl $path/fq2fa.pl ${short_prefix}_1 > ${short_prefix}_1.fa`;
        `perl $path/fq2fa.pl ${short_prefix}_2 > ${short_prefix}_2.fa`;
        $short = $short_prefix . "_1.fa";
    }
    elsif(-e "${short_prefix}.sam"){
        `perl $path/sam2fa_pair.pl ${short_prefix}`;
        $short = $short_prefix . "_1.fa";
    }
    if(! -e "${short_prefix}_1.fa" && ! -e "${short_prefix}.1.fa"){
        die "Not any file related with $short_prefix existed.\n";
    }

    my $outputPrefix = $output_dir . "/blasr_1";
    #`gen_jb_sub.pl "perl $path/runBlasr.pl $short $PB_fa $blasr $outputPrefix $suffix $nproc $option" blasr.1 $nproc 48:00 40gb`;
    `perl $path/runBlasr.pl $short $PB_fa $blasr $outputPrefix $suffix $nproc $option`;
    print "perl $path/runBlasr.pl $short $PB_fa $blasr $outputPrefix $suffix $nproc $option\n";
    #my $tmp1 = `msub blasr.1.PBS`;
    if(-e "${short_prefix}_2.fa"){
        $short = $short_prefix . "_2.fa";
    }
    elsif(-e "${short_prefix}.2.fa"){
        $short = $short_prefix . ".2.fa";
    }
    else{
        die "$short_prefix on the second read does not exist.\n";
    }
    $outputPrefix = $output_dir . "/blasr_2";
    #`gen_jb_sub.pl "perl $path/runBlasr.pl $short $PB_fa $blasr $outputPrefix $suffix $nproc $option" blasr.2 $nproc 48:00 40gb`;
    `perl $path/runBlasr.pl $short $PB_fa $blasr $outputPrefix $suffix $nproc $option`;
    print "perl $path/runBlasr.pl $short $PB_fa $blasr $outputPrefix $suffix $nproc $option\n";
    #my $tmp2 = `msub blasr.2.PBS`;
    #my @tmps = split(/\n/, $tmp1 . $tmp2);
    #my @tmps2;
    #foreach (@tmps){
    #    next if($_ =~ /^$/);
    #    chomp;
    #    push @tmps2, $_;
    #}
    #my $tmp = join(":", @tmps2);

    #return $tmp;
    return 1;
}

sub runFilter_single{
    # now only supports m4
    my ($output_dir, $blasr_jobIDs, $suffix) = @_;
    my $dep = "";
    if($blasr_jobIDs =~ /\d+:\d+/){
        $dep = "-W x=depend=afterok:$blasr_jobIDs";
    }
    #`gen_jb_sub.pl "perl $path/runFilter_single.pl $output_dir $suffix $perc_t $alignLen_t $short_prefix $PB_prefix.fa" filter 1 36:00 40gb`;
    `perl $path/runFilter_single.pl $output_dir $suffix $perc_t $alignLen_t $short_prefix $PB_prefix.fa`;
    #my $match_jobID = `msub $dep filter.PBS`;
    #return $match_jobID;
    return 1;
}
sub runBLASR_single{
    my ($short, $PB_prefix, $output_dir) = @_;
    my $PB_fa = $PB_prefix . ".fa";
    my $outputPrefix = $output_dir . "/blasr";
    #`gen_jb_sub.pl "perl $path/runBlasr.pl $short $PB_fa $blasr $outputPrefix $suffix $nproc $option" blasr $nproc 48:00 40gb`;
    `perl $path/runBlasr.pl $short $PB_fa $blasr $outputPrefix $suffix $nproc $option`;
    #my $tmp = `msub blasr.PBS`;
    #my @tmps = split(/\n/, $tmp);
    #my @tmps2;
    #foreach (@tmps){
    #    next if($_ =~ /^$/);
    #    chomp;
    #    return $_;
    #}
    return 1;
}


