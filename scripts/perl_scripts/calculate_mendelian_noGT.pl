use warnings;
use strict;

if(scalar(@ARGV) == 0){
    die "This calculates venn diagram for trio data given father (the first bed file given), mother (the second bed file given) and child (the third bed file given) bed files. If mode is touch, it use 1bp overlap scheme, otherwise 50% reciprocal overlap\nUsage: $0 <f_bed> <m_bed> <c_bed> <mode>\n";
}

my ($f, $m, $c, $mode) = @ARGV;

my $com = "intersectBed.reciprocal50.sh";
my $com_v = "intersectBed.v.sh";

if($mode eq "touch"){
    $com = "intersectBed.touch.sh";
    $com_v = "intersectBed.v1.sh";
}

my $verbose = 0;
print "father and mother:\n" if($verbose == 1);
my $fm = `$com $f $m 1`; #| sort -u -k 1,1 -k 2,2 -k 3,3 | wc -l`;
chomp $fm;
print "father and mother union:\n" if($verbose == 1);
`$com $f $m > $f.parents`;
print "father not mother:\n" if($verbose == 1);
my $f_m = `$com_v $f $m 1`;
chomp $f_m;
print "mother not father:\n" if($verbose == 1);
my $m_f = `$com_v $m $f 1`;
chomp $m_f;
print "child and parent union:\n" if($verbose == 1);
my $fmc = `$com $c $f.parents 1`;# | sort -u -k 1,1 -k 2,2 -k 3,3 | wc -l`;
chomp $fmc;
print "child and father:\n" if($verbose == 1);
my $fc = `$com $c $f 1`;# | sort -u -k 1,1 -k 2,2 -k 3,3 | wc -l`;
chomp $fc;
print "child and mother:\n" if($verbose == 1);
my $mc = `$com $c $m 1`;#| sort -u -k 1,1 -k 2,2 -k 3,3 | wc -l`;
chomp $mc;
print "child not father:\n" if($verbose == 1);
my $c_f = `$com_v $c $f 1`;
chomp $c_f;

`cat $f $m | cut -f1,2,3 > $f.unionParents`;
print "child not parent union:\n" if($verbose == 1);
`$com_v $c $f.unionParents > $c.only`;
my $c_fm = `$com_v $c $f.unionParents 1`;
my $fc_m = $fc - $fmc;
my $f_mc = $f_m - $fc_m;
my $fm_c = $fm - $fmc;
my $mc_f = $mc - $fmc;
my $m_fc = $m_f - $mc_f;
#my $c_fm = $c_f - $mc_f;

print join("\n", "F_only: $f_mc", "M_only: $m_fc", "C_only: $c_fm", "FM_only: $fm_c", "FC_only: $fc_m", "MC_only: $mc_f", "all: $fmc") ."\n";
    

