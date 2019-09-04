#!/usr/bin/env bash
if [ $# -lt 2 ]; then
    echo "Given two files, one from delly and PB's curated calls, the other from hysa's unique validated calls, prepare files for comparison of the # of supporting reads and PB gaps' distribution in R. \nUsage: $0 <curate_call> <hysa_uniq_call>\n";
    exit
fi


f=$1
#unique.bed.wSRGT.wPBGT.summaryAll
perl -ane 'print $F[4] . "\n"' $f > ILsupport_single.txt
perl -ane '$s = 0; @a = split(/\;/, $F[5]); @b = split(/\,/, $a[0]); foreach $b_ (@b){$s+= abs($b_)}; $s/=scalar(@b); print $s . "\n"' $f > PBstart_mean_single.txt
perl -ane '$s = 0; @a = split(/\;/, $F[5]); @b = split(/\,/, $a[0]); foreach $b_ (@b){$s+= abs($b_)}; $s/=scalar(@b); foreach $b_ (@b){$d += ($s - abs($b_))**2} $d/=scalar(@b); print sqrt($d) . "\n"' $f > PBstart_sd_single.txt

f=$2
#unique.bed.wSRGT.wPBGT.summaryAll
perl -ane 'print $F[4] . "\n"' $f > ILsupport_unique.txt
perl -ane '$s = 0; @a = split(/\;/, $F[5]); @b = split(/\,/, $a[0]); foreach $b_ (@b){$s+= abs($b_)}; $s/=scalar(@b); print $s . "\n"' $f > PBstart_mean_unique.txt
perl -ane '$s = 0; @a = split(/\;/, $F[5]); @b = split(/\,/, $a[0]); foreach $b_ (@b){$s+= abs($b_)}; $s/=scalar(@b); foreach $b_ (@b){$d += ($s - abs($b_))**2} $d/=scalar(@b); print sqrt($d) . "\n"' $f > PBstart_sd_unique.txt
