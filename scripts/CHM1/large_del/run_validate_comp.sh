#!/usr/bin/env bash
if [ $# -lt 2 ]; then
    echo "This script takes the large deletion calls (in bed) as the input, use the WG assembly as the reference for validation (additional one columns, 1 as validated and 0 as not validated), use the short read bam for genotyping (additional one column, 1 as genotyped, 0 no), overlap with curated set by reciprocal 50%. Note: bwa, samtools and intersectBed in bedtools should be made executable to run this script. The output is \$del.wValidByCHMASM.wSRGT.overlapWGT. \nUsage: $0 <ref which the short reads were aligned against for bam> <del_bed> <whole genome assembly in fasta> <curated bed> <bam>"
    exit
fi

# reference fasta file (~/reference/human_g1k.v37.fasta)
ref=$1

# my output deletion file (out.del.GT50.rmdup.bed)
del=$2

# Phillippy's assembly (/scratch/bcb/xfan3/Database/CHM1/Philippy/GCA_000772585.3_ASM77258v3_genomic.fna)
phi=$3

# curated bed (delly_Chaisson.GT50.bed.ILH.GET2.rmdup)
curate=$4

# bam file (~/combinePBIL/short.bam)
bam=$5

# reconstruct the alternative allele
make_config_from_bed_file.pl $del $ref $del.fa 1000 del

# align to Phillippy
bwa mem $phi $del.fa > $del.fa.toPhillippy.sam

# analyze alignment and come up with validated set
analyze_deleted_sequence.pl $del.fa.toPhillippy.sam 0.9 50 50 > bwa_allctg2asm.analysis

perl -ane 'BEGIN{open fh_, "<bwa_allctg2asm.analysis"; while(<fh_>){next if($_ !~ /well_aligned/); @a = split(/\t/, $_); foreach $i (@a[1 .. $#a]){($c1, $b1, $e1, $c2, $b2, $e2) = ($i =~ /^(\S+):(\d+)-(\d+);(\S+):(\d+)-(\d+)$/); $h->{"$c1:$e1-$b2"} = 1; }} close fh_; } if(defined $h->{"$F[0]:$F[1]-$F[2]"}){print join("\t", @F, 1) . "\n"}else{print join("\t", @F, 0) . "\n"}' $del > $del.wValidByCHMASM

find_del_sr_unmapped.pl $del $bam 500 50 $ref $del.wValidByCHMASM.wSRGT short

#compare with those from one signal
overlap_add_one_col.pl $del.wValidByCHMASM.wSRGT $curate $del.wValidByCHMASM.wSRGT.overlapWGT reciprocal50
