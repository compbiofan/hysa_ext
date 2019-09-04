#!/usr/bin/env bash
if [ $# -lt 2 ]; then
    echo "This script takes the large insertion calls (in bed) as the input, use the WG assembly as the reference for validation (additional one columns, 1 as validated and 0 as not validated). Note: bwa, samtools and intersectBed in bedtools should be made executable to run this script. The output is \$ins.wValidation. \nUsage: $0 <ref which the short reads were aligned against for bam> <ins_bed> <whole genome assembly in fasta>"
    exit
fi

# reference fasta file (~/reference/human_g1k.v37.fasta)
ref=$1

# my output insertion file (out.ins.GT50.rmdup.bed)
ins=$2

# Phillippy's assembly (/scratch/bcb/xfan3/Database/CHM1/Philippy/GCA_000772585.3_ASM77258v3_genomic.fna)
phi=$3

# reconstruct the alternative allele
make_config_from_bed_file.pl $ins $ref $ins.fa 1000 ins

# align to Phillippy
bwa mem $phi $ins.fa > $ins.fa.toPhillippy.sam

# analyze alignment and come up with validated set
analyze_inseted_sequence.pl $ins.fa.toPhillippy.sam 0.9 50 50 > bwa_allctg2asm.analysis

perl -ane 'BEGIN{open fh_, "<bwa_allctg2asm.analysis"; while(<fh_>){next if($_ !~ /well_aligned/); @a = split(/\t/, $_); foreach $i (@a[1 .. $#a]){($c1, $b1, $e1, $c2, $b2, $e2) = ($i =~ /^(\S+):(\d+)-(\d+);(\S+):(\d+)-(\d+)$/); $h->{"$c1:$e1-$b2"} = 1; }} close fh_; } if(defined $h->{"$F[0]:$F[1]-$F[2]"}){print join("\t", @F, 1) . "\n"}else{print join("\t", @F, 0) . "\n"}' $ins > $ins.wValidByCHMASM


