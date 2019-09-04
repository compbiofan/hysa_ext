#!/usr/bin/env bash
if [ $# -lt 2 ]; then
    echo "This script takes a deletion calls (in bed) as the input, use reads in Illumina and Pacbio bam file to genotype the calls. It adds extra columns, in the order of supporting IL # (reported as significant), PB genotype details, supporting PB # (reported, though not significant), mean dist of PB start pos to the call (reported as significant), std dev of the dist of PB start pos to the call (reported as significant), mean difference of PB INDEL size to the call, std dev of the PB INDEL size to the call. Reference is the one against which the bam files were obtained. Result is in \$bed.wSRGT.wPBGT.summaryAll. Note: samtools will be executed in this script.\nUsage: $0 <bam_IL> <bam_PB> <bed> <ref>\n"; 
    exit
fi

# short and long bam files. (~/combinePBIL/short.bam and /scratch/1KGENOME/Eichler/rawPacbio/rawH5.bam)
bam_IL=$1
bam_PB=$2

# either single or unique bed files
bed=$3

# build37
ref=$4

# add the first column: # of supporting IL
find_del_sr_unmapped.pl $bed $bam_IL 500 50 $ref $bed.wSRGT short

# add the second column: details of PB supporting
find_del_pb_aln.pl $bam_PB $bed.wSRGT $bed.wSRGT.wPBGT 300 30 500 100 100 0.2

# analyze PB genotype results to create the final file
analyze_std.pl $bed.wSRGT.wPBGT 5 > $bed.wSRGT.wPBGT.summaryAll
