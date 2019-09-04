#!/usr/bin/env bash
if [ $# -lt 1 ]; then
    echo "This take a sam file and make the sorted bam from it. \nUsage: $0 <sam_prefix>\n";
    exit
fi
# given a sam, make a sorted bam
samtools view -bhS $1.sam > $1.bam
#if [ -z "$2" ]; 
#then
#    samtools sort -@ $2 $1.bam $1.sorted
#else
samtools sort $1.bam $1.sorted
#fi

samtools index $1.sorted.bam
