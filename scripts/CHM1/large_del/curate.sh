#!/usr/bin/env bash
if [ $# -lt 2 ]; then
    echo "This script takes the delly and Pacbio (generaetd in Chaisson's paper) large (>50bp) deletion calls and use Illumina bam for genotyping, generating a high confidence curated call set. The calls have been duplicate removed. The output is both.ILH.get2.rmdup. Note: samtools will be executed in this script. \nUsage: $0 <bam_IL> <delly_call> <Chaisson_call> <ref>\n";
    exit
fi

#illumina bam
bam_IL=$1

# delly calls (delly.GT50.bed)
delly=$2

# Chaisson calls
PBcall=$3

# reference
ref=$4

# curate delly calls
find_del_sr.highQualReads.pl $delly $bam_IL 500 50 $ref $delly.ILH short

# curate Chaisson calls
find_del_sr.highQualReads.pl $PBcall $bam_IL 500 50 $ref $PBcall.ILH short

# put the two together and select those with genotyping Illumina reads >= 2 (delly_Chaisson.GT50.bed.ILH.GET2)
cat $delly.ILH $PBcall.ILH | perl -ane 'print $_ if($F[$#F] >= 2)' > both.ILH.GET2

# remove duplications
rm_dup_bed1.pl both.ILH.get2 reciprocal50 > both.ILH.get2.rmdup
