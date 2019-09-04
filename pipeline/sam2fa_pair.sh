#!/usr/bin/env bash
if [ $# -lt 1 ]; then
    echo "usage: $0 <sam_prefix>" 
    exit
fi
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
perl $DIR/sam2fa_pair.pl $1.sam
#perl $DIR/fq2fa.pl ${1}_1 > ${1}_1.fa
#perl $DIR/fq2fa.pl ${1}_2 > ${1}_2.fa
