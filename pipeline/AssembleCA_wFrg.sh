#!/usr/bin/env bash
if [ $# -lt 2 ]; then
    echo "usage: $0 prefix threads frgs"
  exit
fi
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#source ~/.bashrc
#source ~/.profile


region=$1
threads=$2

if [  $# = 3 ]; then
    runCA -p $region -d $region ovlErrorRate=0.40 utgGraphErrorRate=0.40 cnsErrorRate=0.40 cgwErrorRate=0.40 unitigger=bogart obtErrorRate=0.30 ovlThreads=$2 $3;
elif [ $# = 4 ]; then
    runCA -p $region -d $region ovlErrorRate=0.40 utgGraphErrorRate=0.40 cnsErrorRate=0.40 cgwErrorRate=0.40 unitigger=bogart obtErrorRate=0.30 ovlThreads=$2 $3 $4;
elif [ $# = 5 ]; then
    $CA/runCA -p $region -d $region ovlErrorRate=0.40 utgGraphErrorRate=0.40 cnsErrorRate=0.40 cgwErrorRate=0.40 unitigger=bogart obtErrorRate=0.30 ovlThreads=$2 $3 $4 $5;
fi
cp $region/9-terminator/$region.ctg.fasta $region/
#~/find_putative/scripts/runBWA.sh $region/$region.ctg.fasta $threads

