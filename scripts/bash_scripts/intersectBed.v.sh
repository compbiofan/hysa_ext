#!/usr/bin/env bash
if [ $# -lt 2 ]; then
    echo "Launch intersectBed (entries in file A but not in file B) with reciprocal50 overlap in bedtools. If the 3rd argument (1-based) is 1, just output the number of quanlifying entries. Otherwise output the entries. Note: intersectBed will be executed in this script. \nUsage: $0 <fileA> <fileB> <print_mode>\n";
    exit
fi


a=$1
b=$2
count=$3
c=0.5
if [ $# -gt 3 ]; then
    c=$4
fi

if [ $count = 1 ]; then
    intersectBed -a $a -b $b -v -r -f $c | sort -n -k 2 | uniq | wc -l
else
    intersectBed -a $a -b $b -v -r -f $c 
fi
