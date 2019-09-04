#!/usr/bin/env bash
if [ $# -lt 2 ]; then
    echo "Launch intersectBed (entries in file A but not in file B) with 1bp overlap in bedtools. If the 3rd argument (1-based) is 1, just output the number of entries in file A overlapping with file B. Otherwise output the entries. Note: intersectBed will be executed in this script. \nUsage: $0 <fileA> <fileB> <print_mode>\n";
    exit
fi


a=$1
b=$2
count=$3

if [ $count = 1 ]; then
    intersectBed -a $a -b $b -v | wc -l
else
    intersectBed -a $a -b $b -v 
fi
