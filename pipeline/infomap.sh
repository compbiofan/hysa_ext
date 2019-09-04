#!/usr/bin/env bash
if [ $# -lt 2 ]; then
    echo "usage: $0 <m4> <unionfind_prefix> <pb_threshold> <infomap_prefix>"
    exit
fi
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
m4=$1
unionfind_prefix=$2
max=$3
infomap_prefix=$4
infomap_dir="$( cd "$( dirname "$infomap_prefix" )" && pwd )"

# get the reads of the group out
out_txt=$infomap_dir/out.large.txt
perl $DIR/generate_group_file_max.pl $unionfind_prefix $max $out_txt
# get the m4 out
large_m4=$infomap_dir/out.m4.large
perl $DIR/get_new_m4_from_out.pl $m4 $large_m4 $out_txt 
# renumber the nodes
renumbered=$large_m4.renumbered
perl $DIR/renumber_nodes.pl $large_m4 $renumbered 
# run infomap
link=$renumbered.link
tree=$renumbered.tree
node=$renumbered.node
Infomap -i link-list $link ./ > $tree.out 
mv out.m4.large.renumbered.tree $tree 

# turns back to original node index and convert from tree to union find compatible files
if [ -e $infomap_prefix.txt ]; then
    rm $infomap_prefix.txt
fi
if [ -e $infomap_prefix.all ]; then
    rm $infomap_prefix.all
fi
perl $DIR/analyze_infomap_tree.pl $tree $node $infomap_prefix
