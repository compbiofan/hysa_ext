#!/usr/bin/env bash
if [ $# -lt 2 ]; then
    echo "usage: $0 <m4> <unionfind_prefix> <group>"
    exit
fi
m4=$1
unionfind_prefix=$2
group=$3

# get the reads of the group out
#perl ~/u/generate_group_file.pl $unionfind_prefix $group out.$group.txt
# get the m4 out
#perl ~/u/get_new_m4_from_out.pl $m4 $m4.$group out.$group.txt 
# renumber the nodes
#perl ~/u/renumber_nodes.pl $m4.$group $m4.$group.renumbered 
# run infomap
#~/pkg/Infomap/Infomap -i link-list $m4.$group.renumbered.link ./ > $m4.$group.renumbered.tree 
# turns back to original node index and convert from tree to union find compatible files
perl ~/u/analyze_infomap_tree.pl out.all.m4.$group.renumbered.tree $m4.$group.renumbered.node out.all.m4.$group.renumbered.tree
