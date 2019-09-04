#!/usr/bin/perl
use warnings;
use strict;
require cluster;

if($#ARGV == -1){
    die "This checks the ils on pb, if the reads overlap with each other, and size is greater than cluster_s (it could be a range min:max), then report this relation in m4. Size_t is the maximum cluster length from beginning to end of the read.\nUsage: $0 <m4> <cluster_s> <out>\n";
}
cluster::smart_cluster($ARGV[0], 0, @ARGV[1 .. 2]);
