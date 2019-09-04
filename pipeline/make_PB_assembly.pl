#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;

if(scalar(@ARGV) == 0){
    die "This script generates the assembly commands, one for each cluster, with a separator \"#\" for every 500 lines. \nUsage: $0 <cluster_prefix> <output_dir> <pb_prefix> <tmp_dir> <time_limit> <cov: min:max> <cmd_output_file> <thread>\n";
}

my $path = dirname($0);
# change.log (11252015): add one more input argument at the end specifying the analysis directory of the trio. For YRI it is ~/d
# 03222016 change from number of core form 2 to 4, since there are many clusters need large number of core to finish within time
my ($cluster_prefix, $output_dir, $pb_prefix, $tmp_dir, $time_limit, $cov, $cmd_output_file, $thread) = @ARGV; 
my $num_per_job = 100;
my $txt = $cluster_prefix . ".all";
my $cluster_file = $cluster_prefix . ".txt";
my $bin = "$path/asm_pb.sh";
my ($cov_min, $cov_max) = split(/:/, $cov);
my $id = 0;
my $num = 0;
open OUT_fh, ">$cmd_output_file" or die $!;
open TXT, "<$txt" or die $!;
while(<TXT>){
    chomp;
    my @a = split(/\t/, $_);
    # change PB from 3 to 5, since CA can assemble with >= 5 PBs. 03222016
    #if($a[1] >= 3 && $a[2] >=3 && $a[2] < 50000){
    if($a[1] >= $cov_min && $a[2] >=5 && $a[2] < $cov_max){
        print OUT_fh "$bin $pb_prefix $tmp_dir/g$a[0] $cluster_prefix $output_dir/$id $a[0] $time_limit $thread\n";
        $num ++;
        if($num >= $num_per_job){
            $num = 0;
            $id ++;
            print OUT_fh "#\n";
        }
    }
}
close TXT;
close OUT_fh;
