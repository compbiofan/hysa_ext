#!/usr/bin/perl
use warnings;
use strict;
require scaffold;

if(@ARGV == 0){
    die "If a Pacbio read appears in two contigs, these two contigs are put together for assembly. Similarly, if other contigs are put together with one of these contigs, they are all put together for assembly. After assembly, align all ILs belonging to related clusters to the contig. Cluster ILs. Align the ctg to the reference. Use infer_SV to check if all IL clusters are covered (next step: use IL to refine breakpoints). If not, extract the two segments on the left and right of the IL cluster (how to define the breakpoint is tricky, may be overlap a little), and realign them to the reference. Use infer_SV to check any SV is called with the new alignments added. Also, for smaller events, extract the segment (since it has a certain range with IL covers) with some flanking regions and align it to the reference (small INV should show up as three alignments in this way). \nUsage: $0 <ctg> <cluster_file> <IL_fa_prefix> <out>\n";
}

my ($ctg, $clu, $ordered_IL, $out) = @ARGV;
my $ref = "~/reference/build38.fa";
my $blasr_binary_new = "~/blasrn";
#    my $velveth = "~/pkg/velvet/velveth";
#    my $velvetg = "~/pkg/velvet/velvetg";

my $h;
my $h_pb;
my $cluster;
# read cluster file and see if any cluster can be merged
open fh_, "<$clu" or die $!;
while(<fh_>){
    chomp;
    if($_ =~ /^#(\d+)/){
        $cluster = $1;
    }
    elsif($_ !~ /\./){
        $h->{$cluster}->{IL}->{$_} = 1;
    }
    else{
        my @a = split(/\./, $_);
        $h_pb->{join(".", @a[0 .. $#a - 1])}->{$cluster} = 1;
        $h->{$cluster}->{PB}->{$_} = 1;
    }
}
close fh_;

my $h_c;
# cluster the cluster together
foreach my $c (keys %$h){
    $h_c->{$c} = $c if(!defined $h_c->{$c});
    foreach my $pb (keys %{$h->{$c}->{PB}}){
        my @a = split(/\./, $pb);
        my $p = join(".", @a[0 .. $#a - 1]);
        if(scalar(keys %{$h_pb->{$p}}) > 1){
            # see if this read related with other clusters
            foreach my $c_ (keys %{$h_pb->{$p}}){
                if($c_ != $c){
                    # not including myself
                    my $c_1 = scaffold::root($c_, $h_c);
                    my $c1 = scaffold::root($c, $h_c);
                    $h_c = scaffold::connect($c_1, $c1, $h_c);
                }
            }
        }
    }
}

# get the clusters from root of each node
my $h_cluster = scaffold::tree($h_c);

# read contig file
my $g;
my $ctg_n;
my $h_ctg;
open fh_, "<$ctg" or die $!;
while(<fh_>){
    chomp;
    if($_ =~ /^>g(\d+).ctg(\d+)/){
        $g = $1;
        $ctg_n = $2;
        $h_ctg->{$g}->{$ctg_n} = "";
    }
    else{
        $h_ctg->{$g}->{$ctg_n} .= $_; 
    }
}
close fh_;

# assemble big cluster, and align short reads to the assembled contig
my $num = 0;
`rm $out` if(-e $out);
foreach my $c (keys %$h_cluster){
    if(! -d "clu$num"){
        `mkdir clu$num`;
    }
    my $dir = "clu$num";
    my $file_fa = "$dir/seq.fa";
    my $file_clu = "$dir/clu.txt";
    my $assembled = 0;
    open fh_clu, ">$file_clu" or die $!;
    open fh_fa, ">$file_fa" or die $!;
    # for each contig, get all the sequences
    foreach my $ctg (keys %{$h_cluster->{$c}}){
        print fh_clu $ctg . "\n";
        # read contig fa and write to the big cluster file
        # it is possible nothing is assembled 
        foreach my $ctg_sub (keys %{$h_ctg->{$ctg}}){
            print fh_fa ">g$ctg.ctg$ctg_sub\n";
            print fh_fa $h_ctg->{$ctg}->{$ctg_sub} . "\n";
            $assembled = 1;
        }
    }
    close fh_clu;
    close fh_fa;
    if($assembled == 0){
        &clean();
        next;
    }
    my $new_file_fa = scaffold::assemble($file_fa, $dir);
    print $new_file_fa . "\n";
    # check if merged
    $DB::single = 1;
    my $ctg_names = scaffold::if_merged($new_file_fa);
    #if(scalar(@$ctg_names) != 1){
        #next;
     #}
    my $IL_align = scaffold::align_IL2ctg($new_file_fa, $h_cluster->{$c}, $h, $ordered_IL, $num, $ctg_names);
    print $IL_align . "\n";
    # cluster ILs aligned to ctg
    my $IL_clusters = scaffold::cluster($IL_align);
        $DB::single = 1;

    # check overlap
    if(scaffold::overlap($IL_clusters) == 0){
        # segment into three parts: 1) no IL at all; 2) IL with flanking region on the left and right; 3) IL with flanking region on one side and no IL sequence on the other side
        # for two cluster ctg, the following 9 segments will be taken
        #   ---------rrrrrrrrr---------rrrrrrrrrr---------
        # 1.---------
        # 2.--------------------
        # 3.       -------------
        # 4.       --------------------
        # 5.                  ---------
        # 6.                  ---------------------
        # 7.                         --------------
        # 8.                         ---------------------
        # 9.                                     ---------
        # 1, 5 and 9 are in the first set of segments; 2, 4, 6 and 8 are in teh second set of segments; 3 and 7 are in the third segments 
        # this is for all SV types
        # NOTE: on Sep. 27, 2016, changed to the following
        # couple each neighboring pair of clusters
        #   ---------rrrrrrrrr---------rrrrrrrrrr---------rrrrrrrrrrr----------
        # 1.       --------------------------------
        # 2.                         ----------------------------------
        # In case 1. is a large insertion, the breakpoint would be at
        # 1.           *--inserted sequence--*   (if novel, inserted tag; otherwise, two breakpoints)
        # In case 1. is a small event + a small event, should be aligned from end to end
        # In case 1. is a small event + a translocation (one breakpoint of insertion), call the translocation from alignment
        # In case 1. is a translocation + a small event, call the translocation
        # In case 1. is a translocation + a translocation, call the two translocations
        # Note: I only deal with translocation here. No deletion or insertion.
        foreach my $rd_nm (keys %$IL_clusters){
            my $segment_file = scaffold::segment($new_file_fa, $IL_clusters->{$rd_nm}, $rd_nm);
            my $seg_align = scaffold::align_seg2ctg($segment_file, $ref);
            # this is for all SV types
            #&analyze_align($seg_align);
            # this is only for small INV
            scaffold::analyze_align($seg_align, $out);
        }
    }


    &clean();

    $num ++;

}
1;

sub clean{
    #`rm *summary*`;
    #`rm *.*.*.*.ref.fa`;
    `rm tmp_seg.fa`;
    `rm *.IL.?.fa`;
    #`rm -r clu0`;
    #`rm tmp_seg.fa.sam`;
}
