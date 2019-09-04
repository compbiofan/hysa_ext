use warnings;
use strict;
package union_find;
require graph;

# to deal with reads from different groups, $cfg is in a format of each line having "m4 file, suffix to add to readname"
sub union_find_group{
    my ($cfg, $out) = @_;
    # prepare for the IDs, the first col stored in hash, the second in hash1, with the value in order not duplicated
    my $files;
    open fh_, "<$cfg" or die $!;
    while(<fh_>){
        chomp;
        my @tmp = split(/\s+/, $_);
        $files->{$tmp[0]} = $tmp[1];
    }
    close fh_;
    my $hash;
    my $hash1;
    my $hash_rev;
    my $hash1_rev;
    my $n = 0;
    my $n1 = 0;
    my $tr;
    my $sz;
    foreach my $f (keys %$files){
        open fh_, "<$f" or die $!;
        while(<fh_>){
            my @a = split(/\s+/, $_);
            # standardize readname for IL, no interval, no 1 or 2
            ($a[0]) = &standardize_rn($a[0]);

            # this is the suffix to add to the current m4 file
            $a[0] .= $files->{$f};
            if(!defined $hash->{$a[0]}){
                $hash->{$a[0]} = $n;
                $hash_rev->{$n} = $a[0];
                $tr->[$n] = $n;
                $sz->[$n] = 1;
                $n ++;
            }
            if(!defined $hash1->{$a[1]}){
                $hash1->{$a[1]} = $n1;
                $hash1_rev->{$n1} = $a[1];
                $n1 ++;
            }
        }
        close fh_;
    }
    foreach my $key (keys %$hash1){
        $hash1->{$key} += $n;
        my $tmp = $hash1->{$key};
        $tr->[$tmp] = $tmp;
        $sz->[$tmp] = 1;
    }

    # now read the connection
    foreach my $f(keys %$files){
        open fh_, "<$f" or die $!;
        while(<fh_>){
            my @a = split(/\s+/, $_);
            ($a[0]) = &standardize_rn($a[0]);
            $a[0] .= $files->{$f};
            my ($p, $q) = ($hash->{$a[0]}, $hash1->{$a[1]});
            my $p_r = &root($p, $tr);
            my $q_r = &root($q, $tr);
            if($p_r != $q_r){
                # connect the two
                ($tr, $sz) = &connect($p_r, $q_r, $tr, $sz);
            }
        }
        close fh_;
    }
    # cluster the tree to get connected components in ids
    my ($h, $h_s) = &tree_cluster($tr, $n, $hash1_rev, $hash_rev);
    # print the supporting number of reads for each cluster
    &print_stats($h_s, $out);
    # print names
    &print_names_one_f($h, $out);
}

# from a hash of PB reads with the IL reads that they contain in a bipartite graph, get the complement edges, store in a hash and retrun
# n is the number of IL pairs in this graph
sub get_complement_bipartite{
    my ($h, $ILs) = @_;
    my $h_ret;
    my $n = scalar(@$ILs);
    foreach my $PB (keys %$h){
        foreach my $il_id ( 0 .. $n - 1 ){
            my $il = $ILs->[$il_id];
            if(!defined $h->{$PB}->{$il}){
                $h_ret->{$PB}->{$il} = 1;
            }
        }
    }
    return $h_ret;
}

# get the score for each pair of PB
sub get_score_pair{
    my ($h, $p1, $p2) = @_;
    my @il1 = sort keys %{$h->{$p1}};
    my @il2 = sort keys %{$h->{$p2}};
    my ($joint_cardinal) = &get_joint_c(\@il1, \@il2);
    return $joint_cardinal/(scalar(@il1) + scalar(@il2) - $joint_cardinal);
}

# get the overlapping components from two sorted array
sub get_joint_c{
    my ($a1, $a2) = @_;
    my $i = 0;
    my $j = 0;
    my $n = 0;
    while($i < scalar(@$a1) && $j < scalar(@$a2)){
        if($a1->[$i] == $a2->[$j]){
            $n ++;
            $i ++;
            $j ++;
        }
        elsif($a1->[$i] < $a2->[$j]){
            $i ++;
        }
        else{
            $j ++;
        }
    }
    return $n;
}


# get the total score for all pairs of PB
sub get_total_score{
    my ($h, $theta) = @_;
    my $h_s;
    my @ps = sort keys %$h;
    my $p_num = scalar(@ps);
    my @scores;
    my $s = 0;
    for(my $i = 0; $i < $p_num - 1; $i ++){
        for(my $j = 0; $j < $p_num; $j ++){
            my $ss = &get_score_pair($h, $ps[$i], $ps[$j]);
            $h_s->{$i}->{$j} = 1 if($ss > $theta);
            $s += $ss;
        }
    }
    return ($h_s, 0) if($p_num <= 1);
    return ($h_s, $s/($p_num * ($p_num - 1) / 2));
}

sub read_group_file{
    my ($group_file) = @_;
    my $hash;
    my $g = "NA";
    open fh_, "<$group_file" or die $!;
    while(<fh_>){
        chomp;
        if($_ =~ /^#(\d+)/){
            $g = $1;
        }
        else{
            if($_ =~ /\./){
                push @{$hash->{$g}->{PB}}, $_;
            }
            else{
                push @{$hash->{$g}->{IL}}, $_;
            }
        }
    }
    close fh_;
    return $hash;
}



# given a m4 and a union_find out.txt file, read through them, and for each group, report 1) If the connection pass the score; 2) If not, what are the connected components, by separating them into smaller groups. Thus the output is a new out.txt and out.all file. Note 1) if group g9 needs to be separated, the new groups are named g9.1, g9.2 ... 2) Illumina reads might overlap in some of the subgroups. 
sub cluster_analysis{
    my ($m4_file, $group_file, $subset_group, $theta, $prefix) = @_;
    my $hash = &read_m4_record_relation($m4_file, "", "no_print", "PB_only");
    my $g_hash = &read_group_file($group_file);
    # read in the group index that are of interest
    my $sg_hash; 
    open sub_fh, "<$subset_group" or die $!;
    while(<sub_fh>){
        chomp;
        $sg_hash->{$_} = 1;
    }
    close sub_fh;

    # prepare for writing to group files
    open OUT_all, ">$prefix.all" or die $!;
    open OUT_txt, ">$prefix.txt" or die $!;
    open OUT_group, ">$prefix.group" or die $!;
    # do analysis on each qualified group
    foreach my $g (sort keys %$sg_hash){
        my $PBs = $g_hash->{$g}->{PB};
        my $ILs = $g_hash->{$g}->{IL};
        # This h is to be used later on
        my $h;
        foreach my $PB (@$PBs){
            my @this_ILs = split(/,/, $hash->{$PB});
            foreach my $IL (@this_ILs){
                $h->{$PB}->{$IL} = 1;
            }
        }
        my $h_ = &get_complement_bipartite($h, $ILs);
        $h = $h_;

        #h_s is the hash containing PB pairs with score greater than theta
        my ($h_s, $total_score) = &get_total_score($h, $theta);
        if($total_score > $theta){
            print "Group $g does not pass the filter.\n";
            my $cc = &get_cc($h_s);
            foreach my $index (keys %$cc){
                my $g_name = "$g.$index";
                my @tmp_PBs = keys %{$cc->{$index}};
                my $h_IL_tmp;
                foreach my $pb (@tmp_PBs){
                    my @tmp_ILs = split(/,/, $hash->{$pb});
                    foreach my $tmp_IL (@tmp_ILs){
                        $h_IL_tmp->{$tmp_IL} = 1;
                    }
                }
                my $line = join("\t", $g_name, scalar(keys %$h_IL_tmp), scalar(@tmp_PBs)) . "\n";
                print OUT_all $line;
                print OUT_txt join("\n", "#$line", $g_name, scalar(keys %$h_IL_tmp), @tmp_PBs) . "\n";
                print OUT_group $g_name . "\n";
            }
#            return (0, $cc);
        }
        else{
            print "Group $g pass the filter.\n";
            my $line = join("\t", $g, scalar(@$ILs), scalar(@$PBs)) . "\n";
            print OUT_all $line;
            print OUT_txt join("\n", "#$line", @$ILs, @$PBs) . "\n";
            print OUT_group $g . "\n";
#        return (1);
        }
    }
    close OUT_all;
    close OUT_txt;
    close OUT_group;
}

# given a hash representing edges in a graph, do a width first search to get connected component
sub get_cc{
    my ($h) = @_;
    # a hash to mark which has been searched through
    my $h_m;
    my $total_set;
    my $i = 1;
    foreach my $s (sort keys %$h){
        if(! defined $h_m->{$s}){
            $h_m->{$s} = 1;
            my @set;
            my $set1;
            ($set1, $h_m) = &depth_first_search($h, $s, $h_m, \@set);
            $total_set->{$i} = $set1;
            $i ++;
        }
    }
    return $total_set;
}

sub depth_first_search{
    my ($h, $s, $h_m, $set) = @_;
    foreach my $c (keys %{$h->{$s}}){
        if(!defined $h_m->{$c}){
            $h_m->{$c} = 1;
            ($set, $h_m) = &depth_first_search($h, $c, $h_m, $set);
            push @$set, $c;
        }
    }
    return($set, $h_m);
}


sub union_find_from_tree{
    my ($g, $out) = @_;
    # prepare for the IDs, the first col stored in hash, the second in hash1, with the value in order not duplicated
    my $hash;
    my $hash1;
    my $hash_rev;
    my $hash1_rev;
    my $n = 0;
    my $n1 = 0;
    my $tr;
    my $sz;

    my $repeat_max = 80;
    my $count_h;
    foreach my $k1 (keys %$g){
        foreach my $k2 (keys %{$g->{$k1}}){
            next if($k1 =~ /\./);
            my $il = $k1;
            my $pb = $k2;
            if(!defined $hash->{$il}){
                $hash->{$il} = $n;
                $hash_rev->{$n} = $il;
                $tr->[$n] = $n;
                $sz->[$n] = 1;
                $n ++;
                $count_h->{$il} = 0;
            }
            else{
                $count_h->{$il} ++;
                print "$il\n" if($count_h->{$il} > $repeat_max);
            }
            if(!defined $hash1->{$pb}){
                $hash1->{$pb} = $n1;
                $hash1_rev->{$n1} = $pb;
                $n1 ++;
            }
        }
    }
    foreach my $key (keys %$hash1){
        $hash1->{$key} += $n;
        my $tmp = $hash1->{$key};
        $tr->[$tmp] = $tmp;
        $sz->[$tmp] = 1;
    }


    foreach my $k1 (keys %$g){
        foreach my $k2 (keys %{$g->{$k1}}){
            my $il = $k1;
            my $pb = $k2;
            if($k1 =~ /\./){
                $il = $k2;
                $pb = $k1;
            }
            my ($p, $q) = ($hash->{$il}, $hash1->{$pb});
            my $p_r = &root($p, $tr);
            my $q_r = &root($q, $tr);
            if($p_r != $q_r){
                # connect the two
                ($tr, $sz) = &connect($p_r, $q_r, $tr, $sz);
            }
        }
    }
# cluster the tree to get connected components in ids
    my ($h, $h_s) = &tree_cluster($tr, $n, $hash1_rev, $hash_rev);
# print the supporting number of reads for each cluster
#&print_stats($h_s, $out);
# print names
    &print_names_one_f($h, $out);
}



# To implement fast union find based on a file with the first column connecting with the second column

sub union_find{
    my ($f, $out) = @_;
    # prepare for the IDs, the first col stored in hash, the second in hash1, with the value in order not duplicated
    my $hash;
    my $hash1;
    my $hash_rev;
    my $hash1_rev;
    my $n = 0;
    my $n1 = 0;
    my $tr;
    my $sz;

    my $repeat_max = 80;
    my $count_h;
    open fh_, "<$f" or die $!;
    while(<fh_>){
        my @a = split(/\s+/, $_);
        # standardize readname for IL, no interval, no 1 or 2
        ($a[0]) = &standardize_rn($a[0]);

        if(!defined $hash->{$a[0]}){
            $hash->{$a[0]} = $n;
            $hash_rev->{$n} = $a[0];
            $tr->[$n] = $n;
            $sz->[$n] = 1;
            $n ++;
            $count_h->{$a[0]} = 0;
        }
        else{
            # already in
            $count_h->{$a[0]} ++;
            print "$a[0]\n" if($count_h->{$a[0]} > $repeat_max);
        }
        if(!defined $hash1->{$a[1]}){
            $hash1->{$a[1]} = $n1;
            $hash1_rev->{$n1} = $a[1];
            $n1 ++;
        }
    }
    close fh_;
    foreach my $key (keys %$hash1){
        $hash1->{$key} += $n;
        my $tmp = $hash1->{$key};
        $tr->[$tmp] = $tmp;
        $sz->[$tmp] = 1;
    }

    # now read the connection
    open fh_, "<$f" or die $!;
    while(<fh_>){
        my @a = split(/\s+/, $_);
        ($a[0]) = &standardize_rn($a[0]);
        my ($p, $q) = ($hash->{$a[0]}, $hash1->{$a[1]});
        my $p_r = &root($p, $tr);
        my $q_r = &root($q, $tr);
        if($p_r != $q_r){
            # connect the two
            ($tr, $sz) = &connect($p_r, $q_r, $tr, $sz);
        }
    }
    close fh_;
    # cluster the tree to get connected components in ids
    my ($h, $h_s) = &tree_cluster($tr, $n, $hash1_rev, $hash_rev);
    # print the supporting number of reads for each cluster
    #&print_stats($h_s, $out);
    # print names
    &print_names_one_f($h, $out);
}

# standardize readname
sub standardize_rn{
    my ($x) = @_;
    if($x =~ /^(.+)\/\d+_\d+[_\/][12]$/){
        $x = $1;
    }
    elsif($x =~ /^(.+)[_\/][12]$/){
        $x = $1;
    }
    return $x;
}
sub print_stats {
    my ($h, $out) = @_;
    my $new_f = $out . ".statPB.txt";
    open fh_, ">$new_f" or die $!;
    foreach my $i (sort {$a <=> $b } keys %{$h->{P}}){
        print fh_ $h->{P}->{$i} . "\n";
    }
    close fh_;
    $new_f = $out . ".statIL.txt";
    open fh_, ">$new_f" or die $!;
    foreach my $i (sort {$a <=> $b } keys %{$h->{I}}){
        print fh_ $h->{I}->{$i} . "\n";
    }
    close fh_;
    return;
}

# print the names for each cluster in one big file out.txt, and the number of supporting IL and PB reads to out.all
sub print_names_one_f {
    my ($h, $out) = @_;
    my $num = 0;
    open fh_all, ">$out.all" or die $!; 
    open g_all, ">$out.txt" or die $!;
    foreach my $i (sort {$a <=> $b } keys %$h){
        my $str = $h->{$i};
        my @strs = split(/\n/, $str);
        my $il = 0;
        my $pb = 0;
        foreach my $s (@strs){
            # for trios: IL: 38_0
            if($s =~ /^\d+$/ || $s =~ /^\d+\_\d+$/){
                $il ++;
            }
            elsif($s =~ /\./){
                $pb ++;
            }
        }
        print fh_all join("\t", $num, $il, $pb) . "\n"; 
        print g_all "#$num\t$il\t$pb\n";
        print g_all $str . "\n";
        $num ++;
    }
    close fh_all;
    close g_all;
}

# print the names for each cluster, *.g$num.txt, and the number of supporting IL and PB reads to out_all
sub print_names {
    my ($h, $out) = @_;
    my $num = 0;
    open fh_all, ">$out.all" or die $!; 
    foreach my $i (sort {$a <=> $b } keys %$h){
        my $str = $h->{$i};
        my @strs = split(/\n/, $str);
        my $il = 0;
        my $pb = 0;
        foreach my $s (@strs){
            if($s =~ /^\d+$/){
                $il ++;
            }
            elsif($s =~ /\./){
                $pb ++;
            }
        }
        my $new_f = $out . ".g$num.txt";
        print fh_all join("\t", $new_f, $il, $pb) . "\n"; 
        $num ++;
        open fh_, ">$new_f" or die $!;
        print fh_ $str;
        close fh_;
    }
    close fh_all;
}
sub connect_hash{
    my ($p_r, $q_r, $tr, $sz, $id) = @_;
    if($sz->{$p_r} > $sz->{$q_r}){
        # q_r connected to p_r
        $tr->{$q_r} = $p_r;
        $sz->{$p_r} += $sz->{$q_r};
    }
    else{
        $tr->{$p_r} = $q_r;
        $sz->{$q_r} += $sz->{$p_r};
    }
    return ($tr, $sz);
}
sub connect{
    my ($p_r, $q_r, $tr, $sz, $id) = @_;
    if($sz->[$p_r] > $sz->[$q_r]){
        # q_r connected to p_r
        $tr->[$q_r] = $p_r;
        $sz->[$p_r] += $sz->[$q_r];
    }
    else{
        $tr->[$p_r] = $q_r;
        $sz->[$q_r] += $sz->[$p_r];
    }
    return ($tr, $sz);
}
sub root_hash{
    my ($p, $tr) = @_;
    while($tr->{$p} ne $p){
        #$tr->[$p] = $tr->[$tr->[$p]];
        $p = $tr->{$p};
    }
    return $p;
}


sub root{
    my ($p, $tr) = @_;
    while($tr->[$p] != $p){
        #$tr->[$p] = $tr->[$tr->[$p]];
        $p = $tr->[$p];
    }
    return $p;
}
sub tree_cluster_ori_hash{
    my ($tr) = @_;
    my $h;
    my $h_s;
    foreach my $t (keys %$tr){
        my $tmp = &root_hash($t, $tr);
        $h->{$tmp}->{$t} = 1;
        if(!defined $h->{$tmp}){
            $h_s->{$tmp} = 1;
        }
        else{
            $h_s->{$tmp} ++;
        }
    }
    return ($h, $h_s);
}


sub tree_cluster_ori{
    my ($tr) = @_;
    my $h;
    my $h_s;
    foreach my $t (0 .. scalar(@$tr) - 1){
        my $tmp = &root($t, $tr);
        $h->{$tmp}->{$t} = 1;
        if(!defined $h->{$tmp}){
            $h_s->{$tmp} = 1;
        }
        else{
            $h_s->{$tmp} ++;
        }
    }
    return ($h, $h_s);
}

sub tree_cluster{
    my ($tr, $n, $hash1_rev, $hash_rev) = @_;
    # hash to record the number of the roots
    my $h;
    my $h_s;
    my $name;
    foreach my $t (0 .. scalar(@$tr) - 1){
        if($t >= $n){
            my $tt = $t - $n;
            $name = $hash1_rev->{$tt};
        }
        else{
            $name = $hash_rev->{$t};
        }
        my $tmp = &root($t, $tr);
        # record the number of supporting IL and PB for furture filtering or manipulation
        if($name =~ /\./){
            $h_s->{P}->{$tmp} ++;
        }
        else{
            $h_s->{I}->{$tmp} ++;
        }
        # record the read name for printing
        if(!defined $h->{$tmp}){
            $h->{$tmp} = $name;
        }
        else{
            $h->{$tmp} .= "\n$name";
        }
    };
    return ($h, $h_s);
}
# This functions takes a file containing a set of read names, find the corresponding reads an save them to fasta files. It works for trio, so that the Pacbio reads are with the names trio_index.read_index.chomped_read_index. For example, if father, mother and child has 22, 20 and 47 tmp files, a Pacbio read from child tmp3 and read inex 4000 will have the name 44.4000.*, 44 is calculated by adding 22, 20 and 3, and -1 since it counts from 0.
# The trio_str is in the format of num1_num2_num3:str1_str2_str3. num1-3 are the numbers of the tmp files for the three samples, and the str1-3 are the string of the sample names. For example, the YRI trio's str should be 22_20_47:NA19238_NA19239_NA19240.
# Note that the long_fa_dir now should contain an unknown substring: TOREPLACE, so that according to the trio_index, the sample will be found and replace this substring.
# two frg file for short and long
sub readname2fa_dir_for_trio{
    my ($group, $short_fa_prefix, $long_fa_dir, $long_fa_name, $trio_str, $tag, $sam) = @_;
    # analyze the string
    my ($num, $str) = split(/:/, $trio_str);
    my @nums = split(/\_/, $num);
    my @strs = split(/\_/, $str);
    # this takes the readnames in group, make a new fa with both short and long reads, long read name has been hashed in long_name
    my $name;
    my $skip_il = 0;
    my $skip_pb = 0;

    if(defined $tag && $tag eq "PB_only"){
        $skip_il = 1;
    }
    elsif(defined $tag && $tag eq "IL_only"){
        $skip_pb = 1;
    }
#    my ($g) = ($group =~ /(out.g\d+\.txt)/);
    # make the output file in the same directory as out.g*.txt
    my $g = $group;
    my $out1 = $g . ".PB.fa";
    my $out2 = $g . ".IL.1.fa";
    my $out3 = $g . ".IL.2.fa";
    my $out4 = $g . ".IL.fa";
    my $out5 = $g . ".sam";
    if(defined $sam){
        open OUT_sam, ">$out5" or die $!;
    }
    if(defined $tag && $tag eq "PB_only" || !defined $tag || defined $tag && $tag eq "NA"){
        open OUT1, ">$out1" or die $!;
    }
    if(defined $tag && $tag eq "IL_only" || !defined $tag || defined $tag && $tag eq "NA"){
        open OUT2, ">$out2" or die $!;
        open OUT3, ">$out3" or die $!;
        open OUT4, ">$out4" or die $!;
    }
    open g_h, "<$group" or die $!;
    my $short_fa;
    my $sam_file;
    my $h;
    # the hash to indicate which short read file to grep
    my $short_fa_1;
    my $short_fa_2;
    my $readnum = 0;
    if($skip_il == 0){
        if($short_fa_prefix =~ /;/){
            # two source files of short reads, in the format of short1;str:short2
            my @a = split(/;/, $short_fa_prefix);
            foreach my $b (@a){
                if($b =~ /:/){
                    my @c = split(/:/, $b);
                    $h->{$c[0]} = $c[1];
                }
                else{
                    $h->{"NA"} = $b;
                    $short_fa_1 = $b . "_1.fa";
                    $short_fa_2 = $b . "_2.fa";
                    if(defined $sam){
                        $sam_file = $short_fa_prefix . ".sam";
                    }
                }
            }
        }
        else{
            $short_fa_1 = $short_fa_prefix . "_1.fa";
            $short_fa_2 = $short_fa_prefix . "_2.fa";
            if(defined $sam){
                $sam_file = $short_fa_prefix . ".sam";
            }
        }
    } # end of if of skip_il
    while(<g_h>){
        chomp;
        next if($_ =~ /^#/);
        if($_ =~ /^(\d+)\.(\d+)\.(\d+)/){
            # tmp index in trio . tmp index in sample . read index in tmp
            if($skip_pb == 0){
                my $long_fa_dir_ = $long_fa_dir;
                # PB read
                if($1 < $nums[0]){
                    # the first sample
                    $long_fa_dir_ =~ s/TOREPLACE/$strs[0]/;
                }
                elsif($1 < $nums[0] + $nums[1]){
                    # the second sample
                    $long_fa_dir_ =~ s/TOREPLACE/$strs[1]/;
                }
                else{
                    # the third sample
                    $long_fa_dir_ =~ s/TOREPLACE/$strs[2]/;
                }
                my $seq = `samtools faidx ${long_fa_dir_}$1/$long_fa_name $1.$2`;
                print OUT1 $seq;
            }
        }
        elsif($skip_il == 0){
            # short read
            if($_ =~ /^\d+$/){
                # indexed and in mapped
                chomp;
                my $line = ($_+1)*2;
                # extract only the sequence, leaving the readname off for not traversing the whole file (;q will stop the travering at the first hit)
                # fake read name
                print OUT2 ">$readnum\n";
                print OUT2 `sed \'$line!d;q\' $short_fa_1`;
                print OUT3 ">$readnum\n";
                print OUT3 `sed \'$line!d;q\' $short_fa_2`;
                if(defined $sam){
                    # output sam format as well
                    my $prev_line = $line - 1;
                    print OUT_sam `sed \'$prev_line!d;q\' $sam_file`;
                    print OUT_sam `sed \'$line!d;q\' $sam_file`;
                }
                $readnum ++;
            }
            elsif($_ =~ /^(\d+)\.un$/){
                # indexed and in both unmapped
                my $prefix = $h->{un};
                my $f1 = $prefix . "_1.fa";
                my $f2 = $prefix . "_2.fa";
                my $line = ($1+1)*2;
                print OUT2 ">$readnum\n";
                print OUT2 `seq \'$line!d;q\' $short_fa_1`;
                print OUT2 ">$readnum\n";
                print OUT3 `sed \'$line!d;q\' $short_fa_2`;
                $readnum ++;
            }
            else{
                # enter here only if not in the previous two, and previous two only if not here

                my $t;
                if($_ =~ /^(.+)\/\d+_\d+(\/[12])$/){
                    $t = $1 . $2;
                }
                else{
                    $t = $_;
                }
                my $seq;
                if($t =~ /^(.+)[\/_]1$/){
                    $short_fa = $short_fa_prefix . "_1.fa";
                    $seq = `samtools faidx $short_fa $1`;
                    print OUT2 $seq;
                }
                elsif($t =~ /^(.+)[\/_]2$/){
                    $short_fa = $short_fa_prefix . "_2.fa";
                    $seq = `samtools faidx $short_fa $1`;
                    print OUT3 $seq;
                }
                else{
                    $short_fa = $short_fa_prefix . "_1.fa";
                    $seq = `samtools faidx $short_fa $t`;
                    print OUT2 $seq;
                    $short_fa = $short_fa_prefix . "_2.fa";
                    $seq = `samtools faidx $short_fa $t`;
                    print OUT3 $seq;
                }
            }
        } # end of if of skip_il
    }
    close g_h;
    if(defined $tag && $tag eq "PB_only" || !defined $tag || defined $tag && $tag eq "NA"){
        close OUT1;
    }
    if(defined $tag && $tag eq "IL_only" || !defined $tag || defined $tag && $tag eq "NA"){
        close OUT2;
        close OUT3;
        close OUT4;
    }
    if(defined $sam){
        close OUT_sam;
    }
}
# two frg file for short and long
sub readname2fa_dir{
    my ($group, $short_fa_prefix, $long_fa_dir, $long_fa_name, $tag, $sam) = @_;
    # this takes the readnames in group, make a new fa with both short and long reads, long read name has been hashed in long_name
    my $name;
    my $skip_il = 0;
    my $skip_pb = 0;

    if(defined $tag && $tag eq "PB_only"){
        $skip_il = 1;
    }
    elsif(defined $tag && $tag eq "IL_only"){
        $skip_pb = 1;
    }
#    my ($g) = ($group =~ /(out.g\d+\.txt)/);
    # make the output file in the same directory as out.g*.txt
    my $g = $group;
    my $out1 = $g . ".PB.fa";
    my $out2 = $g . ".IL.1.fa";
    my $out3 = $g . ".IL.2.fa";
    my $out4 = $g . ".IL.fa";
    my $out5 = $g . ".sam";
    if(defined $sam){
        open OUT_sam, ">$out5" or die $!;
    }
    if(defined $tag && $tag eq "PB_only" || !defined $tag || defined $tag && $tag eq "NA"){
        open OUT1, ">$out1" or die $!;
    }
    if(defined $tag && $tag eq "IL_only" || !defined $tag || defined $tag && $tag eq "NA"){
        open OUT2, ">$out2" or die $!;
        open OUT3, ">$out3" or die $!;
        open OUT4, ">$out4" or die $!;
    }
    open g_h, "<$group" or die $!;
    my $short_fa;
    my $sam_file;
    my $h;
    # the hash to indicate which short read file to grep
    my $short_fa_1;
    my $short_fa_2;
    my $readnum = 0;
    if($skip_il == 0){
        if($short_fa_prefix =~ /;/){
            # two source files of short reads, in the format of short1;str:short2
            my @a = split(/;/, $short_fa_prefix);
            foreach my $b (@a){
                if($b =~ /:/){
                    my @c = split(/:/, $b);
                    $h->{$c[0]} = $c[1];
                }
                else{
                    $h->{"NA"} = $b;
                    $short_fa_1 = $b . "_1.fa";
                    $short_fa_2 = $b . "_2.fa";
                    if(defined $sam){
                        $sam_file = $short_fa_prefix . ".sam";
                    }
                }
            }
        }
        else{
            $short_fa_1 = $short_fa_prefix . "_1.fa";
            $short_fa_2 = $short_fa_prefix . "_2.fa";
            if(defined $sam){
                $sam_file = $short_fa_prefix . ".sam";
            }
        }
    } # end of if of skip_il
    my $h_pb;
    while(<g_h>){
        chomp;
        next if($_ =~ /^#/);
        if($_ =~ /^(\d+)\.(\d+)/){
            if($skip_pb == 0){
                # PB read
                # remove duplicates
                if(!defined $h_pb->{"$1.$2"}){ 
                    my $seq = `samtools faidx ${long_fa_dir}$1/$long_fa_name $1.$2`;
                    print OUT1 $seq;
                    $h_pb->{"$1.$2"} = 1;
                }
            }
        }
        elsif($skip_il == 0){
            # short read
            if($_ =~ /^\d+$/){
                # indexed and in mapped
                chomp;
                my $line = ($_+1)*2;
                # extract only the sequence, leaving the readname off for not traversing the whole file (;q will stop the travering at the first hit)
                # fake read name
                print OUT2 ">$readnum\n";
                print OUT2 `sed \'$line!d;q\' $short_fa_1`;
                print OUT3 ">$readnum\n";
                print OUT3 `sed \'$line!d;q\' $short_fa_2`;
                if(defined $sam){
                    # output sam format as well
                    my $prev_line = $line - 1;
                    print OUT_sam `sed \'$prev_line!d;q\' $sam_file`;
                    print OUT_sam `sed \'$line!d;q\' $sam_file`;
                }
                $readnum ++;
            }
            elsif($_ =~ /^\d+\_\d+$/){
                # now accomodate the Illumina reads with "_" as separator
                my $tmp = `samtools faidx $short_fa_1:$_`;
                print OUT2 $tmp;
                $tmp = `samtools faidx $short_fa_2:$_`;
                print OUT2 $tmp;
            }
            elsif($_ =~ /^(\d+)\.un$/){
                # indexed and in both unmapped
                my $prefix = $h->{un};
                my $f1 = $prefix . "_1.fa";
                my $f2 = $prefix . "_2.fa";
                my $line = ($1+1)*2;
                print OUT2 ">$readnum\n";
                print OUT2 `seq \'$line!d;q\' $short_fa_1`;
                print OUT2 ">$readnum\n";
                print OUT3 `sed \'$line!d;q\' $short_fa_2`;
                $readnum ++;
            }
            else{
                # enter here only if not in the previous two, and previous two only if not here

                my $t;
                if($_ =~ /^(.+)\/\d+_\d+(\/[12])$/){
                    $t = $1 . $2;
                }
                else{
                    $t = $_;
                }
                my $seq;
                if($t =~ /^(.+)[\/_]1$/){
                    $short_fa = $short_fa_prefix . "_1.fa";
                    $seq = `samtools faidx $short_fa $1`;
                    print OUT2 $seq;
                }
                elsif($t =~ /^(.+)[\/_]2$/){
                    $short_fa = $short_fa_prefix . "_2.fa";
                    $seq = `samtools faidx $short_fa $1`;
                    print OUT3 $seq;
                }
                else{
                    $short_fa = $short_fa_prefix . "_1.fa";
                    $seq = `samtools faidx $short_fa $t`;
                    print OUT2 $seq;
                    $short_fa = $short_fa_prefix . "_2.fa";
                    $seq = `samtools faidx $short_fa $t`;
                    print OUT3 $seq;
                }
            }
        } # end of if of skip_il
    }
    close g_h;
    if(defined $tag && $tag eq "PB_only" || !defined $tag || defined $tag && $tag eq "NA"){
        close OUT1;
    }
    if(defined $tag && $tag eq "IL_only" || !defined $tag || defined $tag && $tag eq "NA"){
        close OUT2;
        close OUT3;
        close OUT4;
    }
    if(defined $sam){
        close OUT_sam;
    }
}


sub readname2fa{
    my ($group, $short_fa, $long_fa) = @_;
    # this takes the readnames in group, make a new fa with both short and long reads, long read name has been hashed in long_name
    my $name;
    my $out = $group . ".fa";
    open OUT, ">$out" or die $!;
    open g_h, "<$group" or die $!;
    while(<g_h>){
        chomp;
        if($_ =~ /^\d+$/){
            # PB read
            my $seq = `samtools faidx $long_fa $_`;
            print OUT $seq;
        }
        else{
            # short read
            my ($t) = ($_ =~ /^(.+)\/\d+_\d+$/);
            my $seq = `samtools faidx $short_fa $t`;
            print OUT $seq;
        }
    }
    close g_h;
    close OUT;
}

# write a txt file, with each row the first column a il or pb read, the second column the pb or ils it relates to, with comma separating them
sub read_m4_record_relation{
    my ($m4, $out, $p_state, $read_state) = @_;
    my $hash;
    open M4, "<$m4" or die $!;
    while(<M4>){
        my @a = split(/\s+/, $_);
        $a[0] = &standardize_rn($a[0]);
        if(!defined $read_state || $read_state eq "IL_only"){
            if(!defined $hash->{IL}->{$a[0]}){
                $hash->{IL}->{$a[0]} = $a[1];
            }
            else{
                $hash->{IL}->{$a[0]} .= ",$a[1]";
            }
        }
        if(!defined $read_state || $read_state eq "PB_only"){
            if(!defined $hash->{PB}->{$a[1]}){
                $hash->{PB}->{$a[1]} = $a[0];
            }
            else{
                $hash->{PB}->{$a[1]} .= ",$a[0]";
            }
        }
    }
    close M4;
    if($p_state eq "print"){
        # write to a file
        open OUT, ">$out" or die $!;
        foreach my $k (sort keys %$hash){
            foreach my $k1 (sort keys %{$hash->{$k}}){
                print OUT join("\t", $k1, $hash->{$k}->{$k1}) . "\n";
            }
        }
        close OUT;
    }
    else{
        if(!defined $read_state){
            return $hash;
        }
        elsif($read_state eq "PB_only"){
            return $hash->{PB};
        }
        else{
            return $hash->{IL};
        }
    }
}

# start from an il, use breadth first search to get PBs; then start from PBs, use breadth first search, and get ils with >= m pbs, do it again with >= m ils to catch PBs. Input is a hash like txt file, each row corresponding to a read (il, or pb), with the second column the reads it connects to , separated by comma
sub sub_clustering{
    my ($txt_f, $m, $out_prefix) = @_;
    my $hash;
    my $max_il = 0;
    my $max_pb = 0;
    my $max_id = -1;
    open fh_, "<$txt_f" or die $!;
    while(<fh_>){
        chomp;
        my @a = split(/\t/, $_);
        my $key = "IL";
        $key = "PB" if($a[0] =~ /\./);
        my @b = split(/,/, $a[1]);
        foreach my $b_ (@b){
            $hash->{$key}->{$a[0]}->{$b_} = 1;
        }
    }
    close fh_;
    my $n = 0;
    while(%{$hash->{IL}}){
        my @keys = keys %{$hash->{IL}};
        # step 1, use this il to catch pbs
        my @pbs = keys %{$hash->{IL}->{$keys[0]}};
        # step 2, use each pb to catch ils
        my $count;
        foreach my $pb_ (@pbs){
            my @ils = keys %{$hash->{PB}->{$pb_}};
            foreach my $i (@ils){
                $count->{$i} ++ if(defined $count->{$i});
                $count->{$i} = 1 if(!defined $count->{$i});
            }
        }
        my $enroll;
        foreach my $i (keys %$count){
            if($count->{$i} >= $m){
                # enroll this il
                $enroll->{IL}->{$i} = 1;
            }
        }
        # step 3, use each enrolled il to catch pbs
        %$count = ();
        foreach my $i (keys %{$enroll->{IL}}){
            my @ps = keys %{$hash->{IL}->{$i}};
            foreach my $p (@ps){
                $count->{$p} ++ if(defined $count->{$p});
                $count->{$p} = 1 if(!defined $count->{$p});
            }
        }
        foreach my $p (keys %$count){
            if($count->{$p} >= $m){
                # enroll this pb
                $enroll->{PB}->{$p} = 1;
            }
        }
        # delete relations if enrolled
        foreach my $p (keys %{$enroll->{PB}}){
            # delete the second layer first
            foreach my $i (keys %{$hash->{PB}->{$p}}){
                delete $hash->{IL}->{$i}->{$p};
            }
            delete $hash->{PB}->{$p};
        }
        foreach my $i (keys %{$enroll->{IL}}){
            # delete the second layer first
            foreach my $p (keys %{$hash->{IL}->{$i}}){
                delete $hash->{PB}->{$p}->{$i};
            }
            delete $hash->{IL}->{$i};
        }
        # delete the first lead il no matter it was enrolled or not for convergence
        foreach my $p (keys %{$hash->{IL}->{$keys[0]}}){
            delete $hash->{PB}->{$p}->{$keys[0]};
        }
        delete $hash->{IL}->{$keys[0]};
        if(scalar(keys %{$enroll->{PB}}) != 0 && scalar(keys %{$enroll->{IL}}) != 0){
            # save this enroll to a file
            my $out = "$out_prefix.g$n.txt";
            $n ++;
            open out_fh, ">$out" or die $!;
            foreach my $keys (sort keys %$enroll){
                foreach my $ks (sort {$a <=> $b} keys %{$enroll->{$keys}}){
                    print out_fh $ks . "\n";
                }
            }
            close out_fh;
            if(scalar keys %{$enroll->{IL}} > $max_il){
                $max_il = scalar keys %{$enroll->{IL}};
                $max_pb = scalar keys %{$enroll->{PB}};
                $max_id = $n-1;
            }
        }
    }
    print "There are in all $n clusters in $out_prefix.g*.txt. The maximum cluster is in $out_prefix.g$max_id.txt, and there are $max_il, $max_pb il and pb in it.\n";
}


sub find_group{
    my ($read_id, $group_file) = @_;
    open g_fh, "<$group_file" or die $!;
    my $l = "";
    while(<g_fh>){
        if($_ =~ /^#/){
            $l = $_;
        }
        elsif($_ !~ /\d+\.\d+/){
            # illumina reads, compare
            chomp;
            if($_ eq $read_id){
                $l =~ s/^#//;
                my @ll = split(/\t/, $l);
                return @ll;
            }
        }
    }
    close g_fh;
    return ("NA", "NA", "NA");
}

sub find_status_group{
    my ($group, $ctg_file, $called_file) = @_;
    my $assembled_ = `grep g$group.ctg $ctg_file`;
    my $called_ = "NA";
    $called_ = `grep g$group.ctg $called_file`;
    my $assembled = 0;
    my $called = 0;
    if($assembled_ ne ""){
        $assembled = 1;
    }
    if($called_ ne ""){
        $called = 1;
    }
    return ($assembled, $called, $called_);
}
sub grep_clusters_reads_print_to_file{
    my ($output_txt, $cluster_file, $file) = @_;
    my $tag = 0;
    open c_fh, "<$cluster_file" or die $!;
    my $h;
    while(<c_fh>){
        my @a = split(/\t/, $_);
        $h->{$a[0]} = 1;
    }
    close c_fh;
    open out_fh, ">$file" or die $!;
    open out_txt, "<$output_txt" or die $!;
    while(<out_txt>){
        if($_ =~ /^#(\d+)/){
            if(defined $h->{$1}){
                print out_fh $_;
                $tag = 1;
            }
            else{
                $tag = 0;
            }
        }
        elsif($_ !~ /^#/ && $tag == 1){
            print out_fh $_;
        }
    }
    close out_txt;
    close out_fh;
}


sub grep_a_single_cluster_reads_print_to_file{
    my ($output_txt, $cluster_num, $file) = @_;
    my $tag = 0;
    open out_fh, ">$file" or die $!;
    open out_txt, "<$output_txt" or die $!;
    while(<out_txt>){
        if($_ =~ /^#(\d+)/ && $tag == 0){
            if($1 == $cluster_num){
                $tag = 1;
            }
        }
        elsif($_ !~ /^#/ && $tag == 1){
            print out_fh $_;
        }
        elsif($_ =~ /^#/ && $tag == 1){
            last;
        }
    }
    close out_txt;
    close out_fh;
}


sub grep_a_single_cluster_reads{
    my ($output_txt, $cluster_num) = @_;
    my $tag = 0;
    open out_txt, "<$output_txt" or die $!;
    while(<out_txt>){
        if($_ =~ /^#(\d+)/ && $tag == 0){
            if($1 == $cluster_num){
                $tag = 1;
            }
        }
        elsif($_ !~ /^#/ && $tag == 1){
            print $_;
        }
        elsif($_ =~ /^#/ && $tag == 1){
            last;
        }
    }
    close out_txt;
}

# given a graph, find the cut vertex without which the graph will be broken into two connected components
# only report those nodes with node number > node_t
sub get_cut_vertex{
    my ($g, $source, $node_t, $cut_vertices_ori) = @_;
    my $debug = 0;
    # get dfs of g => gg
    my ($gg, $node_cor, $node_cor_rev) = graph::dfs_m($g, $source);
    # find missing edges (back traversal edges)
    my $edges = graph::find_missing_edge($g, $gg, $node_cor);
    # key part
    # from small to large nodes, traversal from back edges, tag as visited, stop at visited, find all cycles
    # the cut vertex would be all first vertices of all cycles except the first one
    # k_node is the hash recording which has been visited
    # k is the new graph recording which edges have been traversed
    my $k;
    my $k_node;
    my $tag = 0;
    my @cut_vertices;
    foreach my $node (sort {$a <=> $b} keys %$edges){
        foreach my $n (sort {$a <=> $b} keys %{$edges->{$node}}){
            $k_node->{$node} = 1;
            $k->{$node}->{$n} = 1;
            while(!defined $k_node->{$n}){
                $k_node->{$n} = 1;
                my $x = (keys %{$gg->{$n}})[0];
                $k->{$n}->{$x} = 1;
                $n = $x;
            }
            # difference from finding bridges here
            # check if a cycle
            if($n == $node){
                # a cycle
                if($tag == 1){
                    # not the first cycle
                    if(scalar(keys %{$g->{$node_cor_rev->{$n}}}) > $node_t){
                        push @cut_vertices, $n;
                    }
                }
                else{
                    $tag = 1;
                }
            }
        }
    }
    # convert to original coordinate
    foreach my $node (@cut_vertices){
        my $node_ori = $node_cor_rev->{$node};
        $cut_vertices_ori->{$node_ori} = 1;
    }
    #print join("\n", @cut_vertices_ori) . "\n";
    return $cut_vertices_ori;
}

# run union_find on a tree and return the tree of the largest in tree with the source
sub get_cc_graph{
    my ($g, $size_t, $rest_g) = @_;
    # the graph with trees > size_t
    my $gg;
    my $tr;
    my $sz;
    # the key is the source, then followed by each connected component
    my $gg_cc;
    foreach my $e (keys %$g){
        $tr->{$e} = $e;
        $sz->{$e} = 1;
        foreach my $e1 (keys %{$g->{$e}}){
            $tr->{$e1} = $e1;
            $sz->{$e1} = 1;
        }
    }
    foreach my $e (keys %$g){
        foreach my $e1 (keys %{$g->{$e}}){
            my $e_r = &root_hash($e, $tr);
            my $e1_r = &root_hash($e1, $tr);
            if($e_r ne $e1_r){
                ($tr, $sz) = &connect_hash($e_r, $e1_r, $tr, $sz);
            }
        }
    }
    my ($h, $h_s) = &tree_cluster_ori_hash($tr);
    my $hh;
    # get the trees with size > size_t, save node to hh
    foreach my $k (keys %$h){
        if($h_s->{$k} > $size_t){
            $hh->{$k} = 1;
            foreach my $kk (keys %{$h->{$k}}){
                $hh->{$kk} = 1;
            }
        }
    }
    # source on connected component > size_t only
    # make a new graph only with the nodes in the tree > size_t
    foreach my $e (keys %$g){
        foreach my $e1 (keys %{$g->{$e}}){
            if(defined $hh->{$e} || defined $hh->{$e1}){
                $gg->{$e}->{$e1} = 1;
                $gg->{$e1}->{$e} = 1;
            }
            else{
                $rest_g->{$e}->{$e1} = 1;
                $rest_g->{$e1}->{$e} = 1;
            }
        }
    }
    foreach my $k (keys %$h){
        # all the sources
        if(defined $hh->{$k}){
            # satisfy the size requirement
            foreach my $kk (keys %{$h->{$k}}){
                # all the nodes in this connected component
                $gg_cc->{$k}->{$kk} = $gg->{$kk};
            }
        }
    }


    return ($gg, $h_s, $rest_g, $gg_cc);
}

# given a graph, find the bridge without which the graph will be broken into two connected components
sub get_bridge{
    my ($g, $source) = @_;
    my $debug = 0;
    # g is a hash {node1}->{node2} = 1; (node1 < node2)
    my ($gg, $node_cor, $node_cor_rev) = graph::dfs_m($g, $source);
    if($debug == 1){
        print "Conversion:\n";
        foreach my $t (keys %$node_cor){
            print join("\t", $t, $node_cor->{$t}) . "\n";
        }
        print "Tree:\n";
        graph::print_tree($gg);
        print "End of tree:\n";
    }
    # gg is a tree, each edge pointing towards the root direction, with the nodes renumbered in the order traversed in dfs, node is the correspondence from nodes in g to nodes in gg
    my $edges = graph::find_missing_edge($g, $gg, $node_cor);
    if($debug == 1){
        print "Missing edge 1:\n";
        foreach my $x (keys %$edges){
            foreach my $y (keys %{$edges->{$x}}){
                print join("\t", $x, $y) . "\n";
            }
            print "End of missing edge 1\n";
        }
    }
    # find the missing edge in g that is not found in gg, in gg coordinate, pointing backward from root
    # traverse each edge
    # k_node is the hash recording which has been visited
    # k is the new graph recording which edges have been traversed
    my $k;
    my $k_node;
    foreach my $node (sort {$a <=> $b} keys %$edges){
        foreach my $n (sort {$a <=> $b} keys %{$edges->{$node}}){
            $k_node->{$node} = 1;
            $k->{$node}->{$n} = 1;
            while(!defined $k_node->{$n}){
                $k_node->{$n} = 1;
                my $x = (keys %{$gg->{$n}})[0];
                $k->{$n}->{$x} = 1;
                $n = $x;
            }
        }
    }
    if($debug == 1){
        print "New tree:\n";
        graph::print_tree($k);
        print "End of new tree with missing edges\n";
        print "Bridges:\n";
    }
    my $bridges = graph::find_missing_edge($g, $k, $node_cor);
    # hash to be returned as bridge
    my $bd;
    foreach my $node (keys %$bridges){
        my $node_ori = $node_cor_rev->{$node};
        foreach my $nb (keys %{$bridges->{$node}}){
            my $nb_ori = $node_cor_rev->{$nb};
            if($nb_ori =~ /\./){
                $bd->{$node_ori}->{$nb_ori} = 1;
                #print join("\t", $node_ori, $nb_ori) . "\n";
            }
            else{
                $bd->{$nb_ori}->{$node_ori} = 1;
                #print join("\t", $nb_ori, $node_ori) . "\n";
            }
        }
    }
    return $bd;
}

sub make_graph_from_m4{
    my ($m4) = @_;
    my $g;
    my $s;
    open m4_fh, "<$m4" or die $!;
    while(<m4_fh>){
        my @a = split(/\s+/, $_);
        my @b = split(/\//, $a[0]);
        $g->{$b[0]}->{$a[1]} = 1;
        $g->{$a[1]}->{$b[0]} = 1;
        $s = $b[0];
    }
    close m4_fh;
    my $source;
    $source->{$s} = 1;
    return ($g, $source);
}




1;
