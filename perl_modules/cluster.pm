use warnings;
use strict;
package cluster;

# given a sam from IL to ctg, cluster the reads and return the big clusters in a hash
sub cluster_read{
    my ($sam, $size_t) = @_;
    my $h;
    open fh_, "<$sam" or die $!;
    while(<fh_>){
        next if($_ =~ /^@/);
        my @aa = split(/\s+/, $_);
        $h->{$aa[3]} = 1;
    }
    close fh_;
    my $dist_t = 60;
    my $h_clu;
    my $n = 1;
    my $end;
    my $start;
    my $size;
    my $pre = -1;
    foreach my $key (sort {$a <=> $b} keys %$h){
        if($pre != -1 && $key - $pre < $size_t){
            $size ++;
            $end = $key;
            $pre = $key;
        }
        elsif($pre == -1){
            $pre = $key;
            $size = 1;
            $start = $key;
            $end = $key;
        }
        elsif($key - $pre >= $size_t){
            if($size > $size_t){
                $h_clu->{$n}->{size} = $size;
                $h_clu->{$n}->{start} = $start;
                $h_clu->{$n}->{end} = $end;
                $n ++;
            }
            $pre = $key;
            # start over
            $start = $key;
            $end = $key;
            $size = 1;
        }
    }
    return $h_clu;
}


# This looks into the start and end from all ils on a PB. If the length > a threshold (no both unmapped), then either the PB is repetitive, or the ils have errors.
# delete this PB and corresponding alignments
sub PB_range{
    my ($m4, $length_t, $out) = @_;
    my ($h_pb) = &read_pb_range($m4);
    my $h_pb_d;
    foreach my $pb (keys %$h_pb){
        my $range = $h_pb->{$pb}->{max} - $h_pb->{$pb}->{min};
        if( $range > $length_t ){
            $h_pb_d->{$pb} = 1;
        }
    }
    &write_m4_by_pb($m4, $h_pb_d, $out);
}

# this reads a m4, check if the pb is in the hash, print if not
sub write_m4_by_pb{
    my ($m4, $h, $out) = @_;
    open m4, "<$m4" or die $!;
    open m4_out, ">$out" or die $!;
    while(<m4>){
        my @a = split(/\s+/, $_);
        if(!defined $h->{$a[1]}){
            print m4_out $_;
        }
    }
    close m4;
    close m4_out;
}

sub read_pb_range{
    my ($m4) = @_;
    my $h;
    open m4, "<$m4" or die $!;
    while(<m4>){
        my @a = split(/\s+/, $_);
        my $pb = $a[1];
        my $s = $a[9];
        if($a[8] == 1){
            $s = $a[11] - $s;
        }
        if(!defined $h->{$pb}->{min} || $s < $h->{$pb}->{min}){
            $h->{$pb}->{min} = $s;
        }
        if(!defined $h->{$pb}->{max} || $s > $h->{$pb}->{max}){
            $h->{$pb}->{max} = $s;
        }
    }
    return $h;
}
        



# This throw a ball from an il, let it bounce back and forth for a maximum times, or stop when it hits the start again, record the ids it hits, and let those that it hit frequently be a cluster, iteratively do it on the rest of the alignments
sub bounce{
    my ($m4, $max_bounce, $num_try, $min_hit, $out_prefix) = @_;
    my ($h_il, $h_pb) = &read_cluster_pairs($m4);
    my $n = 0;
    while(%$h_il){
        # there are still ils in the hash
        # get the first key (il)
        my @k = keys %$h_il;
        $DB::single = 1;
        my $il = $k[0];
        my ($i, $p) = &walk($max_bounce, $num_try, $min_hit, $il, $h_il, $h_pb);
        print "Finish walking\n";
        if(scalar(@$i) != 0){
            print "Writing to $n\n";
            &print_cluster($i, $p, $out_prefix, $n);
            $n ++;
            # remove the keys
            foreach my $i_ (@$i){
                # remove second layer from pb first
                foreach my $pb (keys %{$h_il->{$i_}}){
                    delete $h_pb->{$pb}->{$i_};
                    delete $h_pb->{$pb} if(scalar(keys %{$h_pb->{$pb}}) == 0);
                }
                delete $h_il->{$i_} if(defined $h_il->{$i_});
                print "Deleting il $i_\n";
            }
            foreach my $p_ (@$p){
                # remove second layer from il first
                foreach my $il (keys %{$h_pb->{$p_}}){
                    delete $h_il->{$il}->{$p_};
                    delete $h_il->{$il} if(scalar(keys %{$h_il->{$il}}) == 0);
                }
                delete $h_pb->{$p_} if(defined $h_pb->{$p_});
                print "Deleting pb $p_\n";
            }
        }
    }
}

# print a cluster assisting bounce
sub print_cluster{
    my ($i, $p, $out_prefix, $n) = @_;
    my $f = "$out_prefix.g$n.txt";
    open fh_, ">$f" or die $!;
    foreach my $i_ (@$i){
        print fh_ $i_ . "\n";
    }
    foreach my $p_ (@$p){
        print fh_ $p_ . "\n";
    }
    close fh_;
}
        
# take the random number
sub random{
    my ($num) = @_;
    return int(rand($num));
}

# This is a funciton to assist bounce
sub walk{
    my ($max_b, $num_t, $min_h, $s, $h_il, $h_pb) = @_;
    # for record to return
    my @i;
    my @p;
    my $c = $s;
    my $p_h;
    my $i_h;
    foreach my $t (1 .. $num_t){
        my $count_p;
        my $count_i;
        my $it = 0;
        # count the times it hits start again
        my $n = 0;
        while($it < $max_b){
            # forward
            my @keys = keys %{$h_il->{$c}};
            print "IL key: $c\t" . scalar(@keys) . "\n";
            my $random_num = &random(scalar(@keys));
            my $cp = $keys[$random_num];
            $DB::single = 1;
            $count_p->{$cp} ++ if(defined $count_p->{$cp});
            $count_p->{$cp} = 1 if(!defined $count_p->{$cp});
            # backward
            @keys = keys %{$h_pb->{$cp}};
            print "PB key: $cp\t" . scalar(@keys) . "\n";
            $random_num = &random(scalar(@keys));
            $c = $keys[$random_num];
            $count_i->{$c} ++ if(defined $count_i->{$c});
            $count_i->{$c} = 1 if(!defined $count_i->{$c});
            $it ++;
            # check if hit start
            if($c eq $s){
                $n ++;
            }
            if($c eq $s && $n == 3){
                print "Hit start\n";
                last;
            }
        }
        if($it < $max_b){
            # previous run hit start
            # record their counts
            foreach my $key (keys %$count_p){
                # every time should be with hit > min_h
                if($count_p->{$key} > $min_h){
                    $p_h->{$key} = 1;
                }
            }
            foreach my $key (keys %$count_i){
                if($count_i->{$key} > $min_h){
                    $i_h->{$key} = 1;
                }
            }
        }
        else{
            print "Run starting from $c does not end within $max_b times\n";
        }
    }
    # now return the recorded ids
    @i = keys %$i_h;
    @p = keys %$p_h;
    return (\@i, \@p);
}

sub read_cluster_pairs{
    my ($c) = @_;
    my $h_il;
    my $h_pb;
    open c_h, "<$c" or die $!;
    while(<c_h>){
        chomp;
        my @a = split(/\s+/, $_);
        my ($il) = &standardize_rn($a[0]);
        my $pb = $a[1];
        $h_il->{$il}->{$pb} = 1;
        $h_pb->{$pb}->{$il} = 1;
    }
    return ($h_il, $h_pb);
}

sub read_cluster{
    my ($c) = @_;
    my $h_il;
    my $h_pb;
    open c_h, "<$c" or die $!;
    while(<c_h>){
        chomp;
        if($_ =~ /\d+\.\d+/){
            $h_pb->{$_} = 1;
        }
        elsif($_ =~ /^\d+$/){
            $h_il->{$_} = 1;
        }
    }
    close c_h;
    return ($h_il, $h_pb);
}

# simply record each line into a hash (each line is a index of a il read that has more than 80 alignments to pb)
sub read_f{
    my ($f) = @_;
    my $h;
    open f_, "<$f" or die $!;
    while(<f_>){
        chomp;
        $h->{$_} = 1;
    }
    close f_;
    return $h;
}


# read the reads in a cluster, output to a m4 all the related alignments given a m4 file
sub find_alignments{
    my ($c, $m4, $out) = @_;
    my ($h_il, $h_pb) = &read_cluster($c);
    open m4_fh, "<$m4" or die $!;
    open OUT, ">$out" or die $!;
    while(<m4_fh>){
        my @a = split(/\s+/, $_);
        # standardize readname for IL, no interval, no 1 or 2
        my ($il) = &standardize_rn($a[0]);
        if(defined $h_il->{$il} && defined $h_pb->{$a[1]}){
            print OUT $_;
        }
    }
    close m4_fh;
    close OUT;
}

# parallelization, since each line is independent
sub parallel_split{
    my ($f, $out, $line) = @_;
    open out_fh, ">$out" or die $!;
    open f_fh, "<$f" or die $!;
    my $n = 0;
    while(<f_fh>){
        $n ++;
        if($n % $line == 0){
            my $p = tell(f_fh);
            print out_fh $p . "\n";
        }
    }
    close f_fh;
    close out_fh;
}
# an easier way when given a file with a list of repeated il reads
sub split_cluster_quick{
    my ($f, $f_r, $c, $sub_file, $out, $max_pb) = @_;
    my ($h_il, $h_pb) = &read_cluster($c);
    my $h = &read_f($f_r);
    foreach my $k (keys %$h){
        if(defined $h_il->{$k}){
            delete $h_il->{$k};
        }
    }

    open fh_, "<$f" or die $!;
    open OUT_fh, ">$sub_file" or die $!;
    while(<fh_>){
        my @a = split(/\s+/, $_);
        # standardize readname for IL, no interval, no 1 or 2
        my ($il) = &standardize_rn($a[0]);

        if(defined $h_il->{$il}){
            my $pb = $a[1];
            print OUT_fh join("\t", $il, $pb)."\n";
        }
    }
    close fh_;
    close OUT_fh;
}

# remove the alignments with illumina reads appearing in f2, output a new m4
sub remove_aln{
    my ($f2, $m4, $out) = @_;
    my ($h_il, $h_pb) = &read_cluster($f2);
    open m4_fh, "<$m4" or die $!;
    open OUT, ">$out" or die $!;
    while(<m4_fh>){
        my @a = split(/\s+/, $_);
        # standardize readname for IL, no interval, no 1 or 2
        my ($il) = &standardize_rn($a[0]);
        if(!defined $h_il->{$il}){
            print OUT $_;
        }
    }
    close m4_fh;
    close OUT;
}

    
# this is to split large cluster with repeat regions
sub split_cluster{
    my ($f, $c, $sub_file, $out, $max_pb) = @_;
    # f: original m4; c: cluster file containing the read names; sub_file: a temporary file to store the rest of the relations of pb and il; out: cluster after split
    my ($h_il, $h_pb) = &read_cluster($c);
    my $h_pb_this;
    my $h_il_this;
    open fh_, "<$f" or die $!;
    while(<fh_>){
        my @a = split(/\s+/, $_);
        # standardize readname for IL, no interval, no 1 or 2
        my ($il) = &standardize_rn($a[0]);
        my $pb = $a[1];

        if(defined $h_il->{$il}){
            # record the pb it aligns to
            $h_il_this->{$il} .= ",$pb" if(defined $h_il_this->{$il});
            $h_il_this->{$il} = $pb if(!defined $h_il_this->{$il});
        }
        
        if(defined $h_pb->{$pb}){
            # record the il name and il start position
            my $ori = $a[8];
            my $s = $a[9];
            my $size = $a[11];
            $s = $size - $s if($ori == 1);
            $h_pb_this->{$pb} .= ",$il:$s" if(defined $h_pb_this->{$pb});
            $h_pb_this->{$pb} = "$il:$s" if(!defined $h_pb_this->{$pb});
        }
    }
    close fh_;

    open OUT, ">$sub_file" or die $!;
    my $tag_il;
    # look at il
    foreach my $il (keys %$h_il_this){
        my @a = split(/,/, $h_il_this->{$il});
        if(scalar(@a) > $max_pb){
            $tag_il->{$il} = 1;
            $DB::single = 1;
            print join("\t", $il, scalar(@a))."\n";
        }
    }
    # look at pb
    foreach my $pb (keys %$h_pb_this){
        my @a = split(/,/, $h_pb_this->{$pb});
        my $flag = 0;
        foreach my $a_ (@a){
            my @b = split(/:/, $a_);
            if(!defined $tag_il->{$b[0]}){
                $DB::single = 1;
                print OUT join("\t", $b[0], $pb) . "\n";
            }
        }
    }

    close OUT;
    &union_find($sub_file, $out);


}


sub stats{
    # count the nubmer of pbs this il aligns to, and vice versa, print to a file
    my ($f) = @_;
    my $h_p;
    my $h_i;
    open fh_, "<$f" or die $!;
    while(<fh_>){
        my @a = split(/\s+/, $_);
        ($a[0]) = &standardize_rn($a[0]);
        $h_i->{$a[0]} ++;
        $h_p->{$a[1]} ++;
    }
    close fh_;
    open OUT1, ">stat.i.txt";
    open OUT2, ">stat.p.txt";
    foreach my $i (sort {$a <=> $b} keys %$h_i){
        print OUT1 $i . "\t" . $h_i->{$i} . "\n";
    }
    foreach my $i (sort {$a <=> $b} keys %$h_p){
        print OUT2 $i . "\t" . $h_p->{$i} . "\n";
    }
    close OUT1;
    close OUT2;
    print "There are in all " . scalar(keys %$h_i) . " ils and " . scalar(keys %$h_p) . " pbs.\n";
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
    &print_stats($h_s, $out);
    # print names
    &print_names($h, $out);
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

sub print_names {
    my ($h, $out) = @_;
    my $num = 0;
    foreach my $i (sort {$a <=> $b } keys %$h){
        my $str = $h->{$i};
        my $new_f = $out . ".g$num.txt";
        $num ++;
        open fh_, ">$new_f" or die $!;
        print fh_ $str;
        close fh_;
    }
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

sub root{
    my ($p, $tr) = @_;
    while($tr->[$p] != $p){
        #$tr->[$p] = $tr->[$tr->[$p]];
        $p = $tr->[$p];
    }
    return $p;
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

# now take into account a pos corresponding to multiple reads
sub smart_cluster{
    # This cluster looks at a pb, foreach il, count the ils that share the same breakpoint, for those above a threshold, say 3, output the relation. For sharing teh same breakpoint, the ils should overlap with each other, or considering novel micro-insertions, a gap smaller than that threshold. Output is a file with each row the il and pb whose relation is preserved, in two columns.
    # 10202015. now changed so that insertion can be taken in. THere is no size_t limit, but the cluster_s is the cluster length divided by the number of reads. THus allow multiple breakpoints.
    # size_t is according to the largest range a true SV can have, which is novel insertion. Since both unmapped are not here, the largest distance between two ils in one cluster should be < 3*insert size, plus 3*standard_deviation. Give it 1500 seems reasonable.
    my $version = "#V_0.0 (10202015) Allows novel insertions. Cluster size max limit is the average per base coverage. \n";
    print $version;
    my ($m4, $micro_ins_t, $cluster_s, $out) = @_;
    # the hash to record the positions of ils on pb (orientation normalized)
    my $h;
    my $start = "NA";
    my $end = "NA";
    my $pre_qname = "NA";
    $DB::single = 1;
    open m4, "<$m4" or die $!;
    while(<m4>){
        my ($rn_q, $rn_r, $as, $a_per, $ori_q, $s_q, $e_q, $len_q, $ori_r, $s_r, $e_r, $len_r, ) = split(/\s+/, $_); 
        my $qname = &standardize_rn($rn_q);
        # calculate the other end
        if($ori_r == 1){
            # the end to be calculated
            $end = $len_r - ($s_r - $s_q);
        }
        elsif($ori_r == 0){
            # the start to be calculated
            $start = $s_r - $s_q;
        }

        if($qname eq $pre_qname && ($start ne "NA" && $end ne "NA")){
            # s_r is the start on pb for this il
            # record it for this pb, with keys the start, values the end, note it does not count for two reads starting at the same position, ignore it for now
            push @{$h->{$rn_r}->{$start}}, $end;
            $start = "NA";
            $end = "NA";
        }
        $pre_qname = $qname;
    }
    close m4;
    # the range of the ils start that should be included for this pb relation, inclusive
    my $h_r;
    foreach my $pb (keys %$h){
        $h_r->{$pb} = &range_bp_ovl($h->{$pb}, $micro_ins_t, $cluster_s);
    }
    # now read the m4 again and examine if the start is within any range
    # a tag to see the previous read printed or not
    my $pre_print = 0;
    $pre_qname = "NA";
    open OUT, ">$out" or die $!;
    open m4, "<$m4" or die $!;
    while(<m4>){
        my $l1 = $_;
        my ($rn_q, $rn_r, $as, $a_per, $ori_q, $s_q, $e_q, $len_q, $ori_r, $s_r, $e_r, $len_r, ) = split(/\s+/, $l1); 
        my $l2 = <m4>;
        if($ori_r == 1){
            ($rn_q, $rn_r, $as, $a_per, $ori_q, $s_q, $e_q, $len_q, $ori_r, $s_r, $e_r, $len_r, ) = split(/\s+/, $l2); 
        }
        $s_r -= $s_q;
        my $c_id = -1;
        foreach my $range (@{$h_r->{$rn_r}}){
            $c_id ++;
            my ($s, $e) = split(/:/, $range);
            if($s_r >= $s && $s_r <= $e){
                # a little modification to distinguish the PB that supports this cluster and that cluster, add a suffix to PB name
                chomp $l1;
                chomp $l2;
                my @a = split(/\s+/, $l1);
                $l1 = join(" ", $a[0], $a[1] . ".$c_id",  @a[2 .. $#a]). "\n";
                @a = split(/\s+/, $l2);
                $l2 = join(" ", $a[0], $a[1] . ".$c_id", @a[2 .. $#a]) . "\n";
                print OUT $l1, $l2;
                last;
            }
        }
    }
    close m4;
    close OUT;
}

# this function looks at the il pos on pb, find the most possible clusters by share the same bp (or with a gap < micro insertion threshold), see if the cluster size >= cluster_s, if yes, output the range. It looks at each il to see the biggest extension possible. Return is an array reference of the ranges, each range with start-end format. Input is a hash with the keys the starts, the values the end
# 10202015: for INS: modify so that 1). If a cluster (continuous part with overlapping short reads) has < cluster_s_min short reads, skip it. 2). If a cluster has average short read (# short read pair * 200 / total length < $cluster_s_max), suppose read length is 100. Can set cluster_s_max at a relatively high number to allow some fluctuation. 
sub range_bp_ovl{
    my ($h_p, $micro_in_s, $cluster_s) = @_;
    my ($cluster_s_min, $cluster_s_max) = ($cluster_s, 50);
    ($cluster_s_min, $cluster_s_max) = split(/:/, $cluster_s) if($cluster_s =~ /:/);
    # a hash recording the number of reads following up that overlap with this reads, or gap with micro insertion size; and the last read start.
    my $h_c;
    # the p of the last il that overlap with the key of this hash
    my $h_l;
    foreach my $s (sort {$a <=> $b} keys %$h_p){
        foreach my $pre_p (sort {$a <=> $b} keys %$h_l){
            foreach my $end (@{$h_p->{$pre_p}}){
                if($end > $s - $micro_in_s){
                    # there is overlap
                    $h_l->{$pre_p} = $s;
                    $h_c->{$pre_p} ++;
                }
            }
        }
        $h_c->{$s} = 1;
        $h_l->{$s} = $s;
    }
    # select the largest one, second largest, until there is no cluster with size >= cluster_s
    my @keys_by_value = sort {$h_c->{$b} <=> $h_c->{$a}} keys %$h_c;
    my @values = sort {$b <=> $a} values %$h_c;
    my @ret;
    my @dumped;
    $DB::single = 1;
    foreach my $i (0 .. scalar(@values) - 1){
        my $count = $values[$i];
        my $p = $keys_by_value[$i];
        my $l = $h_l->{$p};
        # check min, check max by looking at average sequence coverage
        if($count >= $cluster_s_min && $count / ($l - $p) * 200 <= $cluster_s_max){
            # check if this read has been included and reported before
            my $tag = 1;
            foreach my $r (@ret){
                my ($s, $e) = split(/:/, $r);
                if($s <= $p && $e >= $p || $s <= $l && $e >= $l){
                    $tag = 0;
                }
            }
            foreach my $r (@dumped){
                my ($s, $e) = split(/:/, $r);
                if($s <= $p && $e >= $p || $s <= $l && $e >= $l){
                    $tag = 0;
                }
            }
            if($tag == 1){
                # add this range
                # note: not checking the length of the cluster any more (10202015)
                # next if($l - $p > $size_t);
                # thus allowing overlap of ranges, but not repetitively report the same bp because of some of the ils in the range
                push @ret, "$p:$l";
            }
        }
        else{
            # record the dumped ones in case the subsets will be recorded
            push @dumped, "$p:$l";
        }
    }
    return \@ret;
}
# now take into account a pos corresponding to multiple reads
sub smart_cluster_singleBP{
    # This cluster looks at a pb, foreach il, count the ils that share the same breakpoint, for those above a threshold, say 3, output the relation. For sharing teh same breakpoint, the ils should overlap with each other, or considering novel micro-insertions, a gap smaller than that threshold. Output is a file with each row the il and pb whose relation is preserved, in two columns.
    # size_t is according to the largest range a true SV can have, which is novel insertion. Since both unmapped are not here, the largest distance between two ils in one cluster should be < 3*insert size, plus 3*standard_deviation. Give it 1500 seems reasonable.
    my ($m4, $micro_ins_t, $cluster_s, $out, $size_t) = @_;
    # the hash to record the positions of ils on pb (orientation normalized)
    my $h;
    my $start = "NA";
    my $end = "NA";
    my $pre_qname = "NA";
    $DB::single = 1;
    open m4, "<$m4" or die $!;
    while(<m4>){
        my ($rn_q, $rn_r, $as, $a_per, $ori_q, $s_q, $e_q, $len_q, $ori_r, $s_r, $e_r, $len_r, ) = split(/\s+/, $_); 
        my $qname = &standardize_rn($rn_q);
        # calculate the other end
        if($ori_r == 1){
            # the end to be calculated
            $end = $len_r - ($s_r - $s_q);
        }
        elsif($ori_r == 0){
            # the start to be calculated
            $start = $s_r - $s_q;
        }

        if($qname eq $pre_qname && ($start ne "NA" && $end ne "NA")){
            # s_r is the start on pb for this il
            # record it for this pb, with keys the start, values the end, note it does not count for two reads starting at the same position, ignore it for now
            push @{$h->{$rn_r}->{$start}}, $end;
            $start = "NA";
            $end = "NA";
        }
        $pre_qname = $qname;
    }
    close m4;
    # the range of the ils start that should be included for this pb relation, inclusive
    my $h_r;
    foreach my $pb (keys %$h){
        $h_r->{$pb} = &range_bp_ovl_singleBP($h->{$pb}, $micro_ins_t, $cluster_s, $size_t);
    }
    # now read the m4 again and examine if the start is within any range
    # a tag to see the previous read printed or not
    my $pre_print = 0;
    $pre_qname = "NA";
    open OUT, ">$out" or die $!;
    open m4, "<$m4" or die $!;
    while(<m4>){
        my $l1 = $_;
        my ($rn_q, $rn_r, $as, $a_per, $ori_q, $s_q, $e_q, $len_q, $ori_r, $s_r, $e_r, $len_r, ) = split(/\s+/, $l1); 
        my $l2 = <m4>;
        if($ori_r == 1){
            ($rn_q, $rn_r, $as, $a_per, $ori_q, $s_q, $e_q, $len_q, $ori_r, $s_r, $e_r, $len_r, ) = split(/\s+/, $l2); 
        }
        $s_r -= $s_q;
        my $c_id = -1;
        foreach my $range (@{$h_r->{$rn_r}}){
            $c_id ++;
            my ($s, $e) = split(/:/, $range);
            if($s_r >= $s && $s_r <= $e){
                # a little modification to distinguish the PB that supports this cluster and that cluster, add a suffix to PB name
                chomp $l1;
                chomp $l2;
                my @a = split(/\s+/, $l1);
                $l1 = join(" ", $a[0], $a[1] . ".$c_id",  @a[2 .. $#a]). "\n";
                @a = split(/\s+/, $l2);
                $l2 = join(" ", $a[0], $a[1] . ".$c_id", @a[2 .. $#a]) . "\n";
                print OUT $l1, $l2;
                last;
            }
        }
    }
    close m4;
    close OUT;
}

# this function looks at the il pos on pb, find the most possible clusters by share the same bp (or with a gap < micro insertion threshold), see if the cluster size >= cluster_s, if yes, output the range. It looks at each il to see the biggest extension possible. Return is an array reference of the ranges, each range with start-end format. Input is a hash with the keys the starts, the values the end
sub range_bp_ovl_singleBP{
    my ($h_p, $micro_in_s, $cluster_s, $size_t) = @_;
    my ($cluster_s_min, $cluster_s_max) = ($cluster_s, 50);
    ($cluster_s_min, $cluster_s_max) = split(/:/, $cluster_s) if($cluster_s =~ /:/);
    # a hash recording the number of reads following up that overlap with this reads, or gap with micro insertion size; and the last read start.
    my $h_c;
    # the p of the last il that overlap with the key of this hash
    my $h_l;
    foreach my $s (sort {$a <=> $b} keys %$h_p){
        foreach my $pre_p (sort {$a <=> $b} keys %$h_l){
            foreach my $end (@{$h_p->{$pre_p}}){
                if($end > $s - $micro_in_s){
                    # there is overlap
                    $h_l->{$pre_p} = $s;
                    $h_c->{$pre_p} ++;
                }
            }
        }
        $h_c->{$s} = 1;
        $h_l->{$s} = $s;
    }
    # select the largest one, second largest, until there is no cluster with size >= cluster_s
    my @keys_by_value = sort {$h_c->{$b} <=> $h_c->{$a}} keys %$h_c;
    my @values = sort {$b <=> $a} values %$h_c;
    my @ret;
    my @dumped;
    $DB::single = 1;
    foreach my $i (0 .. scalar(@values) - 1){
        my $count = $values[$i];
        my $p = $keys_by_value[$i];
        my $l = $h_l->{$p};
        if($count >= $cluster_s_min && $count <= $cluster_s_max){
            # check if this read has been included and reported before
            my $tag = 1;
            foreach my $r (@ret){
                my ($s, $e) = split(/:/, $r);
                if($s <= $p && $e >= $p || $s <= $l && $e >= $l){
                    $tag = 0;
                }
            }
            foreach my $r (@dumped){
                my ($s, $e) = split(/:/, $r);
                if($s <= $p && $e >= $p || $s <= $l && $e >= $l){
                    $tag = 0;
                }
            }
            if($tag == 1){
                # add this range
                next if($l - $p > $size_t);
                # thus allowing overlap of ranges, but not repetitively report the same bp because of some of the ils in the range
                push @ret, "$p:$l";
            }
        }
        else{
            # record the dumped ones in case the subsets will be recorded
            push @dumped, "$p:$l";
        }
    }
    return \@ret;
}


# In a big tree, start from the largest trunk (pb), pull all the ils out, pull al the pbs that have at least n ils out, iteratively until nothing can be pulled out. Make it a group. Start from the largest trunk (pb) in the rest of the tree. Iterate until no tree can be pulled out. Tree is a group file with il and pb ids.
sub pull_out_trunk{
    my ($m4, $tree, $n, $out_prefix) = @_;
    # record the number of reads in a big file
    open OUT_all, ">$out_prefix.all" or die $!;
    # record the relation
    my ($h_il, $h_pb) = &read_cluster($tree);
    # read m4, and record the relations
    my ($h_i, $h_p);
    open m4, "<$m4" or die $!;
    while(<m4>){
        my @a = split(/\s+/, $_);
        my $il_nm = &standardize_rn($a[0]);
        if(defined $h_il->{$il_nm} && defined $h_pb->{$a[1]}){
            $h_i->{$il_nm}->{$a[1]} = 1;
            $h_p->{$a[1]}->{$il_nm} = 1;
        }
    }
    close m4;

    my $id = 0;
    # for print for each group;
    my $print_pbs;
    my $print_ils;
    while(scalar(keys %$h_p) != 0){
        $DB::single = 1;

        # print the previous one
        if(scalar(keys %$print_pbs) != 0 && scalar(keys %$print_ils) != 0){
            my $out = $out_prefix . ".g$id.txt";
            $id ++;
            open OUT, ">$out" or die $!;
            print OUT_all join("\t", $out, scalar(keys %$print_ils), scalar(keys %$print_pbs)). "\n";
            foreach my $i (sort {$a <=> $b} keys %$print_ils){
                print OUT $i . "\n";
                foreach my $p (keys %{$h_i->{$i}}){
                    delete $h_p->{$p}->{$i};
                    delete $h_p->{$p} if(scalar(keys %{$h_p->{$p}}) == 0);
                }
                delete $h_i->{$i};
            }
            foreach my $p (sort keys %$print_pbs){
                print OUT $p . "\n";
                foreach my $i (keys %{$h_p->{$p}}){
                    delete $h_i->{$i}->{$p};
                    delete $h_i->{$i} if(scalar(keys %{$h_i->{$i}}) == 0);
                }
                delete $h_p->{$p};
            }
            close OUT;
        }
        $print_pbs = {};
        $print_ils = {};
        # start from the largest trunk (pb)
        my $tr_sz = 0;
        my $select_pb;
        foreach my $pb (keys %$h_p){
            my $tmp = scalar(keys %{$h_p->{$pb}});
            if($tmp > $tr_sz){
                $tr_sz = $tmp;
                $select_pb = $pb;
            }
        }

        if($tr_sz >= $n){
            # pull
            my $pull = 1;
            my $pull_pbs;
            my $pull_ils;
            $print_pbs->{$select_pb} = 1;
            foreach my $tmp (keys %{$h_p->{$select_pb}}){
                $pull_ils->{$tmp} = 1;
                $print_ils->{$tmp} = 1;
            }
            # begin pulling
            while($pull == 1){
                # count the pull horse power for each pulled pb
                $pull = 0;
                my $hash_pb;
                my $hash_il;
                foreach my $il (keys %$pull_ils){
                    foreach my $pulled_pb (keys %{$h_i->{$il}}){
                        # only count those haven't been added to this cluster yet (those added to previous clusters should have been deleted before the start of this round)
                        if(!defined $print_pbs->{$pulled_pb}){
                            if(!defined $hash_pb->{$pulled_pb}){
                                $hash_pb->{$pulled_pb} = 1;
                            }
                            else{
                                $hash_pb->{$pulled_pb} ++;
                            }
                        }
                    }
                }
                $pull_ils = {};
                $pull_pbs = {};
                foreach my $pulled_pb (keys %$hash_pb){
                    if($hash_pb->{$pulled_pb} >= $n){
                        # this pb has been successfully pulled out
                        $pull_pbs->{$pulled_pb} = 1;
                        $print_pbs->{$pulled_pb} = 1;
                        $pull = 1;
                    }
                }
                if($pull == 1){
                    $pull = 0;
                    # use pb to pull il
                    foreach my $pulled_pb (keys %$pull_pbs){
                        foreach my $pulled_il (keys %{$h_p->{$pulled_pb}}){
                            if(!defined $print_ils->{$pulled_il}){
                                if(!defined $hash_il->{$pulled_il}){
                                    $hash_il->{$pulled_il} = 1;
                                }
                                else{
                                    $hash_il->{$pulled_il} ++;
                                }
                            }
                        }
                    }
                    foreach my $pulled_il (keys %$hash_il){
                        if($hash_il->{$pulled_il} >= $n){
                            $pull_ils->{$pulled_il} = 1;
                            $print_ils->{$pulled_il} = 1;
                            $pull = 1;
                        }
                    }
                }
            } # end of while
        } # end of if
        else{
            last;
        }
    } # end of while
    # print the previous one
    $DB::single = 1;
    if(scalar(keys %$print_pbs) != 0 && scalar(keys %$print_ils) != 0){
        my $out = $out_prefix . ".g$id.txt";
        open OUT, ">$out" or die $!;
        print OUT_all join("\t", $out, scalar(keys %$print_ils), scalar(keys %$print_pbs)). "\n";
        foreach my $i (keys %$print_ils){
            print OUT $i . "\n";
            foreach my $p (keys %{$h_i->{$i}}){
                delete $h_p->{$p}->{$i};
                delete $h_p->{$p} if(scalar(keys %{$h_p->{$p}}) == 0);
            }
            delete $h_i->{$i};

        }
        foreach my $p (keys %$print_pbs){
            print OUT $p . "\n";
            foreach my $i (keys %{$h_p->{$p}}){
                delete $h_i->{$i}->{$p};
                delete $h_i->{$i} if(scalar(keys %{$h_i->{$i}}) == 0);
            }
            delete $h_p->{$p};

        }
        close OUT;
    }
    close OUT_all;
}






1;


