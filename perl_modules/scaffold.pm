use warnings;
use strict;
package scaffold;

my $blasr_binary_new = "~/blasrn";
my $bwa = "~/pkg/bwa/bwa";
# check if cigar is good enough: if sum of INDELs > this number, it is counted as not good
my $cigar_threshold = 5;
# flanking size of the segment of the contig
my $segment_flank = 500;
# gap for discerning if two breakpoints belong to one on the contig (in analyze sam)
my $gap = 20;
# maximum gap between two clusters of ILs on contig before they can be merged
my $micro_in_s = 1500;
# in my merge, the minimum overlaplength to tolerate for mergin
my $overlap_length = 300;
# in my merge, the maximum tail gap to tolerate for it counted as a tail concatenation
my $tail_gap = 100;

# check if bwa mem is good enough for local alignment
sub examine_mem{
    my ($aln) = @_;
    foreach my $aln (@$aln){
        next if($aln =~ /^@/);
        my @a = split(/\t/, $aln);
        my $check_cigar_result = &check_cigar_result($a[5], $cigar_threshold);
        # when there is no good alignment, return 0
        if($check_cigar_result == 1){
            return 1;
        }
    }
    return 0;
}

# given a cigar, see if it is good enough
# simply check the sum of D and I
sub check_cigar_result{
    my ($cigar, $threshold) = @_;
    my $num;
    my $tag;
    my $e = 0;
    while($cigar ne ""){
        ($num, $tag, $cigar) = ($cigar =~ /^(\d+)([MISHD])(.*)$/);
        if($tag =~ /[ID]/){
            $e += $num;
            if($e > $threshold){
                return 0;
            }
        }
    }
    if($e > $threshold){
        return 0;
    }
    return 1;
}

# align segments in segment fa file to the reference
# apply INV's alignment process to other types
sub align_seg2ctg{
    my ($segment_file, $ref, $type) = @_;
    my $flank = 500;
    my $h;
    my $seg_align = "$segment_file.sam";
    #if($type eq "INV"){
    # use bwa mem to align the segments
    my @alns = split(/\n/, `$bwa mem $ref $segment_file`);
    open seg_sam_fh, ">$seg_align" or die $!;
    print seg_sam_fh join("\n", @alns) . "\n";
    close seg_sam_fh;
    # see bwa is good enough for detecting INV, look at cigar string
    my $examine_mem_result = &examine_mem(\@alns);
    if($examine_mem_result == 1){
        return $seg_align;
    }

    foreach my $aln (@alns){
        next if($aln =~ /^@/);
        my @a = split(/\t/, $aln);
        $h->{$a[0]}->{"$a[2]:$a[3]:$a[5]"} = 1;
    }
    foreach my $rd_name (keys %$h){
        # remade the test fa file for itself
        #`samtools faidx $segment_file $rd_name > $rd_name.fa`;
        `rm $rd_name.ref.fa`;
        foreach my $key (keys %{$h->{$rd_name}}){
            my ($chr, $p, $cigar) = split(/:/, $key);
            my $num;
            my $tag;
            my $length = 0;
            while($cigar ne ""){
                ($num, $tag, $cigar) = ($cigar =~ /^(\d+)([MDIHS])(.*)$/);
                if($tag =~ /[DM]/){
                    $length += $num;
                }
            }
            my $x = $p - $flank;
            my $y = $p + $length + $flank;
            # extract reference
            `samtools faidx $ref $chr:$x-$y >> $rd_name.ref.fa`;
        }
        # realign the rd with the reconstructed reference with nucmer
        `~/pkg/mummer/MUMmer3.23/nucmer $rd_name.ref.fa $segment_file`; 
        `~/pkg/mummer/MUMmer3.23/delta-filter out.delta > out.delta.filtered`;
        &convert_delta2sam("out.delta.filtered", "out.sam");
        # note: seg_align now should only contain three alignments, and no order disturb is allowed
        `cat out.sam >> $seg_align`;
    }

    #}
    return $seg_align;
}

sub convert_delta2sam{
    my ($delta, $sam) = @_;
    # delta file in the format of 
    # argument
    # command
    # >ref_region read_region ref_length read_length
    # ref_start ref_end read_start read_end
    # a number (0)
    # second alignment of this read corresponding to this reference
    # third ...
    open fh_, "<$delta" or die $!;
    open sam_fh, ">$sam" or die $!;
    my ($ref_name, $read_name, $ref_length, $read_length);
    # read is a subsequence extracted from contig
    my ($ctg_name, $read_start, $read_end, $ctg_length);
    my $tag = 0;
    my $tmp;
    while(<fh_>){
        if($tag != 1){
            next if($_ !~ /^>/);
        }
        my @a = split(/\s+/, $_);
        if($_ =~ /^>/){
            $tag = 1;
            ($ref_name, $read_name, $ref_length, $read_length) = split(/\s+/, $_);
            ($tmp, $ctg_name, $read_start, $read_end, $ctg_length) = split(/\./, $read_name);
            $ref_name =~ s/^>//;
        }
        elsif(scalar(@a) > 4){
            # alignment row
            my $r = -1;
            my ($ref_s, $ref_e, $read_s, $read_e, ) = split(/\s+/, $_);
            if($read_s < $read_e){
                $r = 1;
            }
            my $flag;
            my $head_clip = 0;
            my $tail_clip = 0;
            if($r == -1){
                # reverse
                $flag = 16;
                # note that read_s and read_e are in the reverse order (read_s > $read_e) in delta_filter
                $head_clip = ($read_length - $read_s) + ($ctg_length - $read_length - $read_start) - 1;
                $tail_clip = $read_start + $read_e - 1; 
            }
            else{
                # forward
                $flag = 0;
                $head_clip = $read_s + $read_start - 1;
                $tail_clip = $read_length - $read_e + ($ctg_length - $read_length - $read_start) - 1;
            }
            my ($chr, $s, $e) = ($ref_name =~ /^(.+):(\d+)-(\d+)$/);
            # all are M
            print sam_fh join("\t", $read_name, $flag, $chr, ($ref_s + $s), 99, $head_clip . "S" . ($ref_e - $ref_s + 1) . "M" . $tail_clip . "S") . "\n";
        }
    } # end of while loop
    close sam_fh;
    close fh_;
} # end of convert_delta2sam function

# for strategy, see ~/combinePBIL/scripts/pp/pipe/scaffold_ctgs.pl
sub analyze_align{
    my ($seg_align, $out) = @_;
    my $h;
    my $h_cigar;
    open OUT, ">>$out" or die $!;
    open fh_, "<$seg_align" or die $!;
    while(<fh_>){
        next if($_ =~ /^@/);
        my @a = split(/\t/, $_);
        my @b = split(/\//, $a[0]);
        my @c = split(/\./, $b[0]);
        #my ($ctg_s, $ctg_e) = ($c[$#c - 1], $c[$#c]);
        my $ctg_s = $c[$#c-2];
        my $ctg_e = $c[$#c-1];
        my $chr = $a[2];
        my $pos1 = $a[3];
        my $cigar = $a[5];
        my $length = 0;
        my $num;
        my $tag;
        while($cigar ne ""){
            ($num, $tag, $cigar) = ($cigar =~ /^(\d+)([SHMID])(.*)$/);
            $length += $num if($tag =~ /[MD]/);
        }
        $cigar = $a[5];
        my $pos2 = $pos1 + $length;
        # get orientation (-1: forward; 1: backforward);
        my $r = -1;
        my $r_s = $pos1;
        my $r_e = $pos2;
        if($a[1] & 0x0010){
            $r = 1;
            # do not need to reverse the two since they are already in the increasing order
            #$r_s = $pos2;
            #$r_e = $pos1;
        }
        # get the actual alignment region
        if($cigar =~ /^(\d+)[SH]/){
            if($r == -1){
                $ctg_s += $1;
            }
            else{
                $ctg_e -= $1;
            }
        }
        if($cigar =~ /(\d+)[SH]$/){
            if($r == -1){
                $ctg_e -= $1;
            }
            else{
                $ctg_s += $1;
            }
        }
        # note that reverse need to be reversed again, head -> tail, tail -> head
        #if($r == 1){
        #    my $tmp = $ctg_e;
        #    $ctg_e = $ctg_s;
        #    $ctg_s = $tmp;
        #}
        # only save the clipped end
        if($cigar =~ /^\d+[SH]/){
            if($r == -1){
                # 5' end of contig was clipped, 5' end of reference was clipped
                $h->{$b[0]}->{$ctg_s} = "$chr:$pos1:5";
            }
            else{
                # 3' end of contig was clipped, 5' end of reference was clipped
                $h->{$b[0]}->{-$ctg_e} = "$chr:$pos1:5";
            }
        }
        if($cigar =~ /\d+[SH]$/){
            if($r == -1){
                # 3' end of contig was clipped, 3' end of reference was clipped
                $h->{$b[0]}->{-$ctg_e} = "$chr:$pos2:3";
            }
            else{
                # 5' end of contig was clipped, 3' end of reference was clipped
                $h->{$b[0]}->{$ctg_s} = "$chr:$pos2:3";
            }
        }
    } # end of reading the alignment

    # begin analyzing the alignment
    foreach my $ctg (keys %$h){
        # look at any two contig breakpoints are close to each other (for translocation)
        my @p = sort{$a <=> $b} keys %{$h->{$ctg}};

        for(my $i = 0; $i < scalar(@p) - 1; $i ++){
            for(my $j = $i + 1; $j < scalar(@p); $j ++){
                    $DB::single = 1;
                if(abs(abs($p[$j]) - abs($p[$i])) < $gap && $p[$j] * $p[$i] < 0){
                    # a translocation is confirmed
                    # TODO make it more readable
                    # now in the format of (readname, bp1_on_ctg (-1 means 3'end of ctg clipped), chr1:pos1:[35] (3 means 3' end of reference clipped), bp2_on-ctg, chr2:pos2:[35])
                    print OUT join("\t", $ctg, $p[$i], $h->{$ctg}->{$p[$i]}, $p[$j], $h->{$ctg}->{$p[$j]}) . "\n";
                }
            }
        }
    }
    close OUT;
}





sub analyze_align_INV{
    my ($seg_align, $out) = @_;
    my $h;
    my $h_cigar;
    open OUT, ">>$out" or die $!;
    open fh_, "<$seg_align" or die $!;
    while(<fh_>){
        next if($_ =~ /^@/);
        my @a = split(/\t/, $_);
        my @b = split(/\//, $a[0]);
        my @c = split(/\./, $b[0]);
        #my ($ctg_s, $ctg_e) = ($c[$#c - 1], $c[$#c]);
        my $ctg_s = 0;
        my $ctg_e = $c[$#c];
        my $chr = $a[2];
        my $pos1 = $a[3];
        my $ctg_name = join(".", @c[0 .. $#c - 2]);
        my $cigar = $a[5];
        my $length = 0;
        my $num;
        my $tag;
        while($cigar ne ""){
            ($num, $tag, $cigar) = ($cigar =~ /^(\d+)([SHMID])(.*)$/);
            $length += $num if($tag =~ /[MD]/);
        }
        $cigar = $a[5];
        my $pos2 = $pos1 + $length;
        # get orientation (-1: forward; 1: backforward);
        my $r = -1;
        my $r_s = $pos1;
        my $r_e = $pos2;
        if($a[1] & 0x0010){
            $r = 1;
            # do not need to reverse the two since they are already in the increasing order
            #$r_s = $pos2;
            #$r_e = $pos1;
        }
        # get the actual alignment region
        if($cigar =~ /^(\d+)[SH]/){
            if($r == -1){
                $ctg_s += $1;
            }
            else{
                $ctg_e -= $1;
            }
        }
        if($cigar =~ /(\d+)[SH]$/){
            if($r == -1){
                $ctg_e -= $1;
            }
            else{
                $ctg_s += $1;
            }
        }
        # note that reverse need to be reversed again, head -> tail, tail -> head
        #if($r == 1){
        #    my $tmp = $ctg_e;
        #    $ctg_e = $ctg_s;
        #    $ctg_s = $tmp;
        #}
        $h->{$c[0]}->{"$ctg_s.$ctg_e"} = "$chr:$r_s-$r_e:$r";
        $h_cigar->{$c[0]}->{"$ctg_s.$ctg_e"} = $a[5];
    } # end of while fh_
    close fh_;
    # for each segment, look at any alignments corresponding to INV model
    my $min_gap = 50;;
    foreach my $rd (keys %$h){
        my $n = 0;
        my $pre_chr;
        my $pre_r_s;
        my $pre_r_e;
        my $pre_r;
        my $s;
        my $e;
        my $c;
        my $success = 0;
        # simply assume three alignments, with orientaiton +, -, and +
        foreach my $ctg_region (sort {$a <=> $b} keys %{$h->{$rd}}){
            my ($chr, $r_s, $r_e, $r) = ($h->{$rd}->{$ctg_region} =~ /^(.+):(\d+)-(\d+):(.*\d+)$/);
            if($n == 0){
                $pre_chr = $chr;
                $pre_r_e = $r_e;
                $pre_r = $r;
            }
            elsif($n == 1 || $n == 2){
                if($r * $pre_r != -1 || abs($pre_r_e - $r_s) > $min_gap || $pre_chr ne $chr){
                    $success = 0;
                    last;
                }
                else{
                    $success = 1;
                    $pre_chr = $chr;
                    $pre_r_e = $r_e;
                    $pre_r = $r;
                }
                if($n == 1){
                    $s = $r_s;
                    $e = $r_e;
                    $c = $chr;
                }
            }
            $n ++;
        }
        if($success == 1){
            print OUT join("\t", $c, $s, $e) . "\n";
        }
    }
    close OUT;
}



# look at the segments alignment and infer SV
# 1. Resolve duplicated or overlapping alignments without any conflict
# 2. Build the DAG, find the most complete path (with scores)
# 3. Infer SV from the most complete path. (similar to gustaf, or shaffle-Lagan, but gustaf does not deal with long contig, and they do not deal with small events between 10-30bp and large novel insertion; shuffle-LAGAN does not deal with invertion.)
# Goal is to deal with any type of SVs with size > 10bp. 
# segments
sub analyze_align_deprecated{
    my ($seg_align) = @_;
    my $h;
    my $h_cigar;
    open fh_, "<$seg_align" or die $!;
    while(<fh_>){
        next if($_ =~ /^@/);
        my @a = split(/\t/, $_);
        my @b = split(/\//, $a[0]);
        my @c = split(/\./, $b[0]);
        my ($ctg_s, $ctg_e) = ($c[$#c - 1], $c[$#c]);
        my $chr = $a[2];
        my $pos1 = $a[3];
        my $ctg_name = join(".", @c[0 .. $#c - 2]);
        my $cigar = $a[5];
        my $length = 0;
        while($cigar ne ""){
            my ($num, $tag, $cigar) = ($cigar =~ /^(\d+)([SHMID])(.*)$/);
            $length += $num if($tag =~ /[MD]/);
        }
        my $pos2 = $pos1 + $length;
        # get orientation
        my $r = 0;
        my $r_s = $pos1;
        my $r_e = $pos2;
        if($a[1] & 0x0010){
            $r = 1;
            $r_s = $pos2;
            $r_e = $pos1;
        }
        # get the actual alignment region
        if($cigar =~ /^(\d+)[SH]/){
            if($r == 0){
                $ctg_s += $1;
            }
            else{
                $ctg_e -= $1;
            }
        }
        if($cigar =~ /(\d+)[SH]$/){
            if($r == 0){
                $ctg_e -= $1;
            }
            else{
                $ctg_s += $1;
            }
        }
        $h->{"$ctg_s.$ctg_e"} = "$chr:$r_s-$r_e";
        $h_cigar->{"$ctg_s.$ctg_e"} = $a[5];
    } # end of while fh_

    # resolve duplicated alignments or partial alignments that do not have any conflict, merge them
    $h = &merge_alignments($h);

    # build DAG
} # end of analyze_align

# only to resolve overlapping alignments that do not have any conflict
# connect the two alignments if no conflict
sub merge_alignments{
    my ($h) = @_;
    # h with the keys contig_s:contig_e, the values chr:s-e, s to the point on contig (if reverse, contig_e, if forward, contig_s)

    my $tag = 1;
    while($tag != 0){
        # until there is no more to merge
        $tag = 0;
        my $h_replace;
        foreach my $k (sort {$a <=> $b} keys %$h){
            my ($s, $e) = split(/\./, $k);
            foreach my $k1 (sort {$a <=> $b} keys %$h){
                next if($k1 <= $k);
                my ($s1, $e1) = split(/\./, $k1);
                if($s < $e1 && $s >= $s1 || $s1 < $e && $s1 > $s){
                    # k1 and k overlap, look at their alignments
                    my $merged = &aln_compatible($h, $k, $k1);
                    if($merged ne "NA"){
                        $tag = 1;
                        $h_replace->{$k} = $merged;
                        $h_replace->{$k1} = $merged;
                        last;
                    }
                }
            }
        } # end of look at each key of h
        # replace with h_replace
        if($tag == 1){
            foreach my $k (keys %$h){
                if(defined $h_replace->{$k}){
                    my @x = split(/;/, $h_replace->{$k});
                    delete $h->{$k};
                    if(!defined $h->{$x[0]}){
                        $h->{$x[0]} = $x[1];
                    }
                }
            }
        }
    } # end of while tag
    return $h;
} # end of merge alignment sub

sub aln_compatible{
    my ($h, $k, $k1) = @_;
    my $max_dist = 100;
    my ($ctg_s, $ctg_e) = split(/\./, $k);
    my ($ctg_s1, $ctg_e1) = split(/\./, $k1);
    my $p = $h->{$k};
    my $p1 = $h->{$k1};
    my ($chr, $s, $e) = ($p =~ /^(.+):(\d+)-(\d+)$/);
    my ($chr1, $s1, $e1) = ($p1 =~ /^(.+):(\d+)-(\d+)$/);
    if($chr eq $chr1){
        if($s > $e && $s1 > $e1){
            # both are reverse
            if($s < $s1 && $s > $e1){
                if($e < $e1){
                    # |-------|
                    # e       s
                    #     |--------|
                    #     e1       s1
                    if(abs(($ctg_s - $ctg_e1) - ($s - $e1)) < $max_dist){
                        # merge
                        return "$ctg_s1.$ctg_e;$chr:$s1-$e";
                    }
                }
                else{
                    #       |--|
                    #       e  s
                    #     |--------|
                    #     e1       s1
                    if(abs(($ctg_s1 - $ctg_s) - ($s1 - $s)) < $max_dist && abs(($ctg_e - $ctg_e1) - ($e - $e1)) < $max_dist){
                        # merge
                        return "$ctg_s1.$ctg_e1;$chr:$s1-$e1";
                    }
                }
            } # end of if s in between s1 and e1
            elsif($s1 < $s && $s1 > $e){
                if($e1 < $e){
                    if(abs(($ctg_s1 - $ctg_e) - ($s1 - $e)) < $max_dist){
                        return "$ctg_s.$ctg_e1;$chr:$s-$e1";
                    }
                }
                else{
                    if(abs(($ctg_s - $ctg_s1) - ($s - $s1)) < $max_dist && abs(($ctg_e1 - $ctg_e) - ($e1 - $e)) < $max_dist){
                        return "$ctg_s.$ctg_e;$chr:$s-$e";
                    }
                }
            }
        } # end of if($s > $e && $s1 > $e1)
        elsif($s < $e && $s1 < $e1){
            # both are forward
            if($e < $e1 && $e > $s1){
                if($s < $s1){
                    # |-------|
                    # s       e
                    #     |--------|
                    #     s1       e1
                    if(abs(($ctg_s1 - $ctg_e) - ($s1 - $e)) < $max_dist){
                        # merge
                        return "$ctg_s.$ctg_e1;$chr:$s-$e1";
                    }
                }
                else{
                    #       |--|
                    #       s  e
                    #     |--------|
                    #     s1       e1
                    if(abs(($ctg_s1 - $ctg_s) - ($s1 - $s)) < $max_dist && abs(($ctg_e - $ctg_e1) - ($e - $e1)) < $max_dist){
                        # merge
                        return "$ctg_s1.$ctg_e1;$chr:$s1-$e1";
                    }
                }
            } # end of if s in between s1 and e1
            elsif($e1 < $e && $e1 > $s){
                if($s1 < $s){
                    if(abs(($ctg_s - $ctg_e1) - ($s - $e1)) < $max_dist){
                        return "$ctg_s1.$ctg_e;$chr:$s1-$e";
                    }
                }
                else{
                    if(abs(($ctg_s - $ctg_s1) - ($s - $s1)) < $max_dist && abs(($ctg_e1 - $ctg_e) - ($e1 - $e)) < $max_dist){
                        return "$ctg_s.$ctg_e;$chr:$s-$e";
                    }
                }
            }
        } # end of if($s < $e && $s1 < $e1)
    } #end of if($chr eq $chr1)
    return "NA";
}

# check if only one contig remain in the merged fa file
# TODO deal with 2 or >2 cases properly. Now they are not dealt with at all.
sub if_merged{
    my ($IL_align) = @_;
    my @names;
    foreach my $x (`grep "^>" $IL_align`){
        chomp $x;
        $x =~ s/^>//;
        push @names, $x;
    }
    return \@names;
}


sub segment{
    my ($ctg, $clusters, $key, $type) = @_;
    $DB::single = 1;
    $type = "NA" if(!defined $type || $type eq "");
    my $flank = 500;
    my $segments;
    my $rd_name = "NA";
    my $seq = "";
    my $length = 0;
    my $tag = 0;
    # suppose there is one contig in the ctg file
    open fh_, "<$ctg" or die $!;
    while(<fh_>){
        chomp;
        if($_ =~ /^>(\S+)$/){
            $rd_name = $1;
            if($1 eq $key){
                $tag = 1;
            }
            else{
                if($tag == 1){
                    last;
                }
            }
        }
        elsif($tag == 1){
            $seq .= $_;
        }
    }
    close fh_;
    $length = length($seq);
    my $coords;
    if($type eq "INV"){
        $coords = &get_segments_INV($clusters, $length, $flank);
    }
    else{
        $coords = &get_segments($clusters, $length, $flank);
    }
    my $n = 0;
    my @seqs;
    my $out_file = "tmp_seg.fa";
    open OUT_fh, ">$out_file";
    foreach my $coord (@$coords){
        my ($s, $e) = split(/\./, $coord);
        my $seq_ = substr($seq, $s, ($e - $s + 1));
        my $len = length($seq_);
        # add length of the extracted sequence to the readname
        print OUT_fh ">$rd_name.$s.$e.$len\n$seq_\n";
        push @seqs, $seq_;
    }
    close OUT_fh;
    return $out_file;
    #return \@seqs;
}

sub get_segments_INV{
    my ($clusters, $length, $flank) = @_;
    my $n = 0;
    my $x;
    my $y;
    my @coords;
    foreach my $c (sort {$a <=> $b} @$clusters){
        my ($s, $e) = split(/\./, $c);
        if($s - $flank < 0){
            $x = 0;
        }
        else{
            $x = $s - $flank;
        }
        if($e + $flank > $length - 1){
            $y = $length - 1;
        }
        else{
            $y = $e + $flank;
        }
        # get the condition of 3 and 7 
        push @coords, "$x.$y";

    }
    return \@coords;
}



# couple the neighboring two clusters together, and get the the pairs' coordinates
sub get_segments{
    my ($clusters, $length, $flank) = @_;
    my $n = 0;
    my $x;
    my $y;
    my @coords;
    my @sorted_clusters = sort {$a <=> $b} @$clusters;
    # in case there is only one cluster
    if(scalar(@sorted_clusters) == 1){
        my ($s, $e) = split(/\./, $sorted_clusters[0]);
        if($s - $flank < 0){
            $x = 0;
        }
        else{
            $x = $s - $flank;
        }
        if($e + $flank > $length - 1){
            $y = $length - 1;
        }
        else{
            $y = $e + $flank;
        }

        push @coords, "$x.$y";
    }
    else{
        # 2 or > 2 clusters
        # cut from the left of the first cluster to the right of the second cluster
        for(my $i = 0; $i < scalar(@sorted_clusters) - 1; $i ++){
            my $c = $sorted_clusters[$i];
            my $c1 = $sorted_clusters[$i + 1];
            my ($s, ) = split(/\./, $c);
            my ($tmp, $e) = split(/\./, $c1);


            if($s - $flank < 0){
                $x = 0;
            }
            else{
                $x = $s - $flank;
            }
            if($e + $flank > $length - 1){
                $y = $length - 1;
            }
            else{
                $y = $e + $flank;
            }
            push @coords, "$x.$y";
        }
    } # end of else
    return \@coords;
}







# now compatible with multiple contigs 
sub overlap{
    my ($cs) = @_;
    my $i;
    foreach my $rd_nm (keys %$cs){
        my $cs_a = $cs->{$rd_nm};
        for($i = 0; $i < scalar(@$cs_a) - 1; $i ++){
            my $c = $cs_a->[$i];
            my ($s, $e) = split(/\./, $c);
            my $j;
            for($j = $i + 1; $j < scalar(@$cs_a); $j ++){
                my $c1 = $cs_a->[$j];
                my ($s1, $e1) = split(/\./, $c1);
                if($s1 > $s && $s1 < $e || $s > $s1 && $s < $e1){
                    return 1;
                }
            }
        }
    }
    return 0;
}

# instead of corresponding to one breakpoint, now cluster reads as long as there are any overlap, ignore pair information for now
# for a breakpoint, all others are compared with one read; for continuous, it might be the case that others do not overlap with this read, but overlaped the second read that overlapped with this read
sub cluster_continuous{
    my ($h) = @_;
    # h is in the format of h->{chr}->{start} = @end
    my $start = 0;
    my $end = 0;
    my $ret;
    foreach my $chr (keys %$h){
        foreach my $p (sort {$a <=> $b} keys %{$h->{$chr}}){
            if($start == 0){
                $start = $p;
                foreach my $e (@{$h->{$chr}->{$p}}){
                    if($end < $e){
                        $end = $e;
                    }
                }
            }
            else{
                # already have start and end
                if($end < $p - $micro_in_s){
                    # new cluster
                    push @{$ret->{$chr}}, "$start.$end";
                    $start = $p;
                    $end = 0;
                    foreach my $e (@{$h->{$chr}->{$p}}){
                        if($end < $e){
                            $end = $e;
                        }
                    }
                }
                else{
                    foreach my $e (@{$h->{$chr}->{$p}}){
                        if($end < $e){
                            $end = $e;
                        }
                    }
                }
            }
        }
        push @{$ret->{$chr}}, "$start.$end";
    }
    return $ret;
}



# now can deal with multiple contigs
sub cluster{
    my ($f) = @_;
    my $h_p;
    open fh_, "<$f" or die $!;
    while(<fh_>){
        next if($_ =~ /^@/);
        my @a = split(/\s+/, $_);
        my $cigar = $a[5];
        my ($num, $tag);
        my $aln_len = 0;
        while($cigar ne ""){
            ($num, $tag, $cigar) = ($cigar =~ /^(\d+)([MISHD])(.*)$/);
            if($tag =~ /[MD]/){
                $aln_len += $num;
            }
        }
        my $end = $a[3] + $aln_len;
        push @{$h_p->{$a[2]}->{$a[3]}}, $end;
    }
    close fh_;
    my $ret = &cluster_continuous($h_p);
    # now cluster all the reads together, no mater whether they may correspond to a brekapoint
=cut
    # the p of the last il that overlap with the key of this hash
    # Xian changed two points here. 09262016
    # 1. $h_c ++ not redundant with multiple ends correspondign to the same start
    # 2. a pair of start and end corresonding to the read pair, instead of reads. Note: there should not be discordant (large insert size) read pair, since its within the genome. 
    my $ret;
    my $dumped;
    my $cluster_s_min = 3;
    foreach my $rd_nm (keys %$h_p){
        my $h_l;
        my $h_c;
        foreach my $s (sort {$a <=> $b} keys %{$h_p->{$rd_nm}}){
            foreach my $pre_p (sort {$a <=> $b} keys %$h_l){
                foreach my $end (@{$h_p->{$rd_nm}->{$pre_p}}){
                    if($end > $s - $micro_in_s){
                        # there is overlap
                        $h_l->{$pre_p} = $s;
                        $h_c->{$pre_p} += scalar(@{$h_p->{$rd_nm}->{$s}});
                        last;
                    }
                }
            }
            $h_c->{$s} = 1;
            $h_l->{$s} = $s;
        }
        # select the largest one, second largest, until there is no cluster with size >= cluster_s
        my @keys_by_value = sort {$h_c->{$b} <=> $h_c->{$a}} keys %$h_c;
        my @values = sort {$b <=> $a} values %$h_c;
        foreach my $i (0 .. scalar(@values) - 1){
            my $count = $values[$i];
            my $p = $keys_by_value[$i];
            my $l = $h_l->{$p};
            # check min, check max by looking at average sequence coverage
            if($count >= $cluster_s_min){
                # check if this read has been included and reported before
                my $tag = 1;
                foreach my $r (@{$ret->{$rd_nm}}){
                    my ($s, $e) = split(/\./, $r);
                    if($s <= $p && $e >= $p || $s <= $l && $e >= $l){
                        $tag = 0;
                    }
                }
                foreach my $r (@{$dumped->{$rd_nm}}){
                    my ($s, $e) = split(/./, $r);
                    if($s <= $p && $e >= $p || $s <= $l && $e >= $l){
                        $tag = 0;
                    }
                }
                if($tag == 1){
                    # add this range
                    # note: not checking the length of the cluster any more (10202015)
                    # next if($l - $p > $size_t);
                    # thus allowing overlap of ranges, but not repetitively report the same bp because of some of the ils in the range
                    push @{$ret->{$rd_nm}}, "$p.$l";
                }
            }
            else{
                # record the dumped ones in case the subsets will be recorded
                push @{$dumped->{$rd_nm}}, "$p.$l";
            }
        }
    }
=cut
return $ret;
}

# merge contigs by end to end overlapping (contigs completely encompassed by others will not be output in merged.fasta)
sub merge_contigs{
    my ($delta, $out, $ctg) = @_;
    my $tag = 0;
    # read the line right after >
    open fh_, "<$delta" or die $!;
    my $ctg1;
    my $ctg2;
    my $ctg1_len;
    my $ctg2_len;
    my $ctg_len_h;
    my $merged;
    my $delete;
    while(<fh_>){
        if($_ =~ />(\S+)\s+(\S+)\s+(\d+)\s+(\d+)/){
            if($tag == 0){
                $tag = 1;
            }
            $ctg1 = $1;
            $ctg2 = $2;
            $ctg1_len = $3;
            $ctg2_len = $4;
            $ctg_len_h->{$ctg1} = $ctg1_len;
            $ctg_len_h->{$ctg2} = $ctg2_len;
        }
        elsif($tag == 1){
            $tag = 0;
            next if($ctg1 eq $ctg2);
            my ($s1, $e1, $s2, $e2, ) = split(/\s+/, $_);
            next if($e1 - $s1 < $overlap_length);
            if($s1 < $tail_gap && $ctg2_len - $e2 < $tail_gap){
                # ctg1 head, ctg2 tail
                # head to tail, both are +
                $merged->{"${ctg2}_+:${ctg1}_${e1}_+"} = 1;
            }
            elsif($ctg1_len - $e1 < $tail_gap && $ctg2_len - $s2 < $tail_gap){
                # ctg1 tail, ctg2 tail
                $merged->{"${ctg1}_+:${ctg2}_${e2}_-"} = 1;
            }
            elsif($s1 < $tail_gap && $e2 < $tail_gap){
                # ctg1 head, ctg2 head
                $merged->{"${ctg1}_-:${ctg2}_${s2}_+"} = 1;
            }
            elsif($ctg1_len - $e1 < $tail_gap && $s2 < $tail_gap){
                # ctg1 tail, ctg2 head
                $merged->{"${ctg1}_+:${ctg2}_${e2}_+"} = 1;
            }
            elsif($s1 < $tail_gap && $ctg1_len - $e1 < $tail_gap){
                # ctg1 completely encompassed in ctg2, delete ctg1
                $delete->{"$ctg1:$ctg2"} = 1;
            }
            elsif($s2 < $tail_gap && $ctg2_len - $e2 < $tail_gap || $ctg2_len - $s2 < $tail_gap && $e2 < $tail_gap){
                # ctg2 compltely encompassed in ctg1, delete ctg2, be careful not to delete both
                $delete->{"$ctg2:$ctg1"} = 1 if(!defined $delete->{"$ctg1:$ctg2"});
            }
        }
    }
    close fh_;

    # analyze merged
    $tag = 0;
    my $h_forward;
    my $h_reverse;

    foreach my $k (keys %$merged){
        $tag = 0;
        my ($ctg1, $ori1, $ctg2, $bp2, $ori2) = ($k =~ /^(.+)_([+-]):(.+)_(\d+)_([+-])$/);
        # check if one of the contig is deleted, there should be better overlap
        foreach my $d (keys %$delete){
            my ($c1, $c2) = split(/:/, $d);
            if($c1 eq $ctg1 || $c1 eq $ctg2 || $c2 eq $ctg1 || $c2 eq $ctg2){
                $tag = 1;
                last;
            }
        }
        next if($tag == 1);
        # for one end of one ctg, if there is only one end for concatenating, concatenate it. Otherwise, try both, and see which one yield the longest. Record whatever is in the branches, so that those not involved will be output later.
        $h_forward->{"${ctg1}_$ori1"}->{"${ctg2}_${bp2}_$ori2"} = 1;
        # for checking if this is the head node
        $h_reverse->{"${ctg2}_$ori2"}->{$bp2}->{"${ctg1}_$ori1"} = 1;
    }

    # check all connections used
    my $checked;
    my @paths;
    foreach my $k (keys %$h_forward){
        # start with something on the end 
        my ($ctg, $ori) = split(/\_/, $k);
        my $ori_r = "-";
        $ori_r = "+" if($ori eq "-");
        if(!defined $h_reverse->{"${ctg}_$ori_r"}){
            my $path;
            # begin self-calling
            ($checked, $path) = &path_searching($k, $h_forward, $checked);
            push @paths, $path;
        }
    }

    # output the contig conrresponding to all path


    # check any contig not covered in checked
    my $paths_new = &check_missed(\@paths, $ctg_len_h);

    `rm $out` if(-e $out);
    # output these contigs
    #foreach my $p (@$paths_new, @paths){
    foreach my $p (@paths){
        foreach my $key (keys %$p){
            my $seq = &get_seq($key, $ctg);
            my $name = $key;
            my ($seqs, $names);
            ($seqs, $names) = &output_ctgs($seq, $p->{$key}, $ctg, $name, $seqs, $names); 
            &write_to_fa($seqs, $out, $names);
        }
    }
}

# write the array of sequences to an output file
sub write_to_fa{
    my ($seqs, $out, $names) = @_;
    my $n = scalar(@$seqs);
    open OUT_fh, ">>$out" or die $!;
    for(my $i = 0; $i < $n; $i ++){
        print OUT_fh join("\n", ">" . $names->[$i], $seqs->[$i]) . "\n";
    }
    close OUT_fh;
    return;
}
# output contigs according to the paths
sub output_ctgs{
    my ($seq, $p, $ctg, $name, $seqs, $names) = @_;

    if(ref($p) eq "HASH"){
        foreach my $key (keys %$p){
            print join("\t", "J", $key) . "\n";
            #print "$key\t$ctg\n";
            my $seq1 = &get_seq($key, $ctg);
            $seq .= $seq1;
            $name .= "#$key";
            ($seqs, $names) = &output_ctgs($seq, $p->{$key}, $ctg, $name, $seqs, $names);
        }
    }
    else{
        #print join("\t", "H" , $name) . "\n";
        push @$seqs, $seq;
        push @$names, $name;
    }
    return ($seqs, $names);
}

# get the sequence from the key
sub get_seq{
    my ($key, $ctg_f) = @_;
    my $seq;
    # TODO need to take care of the orientation
    # done
    # output both the first and next contigs
    my ($ctg, $p, $ori);
    if($key =~ /^(.+)_(.+)_(.+)$/){
        $ctg = $1;
        $p = $2;
        $ori = $3;
        if($ori eq "+"){
            # TODO: check if missing end is good
            my @tmps = split(/\n/, `samtools faidx $ctg_f $ctg:$p`);
            $seq = join("", @tmps[1 .. $#tmps]);
        }
        else{
            # TODO need to reverse complement
            # done
            my @tmps = split(/\n/, `samtools faidx $ctg_f $ctg:0-$p`);
            $seq = join("", @tmps[1 .. $#tmps]);
            $seq = seq::revcom($seq);
        }
    }
    elsif($key =~ /^(.+)_.+$/){
        # if the first is reverse, it does not have to be reversed
        $ctg = $1;
        my @tmps = split(/\n/, `samtools faidx $ctg_f $ctg`);
        $seq = join("", @tmps[1 .. $#tmps]);
    }

    return $seq;
}

# find if any contig does not appear in the path, add them to the path, single
sub check_missed{
    my ($paths, $ctg_h) = @_;
    #my $paths_new = $paths;
    my $paths_new;
    foreach my $path (@$paths){
        # path is a list in hash
        # traverse every key in this hash tree
        # again, recursively call itself
        $ctg_h = &check_path_missed($path, $ctg_h);
    }
    foreach my $key (keys %$ctg_h){
        if($ctg_h->{$key} != -1){
            $paths_new->{"${key}_+"} = 1;
        }
    }
    # to be in the same format as paths (array of reference)
    my @paths_return;
    push @paths_return, $paths_new;
    # return missed
    return \@paths_return;
}

# a recursive function used to augment check_missed 
sub check_path_missed{
    my ($path, $ctg_h) = @_;
    my $new_ctg_h;

    if(ref($path) eq "HASH"){
        foreach my $key (keys %$path){
            if($key !~ /\:/){
                if($key =~ /^(.+)_.+/){
                    $ctg_h->{$1} = -1;
                }
            }
            else{
                print "Error in line 1039\n";
            }
            my $path1 = $path->{$key};
            ($ctg_h) = &check_path_missed($path1, $ctg_h);
        }
    }
    return $ctg_h;
}

# traversal all possible paths
# praise the Lord!
sub path_searching{
    my ($k, $h, $checked) = @_;
    my $kk = $k;
    if($k =~ /^(.+)_.+_(.+)$/){
        $kk = join("_", $1, $2);
    }
    my $next;
    my $path;
    # use a temporary hash reference 
    my $tmp;
    ($next, $checked) = &find_next($kk, $h, $checked);
    foreach my $n (@$next){
        # the first node is in the format of ctg:ori; the next ones are in the format of ctg:pos:ori
        #my $key1;
        #if($n =~ /^(.+)\|.+\|(.+)$/){
        #    $key1 = join("|", $1, $2);
        #}
        my $path1;
        ($checked, $path1) = &path_searching($n, $h, $checked);
        # this is just for protection, should have keys always
        if(ref($path1) eq "HASH"){
            @$tmp{keys %$path1} = values %$path1;
        }
    }
    $path->{$k} = $tmp;
    if(scalar(@$next) == 0){
        $path->{$k} = 1;
    }
    return ($checked, $path);
}

# find the next node in merge
sub find_next{
    my ($key, $h, $checked) = @_;
    my @ks;
    my $key1 = $key;
    if($key =~ /^(.+)_.+_(.+)$/){
        $key1 = join("_", $1, $2);
    }
    if(defined $h->{$key1}){
        foreach my $k (keys %{$h->{$key1}}){
            # this actually should never been entered before, therefore if is only for checking
            if(!defined $checked->{"$key1:$k"}){
                push @ks, $k;
                $checked->{"$key1:$k"} = 1;
            }
        }
    }
    return (\@ks, $checked);
}

sub assemble{
    my ($fa, $dir) = @_;

    `nucmer -maxmatch -prefix $dir/out $fa $fa`;
    `delta-filter $dir/out.delta > $dir/out.rq.delta`;
    &merge_contigs("$dir/out.rq.delta", "$dir/merged.fasta", $fa);
    #`~/pkg/quickmerge/merger/quickmerge -d $dir/out.rq.delta -q $fa -r $fa -hco 5.0 -c 1.5 -l n -ml m`;
    #`mv merged.fasta $dir/`;

    #`$velveth $dir 21 -fasta -long $fa`;
    #`$velvetg $dir`;
    return "$dir/merged.fasta";
}

# only one cluster in the big cluster
sub align_IL2ctg_noCluster{
    my ($ctg_fa, $c, $ordered_IL) = @_;
    `rm $c.IL.1.fa`;
    `rm $c.IL.2.fa`;
    # get the Illumina reads
    print "get_IL_by_group.pl $c $ordered_IL $c.IL\n";
    `get_IL_by_group.pl $c $ordered_IL $c.IL`;
    my $nproc = 1;
    my $insert_size = 1000;
    my $output_dir = "clu$c";
    my ($ctg_prefix) = ($ctg_fa =~ /^(.+)\.fa/);

    `runBLASR.bigPipe.lite.pl -b $blasr_binary_new -c 1 -p $nproc -f sam -n 1 -i $insert_size -o $output_dir $c.IL $ctg_prefix`;
    print "runBLASR.bigPipe.lite.pl -b $blasr_binary_new -c 1 -p $nproc -f sam -n 1 -i $insert_size -o $output_dir $c.IL $ctg_prefix\n";
    # out.sam is the short read's alignment sam file
    return "$output_dir/out.sam";
}

# given a contig in ctg_fa, align the IL reads corresponding to the cluster listed in h_cluster ($h->{cluster}->{IL} are all the reads) to the contig, return the alignment. $c_large is the index of the large cluster
# now deal with multiple contigs in one contig file, in case not merged
sub align_IL2ctg{
    my ($ctg_fa, $h_cluster, $h, $ordered_IL, $c_large, $ctg_names) = @_;
    my $nproc = 1;
    my $insert_size = 1000;
    my $output_dir = "clu$c_large";
    `rm $c_large.IL.1.fa` if(-e "$c_large.IL.1.fa");
    `rm $c_large.IL.2.fa` if(-e "$c_large.IL.2.fa");
    `rm $output_dir/out.sam` if(-e "$output_dir/out.sam");
    foreach my $c (keys %$h_cluster){
        # get the Illumina reads
        print "get_IL_by_group.pl $c $ordered_IL $c_large.IL a\n";
        `get_IL_by_group.pl $c $ordered_IL $c_large.IL a`;
    }

    foreach my $ctg_nm (@$ctg_names){
        $DB::single = 1;
        `samtools faidx $ctg_fa $ctg_nm > $output_dir/$ctg_nm.fa`;
        `mv $output_dir/out.sam $output_dir/out.sam.bak` if(-e "$output_dir/out.sam");
        `runBLASR.bigPipe.lite.pl -b $blasr_binary_new -c 1 -p $nproc -f sam -n 1 -i $insert_size -o $output_dir $c_large.IL $output_dir/$ctg_nm`;
        # make the head
        `grep "^@" $output_dir/blasr_1.sam > $output_dir/head`;
        `cat $output_dir/head $output_dir/out.sam > $output_dir/out.sam.1`;
        `mv $output_dir/out.sam.1 $output_dir/out.sam`;
        `cat $output_dir/out.sam.bak >> $output_dir/out.sam` if(-e "$output_dir/out.sam.bak");
        print "runBLASR.bigPipe.lite.pl -b $blasr_binary_new -c 1 -p $nproc -f sam -n 1 -i $insert_size -o $output_dir $c_large.IL $ctg_nm\n";
        `rm $output_dir/out.sam.bak` if(-e "$output_dir/out.sam.bak");
        `rm $output_dir/$ctg_nm.fa`;
    }
    # out.sam is the short read's alignment sam file
    return "$output_dir/out.sam";
}

sub tree{
    my ($h) = @_;
    my $h_;
    foreach my $k (keys %$h){
        my $r = &root($k, $h);
        $h_->{$r}->{$k} = 1;
    }
    return $h_;
}

sub connect{
    my ($p, $q, $h) = @_;
    $h->{$p} = $q;
    return $h;
}

sub root{
    my ($p, $h) = @_;
    if(!defined $h->{$p}){
        $h->{$p} = $p;
        return $p;
    }
    while($h->{$p} != $p){
        $p = $h->{$p};
    }
    return $p;
}
1;
