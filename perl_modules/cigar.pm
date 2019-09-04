use warnings;
use strict;
package cigar;

# analyze cigar string

sub reverse_cigar{
    my ($cigar) = @_;
    my ($tag, $num) = &analyze_cigar($cigar);
    my $str = "";
    for(my $i = scalar(@$tag) - 1; $i >= 0; $i --){
        $str .= $num->[$i] . $tag->[$i];
    }
    return $str;
}
        
sub check_cigar{
    my ($cigar, $str, $n) = @_;
    # This is the main entrance of this package. It calls numerous small funcs. It serves for other programs.
    # The string can be "INS", "DEL" followed by the size, returns the number of contigs with insertion/deletion above that size. It also can be ref_len, returning the number of bases covered by the read on the reference.
    my ($tag, $num) = &analyze_cigar($cigar);
    if($str eq "INS"){
        my $ins = 0;
        for(my $i = 0; $i < scalar(@$tag); $i ++){
            if($num->[$i] > $n && $tag->[$i] eq "I"){
                $ins ++;
            }
        }
        return $ins;
    }
    if($str eq "DEL"){
         my $del = 0;
        for(my $i = 0; $i < scalar(@$tag); $i ++){
            if($num->[$i] > $n && $tag->[$i] eq "D"){
                $del ++;
            }
        }
        return $del;
    }
    if($str eq "ref_len"){
        my $ref_len = 0;
        for(my $i = 0; $i < scalar(@$tag); $i ++){
            if($tag->[$i] =~ /D|M/){
                $ref_len += $num->[$i];
            }
        }
        return $ref_len;
    }
}
# given a line of sam with long cigar string, print an array, each line with ref_pos, cigar_tag, cigar_len
# taking into actual read pos (revcom)
# note: if read is SH, ref does not account for that, meaning, the same position for the region of S/H
sub get_pos_from_cigar_input_consider_rev{
    my ($line) = @_;
    my @a = split(/\t/, $line);
    my $cigar = $a[5];
    my $chr = $a[2];
    my $pos = $a[3];
    my $r = 0;
    if($a[1] & 0x0010){
        $r = 1;
    }
    my $total_len = &get_total_len($cigar);
    my $p = $pos;
    my $p_l = 0;
    # return four pointers to hashes
    my $p_ref;
    my $p_read;
    my ($tag, $num) = &analyze_cigar($cigar);
    for(my $i = 0; $i < scalar(@$tag); $i ++){
        if($i == 0){
            push @$p_ref, $p;
            if($r == 1){
                push @$p_read, ($total_len - $p_l);
            }
            else{
                push @$p_read, $p_l;
            }
        }
        if($tag->[$i] =~ /M/){
            $p += $num->[$i];
            $p_l += $num->[$i];
        }
        if($tag->[$i] =~ /D/){
            $p += $num->[$i];
        }
        if($tag->[$i] =~ /[SHI]/){
            $p_l += $num->[$i];
        }
        # push the positions to array
        push @$p_ref, $p;
        if($r == 1){
            push @$p_read, ($total_len - $p_l);
        }
        else{
            push @$p_read, $p_l;
        }
    }
    return ($p_ref, $p_read, $tag, $num);
}

# get a substring of cigar with the index region
sub get_substr{
    my ($cigar, $a, $b) = @_;
    my ($tag, $num) = &analyze_cigar($cigar);
    my $str = "";
    for(my $i = 0; $i < scalar(@$tag); $i ++){
        if($i <= $b && $i >= $a){
            $str .= $num->[$i] . $tag->[$i];
        }
    }
    return $str;
}

# given a line of sam with long cigar string, print an array, each line with ref_pos, cigar_tag, cigar_len
sub get_pos_from_cigar_input{
    my ($line) = @_;
    my @a = split(/\t/, $line);
    my $cigar = $a[5];
    my $chr = $a[2];
    my $pos = $a[3];
    my $p = $pos;
    my $p_l = 0;
    # return four pointers to hashes
    my $p_ref;
    my $p_read;
    my ($tag, $num) = &analyze_cigar($cigar);
    for(my $i = 0; $i < scalar(@$tag); $i ++){
        if($tag->[$i] =~ /SH/ && $i == 0){
            push @$p_ref, $p - $num->[$i];
            push @$p_read, $p_l;
            #print join("\t", $p - $num->[$i], $p_l, $tag->[$i], $num->[$i]). "\n";
        }
        else{
            push @$p_ref, $p;
            push @$p_read, $p_l;
            #print join("\t", $p, $p_l, $tag->[$i], $num->[$i]) . "\n";
        }
        if($tag->[$i] =~ /M/){
            $p += $num->[$i];
            $p_l += $num->[$i];
        }
        if($tag->[$i] =~ /D/){
            $p += $num->[$i];
        }
        if($tag->[$i] =~ /[SHI]/){
            $p_l += $num->[$i];
        }
    }
    return ($p_ref, $p_read, $tag, $num);
}



# given a line of sam with long cigar string, print an array, each line with ref_pos, cigar_tag, cigar_len
sub get_pos_from_cigar{
    my ($line) = @_;
    my @a = split(/\t/, $line);
    my $cigar = $a[5];
    my $chr = $a[2];
    my $pos = $a[3];
    my $p = $pos;
    my $p_l = 0;
    my ($tag, $num) = &analyze_cigar($cigar);
    for(my $i = 0; $i < scalar(@$tag); $i ++){
        if($tag->[$i] =~ /SH/ && $i == 0){
            print join("\t", $p - $num->[$i], $p_l, $tag->[$i], $num->[$i]). "\n";
        }
        else{
            print join("\t", $p, $p_l, $tag->[$i], $num->[$i]) . "\n";
        }
        if($tag->[$i] =~ /M/){
            $p += $num->[$i];
            $p_l += $num->[$i];
        }
        if($tag->[$i] =~ /D/){
            $p += $num->[$i];
        }
        if($tag->[$i] =~ /[SHI]/){
            $p_l += $num->[$i];
        }
    }
}


# This function analyzes the cigar string of an IL to a PB ctg on the putative INDEL areas of PB ctg. If the comparison between IL and PB reports no INDEL, then the putative INDEL is likely real, otherwise, the PB ctg has indel errors. 
# h is a hash, with keys the integer counting from 0, each pointing to a data structure with of the putative INDEL: chr, pos, len, tp (I/D), ctg_name, ctg_start, ctg_end (INDEL position on the ctg), spt_short (to be added one if this function verify the INDEL with the IL read), line (the putative INDEL line).
# ctg is the contig of where the IL aligns, ctg_p is the pos of the actual position the IL aligns to. cigar is the one to be processed, the actual alignment of the IL to the ctg. Note that due to blasr clipping, the total length of the cigar string might be well below the actual length of the read. 
# perc is the required percentage of alignment (M only) from IL to the surrounding area of the INDEL of ctg to count this IL as supporting the INDEL. The surrounding area is just the inserted sequence if I for tp, or [p-s, p+s] if tp is D, where p is the position of the deletion on ctg, s is given.
sub test_indel{
    my ($h, $ctg, $ctg_p, $cigar, $perc, $s) = @_;
    # c is returned to report whether this read is counted as supportive or not, for a paired end analysis. 0: no; 1: yes
    my $c = 0;

    if($h->{ctg_name} eq $ctg){
        my $begin = $h->{ctg_start};
        my $end = $h->{ctg_end};
        my $tp = $h->{tp};
        if($tp eq "I"){
            # insertion
            my $range = "$begin:$end";
            $c = &test_perc_identity($range, $ctg_p, $cigar, $perc, $tp);
            return $c;
        }
        else{
            # including small deletion, deletion, inversion, duplication, CTX, ITX
            # find the corresponding position on the read
            my $tmp = &find_pos_on_read($ctg_p, $begin, $cigar);
            return 0 if($tmp < 0);
            my $tmp1 = $tmp - $s;
            my $tmp2 = $tmp + $s;
            my $range = "$tmp1:$tmp2";
            $c = &test_perc_identity($range, $ctg_p, $cigar, $perc, $tp);
            return $c;
        }
    }
    return $c;

}

# This function find the corresonding position on the read given a position on the reference, and the cigar string, only the matched part (M, I, D), thus not including the clipped part (H, S)
# p_ref is the starting position of the read on the ctg, and the start is the SV start on the ctg
# return the breakpoint position on the read, excluding the head S/H part
sub find_pos_on_read{
    my ($p_ref, $start, $cigar) = @_;
    # add a protection and saves computation time when the read starts from the coordiante greater than the breakpoint 09142015
    return -1 if($p_ref > $start);
    # p_read is the position of the read base on the reference
    my $p_read = -1;
    # s_read is the position of the read base on the reference, start
    my $s_read = $p_ref;
    # p is the corresonding position on the read to p_read
    my $p = 0;
    my ($tag, $num) = &analyze_cigar($cigar);
    foreach my $i (0 .. scalar(@$tag) - 1){
        next if($tag->[$i] =~ /[SH]/);
        if($tag->[$i] =~ /M/){
            $p_read = $s_read + $num->[$i];
            if($p_read > $start && $start >= $s_read){
                return $p + $start - $s_read;
            }
            else{
                # not yet hit
                $p += $num->[$i];
                $s_read = $p_read;
            }
        }
        elsif($tag->[$i] =~ /D/){
            if($p_read != $start){
                # go through the reference
                $p_read = $s_read + $num->[$i];
                $s_read = $p_read;
            }
            else{
                return $p;
            }
        }
        elsif($tag->[$i] =~ /I/){
            if($p_read != $start){
                # go through the read
                $p += $num->[$i];
            }
            else{
                return $p;
            }
        }
    }
    return -1;
}



# This function tests whether the percentage of identity from a cigar string is over a threshold. Return 0 if no, 1 if yes. The range is where the test is performed. p is the starting position of the cigar string (not including ^S|H), tp is the type of indel (I/D). If I, the range is w.r.t. ctg, if D, the range is w.r.t. IL. 
sub test_perc_identity{
    my ($range, $p, $cigar, $perc, $tp) = @_;
    my $min_req = 3;
    my $c = 0;
    my ($tag, $num) = &analyze_cigar($cigar);
    # look at the coverage of the read on the ctg
    my $len = 0;
    for(my $i = 0; $i < scalar(@$tag); $i ++){
        if($tag->[$i] =~ /[MD]/){
            $len += $num->[$i];
        }
    }
    my $p_end = $p + $len;

    my ($b, $e) = split(/:/, $range);

    # for insertion, only look at the range where the read covers
    if($tp eq "I"){
        # ignore those don't overlap
        if($p_end < $b || $e < $p){
            return $c;
        }
        # begin from whichever comes later, and ends at whichever comes first
        my $b_tmp = $b;
        my $e_tmp = $e;
        $b = $p if($b_tmp < $p);
        $e = $p_end if($e_tmp > $p_end);
        return $c if($b >= $e);
    }
    my $matchLen = 0;
    my $hh;
    if($tp eq "I"){
        # look at the range of the ctg
        $hh = &matching_status_on_ref($p, $cigar, $min_req);
    }
    elsif($tp eq "D"){
        # look at the range of the IL
        $hh = &matching_status_on_read($cigar, $min_req);
    }
    foreach my $pos ($b .. $e){
        if(defined $hh->{$pos} && $hh->{$pos} eq "M"){
            $matchLen ++;
        }
    }
    if($matchLen * 100/($e - $b) > $perc){
        $c = 1;
    }
    return $c;
}
# record the matching status (M, I, S, H, D) on every base on the read given starting point p (excluding H/S) and cigar string
sub matching_status_on_read{
    my ($cigar, $min_req) = @_;
    my ($tag, $num) = &analyze_cigar($cigar);
    my $h;
    my $p = 0;
    foreach my $i (0 .. scalar(@$tag) - 1){
        if($tag->[$i] =~ /[IM]/){
            foreach my $t ($p .. $p + $num->[$i] - 1){
                # require the size to be greater than min_req, otherwise status be NA. 09142015
                if($num->[$i] > $min_req){
                    $h->{$t} = $tag->[$i];
                }
                else{
                    $h->{$t} = "N";
                }
            }
            $p += $num->[$i];
        }
    }
    return $h;
}


# record the matching status (M, I, S, H, D) on every base on the reference given starting point p (excluding H/S) and cigar string
sub matching_status_on_ref{
    my ($p, $cigar, $min_req) = @_;
    my ($tag, $num) = &analyze_cigar($cigar);
    my $h;
    foreach my $i (0 .. scalar(@$tag) - 1){
        if($i == 0 && $tag->[$i] =~ /[SH]/){
            my $pos = $p - $num->[$i];
            foreach my $t ($pos .. $p - 1){
                $h->{$t} = $tag->[$i];
            }
        }
        else{
            if($tag->[$i] =~ /[DSHM]/){
                foreach my $t ($p .. $p + $num->[$i] - 1){
                    # require the size to be greater than min_req, otherwise status be NA. 09142015
                    if($num->[$i] > $min_req){
                        $h->{$t} = $tag->[$i];
                    }
                    else{
                        $h->{$t} = "N";
                    }
                }
                $p += $num->[$i];
            }
        }
    }
    return $h;
}

sub analyze_cigar{
    # return two arrays with tag and length for each, corresponding to each other in the order.
    my ($cigar) = @_;
    my (@tag, @num);
    return (\@tag, \@num) if (!defined $cigar || $cigar eq "" || $cigar eq "*");
    while($cigar ne ""){
        my ($t, $n);
        ($n, $t, $cigar) = ($cigar =~ /^(\d+)(\S)(.*)$/);
        push @tag, $t;
        push @num, $n;
    }
    return (\@tag, \@num);
}

# for each read in a sam, compute the total length on the reference, print the start and end of the ref
sub sum_on_ref{
    my ($sam, $readname) = @_;
    open S, "<$sam" or die $!;
    while(<S>){
        next if($_ =~ /^@/);
        my $len = 0;
        my @a = split(/\t/, $_);
        my ($tag, $num) = &analyze_cigar($a[5]);
        for(my $i = 0; $i < scalar(@$tag); $i ++){
            if($tag->[$i] =~ /M|D/){
                $len += $num->[$i];
            }
        }
        print join("\t",  $a[0], $len, $a[2], $a[3], $a[3] + $len) ."\n";
    }
    close S;
}
# add SV flavor to this package.
# given a cigar, position on the reference, return the breakpoints at the beginning and the end if the soft clipped end length is greater the threshold. return -1 if none.
sub bp_from_cigar{
    # This is mainly for Pacbio related sequence
    my ($cigar, $pos, $clip_t) = @_;
    my ($tag, $num) = &analyze_cigar($cigar);
    my $t = scalar(@$tag);
    my $bp1 = -1;
    my $bp2 = -1;
    if($tag->[0] =~ /[SH]/ && $num->[0] > $clip_t){
        $bp1 = $pos;
    }
    if($tag->[$t - 1] =~ /[SH]/ && $num->[$t - 1] > $clip_t){
        $bp2 = $pos;
        foreach (0 .. $t - 1 - 1){
            $bp2 += $num->[$_] if($tag->[$_] =~ /[MD]/);
        }
    }
    return ($bp1, $bp2);
}

# given a cigar, position on the reference, return the indel information (in the format of type:pos:size, for example, for insertion of size 10 at 500, it is I:500:10). If there are multiple INDELs, report the one with the largest size. Return NA if none.
sub indel_from_cigar{
    # This is mainly for Illumina sequence, i.e., not designed for reads containing many indels.
    my ($cigar, $pos, $indel_t) = @_;
    my ($tag, $num) = &analyze_cigar($cigar);
    my $ret = "NA";
    my $b = $pos;
    my $max_size = -1;
    my $t = scalar(@$tag);
    foreach my $i (0 .. $t - 1){
        if($tag->[$i] =~ /[ID]/ && $num->[$i] > $indel_t){
            if($num->[$i] > $max_size){
                $ret = join(":", $tag->[$i], $b, $num->[$i]);
                $max_size = $num->[$i];
            }
        }
        # update position
        if($tag->[$i] =~ /[MD]/){
            $b += $num->[$i];
        }
    }
    return $ret;
}

# check if the bp on ref is encompassed by the matched read
sub check_crossing_bp{
    my ($line, $bp) = @_;
    my @a = split(/\t/, $line);
    my $len_ref = &get_matchlen_on_ref($a[5]);
    if($bp > $a[3] && $bp < $a[3] + $len_ref){
        return 1;
    }
    return 0;
}


# check if this cigar support the soft clip on the breakpoint on ori (-: 3' clipping, +: 5' clipping)
sub check_breakpoint{
    my ($cigar, $pos, $bp, $ori, $off) = @_;
    return 0 if($cigar !~ /S/);
    if($cigar =~ /S$/ && $ori eq "+" || $cigar =~ /^\d+S/ && $ori eq "-"){
        return 0;
    }
    # check position
    my $bp1;
    my $bp2;
    ($bp1, $bp2) = &bp_from_cigar($cigar, $pos, 0);
    if($ori eq "+" && $bp1 != -1){
        if(abs($bp1 - $bp) < $off){
            return $bp1;
        }
        else{
            return 0;
        }
    }
    elsif($ori eq "-" && $bp2 != -1){
        if(abs($bp2 - $bp) < $off){
            return $bp2;
        }
        else{
            return 0;
        }
    }
    return 0;
}

sub get_total_len{
    my ($cigar) = @_;
    my ($tag, $num) = &analyze_cigar($cigar);
    my $len = 0;
    for(my $i = 0; $i < scalar(@$tag); $i ++){
        if($tag->[$i] =~ /[IMSH]/){
            $len += $num->[$i];
        }
    }
    return $len;
}

# get the match len on read and the total length of the reads, not counting hard clip (works for Pacbio with clipping soft)
sub get_matchlen_on_read_of_all{
    my ($cigar) = @_;
    my ($tag, $num) = &analyze_cigar($cigar);
    my $len = 0;
    my $totallen = 0;
    for(my $i = 0; $i < scalar(@$tag); $i ++){
        if($tag->[$i] =~ /[IM]/){
            $len += $num->[$i];
        }
        if($tag->[$i] =~ /[IMSH]/){
            $totallen += $num->[$i];
        }
    }
    return join("\t", $len, $totallen);
}


sub get_matchlen_on_read{
    my ($cigar) = @_;
    my ($tag, $num) = &analyze_cigar($cigar);
    my $len = 0;
    for(my $i = 0; $i < scalar(@$tag); $i ++){
        if($tag->[$i] =~ /[IM]/){
            $len += $num->[$i];
        }
    }
    return $len;
}

sub get_matchlen_on_ref{
    my ($cigar) = @_;
    my ($tag, $num) = &analyze_cigar($cigar);
    my $len = 0;
    for(my $i = 0; $i < scalar(@$tag); $i ++){
        if($tag->[$i] =~ /[DM]/){
            $len += $num->[$i];
        }
    }
    return $len;
}

# relative to the strand of the read
# if cigar is obtained when read is reversed, start is from the end point
# rev is -1 if reverse, 1 if forward
sub get_startpos_on_read{
    my ($cigar, $rev) = @_;
    my $pos = 0;
    my ($tag, $num) = &analyze_cigar($cigar);
    if($rev == 1){
        for(my $i = 0; $i < scalar(@$tag); $i ++){
            if($tag->[$i] =~ /[SH]/){
                $pos += $num->[$i];
                next;
            }
            else{
                return $pos;
            }
        }
    }
    else{
        for(my $i = scalar(@$tag) - 1; $i >= 0; $i --){
            if($tag->[$i] =~ /[SH]/){
                $pos += $num->[$i];
                next;
            }
            else{
                return $pos;
            }
        }
    }
}

sub get_startpos_on_ref{
    my ($cigar, $pos) = @_;
    return $pos;
}

# get the end position (hard stopped) on reference
sub get_endpos_on_ref{
    my ($cigar, $pos) = @_;
    my ($tag, $num) = &analyze_cigar($cigar);
    for(my $i = scalar(@$tag) - 1; $i >= 0; $i --){
        if($tag->[$i] =~ /[MD]/){
            $pos += $num->[$i];
        }
    }
    return $pos;
}


sub get_endpos_on_read{
    my ($cigar, $rev) = @_;
    return &get_startpos_on_read($cigar, $rev) + &get_matchlen_on_read($cigar);
}

# Given the start (from soft-clipped or hard-clipped position) on ref, the cigar string and a region of interest, output the corresponding cigar string
# return a string with two lines separated by ;, the first line is the ins sizes (>$min_size) separated by :, the second line with corresponding deletion sizes, with each element the size of the indel followed by the position on the reference (start) separated by "|".
sub get_cigar_from_region_on_ref{
    my ($ref_region, $read_start, $cigar, $min_size) = @_;
    my ($start, $end) = split(/-/, $ref_region);
    return -1 if($read_start > $end);
    my @del;
    my @ins;
    # p_read is the position of the read base on the reference
    my $p_read = $read_start;
    # p is the corresonding position on the read to p_read
    my $p = 0;
    my ($tag, $num) = &analyze_cigar($cigar);
    # tag = 1 when starting to count
    foreach my $i (0 .. scalar(@$tag) - 1){
        next if($tag->[$i] =~ /[SH]/);
        if($tag->[$i] =~ /M/){
            $p_read += $num->[$i];
            if($p_read > $end){
                last;
            }
        }
        elsif($tag->[$i] =~ /D/){
            if($p_read + $num->[$i] < $end && $p_read + $num->[$i] > $start || $p_read > $start && $p_read < $end){
                if($num->[$i] > $min_size){
                    push @del, $num->[$i] . "|" . $p_read;
                }
            }
            $p_read += $num->[$i];
            if($p_read > $end){
                last;
            }
        }
        elsif($tag->[$i] =~ /I/){
            if($p_read > $start && $p_read < $end){
                if($num->[$i] > $min_size){
                    push @ins, $num->[$i] . "|" . $p_read;
                }
            }
        }
    }
    if($p_read < $start){
        return -1;
    }
    my $str = "NA";
    if(@del == 0){
        push @del, "NA";
    }
    if(@ins == 0){
        push @ins, "NA";
    }
    $str = join(";", join(":", @ins), join(":", @del));
    return $str;
}

# this function looks at the cigar string, return the left and right clip (both hard and soft clip), 0 if none
sub get_clip{
    my ($cigar) = @_;
    my $head = 0;
    my $tail = 0;
    if($cigar =~ /^(\d+)[SH]/){
        $head = $1;
    }
    if($cigar =~ /(\d+)[SH]$/){
        $tail = $1;
    }
    return ($head, $tail);
}

1;
