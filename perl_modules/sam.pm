use warnings;
use strict;
package sam;
require cigar;
require asm;
require seq;
require hmm;
require math;
require fa;

# note: M_threshold is disabled. As long as the indel is within the cigar (not the first or the alst), and > indel_threshold, it is counted
sub analyze_simple_sam{
    my ($sam, $M_threshold, $indel_threshold) = @_;

    my @b;
    open fh_, "<$sam" or die $!;
    while(<fh_>){
        my $a = $_;
        if($a =~ /^@/){
            next;
        }
        # only analyze the first line
        my @items = split(/\t/, $a);
        my $read_name = $items[0];
        my $ref_name = $items[2];
        my $r = 0;
        if($items[1] & 0x0010){
            $r = 1;
        }

        my ($chr, $start, $end) = ($ref_name =~ /^(\S+)\:(\d+)\-(\d+)$/);
        next if(!defined $start || !defined $end); 
        my ($ctg_name, $ctg_start, $ctg_end) = ($read_name =~ /^(\S+)\:(\d+)\-(\d+)$/);
        #print $_ . "\t" . $ref_name . "\t" . $read_name . "\n";
        my $cigar = $items[5];
        my ($alpha, $num) = cigar::analyze_cigar($cigar);
        my $total = scalar(@$num);
        my $ctg_length = $ctg_end - $ctg_start;
        my $ref_pos;
        my $ctg_pos1;
        my $ctg_pos2;
        my $ref_MD = 0;
        my $ctg_MI = 0;
        # suppose there is only one large INS or DEL, and the two flanking regions (M on the two sides of the SV) are long enough, call the variant
        for(my $i = 0; $i < $total; $i ++){
            if($alpha->[$i] =~ /[ID]/ && $num->[$i] > $indel_threshold){
                if($i != 0 && $i != $total - 1){
                    #if($alpha->[$i-1] eq "M" && $num->[$i-1] > $M_threshold && $alpha->[$i+1] eq "M" && $num->[$i+1] > $M_threshold){
                    $ref_pos = $start + $ref_MD;
                    if($r == 0){
                        $ctg_pos1 = $ctg_start + $ctg_MI;
                        if($alpha->[$i] eq "D"){
                            $ctg_pos2 = $ctg_pos1 + 1; 
                        }
                        elsif($alpha->[$i] eq "I"){
                            $ctg_pos2 = $ctg_pos1 + $num->[$i];
                        }
                    }
                    else{
                        if($alpha->[$i] eq "D"){
                            $ctg_pos1 = $ctg_length - $ctg_MI + $ctg_start; 
                            $ctg_pos2 = $ctg_pos1 + 1;
                        }
                        elsif($alpha->[$i] eq "I"){
                            $ctg_pos1 = $ctg_length - $ctg_MI - $num->[$i] + $ctg_start;
                            $ctg_pos2 = $ctg_length - $ctg_MI + $ctg_start;
                        }
                    }
                    my $b_ = join("\t", $chr, $ref_pos, $num->[$i], $alpha->[$i], $ctg_name, $ctg_pos1, $ctg_pos2);
                    push @b, $b_;
                    #}
                }
            }
            if($alpha->[$i] =~ /[MD]/){
                $ref_MD += $num->[$i];
            }
            if($alpha->[$i] =~ /[MI]/){
                $ctg_MI += $num->[$i];
            }
        }
    }
    close fh_;

    return \@b;
}


# analyze simple cigar string, should be large matching followed by large deletion or insertion, followed by large matching
# input is a string in sam format, may contain header
sub analyze_simple_sam_deprecate{
    my ($sam, $sz_M, $length_diff_t, $r) = @_;
    my @as = split(/\n/, $sam); 
    my $b = "";
    foreach my $a (@as){
        if($a =~ /^@/){
            next;
        }
        # only analyze the first line
        my @items = split(/\t/, $a);
        my $read_name = $items[0];
        my $ref_name = $items[2];
        my ($chr, $start, $end) = ($ref_name =~ /^(\S+)\:(\d+)\-(\d+)$/);
        my ($ctg_name, $ctg_start, $ctg_end) = ($read_name =~ /^(\S+)\:(\d+)\-(\d+)$/);
        my $cigar = $items[5];
        my ($alpha, $num) = cigar::analyze_cigar($cigar);
        my $total = scalar(@$num);
        my $max = 0;
        my $max_alpha = "";
        my $max_ref_pos;
        my $max_ctg_pos;
        my $ref_middle_len = 0;
        my $ctg_middle_len = 0;
        if($num->[0] > $sz_M && $num->[$total-1] > $sz_M && $alpha->[0] eq "M" && $alpha->[$total-1] eq "M"){
            # check the largest in between the two large matching
            for(my $i = 1; $i < $total - 1; $i ++){
                if($num->[$i] > $max){
                    $max = $num->[$i];
                    $max_alpha = $alpha->[$i];
                    $max_ref_pos = $ref_middle_len + $num->[0];
                    $max_ctg_pos = $ctg_middle_len + $num->[0];
                }
                if($alpha->[$i] =~ /[MD]/){
                    $ref_middle_len += $num->[$i];
                }
                if($alpha->[$i] =~ /[MI]/){
                    $ctg_middle_len += $num->[$i];
                }
            }
            # the middle length should not be too far from the total length between the two large matches
            $max_ref_pos += $start;
            if($r == 0){
                $max_ctg_pos += $ctg_start;
            }
            else{
                $max_ctg_pos = $ctg_end - $max_ctg_pos;
            }
            if($max_alpha eq "D" && $max > $ref_middle_len - $length_diff_t){
                if($r == 1){
                    $max_ctg_pos -= 1;
                }
                $b = join("\t", $chr, $max_ref_pos, $max, "D", $ctg_name, $max_ctg_pos, $max_ctg_pos+1);
            }
            elsif($max_alpha eq "I" && $max > $ctg_middle_len - $length_diff_t){
                if($r == 1){
                    $max_ctg_pos -= $max;
                }
                $b = join("\t", $chr, $max_ref_pos, $max, "I", $ctg_name, $max_ctg_pos, $max_ctg_pos+$max);
            }
        }
        last;
    }
    return $b;
}


# make a sam file from the paired HMM cigar (I and D might appear at the end and beginning of a string)
sub make_sam_from_pairHMMcigar{
    my ($ctg_str, $ref_str, $revcom, $cigar, $sam_file, $ctg_name_full) = @_;
    my ($ctg_name, $ctg_start, $ctg_end) = split(/_/, $ctg_str);
    my ($ref_name, $ref_start, $ref_end) = split(/_/, $ref_str);
    my $pos;
    if($revcom == 1){
        $cigar = cigar::reverse_cigar($cigar);
    }
    if($cigar =~ /^(\d+)I(.+)$/){
        # use the second part for cigar
        $cigar = $2;
        if($revcom == 0){
            $ctg_start += $1;
        }
        else{
            $ctg_end -= $1;
        }
    }
    elsif($cigar =~ /^(\d+)D(.+)$/){
        $cigar = $2;
        $ref_start += $1;
    }
    if($cigar =~ /^(.+[MIDSH])(\d+)I$/){
        $cigar = $1;
        if($revcom == 0){
            $ctg_end -= $2;
        }
        else{
            $ctg_start += $2;
        }
    }
    elsif($cigar =~ /^(.+[MIDSH])(\d+)D$/){
        $cigar = $1;
        $ref_end -= $2;
    }
    $ref_str = "$ref_name:$ref_start-$ref_end";
    $ctg_str = "$ctg_name_full:$ctg_start-$ctg_end";
    my $str = join("\t", $ctg_str, $revcom == 1 ? 16 : 8, $ref_str, 0, "60", $cigar, "*", "*", "*", "*", "*", "*") . "\n";
    open fh_, ">$sam_file" or die $!;
    print fh_ $str;
    close fh_;
}


# some tools related with sam file format
# analyze a sam, find the best match (alignment with the highest match bases) and return it or NA if all alignments are bad (at least > threshold of bases are not matched)
sub analyze_sam{
    my ($sam, $t) = @_;
    open SAM, "<$sam" or die $!;
    while(<SAM>){
        next if($_ =~ /^@/);
        my @a = split(/\t/, $_);
        my ($alpha, $num) = cigar::analyze_cigar($a[5]);
        my $match_len = 0;
        my $total_len = 0;
        my $ref_len = 0;
        foreach my $i (0 .. scalar(@$alpha) - 1){
            if($alpha->[$i] eq "M"){
                $match_len += $num->[$i];
            }
            $total_len += $num->[$i] if($alpha->[$i] =~ /[MISH]/);
            $ref_len += $num->[$i] if($alpha->[$i] =~ /[MD]/);
        }
        if($match_len/$total_len > $t){
            close SAM;
            return $a[2] . ":" .  $a[3] . "-" . ($a[3] + $ref_len);
        }
    }
    close SAM;
    return "NA";
}

# filter out short read pairs that have 83, 163 , or 147, 99 flags, and their start distance < 50
sub filter_reads_combinePBIL{
    my ($sam, $dist_t) = @_;
    my $readname = "NA";
    my $flag;
    my $chr;
    my $pos;
    my $line = "";
    open SAM, "<$sam" or die $!;
    #$DB::single = 1;
    while(<SAM>){
        next if($_ =~ /^@/);
        my @a = split(/\t/, $_);
        if($readname eq $a[0]){
            # a pair
            if($flag == 99 && $a[1] == 147 || $flag == 147 && $a[1] == 99 || $flag == 163 && $a[1] == 83 || $a[1] == 163 && $flag == 83){
                # check distance
                if($a[2] eq $chr && abs($a[3] -$pos) < $dist_t){
                    next;
                }
                else{
                    print $line, $_;
                }
            }
            else{
                #    print $line, $_;
            }
        }
        else{
            $readname = $a[0];
            $chr = $a[2];
            $flag = $a[1];
            $pos = $a[3];
            $line = $_;
        }
    }
    close SAM;
}



sub sam2fq_noPair{
    # convert sam to fq. If the quality col is *, then make it to be a high qual, like I for all bases.
    my ($sam) = @_;
    my ($fq) = ($sam =~ /^(.+)\.sam/);
    $fq .= ".fq";
    open sam_fh, "<$sam" or die $!;
    open fq_out, ">$fq" or die $!;
    my $n = 0;
    while(<sam_fh>){
        next if($_ =~ /^@/);
        my @a = split(/\t/, $_);
        my $qual = "I" x length($a[9]);
        if($a[10] ne "*"){
            $qual = $a[10];
        }
        print fq_out join("\n", "@" . $n, $a[9], "+", $qual) . "\n";
        $n ++;
    }
    close fq_out;
    close sam_fh;
}
sub sam2fq{
    # convert sam to fq. If the quality col is *, then make it to be a high qual, like I for all bases.
    my ($sam) = @_;
    my ($fq) = ($sam =~ /^(.+)\.sam/);
    $fq .= ".fq";
    open sam_fh, "<$sam" or die $!;
    open fq_out, ">$fq" or die $!;
    while(<sam_fh>){
        next if($_ =~ /^@/);
        my @a = split(/\t/, $_);
        my $qual = "I" x length($a[9]);
        if($a[10] ne "*"){
            $qual = $a[10];
        }
        print fq_out join("\n", "@" . $a[0], $a[9], "+", $qual) . "\n";
    }
    close fq_out;
    close sam_fh;
}
sub sam2fq_pair{
    # convert a sam to a pair of fqstq files, leave the singles not printed. fq name would be fq_1.fqstq, fq_2.fqstq
    my ($sam) = @_;
    my ($fq) = ($sam =~ /^(.+)\.sam/);
    my $fq1 = $fq. "_1.fq";
    my $fq2 = $fq. "_2.fq";
    my $h;
    open sam_fh, "<$sam" or die $!;
    open fq_out1, ">$fq1" or die $!;
    open fq_out2, ">$fq2" or die $!;
    while(<sam_fh>){
        next if($_ =~ /^@/);
        my @a = split(/\t/, $_);
        print $_ if(!defined $a[9] || $a[9] eq "");
        # remove those that are not primary alignment
        next if($a[1] & 0x0800);
        # convert from reverse to forward strand
        if($a[1] & 0x0010){
            $a[9] = seq::revcom($a[9]);
            $a[10] = scalar reverse($a[10]);
        }
        if(!defined $h->{$a[0]}){
            $h->{$a[0]} = join("\t", @a);
            next;
        }
        # print only in pair
        if($a[1] & 0x0040){
            # this is the first, save to fq1
            print fq_out1 join("\n", "@" . $a[0], $a[9], "+", $a[10])."\n";
            my @b = split(/\t/, $h->{$a[0]});
            print fq_out2 join("\n", "@" . $b[0], $b[9], "+", $b[10])."\n";
            delete $h->{$a[0]};
        }
        elsif($a[1] & 0x0080){
            # this it the second
            print fq_out2 join("\n", "@" . $a[0], $a[9], "+", $a[10])."\n";
            my @b = split(/\t/, $h->{$a[0]});
            print fq_out1 join("\n", "@" . $b[0], $b[9], "+", $b[10])."\n";
            delete $h->{$a[0]};
        }
    }
    close fq_out1;
    close fq_out2;
    close sam_fh;
}
sub sam2fa_pair_wHeader{
    # convert a sam to a pair of fastq files, leave the singles not printed. fa name would be fa_1.fastq, fa_2.fastq
    my ($sam) = @_;
    my ($fa) = ($sam =~ /^(.+)\.sam/);
    my $fa1 = $fa. "_1.fa";
    my $fa2 = $fa. "_2.fa";
    my $h;
    open sam_fh, "<$sam" or die $!;
    open fa_out1, ">$fa1" or die $!;
    open fa_out2, ">$fa2" or die $!;
    while(<sam_fh>){
        next if($_ =~ /^@/);
        my @a = split(/\t/, $_);
        # remove those that are not primary alignment
        next if($a[1] & 0x0800);
        # convert from reverse to forward strand
        if($a[1] & 0x0010){
            $a[9] = seq::revcom($a[9]);
            $a[10] = scalar reverse($a[10]);
        }
        if(!defined $h->{$a[0]}){
            $h->{$a[0]} = join("\t", @a);
            next;
        }
        # print only in pair
        if($a[1] & 0x0040){
            # this is the first, save to fa1
            print fa_out1 join("\n", ">" . join(" ", $a[0], join(",", @a[1 .. 8])), $a[9])."\n";
            my @b = split(/\t/, $h->{$a[0]});
            print fa_out2 join("\n", ">" . join(" ", $b[0], join(",", @b[1 .. 8])), $b[9])."\n";
        }
        elsif($a[1] & 0x0080){
            # this it the second
            print fa_out2 join("\n", ">" . join(" ", $a[0], join(",", @a[1 .. 8])), $a[9])."\n";
            my @b = split(/\t/, $h->{$a[0]});
            print fa_out1 join("\n", ">" . join(" ", $b[0], join(",", @b[1 .. 8])), $b[9])."\n";
        }
    }
    close fa_out1;
    close fa_out2;
    close sam_fh;
}

sub sam2fa_pair{
    # convert a sam to a pair of fastq files, leave the singles not printed. fa name would be fa_1.fastq, fa_2.fastq
    my ($sam) = @_;
    my ($fa) = ($sam =~ /^(.+)\.sam/);
    my $fa1 = $fa. "_1.fa";
    my $fa2 = $fa. "_2.fa";
    my $h;
    open sam_fh, "<$sam" or die $!;
    open fa_out1, ">$fa1" or die $!;
    open fa_out2, ">$fa2" or die $!;
    while(<sam_fh>){
        next if($_ =~ /^@/);
        my @a = split(/\t/, $_);
        # remove those that are not primary alignment
        next if($a[1] & 0x0800);
        # convert from reverse to forward strand
        if($a[1] & 0x0010){
            $a[9] = seq::revcom($a[9]);
            $a[10] = scalar reverse($a[10]);
        }
        if(!defined $h->{$a[0]}){
            $h->{$a[0]} = join("\t", @a);
            next;
        }
        # print only in pair
        if($a[1] & 0x0040){
            # this is the first, save to fa1
            print fa_out1 join("\n", ">" . $a[0], $a[9])."\n";
            my @b = split(/\t/, $h->{$a[0]});
            delete $h->{$a[0]};
            print fa_out2 join("\n", ">" . $b[0], $b[9])."\n";
        }
        elsif($a[1] & 0x0080){
            # this it the second
            print fa_out2 join("\n", ">" . $a[0], $a[9])."\n";
            my @b = split(/\t/, $h->{$a[0]});
            delete $h->{$a[0]};
            print fa_out1 join("\n", ">" . $b[0], $b[9])."\n";
        }
    }
    close fa_out1;
    close fa_out2;
    close sam_fh;
}

# read sam and store to a hash
sub read_sam{
    my ($sam, $mode) = @_;
    my ($header_h, $sam_h);
    open sam_fh, "<$sam" or die $!;
    while(<sam_fh>){
        if($_ =~ /^\@/){
            if($_ =~ /^\@SQ/){
                if($_ =~ /SN:(\S+).+LN:(\S+)/){
                    $header_h->{SQ}->{$1} = $2;
                }
            }
        }
        else{
            my @a = split(/\t/, $_);
            my ($readname, $flag, $chr, $pos, $mapq, $cigar, $mchr, $mpos, $insertsize, $seq, $qual, ) = @a;
            my $h;
            $h->{chr} = $chr;
            $h->{pos} = $pos;
            $h->{flag} = $flag;
            $h->{mapq} = $mapq;
            $h->{mchr} = $mchr;
            $h->{mpos} = $mpos;
            $h->{insertsize} = $insertsize;
            $h->{seq} = $seq;
            $h->{qual} = $qual;
            $h->{readname} = $readname;
            $h->{cigar} = $cigar;
            if($mode =~ /by_readname/){
                $sam_h->{$readname}->{$cigar} = $h;
            }
            elsif($mode =~ /by_position/){
                $sam_h->{$chr}->{$pos}->{$readname} = $h;
            }
        }

    }
    close sam_fh;

    return ($header_h, $sam_h);
}

# check if any alignment in the sam
sub check_sam{
    my ($sam) = @_;
    my ($header_h, $sam_h) = &read_sam($sam, "by_position");
    my @chrs = keys %$sam_h;
    if(scalar(@chrs) == 1 && $chrs[0] eq "*"){
        return 0;
    }
    return 1;
}

# inspect assembly quality
sub asm_qual{
    # report the N50, number of contigs, number of insertion with size > a number, that of deletion, coverage of the genome
    my ($sam, $ins_size, $del_size) = @_;
    # for calculating genome coverage
    my $h;
    # stat of all ctgs
    my $g;
    my @seq_lens = (0);
    my $n = 0;
    my $ins_num = 0; 
    my $del_num = 0;
    my $ref_len;
    #my ($n, $ins_num, $del_num) = (0, 0, 0);
    open SAM, "<$sam" or die $!;
    while(<SAM>){
        chomp;
        # read header to get the total length of the genome
        if($_ =~ /^\@SQ/){
            my ($ctg_name) = ($_ =~ /SN:(\S+)/);
            my ($ctg_len) = ($_ =~ /LN:(\d+)/);
            $g->{$ctg_name} = $ctg_len;
        }
        else{
            my @a = split(/\t/, $_);
            my $ctg_name = $a[2];
            my $seq = $a[9];
            my $seq_len = length($seq);
            push @seq_lens, $seq_len;
            my $cigar = $a[5];
            $ins_num += cigar::check_cigar($cigar, "INS", $ins_size) if($a[3] > 2100 && $a[3] < 2900);
            $del_num += cigar::check_cigar($cigar, "DEL", $del_size) if($a[3] > 2100 && $a[3] < 2900);
            $ref_len = cigar::check_cigar($cigar, "ref_len");
            # cover the reference if matched
            foreach my $i ($a[3] .. $a[3] + $ref_len){
                $h->{$ctg_name}->{$i} = 1;
            }
            $n ++;
        }
    }
    close SAM;
    # check if the repetitive region is covered
    my $tag;
    foreach my $ctg_name (keys %$h){
        foreach my $i (2100 .. 2900){
            if(!defined $h->{$ctg_name}->{$i}){
                $tag->{$ctg_name} = 1; # not covered
                last;
            }
        }
        if(!defined $tag->{$ctg_name}){
            $tag->{$ctg_name} = 0;
        }
    }
    my $N50 = asm::calculate_N50(\@seq_lens);
    # calculate coverage of reference
    print "\n";
    #$DB::single = 1;
    print join("\t", $N50, $ins_num, $del_num);
    foreach my $ctg (keys %$g){
        print "\t" . join("\t", $ctg . ":" . scalar(keys %{$h->{$ctg}})/$g->{$ctg}) . "\t" . join("\t", $ctg . ":" . $tag->{$ctg});
    }
    print "\n";
}

# To deal with fragmented alignment of contigs, if an INDEL (size > 10bp) is bounded by a matching with size at least sz_M, not necessarily next to it, chop the sequence of the contig and the reference, and run pairwise alignment (blasrn with alignOpen 30 deletion 30 insertion 30), and then infer the breakpoint directly from the alignment result. 
# sz_t: first round of checking, if INDEL with size greater than this
# sz_M: minimum matching anchor size, also used in the second round of checking, pairHMM
# ref: reference file
# min_clip: max clip size this program will tolerate to analyze this sam line
# length_diff_t: in second round of checking, the size of the indel inside the two large matching could have maximum difference from that stays inside the two matching
# max_match: the max matching the program will get in pairwise alignment. If greater than this, it will cut and get those within the max_match size. 
# 01142016 Add Illumina_covered_regions, an array whose elements are regions of the contig that are covered by Illumina reads, with format start-end
sub infer_INDEL_PB_wRealign{
    my ($sam, $min_clip, $sz_t, $sz_M, $ref, $length_diff_t, $max_match, $IL_covered_regions) = @_;
    # use pairHMM for alignment
    my $mode = 2;
    my @bs;
    open SAM, "<$sam" or die $!;
    while(<SAM>){
        next if($_ =~ /^@/);
        my @a = split(/\t/, $_);
        my ($ctg_name, $tag, $chr, $pos, $qual, $cigar) = @a[0 .. 5];
        next if($cigar eq "*");

        my ($head_clip, $tail_clip) = (0, 0);
        if($cigar =~ /^(\d+)[SH]/){
            if($1 > $min_clip){
                $head_clip = $1;
                # clipped at 3' to the contig if reverse; 5' if forward
            }
        }

        if($cigar =~ /(\d+)[SH]$/){
            if($1 > $min_clip){
                $tail_clip = $1;
            }
        }

        next if($head_clip != 0 || $tail_clip != 0);
        my ($ref_pos, $read_pos, $alpha, $num) = cigar::get_pos_from_cigar_input($_);
        # read each line of these
        my $t_num = scalar(@$ref_pos);
        my $previous = 0;
        my $previous_line = "";
        my $this_line;
        for(my $i = 0; $i < $t_num; $i ++){
            if($alpha->[$i] eq "M" && $num->[$i] > $sz_M){
                $this_line = join("\t", $ref_pos->[$i], $read_pos->[$i], $alpha->[$i], $num->[$i]);
                # conclude previous
                if($previous == 1){
                    my @tmp = split(/\t/, $previous_line);
                    if(&check_overlap($IL_covered_regions, "$tmp[1]:$read_pos->[$i]")){ 
                        # mode = 1: blasr with high indel penalty
                        # mode = 2: pairHMM

                        my ($sam, $r) = &realign($previous_line, $this_line, $_, $ref, $max_match, $mode);
                        #my $b = &get_combinations_indel($sam, $sz_M, $sz_t);

                        my $b = &analyze_simple_sam($sam, $sz_M, $length_diff_t, $r);
                        push @bs, $b if($b ne "");
                    }
                }
                # start a new session anyway
                $previous_line = $this_line;
                $previous = 0;
            }
            elsif($alpha->[$i] =~ /[ID]/ && $num->[$i] > $sz_t && $previous_line ne ""){
                $previous = 1;
            }
        }
    } # end of sam
    return \@bs;
}
=cut
# given a sam with slightly fragmented alignment with two ends anchored with long matches, construct underlying allele for all combinations, return the alleles in a hash, with keys the putative indel, the value the underlying allele
sub get_combinations_indel{
    my ($sam, $sz_M, $sz_t) = @_;
    # analyze sam cigar, and the sam does not allow hard clip
    my @as = split(/\n/, $sam); 
    my $b = "";
    foreach my $a (@as){
        if($a =~ /^@/){
            next;
        }
        # only analyze the first line
        my @items = split(/\t/, $a);
        my $read_name = $items[0];
        my $ref_name = $items[2];
        my ($chr, $start, $end) = ($ref_name =~ /^(\S+)\:(\d+)\-(\d+)$/);
        my ($ctg_name, $ctg_start, $ctg_end) = ($read_name =~ /^(\S+)\:(\d+)\-(\d+)\/\d+_\d+$/);
        my $cigar = $items[5];
        my ($alpha, $num) = cigar::analyze_cigar($cigar);
        my $total = scalar(@$num);
        my $max = 0;
        my $max_alpha = "";
        my $max_ref_pos;
        my $max_ctg_pos;
        my $ref_middle_len = 0;
        my $ctg_middle_len = 0;
        if($num->[0] > $sz_M && $num->[$total-1] > $sz_M && $alpha->[0] eq "M" && $alpha->[$total-1] eq "M"){
            $ref_middle_len += $num->[0];
            $ctg_middle_len += $num->[0];
            for(my $i = 1; $i < $total - 1; $i ++){
                # check all indels greater than sz_t
                if($num->[$i] > $sz_t){
                    if($alpha->[$i] eq "D"){
                        push @poss, $ctg_middle_len . "-" . ($ctg_middle_len + 1) . ";" . $read_name . ":" . $ref_middle_len . "-" . ($ref_middle_len + $num->[$i]);
                    }
                    elsif($alpha->[$i] eq "I"){
                        push @poss, $ctg_middle_len . "-" . ($ctg_middle_len + $num->[$i]) . ";" . $read_name . ":" . $ref_middle_len . "-" . ($ref_middle_len + 1);
                    }
                }
                if($alpha->[$i] =~ /[MD]/){
                    $ref_middle_len += $num->[$i];
                }
                if($alpha->[$i] =~ /[MI]/){
                    $ctg_middle_len += $num->[$i];
                }
            }
        } # end of checking this sam line's head and tail
    } # end of reading sam lines
    return @$b;
}



=cut

# check if overlap between IL covered regions (in an array with each element start-end) and the str (start-end) on a contig. Return 1 if yes, 0 if no, or if no IL regions given, return 1.
sub check_overlap{
    my ($IL_covered_regions, $str) = @_;
    if(scalar(@$IL_covered_regions) == 0){
        return 1;
    }
    my ($start, $end) = split(/-/, $str);
    foreach my $i (0 .. scalar(@$IL_covered_regions) - 1){
        my ($s, $e) = split(/-/, $IL_covered_regions->[$i]);
        if($s < $end && $s > $start || $start < $e && $start > $s){
            return 1;
        }
    }
    return 0;
}







# given two lines of matching anchors with positions of both reference and contig, the sam line and reference, do pairwise alignment around this area and return the new sam line
sub realign{
    my ($previous, $this, $sam, $ref, $max_match, $mode) = @_;
    my ($ref_pos1, $ctg_pos1, $tmp1, $len1) = split(/\s+/, $previous);
    my ($ref_pos2, $ctg_pos2, $tmps, $len2) = split(/\s+/, $this);
    # forward
    my $r = 0;
    if($ctg_pos1 > $ctg_pos2){
        $r = 1;
    }
    # to limit the running time by avoiding too long matching 
    if($len1 > $max_match){
        $ref_pos1 += $len1 - $max_match;
        if($r == 0){
            $ctg_pos1 += $len1 - $max_match;
        }
        else{
            $ctg_pos1 -= $len1 - $max_match;
        }
    }
    $ref_pos2 += $len2;
    if($r == 0){
        $ctg_pos2 += $len2;
    }
    else{
        $ctg_pos2 -= $len2;
    }
    if($len2 > $max_match){
        $ref_pos2 -= $len2 - $max_match;
        if($r == 0){
            $ctg_pos2 -= $len2 - $max_match;
        }
        else{
            $ctg_pos2 += $len2 - $max_match;
        }
    }
    # get info from sam
    my @a = split(/\t/, $sam);
    my $chr = $a[2];
    my $str = $a[9];
    my $ref_tmp = "ref_$chr.$ref_pos1.$ref_pos2.fasta";
    my $ctg_tmp = "read_$ctg_pos1.$ctg_pos2.fasta";
    # extract ref
    #$DB::single = 1;
    `samtools faidx $ref $chr:$ref_pos1-$ref_pos2 > $ref_tmp`;
    # extract contig
    my $substr;
    if($r == 0){
        $substr = substr($str, $ctg_pos1, $ctg_pos2 - $ctg_pos1);
        `echo ">$a[0]:$ctg_pos1-$ctg_pos2\n$substr\n" > $ctg_tmp`;
    }
    else{
        $substr = reverse(substr($str, $ctg_pos2, $ctg_pos1 - $ctg_pos2));
        `echo ">$a[0]:$ctg_pos2-$ctg_pos1\n$substr\n" > $ctg_tmp`;
    }
    my $ret = "";
    if($mode == 1){
        $ret = `~/blasrn $ctg_tmp $ref_tmp -affineOpen 30 -deletion 30 -insertion 30 -affineExtend 0 -sam`;
    }
    elsif($mode == 2){
        $ret = `perl ~/pairHMM/globalHMM_freeze/testHMM.pl -V $ref_tmp $ctg_tmp ~/pairHMM/globalHMM_freeze/hmm_trans.txt`;
    }
    print STDERR "#perl ~/pairHMM/globalHMM_freeze/testHMM.pl $ref_tmp $ctg_tmp ~/pairHMM/globalHMM_freeze/hmm_trans.txt";
    print STDERR "#~/blasrn $ctg_tmp $ref_tmp -affineOpen 30 -deletion 30 -insertion 30 -affineExtend 0 -sam";
    #`rm $ref_tmp`;
    #`rm $ctg_tmp`;
    return ($ret, $r);
}

# local bwa to resolve complex events such as DEL + INX
sub bwa_local{
    my ($ctg_fa, $ref_fa, $ctg_str, $ref_str, $M_threshold, $indel_threshold, $prefix, $bwa_M_threshold) = @_;

    # changed 07312016
    #my $pairHMM = "pairHMM";
    #my $pairHMM_trans_matrix = "/scratch/bcb/xfan3/p/pairHMM/hmm_trans.ctg.tab.txt";
    #my $pairHMM_trans_matrix = "./pairHMM/hmm_trans.ctg.tab.txt";

    my $ctg_file_name = &analyze_str($ctg_str);
    my $ref_file_name = &analyze_str($ref_str);

    if(! -e "$ctg_fa.fai"){
        `samtools faidx $ctg_fa`;
    }
    `samtools faidx $ctg_fa $ctg_str > $prefix/$ctg_file_name.fa`;
    `samtools faidx $ref_fa $ref_str > $prefix/$ref_file_name.fa`;
    #if($revcom == 1){
    # changed 08282016
    #seq::revcom_file("$prefix/$ref_file_name.fa");
    # changed 07312016
    #my ($rev_fa) = seq::revcom_file("$prefix/$ref_file_name.fa");
    #open fh_rev_out, ">$prefix/$ref_file_name.fa";
    #print fh_rev_out $rev_fa;
    #close fh_rev_out;
    #`~/u/revcom_file.pl $prefix/$ref_file_name.fa > $prefix/$ref_file_name.rev.fa`;
    #`mv $prefix/$ref_file_name.rev.fa $prefix/$ref_file_name.fa`;
    #}
    # here use modification error 0.1 as a constant
    #`$pairHMM -m $modification_error $prefix/$ctg_file_name.fa $prefix/$ref_file_name.fa $pairHMM_trans_matrix > $prefix/$ctg_file_name.aln.csv`;
    #`perl $pairHMM $prefix/$ctg_file_name.fa $prefix/$ref_file_name.fa $pairHMM_trans_matrix $modification_error > $prefix/$ctg_file_name.aln.csv`;
    `bwa index $prefix/$ref_file_name.fa`;
    `bwa mem $prefix/$ref_file_name.fa $prefix/$ctg_file_name.fa > $prefix/$ctg_file_name.sam`; 
# turn 01s to a cigar string
    #my $cigar = hmm::analyze_csv("$prefix/$ctg_file_name.aln.csv");
    #my ($cigar, $matches) = hmm::analyze_csv_w_matches("$prefix/$ctg_file_name.aln.csv");
#print $cigar . "\n";
    #$DB::single = 1;


# make a sam file out of cigar
    #&make_sam_from_pairHMMcigar($ctg_file_name, $ref_file_name, $revcom, $cigar, "$prefix/$ctg_file_name.sam", $ctg_name);
    #my $b = &analyze_simple_sam("$prefix/$ctg_file_name.sam", $M_threshold, $indel_threshold);
    $DB::single = 1;
    my $b = &analyze_complex_sam("$prefix/$ctg_file_name.sam", $bwa_M_threshold, $indel_threshold, $ctg_fa);
    # to ignore this contig if all of the followings are true 
    # there are >= 2 matches, whose left and right neighbors are matches, and either the match length < 100 or the mismatch num > 3
    # First, check if large gaps are separated by matches > 3
    # If yes, check if these gaps can be merged (close enough). 
    #   If yes, output the merged result (dump the individual small deletions). Otherwise, dump all calls from this contig.
    # If no, check if any gaps can be merged (close enough). Outpu both the merged and the individual results for further check of IL anyway.  
    # Xian @ 10032016: now since the DEL + INS is possible, only look for this type of complex SVs. 
    #if(&check_gap_w_match($cigar, $matches, $M_threshold, $indel_threshold) == 1){
    #    print STDERR "#Multiple gaps with erroneous gaps in between: $ctg_name\n";
    #    $b = &take_microinsertion($b, $dist);
    #    if(@$b >= 1){
    #        print STDERR "#Being able to merge the gaps into large ones: \n";
    #        print STDERR "#" . join("\n", @$b) . "\n";
    #    }
    #    else{
    #        print STDERR "#After merging, nothing left. \n";
    #    }
    # now instead of dumping this call, check if can merge them into large deletions and send it for IL confirmation
    #}
    #else{
    #my $c = &take_microinsertion($b, $dist);
    #if(scalar(@$c) >= 1){
    #    print STDERR "#Being able to merge the gaps into large ones: \n";
    #    print STDERR "#" . join("\n", @$c) . "\n";
    #}
    #push @$b, @$c;
    #} 
    # microinsertion in deletion leads to fragmented deletions, separated by small matches, which should have been microinsertions on the ctg and deletion on the reference. For these cases, propose a longer deletion, and leave it to be confirmed by IL. If confirmed by >= ILs for smaller ones, take the large one instead. 
    #$DB::single = 1;
    return $b;
}


# allow complex events such as DEL + INS, or double DELs with spacer
sub analyze_complex_sam{ 
    my ($sam, $M_threshold, $indel_threshold, $ctg_file) = @_;
    $DB::single = 1;
    my $gap = $indel_threshold;
    my $min_SV = $indel_threshold;
    open fh_, "<$sam" or die $!;
    my $ctg_name;
    my @b;
    my $h;
    my $chr;
    while(<fh_>){
        next if($_ =~ /^@/);
        my @x = split(/\t/, $_);
        next if($x[2] eq "*");
        my ($ctg, $ctg_start, $ctg_end) = ($x[0] =~ /^(.+):(\d+)-(\d+)/);
        my ($ref, $ref_start, $ref_end) = ($x[2] =~ /^(.+):(\d+)-(\d+)/);
        $chr = $ref;
        $ctg_name = $ctg;
        my $pos1 = $x[3] + $ref_start;
        my $cigar = $x[5];
        my $length = 0;
        my $num;
        my $tag;
        while($cigar ne ""){
            ($num, $tag, $cigar) = ($cigar =~ /^(\d+)([SHMID])(.*)$/);
            $length += $num if($tag =~ /[MD]/);
        }
        $cigar = $x[5];
        my $pos2 = $pos1 + $length;
        # get orientation (-1: forward; 1: backforward);
        my $r = -1;
        my $r_s = $pos1;
        my $r_e = $pos2;
        if($x[1] & 0x0010){
            $r = 1;
            # do not need to reverse the two since they are already in the increasing order
            #$r_s = $pos2;
            #$r_e = $pos1;
        }
        # the start and end position on ctg that have alignment
        my $ctg_s;
        my $ctg_e;
        # get the actual alignment region
        if($cigar =~ /^(\d+)[SH]/){
            if($r == -1){
                $ctg_s = $ctg_start + $1;
            }
            else{
                $ctg_e = $ctg_end - $1;
            }
        }
        if($cigar =~ /(\d+)[SH]$/){
            if($r == -1){
                $ctg_e = $ctg_end - $1;
            }
            else{
                $ctg_s = $ctg_start + $1;
            }
        }
        # note that reverse need to be reversed again, head -> tail, tail -> head
        #if($r == 1){
        #    my $tmp = $ctg_e;
        #    $ctg_e = $ctg_s;
        #    $ctg_s = $tmp;
        #}
        # only save the clipped end
        # 10032016 (scaffold.pm does not have this change): for the easiness of determining deletion, add one more key, and remove the readname since they are in one sam from bwa run, and all are wiht the same readname
        # determine cigar string pattern
        # if ^SM$ => 1
        # if ^MS$ => 2
        # if ^SMS$ => 3
        #my $cigar_pattern = &determine_cigar_pattern($cigar);
        my $cigar_pattern;
        if($cigar =~ /^\d+[SH].+M$/){
            if($cigar =~ /(\d+)M$/ && $1 > $M_threshold){
                $cigar_pattern = 1;
                if($r == -1){
                    # 5' end of contig was clipped, 5' end of reference was clipped
                    $h->{$cigar_pattern}->{$ctg_s} = "$chr:$pos1:5";
                }
                else{
                    # 3' end of contig was clipped, 5' end of reference was clipped
                    $h->{$cigar_pattern}->{-$ctg_e} = "$chr:$pos1:5";
                }
            }
        }
        elsif($cigar =~ /^(\d+)M.*[SH]$/ && $1 > $M_threshold){
            $cigar_pattern = 2;
            if($r == -1){
                # 3' end of contig was clipped, 3' end of reference was clipped
                $h->{$cigar_pattern}->{-$ctg_e} = "$chr:$pos2:3";
            }
            else{
                # 5' end of contig was clipped, 3' end of reference was clipped
                $h->{$cigar_pattern}->{$ctg_s} = "$chr:$pos2:3";
            }
        }
        elsif($cigar =~ /^\d+[SH].+\d+[SH]$/){
            # for those that align on the middle, only report when the aligned length is > 10
            next if(abs($ctg_s - $ctg_e) < $gap);
            $cigar_pattern = 3;
            # contig position in the order of 5' -> 3' on the reference
            if($r == -1){
                $h->{$cigar_pattern}->{"$ctg_s:$ctg_e"} = "$chr:$pos1:$pos2";
            }
            else{
                # key always in the format of left_breakpoint:right_breakpoint
                $h->{$cigar_pattern}->{"$ctg_e:$ctg_s"} = "$chr:$pos1:$pos2";
            }
        }

    } # end of reading the alignment
    close fh_;

    # begin analyzing the alignment
    # for a deletion, only look at cigar pattern 1 and 2, and see if position difference on the reference is at least 10bp larger than contig breakpoint difference, and then determine if this is a complex event (contig breakpoint difference > 10bp => INS). If yes, further look at cigar pattern 3, and see if any inserted sequence align to the reference, report it.  
    # TODO: implement the above (10032016)
    # done
    if(defined $h->{1} && defined $h->{2}){
        my @ctgs_1 = keys %{$h->{1}};
        my @ctgs_2 = keys %{$h->{2}};
        foreach my $ctg1_p (@ctgs_1){
            foreach my $ctg2_p (@ctgs_2){
                # check if deletion is formed directly
                # gap is the distance between the two breakpoints on the ctg, including overlapping gap
                my $ctg2_p_abs = abs($ctg2_p);
                my $ctg1_p_abs = abs($ctg1_p);
                my $dist = $ctg1_p_abs - $ctg2_p_abs;
                # augment sequence might be micro_ins, or microhology
                my $augment_seq = "";
                if($ctg1_p_abs < $ctg2_p_abs){
                    my @augment_seqs = split(/\n/, `samtools faidx $ctg_file $ctg_name:$ctg1_p_abs-$ctg2_p_abs`);
                    $augment_seq = join("", @augment_seqs[1 .. $#augment_seqs]);
                }
                elsif($ctg1_p_abs > $ctg2_p_abs){
                    my @augment_seqs = split(/\n/, `samtools faidx $ctg_file $ctg_name:$ctg2_p_abs-$ctg1_p_abs`);
                    $augment_seq = join("", @augment_seqs[1 .. $#augment_seqs]);
                }

                my ($chr1, $pos1, $ori1) = split(/:/, $h->{1}->{$ctg1_p});
                my ($chr2, $pos2, $ori2) = split(/:/, $h->{2}->{$ctg2_p});

                # require orientation here
                # inv can be modified here
                # they should be corresponding to 5' clip and 3' clip on the contig
                if($ctg2_p * $ctg1_p > 0){
                    next;
                }
                if($ctg1_p < 0 && length($augment_seq) > 0){
                    $augment_seq = seq::revcom($augment_seq);
                }
                if($chr1 eq $chr2 && $pos1 - $pos2 - abs($dist) > $min_SV){
                    # TODO now didn't consier microhomology > gap (need to take this)
                    if(abs($dist) <= $gap){
                        # check if the corresponding position on the reference are at least 10bp apart
                        # min_SV is the minimum SV we start to report
                        # report this deletion
                        my $b_ = join("\t", $chr1, $pos2, ($pos1 - $pos2), "D", $ctg_name, $ctg2_p_abs, $ctg1_p_abs);
                        # check if microhomology
                        if($dist < 0){
                            $b_ = join("\t", $b_, "homology:$augment_seq:$ctg1_p_abs:$ctg2_p_abs");
                        }
                        elsif($dist == 0){
                            $b_ = join("\t", $b_, "endj");
                        }
                        else{
                            $b_ = join("\t", $b_, "microi:$augment_seq:$ctg2_p_abs:$ctg1_p_abs");
                        }

                        push @b, $b_;
                    }
                    else{
                        # large insertion within the deletion
                        # check if any alignment for pattern 3
                        if(defined $h->{3}){
                            # check the positions on the contig, see if they are within the inserted region
                            my $b_ = join("\t", $chr1, $pos2, ($pos1 - $pos2), "D", $ctg_name, $ctg2_p_abs, $ctg1_p_abs, "largei:$augment_seq");
                            my @tmps;
                            foreach my $key_3 (keys %{$h->{3}}){
                                # match position on the contig
                                my ($ctg_p2, $ctg_p1) = split(/:/, $key_3);
                                # make sure they are mostly within the inserted region
                                # that they don't overlap with the flanking region too much
                                # check orientation by ctg1_p_abs and ctg2_p_abs
                                if($ctg1_p_abs > $ctg2_p_abs && $ctg_p1 - $ctg1_p_abs < $gap && $ctg_p2 - $ctg2_p_abs > -$gap || $ctg1_p_abs < $ctg2_p_abs && $ctg_p1 - $ctg1_p_abs < $gap && $ctg_p2 - $ctg2_p_abs > -$gap){
                                    my @aug_seq_tmps;
                                    # check orientation
                                    if($ctg_p2 > $ctg_p1){
                                        @aug_seq_tmps = split(/\n/, `samtools faidx $ctg_file $ctg_name:$ctg_p1-$ctg_p2`);
                                    }
                                    else{
                                        @aug_seq_tmps = split(/\n/, `samtools faidx $ctg_file $ctg_name:$ctg_p2-$ctg_p1`);
                                    }

                                    my $aug_seq_tmp = join("", @aug_seq_tmps[1 .. $#aug_seq_tmps]);
                                    if($ctg_p2 > $ctg_p1){
                                        $aug_seq_tmp = seq::revcom($aug_seq_tmp);
                                    }
                                    my $tmp = "$key_3:$h->{3}->{$key_3}:$aug_seq_tmp";

                                    push @tmps, $tmp;
                                }
                            }
                            $b_ .= ";" . join(";", @tmps);
                            push @b, $b_;
                        }
                        else{
                            # novel insertion, at least not overlappign with the flanking regino
                            my $b_ = join("\t", $chr1, $pos2, ($pos1 - $pos2), "D", $ctg_name, $ctg2_p_abs, $ctg1_p_abs, "largei:$augment_seq");
                            push @b, $b_;
                        }
                    } # end of if simple deletion or deletion + insertion
                } # end of if a deletion exists
            } # end of foreach {2]
        } # end of foreach {1}
    } # end of if defined {1} and {2}
    return \@b;
}




# memroy efficieint for pairHMM
sub refind_aln_hmm_memEffi{
    my ($ctg_fa, $ref_fa, $ctg_str, $ref_str, $revcom, $skip_aln, $M_threshold, $indel_threshold, $ctg_name, $prefix, $modification_error, $dist, $pairHMM_trans_matrix) = @_;

    # changed 07312016
    my $pairHMM = "pairHMM";
    #my $pairHMM_trans_matrix = "/scratch/bcb/xfan3/p/pairHMM/hmm_trans.ctg.tab.txt";
    #my $pairHMM_trans_matrix = "./pairHMM/hmm_trans.ctg.tab.txt";

    my $ctg_file_name = &analyze_str($ctg_str);
    my $ref_file_name = &analyze_str($ref_str);

    if($skip_aln == 0){
        if(! -e "$ctg_fa.fai"){
            `samtools faidx $ctg_fa`;
        }
        `samtools faidx $ctg_fa $ctg_str > $prefix/$ctg_file_name.fa`;
        `samtools faidx $ref_fa $ref_str > $prefix/$ref_file_name.fa`;
        if($revcom == 1){
            # changed 08282016
            seq::revcom_file("$prefix/$ref_file_name.fa");
            # changed 07312016
            #my ($rev_fa) = seq::revcom_file("$prefix/$ref_file_name.fa");
            #open fh_rev_out, ">$prefix/$ref_file_name.fa";
            #print fh_rev_out $rev_fa;
            #close fh_rev_out;
            #`~/u/revcom_file.pl $prefix/$ref_file_name.fa > $prefix/$ref_file_name.rev.fa`;
            #`mv $prefix/$ref_file_name.rev.fa $prefix/$ref_file_name.fa`;
        }
        # here use modification error 0.1 as a constant
        `$pairHMM -m $modification_error $prefix/$ctg_file_name.fa $prefix/$ref_file_name.fa $pairHMM_trans_matrix > $prefix/$ctg_file_name.aln.csv`;
        #`perl $pairHMM $prefix/$ctg_file_name.fa $prefix/$ref_file_name.fa $pairHMM_trans_matrix $modification_error > $prefix/$ctg_file_name.aln.csv`;
    }
# turn 01s to a cigar string
    #my $cigar = hmm::analyze_csv("$prefix/$ctg_file_name.aln.csv");
    my ($cigar, $matches) = hmm::analyze_csv_w_matches("$prefix/$ctg_file_name.aln.csv");
#print $cigar . "\n";
    #$DB::single = 1;


# make a sam file out of cigar
    &make_sam_from_pairHMMcigar($ctg_file_name, $ref_file_name, $revcom, $cigar, "$prefix/$ctg_file_name.sam", $ctg_name);
    my $b = &analyze_simple_sam("$prefix/$ctg_file_name.sam", $M_threshold, $indel_threshold);
    # to ignore this contig if all of the followings are true 
    # there are >= 2 matches, whose left and right neighbors are matches, and either the match length < 100 or the mismatch num > 3
    # First, check if large gaps are separated by matches > 3
    # If yes, check if these gaps can be merged (close enough). 
    #   If yes, output the merged result (dump the individual small deletions). Otherwise, dump all calls from this contig.
    # If no, check if any gaps can be merged (close enough). Outpu both the merged and the individual results for further check of IL anyway.  
    if(&check_gap_w_match($cigar, $matches, $M_threshold, $indel_threshold) == 1){
        print STDERR "#Multiple gaps with erroneous gaps in between: $ctg_name\n";
        $b = &take_microinsertion($b, $dist);
        if(@$b >= 1){
            print STDERR "#Being able to merge the gaps into large ones: \n";
            print STDERR "#" . join("\n", @$b) . "\n";
        }
        else{
            print STDERR "#After merging, nothing left. \n";
        }
        # now instead of dumping this call, check if can merge them into large deletions and send it for IL confirmation
    }
    else{
        my $c = &take_microinsertion($b, $dist);
        if(scalar(@$c) >= 1){
            print STDERR "#Being able to merge the gaps into large ones: \n";
            print STDERR "#" . join("\n", @$c) . "\n";
        }
        push @$b, @$c;
    }
    # microinsertion in deletion leads to fragmented deletions, separated by small matches, which should have been microinsertions on the ctg and deletion on the reference. For these cases, propose a longer deletion, and leave it to be confirmed by IL. If confirmed by >= ILs for smaller ones, take the large one instead. 
    #$DB::single = 1;
    return $b;
}

# to merge deletions separated by match (which might be caused by micro-insertions) into large deletions
sub take_microinsertion{
    my ($bb, $dist) = @_; 
    my @aa = @$bb;
    my @ret;
    my $h;
    foreach my $a_ (@aa){
        my @b_ = split(/\t/, $a_);
        if($b_[3] eq "D"){
            $h->{$b_[1]} = $a_;
        }
    }
    my $pre = -1;
    my $pre_del_len = 0;
    my $ref_start = -1;
    my $ctg_start = -1;
    my $ctg_end = -1;
    my $length = -1;
    my $chr = "NA";
    my $str = "NA";
    my $ctg = "NA";
    foreach my $p (sort {$a <=> $b} keys %$h){
        if($pre == -1){
            my @x = split(/\t/, $h->{$p});
            $ref_start = $x[1];
            $ctg_start = $x[5];
            $length = $x[2];
            $chr = $x[0];
            $ctg = $x[4];
        }
        elsif($p - ($pre + $pre_del_len) < $dist){
            my @x = split(/\t/, $h->{$p});
            $ctg_end = $x[5];
            $length += $x[2] + ($p - ($pre + $pre_del_len)); 
            $str = join("\t", $chr, $ref_start, $length, "D", $ctg, $ctg_start, $ctg_end); 
        }
        else{
            # conclude
            my @x = split(/\t/, $h->{$p});
            if($str ne "NA"){
                push @ret,$str;
                $str = "NA";
            }
            $ref_start = $x[1];
            $ctg_start = $x[5];
            $length = $x[2];
            $chr = $x[0];
            $ctg = $x[4];
        }
        $pre = $p;
        $pre_del_len = $length;
    }
    if($str ne "NA"){
        push @ret, $str;
    }
    return \@ret;
}


sub check_gap_w_match{
    my ($cigar, $match_str, $M_threshold, $indel_threshold) = @_;
    my ($tag, $num) = cigar::analyze_cigar($cigar);
    my $tag1 = 0;
    my @matches = split(/:/, $match_str);
    for(my $i = 0; $i < scalar(@$tag) - 2; $i ++){
        if($tag->[$i] =~ /[ID]/ && $tag->[$i + 1] eq "M" && $tag->[$i + 2] =~ /[ID]/){
            next if($num->[$i] < $indel_threshold || $num->[$i + 2] < $indel_threshold); 
            if($matches[$i + 1] > 3){
                #|| $num->[$i + 1] < $M_threshold){
                $tag1 ++;
            }
        }
    }
    if($tag1 >= 2){
        return 1;
    }
    else{
        return 0;
    }
}



sub refind_aln_hmm{
    my ($ctg_fa, $ref_fa, $ctg_str, $ref_str, $revcom, $skip_aln, $M_threshold, $indel_threshold, $ctg_name, $prefix, $modification_error) = @_;

    my $pairHMM = "/scratch/bcb/xfan3/p/pairHMM/testHMM.pl";
    my $pairHMM_trans_matrix = "/scratch/bcb/xfan3/p/pairHMM/hmm_trans.ctg.txt";

    my $ctg_file_name = &analyze_str($ctg_str);
    my $ref_file_name = &analyze_str($ref_str);

    if($skip_aln == 0){
        `samtools faidx $ctg_fa $ctg_str > $prefix/$ctg_file_name.fa`;
        `samtools faidx $ref_fa $ref_str > $prefix/$ref_file_name.fa`;
        if($revcom == 1){
            `~/u/revcom_file.pl $prefix/$ref_file_name.fa > $prefix/$ref_file_name.rev.fa`;
            `mv $prefix/$ref_file_name.rev.fa $prefix/$ref_file_name.fa`;
        }
        `perl $pairHMM $prefix/$ctg_file_name.fa $prefix/$ref_file_name.fa $pairHMM_trans_matrix $modification_error > $prefix/$ctg_file_name.aln.csv`;
    }
# turn 01s to a cigar string
    my $cigar = hmm::analyze_csv("$prefix/$ctg_file_name.aln.csv");
#print $cigar . "\n";


# make a sam file out of cigar
    &make_sam_from_pairHMMcigar($ctg_file_name, $ref_file_name, $revcom, $cigar, "$prefix/$ctg_file_name.sam", $ctg_name);
    my $b = &analyze_simple_sam("$prefix/$ctg_file_name.sam", $M_threshold, $indel_threshold);
    return $b;
}

sub analyze_str{
    my ($str) = @_;
    my $ret = "NA";
    if($str =~ /(\S+):(\d+)-(\d+)/){
        $ret = join("_", $1, $2, $3);
        return $ret;
    }
    return $ret;
}

# intersect with cluster position
sub infer_INDEL_PB_wPairHMM_memEffi_wClusterPos{
    my ($sam, $min_clip, $sz_t, $sz_M, $sz_m, $ctg_fa, $ref_fa, $prefix, $modification_error, $microins_dist, $cluster_pos, $pairHMM_trans_matrix) = @_;
    my @lines;
    open SAM, "<$sam" or die $!;
    while(<SAM>){
        next if($_ =~ /^@/);
        my @a = split(/\t/, $_);
        my ($ctg_name, $tag, $chr, $pos, $qual, $cigar) = @a[0 .. 5];
        #$DB::single = 1;
        my ($chr1, $chrs, $chre) = ($chr =~ /^(.+)\:(\d+)-(\d+)/);
        my $ctg_name_for_fa = $ctg_name;
        while($ctg_name_for_fa =~ /(.+)\/\d+_\d+/){
            $ctg_name_for_fa = $1;
        }
        $ctg_name = $ctg_name_for_fa;
        next if($cigar eq "*");

        my $r = 0;
        if($tag & 0x0010){
            # reverse
            $r = 1;
        }

        my ($head_clip, $tail_clip) = (0, 0);
        my ($head_clip_t, $tail_clip_t) = (0, 0);
        if($cigar =~ /^(\d+)[SH]/){
            $head_clip_t = $1;
            if($1 > $min_clip){
                $head_clip = $1;
                # clipped at 3' to the contig if reverse; 5' if forward
            }
        }

        if($cigar =~ /(\d+)[SH]$/){
            $tail_clip_t = $1;
            if($1 > $min_clip){
                $tail_clip = $1;
            }
        }

        my $cigar_tmp = $cigar;
        my $ctg_len = 0;
        my $num;
        my $alpha;
        while($cigar_tmp ne ""){
            ($num, $alpha, $cigar_tmp) = ($cigar_tmp =~ /^(\d+)([SMDIH])(.*)$/);
            $ctg_len += $num if($alpha =~ /[SHMI]/);
        }


        # should be both sides do not have large clip
        next if($head_clip != 0 || $tail_clip != 0);
        # add this line when needing to set -q (max clip) very large to allow contig end to be misassembled
        next if($ctg_len - $head_clip - $tail_clip < $head_clip + $tail_clip);
        my $bp_pos = $pos;
        my $pb_s = 0;
        my $pb_e = $pb_s + 1;
        my $line = "";
        my $check_next = 0;
        my $prev = 0;
        my $con_M = 0;
        # a tag to show that the current status is combined with the previous call
        my $combine = 0;

        # find all M's with size > sz_M. For each gap within two consecutive such M's, qualify it for pairHMM if
        # 1. There is a gap > sz_t; or
        # 2. The gap sum of one kind > sz_t, and the Ms inside the region (between two large M's) are < sz_m
        # Qualify the gap directly if it is the only one in the region and it is > sz_t. Do not need pairHMM. 
        # This is to 1) refine breakpoint; 2) avoid missing ones

        my $h_M;
        my $previous = -1;
        my ($ref_p, $ctg_p);
        ($ref_p, $ctg_p, $alpha, $num) = cigar::get_pos_from_cigar_input_consider_rev($_);
        my $total = scalar(@$alpha);
        for(my $i = 0; $i < $total; $i ++){
            if($num->[$i] > $sz_M && $alpha->[$i] eq "M"){
                $h_M->{$i} = 1;
            }
        }
        for(my $i = 0; $i < $total; $i ++){
            if($previous == -1 && defined $h_M->{$i}){
                $previous = $i;
            }
            elsif(defined $h_M->{$i}){
                # a gap formed
                # check if a single one, no need for pairHMM
                if($i - $previous == 2 && $alpha->[$i-1] =~ /[DI]/ && $num->[$i-1] > $sz_t){
                    if($r == 0){
                        $line = join("\t", $chr1, $chrs + $ref_p->[$i-1], $num->[$i-1], $alpha->[$i-1], $ctg_name, $ctg_p->[$i-1], $ctg_p->[$i]);
                    }
                    else{
                        $line = join("\t", $chr1, $chrs + $ref_p->[$i-1], $num->[$i-1], $alpha->[$i-1], $ctg_name, $ctg_p->[$i], $ctg_p->[$i-1]);
                    }

=cut
my $ctg_pos1;
                    my $ctg_pos2;
                    if($r == 0){
                        $ctg_pos1 = $ctg_p->[$i-1];
                        if($alpha->[$i] eq "D"){
                            $ctg_pos2 = $ctg_pos1 + 1; 
                        }
                        elsif($alpha->[$i] eq "I"){
                            $ctg_pos2 = $ctg_pos1 + $num->[$i-1];
                        }
                    }
                    else{
                        if($alpha->[$i] eq "D"){
                            $ctg_pos1 = $ctg_length - $ctg_p->[$i-1];
                            $ctg_pos2 = $ctg_pos1 + 1;
                        }
                        elsif($alpha->[$i] eq "I"){
                            $ctg_pos1 = $ctg_length - $ctg_p-> 
                            $ctg_MI - $num->[$i] + $ctg_start;
                            $ctg_pos2 = $ctg_length - $ctg_MI + $ctg_start;
                        }
                    }
=cut

if(&check_overlap_cluster($line, $cluster_pos) == 1){
    push @lines, $line;
}
$previous = $i;
next;
                }
                elsif($i - $previous == 2 && ($num->[$i-1] <= $sz_t || $alpha->[$i-1] !~ /[DI]/)){
                    $previous = $i;
                    next;
                }
                else{
                    # see whether this region qualifies for pairHMM
                    my $cigar_substr = cigar::get_substr($cigar, $previous, $i);
                    my $select = 0;
                    my $deselect = 0;
                    my $indel_sum = 0;
                    for(my $j = $previous + 1; $j < $i; $j ++){
                        if($alpha->[$j] =~ /[ID]/ && $num->[$j] > $sz_t){
                            $select = 1;
                            last;
                        }
                        elsif($alpha->[$j] =~ /[ID]/){
                            $indel_sum += $num->[$j];
                        }
                        elsif($alpha->[$j] eq "M"){
                            if($num->[$j] > $sz_m){
                                # not select it if all gap < sz_t, sum > sz_t because there is at least one M > sz_m
                                $deselect = 1;
                            }
                        }
                    }
                    if($select == 1 || $indel_sum > $sz_t && $deselect != 1){ 
                        # qualify it for pairHMM
                        my $ref_str;
                        my $ctg_str;
                        $ref_str = "$chr1:" . ($ref_p->[$previous + 1] - $sz_M + $chrs) . "-" . ($ref_p->[$i] + $sz_M + $chrs);
                        if($r == 1){
                            $ctg_str = "$ctg_name_for_fa:" . ($ctg_p->[$i] - $sz_M) . "-" . ($ctg_p->[$previous + 1] + $sz_M);
                        }
                        else{
                            $ctg_str = "$ctg_name_for_fa:" . ($ctg_p->[$previous + 1] - $sz_M) . "-" . ($ctg_p->[$i] + $sz_M);
                        }
                        my $bs = &refind_aln_hmm_memEffi($ctg_fa, $ref_fa, $ctg_str, $ref_str, $r, 0, $sz_M, $sz_t, $ctg_name, $prefix, $modification_error, $microins_dist, $pairHMM_trans_matrix);
                        foreach my $line (@$bs){
                            if(&check_overlap_cluster($line, $cluster_pos) == 1){
                                push @lines, $line;
                            }
                        }
                    }
                    $previous = $i;
                }
            }
        }
    }
    close SAM;
    return \@lines;
}

# check if the pos in line overlap with the cluster_pos
sub check_overlap_cluster{
    my ($line, $h_clu) = @_;
    my @a = split(/\t/, $line);
    foreach (keys %$h_clu){
        my $s = $h_clu->{$_}->{start};
        my $e = $h_clu->{$_}->{end};
        if($a[$#a - 1] <= $e && $a[$#a - 1] >= $s || $a[$#a] <= $e && $a[$#a] >= $s){
            return 1;
        }
    }
    return 0;
}

# instead of pairHMM, this module takes the local reference and loca contig and use BWA to resolve complex SV, mainly DEL + INS
# 10032016 by Xian
sub infer_INDEL_PB_wLocalBWA{
    my ($sam, $min_clip, $sz_t, $sz_M, $sz_m, $ctg_fa, $ref_fa, $prefix, $bwa_M_threshold, $flank) = @_;
    my @lines;
    open SAM, "<$sam" or die $!;
    while(<SAM>){
        next if($_ =~ /^@/);
        my @a = split(/\t/, $_);
        my ($ctg_name, $tag, $chr, $pos, $qual, $cigar) = @a[0 .. 5];
        my $ctg_name_for_fa = $ctg_name;
        while($ctg_name_for_fa =~ /(.+)\/\d+_\d+/){
            $ctg_name_for_fa = $1;
        }
        $ctg_name = $ctg_name_for_fa;
        next if($cigar eq "*");

        my $r = 0;
        if($tag & 0x0010){
            # reverse
            $r = 1;
        }

        my ($head_clip, $tail_clip) = (0, 0);
        my ($head_clip_t, $tail_clip_t) = (0, 0);
        if($cigar =~ /^(\d+)[SH]/){
            $head_clip_t = $1;
            if($1 > $min_clip){
                $head_clip = $1;
                # clipped at 3' to the contig if reverse; 5' if forward
            }
        }

        if($cigar =~ /(\d+)[SH]$/){
            $tail_clip_t = $1;
            if($1 > $min_clip){
                $tail_clip = $1;
            }
        }

        my $cigar_tmp = $cigar;
        my $ctg_len = 0;
        my $num;
        my $alpha;
        while($cigar_tmp ne ""){
            ($num, $alpha, $cigar_tmp) = ($cigar_tmp =~ /^(\d+)([SMDIH])(.*)$/);
            $ctg_len += $num if($alpha =~ /[SHMI]/);
        }


        # should be both sides do not have large clip
        next if($head_clip != 0 || $tail_clip != 0);
        # add this line when needing to set -q (max clip) very large to allow contig end to be misassembled
        next if($ctg_len - $head_clip - $tail_clip < $head_clip + $tail_clip);
        my $bp_pos = $pos;
        my $pb_s = 0;
        my $pb_e = $pb_s + 1;
        my $line = "";
        my $check_next = 0;
        my $prev = 0;
        my $con_M = 0;
        # a tag to show that the current status is combined with the previous call
        my $combine = 0;

        # find all M's with size > sz_M. For each gap within two consecutive such M's, qualify it for pairHMM if
        # 1. There is a gap > sz_t; or
        # 2. The gap sum of one kind > sz_t, and the Ms inside the region (between two large M's) are < sz_m
        # Qualify the gap directly if it is the only one in the region and it is > sz_t. Do not need pairHMM. 
        # This is to 1) refine breakpoint; 2) avoid missing ones

        my $h_M;
        my $previous = -1;
        my ($ref_p, $ctg_p);
        ($ref_p, $ctg_p, $alpha, $num) = cigar::get_pos_from_cigar_input_consider_rev($_);
        my $total = scalar(@$alpha);
        for(my $i = 0; $i < $total; $i ++){
            if($num->[$i] > $sz_M && $alpha->[$i] eq "M"){
                $h_M->{$i} = 1;
            }
        }
        for(my $i = 0; $i < $total; $i ++){
            if($previous == -1 && defined $h_M->{$i}){
                $previous = $i;
            }
            elsif(defined $h_M->{$i}){
                # a gap formed
                # check if a single one, no need for pairHMM
                if($i - $previous == 2 && $alpha->[$i-1] =~ /[DI]/ && $num->[$i-1] > $sz_t){
                    if($r == 0){
                        $line = join("\t", $chr, $ref_p->[$i-1], $num->[$i-1], $alpha->[$i-1], $ctg_name, $ctg_p->[$i-1], $ctg_p->[$i]);
                    }
                    else{
                        $line = join("\t", $chr, $ref_p->[$i-1], $num->[$i-1], $alpha->[$i-1], $ctg_name, $ctg_p->[$i], $ctg_p->[$i-1]);
                    }

=cut
my $ctg_pos1;
                    my $ctg_pos2;
                    if($r == 0){
                        $ctg_pos1 = $ctg_p->[$i-1];
                        if($alpha->[$i] eq "D"){
                            $ctg_pos2 = $ctg_pos1 + 1; 
                        }
                        elsif($alpha->[$i] eq "I"){
                            $ctg_pos2 = $ctg_pos1 + $num->[$i-1];
                        }
                    }
                    else{
                        if($alpha->[$i] eq "D"){
                            $ctg_pos1 = $ctg_length - $ctg_p->[$i-1];
                            $ctg_pos2 = $ctg_pos1 + 1;
                        }
                        elsif($alpha->[$i] eq "I"){
                            $ctg_pos1 = $ctg_length - $ctg_p-> 
                            $ctg_MI - $num->[$i] + $ctg_start;
                            $ctg_pos2 = $ctg_length - $ctg_MI + $ctg_start;
                        }
                    }
=cut

push @lines, $line;
$previous = $i;
next;
                }
                elsif($i - $previous == 2 && ($num->[$i-1] <= $sz_t || $alpha->[$i-1] !~ /[DI]/)){
                    $previous = $i;
                    next;
                }
                else{
                    # see whether this region qualifies for pairHMM
                    my $cigar_substr = cigar::get_substr($cigar, $previous, $i);
                    my $select = 0;
                    my $deselect = 0;
                    my $indel_sum = 0;
                    for(my $j = $previous + 1; $j < $i; $j ++){
                        if($alpha->[$j] =~ /[ID]/ && $num->[$j] > $sz_t){
                            $select = 1;
                            last;
                        }
                        elsif($alpha->[$j] =~ /[ID]/){
                            $indel_sum += $num->[$j];
                        }
                        elsif($alpha->[$j] eq "M"){
                            if($num->[$j] > $sz_m){
                                # not select it if all gap < sz_t, sum > sz_t because there is at least one M > sz_m
                                $deselect = 1;
                            }
                        }
                    }
                    if($select == 1 || $indel_sum > $sz_t && $deselect != 1){ 
                        $DB::single = 1;
                        # qualify it for pairHMM
                        my $ref_str;
                        my $ctg_str;
                        # Note: flank is to ensure there are enough flank region (this number cannot be too large, because bwa mem will flank the total region with many INDELs if there are both insertion and deletion (complex)). 50bp seem to be good. However, flank is not good enough for reference. In forward strand, ther are cases where ref string is not long enough to cover the whole contig. This would cause forward strand does not have many complex INDELs. Putting sz_M together with flank as the padding, now it is good for the second example in SV companion 2015 Fig. 3b.  
                        $ref_str = "$chr:" . ($ref_p->[$previous + 1] - $flank - $sz_M) . "-" . ($ref_p->[$i] + $flank + $sz_M);
                        if($r == 1){
                            $ctg_str = "$ctg_name_for_fa:" . ($ctg_p->[$i] - $flank) . "-" . ($ctg_p->[$previous + 1] + $flank);
                        }
                        else{
                            $ctg_str = "$ctg_name_for_fa:" . ($ctg_p->[$previous + 1] - $flank) . "-" . ($ctg_p->[$i] + $flank);
                        }
                        $DB::single = 1;
                        #my $bs = &refind_aln_hmm_memEffi($ctg_fa, $ref_fa, $ctg_str, $ref_str, $r, 0, $sz_M, $sz_t, $ctg_name, $prefix, $modification_error, $microins_dist, $pairHMM_trans_matrix);
                        my $bs = &bwa_local($ctg_fa, $ref_fa, $ctg_str, $ref_str, $sz_M, $sz_t, $prefix, $bwa_M_threshold);
                        foreach my $line (@$bs){
                            push @lines, $line;
                        }
                    }
                    $previous = $i;
                }
            }
        }
    }
    close SAM;
    return \@lines;
}

# memory efficient way for pair HMM
sub infer_INDEL_PB_wPairHMM_memEffi{
    my ($sam, $min_clip, $sz_t, $sz_M, $sz_m, $ctg_fa, $ref_fa, $prefix, $modification_error, $microins_dist, $pairHMM_trans_matrix) = @_;
    my @lines;
    open SAM, "<$sam" or die $!;
    while(<SAM>){
        next if($_ =~ /^@/);
        my @a = split(/\t/, $_);
        my ($ctg_name, $tag, $chr, $pos, $qual, $cigar) = @a[0 .. 5];
        my $ctg_name_for_fa = $ctg_name;
        while($ctg_name_for_fa =~ /(.+)\/\d+_\d+/){
            $ctg_name_for_fa = $1;
        }
        $ctg_name = $ctg_name_for_fa;
        next if($cigar eq "*");

        my $r = 0;
        if($tag & 0x0010){
            # reverse
            $r = 1;
        }

        my ($head_clip, $tail_clip) = (0, 0);
        my ($head_clip_t, $tail_clip_t) = (0, 0);
        if($cigar =~ /^(\d+)[SH]/){
            $head_clip_t = $1;
            if($1 > $min_clip){
                $head_clip = $1;
                # clipped at 3' to the contig if reverse; 5' if forward
            }
        }

        if($cigar =~ /(\d+)[SH]$/){
            $tail_clip_t = $1;
            if($1 > $min_clip){
                $tail_clip = $1;
            }
        }

        my $cigar_tmp = $cigar;
        my $ctg_len = 0;
        my $num;
        my $alpha;
        while($cigar_tmp ne ""){
            ($num, $alpha, $cigar_tmp) = ($cigar_tmp =~ /^(\d+)([SMDIH])(.*)$/);
            $ctg_len += $num if($alpha =~ /[SHMI]/);
        }


        # should be both sides do not have large clip
        next if($head_clip != 0 || $tail_clip != 0);
        # add this line when needing to set -q (max clip) very large to allow contig end to be misassembled
        next if($ctg_len - $head_clip - $tail_clip < $head_clip + $tail_clip);
        my $bp_pos = $pos;
        my $pb_s = 0;
        my $pb_e = $pb_s + 1;
        my $line = "";
        my $check_next = 0;
        my $prev = 0;
        my $con_M = 0;
        # a tag to show that the current status is combined with the previous call
        my $combine = 0;

        # find all M's with size > sz_M. For each gap within two consecutive such M's, qualify it for pairHMM if
        # 1. There is a gap > sz_t; or
        # 2. The gap sum of one kind > sz_t, and the Ms inside the region (between two large M's) are < sz_m
        # Qualify the gap directly if it is the only one in the region and it is > sz_t. Do not need pairHMM. 
        # This is to 1) refine breakpoint; 2) avoid missing ones

        my $h_M;
        my $previous = -1;
        my ($ref_p, $ctg_p);
        ($ref_p, $ctg_p, $alpha, $num) = cigar::get_pos_from_cigar_input_consider_rev($_);
        my $total = scalar(@$alpha);
        for(my $i = 0; $i < $total; $i ++){
            if($num->[$i] > $sz_M && $alpha->[$i] eq "M"){
                $h_M->{$i} = 1;
            }
        }
        for(my $i = 0; $i < $total; $i ++){
            if($previous == -1 && defined $h_M->{$i}){
                $previous = $i;
            }
            elsif(defined $h_M->{$i}){
                # a gap formed
                # check if a single one, no need for pairHMM
                if($i - $previous == 2 && $alpha->[$i-1] =~ /[DI]/ && $num->[$i-1] > $sz_t){
                    if($r == 0){
                        $line = join("\t", $chr, $ref_p->[$i-1], $num->[$i-1], $alpha->[$i-1], $ctg_name, $ctg_p->[$i-1], $ctg_p->[$i]);
                    }
                    else{
                        $line = join("\t", $chr, $ref_p->[$i-1], $num->[$i-1], $alpha->[$i-1], $ctg_name, $ctg_p->[$i], $ctg_p->[$i-1]);
                    }

=cut
my $ctg_pos1;
                    my $ctg_pos2;
                    if($r == 0){
                        $ctg_pos1 = $ctg_p->[$i-1];
                        if($alpha->[$i] eq "D"){
                            $ctg_pos2 = $ctg_pos1 + 1; 
                        }
                        elsif($alpha->[$i] eq "I"){
                            $ctg_pos2 = $ctg_pos1 + $num->[$i-1];
                        }
                    }
                    else{
                        if($alpha->[$i] eq "D"){
                            $ctg_pos1 = $ctg_length - $ctg_p->[$i-1];
                            $ctg_pos2 = $ctg_pos1 + 1;
                        }
                        elsif($alpha->[$i] eq "I"){
                            $ctg_pos1 = $ctg_length - $ctg_p-> 
                            $ctg_MI - $num->[$i] + $ctg_start;
                            $ctg_pos2 = $ctg_length - $ctg_MI + $ctg_start;
                        }
                    }
=cut

push @lines, $line;
$previous = $i;
next;
                }
                elsif($i - $previous == 2 && ($num->[$i-1] <= $sz_t || $alpha->[$i-1] !~ /[DI]/)){
                    $previous = $i;
                    next;
                }
                else{
                    # see whether this region qualifies for pairHMM
                    my $cigar_substr = cigar::get_substr($cigar, $previous, $i);
                    my $select = 0;
                    my $deselect = 0;
                    my $indel_sum = 0;
                    for(my $j = $previous + 1; $j < $i; $j ++){
                        if($alpha->[$j] =~ /[ID]/ && $num->[$j] > $sz_t){
                            $select = 1;
                            last;
                        }
                        elsif($alpha->[$j] =~ /[ID]/){
                            $indel_sum += $num->[$j];
                        }
                        elsif($alpha->[$j] eq "M"){
                            if($num->[$j] > $sz_m){
                                # not select it if all gap < sz_t, sum > sz_t because there is at least one M > sz_m
                                $deselect = 1;
                            }
                        }
                    }
                    if($select == 1 || $indel_sum > $sz_t && $deselect != 1){ 
                        # qualify it for pairHMM
                        my $ref_str;
                        my $ctg_str;
                        $ref_str = "$chr:" . ($ref_p->[$previous + 1] - $sz_M) . "-" . ($ref_p->[$i] + $sz_M);
                        if($r == 1){
                            $ctg_str = "$ctg_name_for_fa:" . ($ctg_p->[$i] - $sz_M) . "-" . ($ctg_p->[$previous + 1] + $sz_M);
                        }
                        else{
                            $ctg_str = "$ctg_name_for_fa:" . ($ctg_p->[$previous + 1] - $sz_M) . "-" . ($ctg_p->[$i] + $sz_M);
                        }
                        #$DB::single = 1;
                        my $bs = &refind_aln_hmm_memEffi($ctg_fa, $ref_fa, $ctg_str, $ref_str, $r, 0, $sz_M, $sz_t, $ctg_name, $prefix, $modification_error, $microins_dist, $pairHMM_trans_matrix);
                        foreach my $line (@$bs){
                            push @lines, $line;
                        }
                    }
                    $previous = $i;
                }
            }
        }
    }
    close SAM;
    return \@lines;
}


# infer INDELs, this assumes the contig is with a relatively higher quality than PB, so that the large INDELs (size > 10bp) are real INDEL instead of INDEL error. sz_M is the left and right anchor (M) length for a I/D to qualify as a candidate
# for those indels with size > sz_t but not selected, find the two anchors of M with size > sz_M, extract the sequence from the ctg, run refind_aln_hmm.pl to perform pairHMM
sub infer_INDEL_PB_wPairHMM{
    my ($sam, $min_clip, $sz_t, $sz_M, $sz_m, $ctg_fa, $ref_fa, $prefix, $modification_error) = @_;
    my @lines;
    open SAM, "<$sam" or die $!;
    while(<SAM>){
        next if($_ =~ /^@/);
        my @a = split(/\t/, $_);
        my ($ctg_name, $tag, $chr, $pos, $qual, $cigar) = @a[0 .. 5];
        my $ctg_name_for_fa = $ctg_name;
        if($ctg_name =~ /(.+)\/\d+_\d+/){
            $ctg_name_for_fa = $1;
        }
        next if($cigar eq "*");

        my $r = 0;
        if($tag & 0x0010){
            # reverse
            $r = 1;
        }

        my ($head_clip, $tail_clip) = (0, 0);
        my ($head_clip_t, $tail_clip_t) = (0, 0);
        if($cigar =~ /^(\d+)[SH]/){
            $head_clip_t = $1;
            if($1 > $min_clip){
                $head_clip = $1;
                # clipped at 3' to the contig if reverse; 5' if forward
            }
        }

        if($cigar =~ /(\d+)[SH]$/){
            $tail_clip_t = $1;
            if($1 > $min_clip){
                $tail_clip = $1;
            }
        }

        my $cigar_tmp = $cigar;
        my $ctg_len = 0;
        my $num;
        my $alpha;
        while($cigar_tmp ne ""){
            ($num, $alpha, $cigar_tmp) = ($cigar_tmp =~ /^(\d+)([SMDIH])(.*)$/);
            $ctg_len += $num if($alpha =~ /[SHMI]/);
        }


        # should be both sides do not have large clip
        next if($head_clip != 0 || $tail_clip != 0);
        # add this line when needing to set -q (max clip) very large to allow contig end to be misassembled
        next if($ctg_len - $head_clip - $tail_clip < $head_clip + $tail_clip);
        my $bp_pos = $pos;
        my $pb_s = 0;
        my $pb_e = $pb_s + 1;
        my $line = "";
        my $check_next = 0;
        my $prev = 0;
        my $con_M = 0;
        # a tag to show that the current status is combined with the previous call
        my $combine = 0;

        # find all M's with size > sz_M. For each gap within two consecutive such M's, qualify it for pairHMM if
        # 1. There is a gap > sz_t; or
        # 2. The gap sum of one kind > sz_t, and the Ms inside the region (between two large M's) are < sz_m
        # Qualify the gap directly if it is the only one in the region and it is > sz_t. Do not need pairHMM. 
        # This is to 1) refine breakpoint; 2) avoid missing ones

        my $h_M;
        my $previous = -1;
        my ($ref_p, $ctg_p);
        ($ref_p, $ctg_p, $alpha, $num) = cigar::get_pos_from_cigar_input_consider_rev($_);
        my $total = scalar(@$alpha);
        for(my $i = 0; $i < $total; $i ++){
            if($num->[$i] > $sz_M && $alpha->[$i] eq "M"){
                $h_M->{$i} = 1;
            }
        }
        for(my $i = 0; $i < $total; $i ++){
            if($previous == -1 && defined $h_M->{$i}){
                $previous = $i;
            }
            elsif(defined $h_M->{$i}){
                # a gap formed
                # check if a single one, no need for pairHMM
                if($i - $previous == 2 && $alpha->[$i-1] =~ /[DI]/ && $num->[$i-1] > $sz_t){
                    if($r == 0){
                        $line = join("\t", $chr, $ref_p->[$i-1], $num->[$i-1], $alpha->[$i-1], $ctg_name, $ctg_p->[$i-1], $ctg_p->[$i]);
                    }
                    else{
                        $line = join("\t", $chr, $ref_p->[$i-1], $num->[$i-1], $alpha->[$i-1], $ctg_name, $ctg_p->[$i], $ctg_p->[$i-1]);
                    }

=cut
my $ctg_pos1;
                    my $ctg_pos2;
                    if($r == 0){
                        $ctg_pos1 = $ctg_p->[$i-1];
                        if($alpha->[$i] eq "D"){
                            $ctg_pos2 = $ctg_pos1 + 1; 
                        }
                        elsif($alpha->[$i] eq "I"){
                            $ctg_pos2 = $ctg_pos1 + $num->[$i-1];
                        }
                    }
                    else{
                        if($alpha->[$i] eq "D"){
                            $ctg_pos1 = $ctg_length - $ctg_p->[$i-1];
                            $ctg_pos2 = $ctg_pos1 + 1;
                        }
                        elsif($alpha->[$i] eq "I"){
                            $ctg_pos1 = $ctg_length - $ctg_p-> 
                            $ctg_MI - $num->[$i] + $ctg_start;
                            $ctg_pos2 = $ctg_length - $ctg_MI + $ctg_start;
                        }
                    }
=cut

push @lines, $line;
$previous = $i;
next;
                }
                elsif($i - $previous == 2 && ($num->[$i-1] <= $sz_t || $alpha->[$i-1] !~ /[DI]/)){
                    $previous = $i;
                    next;
                }
                else{
                    # see whether this region qualifies for pairHMM
                    my $cigar_substr = cigar::get_substr($cigar, $previous, $i);
                    my $select = 0;
                    my $deselect = 0;
                    my $indel_sum = 0;
                    for(my $j = $previous + 1; $j < $i; $j ++){
                        if($alpha->[$j] =~ /[ID]/ && $num->[$j] > $sz_t){
                            $select = 1;
                            last;
                        }
                        elsif($alpha->[$j] =~ /[ID]/){
                            $indel_sum += $num->[$j];
                        }
                        elsif($alpha->[$j] eq "M"){
                            if($num->[$j] > $sz_m){
                                # not select it if all gap < sz_t, sum > sz_t because there is at least one M > sz_m
                                $deselect = 1;
                            }
                        }
                    }
                    if($select == 1 || $indel_sum > $sz_t && $deselect != 1){ 
                        # qualify it for pairHMM
                        my $ref_str;
                        my $ctg_str;
                        $ref_str = "$chr:" . ($ref_p->[$previous + 1] - $sz_M) . "-" . ($ref_p->[$i] + $sz_M);
                        if($r == 1){
                            $ctg_str = "$ctg_name_for_fa:" . ($ctg_p->[$i] - $sz_M) . "-" . ($ctg_p->[$previous + 1] + $sz_M);
                        }
                        else{
                            $ctg_str = "$ctg_name_for_fa:" . ($ctg_p->[$previous + 1] - $sz_M) . "-" . ($ctg_p->[$i] + $sz_M);
                        }
                        my $bs = &refind_aln_hmm($ctg_fa, $ref_fa, $ctg_str, $ref_str, $r, 0, $sz_M, $sz_t, $ctg_name, $prefix, $modification_error);
                        foreach my $line (@$bs){
                            push @lines, $line;
                        }
                    }
                    $previous = $i;
                }
            }
        }
    }
    close SAM;
    return \@lines;
}


# infer INDELs, this assumes the contig is with a relatively higher quality than PB, so that the large INDELs (size > 10bp) are real INDEL instead of INDEL error. sz_M is the left and right anchor (M) length for a I/D to qualify as a candidate
sub infer_INDEL_PB{
    my ($sam, $min_clip, $sz_t, $sz_M) = @_;
    my @lines;
    open SAM, "<$sam" or die $!;
    while(<SAM>){
        next if($_ =~ /^@/);
        my @a = split(/\t/, $_);
        my ($ctg_name, $tag, $chr, $pos, $qual, $cigar) = @a[0 .. 5];
        next if($cigar eq "*");

        my $r = 0;
        if($tag & 0x0010){
            # reverse
            $r = 1;
        }

        my ($head_clip, $tail_clip) = (0, 0);
        my ($head_clip_t, $tail_clip_t) = (0, 0);
        if($cigar =~ /^(\d+)[SH]/){
            $head_clip_t = $1;
            if($1 > $min_clip){
                $head_clip = $1;
                # clipped at 3' to the contig if reverse; 5' if forward
            }
        }

        if($cigar =~ /(\d+)[SH]$/){
            $tail_clip_t = $1;
            if($1 > $min_clip){
                $tail_clip = $1;
            }
        }

        my $cigar_tmp = $cigar;
        my $ctg_len = 0;
        my $num;
        my $alpha;
        while($cigar_tmp ne ""){
            ($num, $alpha, $cigar_tmp) = ($cigar_tmp =~ /^(\d+)([SMDIH])(.*)$/);
            $ctg_len += $num if($alpha =~ /[SHMI]/);
        }


        # should be both sides do not have large clip
        next if($head_clip != 0 || $tail_clip != 0);
        # add this line when needing to set -q (max clip) very large to allow contig end to be misassembled
        next if($ctg_len - $head_clip - $tail_clip < $head_clip + $tail_clip);
        my $bp_pos = $pos;
        my $pb_s = 0;
        my $pb_e = $pb_s + 1;
        my $line = "";
        my $check_next = 0;
        my $prev = 0;
        my $con_M = 0;
        # a tag to show that the current status is combined with the previous call
        my $combine = 0;
        while($cigar ne ""){
            ($num, $alpha, $cigar) = ($cigar =~ /^(\d+)([SMDIH])(.*)$/);
            if($check_next == 1){
                if($con_M > $sz_M){
                    # left and right qualify M length for this I/D, add it
                    push @lines, $line;
                    $line = "";
                    $check_next = 0;
                    $prev = 0;
                }
                else{
                    if($alpha ne "M" && $num > $sz_t && $con_M < 10){
                        # encounter a big indel, the previous one hasn't succeeded with a large enough anchor, give up the previous one
                        # now combine this with the previous one if they are the same type
                        my @aa = split(/\t/, $line);
                        if($aa[3] eq $alpha){
                            if($alpha eq "I"){
                                if($r == 1){
                                    $aa[5] = $aa[5] - $num;
                                }
                                else{
                                    $aa[6] = $aa[6] + $num;
                                }
                            }
                            $line = join("\t", @aa[0 .. 1], ($aa[2] + $num), @aa[3 .. 6]);
                            $combine = 1;
                        }
                        else{
                            $line = "";
                            $check_next = 0;
                            $prev = 0;
                            $con_M = 0;
                        }
                    }
                    else{
                        # hasn't arrived, but still going
                        $con_M += $num if($alpha eq "M");
                    }
                }
            }
            if($num > $sz_t && $alpha =~ /[ID]/ && $prev > $sz_M && $combine != 1){
                if($alpha =~ /I/){
                    $pb_e = $pb_s + $num;
                }
                elsif($alpha =~ /D/){
                    $pb_e = $pb_s + 1;
                }
                # since pb_s should be used for the next infererence, discuss the cases directly to output line
                if($r == 1){
                    $line = join("\t", $chr, $bp_pos, $num, $alpha, $ctg_name, $ctg_len - $pb_e, $ctg_len - $pb_s);
                }
                else{
                    $line = join("\t", $chr, $bp_pos, $num, $alpha, $ctg_name, $pb_s, $pb_e);
                }
                #push @lines, $line;
                $bp_pos += $num if($alpha =~ /D/);
                $pb_s += $num if($alpha =~ /I/);
                $check_next = 1;
                $con_M = 0;
            }
            else{
                $bp_pos += $num if($alpha =~ /[MD]/);
                $pb_s += $num if($alpha =~ /[SMIH]/);
            }
            # accumulate from the first M after a large indel
            if($alpha eq "M"){
                $prev += $num;
            }
            elsif($alpha =~ /[ID]/){
                if($num > $sz_t && $combine != 1){
                    $prev = 0;
                }
            }
            $combine = 0;
        }
    }
    close SAM;
    return \@lines;
}

# 06052016: Last line: added remove conflict (chr, gap). For chr, always choose the one with the larger avg. supporting IL reads. For gap, always choose the larger one. 
# This function takes a sam, a reported breakpoint (for INDEL, in format of chr, pos, length, I/D, contig_name, starting position of this indel on the contig, and the end of this indel on the contig), and report the line if there are at least m short read pairs intersecting this breakpoint. $perc is the percentage of the bases that needs to match to the ctg from IL: for insertion, it's simply the insertion sequenced analyzed; for deletion, it is the [-s,s] surround the deleted breakpoint analyzed.
sub check_INDEL_bp{
    my ($cSAM, $line, $m, $perc, $s, $trio) = @_;
    my @r;
    foreach my $l (@$line){
        my $hh;
        chomp $l;
        my @a = split(/\t/, $l);
        if($a[3] =~ /[DI]/){
            # indel
            $hh->{chr} = $a[0];
            if($hh->{chr} =~ /\//){
                my @tmp = split(/\//, $hh->{chr});
                $hh->{chr} = join("/", @tmp[0 .. $#tmp - 1]);
            }
            $hh->{pos} = $a[1];
            $hh->{len} = $a[2];
            $hh->{tp} = $a[3];
            $hh->{ctg_name} = $a[4];
            if($hh->{ctg_name} =~ /\//){
                my @tmp = split(/\//, $hh->{ctg_name});
                $hh->{ctg_name} = join("/", @tmp[0 .. $#tmp - 1]);
            }
            #ctg_start is the coordinate on the contig, in the orientaiton of the contig
            $hh->{ctg_start} = $a[5];
            $hh->{ctg_end} = $a[6];
        }
        else{
            # SV
            $hh->{ctg_name} = $a[$#a-2];
            $hh->{ctg_start} = $a[$#a-1];
            $hh->{ctg_end} = $a[$#a];
            $hh->{tp} = $a[6];
        }
        my $c = 0;
        # count by hash sample name
        my $cc;
        #$DB::single = 1;
        if(defined $trio && $trio ne "NA"){
            my @samples = split(/:/, $trio);
            foreach my $sample (@samples){
                $cc->{$sample} = 0;
            }
        }

        # check the file
        my @f = split(/\./, $hh->{ctg_name});
        $hh->{ctg_name} = $f[0];
        # a hash to record the read pair number that support this INDEL
        my $hc;
        next if(! -e $cSAM);
        open cov_fh, "<$cSAM" or die $!;
        while(<cov_fh>){
            next if($_ =~ /^@/);
            my @a = split(/\s+/, $_);
            my ($ctg_name, $tag, $chr, $pos, $qual, $cigar) = @a[0 .. 5];
            next if($chr eq "*" || $cigar eq "*");
            ($chr) = ($chr =~ /^(g\d+)/);
            my $rp_name;
            if($ctg_name =~ /^(.+)\/[12]/){
                $rp_name = $1;
            }
            else{
                $rp_name = $ctg_name;
            }
            #my ($rp_name) = ($ctg_name =~ /^(.+)\/[12]/);
            $hc->{$rp_name} = cigar::test_indel($hh, $chr, $pos, $cigar, $perc, $s) if(!defined $hc->{$rp_name} || $hc->{$rp_name} == 0);
        }
        close cov_fh;
        # separate multi-sample and single sample
        my $flag = 0;
        foreach my $rp (keys %$hc){
            if($rp =~ /^(\d+)_(\d+)/){
                if($hc->{$rp} == 1){
                    $flag = 1;
                    # by sample
                    $cc->{$1} ++ if(defined $cc->{$1});
                }
            }
            else{
                $c ++ if($hc->{$rp} == 1);
            }
        }
        if(scalar(keys %$cc) != 0 && $flag != 0){
            push @r, "$l\t" . join("\t", map {$cc->{$_}} sort {$a <=> $b} keys %$cc);
        }
        else{
            if($c >= $m){
                push @r, "$l\t$c";
            }
        }
    }
    my $r1 = &remove_conflict(\@r);
    @r = @$r1;
    return \@r;
}

# remove conflict calls
# if two calls are on two chromosomes, but corresponding to the same contig, remove the one with the smaller supporting ILs
# if two calls are overlapping, take the longer one as long as it has at least 1 supporting read
sub remove_conflict{
    my ($r) = @_;
    my $h_overlap;
    my $h_chr;
    my @ret;
    foreach my $rr (@$r){
        my @x = split(/\t/, $rr);
        my $gap;
        if($x[3] eq "I"){
            $gap = "$x[1].$x[1]";
        }
        elsif($x[3] eq "D"){
            my $end = $x[1] + $x[2];
            $gap = "$x[1].$end";
        }
        $h_overlap->{$x[4]}->{$x[0]}->{$gap} = $rr;
        if(!defined $h_chr->{$x[4]}->{$x[0]}){
            $h_chr->{$x[4]}->{$x[0]} = $x[$#x];
        }
        else{
            $h_chr->{$x[4]}->{$x[0]} .= ";$x[$#x]";
        }
    }
    # y is the ctg name
    foreach my $y (keys %$h_chr){
        my $chosen_chr = "NA";
        if(scalar(keys %{$h_chr->{$y}}) >= 2){
            # two chromosomes
            my $max = 0;
            print STDERR "#Conflicting chr: $y\n";
            foreach my $chr (keys %{$h_chr->{$y}}){
                # supporting IL num for each call
                my @sup = split(/;/, $h_chr->{$y}->{$chr});
                my $avg = math::avg(\@sup);
                if($avg > $max){
                    $max = $avg;
                    $chosen_chr = $chr;
                }
            }
            print STDERR "#Choosing $chosen_chr\n";
        }
        else{
            my $tmp = (keys %{$h_chr->{$y}})[0];
            $chosen_chr = $tmp;
        }
        # continue to check overlap
        my @gaps = sort {$a <=> $b} keys %{$h_overlap->{$y}->{$chosen_chr}};
        my $to_ignore;
        for(my $i = 0; $i < scalar(@gaps) - 1; $i ++){
            # the gap that has been deemed ignored will not be able to make others ignored
            next if(defined $to_ignore->{$gaps[$i]});
            my ($start, $end) = split(/\./, $gaps[$i]);
            for(my $j = $i + 1; $j < scalar(@gaps); $j ++){
                my ($start1, $end1) = split(/\./, $gaps[$j]);
                if($end > $start1){
                    # overlap
                    # check which is larger
                    if($end1 - $start1 > $end - $start){
                        print STDERR "#Ignored overlap: " . $h_overlap->{$y}->{$chosen_chr}->{"$start.$end"} . "\n";
                        $to_ignore->{$gaps[$i]} = 1;
                        # ignored cannot ignore others anymore
                        next;
                    }
                    else{
                        print STDERR "#Ignored overlap: " . $h_overlap->{$y}->{$chosen_chr}->{"$start1.$end1"} . "\n";
                        $to_ignore->{$gaps[$j]} = 1;
                    }
                }
            }
        }
        # print only not ignored ones in chosen_chr
        foreach my $gap (@gaps){
            if(!defined $to_ignore->{$gap}){
                push @ret, $h_overlap->{$y}->{$chosen_chr}->{$gap};
            }
        }
    }

    return \@ret;
}







# detect INV, surrounded by anchored regions
sub infer_INV_bwa{
    my ($sam, $gap_t, $min_anchor) = @_;
    my $h;
    open SAM, "<$sam" or die $!;
    while(<SAM>){
        next if($_ =~ /^@/);
        my @a = split(/\t/, $_);
        $h->{$a[0]}->{$a[2]}->{$a[3]} = "$a[1].$a[5]";
    }
    close SAM;
    # output in lines (format: chr pos len type readname read_start_pos read_end_pos)
    my @lines;
    foreach my $rn (keys %$h){
        foreach my $chr (keys %{$h->{$rn}}){
            # inv should have three segments of alignment
            next if(scalar(keys %{$h->{$rn}->{$chr}}) < 3);
            my @poss = ();
            my @cigars = ();
            my @matchlens = ();
            my @startposs = ();
            my @endposs = ();
            my $flip = 0;
            my $pre_ori = "NA";
            foreach my $pos (sort keys %{$h->{$rn}->{$chr}}){
                my ($tag, $cigar) = split(/\./, $h->{$rn}->{$chr}->{$pos});
                # 1 stand for forward
                my $ori = 1;
                $ori = -1 if(($tag & 0x0010) != 0);
                if($pre_ori ne "NA" && $pre_ori * $ori == -1){
                    # count the flip
                    $flip ++;
                    push @poss, $pos;
                    push @cigars, $cigar;
                    push @matchlens, cigar::get_matchlen_on_read($cigar);
                    push @startposs, cigar::get_startpos_on_read($cigar, $ori);
                    push @endposs, cigar::get_endpos_on_read($cigar, $ori);
                }
                else{
                    $flip = 0;
                    @poss = ();
                    push @poss, $pos;
                    @cigars = ();
                    push @cigars, $cigar;
                    push @matchlens, cigar::get_matchlen_on_read($cigar);
                    push @startposs, cigar::get_startpos_on_read($cigar, $ori);
                    push @endposs, cigar::get_endpos_on_read($cigar, $ori);
                }
                $pre_ori = $ori;
                if($flip == 2){
                    #$DB::single = 1;
                    # check match length, overlap/gap length
                    if(!defined $poss[2] || $matchlens[0] < $min_anchor || $matchlens[2] < $min_anchor){
                        last;
                    }
                    # check the gap/overlap on the read
                    if(abs($endposs[0] - $startposs[1]) < $gap_t && abs($endposs[1] - $startposs[2]) < $gap_t && abs($poss[1] - $poss[0] - $matchlens[0]) < $gap_t && abs($poss[2] - $poss[1] - $matchlens[1]) < $gap_t ){
                        my $line = join("\t", $chr, $poss[1], cigar::get_matchlen_on_ref($cigar), "INV", $rn, $startposs[1], $endposs[1]);
                        push @lines, $line;
                    }
                    last;
                }
            } # end of pos
        } # end of chr
    } # end of readname
    return \@lines;
} # end of function

# given a sam file, for each contig, look for the large gaps or insertions not covered in a separate alignment, extract the sequence and realign to the reference
sub fill_in_gaps{
    # max_gap is the minimum gap size for the gap to be realigned
    # insert_off is the maximum distance between either breakpoint of insertion (two alignments, one with I, the other aligned the inserted part) for the insertion to skip realignment
    # min_insert is the minimum insert size for this insertion to be considered
    my ($sam, $fa, $min_gap, $insert_off, $min_insert) = @_;
    open fh_, "<$sam" or die $!;

    # insertion range
    my $x = 0;
    my $y = 0;
    # for recotding insertion coordiantes
    my @inserts;
    # match range (start and stop at soft clip)
    my $s = 0;
    my $e = 0;
    # reverse (1) or not (0)
    my $r = 0;
    # read name
    my $pre = "NA";
    # a running position on cigar 
    my $q_pos;
    # for cigar analysis
    my $num;
    my $alpha;
    my $cigar;
    # where the reads will be written
    my $out = "$sam.fa";
    if(-e $out){
        `rm $out`;
    }
    # readname
    my $name;
    # record start to end intervals
    my @regions;
    # insert position
    my @insert_p;
    while(<fh_>){
        $r = 0;
        $x = 0;
        $y = 0;
        $q_pos = 0;
        next if($_ =~ /^@/);
        my @a = split(/\t/, $_);
        ($name, ) = split(/\//, $a[0]);
        $cigar = $a[5];
        if($pre ne "NA" && $name ne $pre){
            my $gaps = &analyze_regions(\@regions, \@inserts, $min_gap, $insert_off);
            my $h_gaps = &get_seqs($gaps, $pre);
            fa::write_fa_extract($h_gaps, [keys %$h_gaps], $out, $fa, "a");
            foreach my $key (keys %$h_gaps){
                delete $h_gaps->{$key};
            }
            #$DB::single = 1 if(scalar(@$gaps) != 0);
            @regions = ();
            @insert_p = ();
            @inserts = ();
        }
        $pre = $name;
        if($a[1] & 0x0010){
            $r = 1;
        }
        while($cigar ne ""){
            ($num, $alpha, $cigar) = ($cigar =~ /^(\d+)([SMDIH])(.*)$/);
            if($alpha =~ /[SHMI]/){
                $q_pos += $num;
            }
            if($num > $min_insert && $alpha eq "I"){
                my $tmp = $q_pos - $num;
                push @insert_p, "$tmp.$q_pos";
            }
        }
        $cigar = $a[5];
        if($cigar =~ /^(\d+)[SH]/){
            $s = $1;
        }
        else{
            $s = 0;
        }
        if($cigar =~ /(\d+)[SH]$/){
            $e = $1;
        }
        else{
            $e = 0;
        }
        if($r == 0){
            $e = $q_pos - $e;
        }
        else{
            my $tmp = $s;
            $s = $e;
            $e = $q_pos - $tmp;
        }
        if(scalar(@insert_p) != 0){
            foreach my $i (@insert_p){
                ($x, $y) = split(/\./, $i);
                if($r == 1){
                    my $tmp = $x;
                    $x = $q_pos - $y;
                    $y = $q_pos - $tmp;
                }
                my $exist = 0;
                foreach my $tmp (@inserts){
                    if($tmp eq "$x.$y"){
                        $exist = 1;
                    }
                }
                push @inserts, "$x.$y" if($exist == 0);
            }
            @insert_p = ();
        }
        my $exist = 0;
        foreach my $tmp (@regions){
            if($tmp eq "$s.$e"){
                $exist = 1;
            }
        }
        push @regions, "$s.$e" if($exist == 0);
    }
    $DB::single = 1;
    my $gaps = &analyze_regions(\@regions, \@inserts, $min_gap, $insert_off);
    my $h_gaps = &get_seqs($gaps, $name);
    fa::write_fa_extract($h_gaps, [keys %$h_gaps], $out, $fa, "a");
    $DB::single = 1 if(scalar(@$gaps) != 0);
    foreach my $key (keys %$h_gaps){
        delete $h_gaps->{$key};
    }
    return 1;
}

# given the gap name,  
# output is a hash table with the name the contig names (with gap regions), and the value the sequence, extracted from seq 
sub get_seqs{
    my ($gaps, $readnm) = @_;
    my $h;
    $h->{$readnm} = $gaps;
    return $h;
}

# input is the regions that have been covered, and the insert regions
# It looks for the gap that has not been covered, and the inserted sequence that does not have a separate alignment
# output is an array of all the regions that need to be realigned
sub analyze_regions{
    my ($regions, $insert_h, $min_gap, $insert_off) = @_;
    # output: gaps
    my @gaps;
    # look at insertions 
    foreach my $ins (@$insert_h){
        my $exist = 0;
        my ($ins_s, $ins_e) = split(/\./, $ins);
        foreach my $r (@$regions){
            my ($s, $e) = split(/\./, $r);
            if(abs($ins_s - $s) < $insert_off && abs($ins_e - $e) < $insert_off){
                $exist = 1;
                last;
            }
        }
        if($exist == 0){
            push @gaps, "$ins_s.$ins_e";
        }
    }

    my $n = 0;
    my $pre_e;
    foreach my $r (sort {$a <=> $b} @$regions){
        $n ++;
        my ($s, $e) = split(/\./, $r); 
        if($n == 1){
            $pre_e = $e;
            next;
        }

        if($pre_e >= $e){
            next;
        }
        elsif($s <= $pre_e){
            if($pre_e < $e){
                $pre_e = $e;
            }
        }
        elsif($s > $pre_e && $s - $pre_e > $min_gap){
            push @gaps, "$pre_e.$s";
            $pre_e = $e;
        }
    }
    return \@gaps;
}


# different from infer_SV_PB directly, this read out the gaps and insertions that are not covered by separate alignments, align the sequences to the reference, and combine the good alignments with the original ones. For an alignment containing a long insertion, separate it into two with corresponding position and cigar string. In this case, our infer_SV_PB program can work to connect the separately aligned inserted sequence and the two breakpoints.
sub infer_SV_PB_gaps{
    my ($sam, $min_clip, $max_gap, $insert_gap, $fa, $min_gap, $insert_off, $min_insert, $blasr, $ref) = @_;
    $DB::single = 1;
    # find the gaps
    &fill_in_gaps($sam, $fa, $min_gap, $insert_off, $min_insert);
    # align the gaps to the reference
    `$blasr -maxAnchorsPerPosition 100 -advanceExactMatches 10 -affineAlign -affineOpen 100 -affineExtend 0 -insertion 5 -deletion 5 -extend -maxExtendDropoff 20 -clipping subread -bestn 1 -sam -nproc 1 -out $sam.1 $sam.fa $ref`;
    #print "$blasr -maxAnchorsPerPosition 100 -advanceExactMatches 10 -affineAlign -affineOpen 100 -affineExtend 0 -insertion 5 -deletion 5 -extend -maxExtendDropoff 20 -clipping subread -bestn 1 -sam -nproc 1 -out $sam.1 $sam.fa $ref\n";
    # store the length of each read, only once
    my $h_len;
    # store everything in the hash h
    my $h;
    # combine the original sam and the new sam together, sort by readname, break up at the insertion
    my ($num, $alpha, $cigar);
    # a running cigar for breaking at large insertion
    my $this_cigar = "";
    # a running length of each segment that are in between large insertions
    my $this_length = 0;
    # a running length of each segment on the contig
    my $ctg_this_length = 0;
    # a running position of the start of each segment for breakpoint between large insertions
    my $pos;
    open fh_, "<$sam" or die $!;
    while(<fh_>){
        if($_ =~ /^@/){
            next;
        }
        my @a = split(/\t/, $_);
        my ($name, ) = split(/\//, $a[0]);
        my $cigar = $a[5];
        my $pos = $a[3];
        $this_cigar = "";
        $this_length = 0;
        $ctg_this_length = 0;
        my $length = 0;
        while($cigar ne ""){
            ($num, $alpha, $cigar) = ($cigar =~ /^(\d+)([SMDIH])(.*)$/);
            if($alpha =~ /[SMIH]/){
                $length += $num;
            }
        }
        $h_len->{$name} = $length if(!defined $h_len->{$name});
        $cigar = $a[5];

        my $pre_ctg_this_length = 0;

        while($cigar ne ""){
            ($num, $alpha, $cigar) = ($cigar =~ /^(\d+)([SMDIH])(.*)$/);
            if($alpha =~ /[MD]/){
                $this_length += $num;
            }
            if($num > $min_insert && $alpha eq "I"){
                # print the previous ones
                my $hh;
                # see if ever coming to this if condition before, if yes, there should be a head soft clip taking account of the cigars not included in this segment
                if($pre_ctg_this_length != 0){
                    $hh->{cigar} = $pre_ctg_this_length . "S";
                }
                else{
                    $hh->{cigar} = "";
                }
                $hh->{cigar} .= $this_cigar;
                # need to take care of the end by padding the necessary soft clips 
                if($h_len->{$name} - $ctg_this_length != 0){
                    $hh->{cigar} .= ($h_len->{$name} - $ctg_this_length) . "S";
                }
                $hh->{tag} = $a[1];
                $hh->{chr} = $a[2];
                $hh->{pos} = $pos;
                push @{$h->{$name}}, $hh;
                $this_cigar = "";
                # running position on the reference
                $pos += $this_length;
                # total length of a segment on the reference
                $this_length = 0;
                #$ctg_this_length += $num;
                # update the start with the end, but should consider the insertion in the previous line
                $ctg_this_length += $num;
                $pre_ctg_this_length = $ctg_this_length;
            }
            # update the running ctg length, but be careful not to double count
            if($alpha =~ /[MISH]/ && ! ($alpha eq "I" && $num > $min_insert)){
                $ctg_this_length += $num;
            }
            if($alpha =~ /[MISHD]/ && ! ($alpha eq "I" && $num > $min_insert)){
                $this_cigar .= $num . $alpha;
            }
        }
        my $hh;
        # take care of the head padding if ever going into the if condition before
        if($pre_ctg_this_length != 0){
            $hh->{cigar} = $pre_ctg_this_length . "S";
        }
        else{
            $hh->{cigar} = "";
        }
        # no need to take care of the end because it is already at the end
        $hh->{cigar} .= $this_cigar;
        $hh->{tag} = $a[1];
        $hh->{chr} = $a[2];
        $hh->{pos} = $pos;
        push @{$h->{$name}}, $hh;
    }
    close fh_;
    if(-e "$sam.1"){
        open fh_, "<$sam.1" or die $!;
        while(<fh_>){
            if($_ =~ /^@/){
                next;
            }
            my @a = split(/\t/, $_);
            my ($name_, ) = split(/\//, $a[0]);
            my ($name1, $name2, $s, $e) = split(/\./, $name_);
            my $name = join(".", $name1, $name2);
            my $hh;
            # utilize the coordinates on the name to pad the cigar string
            $hh->{cigar} = $s . "S" . $a[5] . ($h_len->{$name} - $e) . "S";
            $hh->{tag} = $a[1];
            $hh->{chr} = $a[2];
            $hh->{pos} = $a[3];
            push @{$h->{$name}}, $hh;
        }
        close fh_;
    }
    my $out_sam = "$sam.2";
    open SAM, ">$out_sam" or die $!;
    foreach my $name (keys %$h){
        foreach my $hh (@{$h->{$name}}){
            print SAM join("\t", $name, $hh->{tag}, $hh->{chr}, $hh->{pos}, "99", $hh->{cigar}) . "\n";
        }
    }
    close SAM;
    my $result = &infer_SV_PB($out_sam, $min_clip, $max_gap, $insert_gap);
    return $result;
}








# infer SV for Pacbio data, starts from DEL
# add ctg_name, bp_start on ctg, bp_end on ctg as the last three columns for all SV types, with the orientation w.r.t. the original sequence of the ctg(read) 
sub infer_SV_PB{
    # min_clip is the minimum size of clip to be counted; max_gap is the maximum gap on the contig between a pair of alignments
    my ($sam, $min_clip, $max_gap, $insert_gap) = @_;
    if(!defined $max_gap){
        $max_gap = 200;
    }
    if(!defined $insert_gap){
        $insert_gap = 100;
    }
    my $min_aln_len = 500;
    my $debug = 0;
    my $h;
    my $hi;
    my $flagged;
    open SAM, "<$sam" or die $!;
    while(<SAM>){
        next if($_ =~ /^@/);
        # convert D to M, ignore I. Generate two breakpoints if both ends are splitted, including both hard and soft.
        my @a = split(/\t/, $_);
        my ($ctg_name, $tag, $chr, $pos, $qual, $cigar) = @a[0 .. 5];
        next if($cigar eq "*" || defined $flagged->{$ctg_name});

        my $r = 0;
        if($tag & 0x0010){
            # reverse
            $r = 1;
        }

        my ($head_clip, $tail_clip) = (0, 0);
        if($cigar =~ /^(\d+)[SH]/){
            if($1 > $min_clip){
                $head_clip = $1;
                # clipped at 3' to the contig if reverse; 5' if forward
            }
        }

        if($cigar =~ /(\d+)[SH]$/){
            if($1 > $min_clip){
                $tail_clip = $1;
            }
        }

        my ($ref_len, $ctg_len) = (0, 0);
        my ($num, $alpha);
        my $pos2 = $pos; # reference position where 3' of ref was not aligned
        my $pos1 = $pos; # reference position where 5' of ref was not aligned
        while($cigar ne ""){
            ($num, $alpha, $cigar) = ($cigar =~ /^(\d+)([SMDIH])(.*)$/);
            #die "alpha not defined in $sam\n" if(!defined $alpha);
            $ref_len += $num if($alpha =~ /[SHMD]/);
            $ctg_len += $num if($alpha =~ /[SHMI]/);
            if($tail_clip != 0){
                $pos2 += $num if($alpha =~ /[MD]/);
            }
        }

        next if($ctg_len - $tail_clip - $head_clip < $min_aln_len);
        # These two are w.r.t. the ref cigar, from 5' on head, and from 3' on tail
        my $tail_clip0 = $tail_clip;
        my $tail_clip1 = $tail_clip;
        my $head_clip1 = $head_clip;

        # This has the orientation w.r.t. the contig, from 5' to 3'
        $tail_clip1 = $ctg_len - $tail_clip if($tail_clip != 0);
        # convert the clip position on the ctg to forward if reverse
        if($r == 1){
            $head_clip = $ctg_len - $tail_clip1 if($tail_clip1 != 0);
            $head_clip = 0 if($tail_clip1 == 0);
            $tail_clip = $ctg_len - $head_clip1 if($head_clip1 != 0);
            $tail_clip = 0 if($head_clip1 == 0);
        }
        # Debugged on 08132015, especially important for DEL
        else{
            $head_clip = $head_clip1;
            $tail_clip = $tail_clip1;
        }

        if($head_clip != 0){
            # clipped on 5' of the ctg
            my $ori = 5;
            my $t_pos = $pos1;
            if($r == 1){
                $ori = 3;
                $t_pos = $pos2;
            }
            push @{$h->{$ctg_name}->{5}->{$head_clip}}, "$chr:$t_pos:$ori";
            #print "Before in head_clip != 0, ori = $ori: " . join("\t", keys %{$hi->{$ctg_name}->{$ori}->{$chr}}) . "\n";
            #print $_;
            push @{$hi->{$ctg_name}->{$ori}->{$chr}->{$t_pos}}, "$head_clip:5";
            #print "After in head_clip != 0, ori = $ori: " . join("\t", keys %{$hi->{$ctg_name}->{$ori}->{$chr}}) . "\n";
        }
        if($tail_clip != 0){
            # clipped on 3' of the ctg
            my $ori = 3;
            my $t_pos = $pos2;
            if($r == 1){
                $ori = 5;
                $t_pos = $pos1;
            }
            push @{$h->{$ctg_name}->{3}->{$tail_clip}}, "$chr:$t_pos:$ori";
            #print "Before in tail_clip != 0, ori = $ori: " .  join("\t", keys %{$hi->{$ctg_name}->{$ori}->{$chr}}) . "\n";
            #print $_;
            push @{$hi->{$ctg_name}->{$ori}->{$chr}->{$t_pos}}, "$tail_clip:3";
            #print "After in tail_clip != 0, ori = $ori: " . join("\t", keys %{$hi->{$ctg_name}->{$ori}->{$chr}}) . "\n";
        }

    }
    my $tp = "NA";
    my @lines = ();
# compare clip position to infer SV
    foreach my $ctg_name (sort keys %$h){

        # report insertion
        my $new_hi = $hi->{$ctg_name};
        foreach my $chr (sort keys %{$hi->{$ctg_name}->{3}}){
            my $hi1 = $hi->{$ctg_name}->{3}->{$chr};
            my $hi2 = $hi->{$ctg_name}->{5}->{$chr};
            foreach my $p1 (sort keys %$hi1){
                foreach my $p2 (sort keys %$hi2){
                    # see the gap size of the inserted positions on the reference
                    if(abs($p2 - $p1) < $insert_gap){
                        my @x = @{$hi1->{$p1}};
                        my @y = @{$hi2->{$p2}};
                        foreach my $x_ (@x){
                            my ($pos_on_ctg1, $ori_on_ctg1) = split(/:/, $x_);
                            foreach my $y_ (@y){
                                my ($pos_on_ctg2, $ori_on_ctg2) = split(/:/, $y_);
                                # look at orientation
                                if(defined $ori_on_ctg1 && defined $ori_on_ctg2 && $ori_on_ctg1 != $ori_on_ctg2){
                                    #($ori_on_ref1 eq "-" && $ori_on_ref2 eq "+" || $ori_on_ref1 eq "+" && $ori_on_ref2 eq "-")){
                                    my $insert_size = abs($pos_on_ctg1 - $pos_on_ctg2);
                                    if($pos_on_ctg2 < $pos_on_ctg1){
                                        my $tmp = $pos_on_ctg1;
                                        $pos_on_ctg1 = $pos_on_ctg2;
                                        $pos_on_ctg2 = $tmp;
                                    }
                                    # now good to call it an insertion
                                    $DB::single = 1;
                                    my $line = join("\t", $chr, $p1, "+", $chr, $p2, "-", "INS", "99", $insert_size, $ctg_name, $pos_on_ctg1, $pos_on_ctg2);
                                    push @lines, $line;
                                }
                            }
                        }
                    }
                }
            }
        }
        # report other types
        my $new_h = $h->{$ctg_name};
        foreach my $clip_pos_5 (sort {$a <=> $b} keys %{$new_h->{5}}){
            foreach my $clip_pos_3 (sort {$a <=> $b} keys %{$new_h->{3}}){
                if(abs($clip_pos_5 - $clip_pos_3) < $max_gap){
                    $DB::single = 1;
                    #print join("\t", $clip_pos_5, $clip_pos_3)."\n";
                    #$DB::single = 1;
                    # this is a pair
                    my $size;
                    foreach my $x_ (@{$new_h->{5}->{$clip_pos_5}}){
                        foreach my $y_ (@{$new_h->{3}->{$clip_pos_3}}){
                            my ($chr1, $pos1, $ori1) = split(/:/, $x_);
                            my ($chr2, $pos2, $ori2) = split(/:/, $y_);
                            if($chr1 eq $chr2){
                                # the same chromosome
                                $size = abs($pos2 - $pos1);
                                if($ori1 == $ori2){
                                    # inversion
                                    $tp = "INV";
                                }
                                else{
                                    if($pos1 < $pos2 && $ori1 == 3 || $pos1 > $pos2 && $ori2 == 3){
                                        $tp = "DEL";
                                    }
                                    else{
                                        $tp = "DUP";
                                    }
                                }
                            }
                            else{
                                $tp = "CTX";
                                # a default inter-chromosomal size
                                $size = 10000;
                            }
                            # change 1 and 2 so that pos1 is always < pos2
                            if($pos1 > $pos2 && $chr1 eq $chr2){
                                my $ori_t = $ori1;
                                $ori1 = $ori2;
                                $ori2 = $ori_t;
                                my $pos_t = $pos1;
                                $pos1 = $pos2;
                                $pos2 = $pos_t;
                            }
                            if($clip_pos_5 > $clip_pos_3){
                                my $tmp = $clip_pos_3;
                                $clip_pos_3 = $clip_pos_5;
                                $clip_pos_5 = $tmp;
                            }
                            $ori1 = $ori1 == 5 ? "+":"-";
                            $ori2 = $ori2 == 5 ? "+":"-";
                            my $line = join("\t", $chr1, $pos1, $ori1, $chr2, $pos2, $ori2, $tp, "99", $size, $ctg_name, $clip_pos_5, $clip_pos_3);
                            push @lines, $line;
                        }
                    }
                }
            }
        }
    }

    close SAM;
    return \@lines;

}

# This is probably a temporary function. Given a small sam, it checks if all the left clipped (including indel cigars) are the same, and does this for the right clipped (including indels). For each left bp, print the number of supporting reads; the same for the right bp.  
sub bp_consistency{
    my ($sam, $indel_t, $clip_t) = @_;
    #$DB::single = 1;
    my $h_indel = &indel_from_cigar($sam, $indel_t);
    my $mute = 1;
    my $h_clip = &bp_from_cigar($sam, $clip_t, $mute);
    # assort left and right breakpoint
    my $h_l;
    my $h_r;
    my $h_l_indel;
    my $h_r_indel;
    my $l_max = -1;
    my $r_max = -1;
    my $l_pos = -1;
    my $r_pos = -1;
    my $debug = 0;
    my $interest_range = 20;
    foreach (keys %$h_indel){
        my ($chr, $indel_type, $start_p, $indel_size) = split(/:/, $_);
        if($indel_type eq "D"){
            if(defined $h_l->{"$chr:$start_p"}){
                $h_l_indel->{"$chr:$start_p"} += $h_indel->{$_};
            }
            else{
                $h_l_indel->{"$chr:$start_p"} = $h_indel->{$_};
            }
            my $right = $start_p + $indel_size;
            if(defined $h_r->{"$chr:$right"}){
                $h_r_indel->{"$chr:$right"} += $h_indel->{$_};
            }
            else{
                $h_r_indel->{"$chr:$right"} = $h_indel->{$_};
            }
            if($l_max < $h_l_indel->{"$chr:$start_p"}){
                $l_max = $h_l_indel->{"$chr:$start_p"};
                $l_pos = $start_p;
            }
            if($r_max < $h_r_indel->{"$chr:$right"}){
                $r_max = $h_r_indel->{"$chr:$right"};
                $r_pos = $right;
            }
        }
        elsif($indel_type eq "I"){
            if(defined $h_l_indel->{"$chr:$start_p"}){
                $h_l_indel->{"$chr:$start_p"} += $h_indel->{$_};
            }
            else{
                $h_l_indel->{"$chr:$start_p"} = $h_indel->{$_};
            }
            if(defined $h_r_indel->{"$chr:$start_p"}){
                $h_r_indel->{"$chr:$start_p"} += $h_indel->{$_};
            }
            else{
                $h_r_indel->{"$chr:$start_p"} = $h_indel->{$_};
            }
            if($l_max < $h_l_indel->{"$chr:$start_p"}){
                $l_max = $h_l_indel->{"$chr:$start_p"};
                $l_pos = $start_p;
            }
            if($r_max < $h_r_indel->{"$chr:$start_p"}){
                $r_max = $h_r_indel->{"$chr:$start_p"};
                $r_pos = $start_p;
            }
        }
    }
    #$DB::single = 1;
    foreach (keys %$h_clip){
        my ($chr, $b1, $b2) = split(/:/, $_);
        #my $q = 0;
        #foreach my $qual (@{$h_clip->{$_}}){
        #    $q += $qual;
        #}
        #if($q/scalar(@{$h_clip->{$_}}) > $mean_qual_t){
        if($b1 != -1){
            if(defined $h_r->{"$chr:$b1"}){
                $h_r->{"$chr:$b1"} += $h_clip->{$_};
            }
            else{
                $h_r->{"$chr:$b1"} = $h_clip->{$_};
            }
            if($debug == 1){
                print "\t" . $h_r->{"$chr:$b1"};
            }
            if($r_max < $h_r->{"$chr:$b1"}){
                $r_max = $h_r->{"$chr:$b1"};
                $r_pos = $b1;
            }
        }
        if($debug == 1){
            print "\n";
        }
        if($b2 != -1){
            if(defined $h_l->{"$chr:$b2"}){
                $h_l->{"$chr:$b2"} += $h_clip->{$_};
            }
            else{
                $h_l->{"$chr:$b2"} = $h_clip->{$_};
            }
            if($debug == 1){
                print "\t" . $h_l->{"$chr:$b2"};
            }
            if($l_max < $h_l->{"$chr:$b2"}){
                $l_max = $h_l->{"$chr:$b2"};
                $l_pos = $b2;
            }
        }
        if($debug == 1){
            print "\n";
        }
        #}
    }
    if($l_max <= 1 || $r_max <= 1){
        print "\t" . join("\t", -1, -1) . "\n";
    }
    elsif(scalar(keys %$h_l_indel) == 0 && scalar(keys %$h_r_indel) == 0){
        print "\t" . join("\t", -2, -2) . "\n";
    }
    else{
        my $count_l = 0;
        my $count_r = 0;
        foreach (keys %$h_l_indel){
            my ($chr_nb, $b_nb) = split(/:/, $_);
            if(abs($b_nb - $l_pos) < $interest_range && $b_nb != $l_pos){
                $count_l ++;
            }
        }
        foreach (keys %$h_r_indel){
            my ($chr_nb, $b_nb) = split(/:/, $_);
            if(abs($b_nb - $r_pos) < $interest_range && $b_nb != $r_pos){
                $count_r ++;
            }
        }
        print "\t" . join("\t", $count_l, $count_r) . "\n";

        #print "\t" . join("\t", scalar(keys %$h_l), scalar(keys %$h_r)) . "\n";
        #if($debug == 1){
        #    print "Left:" . join("\t", keys %$h_l) . "\n" . "Right:" .  join("\t", keys %$h_r) . "\n";
        #}
    }
=cut
print "Left:\n";
if(scalar(keys %$h_l) > 1){
    foreach my $k (keys %$h_l){
        print $k . "\t" . join(":", $h_l->{$k}) . "\n";
    }
}
print "Right:\n";
if(scalar(keys %$h_r) > 1){
    foreach my $k (keys %$h_r){
        print $k . "\t" . join(":", $h_r->{$k}) . "\n";
    }
}
=cut
}

# read a sam, report if any INDEL in the cigar string with size greater than indel_t. If yes, report the number, INDEL positions, INDEL type for each distinct combinations. 
sub indel_from_cigar{
    # This is now designed for Illumina reads, which does not have lots of INDELs in alignment.
    my ($sam, $indel_t) = @_;
    my $h;
    open SAM, "<$sam" or die $!;
    while(<SAM>){
        next if($_ =~ /^@/);
        chomp;
        my @a = split(/\t/, $_);
        next if($a[1] & 0x0800 || $a[5] eq "*" || $a[4] == 0);
        my $b = cigar::indel_from_cigar($a[5], $a[3], $indel_t);
        if($b ne "NA"){
            my $s = "$a[2]:$b";
            if(defined $h->{$s}){
                $h->{$s} ++;
            }
            else{
                $h->{$s} = 1;
            }
        }
    }
    close SAM;
    return $h;
}

# this applies to novel insertions, since it does not try to pair of the cigar
sub bp_from_cigar{
    # This is mainly for Pacbio related sequence
    my ($sam, $clip_t, $mute) = @_;
    my $h;
    open SAM, "<$sam" or die $!;
    #$DB::single = 1;
    while(<SAM>){
        next if($_ =~ /^@/);
        chomp;
        my @a = split(/\t/, $_);
        next if($a[1] & 0x0800 || $a[5] eq "*" || $a[4] == 0);
        my @b = cigar::bp_from_cigar($a[5], $a[3], $clip_t);
        foreach my $bb (@b){
            if($bb != -1 && (!defined $mute || $mute != 1)){
                print join(":", $a[2], $bb) . "\n";
            }
        }
        my $str = join(":", $a[2], @b);
        if(defined $h->{$str}){
            $h->{$str} ++;
        }
        else{
            $h->{$str} = 1;
        }
    }
    close SAM;
    return $h;
}

# look at both match/ref and match/ctg, > per. Good for validating alternative allele constructed with del, mapped to assembled sequence. 
# Output: #well_aligned\tch1:start1-end1;start2-end2\tnext call (end1 and start2 are the two breakpoitns, with the flanking region in start1-end1, and start2-end2)
# classify reads in a sam file into well-aligned and poor-aligned 
sub classify_aln_del{
    my ($sam, $per, $min_clip, $indel_size) = @_;
    my $h_good;
    my $h;
    open sam_fh, "<$sam" or die $!;
    while(<sam_fh>){
        next if($_ =~ /^@/);
        my $last = 0;
        my $acc_ctg_len = 0;
        my $acc_ref_len = 0;
        my $match_len = 0;
        my @a = split(/\t/, $_);
        $h->{$a[0]} = 1;
        my ($tag, $num) = cigar::analyze_cigar($a[5]);
        foreach my $i (0 .. scalar(@$tag) - 1){
            if($tag->[$i] =~ /[DI]/ && $num->[$i] > $indel_size || $tag->[$i] =~ /[HS]/ && $num->[$i] > $min_clip){
                $last = 1;
                last;
            }
            else{
                if($tag->[$i] =~ /M/){
                    $match_len += $num->[$i];
                }
                if($tag->[$i] =~ /[HSMI]/){
                    $acc_ctg_len += $num->[$i];
                }
                if($tag->[$i] =~ /[HSMD]/){
                    $acc_ref_len += $num->[$i];
                }
            }
        }
        if($last != 1 && $match_len/$acc_ctg_len > $per && $match_len/$acc_ref_len > $per){
            $h_good->{$a[0]} = 1;
        }
    }
    close sam_fh;
    print join("\t", "#well_aligned", keys %$h_good) . "\n";
    print "#poor_aligned";
    foreach my $key (keys %$h){
        if(!defined $h_good->{$key}){
            print "\t" . $key;
        }
    }
    print "\n";
}


# classify reads in a sam file into well-aligned and poor-aligned 
sub classify_aln{
    my ($sam, $per, $min_clip, $indel_size) = @_;
    my $h_good;
    my $h;
    open sam_fh, "<$sam" or die $!;
    while(<sam_fh>){
        next if($_ =~ /^@/);
        my $last = 0;
        my $acc_len = 0;
        my $soft_clipped = 0;
        my @a = split(/\t/, $_);
        $h->{$a[0]} = 1;
        my ($tag, $num) = cigar::analyze_cigar($a[5]);
        foreach my $i (0 .. scalar(@$tag) - 1){
            if($tag->[$i] =~ /[DI]/ && $num->[$i] > $indel_size || $tag->[$i] =~ /[HS]/ && $num->[$i] > $min_clip){
                $last = 1;
                last;
            }
            else{
                if($tag->[$i] =~ /[HS]/){
                    $soft_clipped += $num->[$i];
                }
                if($tag->[$i] =~ /[HSMI]/){
                    $acc_len += $num->[$i];
                }
            }
        }
        if($last != 1 && 1 - $soft_clipped/$acc_len > $per){
            $h_good->{$a[0]} = 1;
        }
    }
    close sam_fh;
    print join("\t", "#well_aligned", keys %$h_good) . "\n";
    print "#poor_aligned";
    foreach my $key (keys %$h){
        if(!defined $h_good->{$key}){
            print "\t" . $key;
        }
    }
    print "\n";
}

# This returns a hash table with the region of the reads as values, read names as the keys given a sam file, starting from M and end with M
sub get_coordinates{
    my ($sam, $tag) = @_;
    my $h;
    open fh_, "<$sam" or die $!;
    while(<fh_>){
        next if ($_ =~ /^@/);
        my $str = &get_coordinate($_, $tag);
        my @a = split(/\t/, $_);
        if(!defined $h->{$a[0]}){
            $h->{$a[0]} = $str;
        }
        else{
            $h->{$a[0]} .= ";$str";
        }
    }
    close fh_;
    return $h;
}

# given a line of sam, get the region in a string
sub get_coordinate{
    my ($line, $tag) = @_;
    my @a = split(/\t/, $line);
    return "NA" if($a[2] eq "*" || $a[2] eq "*" || $a[5] eq "*");
    my ($matchlen, $totallen);
    if(defined $tag && $tag =~ /match/){
        ($matchlen, $totallen) = split(/\t/, cigar::get_matchlen_on_read_of_all($a[5]));
    }
    my $matchlen_ref = cigar::get_matchlen_on_ref($a[5]);
    my $e = $a[3] + $matchlen_ref;
    my $str = "$a[2]:$a[3]-$e";
    if(defined $tag && $tag =~ /match/){
        $str .= ",$matchlen/$totallen";
    }
    return $str;
}


1;
