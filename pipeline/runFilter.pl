#!/usr/bin/perl
use warnings;
use strict;
require cigar;

if($#ARGV == -1){
    die "Function: filter/pair up alignments from blasr. Starting 08252015, now accepts blasr alignment result in sam format. Output will be in the same format as the input. \nUsage: $0 <output_dir> <suffix_of_blasr> <perc_t> <alignLen_t> <insert_t> <short_prefix>\n";
}

my ($output_dir, $suffix, $perc_t, $alignLen_t, $insert_t, $short_prefix) = @ARGV;
my $h;
my $h_il;
my $h_pb;
my $out = $output_dir . "/out.$suffix";
open OUT, ">$out" or die $!;
my $mode;
if($suffix =~ /m4/){
    # first time get all those have length > t in blasr_1
    $mode = "all";
    $h = &readAln_m4("$output_dir/blasr_1.$suffix");
    &checkAln_m4($h, "$output_dir/blasr_2.$suffix", $mode);
    $h = {};
    # second time get all those have length > t in blasr_2, and those have length < t in blasr_1, to get rid of both > t, which is duplicated in the first time
    $mode = "exc";
    $h = &readAln_m4("$output_dir/blasr_2.$suffix");
    &checkAln_m4($h, "$output_dir/blasr_1.$suffix", $mode);
}
elsif($suffix =~ /sam/){
    $mode = "all";
    $h =  &readAln_sam("$output_dir/blasr_1.$suffix");
    &checkAln_sam($h, "$output_dir/blasr_2.$suffix", $mode);
    $h = {};
    # second time get all those have length > t in blasr_2, and those have length < t in blasr_1, to get rid of both > t, which is duplicated in the first time
    $mode = "exc";
    $h = &readAln_sam("$output_dir/blasr_2.$suffix");
    &checkAln_sam($h, "$output_dir/blasr_1.$suffix", $mode);
}

# print pbs
my $pb_out = $output_dir . "/pb.txt";
open pb_out, ">$pb_out" or die $!;
foreach my $pb (sort keys %$h_pb){
    print pb_out join("\t", $pb, scalar(keys %{$h_pb->{$pb}})) . "\n";
}
print pb_out "###\n";
foreach my $il (keys %$h_il){
    print pb_out join("\t", $il, scalar(keys %{$h_il->{$il}})) . "\n";
}

my $short_all;
if(-e "${short_prefix}_1.fa"){
    $short_all = `grep -c \'\^>\' ${short_prefix}_1.fa`;
}
elsif(-e "${short_prefix}.1.fa"){
    $short_all = `grep -c \'\^>\' ${short_prefix}.1.fa`;
}
chomp $short_all;
#my $PB_all = `grep -c > $PB_fa`;
#chomp $PB_all;
my $short_improved = scalar(keys %$h_il);
my $PB_involved = scalar(keys %$h_pb);
print pb_out "There are in all $short_all short read pair, of which there are $short_improved that have better alignment if to PB. There are $PB_involved PB involved\n";
close pb_out;
1;
sub checkAln_sam{
    my ($h, $f, $mode) = @_;
    open SAM, "<$f" or die $!;
    while(<SAM>){
        next if($_ =~ /^@/);
        my ($rn_q, $flag, $chr, $pos, $mapq, $cigar, ) = split(/\s+/, $_);
        # deal with readname, _1 and _2 only in the function necessarily
        my @tmp = split(/\//, $rn_q);
        #rname is the name without blasr adding
        my $rname = $rn_q;
        $rname = join("/", @tmp[0 .. $#tmp - 1]) if($rn_q =~ /\//);
        #new_rname is the name wihtout fastq pair adding
        my $new_rname = $rname;
        if($rname =~ /_2$/){
            $new_rname =~ s/_2$//;
        }
        elsif($rname =~ /_1$/){
            $new_rname =~ s/_1$//;
        }

        my $ori_r = 0;
        if($flag & 0x0010){
            $ori_r = 1;
        }
        my $rn_r = $chr;
        my $s_r = $pos;
        if($cigar =~ /^(\d+)[SH]/){
            $s_r -= $1;
        }

        my ($aln_len, $len, $cut_len) = &analyze_cigar($cigar);
        # blasr cigar string does not input the true cigar, need to get the total read length 
        ($len) = ($_ =~ /XQ:i:(\d+)/);
        my $a_per = ($aln_len * 100 )/$len;

        # add filters
        next if($a_per < $perc_t);
        #next if($e_q - $s_q < $alignLen_t);
        if($ori_r == 1){
            # reverse, record the end
            $s_r += $len;
        }
        # get the orientation of the mate
        my $ori = 1;
        # when ori_r is 1, this query is reverse strand to the ref
        if($ori_r == 1){
            # require the mate to be forward
            $ori = 0;
        }
        my $tmp = $rn_r . ":" . $ori;
        $new_rname .= "_1" if($rname =~ /_2$/);
        if(defined $h->{$new_rname}->{$tmp}){
            foreach my $p (keys %{$h->{$new_rname}->{$tmp}}){

                next if(abs($s_r - $p) > $insert_t);
                my $l1 = $h->{$new_rname}->{$tmp}->{$p};
                # require at least one having mapped length > min
                # exclusively only for those with alignment length smaller than the threshold
                if($mode eq "exc"){
                    next if($len - $cut_len >= $alignLen_t);
                }
                #next if($x[6] - $x[5] < $alignLen_t && $e_q - $s_q < $alignLen_t);
                if($ori == 0 && $s_r > $p){
                    print OUT &change_rn($l1, $_, $mode);
                    $h_pb->{$rn_r}->{$rname} = 1;
                    $h_il->{$rname}->{$rn_r} = 1;
                }
                elsif($ori == 1 && $s_r < $p){
                    print OUT &change_rn($l1, $_, $mode);
                    $h_pb->{$rn_r}->{$rname} = 1;
                    $h_il->{$rname}->{$rn_r} = 1;
                }
            }
        }
    }
    close SAM;
}
sub checkAln_m4{
    my ($h, $f, $mode) = @_;
    open M4, "<$f" or die $!;
    while(<M4>){
        next if($_ =~ /^@/);
        # new m4 format
        #my ($rn_q, $rn_r, $ori_q, $ori_r, $as, $a_per, $s_r, $e_r, $len_r, $s_q, $e_q, $len_q, ) = split(/\s+/, $_);
        my ($rn_q, $rn_r, $as, $a_per, $ori_q, $s_q, $e_q, $len_q, $ori_r, $s_r, $e_r, $len_r, ) = split(/\s+/, $_); 
        # deal with readname, _1 and _2 only in the function necessarily
        my @tmp = split(/\//, $rn_q);
        #rname is the name without blasr adding
        my $rname = $rn_q;
        $rname = join("/", @tmp[0 .. $#tmp - 1]) if($rn_q =~ /\//);
        #new_rname is the name wihtout fastq pair adding
        my $new_rname = $rname;
        if($rname =~ /_2$/){
            $new_rname =~ s/_2$//;
        }
        elsif($rname =~ /_1$/){
            $new_rname =~ s/_1$//;
        }
        # add filters
        next if($a_per < $perc_t);
        #next if($e_q - $s_q < $alignLen_t);
        if($ori_r == 1){
            # reverse, the pos is reference len - reference end
            $s_r = $len_r - $e_r;
        }
        # get the orientation of the mate
        my $ori = 1;
        # when ori_r is 1, this query is reverse strand to the ref
        if($ori_r == 1){
            # require the mate to be forward
            $ori = 0;
        }
        my $tmp = $rn_r . ":" . $ori;
        $new_rname .= "_1" if($rname =~ /_2$/);
        if(defined $h->{$new_rname}->{$tmp}){
            foreach my $p (keys %{$h->{$new_rname}->{$tmp}}){

                next if(abs($s_r - $p) > $insert_t);
                my $l1 = $h->{$new_rname}->{$tmp}->{$p};
                # require at least one having mapped length > min
                # exclusively only for those with alignment length smaller than the threshold
                if($mode eq "exc"){
                    next if($e_q - $s_q >= $alignLen_t);
                }
                #next if($x[6] - $x[5] < $alignLen_t && $e_q - $s_q < $alignLen_t);
                if($ori == 0 && $s_r > $p){
                    print OUT &change_rn($l1, $_, $mode);
                    $h_pb->{$rn_r}->{$rname} = 1;
                    $h_il->{$rname}->{$rn_r} = 1;
                }
                elsif($ori == 1 && $s_r < $p){
                    print OUT &change_rn($l1, $_, $mode);
                    $h_pb->{$rn_r}->{$rname} = 1;
                    $h_il->{$rname}->{$rn_r} = 1;
                }
            }
        }
    }
    close M4;
}

# change the read name so that the first and seoncd read can be distinguished
sub change_rn{
    my ($line1, $line2, $mode) = @_;
    chomp $line1;
    chomp $line2;
    my @a = split(/\s+/, $line1);
    my @b = split(/\s+/, $line2);
    return join(" ", $a[0] . "/1", @a[1 .. $#a]) . "\n", join(" ", $b[0] . "/2", @b[1 .. $#b]) . "\n" if($mode eq "all");
    return join(" ", $a[0] . "/2", @a[1 .. $#a]) . "\n", join(" ", $b[0] . "/1", @b[1 .. $#b]) . "\n" if($mode eq "exc");
}

sub analyze_cigar{
    my ($cigar) = @_;
    my ($tag, $num) = cigar::analyze_cigar($cigar);
    my $aln_len = 0;
    my $len = 0;
    my $cut_len = 0;
    for(my $i = 0; $i < scalar(@$tag); $i ++){
        if($tag->[$i] eq "M"){
            $aln_len +=$num->[$i];
        }
        if($tag->[$i] =~ /[MISH]/){
            $len += $num->[$i];
        }
        if($tag->[$i] =~ /[SH]/){
            $cut_len += $num->[$i];
        }
    }
    return ($aln_len, $len, $cut_len);
}

# read the alignment of the first fq to a hash for comparison
sub readAln_sam{
    my ($f) = @_;
    my $h;
    open SAM, "<$f" or die $!;
    while(<SAM>){
        next if($_ =~ /^@/);
        my $ori_r = 0;
        my ($rn_q, $flag, $chr, $pos, $mapq, $cigar, ) = split(/\s+/, $_);
        my @tmp = split(/\//, $rn_q);
        $rn_q = join("/", @tmp[0 .. $#tmp - 1]) if($rn_q =~ /\//);

        if($flag & 0x0010){
            $ori_r = 1;
        }
        my $rn_r = $chr;
        my $s_r = $pos;
        if($cigar =~ /^(\d+)[SH]/){
            $s_r -= $1;
        }

        my ($aln_len, $len, $cut_len) = &analyze_cigar($cigar);
        # blasr cigar string does not input the true cigar, need to get the total read length 
        ($len) = ($_ =~ /XQ:i:(\d+)/);
        my $a_per = ($aln_len*100)/$len;
        next if($a_per < $perc_t);
        next if($len - $cut_len < $alignLen_t);
        $h->{$rn_q}->{"$rn_r:$ori_r"}->{$s_r} = $_;
    }
    close SAM;
    return $h;
}


# read the alignment of the first fq to a hash for comparison
sub readAln_m4{
    my ($f) = @_;
    my $h;
    open M4, "<$f" or die $!;
    while(<M4>){
        next if($_ =~ /^@/);
        my ($rn_q, $rn_r, $as, $a_per, $ori_q, $s_q, $e_q, $len_q, $ori_r, $s_r, $e_r, $len_r, ) = split(/\s+/, $_);
        my @tmp = split(/\//, $rn_q);
        $rn_q = join("/", @tmp[0 .. $#tmp - 1]) if($rn_q =~ /\//);
        next if($a_per < $perc_t);
        # as one_mapped_c, check the insert size and one of the mapped length
        next if($e_q - $s_q < $alignLen_t);

        if($ori_r == 1){
            # reverse, the pos is $len - $end on the reference
            $s_r = $len_r - $e_r;
        }
        $h->{$rn_q}->{"$rn_r:$ori_r"}->{$s_r} = $_;
    }
    close M4;
    return $h;
}



