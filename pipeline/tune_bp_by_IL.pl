#!/usr/bin/perl
use warnings;
use strict;
require cigar;

if(scalar(@ARGV) == 0){
    die "Given a list of putative bps, This script checks with IL in the same group with sam file separated by group number, and flanking size within the left and right of each breakpoint, call the accurate breakpoint from IL if there are >= num_sp IL reads soft clipped at the same position. \nUsage: $0 bp_list IL_sam_sep_by_group flanking_size IL_sp_num\n";
}

my ($b1, $IL_sam, $flank, $sp) = @ARGV;
# by file is only for test
# in /scratch/bcb/xfan3/d/NA19238/pipe/analysis_10202015
# perl ~/combinePBIL/scripts/tune_bp_by_IL.pl debug/test debug/test.sam 50 1
my $b_h = &read_bp_by_file($b1);
#my $b_h = &read_bp($b1);
my $b_IL = &infer_bp_from_sam($IL_sam, $b_h, $flank);
my $accurate_bp = &infer_accurate_bp($b_h, $b_IL, $flank);
print @$accurate_bp if(defined $accurate_bp);

1;

# check if a point is within a range
sub within_range{
    my ($p, $r) = @_;
    my ($a, $b) = split(/:/, $r);
    if($p <= $b && $p >= $a){
        return 1;
    }
    else{
        return 0;
    }
}

# given two hash tables, this compares and report the accurate breakpoint if any of b_IL within flanking size of those called by b
sub infer_accurate_bp{
    my ($b_h, $b_IL, $flank) = @_;
    my $return;
    my ($chr, $pos, $SV_size, $SV_type, $PB, $PB_start, $PB_end, $PB_spt_num) = split(/:/, $b_h);
    if($SV_type !~ /[ID]/){
        die "The SV type is not within I or D: $b_h\n";
    }
    my ($l_range, $r_range);
    if($SV_type eq "I"){
        $l_range = ($pos - $flank) . ":" . ($pos + $flank);
        $r_range = $l_range;
    }
    elsif($SV_type eq "D"){
        $l_range = ($pos - $flank) . ":" . ($pos + $flank);
        $r_range = ($pos + $SV_size - $flank) . ":" . ($pos + $SV_size + $flank);
    }

    my $this_IL = $b_IL;
    my $h_x;
    my $h_num;
    $h_num->{L} = -1;
    $h_num->{R} = -1;
    if(defined $this_IL){
        foreach my $key (keys %$this_IL){
            my $num_sp = $this_IL->{$key};
            my @x = split(/:/, $key);
            next if($x[1] ne $chr);
            if($x[0] eq "L"){
                # left soft clip, right breakpoint
                if(&within_range($x[2], $r_range)){
                    if($h_num->{R} < $num_sp){
                        $h_x->{R} = $x[2];
                        $h_num->{R} = $num_sp;
                    }
                }
            }
            elsif($x[0] eq "R"){
                # right soft clip, left breakpoint
                if(&within_range($x[2], $l_range)){
                    if($h_num->{L} < $num_sp){
                        $h_x->{L} = $x[2];
                        $h_num->{L} = $num_sp;
                    }
                }
            }
        }
    }
    my $output;
    my ($pos_from_IL, $spt_IL, $new_pos, $new_size);
    if(defined $h_x->{L}){
        if(defined $h_x->{R}){
            $pos_from_IL = "$h_x->{L}:$h_x->{R}";
            $spt_IL = "$h_num->{L}:$h_num->{R}";
            $new_pos = $h_x->{L};
            if($SV_type eq "I"){
                # take the left 
                $new_size = $SV_size;
            }
            elsif($SV_type eq "D"){
                $new_size = $h_x->{R} - $new_pos;
            }
        }
        else{
            $pos_from_IL = "$h_x->{L}:";
            $spt_IL = "$h_num->{L}:";
            $new_pos = $h_x->{L};
            $new_size = $SV_size;
        }
    }
    elsif(defined $h_x->{R}){
        # assume the size is correct
        if($SV_type eq "I"){
            $new_pos = $h_x->{R};
            $new_size = $SV_size;
        }
        elsif($SV_type eq "D"){
            $new_pos = $h_x->{R} - $SV_size;
            $new_size = $SV_size;
        }
        $pos_from_IL = ":$h_x->{R}";
        $spt_IL = ":$h_num->{R}";
    }
    else{
        # nothing
        $new_pos = $pos;
        $new_size = $SV_size;
        $pos_from_IL = ":";
        $spt_IL = ":";
    }

    $output = join("\t", $chr, $new_pos, $new_size, $SV_type, $PB, $PB_start, $PB_end, $PB_spt_num, $pos_from_IL, $spt_IL, $pos, $SV_size) . "\n";
    push @$return, $output if($new_size > 0 && $SV_type eq "D" || $SV_type eq "I");

    return $return;

}
# read a file, and put them into a hash by group
sub read_bp_by_file{
    my ($b) = @_;
    my $h;
    open b_fh, "<$b" or die $!;
    while(<b_fh>){
        my $l = $_;
        chomp $l;
        my ($chr, $pos, $SV_size, $SV_type, $ctg, $ctg_start, $ctg_end, $ctg_spt_num ) = split(/\t/, $l);
        my @ctgs = split(/\./, $ctg);
        if($ctgs[0] =~ /^g(\d+)/){
            my $g = $1;
            $h->{$g} = join(":", $chr, $pos, $SV_size, $SV_type, $ctg, $ctg_start, $ctg_end, $ctg_spt_num);
        }
    }
    close b_fh;
    return $h;
}


# read an array, and put them into a hash by group
sub read_bp{
    my ($b) = @_;
    my $h;
    foreach my $l (@$b){
        chomp;
        my ($chr, $pos, $SV_size, $SV_type, $ctg, $ctg_start, $ctg_end, $ctg_spt_num ) = split(/\t/, $l);
        my @ctgs = split(/\./, $ctg);
        if($ctgs[0] =~ /^g(\d+)/){
            my $g = $1;
            $h->{$g} = join(":", $chr, $pos, $SV_size, $SV_type, $ctg, $ctg_start, $ctg_end, $ctg_spt_num);
        }
    }
    return $h;
}

# This reads a sam file containing IL reads separated by group number, infer breakpoint if there are at least num of reads supporting it, save results to a hash table with keys the group
sub infer_bp_from_sam{
    my ($IL_sam, $b_h, $flank) = @_;
    my $h_t;
    my $g = "NA";
    open SAM, "<$IL_sam" or die $!;
    while(<SAM>){
        if($_ =~ /^#(\d+)/){
            # do comparison for every cluster
            if(defined $b_h->{$g}){
                my $accurate_bp = &infer_accurate_bp($b_h->{$g}, $h_t, $flank);
                print @$accurate_bp if(defined $accurate_bp);
            }
            $g = $1;
            $h_t = {};
        }
        else{
            my @a1 = split(/\t/, $_);
            my $cigar = $a1[5];
            my $str = "";
            if($cigar =~ /S/){
                if($cigar =~ /^\d+S/){
                    $str = "L:$a1[2]:$a1[3]";
                }
                if($cigar =~ /\d+S$/){
                    my ($tag, $num) = cigar::analyze_cigar($cigar);
                    my $pos = $a1[3];
                    for(my $i = 0; $i < scalar(@$tag); $i ++){
                        if($tag->[$i] =~ /[MD]/){
                            $pos += $num->[$i];
                        }
                    }
                    $str = "R:$a1[2]:$pos";
                }
                if(defined $h_t->{$str}){
                    $h_t->{$str} ++;
                }
                else{
                    $h_t->{$str} = 1;
                }
            }
        }
    }
    close SAM;
}


