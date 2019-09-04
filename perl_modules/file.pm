use warnings;
use strict;
package file;
require region;
require bam;
# This is a package dealing with files.
#

sub get_max{
    my ($prefix, $threshold, $out) = @_;
    open fh_, ">$out" or die $!;
    open fh_in, "<$prefix.all";
    while(<fh_in>){
        chomp;
        my @aa = split(/\t/, $_);
        if($aa[2] > $threshold){
            print fh_ $aa[0] . "\n";
        }
    }
    close fh_in;
    close fh_;
}

# for multiple groups in a file
# Given the prefix of two files (out.all: containing the stats of all groups; out.txt: containing the ILs and PBs corresponding to a group, with the group number appearing before all following #)
sub generate_groups_file{
    # this is a fasta version by counting the lines
    my ($prefix, $group_file, $out) = @_;
    my $h;
    open fh_in, "<$group_file" or die $!;
    while(<fh_in>){
        chomp;
        $h->{$_} = 1;
    }
    close fh_in;
    if(-e $out){
        `rm $out`;
    }
    my $count_f = "$prefix.all";
    my $read_f = "$prefix.txt";
    my $line_ext = 1;
    my $last = 1;
    my $begin;
    open CF, "<$count_f" or die $!;
    while(<CF>){
        chomp;
        my @a = split(/\t/, $_);
        $begin = $last + 1;
        # duration
        $line_ext = $a[1] + $a[2];
        # last line of this block
        $last = $begin + $line_ext;
        if(defined $h->{$a[0]}){
            my $last1 = $last - 1;
            `sed -n $begin,${last1}p $read_f >> $out`;
            print "Successfully extracted block starting $begin, with $line_ext lines to $out\n";
        }
    }
    close CF;
    return;
}


# Given the prefix of two files (out.all: containing the stats of all groups; out.txt: containing the ILs and PBs corresponding to a group, with the group number appearing before all following #)
sub generate_group_file{
    # this is a fasta version by counting the lines
    my ($prefix, $group, $out) = @_;
    die "Group not an integer number\n" if($group !~ /^\d+$/);
    my $count_f = "$prefix.all";
    my $read_f = "$prefix.txt";
    my $line = 1;
    my $line_ext = 1;
    open CF, "<$count_f" or die $!;
    while(<CF>){
        chomp;
        my @a = split(/\t/, $_);
        if($a[0] ne $group){
            # not yet, add up the lines
            $line += $a[1] + $a[2] + 1;
        }
        else{
            # this is it, record the lines to be extracted and jump out of while
            $line_ext = $a[1] + $a[2];
            last;
        }
    }
    close CF;
    my $last = $line + $line_ext;
    `sed -n $line,${last}p $read_f > $out`;
    print "Successfully extracted block starting $line, with $line_ext lines to $out\n";
    return;
}

# This checks the number of split reads supporting the breakpoint on the two bam files
sub check_ground_truth_split{
    my ($bed, $bam_all, $bam_spe, $flanking, $off, $out_bed) = @_;
    open out_bed, ">$out_bed" or die $!;
    open BED, "<$bed" or die $!;
    while(<BED>){
        chomp;
        my @a = split(/\t/, $_);
        my $region = "$a[0]:" . ($a[1] - $flanking) . "-" . ($a[1] + $flanking);
        my ($r1, $br1) = bam::extract_soft_clip($bam_all, $region, $a[1], "-", $off);
        $DB::single= 1;
        my ($s1, $bs1) = bam::extract_soft_clip($bam_spe, $region, $a[1], "-", $off);
        $region = "$a[0]:" . ($a[2] - $flanking) . "-" . ($a[2] + $flanking);
        my ($r2, $br2) = bam::extract_soft_clip($bam_all, $region, $a[2], "+", $off);
        my ($s2, $bs2) = bam::extract_soft_clip($bam_spe, $region, $a[2], "+", $off);
        my $r = $r1 + $r2;
        my $s = $s1 + $s2;
        print out_bed join("\t", @a, $r, $s, $r1, $br1, $r2, $br2, $s1, $bs1, $s2, $bs2) . "\n";
    }
    close BED;
    close out_bed;
}

# This checks the number of reads on the two bam files, and cat them till the end of the call 
sub check_ground_truth_quick{
    my ($bed, $bam_all, $bam_spe, $flanking, $readlen, $out_bed) = @_;
    my $r = $flanking;
    open out_bed, ">$out_bed" or die $!;
    open BED, "<$bed" or die $!;
    while(<BED>){
        chomp;
        my @a = split(/\t/, $_);
        my $key = "$a[0]:$a[1]-$a[2]";
        my $h = 0;
        my $l = 0;
        if($a[2] - $a[1] > 2 * $readlen){
            my $region = "$a[0]:" . ($a[1] - $r) . "-" . ($a[1] + $r);
            $h = bam::samtools_count($bam_all, $region);
            $l = bam::samtools_count($bam_spe, $region);

            $region = "$a[0]:" . ($a[2] - $r) . "-" . ($a[2] + $r);
            $h += bam::samtools_count($bam_all, $region);
            $l += bam::samtools_count($bam_spe, $region);
        }
        else{
            my $region = "$a[0]:" . ($a[1] - $r) . "-" . ($a[2] + $r);
            $h += bam::samtools_count($bam_all, $region);
            $l += bam::samtools_count($bam_spe, $region);
        }
        print out_bed join("\t", @a, $h, $l) . "\n";
    }
    close BED;
    close out_bed;
}


# This takes a bed file, for each row, take the short read id from the bam file, then check whether at least three read ids appear in one cluster, by checking a big file containing the names of the text files, which contains an IL ID per row that belongs to this cluster. If yes, further check the number of PBs in this cluster, and whether it has assembled contig, and how many of the stated ILs can be mapped normally to the assembled contig, whether it calls an SV and whether it calls it correclty
# Input: big_txt, a big text each row containing a text file, including the path
# bed, gold standard to be checked
# bam, the bam file to be checked
# flanking and readlen is for overlapping metrics
# min_IL is for a minimum bound to take a cluster to belong to this call
# my_call_bed is a bed file that I called for comparison (deleted)
sub combined_stat{
    my ($big_txt, $bed, $bam, $flanking, $readlen, $min_IL) = @_;
    # first read the big txt and record all the ILs corresponding to the cluster
    my $r = $flanking;
    my $h;
    my $h_pb;
    # a hash containing the files of the clusters, with keys the cluster number only
    my $h_files;
    open big_txt, "<$big_txt" or die $!;
    $DB::single = 1;
    while(<big_txt>){
        chomp;
        my @F = split(/\t/, $_);
        next if(! -e $F[0]);
    $DB::single = 1;
        open c_fh, "<$F[0]" or die $!;
        my ($c) = ($F[0] =~ /out\.g(\d+)\./);
        $h_files->{$c} = $F[0];
        # also count the number of pbs belonging to this cluster
        $h_pb->{$c} = 0;
        while(<c_fh>){
            chomp;
            if($_ !~ /\./){
                $h->{$_} = $c;
            }
            else{
                $h_pb->{$c} ++;
            }
        }
        close c_fh;

    }
    close big_txt;

    # second foreach call in bed, take the ILs from bam, check them in $h
    open BED, "<$bed" or die $!;
    while(<BED>){
        chomp;
        my @a = split(/\t/, $_);
        my $key = "$a[0]:$a[1]-$a[2]";
        # an array to record all the read IDs
        my @b = ();
        # a hash to record the cluster
        my $h_tmp;
        # cluster to be selected and the number of the ILs in it
        my $max = 0;
        my $select_cl = -1;
        # a hash to record the readnames, to avoid duplicates of pairs
        my $h_rd;
        if($a[2] - $a[1] > 2 * $readlen){
            my $region = "$a[0]:" . ($a[1] - $r) . "-" . ($a[1] + $r);
            $h_rd = bam::samtools_rdname($bam, $region, $h_rd);

            $region = "$a[0]:" . ($a[2] - $r) . "-" . ($a[2] + $r);
            $h_rd = bam::samtools_rdname($bam, $region, $h_rd);
        }
        else{
            my $region = "$a[0]:" . ($a[1] - $r) . "-" . ($a[2] + $r);
            $h_rd = bam::samtools_rdname($bam, $region, $h_rd);
        }
        $DB::single = 1;
        foreach my $b_ (keys %$h_rd){
            if(defined $h->{$b_}){
                my $cl = $h->{$b_};
                if(!defined $h_tmp->{$cl}){
                    $h_tmp->{$cl} = 1; 
                }
                else{
                    $h_tmp->{$cl} ++;
                }
                if($h_tmp->{$cl} > $max){
                    $max = $h_tmp->{$cl};
                    $select_cl = $cl;
                }
            }
        }
        my $IL = 0;
        my $PB = 0;
        my $AS = 0;
        my $SV = 0;
        my $correct_SV = 0;
        # check if the # of IL satisfies the minimum requirement
        if($max >= $min_IL){
            $DB::single = 1;
            $IL = $max;
            $PB = $h_pb->{$select_cl};
            $AS = &check_assembled($h_files, $select_cl);
            #($SV, $correct_SV) = &check_SV($h_files, $select_cl, $my_call, $key);
        }

        print join("\t", @a, $select_cl, $IL, $PB, $AS, $SV, $correct_SV)."\n";

    } # end of bed file
    close BED;
}


# given a dir/list of error files, each line in the form of "Cluster $num does not assemble successfully within $time seconds.", and an sh file containing all the commands running assembly, each line in the form of "~/c/prepare_reads_spaceEfficient.sh ~/c/analyze_all_3/mapped/rerun_09092015 $num ~/c/analyze_all_3/mapped/rerun_09092015/pipe/union_find/out $thread_num /scratch/bcb/xfan3/c/analyze_all_3/mapped/rerun_09092015/pipe/assembly/$output_prefix $time", output the same command on all those numbers appearing in the error file, with the following modifications: 1) multiple the time by $times, 2) new thread_num, 3) new output_prefix. Put all output in one file. Note: those $num already in $bash_files_exclude (can be multiple files, separated by ;) will not be in the output. $sep is the number of batch before a separator # is inserted.
sub run_assembly_w_moretime{
    my ($dir_errs, $bash_file, $bash_files_exclude, $times, $thread_num, $output_prefix, $sep) = @_;
    my $h_exclude;
    # read the $num in bash_files_exclude so that these num will not be added for rerun
    my @f_exclude = split(/;/, $bash_files_exclude);
    foreach (@f_exclude){
        if(-e $_){
            open fh_exc, "<$_" or die $!;
            while(<fh_exc>){
                chomp;
                next if($_ =~ /^#/);
                my @a = split(/\s+/, $_);
                $h_exclude->{$a[2]} = 1;
            }
            close fh_exc;
        }
    }

    # read the template of bash command and record what are in the bash by $num
    # hash recording the bash command from the first line in bash without # in prefix
    my $c;
    my $n = 0;
    # hash recording what are in the bash by $num
    my $h;
    open b_fh, "<$bash_file" or die $!;
    while(<b_fh>){
        chomp;
        next if($_ =~ /^#/);
        my @a = split(/\s+/, $_);
        if($n == 0){
            $c->{cmd} = $a[0];
            $c->{dir} = $a[1];
            $c->{src} = $a[3];
            if($a[5] =~ /\//){
                my @tmp_out = split(/\//, $a[5]);
                $c->{out} = join("/", @tmp_out[0 .. $#tmp_out - 1], $output_prefix);
            }
            else{
                $c->{out} = $output_prefix;
            }
            $c->{tm} = $a[6] * $times;
        }
        if(!defined $h_exclude->{$a[2]}){
            $h->{$a[2]} = 1;
        }
        $n ++;
    }
    close b_fh;

    my @errs = split(/\n/, `ls $dir_errs`);
    # m controls which batch it is in, output prefix changes correspondingly
    my $m = 0;
    $n = 0;
    foreach my $err_f (@errs){
        open e_fh, "<$err_f" or die $!;
        while(<e_fh>){
            if($_ =~ /^Cluster (\d+) does not/){
                my $num = $1;
                if(defined $h->{$num}){
                    print join(" ", $c->{cmd}, $c->{dir}, $num, $c->{src}, $thread_num, $c->{out}.".$m", $c->{tm}) . "\n";
                    $n ++;
                    if($n >= $sep){
                        print "#\n";
                        $m ++;
                        $n = 0;
                    }
                }
            }
        }
        close e_fh;
    }
}


# given a hash table, retrieve the specific location of the assembly directory, check if the contig is there and with size nonzero
sub check_assembled{
    my ($h_files, $select_cl) = @_;
    my $f = $h_files->{$select_cl};
    my $str = "asm_g$select_cl/asm_g$select_cl.ctg.fasta";
    if($f =~ /\//){
        my @a = split(/\//, $f);
        $DB::single = 1;
        $f = join("/", @a[0 .. $#a - 1]) . "/$str";
    }
    else{
        $f = "./$str";
    }

    if(-e $f && ! -z $f){
        return 1;
    }
    return 0;
}




# This takes a txt file with the readname of each row, and a sam file. Output a sam file with only those readnames appearing in the txt file.
sub take_reads{
    my ($sam, $txt, $out_sam) = @_;
    # read in the short read ids in the large cluster
    my $h_txt;
    open TXT, "<$txt" or die $!;
    while(<TXT>){
        chomp;
        if($_ !~ /\./){
            $h_txt->{$_} = 1;
        }
    }
    close TXT;

    open OUT, ">$out_sam" or die $!;
    open SAM, "<$sam" or die $!;
    while(<SAM>){
        my @a = split(/\t/, $_);
        last if($a[0] =~ /\./);
        if(defined $h_txt->{$a[0]}){
            print OUT $_;
        }
    }
    close SAM;
    close OUT;
}

# This takes a bed file of GT calls (mainly deletions), a sam file of extracted reads and a text file containing the read names corresponding to one cluster. For each GT call, check if any extracted read falls in this range (overlapping with the breakpoints), if yes, further check if they are in the text file. Output is 1). a sam file containing all the reads overlapping with the breakpoints, with the last column the breakpoints in the format of chr:start-end; 2). the original bed file with extra columns: (1) the # of extracted reads it contains; (2) the # of extracted reads it contains falling into the large cluster.
# flanking is the region allowed for inaccurate breakpoint
sub check_ground_truth{
    my ($bed, $sam, $txt, $flanking, $readlen) = @_;
    # first of all, read the bed and save to a hash
    my $h_bed;
    my $r = $flanking;
    open BED, "<$bed" or die $!;
    while(<BED>){
        chomp;
        my @a = split(/\t/, $_);
        my $key = "$a[0]:$a[1]-$a[2]";
        if($a[2] - $a[1] > 2 * $readlen){
            my $region = "$a[0]:" . ($a[1] - $r) . "-" . ($a[1] + $r);
            push @{$h_bed->{$key}}, $region;
            $region = "$a[0]:" . ($a[2] - $r) . "-" . ($a[2] + $r);
            push @{$h_bed->{$key}}, $region;
        }
        else{
            my $region = "$a[0]:" . ($a[1] - $r) . "-" . ($a[2] + $r);
            push @{$h_bed->{$key}}, $region;
        }
    }
    close BED;

    # read in the short read ids in the large cluster
    my $h_txt;
    open TXT, "<$txt" or die $!;
    while(<TXT>){
        chomp;
        if($_ !~ /\./){
            $h_txt->{$_} = 1;
        }
    }
    close TXT;

    # second, read sam file and analyze read one by one
    open OUT_SAM, ">$sam.GT" or die $!;
    open SAM, "<$sam" or die $!;
    my $h_count;
    my $h_count_cluster;
    while(<SAM>){
        chomp;
        my @a = split(/\t/, $_);
        my ($rn, $flag, $chr, $pos, $qual, $cigar) = @a[0 .. 5];
        if($cigar =~ /^(\d+)[SH]/){
            $pos -= $1;
        }
        my $end = $pos + $readlen;
        my $region = "$chr:$pos-$end";
        foreach my $key (keys %$h_bed){
            my $if_overlap = region::overlap1($region, $h_bed->{$key});
            if($if_overlap){
                print OUT_SAM join("\t", $_, $key) . "\n";
                $h_count->{$key} ++;
                if(defined $h_txt->{$rn}){
                    $h_count_cluster->{$key} ++;
                }
            }
        }
    }
    close SAM;

    # output
    open BED, ">$bed.GTanalyze" or die $!;
    foreach my $key (keys %$h_bed){
        my @a = split(/:/, $key);
        my @b = split(/-/, $a[1]);
        my $h1 = 0;
        my $h2 = 0;
        if(defined $h_count->{$key}){
            $h1 = $h_count->{$key};
        }
        if(defined $h_count_cluster->{$key}){
            $h2 = $h_count_cluster->{$key};
        }
        print BED join("\t", $a[0], @b, $h1, $h2);
    }
    close BED;

}


# check if the ith column of fileA contains all of the jth column of fileB
sub A_contains_B{
    my ($fA, $i, $fB, $j) = @_;
    # read in the column i of fA
    my $hA = &read_col_i($fA, $i);
    # check column j of fB
    my @missed = ();
    open fh_, "<$fB" or die $!;
    while(<fh_>){
        chomp;
        my @a = split(/\s+/, $_);
        if(!defined $hA->{$a[$j]}){
            push @missed, $a[$j];
        }
    }
    close fh_;
    if(@missed == 0){
        return 1;
    }
    else{
        return join("\t", @missed)."\n";
    }
}

# read in ith column of file and return a hash with keys of it
sub read_col_i{
    my ($f, $i) = @_;
    my $h;
    open fh_, "<$f" or die $!;
    while(<fh_>){
        chomp;
        my @a = split(/\s+/, $_);
        $h->{$a[$i]} = 1;
    }
    close fh_;
    return $h;
}

sub cat_by_adding_name{
    my ($f, $i, $c) = @_;
    # Given a file $f, this function adds $i to $c th column, as a prefix. 
    open fh_, "<$f" or die $!;
    while(<fh_>){
        chomp;
        my @a = split(/\s+/, $_);
        my $s = "$i.$a[$c]";
        my $prefix = "";
        $prefix = join(" ", @a[0 .. $c - 1]) if($c >= 1);
        my $suffix = "";
        $suffix = join(" ", @a[$c + 1 .. $#a]) if($c < $#a);
        if($prefix ne "" && $suffix ne ""){
            print join(" ", $prefix, $s, $suffix);
        }
        elsif($prefix eq ""){
            print join(" ", $s, $suffix);
        }
        elsif($suffix eq ""){
            print join(" ", $prefix, $s);
        }
        print "\n";
    }
    close fh_;
}

# check if each line in fileB contains any line in fileA, with the prefix and suffix strings.
sub compare_two_files{
    my ($fileA, $colA, $fileB, $colB, $prefix, $suffix) = @_;
    my $h;
    open fA, "<$fileA" or die $!;
    while(<fA>){
        chomp;
        my @a = split(/\s+/, $_);
        my $str = $prefix . $a[$colA] . $suffix;
        $h->{$str} = 1;
    }
    close fA;
    open fB, "<$fileB" or die $!;
    while(<fB>){
        my @a = split(/\s+/, $_);
        foreach my $keys (keys %$h){
            if($a[$colB] =~ /$keys/){
                print $_;
                last;
            }
        }
    }
    close fB;
}


1;
