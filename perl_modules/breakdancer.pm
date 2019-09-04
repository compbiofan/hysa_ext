use warnings;
use strict;

package breakdancer;
require vcf;
require bam;
require overlap;

# Function: This package manipulates the file format of breakdancer. 

# change the type of brekadancer according to string, and write to the same file
# add header
sub format_breakdancer{
    my ($f, $str, $header_str) = @_;
    my $h = &change_type($str);
    my $output = $f . ".tmp";
    open OUT, ">$output" or die $!;
    my $tag = 0;
    open FH, "<$f" or die $!;
    while(<FH>){
        if($_ !~ /^#/ && $tag == 0){
            print OUT join("\t", split(/:/, $header_str))."\n"  if(defined $header_str);
            $tag = 1;
        }
        else{
            my @a = split(/\t/, $_);
            $a[6] = $h->{$a[6]} if(defined $h->{$a[6]});
            #$a[6] = "CTX" if($a[0] ne $a[3]);
            print OUT join("\t", @a);
        }
    }
    close FH;
    `mv $output $f`;
}

sub read_breakdancer{
    # It reads a breakdancer formatted file and record it to a hash in a form of chr1:pos1:chr2:pos2:type -> line
    my ($f, $mode, $spec_hash) = @_;
    my $hash;

    open fh_, "<$f" or die $!;
    while(<fh_>){
        next if($_ =~ /^#/);
        chomp;
        my @a = split(/\s+/, $_);
        if(!defined $mode || $mode eq ""){
            $hash->{join(":", @a[0 .. 1], @a[3 .. 4])} = join(":", @a);
        }
        elsif($mode eq "sort_by_all"){
            $hash->{$_} = 1;
        }
        elsif($mode eq "CTX_overlap"){
            $DB::single = 1;
            my $str = join(":", @a[0 .. 1], @a[3 .. 4]);
            $hash->{join(":", @a[0 .. 1])} = $str;
            $hash->{join(":", @a[3 .. 4])} = $str;
        }
        elsif($mode eq "ITX_overlap"){
            next if($a[4] < $a[1]);
            $hash->{$a[0]}->{"$a[1].$a[4]"} = $a[6];
        }
        elsif($mode eq "SPEC"){
            my $str = &print_spec(\@a, $spec_hash);
            if($a[0] eq $a[3]){
                next if($a[4] < $a[1]);
                $hash->{ITX}->{$a[0]}->{"$a[1].$a[4]"} = $str;
            }
            else{
                $hash->{CTX}->{join(":", @a[0 .. 1])} = $str;
                $hash->{CTX}->{join(":", @a[3 .. 4])} = $str;
            }
        }
    }
    close fh_;

    return $hash;
}

# print a string according to the spec, if col to be print, it is in the key of spec, otherwise not, and data in @a. Output format is col_id:entry:(P|N);col_id:entry:(P|N)
sub print_spec{
    my ($array, $spec_hash) = @_;
    my @strs;
    for(my $i = 0; $i < scalar(@$array); $i ++){
        my $tag = "N";
        $tag = "P" if(defined $spec_hash->{$i});
        push @strs, join(":", $i, $array->[$i], $tag);
    }
    return join(";", @strs);
}

sub write_breakdancer{
    # it reads a hash with $str (chr1:start:chr2:end) -> @samples and write to a breakdancer file
    my ($hash, $sv) = @_;
    open SV, ">$sv" or die $!;
    foreach my $key (keys %$hash){
        my @a = split(/:/, $key);
        # assuming ITX for all
        print SV join("\t", @a[0 .. 1], "+", @a[2 .. 3], "-", "ITX", $a[3] - $a[1], "99", "10");
        foreach my $s (@{$hash->{$key}}){
            print SV "\t" . join("|", $s, 5);
        }
        print SV "\n";
    }
    close SV;
}

# convert a string connected by ":" into a hash
sub analyze_spec{
    my ($spec) = @_;
    my @a = split(/:/, $spec);
    my $hash;
    foreach (@a){
        $hash->{$_} = 1;
    }
    return $hash;
}

# change type; inputs are old:new;old:new; output types are in the form of $hash->{$old} = $new
sub change_type{
    my $h_str = shift;
    my $h;
    my @a = split(/;/, $h_str);
    $DB::single = 1;
    foreach (@a){
        my ($old, $new) = split(/:/, $_);
        $h->{$old} = $new;
    }
    return $h;
}


# this script separate CTX and non-CTX in fileA and fileB in breakdancer format, and add specs for both of them to be a hash, then run overlap or somatic by overlap.pm. spec is in the format of col1:col2 (all those appear will be printed)
# An example is in ~/TIGRA-ext/Stephens*/run_overlap.pl
sub run_overlap_spec{
    my ($fileA, $fileB, $specA, $specB, $flank) = @_;
    my $mode_breakdancer = "SPEC";
    my $specA_hash = &analyze_spec($specA);
    my $specB_hash = &analyze_spec($specB);
    my $hashA = &read_breakdancer($fileA, $mode_breakdancer, $specA_hash);
    my $hashB = &read_breakdancer($fileB, $mode_breakdancer, $specB_hash);
    overlap::run_ctx($hashA->{CTX}, $hashB->{CTX}, $flank);
    #overlap::run_ctx_somatic($hashA->{CTX}, $hashB->{CTX}, $flank) if($mode =~ /somatic/i);
    #$DB::single = 1; 
    overlap::run_itx($hashA->{ITX}, $hashB->{ITX});
    #overlap::run_itx_somatic($hashA->{ITX}, $hashB->{ITX}) if($mode =~ /somatic/i);

}

sub comp_bd_output{
    # it compares two breakdancer files and output those entries of A appearing B, and summarize the breakpoint resolution histogram
    my ($fileA, $fileB, $num_shift) = @_;
    my $bd_hash = &read_breakdancer($fileB, "");
    my $start;
    my $end;
    my $chr2;
    my $hash;
    my $hash_str;
    open CM_FH, "$fileA" or die $!;
    while(<CM_FH>){
        next if($_ =~ /^#/);
        my @a = split(/\t/, $_);
        $start = $a[1];
        $end = $a[4];
        $chr2 = $a[2];
        $chr2 = $a[3];
        my $break = 0;
        # allow for a subtle difference
        foreach my $i (-$num_shift .. $num_shift){
            foreach my $j (-$num_shift .. $num_shift){
                my $key1 = join(":", ($a[0], $start + $i, $chr2, $end + $j));
                my $key2 = join(":", ($chr2, $end + $j, $a[0], $start + $i));
                if(defined $bd_hash->{$key1} || defined $bd_hash->{$key2}){
                    $hash->{$i} ++;
                    $hash->{$j} ++;
                    my $str = $_;
                    $hash_str->{$str} = 1;
                    $break = 1;
                    last;
                }
            }
            last if($break == 1);
        }
    }
    close CM_FH;
    foreach my $key (sort keys %$hash_str){
        print $key;
    }
    foreach my $i (sort {$a <=> $b} keys %$hash){
        print join("\t", $i, $hash->{$i}) . "\n";
    }
}



sub read_breakdancer_for_idup{
    # store the two breakpoints separately in their orientations, and the corresponding samples and the other end
    my ($f) = @_;
    $DB::single = 1;
    my $hash;
    open fh_, "<$f" or die $!;
    while(<fh_>){
        next if($_ =~ /^#/);
        chomp;
        my @a = split(/\s+/, $_);

        my $ori1 = &analyze_ori($a[2]);
        my $ori2 = &analyze_ori($a[5]);
        next if($ori1 eq "" || $ori2 eq "");
        my $hash_sup_sample = &read_sup_sample($_);
        $hash->{$ori1}->{$a[0]}->{$a[1]}->{samples} = $hash_sup_sample;
        $hash->{$ori1}->{$a[0]}->{$a[1]}->{pair} = "$a[3]:$a[4]";
        $hash->{$ori1}->{$a[0]}->{$a[1]}->{pair_ori} = $ori2;
        $hash->{$ori2}->{$a[3]}->{$a[4]}->{samples} = $hash_sup_sample;
        $hash->{$ori2}->{$a[3]}->{$a[4]}->{pair} = "$a[0]:$a[1]";
        $hash->{$ori2}->{$a[3]}->{$a[4]}->{pair_ori} = $ori1;
    }
    close fh_;
    return $hash;
}

sub analyze_ori{
    # given a+b- ori, determine which is the majority
    my $ori = shift;
    my ($x, $y);
    if($ori =~ /^(\d+)\+(\d+)\-$/){
        ($x, $y) = ($1, $2);
    }
    elsif($ori =~ /\+/ && $ori !~ /\-/){
        return '+';
    }
    elsif($ori =~ /\-/ && $ori !~ /\+/){
        return '-';
    }

    if($x eq $y){
        return "";
    }
    elsif($x > $y){
        return '+';
    }
    else{
        return '-';
    }
}

# for 1000G to distinguish config files into populations, one pop per cfg file

sub group_cfg_by_pop{
    my ($sample_list, $dir, $pop_list) = @_;
    # read population list file, one per line
    my $pop_hash;
    open POP, "<$pop_list" or die $!;
    while(<POP>){
        chomp;
        $pop_hash->{$_} = 1;
    };
    close POP;

    # read sample list file, in the format of SAMPLE:LOCATION, in which the full bam location contains the population name
    my $hash;
    open SAMPLE, "<$sample_list" or die $!;
    while(<SAMPLE>){
        chomp;
        my @a = split(/:/, $_);
        foreach my $pop (keys %$pop_hash){
            if($a[1] =~ /\.$pop\./){
                if(-e "$dir/$a[0].cfg"){
                    if(!defined $hash->{$a[0]}){
                        $hash->{$a[0]} = "";
                    }
                    open CFG, "<$dir/$a[0].cfg" or die $!;
                    while(<CFG>){
                        $hash->{$pop} .= $_;
                    }
                    close CFG;
                }
                last;
            }
        }
    }
    foreach my $pop (keys %$hash){
        open CFG, ">$dir/$pop.cfg" or die $!;
        print CFG $hash->{$pop};
        close CFG;
    }

}

# for running breakdancer-max by chromosome
sub run_breakdancer_by_chr{
    # cfg_name could be a population name, such as YRI.cfg with absolute path, or a command such as POP=`head -n ${PBS_ARRAYID} ../population.lst | tail -n 1`
    my $cfg_name = shift;
    my $cfg_dir = "";
    my @chr = 1 .. 22;
    push @chr, 'X';
    push @chr, 'Y';
    foreach (@chr){
        if($cfg_name =~ /head.+tail/){
            # multiple array jobs
            print $cfg_name . "\n";
            $cfg_dir = shift if($cfg_dir eq "");
            print "/scratch/xfan3/pkg/breakdancer/breakdancer-max -o $_ $cfg_dir/\${POP}.cfg > \${POP}.chr$_.sv\n";
        }
        else{
            print "/scratch/xfan3/pkg/breakdancer/breakdancer-max -o $_ $cfg_name > $cfg_name.chr$_.sv\n";
        }
        print "#\n";
    }
}

# for ctx detection, similar with by_chr but without multiple chromosomes
sub run_breakdancer_ctx{
    # cfg_name could be a population name, such as YRI.cfg with absolute path, or a command such as POP=`head -n ${PBS_ARRAYID} ../population.lst | tail -n 1`
    my $cfg_name = shift;
    my $cfg_dir = "";
    if($cfg_name =~ /head.+tail/){
        # multiple array jobs
        print $cfg_name . "\n";
        $cfg_dir = shift if($cfg_dir eq "");
        print "/scratch/xfan3/pkg/breakdancer/breakdancer-max -t $cfg_dir/\${POP}.cfg > \${POP}.ctx.sv\n";
    }
    else{
        print "/scratch/xfan3/pkg/breakdancer/breakdancer-max -t $cfg_name > $cfg_name.ctx.sv\n";
    }
    print "#\n";

}

# read the samples and store to a hash
sub read_sup_sample{
    my $line = shift;
    my @a = split(/:/, $line);
    my $hash;
    foreach (@a){
        if($_ =~ /FULL_1KGENOME\/(.+).mapped/){
            $hash->{$1} = 1;
        }
    }
    return $hash;
}

# detect inverted tandem dup, including the first and second copy inversion
# start from the two bunches of reads pointing to the same direction on the same pos, looking for the opposite direction on the same position, and its pair (the other breakpoint). See my drawing for illustartion.
sub invtdup_detection{
    my ($bd, $a_, $b_, $c_) = @_;
    # a_ is the smallest dup length, b_ is the highest distance between the first same direction pair reads. c_ is the highest distance between the first pair and the second pair for the read closer to the first pair
    my $bd_hash = &read_breakdancer_for_idup($bd);
    # bd_hash->{ori}->{chr}->{pos}->{samples}, {pair}
    my %ori_hash = ('+'=>'-',
        '-'=>'+'
    );
    foreach my $ori (keys %{$bd_hash}){
        my $p_ori = $ori_hash{$ori};
        foreach my $chr (keys %{$bd_hash->{$ori}}){
            foreach my $pos1 (sort {$a <=> $b} keys %{$bd_hash->{$ori}->{$chr}}){
                # remove those that do not have the other end the same ori and same pos
                $DB::single = 1;
                next if($bd_hash->{$ori}->{$chr}->{$pos1}->{pair_ori} eq $p_ori);
                my ($p_chr1, $p_pos1) = split(/:/, $bd_hash->{$ori}->{$chr}->{$pos1}->{pair});
                # remove those that the pair is far away from this pos
                next if($p_chr1 ne $chr || abs($p_pos1 - $pos1) > $b_);
                # now search the other ori
                foreach my $pos2 (sort {$a <=> $b} keys %{$bd_hash->{$p_ori}->{$chr}}){
                    # not yet
                    next if($pos2 - $pos1 < -$c_);
                    # has passed
                    last if($pos2 - $pos1 > $c_);
                    # the other read should be of the same ori as this one
                    my $p_ori2 = $bd_hash->{$p_ori}->{$chr}->{$pos2}->{pair_ori};
                    next if($p_ori2 ne $p_ori);
                    # ignore those with relative positions and orientations not make sense
                    my ($p_chr2, $p_pos2) = split(/:/, $bd_hash->{$p_ori}->{$chr}->{$pos2}->{pair});
                    next if($p_chr2 ne $chr);                
                    if($ori eq "+"){
                        next if($p_pos2 > $pos2);
                    }
                    else{
                        next if($p_pos2 < $pos2);
                    }

                    # dup length restriciton
                    my $length = abs($p_pos2 - $pos2);
                    next if($length < $a_);

                    # for print and depth check
                    my @samples;
                    my @rds;
                    $DB::single = 1;
                    foreach my $sample (keys %{$bd_hash->{$ori}->{$chr}->{$pos1}->{samples}}){
                        if(defined $bd_hash->{$p_ori}->{$chr}->{$pos2}->{samples}->{$sample}){
                            push @samples, $sample;

                            my ($p_pos_1, $p_pos_2) = ($p_pos2 < $pos2 ? $p_pos2:$pos2, $p_pos2 < $pos2 ? $pos2:$p_pos2);
                            my ($p_pos_s, $p_pos_l) = ($p_pos_1 - $length, $p_pos_2 + $length);
                            my $region_s = "$chr:$p_pos_s-$p_pos_1";
                            my $region_l = "$chr:$p_pos_2-$p_pos_l";
                            my $region = "$chr:$p_pos_1-$p_pos_2";
                            #push @rds, bam::count_reads_1KG($sample, $region_s);
                            #push @rds, bam::count_reads_1KG($sample, $region);
                            #push @rds, bam::count_reads_1KG($sample, $region_l);
                        }
                    }
                    my ($p1, $p2) = ($pos1 < $p_pos2 ? $pos1 : $p_pos2, $pos1 < $p_pos2 ? $p_pos2 : $pos1);
                    &print_idup($bd_hash, $chr, $p1, $p2, \@samples, \@rds) if(scalar(@samples) != 0);
                }
            }
        }
    }
}

# start from the third point, looking for interspersed duplication pattern in breakdancer file
sub idup_detection{
    my ($bd, $a_, $b_, $c_) = @_;
    # a_ is the smallest dup length, b_ is the highest distance between the two third points, c_ is how far away the third point thousld be from the toher two points if on the same chr
    my $bd_hash = &read_breakdancer_for_idup($bd);
    # bd_hash->{ori}->{chr}->{pos}->{samples}, {pair}
    foreach my $chr (keys %{$bd_hash->{"+"}}){
        foreach my $pos1 (sort {$a <=> $b} keys %{$bd_hash->{"+"}->{$chr}}){
            foreach my $pos2 (sort {$a <=> $b} keys %{$bd_hash->{"-"}->{$chr}}){
                last if($pos2 - $pos1 > $b_);
                next if($pos1 - $pos2 > $b_);
                # check if these two are the pair
                my ($p_chr1, $p_pos1) = split(/:/, $bd_hash->{"+"}->{$chr}->{$pos1}->{pair});
                my ($p_chr2, $p_pos2) = split(/:/, $bd_hash->{"-"}->{$chr}->{$pos2}->{pair});
                next if($p_chr1 eq $chr && $p_pos1 == $pos2 || $p_chr2 eq $chr && $p_pos2 == $pos1);

                # check to see if the third point is far away from the region
                next if($p_chr1 eq $chr && (abs($p_pos1 - $pos1) < $c_ || abs($p_pos2 - $pos1) < $c_) || $p_chr2 eq $chr && (abs($p_pos1 - $pos2) < $c_ || abs($p_pos2 - $pos2) < $c_));

                # check if the otehr end on the same chromosome and the size
                next if(!($p_chr1 eq $p_chr2 && abs($p_pos1 - $p_pos2) > $a_));
                # remove those with the third point inside the region
                next if($p_chr1 eq $chr && ($pos1 - $p_pos1) * ($pos1 - $p_pos2) < 0);
                next if($p_chr2 eq $chr && ($pos2 - $p_pos1) * ($pos2 - $p_pos2) < 0);
                # check if the supporting samples overlap
                my @samples;
                my @rds;
                foreach my $sample (keys %{$bd_hash->{"+"}->{$chr}->{$pos1}->{samples}}){
                    if(defined $bd_hash->{"-"}->{$chr}->{$pos2}->{samples}->{$sample}){
                        # find it, store the sample
                        push @samples, $sample;
                        # calculate the read depth
                        my $p_pos_1 = $p_pos1 < $p_pos2 ? $p_pos1 : $p_pos2;
                        my $p_pos_2 = $p_pos1 <= $p_pos2 ? $p_pos2 : $p_pos1;
                        my $region = "$p_chr1:$p_pos_1-$p_pos_2";
                        my $len = abs($p_pos2 - $p_pos1);
                        my $p_pos_s = $p_pos_1 - $len;
                        my $p_pos_l = $p_pos_2 + $len;
                        my $region_s = "$p_chr1:$p_pos_s-$p_pos_1";
                        my $region_l = "$p_chr1:$p_pos_2-$p_pos_l";
                        push @rds, bam::count_reads_1KG($sample, $region_s);
                        push @rds, bam::count_reads_1KG($sample, $region);
                        push @rds, bam::count_reads_1KG($sample, $region_l);
                    }
                }

                $DB::single = 1;
                &print_idup($bd_hash, $chr, $pos1, $pos2, \@samples, \@rds) if(scalar(@samples) != 0);
            }
        }
    }
}

# print idup
sub print_idup{
    my ($hash, $chr, $pos1, $pos2, $samples, $rds) = @_;
    my $hash1 = $hash->{"+"}->{$chr}->{$pos1};
    my $hash2 = $hash->{"-"}->{$chr}->{$pos2};
    print join("\t", "$chr:$pos1-$pos2", $hash1->{pair}, $hash2->{pair}, join(":", @$samples));
    if(scalar(@$rds) != 0){
        print "\t";
        for(my $i = 0; $i < scalar(@$rds); $i ++){
            print $rds->[$i];
            if($i%3 == 2){
                print ";";
            }
            else{
                print "|";
            }
        }
    }
    print "\n";
}

# overlap with vcf with copy number gain breakpoints to detect interspersed dup
sub overlap_vcf_idup{
    my ($bd, $vcf, $a_s, $a_l, $b_m) = @_;
    # vcf is in the format of hash->{chr1}->{pos1}->{ID} and hash->{$chr2}->{pos2}->{ID}, and overall_hash->{ID}->{line} = 1
    my ($vcf_hash) = vcf::read_vcf($vcf, 5);
    my $bd_hash = &read_breakdancer($bd);
    my $gold_sand;
    $DB::single = 1;
    foreach my $key (keys %$bd_hash) {
        my @a = split(/:/, $key);
        my @chrs = ($a[0], $a[2]);
        my @poss = ($a[1], $a[3]);
        my $hash_sup_sample = &read_sup_sample($bd_hash->{$key});
        my $tag = 0;
        foreach my $chr(@chrs){

            if(defined $vcf_hash->{$chr}){
                foreach my $pos (sort {$a <=> $b} keys %{$vcf_hash->{$chr}}){
                    if($poss[$tag] > $pos + $a_l){
                        # keep on going
                        next;
                    }
                    elsif($poss[$tag] < $pos - $a_l){
                        # stop searching
                        last;
                    }
                    else{
                        $DB::single = 1;
                        # within the scope
                        my $vcf_ID = $vcf_hash->{$chr}->{$pos}->{ID};
                        my $str1 = $chrs[1-$tag];
                        my $str2 = $poss[1-$tag];
                        # compare the other side to determine if tdup
                        my $other = $vcf_hash->{$chr}->{$pos}->{other};
                        my ($other1, $other2) = split(/:/, $other);
                        if($str1 eq $other1 && abs($str2 - $other2) < $a_s){
                            $DB::single = 1;
                            # tandem duplication
                        }
                        elsif($str1 eq $chr && ($str2 - $pos) * ($str2 - $other2) < 0){
                            # third_pos is inside of vcf two points
                        }
                        else{
                            $gold_sand->{$vcf_ID}->{$str1}->{$str2}->{bd_str} = $bd_hash->{$key};
                            $gold_sand->{$vcf_ID}->{$str1}->{$str2}->{first_pos} = "$chr:$pos";
                            $gold_sand->{$vcf_ID}->{$str1}->{$str2}->{sup_sample} = join(":", keys %$hash_sup_sample);
                        }
                        if(abs($poss[$tag] - $pos) < $a_s){
                            $DB::single = 1;
                            last;
                        }
                    }
                }
            }
            $tag ++;
        }
    }

    $DB::single = 1;
    # analyze $gold_sand and filter out sand
    foreach my $vcf_ID (keys %$gold_sand){
        foreach my $str1 (keys %{$gold_sand->{$vcf_ID}}){
            next if(scalar(keys %{$gold_sand->{$vcf_ID}->{$str1}}) == 1);
            my $third_pos = -1;
            my $first_pos = "";
            my $hash_sup;
            foreach my $str2 (keys %{$gold_sand->{$vcf_ID}->{$str1}}){
                if($third_pos == -1){
                    $third_pos = $str2;
                    $first_pos = $gold_sand->{$vcf_ID}->{$str1}->{$str2}->{first_pos};
                    my @sup_samples = split(/:/, $gold_sand->{$vcf_ID}->{$str1}->{$str2}->{sup_sample});
                    foreach (@sup_samples){
                        $hash_sup->{$_} = 1;
                    }
                    next;
                }
                if(abs($third_pos - $str2) < $b_m && $gold_sand->{$vcf_ID}->{$str1}->{$str2}->{first_pos} ne $first_pos){
                    # check if any sample overlapping
                    my @sup_samples = split(/:/, $gold_sand->{$vcf_ID}->{$str1}->{$str2}->{sup_sample});
                    foreach (@sup_samples){
                        if(defined $hash_sup->{$_}){
                            print join("\t", $vcf_ID, $str1, $third_pos, $str2, $_) . "\n";
                        }
                    }
                }
            }
        }

    }
}

# sort breakdancer by $_
sub sort_breakdancer{
    my ($file) = @_;
    my $hash = &read_breakdancer($file, "sort_by_all");
    foreach my $a (sort keys %$hash){
        print $a . "\n";
    }
}

1;
