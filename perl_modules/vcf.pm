use warnings;
use strict;
package vcf;
require fa;
require string;

# Function: This package manipulates the file format of vcf. 

sub rm_dup{
    # remove duplicates with the same coordinates
    my ($f) = @_;
    my $hash = &read_vcf($f);
    foreach my $key (sort keys %$hash){
        $DB::single = 1;
        my @a = split(/:/, $hash->{$key});

        print join("\t", @a)."\n";
    }
}

sub read_vcf{
    # It reads a vcf formatted file and record it to a hash in a form of chr1:pos1:chr2:pos2:type -> line
    my ($f, $mode) = @_;
    $mode = 0 if(!defined $mode);
    my $hash;
    open fh_, "<$f" or die $!;
    while(<fh_>){
        next if($_ =~ /^#/);
        chomp;        
        my @a = split(/\t/, $_);
        next if(scalar(@a) == 0);
        if($mode == 0){
            if($_ =~ /CHR2=(.+);END=(\d+);/){
                # chr1:start:chr2:end
                my $key = join(":", @a[0 .. 1], $1, $2); 
                $hash->{$key} = join(":", @a) if(!defined $hash->{$key});
            }
        }
        elsif($mode == 1){
            if($_ =~ /CHR2=(.+);END=(\d+);/){
                # chr1:start:chr2:end:type
                my $key = join(":", @a[0 .. 1], $1, $2); 
                if($_ =~ /SVTYPE=(\w)/){
                    $key .= ":$1";
                    $hash->{$key} = join(":", @a) if(!defined $hash->{$key});
                }
            }
        }
        elsif($mode == 2){
            # for /scratch/bcb/kchen3/1000genomes/Phase3/TIGRA/ALLPhase3merged/*.vcf format, no chr2
            # for find_pair_coord func
            if($_ =~ /END=(\d+);/){

                # start or end
                my $str = join(":", @a[0 .. 8]);
                my @t;
                @t = @{$hash->{$a[1]}} if(defined $hash->{$a[1]});
                push @t, $str;
                $hash->{$a[1]} = \@t;

                @t = @{$hash->{$1}} if(defined $hash->{$1});
                push @t, $str;
                $hash->{$1} = \@t;
            }
        }
        elsif($mode == 3){
            # the same as mode = 2, but now with the ID as the key
            my $key = $a[2];
            my $str = join(":", @a[0 .. 8]);
            $hash->{$key}->{str} = $str;
            my $num = 0;
            foreach (@a[9 .. $#a]){
                if($_ =~ /0\/1/ || $_ =~ /^1\/1/ || $_ =~ /^0\|1/ || $_ =~ /^1\|1/){
                    $num ++;
                }
            }
            $hash->{$key}->{num} = $num;
        }
        elsif($mode == 4){
            # chr1.start, and chr2.end as the key, ID as the value
            my $key = join(".", @a[0 .. 1]);
            $hash->{$key}->{ID} = $a[2];
            if($_ =~ /end=(\d+)/){
                my $end = $1;
                if($_ =~ /chr2=(\w+)/i){
                    $key = join(".", $1, $end);
                    $hash->{$key}->{id} = $a[2];
                }
                else{
                    $key = join(".", $a[0], $end);
                    $hash->{$key}->{ID} = $a[2];
                }
            }
        }
        elsif($mode == 5){
            # vcf_hash->{$chr}->{$pos1}->{ID}, vcf_hash->{$chr2}->{$pos2}->{ID}
            if($_ =~ /END=(\d+)/){
                my $end = $1;
                if($_ =~ /CHR2=(\w+)/i){
                    $hash->{$a[0]}->{$a[1]}->{ID} = $a[2];
                    $hash->{$1}->{$end}->{ID} = $a[2];
                    $hash->{$1}->{$end}->{other} = "$a[0]:$a[1]";
                    $hash->{$a[0]}->{$a[1]}->{other} = "$1:$end";
                }
                else{
                    $hash->{$a[0]}->{$a[1]}->{ID} = $a[2];
                    $hash->{$a[0]}->{$end}->{ID} = $a[2];
                    $hash->{$a[0]}->{$end}->{other} = "$a[0]:$a[1]";
                    $hash->{$a[0]}->{$a[1]}->{other} = "$a[0]:$end";
                }
            }
        }
    }
    close fh_;
    return $hash;
}

sub read_vcf_w_spec_sample{
    # It reads a vcf formatted file and record it to a hash in a form of ID -> chr1, pos1, chr2, pos2, at which any of the specified samples have variants. Samples are split by comma (,).
    my ($f, $sample) = @_;
    my $hash;
    my @samples = split(/,/, $sample);
    open fh_, "<$f" or die $!;
    my @array;
    while(<fh_>){
        next if($_ =~ /^##/);
        chomp;        
        my @a = split(/\t/, $_);
        next if(scalar(@a) == 0);
        if($_ =~ /^#/){
            # record the column index of where the samples specified appear
            for(my $i = 9; $i <= $#a; $i ++){
                for(my $j = 0; $j <= $#samples; $j ++){
                    if($a[$i] eq $samples[$j]){
                        push @array, $i;
                    }
                }
            }
            next;
        }
 
        # look at only those with specified samples have variants
        my $tag = 0;
        foreach my $col (@array){
            if($a[$col] =~ /^[01][\/\|]1/){
                $tag = 1;
                last;
            }
        }

        next if($tag == 0);
        # record variant samples
        my $key = $a[2];
        my ($chr1, $start) = @a[0 .. 1];
        my ($chr2, $end);
        if($_ =~ /;end=(\d+);/i){
            $end = $1;
            if($_ =~ /chr2=(\w+)/){
                $chr2 = $1;
            }
            else{
                # assume the same if no chr2
                $chr2 = $chr1;
            }
        }
        else{
            $chr2 = $chr1;
            $end = $start + 1;
        }
        $hash->{$key} = "$chr1:$start:$chr2:$end:$a[4]:$a[7]";

    }
    close fh_;

    return $hash;
}



sub read_vcf_w_sample{
    # It reads a vcf formatted file and record it to a hash in a form of ID -> chr1, pos1, chr2, pos2, @samples with variants
    my ($f, $mode) = @_;
    $mode = 0 if(!defined $mode);
    my $hash;
    my @samples;
    open fh_, "<$f" or die $!;
    while(<fh_>){
        next if($_ =~ /^##/);
        chomp;        
        my @a = split(/\t/, $_);
        next if(scalar(@a) == 0);
        if($_ =~ /^#/){
            @samples = @a[9 .. $#a];
        }


        if($mode == 0){
            # record variant samples
            my $key = $a[2];
            my ($chr1, $start) = @a[0 .. 1];
            my ($chr2, $end);
            if($_ =~ /end=(\d+)/){
                my $end = $1;
                if($_ =~ /chr2=(\w+)/){
                    $chr2 = $1;
                }
                else{
                    # assume the same if no chr2
                    $chr2 = $chr1;
                }
            }
            $hash->{$key}->{str} = "$chr1:$start:$chr2:$end";
            # deal with samples
            my @sample;
            my $tag = 0;
            foreach (@a[9 .. $#a]){
                if($_ =~ /0\/1/ || $_ =~ /^1\/1/ || $_ =~ /^0\|1/ || $_ =~ /^1\|1/){
                    push @sample, $samples[$tag];
                }
                $tag ++;
            }
            $hash->{$key} = \@sample;
        }

    }
    close fh_;

    return $hash;
}

sub aligned_calls{
    # this function reads a vcf file, with the last column containing original call info, and connect with the original vcf hash, which was read by mode = 4 (hash from mode = 4, hash1 from mode = 3);
    my ($hash, $hash1, $vcf_file) = @_;
    open VCF, "<$vcf_file" or die $!;

    while(<VCF>){
        chomp;
        my @a = split(/\t/, $_);
        my @b = split(/\./, $a[$#a - 4]);
        my $key = join(".", @b[2 .. 3]);
        $DB::single = 1;
        if(defined $hash->{$key}->{ID}){
            my $ID = $hash->{$key}->{ID};
            $hash1->{$ID}->{aligned} = 1 if(defined $hash1->{$ID});
        }
    }
    close VCF;
    return $hash1;
}

sub analyze_calls{
    # it divides into Delly vs otheres, precise vs imprecise, DUP vs CNV, common vs rare
    my ($vcf) = @_;
    my $hash = &read_vcf($vcf, 3);
    foreach my $key (sort keys %$hash){
        $hash->{$key}->{delly} = 0;
        $hash->{$key}->{imprecise} = 0;
        $hash->{$key}->{dup} = 0;
        $hash->{$key}->{cnv} = 0;
        if($key =~ /delly/i){
            $hash->{$key}->{delly} = 1;
        }
        if($hash->{$key}->{str} =~ /imprecise/i){
            $hash->{$key}->{imprecise} = 1;
        }
        my @a = split(/:/, $hash->{$key}->{str});
        if($key =~ /DUP/i){
            $hash->{$key}->{dup} = 1;
        }
        elsif($key =~ /CN/i){
            $hash->{$key}->{cnv} = 1;
        }
    }
    return $hash;
}


sub distinguish_interspersed{
    # given a vcf file, find tandum and interspersed dup
    my ($vcf) = @_;
    my $s_size = 1000; 
    my $l_size = 3000;
    my $size = 10000;
    my $third_range = 500; # the distance between two bpts indicated from two alignment (the third hidden bpt)
    my $str = "";
    my $hash;
    open VCF, "<$vcf" or die $!;
    while(<VCF>){
        $str = $_;
        my @a = split(/\s+/, $_);
        # called should be on the same chr
        my ($chr1) = ($_ =~ /CHR2=(\w+)/);
        next if($a[0] ne $chr1);

        my @b = split(":", $a[$#a]);
        my ($chr, $s, $e) = @b[1 .. 3];
        my ($e1) = ($_ =~ /END=(\d+)/);
        my $s1 = $a[1];

        if($e - $s > $size){
            # deal with large ones
            if(abs($e - $e1) < $s_size && abs($s - $s1) < $s_size){
                # tandum
                $str = s/<DUP>/<TANDUP>/;
            }
            else{
                my $str_ = "$chr:$s-$e";
                # putative interspersed
                if(abs($e - $e1) > $s_size){
                    $hash->{$str_}->{$e1}->{pos} = $s1;
                    $hash->{$str_}->{$e1}->{line} = $_;
                }
                if(abs($s - $s1) > $s_size){
                    $hash->{$str_}->{$s1}->{pos} = $e1;
                    $hash->{$str_}->{$s1}->{line} = $_;
                }
                if($s1 !~ /\d+/ || $e1 !~ /\d+/){
                    $DB::single = 1;
                }
            }
        }
        else{ # deal with small ones
            if(abs($e - $e1) < $l_size && abs($s - $s1) < $l_size){
                # tandum
                my @c = split(/\t/, $str);
                print join("\t", @c[0 .. 3], "<TANDUP>", @c[5 .. $#c]);
            }
            else{
                # putative interspersed
                my $str_ = "$chr:$s-$e";
                if(abs($e - $e1) > $l_size){
                    $hash->{$str_}->{$e1}->{pos} = $s1;
                    $hash->{$str_}->{$e1}->{line} = $_;
                }
                if(abs($s - $s1) > $l_size){
                    $hash->{$str_}->{$s1}->{pos} = $e1;
                    $hash->{$str_}->{$s1}->{line} = $_;
                } 
            }

        }
    }
    close VCF;

    foreach my $key (sort keys %$hash){
        my ($chr, $s, $e) = ($key =~ /^(.+):(\d+)-(\d+)/);
        my @pos = sort {$a <=> $b} keys %{$hash->{$key}};
        next if(scalar(@pos) == 1);
        for(my $i = 0; $i < scalar(@pos) - 1; $i ++){
            for(my $j = $i + 1; $j < scalar(@pos); $j ++){
                if($pos[$j] - $pos[$i] < $third_range){
                    # check pattern
                    my $key1 = $pos[$i];
                    my $key2 = $pos[$j];
                    my $pos1 = $hash->{$key}->{$key1}->{pos};
                    my $pos2 = $hash->{$key}->{$key2}->{pos};
                    if($e !~ /^\d+$/ || $s !~ /^\d+$/ || $pos1 !~ /^\d+$/ || $pos2 !~ /^\d+$/){
                        $DB::single = 1;
                    }
                    if(abs($e - $s) < $size && (abs($pos1 - $s) < $s_size && abs($pos2 - $e) < $s_size || abs($pos1 - $e) < $s_size && abs($pos2 - $s) < $s_size)){
                        $str = $hash->{$key}->{$key1}->{line};
                        my @a = split(/\t/, $str);
                        print join("\t", @a[0 .. 3], "<INTDUP>", @a[5 .. $#a]);
                        $str = $hash->{$key}->{$key2}->{line};
                        @a = split(/\t/, $str);
                        print join("\t", @a[0 .. 3], "<INTDUP>", @a[5 .. $#a]);
                    }
                    elsif(abs($e - $s) >= $size && (abs($pos1 - $s) < $l_size && abs($pos2 - $e) < $l_size || abs($pos1 - $e) < $l_size && abs($pos2 - $s) < $l_size)){
                        $str = $hash->{$key}->{$key1}->{line};
                        my @a = split(/\t/, $str);
                        print join("\t", @a[0 .. 3], "<INTDUP>", @a[5 .. $#a]);
                        $str = $hash->{$key}->{$key2}->{line};
                        @a = split(/\t/, $str);
                        print join("\t", @a[0 .. 3], "<INTDUP>", @a[5 .. $#a]);
                    }
                }
                else{
                    last;
                }
            }
        }
    }

}

sub get_var{
    # it detects the variant by 0/1 or 1/1 for each column after 8 and output the bam file name
    my ($f) = @_;
    my @bams;
    my $vars;

    open fh_, "<$f" or die $!;
    while(<fh_>){
        next if($_ =~ /^##/);
        chomp;
        my @a = split(/\t/, $_);
        if($_ =~ /^#CHROM/){
            @bams = @a[9 .. $#a];
        }
        else{
            my $i = 0;
            foreach my $item (@a[9 .. $#a]){
                my @b = split(/:/, $item);
                if($b[0] ne "0" && $b[0] ne "0/0" && $b[0] ne "0|0"){
                    #if($item =~ /^0\/1/ || $item =~ /^1\/1/ || $item =~ /^0\|1/ || $item =~ /^1\|1/){
                    my $key = join(":", @a[0 .. 7]);
                    push @{$vars->{$key}}, "$bams[$i]";
                }
                $i ++;
            }
        }
    }
    foreach my $key (sort keys %$vars){
        print join("\t", $key, @{$vars->{$key}}) . "\n";
    }
    return $vars;
}

sub find_coord{
    # given an array of coordinates, and optionally chromosome if they are in the same chr, output those lines containing any of the coordinate to a new vcf file
    my ($f, $coords, $new_f, $chr) = @_;

    open fh_, "<$f" or die $!;
    open fh_new, ">$new_f" or die $!;
    while(<fh_>){
        if($_ =~ /^#/){
            print fh_new $_;
        }
        else{
            next if($chr ne "" && $_ !~ /^$chr/); 
            foreach my $coord (@$coords){
                if($_ =~ /$coord/){
                    print fh_new $_;
                    last;
                }
            }
        }

    }

    close fh_;
    close fh_new;

}

# Given a vcf file and a fasta file, extract the supporting contig name and output it to a fasta file.
sub extract_sup_allele{
    my ($vcf, $fa, $out_fa) = @_;
    my $hash = &read_vcf($vcf);
    #&write_vcf($hash, $out_vcf);
    my $fa_hash = fa::read_fa($fa);
    my @keys;
    foreach my $key (sort keys %$hash){
        my @a = split(/:/, $hash->{$key});
        push @keys, $a[$#a - 4];
    }

    fa::write_fa($fa_hash, \@keys, $out_fa);
}

# Display a certain line in vcf for the first 9 columns
# It also outputs the header line #CHROM
sub display_vcf{
    my ($sig, $vcf) = @_;
    open VCF, "<$vcf" or die $!;
    while(<VCF>){
        next if($_ =~ /^##/);
        my @a = split(/\t/, $_);
        print join("\t", @a[0 .. 8]) . "\n" if($_ =~ /^#CHROM/ || $_ =~ /$sig/);
    }
    close VCF;
}


sub find_pair_coord{
    # with one of the bpt indicated in contig name of the last column in a vcf, find the other bpt in the other vcf and append the coordinates to the end of the current one (supplement with vcf2 to vcf1)
    my ($vcf1, $vcf2) = @_;

    my $mode = 2;
    my $hash = &read_vcf($vcf2, $mode);
    open FH, "<$vcf1" or die $!;
    while(<FH>){
        next if($_ =~ /^#/);
        chomp;
        my @a = split(/\t/, $_);
        my @t = split(/\./, $a[$#a]);
        my $bpt = $t[3];
        $DB::single = 1;
        if(defined $hash->{$bpt}){
            my $strs = $hash->{$bpt};
            my $str = $strs->[0];
            my @w = split(/:/, $str);
            if($w[7] =~ /SVTYPE=(\w+).+END=(\d+)/){
                my $append = join(":", $w[0], $w[1], $2, $1);
                print $_ . ":" . $append . "\n";
            }
        }
    }
    close FH;
}

sub comp_crossmatch_output{
    # it compare a vcf file from Zechen's and a crossmatch result from Ken's, output in vcf format of those calls overlapping with crossmatch's. Allowing a shift of breakpoint of $num_shift
    my ($crossmatch_file, $vcf_file, $num_shift) = @_;
    my $vcf_hash = &read_vcf($vcf_file);
    my $start;
    my $end;
    open CM_FH, "$crossmatch_file" or die $!;
    while(<CM_FH>){
        next if($_ =~ /^#/);
        my @a = split(/\t/, $_);
        if($a[2] =~ /(\d+)\(/){
            $start = $1;
        }
        if($a[4] =~ /(\d+)\(/){
            $end = $1;
        }
        $DB::single = 1;
        my $break = 0;
        # allow for a subtle difference
        foreach my $i (-$num_shift .. $num_shift){
            foreach my $j (-$num_shift .. $num_shift){
                my $key = join(":", ($a[1], $start + $i, $a[3], $end + $j));
                if(defined $vcf_hash->{$key}){
                    $break = 1;
                    last;
                }
            }
            last if($break == 1);
        }
        print $_ if($break != 1);
    }
    close CM_FH;
}

sub comp_bd_output{
    # it compare a vcf file from Zechen's and a crossmatch result from Ken's, output in vcf format of those calls overlapping with crossmatch's. Allowing a shift of breakpoint of $num_shift
    # breakdancer output formatic: bd = 1
    my ($crossmatch_file, $vcf_file, $num_shift, $bd) = @_;
    my $vcf_hash = &read_vcf($vcf_file);
    my $vcf_hash_id;
    my $hash_miss;
    my $start;
    my $end;
    my $chr1;
    my $chr2;
    my $type;
    my $hash;
    my $hash_str;
    my $header;
    my $stat_h;
    # order of numerating
    my @order = (0);
    foreach my $i (1 .. $num_shift){
        push @order, $i;
        push @order, -$i;
    }
    open CM_FH, "$crossmatch_file" or die $!;
    while(<CM_FH>){
        chomp;
        my @a = split(/\s+/, $_);
        if($_ =~ /^#/){
            my $i = 0;
            foreach my $entry (@a){
                $header->{chr1} = $i if($entry =~ /chr1/i);
                $header->{start} = $i if($entry =~ /start/i);
                $header->{chr2} = $i if($entry =~ /chr2/i);
                $header->{end} = $i if($entry =~ /end/i);
                $header->{type} = $i if($entry =~ /type/i);
                $header->{micro} = $i if($entry =~ /micro/i);
                $header->{template} = $i if($entry =~ /template/i);
                $header->{nontemplate} = $i if($entry =~ /nontemplate/i);
                $i ++;
            }
            next;
        }
        my ($micro, $template, $nontemplate) = ("NA", "NA", "NA");
        $start = $a[$header->{start}];
        $end = $a[$header->{end}];
        $chr1 = $a[$header->{chr1}];
        $chr2 = $a[$header->{chr2}];
        $type = $a[$header->{type}] if(defined $header->{type});

        $micro = $a[$header->{micro}] if(defined $header->{micro});
        $template = $a[$header->{template}] if(defined $header->{template});
        $nontemplate = $a[$header->{nontemplate}] if(defined $header->{nontemplate});


        my $break = 0;
        # allow for a subtle difference
        foreach my $i (@order){
            foreach my $j (@order){
                my $key1 = join(":", ($chr1, $start + $i, $chr2, $end + $j));
                my $key2 = join(":", ($chr2, $end + $j, $chr1, $start + $i));
                if(defined $vcf_hash->{$key1} || defined $vcf_hash->{$key2}){
                    # need to specify those that are not identified
                    $vcf_hash_id->{$key1} = 1 if(defined $vcf_hash->{$key1});
                    $vcf_hash_id->{$key2} = 1 if(defined $vcf_hash->{$key2});
                    # check sv type
                    my $vcf_str;
                    $vcf_str = $vcf_hash->{$key1} if(defined $vcf_hash->{$key1});
                    $vcf_str = $vcf_hash->{$key2} if(defined $vcf_hash->{$key2});
                    my ($type_comp, $micro_comp, $template_comp, $nontemplate_comp);
                    my ($vcf_type) = ($vcf_str =~ /SVTYPE=(\w+)/);
                    ($type_comp, $stat_h) = &check_equal($vcf_type, $type, "Type", $stat_h);
                    my ($vcf_micro) = ($vcf_str =~ /MICRO=(\w+)/);
                    ($micro_comp, $stat_h) = &check_equal($vcf_micro, $micro, "MICRO", $stat_h);
                    my ($vcf_template) = ($vcf_str =~ /;TEMPLATE=(\w+)/);
                    ($template_comp, $stat_h) = &check_equal($vcf_template, $template, "TEMPLATE", $stat_h, $i, $j);
                    my ($vcf_non_template) = ($vcf_str =~ /NONTEMPLATE=(\w+)/);
                    ($nontemplate_comp, $stat_h) = &check_equal($vcf_non_template, $nontemplate, "NONTEMPLATE", $stat_h, $i, $j);

                    my $off_str;
                    if(defined $vcf_micro && abs($i) == abs($j) && abs($i) <= length($vcf_micro)){
                        $off_str = join("\t", "$i(0)", "$j(0)");
                        $hash->{0} += 2;
                    }
                    else{
                        $off_str = join("\t", $i, $j);
                        $hash->{$i} ++;
                        $hash->{$j} ++;
                    }

                    # deal with false positive
                    if($break != 1){
                        my $str = join("\t", $_, $off_str, $type_comp, $micro_comp, $template_comp, $nontemplate_comp) . "\n";
                        $hash_str->{$str} = 1;
                    }
                    $break = 1;

                    # last;
                }
            }
            #last if($break == 1);
        }
        if($break != 1){
            $hash_miss->{"$chr1.$start.$chr2.$end.$type"} = 1;
        }
    }
    close CM_FH;
    print "Hit (TP) : " . scalar(keys %$hash_str) . "\n";
    foreach my $key (sort keys %$hash_str){
        print $key;
    }
    foreach my $key (sort keys %$stat_h){
        print join("\t", $key, $stat_h->{$key}) . "\n";
    }
    foreach my $i (sort {$a <=> $b} keys %$hash){
        print join("\t", $i, $hash->{$i}) . "\n";
    }

    # print those not identified
    # false positives
    my $FP = 0;
    print "False positives\n";
    foreach my $key (keys %$vcf_hash){
        if(!defined $vcf_hash_id->{$key}){
            print $vcf_hash->{$key} . "\n";
            $FP ++;
        }
    }
    print "In all we have $FP false positives.\n";
    print "Missing: " . scalar(keys %$hash_miss) . "\n";
    foreach my $key (keys %$hash_miss){
        print $key . "\n";
    }
}

sub check_equal{
    my ($vcf_entry, $bd_entry, $key, $h, $i, $j) = @_;
    return "No$key" if($bd_entry eq "" || !defined $vcf_entry);
    if($key =~ /type/i){
        if($vcf_entry eq $bd_entry || $vcf_entry eq "DUP" && $bd_entry eq "ITX" || $bd_entry =~ /Amplified/i){
            $h->{"type:Correct"} ++;
            return ("$key:correct", $h);
        }
        else{
            $h->{"type:Incorrect"} ++;
            return ("$key:$vcf_entry", $h);
        }
    }
    else{
        # micro, template, nontemplate
        if($bd_entry eq "NA"){
            # true value is nothing
            if($vcf_entry eq "NA"){
                $h->{"$key:TN"} ++;
                return ("$key:TN", $h);
            }
            else{
                # FP
                $h->{"$key:FP"} ++;
                return ("$key:FP:$vcf_entry", $h);
            }
        }
        else{
            # true value has something
            if(defined $vcf_entry && $vcf_entry eq "NA"){
                $h->{"$key:FN"} ++;
                return ("$key:FN", $h);
            }
            else{
                # both have something, compare
                my $str1 = "$key:TP";
                $h->{"$key:TP"} ++;
                my $ret_str;
                $ret_str = string::check_same($bd_entry, $vcf_entry) if($key =~ /micro/i);
                $ret_str = string::check_same($bd_entry, $vcf_entry, $i, $j) if($key =~ /template/i);
                if($ret_str == 1){
                    $h->{"$key:correctStr"} ++;
                    return ("$str1:correctStr", $h);
                }
                elsif($ret_str == 2){
                    $h->{"$key:correctSubStr"} ++;
                    return ("$str1:correctSubStr:$vcf_entry", $h);
                }
                else{
                    $h->{"$key:incorrectStr"} ++;
                    return ("$str1:$vcf_entry", $h);
                }
            }
        }
    }
}

1;
