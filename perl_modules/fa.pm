use warnings;
use strict;
require vcf;
require breakdancer;
require tmpfile;
#use Cwd 'abs_path';
#use File::Basename;
package fa;

# Function: This package manipulates the file format of fa. 

my $line_length = 60;

# given an fa with readnames _1 and _2, put _1 in fa1 and _2 in fa2
sub separate_fa{
    my ($fa, $fa1, $fa2) = @_;
    my $tag = 0;
    open fh_1, ">$fa1" or die $!;
    open fh_2, ">$fa2" or die $!;
    open fh_, "<$fa" or die $!;
    while(<fh_>){
        if($_ =~ /\_[1|2]$/ && $_ !~ /^>/){
            $_ = ">$_";
        }
        if($_ =~ /^(>.+)\_1$/){
            $tag = 1;
            print fh_1 $1 . "\n";
        }
        elsif($_ =~ /^(>.+)\_2$/){
            $tag = 2;
            print fh_2 $1 . "\n";
        }
        elsif($tag == 1){
            print fh_1 $_;
        }
        elsif($tag == 2){
            print fh_2 $_;
        }
    }
    close fh_;
    close fh_1;
    close fh_2;
}


# find the reads in fa1 but not in fa2, output them to fa3
sub take_diff{
    my ($fa1, $fa2, $fa3) = @_;
    my $h;
    open fa_, "<$fa2" or die $!;
    while(<fa_>){
        if($_ =~ /^>/){
            chomp;
            $h->{$_} = 1;
        }
    }
    close fa_;
    my $tag = 0;
    open fa_, "<$fa1" or die $!;
    open fh_3, ">$fa3" or die $!;
    while(<fa_>){
        if($_ =~ /^>/){
            chomp;
            if(! defined $h->{$_}){
                # print
                print fh_3 $_ . "\n";
                $tag = 1;
            }
            else{
                $tag = 0;
            }
        }
        else{
            if($tag == 1){
                print fh_3 $_;
            }
        }
    }
    close fa_;
    close fh_3;
}




# Mixing the contigs from different groups together, making a big contig fa file, with each read followed by .g$num in its suffix of readname, in which num is the number after "g" in the contig file name.
sub mix_group_fa{
    my ($dir, $out) = @_;
    my @files = split(/\n/, `ls $dir`);
    open OUT, ">$out" or die $!;
    foreach my $f (@files){
        my ($g) = ($f =~ /g(\d+)/);
        open FH, "<$f" or die $!;
        while(<FH>){
            if($_ =~ /^>/){
                chomp;
                print OUT "$_.g$g\n";
            }
            else{
                print OUT $_;
            }
        }
        close FH;
    }
    close OUT;
}

sub index_rn{
    my ($fa) = @_;
    my $h;
    my $n = 0;
    open FA, "<$fa" or die $!;
    while(<FA>){
        if($_ =~ /^>/){
            chomp;
            $_ =~ s/^>//;
            $h->{$_} = $n;
            $n ++;
        }
    }
    close FA;
    return $h;
}

# remove a suffix string
sub rm_str_rn {
    my ($fa, $str) = @_;
    open FA, "<$fa" or die $!;
    while(<FA>){
        if($_ =~ /^(>.+)$str$/){
            print "$1\n";
        }
        else{
            print $_;
        }
    }
    close FA;
}



# add a string at the beginning of all readnames, separated by separation
sub add_str_rn {
    my ($fa, $str, $sep) = @_;
    open FA, "<$fa" or die $!;
    while(<FA>){
        if($_ =~ /^>/){
            chomp;
            $_ =~ s/^>//;
            print ">$str" . $sep . $_ . "\n";
        }
        else{
            print $_;
        }
    }
    close FA;
}

sub readlen{
    my ($fa) = @_;
    # calculate the length of each read, output each line "readname readlength"
    my $readlen = 0;
    my $readname = "NA";
    my @str;
    open FA, "<$fa" or die $!;
    while(<FA>){
        chomp;
        if($_ =~ /^>/){
            # process the previous one if any
            if($readname ne "NA"){
                if($readlen == 0){
                    print "Warning: $readname has read length 0\n";
                }
                else{
                    push @str, join("\t", $readname, $readlen);
                }
            }
            $_ =~ s/^>//;
            # prepare for this read
            $readname = $_;
            $readlen = 0;
        }
        else{
            $readlen += length($_);
        }
    }
    close FA;
    push @str, join("\t", $readname, $readlen);
    #print join("\n", @str);
    return join("\n", @str);

}

sub readlen_fq{
    my ($fq) = @_;
    # calculate the length of each read, output each line "readname readlength"
    my $readlen = 0;
    my $readname = "NA";
    my $tag = 0;
    open FQ, "<$fq" or die $!;
    while(<FQ>){
        chomp;
        if($_ =~ /^\@/){
            # process the previous one if any
            if($readname ne "NA"){
                if($readlen == 0){
                    print "Warning: $readname has read length 0\n";
                }
                else{
                    print join("\t", $readname, $readlen) . "\n";
                }
            }
            $_ =~ s/^\@//;
            # prepare for this read
            $readname = $_;
            $readlen = 0;
            # let it know the next line can be counted in
            $tag = 1;
        }
        elsif($_ !~ /^\+$/ && $tag == 1){
            $readlen += length($_);
        }
        elsif($_ =~ /^\+$/){
            $tag = 0;
        }
        elsif($_ !~ /^\+$/ && $tag == 0){
            next;
        }
    }
    close FQ;
    print join("\t", $readname, $readlen) . "\n";
}


sub fa2fq{
    my ($fa) = @_;
    # convert an fa to fq. 
    my $seq = "";
    open FA, "<$fa" or die $!;
    while(<FA>){
        if($_ =~ /^>/){
            # process the previous read
            if(length($seq) != 0){
                print join("\n", $seq, "+", "I" x length($seq))."\n";
                $seq = "";
            }
            $_ =~ s/^>/@/;
            print $_;
        }
        else{
            chomp;
            # concatenate the sequence together
            $seq .= $_;
        }
    }
    if(length($seq) != 0){
        print join("\n", $seq, "+", "I" x length($seq))."\n";
    }
    close FA;
}

sub fa_by_len{
    my ($fa, $threshold) = @_;
    my $st = "";
    my $nm = "NA";
    open FA, "<$fa" or die $!;
    while(<FA>){
        chomp;
        if($_ =~ /^>/){
            if(length($st) > $threshold){
                print join("\n", $nm, $st) . "\n";
                $st = "";
            }
            $nm = $_;
        }
        else{
            $st .= $_;
        }
    }
    if(length($st) > $threshold){
        print join("\n", $nm, $st);
    }
    close FA;
}

# given an fa, this function index the readname, output a new fa and index file
# the new fa is called prefix.new.fa, index is prefix.index
# note: it dumps the second and following columns
sub fa_index{
    my ($fa_prefix) = @_;
    my $fa;
    my $index = $fa_prefix . ".index";
    if(-e "$fa_prefix.fa"){
        $fa = "$fa_prefix.fa";
    }
    elsif(-e "$fa_prefix.fasta"){
        $fa = "$fa_prefix.fasta";
    }
    else{
        die "There is no $fa.fa or $fa.fasta\n";
    }
    open id_fh, ">$index" or die $!;
    open new_fa_fh, ">$fa_prefix.new.fa" or die $!;
    open old_fa_fh, "<$fa" or die $!;
    my $num1 = 0;
    while(<old_fa_fh>){
        if($_ =~ /^>/){
            chomp;
            $_ =~ s/^>//;
            my @a = split(/\s+/, $_);
            print id_fh join("\t", $a[0], $num1) . "\n";
            print new_fa_fh ">$num1\n";
            $num1 ++;
        }
        else{
            print new_fa_fh $_;
        }
    }
    close id_fh;
    close new_fa_fh;
    close old_fa_fh;
}





# this is a temporary function for indexing read name and convert fq to fa
# write the index to fq_prefix.index
sub fq2fa_index{
    my ($fq_prefix, $suffix) = @_;
    my $fa = $fq_prefix . ".fa";
    my $sf = "";
    $sf = "_$suffix" if(defined $suffix && $suffix ne "" && $suffix ne "NA");
    my $index = $fq_prefix . ".index";
    my $num = 0;
    my $num1 = 0;
    my $hash;
    open id_fh, ">$index" or die $!;
    open fa_fh, ">$fa" or die $!;
    my $fq;
    if(-e "$fq_prefix.fq"){
        $fq = "$fq_prefix.fq";
    }
    elsif(-e "$fq_prefix.fastq"){
        $fq = "$fq_prefix.fastq";
    }
    else{
        die "There is no $fq.fq or $fq.fastq\n";
    }
    open FQ, "<$fq" or die $!;
    while(<FQ>){
        if($num%4 == 0){
            chomp;
            $_ =~ s/^@//;
            my @a = split(/\s+/, $_);
            print id_fh join("\t", $num1, $a[0]) . "\n";
            print fa_fh ">$num1" . "$sf\n";
            $num1 ++;
        }
        elsif($num%4 == 1){
            print fa_fh $_;
        }
        $num ++;
    }
    close FQ;
    close fa_fh;
    close id_fh;
}


sub fq_to_fa{
    my ($fq_prefix, $suffix) = @_;
    my $num = 0;
    my $sf = "";
    $sf = "_$suffix" if(defined $suffix && $suffix ne "" && $suffix ne "NA");
    my $fq;
    if(-e "$fq_prefix.fq"){
        $fq = "$fq_prefix.fq";
    }
    elsif(-e "$fq_prefix.fastq"){
        $fq = "$fq_prefix.fastq";
    }
    elsif(-e $fq_prefix && ($fq_prefix =~ /\.fq/ || $fq_prefix =~ /\.fastq/)){
        $fq = $fq_prefix;
        if($fq_prefix =~ /^(.+)\.fq$/){
            $fq_prefix = $1;
        }
        elsif($fq_prefix =~ /^(.+)\.fastq$/){
            $fq_prefix = $1;
        }
    }
    else{
        die "There is no $fq.fq or $fq.fastq\n";
    }

    open FA, ">$fq_prefix.fa" or die $!;
    open FQ, "<$fq" or die $!;
    while(<FQ>){
        if($num%4 == 0){
            chomp;
            $_ =~ s/^@/>/;
            print FA $_ . "$sf\n";
        }
        elsif($num%4 == 1){
            print FA $_;
        }
        $num ++;
    }
    close FQ;
    close FA;
}

sub display_fa{
# Display a certain line in vcf for the first 9 columns
# It also outputs the header line #CHROM
    my ($sig, $fa) = @_;
    open FA, "<$fa" or die $!;
    while(<FA>){
        print $_ if($_ =~ /^>/ && $_ =~ /$sig/);
    }
    close FA;
}

# given an original big vcf, a TIGRA contig fa file, return if it is assembled in a hash, input hash is from vcf::read_vcf(mode = 4), hash1 from mode = 3;
sub stats_assembled{
    my ($hash, $hash1, $fa_file) = @_;
    my $mode = 1;
    my $fa_hash = &read_fa($fa_file, $mode);
    foreach my $key (sort keys %$fa_hash){
        my @a = split(/\./, $key);
        my $key1 = join(".", @a[2 .. 3]);
        if(defined $hash->{$key1}){
            my $ID = $hash->{$key1}->{ID};
            $hash1->{$ID}->{assembled} = 1 if(defined $hash1->{$ID});
        }
    }
    return $hash1;
}



sub read_fa{
    # It reads a fa formatted file and record it to a hash in a form of readname -> seq
    my ($f, $mode) = @_;
    $mode = 0 if(!defined $mode);
    my $hash;
    my $key = "NA";
    return $hash if(!-e $f);
    open fh_, "<$f";
    while(<fh_>){
        chomp;
        if($_ =~ /^>(.+)$/){
            my @a = split(/\s+/, $1);
            $key = $a[0];
            $hash->{$key} = "" if($mode == 0);
            $hash->{$key} = 1 if($mode == 1);
        }
        else{
            $hash->{$key} .= $_ if($mode == 0);
        } 
    }
    close fh_;

    return $hash;
}

# Write a fasta with one contig, given contig name and contig seq
sub write_fa_simple{
    my ($contig_name, $contig_seq, $contig_file) = @_;
    open fh_, ">$contig_file" or die $!;
    print fh_ join("\n", ">$contig_name", $contig_seq) . "\n";
    close fh_;
}
# write fasta file but with extraction of certain range of sequence
# input is a hash with the key the readname, the value an array of the read coordinates
# output is a file that has all the sequences specified in the key and value, with the readname changed (start and end of the contig append with dot as a separtor
sub write_fa_extract{
    my ($hash, $key, $out, $fa, $append) = @_;
    if(defined $append && $append eq "a"){
        open OUT, ">>$out" or die $!;
    }
    else{
        open OUT, ">$out" or die $!;
    }
#	open OUT, ">>$out" or die $!;
    foreach my $name (@$key){
        if(defined $hash->{$name}){
            my $value = $hash->{$name};
            foreach my $v (@$value){
                my ($s, $e) = split(/\./, $v);
                my $str = `samtools faidx $fa $name:$s-$e`;
                if($str eq ""){
                    warn "samtools faidx $fa $name:$s-$e does not work!\n";
                }
                print "samtools faidx $fa $name:$s-$e\n";
                my @strs = split(/\n/, $str);
                $str = join("", @strs[1 .. $#strs]);
                print OUT ">$name.$s.$e\n";
                $str = &sep_line_length($str);
#            print OUT $hash->{$name};
                print OUT $str;
            }
        }
    }
    close OUT;
}


# Write fasta file
sub write_fa{
    my ($hash, $key, $out, $append) = @_;
    if(defined $append && $append eq "a"){
        open OUT, ">>$out" or die $!;
    }
    else{
        open OUT, ">$out" or die $!;
    }
#	open OUT, ">>$out" or die $!;
    foreach my $name (@$key){
        if(defined $hash->{$name}){
            print OUT ">$name\n";
            my $str = &sep_line_length($hash->{$name});
#            print OUT $hash->{$name};
            print OUT $str;
        }
    }
    close OUT;
}

# mode: del: get the flank size left of the left bp, right of hte right bp, construct the allele, output to the same fa for all dels.
# On Oct.15, add ins mode. The inserted sequence is on the last columns. 
sub config_fa_bed_file{
    my ($bed_file, $reference, $out, $flank, $mode) = @_;
    my $h;
    if($mode =~ /del/i){
        open fh_, "<$bed_file" or die $!;
        while(<fh_>){
            chomp;
            my @a = split(/\t/, $_);
            my $left = $a[1] - $flank;
            my $right = $a[2] + $flank; 
            if($a[1] - $flank < 0){
                $left = 0;
            }
            my $fa1 = &make_contig_simple($reference, $a[0], $left, $a[1]);
            my $fa2 = &make_contig_simple($reference, $a[0], $a[2], $right);
            my $str = "$a[0]:$left-$a[1];$a[0]:$a[2]-$right";
            $h->{$str} = &combine_contig($fa1, $fa2, "+", "-", -1, "NA");
        }
        close fh_;
        &write_fa($h, [keys %$h], $out);
    }
    elsif($mode =~ /ins/i){
        open fh_, "<$bed_file" or die $!;
        while(<fh_>){
            chomp;
            my @a = split(/\t/, $_);
            my $left = $a[1] - $flank;
            my $right = $a[1] + $flank; 
            if($a[1] - $flank < 0){
                $left = 0;
            }
            my $fa1 = &make_contig_simple($reference, $a[0], $left, $a[1]);
            my $fa2 = &make_contig_simple($reference, $a[0], $a[2], $right);
            my $str = "$a[0]:$left-$a[1];$a[0]:$a[2]-$right";
            my $str1 = &combine_contig($fa1, $a[$#a], "+", "-", -1, "NA");
            $h->{$str} = &combine_contig($str1, $fa2, "+", "-", -1, "NA");
        }
        close fh_;
        &write_fa($h, [keys %$h], $out);
    }

    #print "#Error: No config allele other than del type.\n";
    #}
}


# configure contig with breakpoints in breakdancer format (microhomology in the column after score, and templated insertion in the column after), with reference file given
# An example in ~/TIGRA-ext/Stephens_Data/use_config_fa.pl
sub config_fa_bed{
    my ($str, $reference, $out) = @_;
    my $h;
    my ($chr1, $start1, $end1, $chr2, $start2, $end2) = ($str =~ /^(\S+):(\d+)\-(\d+)\;(\S+):(\d+)\-(\d+)$/);
    my $fa1 = &make_contig_simple($reference, $chr1, $start1, $end1);
    my $fa2 = &make_contig_simple($reference, $chr2, $start2, $end2);
    $h->{$str} = &combine_contig($fa1, $fa2, "+", "-", -1, "NA");
    &write_fa($h, [keys %$h], $out);
}


# configure contig with breakpoints in breakdancer format (microhomology in the column after score, and templated insertion in the column after), with reference file given
# An example in ~/TIGRA-ext/Stephens_Data/use_config_fa.pl
sub config_fa{
    my ($breakdancer, $reference, $out, $flank, $micro_col, $templated_col) = @_;
    my $hash = breakdancer::read_breakdancer($breakdancer);
    my $h;
    foreach my $key (sort keys %$hash){
        my $line = $hash->{$key};
        my @lines = split(/:/, $line);
        # get the two breakpoints and their orientations
        my ($chr1, $bp1, $ori1, $chr2, $bp2, $ori2, $type) = @lines[0 .. 6];
        my ($micro, $templated) = ("", "");
        $micro = $lines[$micro_col] if($micro_col != -1);
        $templated = $lines[$templated_col] if($templated_col != -1);
        my $fa1 = &make_contig($reference, $chr1, $bp1, $ori1, $flank);
        my $fa2 = &make_contig($reference, $chr2, $bp2, $ori2, $flank);
        $DB::single = 1;
        $h->{$key} = &combine_contig($fa1, $fa2, $ori1, $ori2, $micro, $templated);
    }
    &write_fa($h, [keys %$hash], $out);
}

# make a contig according to breakpoint position, orientation and flanking region size
sub make_contig_simple{
    my ($ref, $chr, $start, $end) = @_;
    $start = 0 if($start < 0);
    return &faidx($ref, $chr, $start, $end, "contig_only");
}
# make a contig according to breakpoint position, orientation and flanking region size
sub make_contig{
    my ($ref, $chr, $pos, $ori, $flank) = @_;
    $DB::single = 1;
    my $end = $pos;
    $end = $pos + $flank if($ori =~ /\-/);
    $pos = $pos - $flank if($ori =~ /\+/);
    $pos = 0 if($pos < 0);
    return &faidx($ref, $chr, $pos, $end, "contig_only");
}

# make fastq file given sequence, and file_name, if no file_name, then make a temporary file in the current directory, and return the file name. qualities are all 2
# for multiple sequences, they are separated by ";"
sub make_fq{
    my ($seq, $file) = @_;
    if(!defined $file){
        my $time = tmpfile::time();
        $file = "seq.$time.fq";
    }
    my @seqs = split(/;/, $seq);
    my $num = 1;
    open fh_, ">$file" or die $!;
    foreach (@seqs){
        print fh_ join("\n", "@" . "read_$num", $_, "+", "2" x length($_)) . "\n";
        $num ++;
    }
    close fh_;

    return $file;
}

# from bam file, it makes a fastq file 
sub make_fq_from_bam{
    my ($bam, $file) = @_;
    `samtools view $bam | perl -ane \'print join("\n", "@" . \$F[0], \$F[9], "+", \$F[10]) ."\n" if(\$F[5] !~ /^\\d+H/ && \$F[5] !~ /\\d+H\$/)\' > $file`;
}

# from bam file, it makes a fastq file 
sub make_fa_from_bam{
    my ($bam, $file) = @_;
    `samtools view $bam | perl -ane \'print join("\n", ">" . \$F[0], \$F[9]) ."\n" if(\$F[5] !~ /^\\d+H/ && \$F[5] !~ /\\d+H\$/)\' > $file`;
}

# from bam file, make three fq files, two for paired, and one for unpaired, needs large memory
sub make_sorted_fq{
    my ($bam, $prefix) = @_;
    my $h;
    my @as = `samtools view $bam`;
    my $fq_1 = $prefix . "read_1.fq";
    my $fq_2 = $prefix . "read_2.fq";
    my $fq_single = $prefix. "read_single.fq";
    open pa_1, ">$fq_1" or die $!;
    open pa_2, ">$fq_2" or die $!;
    open fh_s, ">$fq_single" or die $!;

    foreach my $a (@as){
        my @b = split(/\t/, $a);
        # note: apply to read_name like ***/1, ***/2
        my $read_seq = substr $b[0], length($b[0])-1;
        my $read_name = substr $b[0], 0, length($b[0]) - 2;
        #$h->{$read_name}->{$read_seq} = "" if(!defined $h->{$read_name}->{$read_seq});
        $h->{$read_name}->{$read_seq} = $b[9] . "\n" . $b[10] . "\n";
    }
    foreach my $pair (keys %$h){
        my @read = sort {$a <=> $b} keys %{$h->{$pair}};
        if(scalar(@read) == 1){
            my @a = split(/\n/, $h->{$pair}->{$read[0]});
            print fh_s join("\n", "@" . $pair, $a[0], "+", $a[1]) . "\n";
        }
        elsif(scalar(@read) == 2){
            my @a = split(/\n/, $h->{$pair}->{$read[0]});
            print pa_1 join("\n", "@" . $pair . "/" . $read[0], $a[0], "+", $a[1]) . "\n";
            @a = split(/\n/, $h->{$pair}->{$read[1]});
            print pa_2 join("\n", "@" . $pair . "/" . $read[1], $a[0], "+", $a[1]) . "\n";
        }
    }
    close pa_1;
    close pa_2;
    close fh_s;
    print "Completed transformation from BAM to FQ: $fq_1, $fq_2, $fq_single\n";
}



# faidx of a pariticular position of reference, with mode contig_only that contains only the ACG string, not header
sub faidx{
    my ($ref, $chr, $pos, $end, $mode) = @_;
    my $result = `samtools faidx $ref $chr:$pos-$end`;
    my @results = split(/\n/, $result);
    my $str = "";
    foreach (@results){
        next if($_ =~ /^>/);
        $str .= $_;
    }
    return $str;
}

# combine the contigs with microhomology or tempalted insertion (non-templated), the orientations of the left and right are as follows: if ori1 is +, nothing to be changed (same for ori2 is -); if ori1 is -, complemented reverse to fa1 (same for ori2 is +)
sub combine_contig{
    my ($fa1, $fa2, $ori1, $ori2, $micro, $template) = @_;
    if($ori1 =~ /\-/){
        $fa1 = &reverse_complement($fa1);
    }
    if($ori2 =~ /\+/){
        $fa2 = &reverse_complement($fa2);
    }
    $template = "" if($template eq "NA");
    return $fa1 . $template . $fa2;
}

# reverse complement
sub reverse_complement{
    my ($str) = shift;
    my $rev_str = scalar reverse("$str");
    $rev_str =~ tr/ACGT/TGCA/;
    return $rev_str;
}

# separate by line_length
sub sep_line_length{
    my ($str) = shift;
    my $ret_str = "";
    while(length($str) > $line_length){
        $ret_str = $ret_str . substr($str, 0, $line_length) . "\n";
        $str = substr($str, $line_length);
    }
    $ret_str = $ret_str . $str . "\n";
    return $ret_str;
}
# extract contigs according to $contig_f, with g ahead, ignoring the string afterwards in fa
sub extract_contigs_first_col{
    my ($fa, $contig_f, $out) = @_;
    my $h;
    open out_fh, ">$out" or die $!;
    open ctg_f, "<$contig_f" or die $!;
    while(<ctg_f>){
        chomp;
        $h->{$_} = 1;
    }
    close ctg_f;
    open fa_, "<$fa" or die $!;
    my $tag = 0;
    while(<fa_>){
        chomp;
        $DB::single = 1;
        if($_ =~ /^>(g\d+)\./ && defined $h->{$1}){
            # print this line and open the printing for the next few lines with the sequence
            $tag = 1;
            print out_fh $_ . "\n";
        }
        elsif($_ !~ /^\>/ && $tag == 1){
            print out_fh $_ . "\n";
        }
        elsif($_ =~ /^\>/ && $tag == 1){
            $tag = 0;
        }
    }
    close fa_;
    close out_fh;
}

# extract multiple contigs with fa name stored in a file and output in fa
# only the first part before / should match the given fa name
sub extract_reads_PB{
    my ($fa, $contig_f, $out) = @_;
    my $h;
    open out_fh, ">$out" or die $!;
    open ctg_f, "<$contig_f" or die $!;
    while(<ctg_f>){
        chomp;
        $h->{">" .$_} = 1;
    }
    close ctg_f;
    open fa_, "<$fa" or die $!;
    my $tag = 0;
    while(<fa_>){
        chomp;
        if($_ =~ /^>/){
            my @a = split(/\//, $_);
            if(defined $h->{$a[0]}){
                # print this line and open the printing for the next few lines with the sequence
                $tag = 1;
                print out_fh $_ . "\n";
            }
            elsif($tag == 1){
                $tag = 0;
            }
        }
        elsif($_ !~ /^\>/ && $tag == 1){
            print out_fh $_ . "\n";
        }
    }
    close fa_;
    close out_fh;
}


# extract multiple contigs with fa name stored in a file and output in fa
# modified on 071315, added the third argument $out
sub extract_contigs{
    my ($fa, $contig_f, $out) = @_;
    my $h;
    open out_fh, ">$out" or die $!;
    open ctg_f, "<$contig_f" or die $!;
    while(<ctg_f>){
        chomp;
        $h->{">" .$_} = 1;
    }
    close ctg_f;
    open fa_, "<$fa" or die $!;
    my $tag = 0;
    while(<fa_>){
        chomp;
        if(defined $h->{$_}){
            # print this line and open the printing for the next few lines with the sequence
            $tag = 1;
            print out_fh $_ . "\n";
        }
        elsif($_ !~ /^\>/ && $tag == 1){
            print out_fh $_ . "\n";
        }
        elsif($_ =~ /^\>/ && $tag == 1){
            $tag = 0;
        }
    }
    close fa_;
    close out_fh;
}


# extract a contig with fa name and output in fa
sub extract_contig{
    my ($fa, $contig) = @_;
    open fa_, "<$fa" or die $!;
    my $tag = 0;
    while(<fa_>){
        if($_ =~ /^\>$contig/){
            # print this line and open the printing for the next few lines with the sequence
            $tag = 1;
            print $_;
        }
        elsif($_ !~ /^\>/ && $tag == 1){
            print $_;
        }
        elsif($_ =~ /^\>/ && $tag == 1){
            last;
        }
    }
    close fa_;
}

# extract a contig with fq name and output in fq
sub extract_contig_fq{
    my ($fq, $contig) = @_;
    open fq_, "<$fq" or die $!;
    my $tag = 0;
    while(<fq_>){
        if($_ =~ /^\@$contig/){
            # print this line and open the printing for the next few lines with the sequence
            $tag = 1;
            print $_;
        }
        elsif($tag != 0 && $tag < 4){
            print $_;
            $tag ++;
        }
        elsif($tag == 4){
            last;
        }
    }
    close fq_;
}

# remove duplicated fa name and sequence
sub dedup {
    my ($fa) = @_;
    my $fh;
    my $tag = 1;
    open fa_, "<$fa" or die $!;
    while(<fa_>){
        if($_ =~ /^>/){
            if(defined $fh->{$_}){
                $tag = 1;
                next;
            }
            else{
                $fh->{$_} = 1;
                $tag = 0;
                print $_;
            }
        }
        else{
            if($tag == 0){
                print $_;
            }
        }
    }
    close fa_;
}

# This concatenate a ls fasta files, changing the readnames to be fastafile_name_readname
sub concatenate_fa{
    my ($fas) = @_;
    my @fa = split(/\n/, `ls $fas`);
    foreach my $f(@fa){
        $DB::single = 1;
        my $file = `readlink -f $f`;
        chomp $file;
        open fh_, "<$file" or die $!;
        while(<fh_>){
            if($_ =~ /^>/){
                chomp;
                my ($name) = ($_ =~ /^>(.+)$/);
                print ">$f" . "_$name\n";
            }
            else{
                print $_;
            }
        }
        close fh_;
    }
}
# get reads from read names of IL and PB
# originally in union_find.pm, move it here so that the way to extract reads can be optimized. should use this one in future.
sub readname2fa_dir{
    my ($group, $short_fa_prefix, $long_fa_dir, $long_fa_name, $tag) = @_;
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
    my ($g) = ($group =~ /(out.g\d+\.txt)/);
    my $out1 = $g . ".PB.fa";
    my $out2 = $g . ".IL.1.fa";
    my $out3 = $g . ".IL.2.fa";
    my $h;
    # the hash to indicate which short read file to grep
    my $short_fa_1;
    my $short_fa_2;
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
                }
            }
        }
        else{
            $short_fa_1 = $short_fa_prefix . "_1.fa";;
            $short_fa_2 = $short_fa_prefix . "_2.fa";;
        }
    } # end of if of skip_il
    # read the read names and store to respective hash
    my $h_i;
    my $h_p;
    open g_h, "<$group" or die $!;
    while(<g_h>){
        chomp;
        if($skip_pb == 0 && $_ =~ /\./){
            $h_p->{$_} = 1;
        }
        if($skip_il == 0 && $_ !~ /\./){
            $h_i->{$_} = 1;
        }
    }
    close g_h;
    # assume only pair now
    $DB::single = 1;
    if($skip_il == 0){
        # read fa one by one and output the read when it is hit
        # get the maximum one, so that after it, the loop will jump out
        open OUT2, ">$out2" or die $!;
        open OUT3, ">$out3" or die $!;
        my @keys = sort {$b <=> $a} keys %$h_i;
        my $max_n = $keys[0];
        open fa_1, "<$short_fa_1" or die $!;
        my $n = 0;
        while(<fa_1>){
            if(defined $h_i->{$n}){
                print OUT2 $_;
                my $tmp = <fa_1>;
                print OUT2 $tmp;
                $n ++;
                last if($n > $max_n);
            }
            else{
                my $tmp = <fa_1>;
                $n ++;
            }
        }
        close fa_1;
        close OUT2;
        open fa_2, "<$short_fa_2" or die $!;
        $n = 0;
        while(<fa_2>){
            if(defined $h_i->{$n}){
                print OUT3 $_;
                my $tmp = <fa_2>;
                print OUT3 $tmp;
                $n ++;
                last if($n > $max_n);
            }
            else{
                my $tmp = <fa_2>;
                $n ++;
            }
        }
        close fa_2;
        close OUT3;
    }

    # deal with pb
    if($skip_pb == 0){
        open OUT1, ">$out1" or die $!;
        foreach my $pb (keys %$h_p){
            if($pb =~ /^(\d+)\.(\d+)/){
                my $seq = `samtools faidx ${long_fa_dir}$1/$long_fa_name $1.$2`;
                print OUT1 $seq;
            }
        }
        close OUT1;
    }
}

# insert string at the head/tail of readname, overwrite the original fasta. By default at tail.
sub change_readname{
    my ($fasta, $str, $pos) = @_;
    open fa, "<$fasta" or die $!;
    open fa_out, ">$fasta.tmp10000" or die $!;
    while(<fa>){
        if($_ =~ /^>(.+)$/){
            chomp;
            if($pos eq "head"){
                print fa_out join("", ">", $str, ".", $1) . "\n";
            }
            else{
                print fa_out join("", ">", $1, ".", $str) . "\n";
            }
        }
        else{
            print fa_out $_;
        }
    }
    close fa;
    close fa_out;
    `mv $fasta.tmp10000 $fasta`;
}

# get the substring from an fa according to SV call txt file
sub get_substr_by_call{
    my ($fa, $txt, $subset_txt) = @_;
    my $hh;
    if($subset_txt ne "NA"){
        open sub_fh, "<$subset_txt" or die $!;
        while(<sub_fh>){
            chomp;
            my @a = split(/\//, $_);
            # a[0] in the form of chr:start-end
            $hh->{$a[0]} = 1;
        }
        close sub_fh;
    }
    my $h;
    open txt_fh, "<$txt" or die $!;
    while(<txt_fh>){
        my @a = split(/\t/, $_);
        my @b = split(/\//, $a[4]);
        next if($subset_txt ne "NA" && !defined $hh->{"$b[0]:$a[5]-$a[6]"});
        $h->{"$b[0]:$a[5]-$a[6]"} = `samtools faidx $fa $b[0]:$a[5]-$a[6]`;
        print $h->{"$b[0]:$a[5]-$a[6]"};
    }
    close txt_fh;
}



1;
