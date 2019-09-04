use warnings;
use strict;

# Deal with bam files.
require sample_list;
require cigar;
package bam;

my @chrs = (1 .. 22, 'X');

sub extract_reads{
    my ($bam, $chr, $pos, $flank, $out) = @_;
    my ($start, $end) = ($pos - $flank, $pos + $flank);
    my $str = "$chr:$start-$end";
    `samtools view -bh $bam $str > $out`;
    `samtools index $out`;
}

sub bam2fa{
    my ($bam) = @_;
    my ($fa) = ($bam =~ /^(.+)\.bam/);
    $fa .= ".fa";
    open BAM, "samtools view $bam |" or die $!;
    open fa_out, ">$fa" or die $!;
    while(<BAM>){
        next if($_ =~ /^\@/);
        my @a = split(/\t/, $_);
        print fa_out join("\n", ">" .$a[0], $a[9])."\n";
    }
    close BAM;
    close fa_out;
}
sub sam2bam{
    # write a bam from sam
    my ($sam) = @_;
    my ($bam) = ($sam =~ /^(.+)\.sam/); 
    $bam .= ".bam";
    `samtools view -bhS $sam > $bam`;
    `samtools sort $bam $bam.sorted`;
    `mv $bam.sorted.bam $bam.bam`;
    `samtools index $bam.bam`;
}

# the following is to avoid merging to a big bam, and thus extract reads from small bams
sub extract_reads_mBPs_smallBams{
    my ($bam_list, $file, $flank, $out) = @_;
    open bam_fh, "<$bam_list" or die $!;
    while(<bam_fh>){
        chomp;
        &extract_reads_mBPs($_, $file, $flank, $out, 3);
    }
    close bam_fh;
}


sub extract_reads_mBPs{
    my ($bam, $file, $flank, $out, $out_mode) = @_;
    # $file in the format of chr1.start.chr2.end.TYPE or in bed format if end up with .bed
    # generate head
    # The out_mode is the output file format. If 0, then each line is a read name. If 1, then in sam. If 2, then in bam. If 3, then in fa, but do not erase the previous content in the same file name.
    my $prefix = $out;
    open fh_, "<$file" or die $!;
    # default type if not present in bed
    my $type = "DEL";
    while(<fh_>){
        next if($_ =~ /^#/);
        #|| $_ !~ /\./ && $file !~ /\.bed/);
        # output to each separate file
        chomp;
        my @a;
        if($file =~ /\.bed/){
            my ($chr, $start, $end) = split(/\s+/, $_);
            @a = ($chr, $start, $chr, $end);
        }
        else{
            @a = split(/\./, $_);
            next if($#a < 3);
        }
        $out = "$prefix." . join(".", @a);
        `rm $out.txt` if(-e "$out.txt" && $out_mode == 0);
        `rm $out.sam` if(-e "$out.sam" && $out_mode == 1);
        `rm $out.bam` if(-e "$out.bam" && $out_mode == 2);
        if($out_mode == 1 || $out_mode == 2){
            `samtools view -H $bam > $out.sam`;
        }

        if($a[0] eq $a[2] && abs($a[3] - $a[1]) < 2 * $flank){
            my ($s, $e) = ($a[1] - $flank, $a[3] + $flank);
            my $str = "$a[0]:$s-$e";
            if($out_mode == 0){
                `samtools view $bam $str | perl -ane 'print ">" . \$F[0] . "\n"' >> $out.txt`;
            }
            elsif($out_mode == 1 || $out_mode == 2){
                `samtools view $bam $str >> $out.sam`;
            }
            elsif($out_mode == 3){
                # to extract reads from multiple bams but to avoid duplication, assuming each bam contains different set of reads
                my $r = `samtools view $bam $str`;
                &write_wo_dup($r, $out, $out_mode);
            }
        }
        else{
            my ($s, $e) = ($a[1] - $flank, $a[1] + $flank);
            my $str = "$a[0]:$s-$e";
            if($out_mode == 0){
                `samtools view $bam $str | perl -ane 'print ">" . \$F[0] . "\n"' >> $out.txt`;
            }
            elsif($out_mode == 1 || $out_mode == 2){
                `samtools view $bam $str >> $out.sam`;
            }
            elsif($out_mode == 3){
                my $r = `samtools view $bam $str`;
                &write_wo_dup($r, $out, $out_mode);
            }
            ($s, $e) = ($a[3] - $flank, $a[3] + $flank);
            $str = "$a[2]:$s-$e";
            if($out_mode == 0){
                `samtools view $bam $str | perl -ane 'print ">" . \$F[0] . "\n"' >> $out.txt`;
            }
            elsif($out_mode == 1 || $out_mode == 2){
                `samtools view $bam $str >> $out.sam`;
            }
            elsif($out_mode == 3){
                my $r = `samtools view $bam $str`;
                &write_wo_dup($r, $out, $out_mode);
            }
        }

    }
    close fh_;
    if($out_mode == 2){
        # write the sam to bam
        `samtools view -bhS $out.sam > $out.bam`;
        `samtools sort $out.bam $out.sorted`;
        `mv $out.sorted.bam $out.bam`;
        `samtools index $out.bam`;
        print "Done: extracting reads and make to a bam file: $out.bam\n";
    }
    elsif($out_mode == 1){
        print "Done: extracting reads and make to a sam file: $out.sam\n";
    }
    elsif($out_mode == 0){
        print "Done: extracting reads and make to a file with each line the read name: $out.txt\n";
    }
    elsif($out_mode == 3){
        print "Done: extracting reads without writing over the original file, result in a fasta file.\n";
    }

}

sub write_wo_dup{
    # write the extracted reads with out_mode, to out, without writing the duplicated reads
    my ($r, $out, $out_mode) = @_;
    if($out_mode == 3){
        my @F = split(/\n/, $r);
        my $ff_h;
        open out_fa, ">>$out.fa" or die $!;
        foreach my $f (@F){
            my @ff = split(/\t/, $f);
            if(!defined $ff_h->{$ff[0]}){
                # no duplication
                print out_fa ">" . $ff[0] . "\n$ff[9]\n";
                $ff_h->{$ff[0]} = 1;
            }
        }
        close out_fa;
    }
}

# extract reads in a more comprehensive way
sub extract_reads_comp{
    my ($bam, $pos, $flank, $out) = @_;
    my @poss = split(/;/, $pos);
    if(scalar(@poss) != 0){
        $DB::single = 1;
        `samtools view -H $bam > $out.sam`;
    }
    foreach (@poss){
        if($_ =~ /\-/ && $_ =~ /:/){
            `samtools view $bam $_ >> $out.sam`;
        }
        elsif($_ =~ /:/){
            my @a = split(/:/, $_);
            my ($start, $end) = ($a[1] - $flank, $a[1] + $flank);
            my $str = "$a[0]:$start-$end";
            `samtools view $bam $str >> $out.sam`;
        }
    }
    `samtools view -bS $out.sam > $out.bam`;
    `samtools sort $out.bam $out.sorted`;
    `samtools index $out.sorted.bam`;
}

# convert a bam to a position file for BIC-seq
sub prepare_BIC_seq{
    my ($t_bam, $n_bam, $prefix) = @_;
    `mkdir $prefix/tumor` if(!-d "$prefix/tumor");
    `mkdir $prefix/normal` if(!-d "$prefix/normal");
    foreach my $chr (@chrs){
        `samtools view $t_bam $chr | perl -ane \'print \$F[3] ."\n" if(\$F[4] ne "0")\' > $prefix/tumor/chr$chr.seq`;
        `samtools view $n_bam $chr | perl -ane \'print \$F[3] ."\n" if(\$F[4] ne "0")\' > $prefix/normal/chr$chr.seq`;
    }
    print "Complete preparation of BIC-seq.\n";
}

sub group_reads{
    # group those paired end reads as 0, and the others as 1
    # if insert, input insert size as the last input
    # if CTX, input the two chrs
    # if one breakpoint and manual, input the two connected chrs
    my ($bam, $chr, $pos, $flank, $out, $mode, $insert, $pos_range) = @_;
    my @a;
    if($chr ne "" && $pos != -1){
        my ($start, $end) = ($pos - $flank, $pos + $flank);
        my $str = "$chr:$start-$end";
        @a = `samtools view -h $bam $str`;
    }
    else{
        @a = `samtools view -h $bam`;
    }
    my $sam = $out . ".tmp.sam";
    my @chrs;
    my $hash;
    my $col = 6;
    #limiting the position of higher group reads
    my ($pos1, $pos2);
    if($mode eq "CTX"){
        # $insert as the connected chr for CTX mode
        # separated by :
        @chrs = split(/:/, $insert);
    }
    elsif($mode eq "manual"){
        @chrs = split(/:/, $insert);
        for(my $i = 0; $i < scalar(@chrs); $i ++){
            $hash->{$i} = 0;
        }
        ($pos1, $pos2) = split(/:/, $pos_range) if($pos_range ne "");

    }

    open TMP_SAM, ">$sam" or die $!;
    foreach (@a){
        if($_ =~ /^@/){
            print TMP_SAM $_;
            next;
        }
        my @b = split(/\t/, $_);
        my $RG = 1;
#        $DB::single = 1;
        if($mode eq "insert" && abs($b[8]) > $insert){
            $RG = 0;
        }
        elsif($mode eq "CTX"){
            foreach my $c (@chrs){{
                    if($b[6] eq $c){
                        $RG = 0;
                    }
                }
            }
        }
        elsif($mode eq "manual"){
            #for display a breakpoint connecting to two chromosomes, for IGV, limited number of reads so taht they don't stack deeply
            my $tag = 0;
            for(my $i = 0; $i < scalar(@chrs); $i ++){
                if($b[$col] eq $chrs[$i] && $hash->{$i} < 100 && ($pos_range eq "" || $pos_range ne "" && $b[3] > $pos1 && $b[3] < $pos2)){
                    $DB::single = 1;
                    $RG = 0;
#                    $RG = $i;
                    $tag = 1;
                    $hash->{$i} ++;
                }
            }
#            $RG = scalar(@chrs) if($tag == 0);
        }


        for(my $i = 11; $i < $#b; $i ++){
            if($b[$i] =~ /^RG:Z:/){
                print TMP_SAM join("\t", @b[0 .. $i-1], "RG:Z:$RG", @b[$i + 1 .. $#b]);
            }
        }
    }
    close TMP_SAM;

    `samtools view -bhS $sam > $out`;

    `samtools index $out`;
}


# count reads given a bam and a region for 1000G
sub count_reads_1KG{
    my ($sample, $region) = @_;

    my $file = "/scratch/bcb/xfan3/pkg/use_perl_modules/Sample_bam.lst";
    my $hash_sample;
    $hash_sample->{$sample} = 1;
    my $hash = sample_list::get_samples($hash_sample, $file);

    my $bam = $hash->{$sample};
    my $count = `samtools view -c $bam $region`;
    chomp $count;
    return $count;
}

sub samtools_count{
    my ($bam, $region) = @_;
    my $count = `samtools view -c $bam $region`;
    return $count;
}
sub samtools_rdname{
    my ($bam, $region, $h) = @_;
    my $read = `samtools view $bam $region`;
    my @reads = split(/\n/, $read);
    foreach my $r (@reads){
        my @b = split(/\t/, $r);
        $h->{$b[0]} = 1;
    }
    return $h;
}

# given the breakpoint and the end of the clipping (- for 3' clipping, + for 5' clippign), return the alignments that support such breakpoint, allowing some offset
sub extract_soft_clip{
    my ($bam, $region, $bp, $ori, $off) = @_;
    my @alns = split(/\n/, `samtools view $bam $region`);
    # a hash to store the bps
    my $hash;
    my $max = 0;
    my $ret_bp = -1;
    foreach my $aln (@alns){
        my @a = split(/\t/, $aln);
        next if($a[4] == 0);
        my $new_bp = cigar::check_breakpoint($a[5], $a[3], $bp, $ori, $off);
        next if($new_bp == 0);
        if(!defined $hash->{$new_bp}){
            $hash->{$new_bp} =  1;
        }
        else{
            $hash->{$new_bp} += 1;
        }
        if($max < $hash->{$new_bp}){
            $max = $hash->{$new_bp};
            $ret_bp = $new_bp;
        }
    }
    return ($max, $ret_bp);
}


1;
