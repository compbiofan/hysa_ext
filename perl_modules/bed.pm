use warnings;
use strict;

package bed;

my $home = "/scratch/bcb/xfan3/pkg/bedtools/bedtools-2.17.0/bin";

# Function: This package manipulates bed formatted file.

sub rm_redundancy_l_r{
    my ($file) = @_;
    my $h;
    my $h1;
    open fh_, "<$file" or die $!;
    while(<fh_>){
        my @F = split(/\t/, $_);
        $h->{$F[1]}->{"$F[2].$F[3]"} = $F[0];
        $h1->{$F[0]} = $_; 
    }
    close fh_;

    my $h_union;
    foreach my $chr (keys %$h){
        foreach my $s1 (sort {$a <=> $b} keys %{$h->{$chr}}){
            foreach my $s2 (sort {$a <=> $b} keys %{$h->{$chr}}){
                next if($s1 == $s2);
                if(abs($s1 - $s2) < 500){
                    $DB::single = 1;
                    my $n1 = $h->{$chr}->{$s1};
                    my $n2 = $h->{$chr}->{$s2};
                    $h_union = &check_redundancy($h_union, $h1->{$n1}, $h1->{$n2});
                }
            }
        }
    }

    foreach my $n (sort {$a <=> $b} keys %$h1){
        if(defined $h_union->{$n} && $h_union->{$n} != $n){
            print "#$h_union->{$n}" . "\t" . $h1->{$n};
            next;
        }
        else{
            print $h1->{$n};
        }
    }
}

# if redundant, point from the redundant one to the true one
sub check_redundancy{
    my ($h, $l1, $l2) = @_;
    my @a = split(/\t/, $l1);
    my @b = split(/\t/, $l2);
    if(!defined $a[6] || !defined $a[5]){
        print join("\n", ":$l1", ":$l2");
        return;
    }
    if($a[1] eq $b[1] && $a[4] eq $b[4]){
        if($a[2] <= $b[3] && $a[2] >= $b[2] || $b[2] <= $a[3] && $b[2] >= $a[2]){
            if($a[5] <= $b[6] && $a[5] >= $b[5] || $b[5] <= $a[6] && $b[5] >= $a[5]){
                # get the root
                if(!defined $h->{$a[0]}){
                    $h->{$a[0]} = $a[0];
                }
                $h->{$b[0]} = $h->{$a[0]};
            }
        }
    }
    return $h;
}





sub intersect_bed{
    # intersect two bed files with option and write to output_file
    my ($a, $b, $of, $option) = @_;
    `$home/intersectBed -a $a -b $b $option > $of`;
} # end of intersect_bed

sub write_bed{
    # Given a hash with "chr1:pos1:chr2:pos2" -> line, write a bed file with a certain range of flanking window. If chr1 eq chr2, then the whole pos1-pos2 is included. Otherwise, to two breakpoints separately as two lines.
    my ($hash, $of, $ob) = @_;
    my ($l, $r);
    print $of;
    open of_fh, ">$of" or die $!;

    foreach my $key (keys %$hash){
        my @a = split(/:/, $key);
        if($a[0] eq $a[2]){
            next if($a[1] < 0 || $a[1] > $a[3]);
            $l = $a[1] - $ob;
            $l = 0 if($l < 0);
            $r = $a[3] + $ob;
            print of_fh join("\t", ($a[0], $l, $r, $hash->{$key}))."\n";
        }
        else{
            next if($a[1] < 0 || $a[3] < 0);
            $l = $a[1] - $ob;
            $l = 0 if($l < 0);
            $r = $a[1] + $ob;
            print of_fh join("\t", ($a[0], $l, $r, $hash->{$key}))."\n";
            $l = $a[3] - $ob;
            $l = 0 if($l < 0);
            $r = $a[3] + $ob;
            print of_fh join("\t", ($a[2], $l, $r, $hash->{$key}))."\n";
        }
    }
    close of_fh;
}

sub bed2breakdancer{
    # convert bed to breakdancer for SVs (use DEL as default type if not otherwise indicated in the column or the 2nd argument (0-based). 
    # If the 2nd argument (0-based) is a number, take that as the column index to find type. The number is 0-based. Otherwise, take that string as the SV type for the whole bed file's SVs.
    my ($bed_file, $bd_file, $type) = @_;
    my $sv_type = "DEL";
    my $sv_type_col = -1;
    if(defined $type){
        if($type =~ /^\d+$/){
            $sv_type_col = $type;
        }
        else{
            $sv_type = $type;
        }
    }
    open out_bd, ">$bd_file" or die $!;
    open out_bed, "<$bed_file" or die $!;
    while(<out_bed>){
        my @a = split(/\t/, $_);
        if($sv_type_col != -1){
            $sv_type = $a[$sv_type_col];
        }
        print out_bd join("\t", ($a[0], $a[1], "+", $a[0], $a[2], "-", $sv_type, $a[2] - $a[1], "99"))."\n"; 
    }
    close out_bd;
    close out_bed;
}


sub bed2vcf_1kg{
    # convert bed to vcf for SVs (use DEL as default type if not indicated in the fourth column (1-based)
    my ($bed_file, $vcf_file, $type) = @_;
    # write header to vcf_file, copied from Zechen's code
    open out_fh, ">$vcf_file" or die $!;
    print out_fh qq(##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion">
##FORMAT=<ID=SPNUM,Number=1,Type=Integer,Description="Number of supporting Illumina reads.">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
);
open bed_fh, "<$bed_file" or die $!;
while(<bed_fh>){
    next if($_ =~ /^#/);
    chomp;
    my @a = split(/\t/, $_);
    $a[0] =~ s/^chr//;
    if($type eq "INS" || $bed_file =~ /flanking50/){
        $a[1] += 50;
        $a[2] -= 50;
    }
    my $SV_len = $a[3];
    print out_fh join("\t", @a[0 .. 1], ".", "N", $type, ".", "PASS", "CIPOS=-50,50;CIEND=-50,50;SVTYPE=$type;END=$a[2];SVLEN=$SV_len", "SPNUM", $a[$#a]) . "\n";
}
close out_fh;
close bed_fh;
}


sub bed2vcf{
    # convert bed to vcf for SVs (use DEL as default type if not indicated in the fourth column (1-based)
    my ($bed_file, $vcf_file) = @_;
    # write header to vcf_file, copied from Zechen's code
    open out_fh, ">$vcf_file" or die $!;
    print out_fh qq(##fileformat=VCFv4.1
##phasing=none
##INDIVIDUAL=TRUTH
##SAMPLE=<ID=TRUTH,Individual="TRUTH",Description="bamsurgeon spike-in">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation in primary">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">
##INFO=<ID=UNKNOWN_LEN,Number=0,Type=Flag,Description="Unknown the length of SV">
##INFO=<ID=SOURCE1,Number=1,Type=String,Description="Paper where this call comes from.">
##INFO=<ID=SOURCE2,Number=1,Type=String,Description="Method that generates this call.">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=CTX,Description="Translocation">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SPIKEIN
);
open bed_fh, "<$bed_file" or die $!;
while(<bed_fh>){
    chomp;
    my @a = split(/\t/, $_);
    $a[0] =~ s/^chr//;
    my $SV_len = $a[2] - $a[1];
    print out_fh join("\t", @a[0 .. 1], "N", ".", "<$a[3]>", 30, "PASS", "PRECISE;CIPOS=-10,10;CIEND=-10,10;SVTYPE=$a[3];CHR2=$a[0];END=$a[2];SVLEN=$SV_len");
    # paper source
    print out_fh ";SOURCE1=$a[5]" if(defined $a[5]);
    # method source
    print out_fh ";SOURCE2=$a[4]" if(defined $a[4]);
    print out_fh "\n";
}
close out_fh;
close bed_fh;
}

1;

