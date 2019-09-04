# Supporting scripts for analysis of NA12878 complex deletions.

## Steps

### Step 1. Run HySA on complex mode, with the output in complex.txt. 

### Step 2. Prepare the bed file with complex deletions whose inserted sequence has the origin. 

```
perl -ane 'chomp; @b = split(/;/, $F[7]); next if(@b < 2); @c = split(/:/, $b[0]); print join("\t", $F[0], $F[1], ($F[1] + $F[2]), $F[2], $c[1]); for(my $i = 1; $i < scalar(@b); $i ++){@a = split(/:/, $b[$i]); $ori = "normal"; $type = "duplication"; if($a[0] > $a[1]){$ori = "inverted";}$del_end = $F[1] + $F[2];if($a[3] <= $del_end && $a[3] >= $F[1] || $a[4] <= $del_end && $a[4] >= $F[1]){$type = "spacer";} $location="$a[2]:$a[3]-$a[4]"; print "\t". join("\t", $location, $type, $ori)}print "\n"' complex.txt | sort -k1,1 > complex_w_origin.bed;
```

This generates Supplemental Table 5 sheet "Complex_DEL (DEL+ INS_w_origin)". 

### Step 3. Prepare the bed file with all complex deletions.

```
perl -ane 'chomp; @b = split(/;/, $F[7]); @c = split(/:/, $b[0]); print join("\t", $F[0], $F[1], ($F[1] + $F[2]), $F[2], $c[1]) . "\n"' complex.txt | sort -k1,1 > complex_all.bed;
```

This generates Supplemental Table 5 sheet "All_Complex_DEL (DEL+INS)". 

### Step 4. Forsmid validation of complex deletions. fosmid.covered.bed and fosmid.sam are the same files generated in ../fosmid_validation/

Prepare bed files that overlap with fosmid covered region.

```
bed=complex_all_overlap_w_fosmid.bed;
~/bin/intersectBed.touch.sh complex_all.bed fosmid.covered.bed > $bed;
~/bin/intersectBed.touch.sh fosmid.covered.bed complex_all.bed > $bed.tmp;
```

Prepare fosmid reads whose alignment overlap with complex deletions.

```
sam=fosmid.sam;
perl -ane 'BEGIN{open fh_, "<complex_all_overlap_w_fosmid.bed.tmp" or die $!;while(<fh_>){@a = split(/\t/, $_); $h->{$a[0]}->{$a[1]} = 1;}close fh_; } next if($_ =~ /^@/); if(defined $h->{$F[2]}->{$F[3]}){print join("\n", ">$F[0]", $F[9]) . "\n"}' $sam > fosmid.fa;
~/pkg/bwa/bwa index fosmid.fa;
```

Prepare the alternative alleles of these complex deletions that overlap with fosmid covered regions and align them to fosmid.fa.

```
rm bwa_allctg2asm.analysis;
perl -ane ' $ins_col = 4; $ref = "~/reference/build38.fa"; $l_s = $F[1] - 50; $r_e = $F[2] + 50; $s = $F[1] - 1; $e = $F[2] + 1; $bed_str = "$F[0]:$l_s-$s"; @a = split(/\n/, `samtools faidx $ref $bed_str`); $str = join("", @a[1 .. $#a]) . $F[$ins_col]; $bed_str = "$F[0]:$e-$r_e"; @a = split(/\n/, `samtools faidx $ref $bed_str`); $str = $str . join("", @a[ 1 .. $#a]); $name = "$F[0]:$l_s-$F[1];$F[0]:$F[2]-$r_e"; open fh_, ">tmp.alt.ctg.fasta"; print fh_ join("\n", ">$name", $str); close fh_; `bwa mem fosmid.fa tmp.alt.ctg.fasta > bwa_ctg2asm.sam;`; `perl ~/combinePBIL/scripts/pp/supporting_scripts/analyze_deleted_sequence.pl bwa_ctg2asm.sam 0.9 9 9 >> bwa_allctg2asm.analysis`;' $bed;
```

Find out which complex deletions have good alignment to fosmid.fa.

```
perl -ane 'BEGIN{open fh_, "<bwa_allctg2asm.analysis"; while(<fh_>){next if($_ !~ /well_aligned/); @a = split(/\t/, $_); foreach $i (@a[1 .. $#a]){($c1, $b1, $e1, $c2, $b2, $e2) = ($i =~ /^(\S+):(\d+)-(\d+);(\S+):(\d+)-(\d+)$/); $h->{"$c1:$e1-$b2"} = 1; }} close fh_; } if(defined $h->{"$F[0]:$F[1]-$F[2]"}){print join("\t", @F, 1) . "\n"}else{print join("\t", @F, 0) . "\n"}' $bed > $bed.validate;
perl -ane 'print $_ if($F[$#F] == 1)' $bed.validate > $bed.validated;
```

Final results are in complex_all_overlap_w_fosmid.bed.validated. 
