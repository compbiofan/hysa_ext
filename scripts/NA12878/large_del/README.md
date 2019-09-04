# Supporting scripts for analysis of NA12878 large deletion.

## Steps

### Step 1. Run or download other programs. 

#### Run delly. 

Prepare exclusion file for delly parallel running by chromosome:

```for i in `seq 1 22` X Y; do env x=$i perl -F: -lane 'print $_ if($_ !~ /^chr$ENV{x}$/)' human.hg38.excl.tsv > excl_for_chr$i.tsv; done```

Generate script to run delly for each chromosome.

```for i in `seq 1 22` X Y; do echo "delly -t DEL -o del.chr$i.vcf.gz -g $ref -x excl_for_chr$i.tsv $bam" >> run_delly.sh; echo "#" >> run_delly.sh; done```

$bam is the short read bam aligned to the reference. $ref is the reference against which $bam is generated. For NA12878 we use build38. human.hg38.excl.tsv is included in delly package. The above script generated a bash script with each line the command to run delly separated by "#".  

Convert results to one bed file.

```for i in `seq 1 22` X; do gunzip del.chr$i.vcf.gz; done```

```cat del.chr*.vcf | grep -v "^#" | perl -ane 'next if($_ =~ /^#/); if($_ =~ /END=(\d+).+SR=(\d+)/){$end = $1; $SR = $2; if($_ =~ /0\/1/ || $_ =~ /1\/1/){print join("\t", $F[0], $F[1], $end, $SR) . "\n" if($end - $F[1] > 50)}}' > del.delly.bed```

#### Download pbhoney and custom from [1].

Download the spreadsheet from http://www.nature.com/nmeth/journal/v12/n8/extref/nmeth.3454-S3.xlsx. Convert to a txt format, and then bed format. 

```perl -ane 'print $_ if($F[3] == 1 && $F[2] - $F[1] > 50)' nmeth.3454-S3.txt | grep deletion > pbhoney.bed```

```perl -ane 'print $_ if($F[3] == 1 && $F[2] - $F[1] > 50)' nmeth.3454-S3.txt | grep custom > custom.bed```

#### Download svclassify.

Download the bed file from ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/technical/svclassify_Manuscript/Supplementary_Information/Personalis_1000_Genomes_deduplicated_deletions.bed. Lift over to build38.

```perl -ane 'print $_ if($_ !~ /^chr/)' Personalis_1000_Genomes_deduplicated_deletions.bed > svclassify.beforeLiftover.bed```
```liftover_hg19_hg38.sh svclassify.beforeLiftover.bed```
```ln -s svclassify.beforeLiftover.bed svclassify.bed```

### Step 2. Generate validated call set for each program by aligning the reconstructed allele to the whole genome assembly. 

Download whole genome assembly from ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_001013985.1_ASM101398v1/GCA_001013985.1_ASM101398v1_genomic.fna.gz [1] and gunzip it. 

```
for i in hysa delly pbhoney custom svclassify; 
do 
    bed=$i.bed
    rm bwa_ctg2asm.analysis;
    perl -M"fa" -ane '$l_s = $F[1] - 1000; $r_e = $F[2] + 1000; $bed_str = "$F[0]:$l_s-$F[1];$F[0]:$F[2]-$r_e"; fa::config_fa_bed($bed_str, $ref, "tmp.alt.ctg.fasta"); `bwa mem GCA_001013985.1_ASM101398v1_genomic.fna tmp.alt.ctg.fasta > bwa_ctg2asm.sam;`; `perl analyze_inserted_sequence.pl bwa_ctg2asm.sam >> bwa_allctg2asm.analysis`;' $bed
    less bwa_allctg2asm.analysis | grep "well_aligned" | perl -ane '@a = split(/;/, $F[1]); @b = split(/-/, $a[0]); @c = split(/-/, $a[1]); @d = split(/:/, $c[0]); print join("\t", $d[0], $b[1], $d[1]) . "\n"' > del.bwa_allctg2asm.analysis.well_aligned.txt
    perl -ane 'BEGIN{open fh_, "<del.bwa_allctg2asm.analysis.well_aligned.txt"; while(<fh_>){chomp; @F = split(/\s+/, $_); $h->{"$F[0]:$F[1]-$F[2]"} = 1}; close fh_} if(defined $h->{"$F[0]:$F[1]-$F[2]"}){print join("\t", @F, 1) . "\n"}else{print join("\t", @F, 0) . "\n"}' $bed > $bed.validate
    perl -ane 'print $_ if($F[$#F] == 1)' $bed.validate > $bed.validated
done
```

Note: HySA calls (hysa.bed) are the same call file namely out.del.GT50.bed from HySA progam. In the above for loop, $ref should be replaced with the reference (build38) fasta. 

Concatenate all validated calls from these five programs, and deduplicate them.

```
cat hysa.bed.validated delly.bed.validated pbhoney.bed.validated svclassify.bed.validated custom.bed.validated > all.bed
perl -ane 'print join("\t", $F[0], $F[1], $F[2]) . "\n"' all.bed > all.bed1
perl rm_dup_bed1.pl all.bed1 reciprocal50 > curated.bed
```

### Step 3. Evaluate the sensitivity and specificity for each program.

For delly and HySA, there are ROC curves on supporting soft clipped reads and supporting Illumina reads, respectively.

```
for p in hysa delly; do  
    for i in `seq 0 100`;
    do
        env x=$i perl -F: -lane '@F = split(/\s+/, $_); if($F[3] > $ENV{x}){print $_}' $p.bed > $p.$i.bed;
        intersectBed.reciprocal50.sh curated.bed $p.$i.bed 1 >> $p.sen.out
    done
    for i in `seq 0 100`;
    do
        intersectBed.v.sh $p.$i.bed curated.bed 1 >> $p.spe.out
    done
done
```

Make final ROC file with each row two columns (sensitivity and specificity) for both delly and HySA.

```
perl -ane 'BEGIN{open fh_, "<$p.sen.out"; while(<fh_>){chomp; push @a, $_;}close fh_; $n = 0}print join("\t", $a[$n], $_); $n ++;' $p.sep.out > sen_spe.csv
less sen_spe.csv | perl -ane '$F[0] /= $n_curate; $F[1] = 1 - $F[1]/$n_p; print join("\t", $F[0], $F[1]) . "\n"' > $p.roc.csv
```

Note: in the above two commands, replace $p with delly and hysa subsequently, and $n_p with total large deletion call from delly and HySA, respectively. $n_curate is the total number of curated calls in curated.bed for both programs.

For pbhoney, custom and svclassify, there are only three dots on the comparison figure. 

```
for p in pbhoney custom svclassify; do
    echo "$p:"
    echo "total $p calls:"
    wc -l $p.bed
    echo "overlap with gs.bed:"
    intersectBed.reciprocal50.sh curated.bed $p.bed 1
    echo "false positive:"
    intersectBed.v.sh $p.bed curated.bed 1
done
```

Suppose $a, $b and $c are the three numbers obtained from the above for loop for pbhoney, custom and svclassify. Then sensitivity = $b/$n_curate; specificity = 1- $c/$a.  

Generate three files (pbhoney.roc.csv, custom.roc.csv and svclassify.roc.csv) for pbhoney, custom and svclassify, each with one row (two columns separated by tab, the first column the sensitivity and the second the specificity of the program). 

### Step 4. Draw ROC curves. 

Launch R. Run

`source $dir/scripts/R_scripts/plot_roc.r`

in which $dir is the directory where this package is installed. 

The figure "roc.png" is produced.

## Reference

1. Pendleton M, Sebra R, Pang AW, Ummat A, Franzen O, Rausch T, Stutz AM, Stedman W, Anantharaman T, Hastie A et al. 2015. Assembly and diploid architecture of an individual human genome via single-molecule technologies. Nat Methods 12: 780-786.




