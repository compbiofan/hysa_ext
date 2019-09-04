# Supporting scripts for analysis of CHM1 small deletion and insertion.

## Steps

### Step 1. Prepare gatk results

Prepare CHM1 bam with read groups.

```samtools reheader build37.head.wRG CHM1.bam > CHM1.wHeader.bam```
```samtools view -h CHM1.wHeader.bam | perl -ane 'chomp; print $_ . "\n" if($_ =~ /^@/); print join("\t", $_, "RG:Z:0001") . "\n" if($_ !~ /^@/);' | samtools view -bhS - > CHM1.withRG.bam```
```samtools index CHM1.withRG.bam```

CHM1.bam is the original bam. "build37.head.wRG" is provided in this folder. This command generates a bam that has the read group for each line.

Generate GATK command on each chromosome:

```for file in `seq 1 22` X Y; do echo "java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R $ref -I $bam  --genotyping_mode DISCOVERY -stand_emit_conf 10  -stand_call_conf 30 -o raw.variants.chr$file.vcf -L $file" >> run_gatk.sh; echo "java -Xmx2g -jar GenomeAnalysisTK.jar -T SelectVariants -R $ref -V raw.variants.chr$file.vcf -selectType INDEL -o raw.indels.chr$file.vcf" >> run_gatk.sh; echo "java -Xmx2g -jar GenomeAnalysisTK.jar -T VariantFiltration -R $ref -V  raw.indels.chr$file.vcf --filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\" --filterName \"indel_filter\" -o filtered.indels.chr$file.vcf" >> run_gatk.sh ; echo "#" >> run_gatk.sh; done```

$bam is CHM1.withRG.bam. $ref is the reference; here we use build37. 

Then run commands in run_gatk.sh. You can parallelize the tasks by chromosome, separated by "#" in this bash script.

Summarize the filtered results in bed format and use only those with size < 50bp:

```cat filtered.indels.chr*.vcf | perl -ane 'next if($_ =~ /^#/); next if(length($F[4] > length($F[3]))); $len = length($F[3]); print join("\t", $F[0], $F[1], ($F[1] + $len), $len) . "\n" if($len < 50)' > gatk.deletion.bed
```cat filtered.indels.chr*.vcf | perl -ane 'next if($_ =~ /^#/); next if(length($F[4] < length($F[3]))); $len = length($F[4]); print join("\t", $F[0], ($F[1] - 50), ($F[1] + 50), $len) . "\n" if($len < 50)' > gatk.insertion.bed

### Step 2. Prepare pindel results. 

Gerenate config file:
```echo "$bam $insert_size    $tag" > $cfg``` 
```for chr in `seq 1 22` X Y; do echo "pindel -f $build37 -i $cfg -c 1 -T 4 -o out_chr$chr" >> run_pindel.sh; echo "#" >> run_pindel.sh; done```

$bam is CHM1.bam. $insert_size is the average insert size. We use breakdancer bam2cfg to calculate it: 

```bam2cfg.pl $bam```

Average insert size is the number following "mean:" in the output of this command. $tag is the name of the sample. We use "short" here. 

Then run commands in run_pindel.sh. You can parallelize the tasks by chromosome, separated by "#" in this bash script.

Summarize the results after filtering and use only those with size (10, 50]bp:

```for i in `seq 1 22` X Y; do grep ChrID pindel0.2.5_default_chr${i}_D | perl -ane 'print join("\t", $F[7], $F[9], $F[10], ($F[10] - $F[9])) . "\n" if($F[2] > 10 && $F[15] >= 5 && $F[10] - $F[9] <= 50)' >> pindel.deletion.bed; done```
```for i in `seq 1 22` X Y; do grep ChrID pindel0.2.5_default_chr${i}_SI | perl -ane 'print join("\t", $F[7], ($F[9] - 50), ($F[10] + 50), $F[2]) . "\n" if($F[2] > 10 && $F[15] >= 5 && $F[2] <= 50)' >> pindel.insertion.bed; done``` 

Note: insertion calls are with 50bp two-sided padding.

### Step 3. Make venn diagrams for small deletion and insertion. 

```calculate_mendelian_noGT.pl pindel.deletion.bed gatk.deletion.bed out.del.SET50.rmdup.bed touch```
```calculate_mendelian_noGT.pl pindel.insertion.bed gatk.insertion.bed out.ins.SET50.flanking50.rmdup.bed touch```

out.del.SET50.rmdup.bed and out.ins.SET50.flanking50.rmdup.bed are the small deletion and insertion files produced by HySA.

### Step 4. Overlap our unique calls with database (dbSNP and short tandem repeat).

In step 3, the files of our unique calls are generated, namely, "out.del.SET50.rmdup.bed.only" for deletion and "out.ins.SET50.flanking50.rmdup.bed.only" for insertion.

Download the dbSNP file from 

http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/snp146.txt.gz

Convert it to bed format:

```gunzip snp146.txt.gz```
```perl -ane 'if($_ =~ /deletion/){print join("\t", @F[1 .. 3]) . "\n"}' snp146.txt > snp146.deletion.bed```
```perl -ane 'if($_ =~ /insertion/){print join("\t", $F[1], ($F[2] - 50), ($F[2] + 50)) . "\n"}' snp146.txt > snp146.insertion.bed```

Download the short tandem repeat file from

http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz

Convert it to bed format:

```gunzip simpleRepeat.txt.gz```
```perl -ane '$F[1] =~ s/chr//; print join("\t", @F[1 .. 3]) . "\n"' simpleRepeat.txt > simpleRepeat.bed```

Overlap our unique small indels with the two databases:

```intersectBed.touch.sh out.del.SET50.rmdup.bed.only $STR 1```
```intersectBed.touch.sh out.del.SET50.rmdup.bed.only $dbSNP_del 1```
```intersectBed.touch.sh out.ins.SET50.flanking50.rmdup.bed.only $STR 1```
```intersectBed.touch.sh out.ins.SET50.flanking50.rmdup.bed.only $dbSNP_ins 1```

$STR is simpleRepeat.bed. $dbSNP_del is snp146.deletion.bed. $dbSNP_ins is snp146.insertion.bed.

The four numbers shown from these four commands are the number of HySA unique small deletion calls overlapping with short tandem repeat database, dbSNP database, and the number of HySA unique small insertion calls overlapping with short tandem repeat database and dbSNP database. 

To calculate the number of deletion and insertion calls overlapping both database, run the following command:

```cat snp146.insertion.bed simpleRepeat.bed | perl -ane 'chomp; print join("\t", @F[0 .. 2]) . "\n"' > database.ins.bed```
```intersectBed.touch.sh out.ins.SET50.flanking50.rmdup.bed.only database.ins.bed 1```

```cat snp146.deletion.bed simpleRepeat.bed | perl -ane 'chomp; print join("\t", @F[0 .. 2]) . "\n"' > database.del.bed```
```intersectBed.touch.sh out.del.SET50.rmdup.bed.only database.del.bed 1```

