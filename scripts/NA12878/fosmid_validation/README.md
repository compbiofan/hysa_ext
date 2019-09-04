# Supporting scripts for validation with fosmid data 

## Steps

### Step 1. Align fosmid reads to NCBI build v38.

```
bwa mem build38.fa NA12878.fosmid.ABC12.Q80.cleaned.fastq > fosmid.sam;
makebam.2.sh fosmid;
```

### Step 2. Get the range of fosmid covered region.

```
cover=fosmid.covered.bed;
samtools view fosmid.sorted.bam | perl -ane 'next if($F[1] & 0x100 || $F[1] & 0x400 || $F[1] & 0x800); $chr = $F[2]; $pos = $F[3]; $end = $pos; $cigar = $F[5]; while($cigar ne ""){($num, $tag, $cigar) = ($cigar =~ /^(\d+)(\S)(.*)/); if($tag =~ /[MD]/){ $end += $num; }} print join("\t", $chr, $pos, $end) . "\n"' > $cover;
```

### Step 3. Generate validated call set from the cigar string in the alignment of fosmid reads 

```
smallINDEL=fosmid.smallINDEL.bed;
samtools view fosmid.sorted.bam | perl -ane 'next if($F[1] & 0x100 || $F[1] & 0x400 || $F[1] & 0x800); $chr=$F[2]; $pos = $F[3]; $cigar=$F[5]; while($cigar ne ""){($num, $tag, $cigar) = ($cigar =~ /^(\d+)(\S)(.*)/); if($tag eq "I"){print join("\t", $chr, $pos, ($pos + 1), "I", $num) . "\n"}elsif($tag eq "D"){print join("\t", $chr, $pos, ($pos + $num), "D", $num) . "\n"}if($tag =~ /[DM]/){$pos += $num;}}; ' > $smallINDEL;
```

small insertion

```
perl -ane 'print $_ if($F[3] eq "I")' $smallINDEL | perl -ane 'print join("\t", $F[0], $F[1], ($F[1] + 1), $F[4]) . "\n" if($F[4] <= 50)' > fosmid.smallINS.bed
```

small deletion

```
perl -ane 'print $_ if($F[3] eq "D")' $smallINDEL | perl -ane 'print join("\t", $F[0], $F[1], $F[2], $F[4]) . "\n" if($F[4] <= 50)' > fosmid.smallDEL.bed
```

large deletion

```
samtools view fosmid.sorted.bam | perl -ane 'next if($F[1] & 0x100 || $F[1] & 0x400 || $F[1] & 0x800); $pos = $F[3]; if($F[5] =~ /[ID]/){$cigar = $F[5]; while($cigar ne ""){($num, $tag, $cigar) = ($cigar =~ /^(\d+)(\S)(.*)/); if($num >= 50 && $tag =~ /[DI]/){print join("\t", $F[2], $F[1], $pos, $tag, $num) . "\n"}; if($tag =~ /[MD]/){$pos += $num;}}}' | grep -v alt | grep D | perl -ane 'print join("\t", $F[0], $F[2], ($F[2] + $F[4]), $F[4]) . "\n"' > fosmid.largeDEL.bed
```

### Step 3. Validation of small deletions (padding 50bp) (hysa, platinum and giab calls are named smallDEL.hysa.bed, smallDEL.platinum.bed and smallDEL.giab.bed, respectively. Each input bed file has no padding, and the fourth column (1-based) is the insertion size). 

```
gs_del=fosmid.smallDEL.bed;
for i in hysa platinum giab; 
do 
    echo $i; 
    file=smallDEL.$i.bed; 
    ~/bin/intersectBed.touch.sh $file $cover | perl -ane 'print join("\t", @F[0 .. 2], ($F[2] - $F[1]) ) . "\n"' > $file.o; 
    file=$file.o; 
    wc -l $file; 
    perl -ane 'print join("\t", $F[0], ($F[1] - 50), ($F[2] + 50), $F[3]) . "\n"' $file > $file.padding50; file=$file.padding50; 
    ~/pkg/bedtools/bedtools-2.17.0/bin/intersectBed -a $file -b $gs_del -wo | perl -ane '$len1=$F[3]; $len2 = $F[$#F - 1]; if(abs($len1 - $len2) <= 10){print join("\t", @F[0 .. 2]) . "\n"}' | sort | uniq | wc -l; 
done
```

### Step 4. Validation of small insertions (padding 50bp) (hysa, platinum and giab calls are named smallINS.hysa.bed, smallINS.platinum.bed and smallINS.giab.bed, respectively. Each input bed file has no padding, and the fourth column (1-based) is the insertion size). 

```
gs_ins=fosmid.smallINS.bed;
for i in hysa platinum giab
do 
    echo $i; 
    file=smallINS.$i.bed; perl -ane 'print join("\t", $F[0], ($F[1] - 50), ($F[2] + 50), @F[ 3 .. $#F ]) . "\n"' $file > $file.padding50; 
    file=$file.padding50; 
    ~/bin/intersectBed.touch.sh $file $cover > $file.o; 
    file=$file.o; 
    wc -l $file; 
    ~/pkg/bedtools/bedtools-2.17.0/bin/intersectBed -a $file -b $gs_ins -wo | perl -ane '$len1=$F[3]; $len2 = $F[$#F - 1]; if(abs($len1 - $len2) <= 10){print join("\t", @F[0 .. 2]) . "\n"}' | sort | uniq | wc -l; 
done
```

### Step 5. Validation of large deletions (no padding, 1bp overlap) (hysa, delly, pbhoney, customized approach and svclassify calls are named largeDEL.hysa.bed, largeDEL.delly.bed, largeDEL.pbhoney.bed, largeDEL.custom.bed, largeDEL.svclassify.bed, respectively. Each input bed file has no padding, and the fourth column (1-based) is the deletion size).

```
gs_l_del=fosmid.largeDEL.bed;
for i in ours delly pbhoney custom svclassify
do 
    echo $i; 
    file=largeDEL.$i.bed; ~/bin/intersectBed.touch.sh $file $cover > $file.o; 
    file=$file.o; 
    wc -l $file; 
    ~/pkg/bedtools/bedtools-2.17.0/bin/intersectBed -a $file -b $gs_l_del -wo | perl -ane 'print join("\t", @F[0 .. 2]) ."\n" if(abs($F[4] - $F[$#F - 1]) < 10)' | sort | uniq | wc -l; 
done
```

