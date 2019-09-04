# Supporting scripts for analysis of CHM1 large deletion.

## Steps

### Step 1. Curate gold standard call set from Delly and Pacbio

```./curate.sh <bam_IL> <delly_call> <PB_call> <reference>```

$bam_IL is the Illumina bam aligned against the reference (build37 or build38). $delly_call and $PB_call are bed files from delly and Chaisson large (>50bp) deletion call set. $reference is either build37 or build38 (in CHM1, we use build37).

### Step 2. Validate HySA callset and compare with curated ones

```./run_validate_comp.sh <ref> <hysa_del> <WGA> <curated> <bam>```

$ref is the reference against which the Illumina bam was obtained. $hysa_del is hysa deletion bed file to be validated and compared (here, it is out.del.GT50.rmdup.bed as described in pipeline). $WGA is the whole genome assembly for validating HySA calls (here, it is all.ctg.fa as described in pipeline). $curated is the curated call set in bed format. $bam is the Illumina bam file aligned to build37 or build38. 

### Step 3. Find hysa validated calls that is not in curate call 

```perl -ane 'print join("\t", @F[0 .. 2], $F[$#F - 1]) . "\n" if($F[$#F] == 0 && $F[$#F - 2] == 1)' $hysa_del.wValidByCHMASM.wSRGT.overlapWGT > unique.bed```

$hysa_del is the same $hysa_del in the previous step.

### Step 4. Calculate the supporting IL and PB # and mean and standard deviation of the distance from the PB gap start position to the start position of the call. 

```./calc_mean_std.sh <bam_IL> <bam_PB> <bed> <ref>```

$bam_IL and $bam_PB are the bams of Illumina and Pacbio reads alignment to the reference: build37 or build38. $bed is the call set upon which these metrics are to be calculated. In our case, they are the output of curate.sh: both.ILH.get2.rmdup (the curated calls from Delly and Pacbio) and unique.bed (HySA validated unique calls). 

### Step 5. Prepare files for comparison in R.

```./get_p_value.sh both.ILH.get2.rmdup.wSRGT.wPBGT.summaryAll unique.bed.wSRGT.wPBGT.summaryAll```

### Step 6. Run R to get the p values. 

Launch R. Then type
```source comp.R```

The results are subsequently the p value for # of supporting Illumina reads, mean of starting position of PB gap to the call, standard deivation of starting position of PB gap to the call.  

### Step 7. Calculate FDR.

Calculate not validated by assembly (col 9 from hysa output format): (output as 1)
```perl -ane 'print $_ if($F[$#F-2] == 0)' $hysa_del.wValidByCHMASM.wSRGT.overlapWGT  | wc -l```

Calculate not validated by assembly but overlap with curated set: (output as 2)
```perl -ane 'print $_ if($F[$#F-2] == 0 && $F[$#F] == 1)' $hysa_del.wValidByCHMASM.wSRGT.overlapWGT  | wc -l```

Calculate not validated by assembly, not overlap with curated set, but can be genotype by >= 2 ILs: (output as 3)
```perl -ane 'print $_ if($F[$#F-2] == 0 && $F[$#F] == 0 && $F[$#F-1] >= 2)' $hysa_del.wValidByCHMASM.wSRGT.overlapWGT  | wc -l```

Those false discoveries are the number of 1 minus the sum of 2 and 3. 

## Note
### The steps cannot be run in parallel due to the dependency of the latter steps on the previous ones, and some temporary file with the same name appearing in more than one step. 
