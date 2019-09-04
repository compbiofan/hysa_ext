# Supporting scripts for analysis of NA12878 small deletion and insertion.

## Steps

### Step 1. Prepare Platinum and GIAB calls.

Download Platinum calls from ftp://platgene_ro@ussd-ftp.illumina.com/hg19/8.0.1/NA12878/NA12878.vcf.gz. Convert from vcf to bed and lift over to build38 with a filtering of size (only select size > 10bp):

```
gunzip NA12878.vcf.gz
mv NA12878.vcf platinum.vcf
./process_platinum.sh platinum.vcf
```

Results are in platinum.del.bed and platinum.ins.bed for deletion and insertion, respectively.

Download GIAB calls from ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.vcf.gz. Convert from vcf to bed and lift over to build38 with a filtering of size (only select size > 10bp):

```
gunzip NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.vcf.gz 
mv NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.vcf giab.vcf
./process_GIAB.sh giab.vcf 
```

Results are in giab.del.bed and giab.ins.bed for deletion and insertion, respectively.

### Step 2

Generate venn diagram among platinum.bed, giab.bed and HySA calls for deletion and insertion, respectively.

```calculate_mendelian_noGT.pl platinum.del.bed giab.del.bed out.del.SET50.rmdup.bed touch```
```calculate_mendelian_noGT.pl platinum.ins.bed giab.ins.bed out.ins.SET50.flanking50.rmdup.bed touch```

out.del.SET50.rmdup.bed and out.ins.SET50.flanking50.rmdup.bed are the small deletion and insertion files produced by HySA.

### Step 3. Overlap our unique calls with database (dbSNP and short tandem repeat).

Prepare database files in the same way as in CHM1/small_indel/README.md. Liftover both database files to build38. 

Overlap our unique small indels with the two databases:

```intersectBed.touch.sh out.del.SET50.rmdup.bed.only $STR 1```
```intersectBed.touch.sh out.del.SET50.rmdup.bed.only $dbSNP_del 1```
```intersectBed.touch.sh out.ins.SET50.flanking50.rmdup.bed.only $STR 1```
```intersectBed.touch.sh out.ins.SET50.flanking50.rmdup.bed.only $dbSNP_ins 1```

$STR is simpleRepeat.build38.bed. $dbSNP_del is snp146.deletion.build38.bed. $dbSNP_ins is snp146.insertion.build38.bed.

The four numbers shown from these four commands are the number of HySA unique small deletion calls overlapping with short tandem repeat database, dbSNP database, and the number of HySA unique small insertion calls overlapping with short tandem repeat database and dbSNP database. 

