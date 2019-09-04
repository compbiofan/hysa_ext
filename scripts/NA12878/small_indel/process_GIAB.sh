#!/usr/bin/env bash
if [ $# -lt 2 ]; then
    echo "This script takes the gunzipped GIAB vcf file and convert the deletion and insertion files to bed, followed by liftover to build38. It generates two files: GIAB.del.bed and GIAB.ins.bed for deletion and insertion. Note: a filtering of size > 10bp is applied.\nUsage: $0 <gunzipped_vcf>\n";
    exit
fi

# input is in the original HySA output format, deduplicated
in=$1
perl -ane '@a = split(/,/, $F[3]); @b = split(/,/, $F[4]); foreach my $aa (@a){if(length($aa) > 10){$tag = 1}} foreach my $bb (@b){if(length($bb) > 10){$tag = 1;}} if($tag == 1){print $_}' $in > $in.GT10.vcf

# get the actual deletion
perl -ane '@a = split(/\,/, $F[3]); @b = split(/\,/, $F[4]); next if(@a > 1 && @b > 1); if(@a == 1){foreach my $c (@b){ if(abs(length($c) - length($F[3])) > 10){if(length($c) < length($F[3])){print join("\t", $F[0], $F[1], ($F[1] + length($F[3]) - length($c))) . "\n"}}}}' $in.GT10.vcf  > $in.GT10.del.bed
perl -ane '@a = split(/\,/, $F[3]); @b = split(/\,/, $F[4]); next if(@a > 1 && @b > 1); if(@b == 1 && @a != 1){foreach my $c (@a){ if(abs(length($c) - length($F[4])) > 10){if(length($c) > length($F[4])){print join("\t", $F[0], $F[1], ($F[1] + length($c) - length($F[3]))) . "\n"}}}}' $in.GT10.vcf  >> $in.GT10.del.bed
liftover_hg19_hg38.sh $in.GT10.del.bed
ln -s $in.GT10.del.bed.wchr.liftOver.build38 ./GIAB.del.bed 

# get the actual insertion
perl -ane '@a = split(/\,/, $F[3]); @b = split(/\,/, $F[4]); next if(@a > 1 && @b > 1); if(@a == 1){foreach my $c (@b){ if(abs(length($c) - length($F[3])) > 10){if(length($c) > length($F[3])){print join("\t", $F[0], $F[1], abs(length($F[3]) - length($c))) . "\n"}}}}' $in.GT10.vcf  > $in.GT10.ins.bed
perl -ane '@a = split(/\,/, $F[3]); @b = split(/\,/, $F[4]); next if(@a > 1 && @b > 1); if(@b == 1 && @a != 1){foreach my $c (@a){ if(abs(length($c) - length($F[4])) > 10){if(length($c) < length($F[4])){print join("\t", $F[0], $F[1], abs(length($c) - length($F[3]))) . "\n"}}}}' $in.GT10.vcf  >> $in.GT10.ins.bed
perl -ane 'print join("\t", $F[0], ($F[1] - 50), ($F[1] + 50)) . "\n"' $in.GT10.ins.bed > $in.GT10.ins.flanking50.bed
liftover_hg19_hg38.sh $in.GT10.ins.flanking50.bed
ln -s $in.GT10.ins.flanking50.bed.wchr.liftOver.build38 ./GIAB.ins.bed 

