#!/usr/bin/env bash
if [ $# -lt 0 ]; then
    echo "This script partition HySA calls by size (50bp) and types (insertion, deletion), followed by deduplication. Note: insertion will be padded 50bp on the left and right. \nUsage: $0\n"
  exit
fi


# large deletion
perl -ane 'chomp;print join("\t", $F[0], $F[1], ($F[1] + $F[2]), @F[2 .. $#F]) . "\n" if($F[2] > 50 && $F[3] eq "D")' out.all.txt > out.del.GT50.bed
rm_dup_bed1.pl out.del.GT50.bed reciprocal50 > out.del.GT50.rmdup.bed

# large insertion
perl -ane 'chomp;print join("\t", $F[0], ($F[1] - 50), ($F[1] + 50), @F[2 .. $#F]) . "\n" if($F[2] > 50 && $F[3] eq "I")' out.all.txt > out.ins.GT50.flanking50.bed
rm_dup_bed1.pl out.ins.GT50.flanking50.bed touch > out.ins.GT50.flanking50.rmdup.bed

# small deletion
perl -ane 'chomp;print join("\t", $F[0], $F[1], ($F[1] + $F[2]), @F[2 .. $#F]) . "\n" if($F[2] <= 50 && $F[3] eq "D")' out.all.txt > out.del.SET50.bed
rm_dup_bed1.pl out.del.SET50.bed touch > out.del.SET50.rmdup.bed

# small insertion
perl -ane 'chomp;print join("\t", $F[0], ($F[1] - 50), ($F[1] + 50), @F[2 .. $#F]) . "\n" if($F[2] <= 50 && $F[3] eq "I")' out.all.txt > out.ins.SET50.flanking50.bed
rm_dup_bed1.pl out.ins.SET50.flanking50.bed touch > out.ins.SET50.flanking50.rmdup.bed

`rm tmp.intersect.bed`;
