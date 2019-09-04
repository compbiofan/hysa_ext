#!/usr/bin/env bash
if [ $# -lt 2 ]; then
    echo "This script is to validate large insertions with size > 500bp by aligning the inserted sequence to build37, build38, assembly 1 (Chaisson MJP, Huddleston J, Dennis MY, Sudmant PH, Malig M, Hormozdiari F, Antonacci F, Surti U, Sandstrom R, Boitano M et al. 2015b. Resolving the complexity of the human genome using single-molecule sequencing. Nature 517: 608-U163) and assembly 2 (Berlin K, Koren S, Chin CS, Drake JP, Landolin JM, Phillippy AM. 2015. Assembling large genomes with single-molecule sequencing and locality-sensitive hashing. Nat Biotechnol 33: 623-630). Note: since assembly 1 has corresponding insertion coordinates, besides comparing sequences, we also compare coordinates. Samtools and blasr will be run in this script.\nUsage: $0 <HySA_output_in_bed_deduplicated> <HySA_assembly_in_fa> <build37.fasta> <build38.fasta> <assembly1.bed> <assembly2.fasta> <processor_num>\n";
    exit
fi

# input is in the original HySA output format, deduplicated
in=$1
# fasta of the assembled contigs
fa=$2
ln -s $fa ./all.ctg.fa
# build37 fasta
build37=$3
# build38 fasta
build38=$4
# Chaisson's insertion in bed
cha_bed=$5
# phillippy's assembly
phi_asm=$6
# number of processor in running this script
nproc=$7

# prepare files
# convert from bed to txt since deduplication was operated on bed format
perl -ane 'print join("\t", $F[0], ($F[1] + 50), @F[3 .. $#F]) . "\n"' $in > $in.out
in_ori=$in
in=$in.out
perl -ane 'if($F[3] eq "I" && $F[2] >= 500){print $_}' $in > $in.selected.txt
perl -ane 'my @a = split(/\//, $F[4]); my $x = `samtools faidx all.ctg.fa $a[0]`; print $x' $in.selected.txt > $in.selected.fa
extract_inserted.pl $in.selected.fa $in.selected.txt 0 > $in.selected.inserted.fa
mkdir out_ctg_aln
# make chaisson insertion fasta file
perl -ane '$F[0] =~ s/^chr//; if(length($F[5])>500){print join("\n", ">$F[0]:$F[1]-$F[2]", $F[5]) . "\n"}' $cha_bed > $cha_bed.fa
perl -ane '$F[0] =~ s/^chr//; print join("\t", @F[0 .. 4]) . "\n"' $cha_bed > $cha_bed.coord
cha_asm=$cha_bed.fa
cha_coord_bed=$cha_bed.coord

# align to build37 
ref=$build37
blasr -maxAnchorsPerPosition 100 -advanceExactMatches 10 -affineAlign -affineOpen 100 -affineExtend 0 -insertion 5 -deletion 5 -extend -maxExtendDropoff 20 -clipping subread -bestn 3 -nproc $nproc -sam -out out_ctg_aln/blasr_ctg2ref.build37.sam $in.selected.inserted.fa $ref
# check how many inserted sequence can be found on the reference
analyze_inserted_sequence.pl out_ctg_aln/blasr_ctg2ref.build37.sam 0.9 100 50 > build37.ins.txt
less build37.ins.txt | grep "well_aligned" | perl -ane 'print join("\n", @F[1 .. $#F]) . "\n"' > build37.well_aligned.ins.txt
less build37.ins.txt | grep "poor_aligned" | perl -ane 'print join("\n", @F[1 .. $#F]) . "\n"' > build37.poor_aligned.ins.txt
# add those that do not appear in the sam file to the poor_aligned pool
get_disappeared.pl out_ctg_aln/blasr_ctg2ref.build37.sam $in.selected.inserted.fa >> build37.poor_aligned.ins.txt
get_txt1.pl build37.poor_aligned.ins.txt $in.selected.txt > build37.poor_aligned.ins.txt.txt
extract_inserted.pl $in.selected.fa build37.poor_aligned.ins.txt.txt 0 > $in.selected.forBuild38.inserted.fa 

# align to build38
ref=$build38
blasr -maxAnchorsPerPosition 100 -advanceExactMatches 10 -affineAlign -affineOpen 100 -affineExtend 0 -insertion 5 -deletion 5 -extend -maxExtendDropoff 20 -clipping subread -bestn 3 -nproc $nproc -sam -out out_ctg_aln/blasr_ctg2ref.build38.sam $in.selected.forBuild38.inserted.fa $ref
analyze_inserted_sequence.pl out_ctg_aln/blasr_ctg2ref.build38.sam 0.9 100 50 > build38.ins.txt
# see how many can be aligned well to Chaisson's assembly for both well and poor aligned build38
less build38.ins.txt | grep "well_aligned" | perl -ane 'print join("\n", @F[1 .. $#F]) . "\n"' > build38.well_aligned.ins.txt
less build38.ins.txt | grep "poor_aligned" | perl -ane 'print join("\n", @F[1 .. $#F]) . "\n"' > build38.poor_aligned.ins.txt
get_disappeared.pl out_ctg_aln/blasr_ctg2ref.build38.sam $in.selected.forBuild38.inserted.fa >> build38.poor_aligned.ins.txt
get_txt1.pl build38.poor_aligned.ins.txt $in.selected.txt > build38.poor_aligned.ins.txt.txt
get_txt1.pl build38.well_aligned.ins.txt $in.selected.txt > build38.well_aligned.ins.txt.txt
extract_inserted.pl $in.selected.fa build38.well_aligned.ins.txt.txt 0 > $in.selected.well_aligned_of_build38.inserted.fa
extract_inserted.pl $in.selected.fa build38.poor_aligned.ins.txt.txt 0 > $in.selected.poor_aligned_of_build38.inserted.fa
r

# align to Chaisson
ref=$cha_asm
blasr -maxAnchorsPerPosition 100 -advanceExactMatches 10 -affineAlign -affineOpen 100 -affineExtend 0 -insertion 5 -deletion 5 -extend -maxExtendDropoff 20 -clipping subread -bestn 3 -nproc $nproc -sam -out out_ctg_aln/blasr_ctg2Chaisson.well.sam $in.selected.well_aligned_of_build38.inserted.fa $ref
analyze_inserted_sequence.pl out_ctg_aln/blasr_ctg2Chaisson.well.sam 0.9 100 50 > well2Chaisson.ins.txt
less well2Chaisson.ins.txt | grep "well_aligned" | perl -ane 'print join("\n", @F[1 .. $#F]) . "\n"' > well2Chaisson.well_aligned.ins.txt
less well2Chaisson.ins.txt | grep "poor_aligned" | perl -ane 'print join("\n", @F[1 .. $#F]) . "\n"' > well2Chaisson.poor_aligned.ins.txt
get_disappeared.pl out_ctg_aln/blasr_ctg2Chaisson.well.sam $in.selected.well_aligned_of_build38.inserted.fa >> well2Chaisson.poor_aligned.ins.txt
# now add coordinates confirmation as well, those in well2Chaisson.poor_aligned.yetCoordOverlap.txt should be added to well2Chaisson.well_aligned.ins.txt and substracted from well2Chaisson.poor_aligned.ins.txt
get_coordinate1.pl well2Chaisson.poor_aligned.ins.txt $in_ori | sort > well2build38_poor2Chaisson.bed
intersectBed.touch.sh well2build38_poor2Chaisson.bed $cha_coord_bed | perl -ane 'print join("\t", @F[0 .. 3], "build38", "Y") . "\n"' > well2Chaisson.poor_aligned.yetCoordOverlap.txt

# align to Phillippy's assembly
# now deal with alignment of the poor2build38 to whole genome assembly of Phillippy's 
ref=$phi_asm
blasr -maxAnchorsPerPosition 100 -advanceExactMatches 10 -affineAlign -affineOpen 100 -affineExtend 0 -insertion 5 -deletion 5 -extend -maxExtendDropoff 20 -clipping subread -bestn 3 -nproc $nproc -sam -out out_ctg_aln/blasr_ctg2Philippy.poor.sam $in.selected.poor_aligned_of_build38.inserted.fa $ref
analyze_inserted_sequence.pl out_ctg_aln/blasr_ctg2Philippy.poor.sam 0.9 100 50 > poor2Philippy.ins.txt
less poor2Philippy.ins.txt | grep "well_aligned" | perl -ane 'print join("\n", @F[1 .. $#F]) . "\n"' > poor2Philippy.well_aligned.ins.txt
less poor2Philippy.ins.txt | grep "poor_aligned" | perl -ane 'print join("\n", @F[1 .. $#F]) . "\n"' > poor2Philippy.poor_aligned.ins.txt
get_disappeared.pl out_ctg_aln/blasr_ctg2Philippy.poor.sam $in.selected.poor_aligned_of_build38.inserted.fa >> poor2Philippy.poor_aligned.ins.txt
get_txt1.pl poor2Philippy.poor_aligned.ins.txt $in.selected.txt > poor2Philippy.poor_aligned.ins.txt.txt
get_txt1.pl poor2Philippy.well_aligned.ins.txt $in.selected.txt > poor2Philippy.well_aligned.ins.txt.txt
extract_inserted.pl $in.selected.fa poor2Philippy.well_aligned.ins.txt.txt 0 > $in.selected.poor_aligned_of_build38.well_aligned_of_Philippy.inserted.fa
extract_inserted.pl $in.selected.fa poor2Philippy.poor_aligned.ins.txt.txt 0 > $in.selected.poor_aligned_of_build38.poor_aligned_of_Philippy.inserted.fa

# check how many well aligned to whole genome but not build38 can be well aligned to Chaisson
ref=$cha_asm
blasr -maxAnchorsPerPosition 100 -advanceExactMatches 10 -affineAlign -affineOpen 100 -affineExtend 0 -insertion 5 -deletion 5 -extend -maxExtendDropoff 20 -clipping subread -bestn 3 -nproc $nproc -sam -out out_ctg_aln/blasr_ctg2Chaisson.well2Philippy.sam $in.selected.poor_aligned_of_build38.well_aligned_of_Philippy.inserted.fa $ref
analyze_inserted_sequence.pl out_ctg_aln/blasr_ctg2Chaisson.well2Philippy.sam 0.9 100 50 > poor2build38.well2Philippy.align2Chaisson.ins.txt
less poor2build38.well2Philippy.align2Chaisson.ins.txt | grep "well_aligned" | perl -ane 'print join("\n", @F[1 .. $#F]) . "\n"' > poor2build38.well2Philippy.align2Chaisson.well_aligned.ins.txt
less poor2build38.well2Philippy.align2Chaisson.ins.txt | grep "poor_aligned" | perl -ane 'print join("\n", @F[1 .. $#F]) . "\n"' > poor2build38.well2Philippy.align2Chaisson.poor_aligned.ins.txt
get_disappeared.pl out_ctg_aln/blasr_ctg2Chaisson.well2Philippy.sam $in.selected.well_aligned_of_build38.well_aligned_of_Philippy.inserted.fa >> poor2build38.well2Philippy.align2Chaisson.poor_aligned.ins.txt
get_coordinate1.pl poor2build38.well2Philippy.align2Chaisson.poor_aligned.ins.txt $in_ori | sort > poor2build38.well2Philippy.poor2Chaisson.bed
intersectBed.touch.sh poor2build38.well2Philippy.poor2Chaisson.bed $cha_coord_bed | perl -ane 'print join("\t", @F[0 .. 3], "build38", "Y") . "\n"' > poor2build38.well2Philippy.poor2Chaisson.yetCoordOverlap.txt

# check how many poor aligned to whole genome and not build38 can be well aligned to Chaisson
blasr -maxAnchorsPerPosition 100 -advanceExactMatches 10 -affineAlign -affineOpen 100 -affineExtend 0 -insertion 5 -deletion 5 -extend -maxExtendDropoff 20 -clipping subread -bestn 3 -nproc $nproc -sam -out out_ctg_aln/blasr_ctg2Chaisson.poor2Philippy.sam $in.selected.poor_aligned_of_build38.poor_aligned_of_Philippy.inserted.fa $ref
analyze_inserted_sequence.pl out_ctg_aln/blasr_ctg2Chaisson.poor2Philippy.sam 0.9 100 50 > poor2build38.poor2Philippy.align2Chaisson.ins.txt
less poor2build38.poor2Philippy.align2Chaisson.ins.txt | grep "well_aligned" | perl -ane 'print join("\n", @F[1 .. $#F]) . "\n"' > poor2build38.poor2Philippy.align2Chaisson.well_aligned.ins.txt
less poor2build38.poor2Philippy.align2Chaisson.ins.txt | grep "poor_aligned" | perl -ane 'print join("\n", @F[1 .. $#F]) . "\n"' > poor2build38.poor2Philippy.align2Chaisson.poor_aligned.ins.txt
get_disappeared.pl out_ctg_aln/blasr_ctg2Chaisson.poor2Philippy.sam $in.selected.well_aligned_of_build38.poor_aligned_of_Philippy.inserted.fa >> poor2build38.poor2Philippy.align2Chaisson.poor_aligned.ins.txt
get_coordinate1.pl poor2build38.poor2Philippy.align2Chaisson.poor_aligned.ins.txt $in_ori | sort > poor2build38.poor2Philippy.poor2Chaisson.bed
intersectBed.touch.sh poor2build38.poor2Philippy.poor2Chaisson.bed $cha_coord_bed | perl -ane 'print join("\t", @F[0 .. 3], "build38", "Y") . "\n"' > poor2build38.poor2Philippy.poor2Chaisson.yetCoordOverlap.txt

echo "In all, there are the following number of insertion:"
perl -ane '@a = split(/\//, $F[4]); print "$a[0]:$F[5]-$F[6]\n"' out.ins.selected.txt | sort | uniq | wc -l 
echo "in which the following number of insertions can be well aligned to build37, thus not novel:"
wc -l build37.well_aligned.ins.txt
echo "and the following number of insertions cannot be well aligned to build37, thus novel:"
wc -l build37.poor_aligned.ins.txt
echo "Of those that cannot be aligned to build37, the following number of insertions can be well aligned to build38:"
wc -l build38.well_aligned.ins.txt
echo "and the following number cannot be aligned to build38:"
wc -l build38.poor_aligned.ins.txt
echo "For those that can be well aligned to build38, the following number can be well aligned to Chaisson's detection:"
wc -l well2Chaisson.well_aligned.ins.txt
echo "and the following number cannot be aligned to Chaisson's detection:"
wc -l well2Chaisson.poor_aligned.ins.txt
echo "yet the following number have coordinates overlapping with Chaisson's detection:"
wc -l well2Chaisson.poor_aligned.yetCoordOverlap.txt
echo "For those that cannot be well aligned to build38, the following number can be well aligned to whole genome assembly:"
wc -l poor2Philippy.well_aligned.ins.txt 
echo "and the following number cannot be aligned to whole genome assembly:"
wc -l poor2Philippy.poor_aligned.ins.txt
echo "For those that cannot be well aligned to build38, can be well aligned to whole genome assembly, the following number can be well aligned to Chaisson's detection:"
wc -l poor2build38.well2Philippy.align2Chaisson.well_aligned.ins.txt
echo "For those that cannot be well aligned to build38, can be well aligned to whole genome assembly, the following number cannot be well aligned to Chaisson's detection:"
wc -l poor2build38.well2Philippy.align2Chaisson.poor_aligned.ins.txt
echo "yet the following number have coordinates overlapping with Chaisson's detection:"
wc -l poor2build38.well2Philippy.poor2Chaisson.yetCoordOverlap.txt
echo "For those that cannot be well aligned to build38, cannot be well aligned to whole genome assembly, the following number can be well aligned to Chaisson's detection:"
wc -l poor2build38.poor2Philippy.align2Chaisson.well_aligned.ins.txt
echo "For those that cannot be well aligned to build38, cannot be well aligned to whole genome assembly, the following number cannot be well aligned to Chaisson's detection:"
wc -l poor2build38.poor2Philippy.align2Chaisson.poor_aligned.ins.txt
echo "yet the following number have coordinates overlapping with Chaisson's detection:"
wc -l poor2build38.poor2Philippy.poor2Chaisson.yetCoordOverlap.txt

