#!/usr/bin/env bash
if [ $# -lt 2 ]; then
    echo "This script is to validate large insertions with size > 500bp by aligning the inserted sequence to build38, pacbio whole genome assembly (Pendleton M, Sebra R, Pang AW, Ummat A, Franzen O, Rausch T, Stutz AM, Stedman W, Anantharaman T, Hastie A et al. 2015. Assembly and diploid architecture of an individual human genome via single-molecule technologies. Nat Methods 12: 780-786.) and fosmid clones (Kidd JM, Cooper GM, Donahue WF, Hayden HS, Sampas N, Graves T, Hansen N, Teague B, Alkan C, Antonacci F et al. 2008. Mapping and sequencing of structural variation from eight human genomes. Nature 453: 56-64.)\nUsage: $0 <<HySA_output_in_bed_deduplicated> <HySA_assembly_in_fa> <build38.fasta> <wg_asm.fasta> <fosmid.fasta> <processor_num>\n";
    exit
fi

# input is in the original HySA output format, deduplicated
in=$1
# fasta of the assembled contigs
fa=$2
ln -s $fa ./all.ctg.fa
# build38 fasta
build38=$3
# whole genome diploid assembly 
wg_asm=$4
# fosmid clones
fosmid=$5
# number of processor in running this script
nproc=$6

# prepare files
# convert from bed to txt since deduplication was operated on bed format, with filtering applied.
perl -ane 'print join("\t", $F[0], ($F[1] + 50), @F[3 .. $#F]) . "\n" if($F[3] > 500)' $in > $in.selected.txt
perl -ane 'my @a = split(/\//, $F[4]); my $x = `samtools faidx all.ctg.fa $a[0]`; print $x' $in.selected.txt > $in.selected.fa
extract_inserted.pl $in.selected.fa $in.selected.txt 0 > $in.selected.inserted.fa
mkdir out_ctg_aln
ref=$build38
blasr -maxAnchorsPerPosition 100 -advanceExactMatches 10 -affineAlign -affineOpen 100 -affineExtend 0 -insertion 5 -deletion 5 -extend -maxExtendDropoff 20 -clipping subread -bestn 3 -nproc $nproc -sam -out out_ctg_aln/blasr_ctg2ref.build38.sam $in.selected.inserted.fa $ref
analyze_inserted_sequence.pl out_ctg_aln/blasr_ctg2ref.build38.sam 0.9 100 50 > build38.ins.txt
# see how many can be aligned well to diploid's assembly for both well and poor aligned build38
less build38.ins.txt | grep "well_aligned" | perl -ane 'print join("\n", @F[1 .. $#F]) . "\n"' > build38.well_aligned.ins.txt
less build38.ins.txt | grep "poor_aligned" | perl -ane 'print join("\n", @F[1 .. $#F]) . "\n"' > build38.poor_aligned.ins.txt
get_disappeared.pl out_ctg_aln/blasr_ctg2ref.build38.sam $in.selected.forBuild38.inserted.fa >> build38.poor_aligned.ins.txt
extract_inserted.pl $in.selected.fa $in.selected.txt 0 build38.poor_aligned.ins.txt > $in.selected.poor_aligned_of_build38.inserted.fa
ref=$wg_asm
blasr -maxAnchorsPerPosition 100 -advanceExactMatches 10 -affineAlign -affineOpen 100 -affineExtend 0 -insertion 5 -deletion 5 -extend -maxExtendDropoff 20 -clipping subread -bestn 3 -nproc $nproc -sam -out out_ctg_aln/blasr_ctg2diploid.sam $in.selected.poor_aligned_of_build38.inserted.fa $ref
analyze_inserted_sequence.pl out_ctg_aln/blasr_ctg2diploid.sam 0.9 100 50 > diploid.ins.txt
less diploid.ins.txt | grep "well_aligned" | perl -ane 'print join("\n", @F[1 .. $#F]) . "\n"' > poor2build38.well_aligned_diploid.ins.txt
less diploid.ins.txt | grep "poor_aligned" | perl -ane 'print join("\n", @F[1 .. $#F]) . "\n"' > poor2build38.poor_aligned_diploid.ins.txt
get_disappeared.pl out_ctg_aln/blasr_ctg2diploid.sam $in.selected.poor_aligned_of_build38.inserted.fa >> poor2build38.poor_aligned_diploid.ins.txt
echo "In all, there are the following number of insertion:"
perl -ane '@a = split(/\//, $F[4]); print "$a[0]:$F[5]-$F[6]\n"' $in.selected.txt | sort | uniq | wc -l
echo "in which the following number of insertions can be well aligned to build38, thus not novel:"
wc -l build38.well_aligned.ins.txt
echo "and the following number of insertions cannot be well aligned to build38, thus novel:"
wc -l build38.poor_aligned.ins.txt
echo "For those that can be well aligned to build38, the following number can be well aligned to diploid's detection:"
wc -l poor2build38.well_aligned_diploid.ins.txt
echo "and the following number cannot be aligned to diploid's detection:"
wc -l poor2build38.poor_aligned_diploid.ins.txt
extract_inserted.pl $in.selected.fa $in.selected.txt 0 poor2build38.poor_aligned_diploid.ins.txt > $in.selected.poor_aligned_of_diploid.inserted.fa
ref=$fosmid
blasr -maxAnchorsPerPosition 100 -advanceExactMatches 10 -affineAlign -affineOpen 100 -affineExtend 0 -insertion 5 -deletion 5 -extend -maxExtendDropoff 20 -clipping subread -bestn 3 -nproc $nproc -sam -out out_ctg_aln/blasr_ctg2fosmid.sam $in.selected.poor_aligned_of_diploid.inserted.fa $ref
analyze_inserted_sequence.pl out_ctg_aln/blasr_ctg2fosmid.sam 0.9 100 50 > fosmid.ins.txt
less fosmid.ins.txt | grep "well_aligned" | perl -ane 'print join("\n", @F[1 .. $#F]) . "\n"' > poor2diploid.well_aligned_fosmid.ins.txt
less fosmid.ins.txt | grep "poor_aligned" | perl -ane 'print join("\n", @F[1 .. $#F]) . "\n"' > poor2diploid.poor_aligned_fosmid.ins.txt
get_disappeared.pl out_ctg_aln/blasr_ctg2fosmid.sam $in.selected.poor_aligned_of_diploid.inserted.fa >> poor2diploid.poor_aligned_fosmid.ins.txt
echo "in which the following number of insertions can be well aligned to fosmid, thus not novel:"
wc -l poor2diploid.well_aligned_fosmid.ins.txt
echo "and the following number of insertions cannot be well aligned to fosmid, thus novel:"
wc -l poor2diploid.poor_aligned_fosmid.ins.txt
