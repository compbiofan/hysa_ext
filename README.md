### This repository is for reproducing results of the paper of 
### HySA: a Hybrid Structural variant Assembly approach using next-generation and single-molecule sequencing technologies Genome Research (19 January 2017), gr.214767.116, doi:10.1101/gr.214767.116 by Xian Fan, Mark Chaisson, Luay Nakhleh, Ken Chen

# Pipeline #

## Software Requirements ##

1. samtools (https://samtools.github.io)
2. perl 5.10.1 or up
3. customized blasr (downloadable from https://github.com/mchaisso/blasr with commit number 82553d6)
4. Infomap (downloadable from http://www.mapequation.org/code.html)
5. bwa 0.7.10 (downloadable from https://github.com/lh3/bwa)
6. g++ 4.4.7 or up
7. pairHMM (included in this package, need to install shown in the next section)
8. bedtools 2.17.0 or up (the binary intersectBed will be used)

## Environment Setup ##

Suppose "$DIR" is the path of this package.

* Add perl_modules to your PERL5LIB. In bash, 

  ```export PERL5LIB=$DIR/HybridAssemblySV/perl_modules```
  

* Add pipeline to your directory. In bash,

```HySA="$DIR/hybridassemblysv/pipeline"```

```if [ -d $HySA ]; then PATH="$HySA:$PATH" fi```

Similarly, make sure pls2fasta, blasr, Infomap and intersectBed are executable by adding the path to your directory. pls2fasta can be found in $blasr_dir/pbihdfutils/bin/. blasr can be found in $blasr_dir/alignment/bin/. Infomap can be found in $infomap_dir/. intersectBed could be found in bedtools package bin directory.

* Make sure all perl scripts and bash scripts in pipeline are executable.

```chmod u+x $DIR/hybridassemblysv/pipeline/*.pl```
```chmod u+x $DIR/hybridassemblysv/pipeline/*.sh```

* Install pairHMM binary. In $DIR/HybridAssemblySV/perl_modules/pairHMM, 

```g++ -g -Wall pairHMM.cpp -o pairHMM -lm -lz```

Add the binary's path to your directory.

```if [ -d $pairHMM_dir ]; then PATH="$pairHMM_dir:$PATH" fi ```

## Note for complex deletions

For inferring complex deletions, all steps are the same as the simple indels except in Step 6 and 7 (complex mode).

### Step 1. Extract Illumina reads given a bam file. ###

1.1 Extract unmapped reads

  ```extract_both_unmapped.sh <bam> <output_dir> <prefix>```

In $output_dir, files "$prefix.filtered.sam" and "index.$prefix.txt" will be generated.

1.2 Extract discordant read pairs, split reads and reads containing large gaps

1.2.1 Get all putative reads. 

```extract_short.pl -r <reference> -c <chromosome> -t <ref_region_for_insert> <bam> <output_dir> 1```

1.2.2 Extract intra-chromosomal reads.

```extract_short.pl -r <reference> -c <chromosome> -t <ref_region_for_insert> <bam> <output_dir> 2```

"$reference" is the fasta file of reference genome used in generating the bam file, e.g., NCBI build 38. "$chromosome" needs to be in the same format of the reference genome. "$ref_region_for_insert" is a region (in the format of chr:start-end) given by the user, the mean and standard deviation of insert size are computed based on reads mapped to this region. In the "$output_dir", files "discordant.chr$chromosome.sam" and "index.chr$chromosome.txt" will be generated.

1.2.3 Extract inter-chromosomal reads.

First, concatenate all candidate inter-chromosomal reads into one file:

```cat $output_dir/tmp.inter.*.txt > $output_dir/tmp.inter.txt```

Then,

```extract_short.pl -r <reference> -t <ref_region_for_insert> <bam> <output_dir> 3```

$reference, $ref_region_for_insert, $bam and $output_dir are all the same with that in 1.2.2.

The output is in $output_dir/discordant.inter.sam and $output_dir/tmp.inter.txt.

1.3 Combine all extracted files

```combine_all.pl <config_file> <new_output_dir> <new_prefix>```

Each row of the "$config_file" contains two columns, the first the extracted reads' sam file (from unmapped in 1.1, or intra-chromosomal ones in 1.2.2 or inter-chromosomal ones in 1.2.3), the second the corresponding index file, all in full path.

In the "$new_output_dir", four files "${new_prefix}_1.fa", "${new_prefix}_2.fa", "$new_prefix.id" and "$new_prefix.drp" will be generated.

### 2. Prepare long reads in proper format.###

2.1 If given a list of HDF5 Pacbio files, run the following command to convert it to fasta.

```pls2fasta <fofn> <fasta> -trimByRegion```

"$fofn" is a file with each row the full path of the HDF5 file. "$fasta" is the output file.

If the input Pacbio read is already in fasta format, concatenate them into one big fasta.

2.2 Partition fasta files into small ones.

```partition_pb.pl <fasta> <output_prefix> <mode>```

"$fasta" is the output of 2.1. "$mode" is either 37 or 38, depending on the reference you are using (NCBI build 37 or NCBI build 38). Give either one of them if the reference is different from these two. 

The output is a set of new folders, each folder with a name "$output_prefix$ID", in which "$ID" is a non-negative integer counting from 0. In each folder, two files "ref.fa" and "id.txt" are generated .

###3. Align extracted short reads to all long reads.###

3.1 Use the customized version of blasr to align "${new_prefix}_1.fa" and "${new_prefix}_2.fa" from 1.3, respectively, to each of "$output_prefix$ID/ref.fa" from 2.2, with the output files "${m4out}_1.m4" and "${m4out}_2.m4" in "$m4_dir".

```blasr -nCandidates 40 -minInterval 40 –maxScore 0 –minMatch 4 –maxMatch 13 -bestn 5 ${new_prefix}_1.fa -m 4 -out ${m4out}_1.m4 -nproc $proc $output_prefix$ID/ref.fa```
```blasr -nCandidates 40 -minInterval 40 –maxScore 0 –minMatch 4 –maxMatch 13 -bestn 5 ${new_prefix}_2.fa -m 4 -out ${m4out}_2.m4 -nproc $proc $output_prefix$ID/ref.fa```

in which "$proc" is the processor you will use for running the program. This process should be done on all "$ID"s generated in step 2.2.

3.2 Use paired end information to filter alignments.

```pairFilter.pl <m4_dir> <percentage_identity_threshold> <align_length_threshold> <minimum_insert> <m4out>```

"${m4out}" is the directory of each pair of m4 files, indicated in 3.1. Suggested "${percentage_identity_threshold}" is between 70 and 90. "${align_length_threshold}" depends on Illumina read length. If Illumina read length is 100, suggested number is between 70 and 90. "$minimum_insert" is the minimum total length of a Illumina read pair's alignment on a Pacbio read. "${m4out}" is the one indicated in 3.1. 

3.3 Use "pseudo" reference read to filter out long reads from reference allele.

First sort the output of 3.2:

```sort -k 1 -h <m4out> > <sorted_m4out>```

Next run the filter of reference allele:

```refFilter.pl <sorted_m4out> <out> <drp_file>```

"$sorted_m4out" is the output of the sort; "$out" is a file name, the output of this step; "$drp_file" is the same file as "$new_prefix.drp" of 1.3.

3.4 Use Illumina reads' coverage to filter out false alignments.

```covFilter.pl <refFiltered_m4out> <cov_threshold> <out>```

"$refFiltered_m4out" is the output of step 3.3 "$out".

"$cov_threshold" is a string concatenated by two integers with separator ":". The integers are the minimum and maximum allowed Illumina coverage, e.g., 3:50. Suggested maximum coverage is between 1 and 1.2 times of the Illumina sequence coverage.

3.5 Concatenate all the output of step 3.4 (in m4 format) together into one m4.

###4. Clustering of reads from m4 file.###

4.1 Union find.

```union_find.pl <m4> <union_find_prefix>```

"$m4" is the output of 3.5. "$union_find_prefix" is the prefix of the two output files ("$union_find_prefix.all" and $union_find_prefix.txt") in this step.

4.2 Infomap.

```infomap.sh <m4> <unionfind_prefix> <pb_threshold> <infomap_prefix>```

"$m4" is the output of 3.5. "$unionfind_prefix" is the same as 4.1. "$pb_threshold" is the minimum number of long reads in a cluster from Union-Find so that the cluster can be considered for community partition with Infomap, e.g., 5000. "$infomap_prefix" is the prefix of the output file in this step.

###5. Assembly of long reads###

```perl make_PB_assembly.pl <cluster_prefix> <ctg_output_dir> <pb_prefix> <tmp_dir> <time_limit> <cov_threshold> <cmd_output_file> <thread>```

The assembly is performed for both clusters of 4.1 and 4.2. When assembling clusters for 4.1, "$cluster_prefix" is "$union_find_prefix" from 4.1. When assembling clusters for 4.2, "$cluster_prefix" is "$infomap_prefix" from 4.2. "$ctg_output_dir" is the directory for the assembled contig, which needs to be separate for Union-Find and Infomap. "$pb_prefix" is "$output_prefix" of step 2.2. "$tmp_dir" is the temporary directory for running Celera Assembly, the content of which will be deleted inside the program. Give separate "$tmp_dir" for assembling clusters in Union-Find and Infomap to avoid potential data race. "$time_limit" is the maximum allowed time (minutes) for running assembly for the long reads in one cluster. A reasonable $time_limit is between 5min to 60min. "$cov_threshold" is the same as that in step 3.4. "$cmd_output_file" contains the assembly commands generated by this program, further to be embedded in jobs submitted to the servers. $thread is the thread number for Celera Assembly (the jobs of assembly submitted to the server have to have the same number of threads claimed). Inside this file, a separator "#" is in place between every 500 assemblies for parallelization.

Note: all prefix, directory and files in the input for this step should be in full path.

With assembly commands generated and saved in "$cmd_output_file", submit jobs to clusters. 

###6. Infer SV.###

6.1 Align contigs to the reference.

Concatenate the assembled contigs in step 5 into two fasta files, "$union_find_ctg" and "infomap_ctg" for the clusters generated in 4.1 and 4.2, respectively, and distinguish contig names among union_find and infomap clusters:

For ctg_output_dir corresponding to union_find cluster:

```cat $ctg_output_dir/*.fa | perl -ane 'if($\_ =~ /^>/){$\_ =~ s/ctg/ctg0/g; print $\_}else{print $\_}' > $union_find_ctg

For ctg_output_dir corresponding to infomap cluster:

```cat $ctg_output_dir/*.fa | perl -ane 'if($\_ =~ /^>/){$\_ =~ s/ctg/ctg1/g; print $\_}else{print $\_}' > $infomap_ctg

Align them to the reference and output in sam format. 

```blasr -maxAnchorsPerPosition 100 –advanceExactMatches 10 –affineAlign –affineOpen 100 –affineExtend 0 –insertion 5 –deletion 5 –extend –maxExtendDropoff 20 –clipping subread –bestn 3 -nproc $nproc -sam -out $union_find_sam $union_find_ctg $ref```

```blasr -maxAnchorsPerPosition 100 –advanceExactMatches 10 –affineAlign –affineOpen 100 –affineExtend 0 –insertion 5 –deletion 5 –extend –maxExtendDropoff 20 –clipping subread –bestn 3 -nproc $nproc -sam -out $infomap_sam $infomap_ctg $ref```

$ref is the reference against which the short read was aligned to create the bam. $nproc is the number of processors to be used. Note that the corresponding jobs submitted to the server should have at least $nproc processors claimed. $union_find_sam and $infomap_sam are the output alignment file in SAM format, to be further analyzed for inferring SVs in step 6.3.

For potential usage of the assembled contig, concatenate the fasta files from union_find and infomap:

```cat $union_find_ctg $infomap_ctg > all.ctg.fa```

6.2 Reorder Illumina reads according to clusters.

```reorder_short.pl <union_find_prefix> <short_prefix> <reordered_union_find>```

and 

```reorder_short.pl <infomap_prefix> <short_prefix> <reodered_infomap>```

"$union_find_prefix" and "$infomap_prefix" are the same arguments as those in 4.1 and 4.2. "$short_prefix" is the same as "$new_prefix" in 1.3. "$reordered_union_find" and "$reordered_infomap" are output files to be used in the next step.

6.3 Infer SV.

Make commands for inferring SV:

```make_inferSV.sh <union_find_dir> <union_find_sam> <reference> <union_find_prefix> <reordered_union_find> <union_find_ctg> <pairHMM_trans_matrix_file>```

```make_inferSV.sh <infomap_dir> <infomap_sam> <reference> <infomap_prefix> <reordered_infomap> <infomap_ctg> <pairHMM_trans_matrix_file>```

"$union_find_dir" and "$infomap_dir" are the output folder storing the SV results for the clusters from 4.1 and 4.2, respectively. "$union_find_sam" and "$infomap_sam" are the output of 6.1. "$union_find_prefix" and "$infomap_prefix" are the same arguments as those in 4.1 and 4.2. "$reordered_union_find" and "$reorderred_infomap" are the output in 6.2. "$union_find_ctg" and "$infomap_ctg" are the fasta files concatenated in 6.1. "$pairHMM_trans_matrix_file" is the location of the transition probability matrix (A file that contains a 3*3 matrix, which represents three states: Match, Deletion and Insertion from left to right and top to bottom. An entry at (i,j) represents a transition probability from state i to state j. ) An example of "$pairHMM_trans_matrix_file" could be found in "$DIR/HybridAssemblySV/perl_modules/pairHMM/hmm_trans.ctg.tab.txt". 

The output is a file containing the commands for interring SVs: the "run_inferSV.sh" inside "$union_find_dir" and "$infomap_dir", respectively. 

With the commands for inferring SVs generated in "$union_find_dir/run_inferSV.sh" and "$infomap_dir/run_inferSV.sh" (a separator "#" is placed for every 100 commands), submit jobs to run these commands. The output are the small SV files in the format "out.\d+_\d+.txt"

###Complex mode: 
Run the same command but with the extra input argument indicating complex mode (value 1):

```make_inferSV.sh <union_find_dir> <union_find_sam> <reference> <union_find_prefix> <reordered_union_find> <union_find_ctg> <pairHMM_trans_matrix_file> 1```

```make_inferSV.sh <infomap_dir> <infomap_sam> <reference> <infomap_prefix> <reordered_infomap> <infomap_ctg> <pairHMM_trans_matrix_file> 1```

###7. Summarize results. ###
Concatenate all the output SV files into one file:

```cat $union_find_dir/out.*_*.txt $infomap_dir/out.*_*.txt > out.all.txt```

###Complex mode

```grep largei out.all.txt > complex.txt```

The output of complex calls is in complex.txt, ready for further analysis (seen in scripts/NA12878/complex_del/).

Convert the output file to VCF file:

```txt2vcf.pl out.all.txt out.all.vcf```

Partition calls by size and type, followed by deduplication:

```get_analysis_files.pl```

Results are:

1. large deletion (>50bp): **out.del.GT50.rmdup.bed**

2. large insertion (>50bp with 50bp two-sided padding): **out.ins.GT50.flanking50.rmdup.bed**

3. small deletion (<=50bp): **out.del.SET50.rmdup.bed**

4. small insertion (<=50bp with 50bp two-sided padding): **out.ins.SET50.flanking50.rmdup.bed**
