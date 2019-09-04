# Supporting scripts for analysis of CHM1 large insertion.

## Steps

```./run_validate.sh <HySA_output> <HySA_assembly> <build37> <build38> <asm1> <asm2> <nproc>```

$HySA_output is the output of HySA in bed format after deduplication (out.ins.GT50.flanking50.rmdup.bed in pipeline). $HySA_assembly is a fasta file with HySA assembled contig (all.ctg.fa in pipeline). $build37, $build38, $asm1 and $asm2 are all fasta files. $asm1 is the insertion file in [1], downloaded from http://eichlerlab.gs.washington.edu/publications/chm1-structural-variation/data/GRCh37/insertions.bed. $asm2 is the whole genome assembly of CHM1 in [2], downloaded from ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000772585.3_ASM77258v3/GCA_000772585.3_ASM77258v3_genomic.fna.gz and gunzipped. 

## Reference

1. Chaisson MJP, Huddleston J, Dennis MY, Sudmant PH, Malig M, Hormozdiari F, Antonacci F, Surti U, Sandstrom R, Boitano M et al. 2015b. Resolving the complexity of the human genome using single-molecule sequencing. Nature 517: 608-U163.
2. Berlin K, Koren S, Chin CS, Drake JP, Landolin JM, Phillippy AM. 2015. Assembling large genomes with single-molecule sequencing and locality-sensitive hashing. Nat Biotechnol 33: 623-630.
