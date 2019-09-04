# Supporting scripts for analysis of NA12878 large insertion.

## Steps

```./run_validate.sh <HySA_output> <HySA_assembly> <build38> <asm> <fosmid> <nproc>```

$HySA_output is the output of HySA in bed format after deduplication (out.ins.GT50.flanking50.rmdup.bed). $HySA_assembly is a fasta file with HySA assembled contig (all.ctg.fa in pipeline). $build38, $asm1 and $fosmid are all fasta files. $asm is the whole genome assembly of NA12878 in [1], downloaded from ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_001013985.1_ASM101398v1/GCA_001013985.1_ASM101398v1_genomic.fna.gz and gunzipped. $fosmid is downloaded from http://www.cell.com/cell/fulltext/S0092-8674(10)01197-9 with GenBank Bioproject 29893 [2]. 

Note: after $fosmid is downloaded, only the fosmids corresponding to NA12878 are selected and used by the following command:

```perl -ane 'BEGIN{$tag = 0}{if($_ =~ /^>/){if($_ =~ /ABC12/){$tag = 1; print $_}else{$tag = 0;}}elsif($tag == 1){print $_}}' $fosmid_downloaded > $fosmid```

in which $fosmid_downloaded is the fasta file downloaded from the link provided above, and $fosmid is the file feeding run_validate.sh as the input.

## Reference

1. Pendleton M, Sebra R, Pang AW, Ummat A, Franzen O, Rausch T, Stutz AM, Stedman W, Anantharaman T, Hastie A et al. 2015. Assembly and diploid architecture of an individual human genome via single-molecule technologies. Nat Methods 12: 780-786.
2. Kidd JM, Cooper GM, Donahue WF, Hayden HS, Sampas N, Graves T, Hansen N, Teague B, Alkan C, Antonacci F et al. 2008. Mapping and sequencing of structural variation from eight human genomes. Nature 453: 56-64.



