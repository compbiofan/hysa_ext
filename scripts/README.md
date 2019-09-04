# This directory contains scripts to repeat the analysis on CHM1 and NA12878.

## Software Requirements ##

1. samtools (https://samtools.github.io)
2. perl 5.10.1 or up
3. customized blasr (downloadable from https://github.com/mchaisso/blasr with commit number 82553d6)
4. bedtools 2.17.0 or up

## Environment Setup ##

Suppose "$DIR" is the path of this package.

* Make sure all perl scripts and bash scripts in pipeline are executable.

```chmod u+x $DIR/hybridassemblysv/scripts/perl_scripts/*.pl```
```chmod u+x $DIR/hybridassemblysv/scripts/bash_scripts/*.sh```

* Add perl_modules to your PERL5LIB. In bash, 

  ```export PERL5LIB=$DIR/hybridassemblysv/perl_modules```

* Add scripts to your directory. In bash,

```HySA_sup_perl="$DIR/hybridassemblysv/scripts/perl_scripts"```

```if [ -d $HySA_sup_perl ]; then PATH="$HySA_sup_perl:$PATH" fi```

```HySA_sup_bash="$DIR/hybridassemblysv/scripts/bash_scripts"```

```if [ -d $HySA_sup_bash ]; then PATH="$HySA_sup_bash:$PATH" fi```

Similarly, make intersectBed in bedtools bin directory, GATK's GenomeAnalysisTK.jar, Pindel and delly binaries executable from any location. Also, add liftover folder (both the folder containing the binary and the folder containing the chain file) in path. 

## Locations

In the directories of CHM1, there are three subfolders. They are large_del, large_ins, small_indel. In the directories of NA12878, in addition to large_del, large_ins, small_indel, there are extra two subfolders, complex_del and fosmid_validation. In each subfolder, the README is where you would start. 
