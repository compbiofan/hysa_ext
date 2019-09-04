#!/risapps/rhel6/perl/5.10.1/bin/perl
use warnings;
use strict;

use bed;

if($#ARGV == -1){
    die "This converts a bed file into a vcf file. \nUsage: bed2vcf.pl <bed_file> <vcf_file>\n";
}

bed::bed2vcf(@ARGV);
