#!/usr/bin/perl
use warnings;
use strict;

require sam;

if($#ARGV == -1){
    die "This converts a sam file to a fa pair. \nUsage:perl $0 <sam_file>\n"; 
}

sam::sam2fa_pair(@ARGV);
