#!/usr/bin/perl
use warnings;
use strict;
if($#ARGV == -1){
    die "This partitions pacbio subreads, either scattered in multiple files or in one big fasta file, into multiple batches, each batch is a directory of out_prefix + id, with the fasta inside with the size of the reference. The mode can be either 37 or 38, representing build37 and build38. \nUsage: $0 <config file || one big fasta file> <out_prefix> <mode>\n";
}

# Now the reference file contains the partition index in readname.
my ($config_file, $out_prefix, $mode) = @ARGV;
my $ref_size;
if($mode == 37){
    $ref_size = 3101804739;
}
elsif($mode == 38){
    $ref_size = 3217346917;
}
elsif($mode == 0){
    # debug mode
    $ref_size = 100000;
}
else{
    die "Mode can only be 37 or 38, representing reference version build37 or build38\n";
}

my @f;
# analyze config file
if($config_file =~ /\.fasta$/ || $config_file =~ /\.fa$/){
    # this is a big fasta file
    push @f, $config_file;
}
else{
    # this is a file with small subread fastas
    open CON, "<$config_file";
    while(<CON>){
        chomp;
        push @f, $_;
    }
    close CON;
}

# control the batch index
my $ID = 0;
my $batch_dir = $out_prefix . $ID;
`mkdir $batch_dir` if(! -d $batch_dir);
my $ID_file = "$batch_dir/id.txt";
open ID_fh, ">$ID_file" or die $!;
my $OUT_file = "$batch_dir/ref.fa";
open OUT_fh, ">$OUT_file" or die $!;
my $sz = 0;
# control the read index in each batch 
my $id = 0;
# now partition the reads according to the read base number
foreach my $ff (@f){
    open fh_, "<$ff" or die $!;
    while(<fh_>){
        my $next = $_;
        if($sz > $ref_size){
            # process the bases of the previous reads, then start a new batch
            if($next !~ /^>/){
                print OUT_fh $next;
                $next = <fh_>;
                while($next !~ /^>/){
                    print OUT_fh $next;
                    $next = <fh_>;
                }
            }


# close the previous batch if there was one
            close ID_fh;
            close OUT_fh;
# begin a new batch
            $ID ++;
            $sz = 0;
            $id = 0;
            $batch_dir = $out_prefix . $ID;
            `mkdir $batch_dir` if(! -d $batch_dir);
            $ID_file = "$batch_dir/id.txt";
            open ID_fh, ">$ID_file" or die $!;
            $OUT_file = "$batch_dir/ref.fa";
            open OUT_fh, ">$OUT_file" or die $!;
        }
        if($next =~ /^>(.+)$/){
            my $rn = $1;
# index it now
            print ID_fh "$ff:$rn\t$id\n";
            print OUT_fh ">$ID.$id\n";
            $id ++;
        }
        else{
# bases: count the size
            chomp $next;
            $sz += length($next);
            print OUT_fh $next . "\n";
        }
    }
    close fh_;
}
close ID_fh;
close OUT_fh;

