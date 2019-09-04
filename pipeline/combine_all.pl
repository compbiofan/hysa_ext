#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;

if(@ARGV == 0){
    die "This combines all sam files into one, reorder the reads, and generate one new sam and index file in the new directory. \nUsage: $0 <contig_file:sam in each row> <new_dir> <prefix>\n";
}

my $path = dirname($0);
my ($config_file, $new_dir, $prefix, $head) = @ARGV;
my $tag = 1;
my $id = -1;
`mkdir $new_dir` if(! -d $new_dir);
my $new_index = "$new_dir/$prefix.id";
my $new_sam = "$new_dir/$prefix.sam";
my $new_drp = "$new_dir/$prefix.drp";
`rm $new_index` if(-e $new_index);
`rm $new_sam` if(-e $new_sam);
`rm $new_drp` if(-e $new_drp);
open index_fh, ">>$new_index" or die $!;
open sam_fh, ">>$new_sam" or die $!;
open drp_fh, ">>$new_drp" or die $!;
my $h = &read_cfg($config_file);
#foreach my $sam (sort keys %$h){
#    my $index = $h->{$sam};
#    &read_index($index);
#}
foreach my $sam (sort keys %$h){
    my $index = $h->{$sam}->{id};
    my $drp = $h->{$sam}->{drp};
    my $id_hash = &read_index($index);
    my $drp_hash = &read_drp($drp) if($drp ne "NA");
    open fh_, "<$sam" or die $!;
    while(<fh_>){
        chomp;
        my @a = split(/\t/, $_);
        if($tag == 0){
            $tag = 1;
            if($a[0] =~ /\d+_ref$/){
                print sam_fh join("\t", "${id}_ref", @a[1 .. $#a]) . "\n";
                next;
            }

        }
        elsif($tag == 1){
            $id ++;
            $tag = 0;
            if(!defined $id_hash->{$a[0]} && $a[0] !~ /_ref$/){
                $id_hash->{$a[0]} = $a[0];
            }
            elsif($a[0] =~ /(\d+)_ref$/){
                $id --;
                print sam_fh join("\t", "${id}_ref", @a[1 .. $#a]) . "\n";
                next;
            }

            my $index_line = join("\t", $id_hash->{$a[0]}, $id) . "\n";
            print index_fh $index_line;
            if($drp ne "NA" && defined $drp_hash->{$a[0]}){
                my $drp_line = join("\t", $id, $drp_hash->{$a[0]}) . "\n";
                print drp_fh $drp_line;
            }
        }
        my $line = join("\t", $id, @a[1 .. $#a]) . "\n";
        print sam_fh $line;
    }
    close fh_;
}
close index_fh;
close sam_fh;
close drp_fh;
`$path/sam2fa_pair.sh $new_dir/$prefix`;
#`~/bin/makebam.3.sh $new_dir/$prefix $head`;

1;

# read the drp for each small sam, save to a hash with the keys the original index in small sam, and the values the insert size
sub read_drp{
    my ($drp_index) = @_;
    my $hash;
    open discordant_fh, "<$drp_index" or die $!;
    while(<discordant_fh>){
        chomp;
        my @a = split(/\t/, $_);
        if(!defined $a[1]){
            print $drp_index . "\t$_\n";
        }
        $hash->{$a[0]} = $a[1];
    }
    close discordant_fh;
    return $hash;
}


# this is to read the index file for small sams, so that in the big index file, the name is readname, not the small index file's id
sub read_index{
    my ($index) = @_;
    my $hash;
    open id_fh, "<$index" or die $!;
    while(<id_fh>){
        chomp;
        my @a = split(/\t/, $_);
        if(!defined $a[1]){
            print $index . "\t$_\n";
        }
        $hash->{$a[1]} = $a[0];
    }
    close id_fh;
    return $hash;
}

sub read_cfg{
    my ($cfg) = @_;
    my $h;
    open cfg_fh, "<$cfg" or die $!;
    while(<cfg_fh>){
        chomp;
        my @a = split(/\s+/, $_);
        $h->{$a[0]}->{id} = $a[1];
        $h->{$a[0]}->{drp} = $a[2];
    }
    close cfg_fh;
    return $h;
}
