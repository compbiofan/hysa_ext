#!/usr/bin/perl
use warnings;
use strict;
require cigar;

if($#ARGV == -1){
    die "Split the file into chunks. Read the paired chunks at the same time. If a read is paired up, output those concordants and dump the disconcordants. Otherwise store it until its pair appears in a chunk. Assume one read occur in one chunk only for one end. Function: filter/pair up alignments from blasr. Starting 08252015, now accepts blasr alignment result in sam format. Output will be in the same format as the input. \nUsage: $0 <output_dir> <perc_t> <alignLen_t> <insert_t> <m4_prefix> <chunk_read_num>\n";
}

my ($output_dir, $perc_t, $alignLen_t, $insert_t, $m4_prefix) = @ARGV;
# by default 10000 reads
my $chunk_size = 10000;
if(defined $ARGV[5]){
    $chunk_size = $ARGV[5];
}
my $h;
# changed on 06292016
my $out = "$output_dir/${m4_prefix}.m4"; 
#my $out = $output_dir . "/out.m4";
open OUT, ">$out" or die $!;
close OUT;

# begin reading by chunk
my $f1 = "$output_dir/${m4_prefix}_1.m4";
my $f2 = "$output_dir/${m4_prefix}_2.m4";
open M4_1, "<$f1" or die $!;
open M4_2, "<$f2" or die $!; 
while(eof(M4_1) != 1 || eof(M4_2) != 1){
    $h = &readAln_m4($h, $chunk_size); 
    open OUT, ">>$out" or die $!;
    $h = &checkAln_m4($h);
    close OUT;
}

close M4_1;
close M4_2;

1;
# read both chunk_size of the following M4 for 1 and 2
sub readAln_m4{
    my ($h, $chunk_size) = @_;
    my $h_num1;
    my $h_num2;
    while(<M4_1>){
        next if($_ =~ /^@/);
        my ($rn_q, $rn_r, $as, $a_per, $ori_q, $s_q, $e_q, $len_q, $ori_r, $s_r, $e_r, $len_r, ) = split(/\s+/, $_);
        my @tmp = split(/\//, $rn_q);
        $rn_q = join("/", @tmp[0 .. $#tmp - 1]) if($rn_q =~ /\//);
        # control the size of chunk by the number of reads
        $h_num1->{$rn_q} = 1;
        if(scalar(keys %$h_num1) >= $chunk_size){
            # retrive back one line
            seek(M4_1, -length($_), 1);
            last;
        }
        next if($a_per < $perc_t);
        # as one_mapped_c, check the insert size and one of the mapped length
        #next if($e_q - $s_q < $alignLen_t);

        if($ori_r == 1){
            # reverse, the pos is $len - $end on the reference
            $s_r = $len_r - $e_r;
        }
        $h->{$rn_q}->{1}->{$rn_r}->{"$ori_r.$s_r"} =  &change_rn($_, 1);
    }
    while(<M4_2>){
        next if($_ =~ /^@/);
        my ($rn_q, $rn_r, $as, $a_per, $ori_q, $s_q, $e_q, $len_q, $ori_r, $s_r, $e_r, $len_r, ) = split(/\s+/, $_);
        my @tmp = split(/\//, $rn_q);
        $rn_q = join("/", @tmp[0 .. $#tmp - 1]) if($rn_q =~ /\//);
        # control the size of chunk by the number of reads
        $h_num2->{$rn_q} = 1;
        if(scalar(keys %$h_num2) >= $chunk_size){
            # retrive back one line
            seek(M4_2, -length($_), 1);
            last;
        }
        next if($a_per < $perc_t);
        # as one_mapped_c, check the insert size and one of the mapped length
        #next if($e_q - $s_q < $alignLen_t);

        if($ori_r == 1){
            # reverse, the pos is $len - $end on the reference
            $s_r = $len_r - $e_r;
        }
        $h->{$rn_q}->{2}->{$rn_r}->{"$ori_r.$s_r"} = &change_rn($_, 2);
    }
    return $h;
}
# change the read name so that the first and seoncd read can be distinguished
sub change_rn{
    my ($line, $tag) = @_;
    chomp $line;
    my @aa = split(/\s+/, $line);
    return join(" ", $aa[0] . "/$tag", @aa[1 .. $#aa]);
}

sub checkAln_m4{
    my ($h) = @_;
    my @to_delete;
    foreach my $read (keys %$h){
        next if(keys %{$h->{$read}} == 1);
        my $h1 = $h->{$read}->{1};
        my $h2 = $h->{$read}->{2};
        foreach my $pb (keys %$h1){
            if(defined $h2->{$pb}){
                # aligned to the same pb
                # check orientation, insert size and aligned length
                my $h1pb = $h1->{$pb};
                my $h2pb = $h2->{$pb};
                foreach my $p1 (keys %$h1pb){
                    my ($ori1, $pos1) = split(/\./, $p1);
                    foreach my $p2 (keys %$h2pb){
                        my ($ori2, $pos2) = split(/\./, $p2);
                        next if($ori1 == $ori2 || abs($pos2 - $pos1) > $insert_t || $ori1 == 0 && $pos1 >= $pos2 || $ori1 == 1 && $pos1 <= $pos2);
                        # now qualify in pair, need to check the alignment length
                        my @x = split(/\s+/, $h1pb->{$p1});
                        my @y = split(/\s+/, $h2pb->{$p2});
                        if(abs($x[6] - $x[5]) >= $alignLen_t || abs($y[6] - $y[5]) >= $alignLen_t){
                            print OUT join("\n", $h1pb->{$p1}, $h2pb->{$p2}) . "\n";
                        }
                    }
                }
            }
        }
        push @to_delete, $read;
    }
    foreach my $read (@to_delete){
        delete $h->{$read};
    }
    return $h;
}



