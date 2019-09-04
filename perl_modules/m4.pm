use warnings;
use strict;
package m4;

# read m4, focus on the relation between IL and PB
# not take care of read name changes, etc
sub read_m4{
    my ($m4) = @_;
    my $h;
    open fh_, "<$m4" or die $!;
    while(<fh_>){
        my @a = split(/\s+/, $_);
        $h->{$a[0]}->{$a[1]} = 1;
    }
    close fh_;
    return $h;
}

# This function takes a m4, a list of IL read names, each at one line, a list of PB read names, each at one line, and output a new m4 which contains lines only when both IL and PB are included in the two files.
sub get_new_m4{
    my ($m4, $new_m4, $IL_f, $PB_f) = @_;
    my $h_il;
    my $h_pb;
    open il_fh, "<$IL_f" or die $!;
    while(<il_fh>){
        chomp;
        $h_il->{$_} = 1;
    }
    close il_fh;
    open pb_fh, "<$PB_f" or die $!;
    while(<pb_fh>){
        chomp;
        $h_pb->{$_} = 1;
    }
    close pb_fh;
    open m4_fh, "<$m4" or die $!;
    open new_m4_fh, ">$new_m4" or die $!;
    while(<m4_fh>){
        my @a = split(/\s+/, $_);
        my @b = split(/\//, $a[0]);
        if(defined $h_il->{$b[0]} && defined $h_pb->{$a[1]}){
            print new_m4_fh $_;
        }
    }
    close m4_fh;
    close new_m4_fh;
}
# This function takes a m4, a list of PB and ILs in one file, output the new m4
sub get_new_m4_from_out{
    my ($m4, $new_m4, $out) = @_;
    my $h_il;
    my $h_pb;
    open fh_, "<$out" or die $!;
    while(<fh_>){
        next if($_ =~ /^#/);
        chomp;
        if($_ !~ /\./){
            $h_il->{$_} = 1;
        }
        else{
            $h_pb->{$_} = 1;
        }
    }
    close fh_;
    open m4_fh, "<$m4" or die $!;
    open new_m4_fh, ">$new_m4" or die $!;
    while(<m4_fh>){
        my @a = split(/\s+/, $_);
        my @b = split(/\//, $a[0]);
        if(defined $h_il->{$b[0]} && defined $h_pb->{$a[1]}){
            print new_m4_fh $_;
        }
    }
    close m4_fh;
    close new_m4_fh;
}

# given a m4, and a file containing the nodes to be removed, output a new m4 that does not contain these nodes
sub remove_node{
    my ($m4, $nodes_to_rm, $new_m4) = @_;
    my $h;
    open node_fh, "<$nodes_to_rm" or die $!;
    while(<node_fh>){
        chomp;
        $h->{$_} = 1;
    }
    close node_fh;
    open m4_fh, "<$m4" or die $!;
    open new_m4_fh, ">$new_m4" or die $!;
    while(<m4_fh>){
        my @a = split(/\s+/, $_);
        my @b = split(/\//, $a[0]);
        if(defined $h->{$b[0]} || defined $h->{$a[1]}){
            next;
        }
        print new_m4_fh $_;
    }
    close m4_fh;
    close new_m4_fh;
}

# given a m4, and a file containing the edges to be removed, output a new m4 that does not contain these edges
sub remove_edge{
    my ($m4, $edges_to_rm, $new_m4) = @_;
    my $h;
    open edge_fh, "<$edges_to_rm" or die $!;
    while(<edge_fh>){
        chomp;
        my @a = split(/\t/, $_);
        $h->{$a[0]}->{$a[1]} = 1;
    }
    close edge_fh;
    open m4_fh, "<$m4" or die $!;
    open new_m4_fh, ">$new_m4" or die $!;
    while(<m4_fh>){
        my @a = split(/\s+/, $_);
        my @b = split(/\//, $a[0]);
        if(defined $h->{$b[0]}->{$a[1]}){
            next;
        }
        print new_m4_fh $_;
    }
    close m4_fh;
    close new_m4_fh;
}

# This function takes a m4, a reported breakpoint (for INDEL, in format of chr, pos, length, I/D, contig_name, starting position of this contig on the reference, and the contig orientation w.r.t. the reference), and report the line if there are at least m short read pairs intersecting this breakpoint
sub check_INDEL_bp{
    my ($m4, $line, $m, $min_no_check) = @_;
    my $h;
    my $n = 0;
    foreach my $l (@$line){
        chomp $l;
        my @a = split(/\t/, $l);
        if($a[2] >= $min_no_check){
            print $l . "\n";
            next;
        }
        my $hh;
        $hh->{chr} = $a[0];
        $hh->{pos} = $a[1];
        $hh->{len} = $a[2];
        $hh->{tp} = $a[3];
        $hh->{ctg_name} = $a[4];
        #ctg_start is the coordinate on the contig, in the orientaiton of the contig
        $hh->{ctg_start} = $a[5];
        $hh->{ctg_end} = $a[6];
        $hh->{spt_short} = 0;
        $hh->{line} = $l;
        $h->{$n} = $hh;
        $n ++;

    }
    open M4, "<$m4" or die $!;
    while(<M4>){
        my ($rn_q, $rn_r, $as, $a_per, $ori_q, $s_q, $e_q, $len_q, $ori_r, $s_r, $e_r, $len_r, ) = split(/\s+/, $_);  
        my $next = <M4>;
        my @b = split(/\s+/, $next);
        my ($s_r1, $s_q1) = ($b[9], $b[5]);

        # compensate the clipped ends
        $s_r -= $s_q;
        $s_r1 -= $s_q1;
        my $start;
        my $end;
        if($ori_r == 0){
            $end = $len_r - $s_r1;
            $start = $s_r;
        }
        elsif($ori_r == 1){
            $end = $len_r - $s_r;
            $start = $s_r1;
        }
        # convert to the coordinate on the reference
        foreach my $i (keys %$h){
            if($h->{$i}->{ctg_name} eq $rn_r){
                if($h->{$i}->{ctg_end} > $start && $h->{$i}->{ctg_start} < $end){
                    # within the short read covered region
                    $h->{$i}->{spt_short} ++;
                }
            }
        }
    }
    close M4;

    # count those with suppporting shorts > m, then report it
    my $p;
    foreach my $i (keys %$h){
        if($h->{$i}->{spt_short} >= $m){
            my $ll = $h->{$i}->{line};
            $p->{$ll} = 1;
        }
    }
    foreach my $key (keys %$p){
        print $key . "\n";
    }
}

1;
