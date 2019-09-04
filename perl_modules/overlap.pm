#!/bin/perl -w
use warnings;
use strict;
require tmpfile;

package overlap;

my $intersectBed = "/scratch/bcb/xfan3/pkg/bedtools/bedtools-2.17.0/bin/intersectBed";

# Function: Compare the overlap between hash A and hash B. Hash in a form of hash->{chr}->{start.end}. Output the FP and FN of B relative to A.

sub run{
    my ($hash_a, $hash_b, $f, $r, $v, $wa, $wb, $wo) = @_;
    $DB::single = 1;
    &write_tmpbed($hash_a, "tmp.a");
    &write_tmpbed($hash_b, "tmp.b");
    my $r_ = "";
    $r_ = "-r" if($r == 1);
    my $f_ = 0.5;
    $f_ = "-f $f" if($f <= 1 && $f >= 0);
    my $v_ = "";
    $v_ = "-v" if($v == 1);
    my $wa_ = "";
    $wa_ = "-wa" if($wa == 1);
    my $wb_ = "";
    $wb_ = "-wb" if($wb == 1);
    my $wo_ = "";
    $wo_ = "-wo" if($wo == 1);

    my $time = tmpfile::time();
    my $outputfile = "overlapped.$time.bed";
    #print "#$intersectBed -a tmp.a -b tmp.b $r_ $f_ $v_ $wa_ $wb_ $wo_ > $outputfile\n";
    `$intersectBed -a tmp.a -b tmp.b $r_ $f_ $v_ $wa_ $wb_ $wo_ > $outputfile`;
    return $outputfile;
}

# overlap itx, hash in a form of hash->{chr}->{s.e} = type
sub run_itx{
    my ($hash_a, $hash_b) = @_;
    my ($f, $r, $v, $wa, $wb, $wo) = (0.8, 1, 0, 0, 0, 1);
    &print_itx($hash_a, $hash_b, $f, $r, $v, $wa, $wb, $wo);
}

# overlap ctx, hash in a form of hash->{chr1:pos1} = chr1:pos1:chr2:pos2, the two ends should be both overlapped, a flanking region is placed
sub run_ctx{
    my ($hash_a, $hash_b, $flank) = @_;
    # prepare hash from hash_a and hash_b
    my ($hash_1, $hash_2);
    $hash_1 = &prepare_ctx($hash_a, $flank);
    $hash_2 = &prepare_ctx($hash_b, $flank);

    # run overlap
    my $somatic = 0;
    my ($f, $r, $v, $wa, $wb, $wo) = (0.5, 1, 2, 2, 2, 1);
    &print_ctx($hash_1, $hash_2, $f, $r, $v, $wa, $wb, $wo, $somatic);
}

sub print_ctx{
    my ($hash_1, $hash_2, $f, $r, $v, $wa, $wb, $wo, $somatic) = @_;

    my $outputfile = &run($hash_1, $hash_2, $f, $r, $v, $wa, $wb, $wo);


    my $hash;
    open fh_, "<$outputfile" or die $!;
    while(<fh_>){
        chomp;
        my @b = split(/\t/, $_);
        if($somatic == 1){
            my @a = split(/\:/, $b[3]);
            print join("\t", @a[0 .. 1], "+", @a[2 .. 3], "-", "CTX", "99", "99")."\n";
        }
        else{
            $hash->{$b[3]}->{$b[7]}->{$b[0]} = 1;
        }
    }
    close fh_;

    if($somatic != 1){
        foreach my $key (keys %$hash){
            foreach my $key2 (keys %{$hash->{$key}}){
                if(scalar(keys %{$hash->{$key}->{$key2}}) == 2){
                    if($key =~ /:P/ || $key =~ /:N/){
                        my $str = &print_specific($key, $key2);
                        print $str;
                        last;
                    }
                    #print $key . "\t" . $key2 . "\n";
                    my @tmp = split(/\:/, $key2);
                    print join("\t", @tmp[0 .. 1], "+", @tmp[2 .. 3], "-", "CTX", "99", "99")  . "\n";
                    last;
                }
            }
        }
        #print "#In all, $num are recognized out of $total_num SVs.\n";
    }
}

# read a key with the format ";", each entry "col:entry:[P|N]"
sub read_key{
    my ($str, $hash) = @_;
    my @a = split(/;/, $str);
    foreach (@a) {
        my @b = split(/:/, $_);
        if($b[2] =~ /P/){
            $hash->{$b[0]} = $b[1];
        }
    }
    return $hash;
}

# print in a user defined format combining two strings, each in a format of col1:$$:P;    col2:$$:N; (P is print, N is not print, $$ stand for the entry, and col1 col2 are the     column index
sub print_specific{
    my ($key1, $key2) = @_;
    my $hash;
    $hash = &read_key($key1, $hash);
    $hash = &read_key($key2, $hash);
#    my $previous = -1;
    my @a;
    foreach my $col_id (sort {$a <=> $b} keys %$hash){
#        if($col_id - $previous > 1){
#            push @a, "NA";
#        }
#        else{
            push @a, $hash->{$col_id};
#        }
#        $previous = $col_id;
    }
    return join("\t", @a) . "\n";
}

# overlap itx, hash in a form of hash->{chr1}->{start.end}, 0.5 criteria is used
sub run_itx_somatic{
    my ($hash_t, $hash_n) = @_;
    my ($f, $r, $v, $wa, $wb, $wo) = (0.5, 1, 1, 0, 0, 1);
    &print_itx($hash_t, $hash_n, $f, $r, $v, $wa, $wb, $wo);
}


sub print_itx{
    my ($hash_t, $hash_n, $f, $r, $v, $wa, $wb, $wo) = @_;
    my $outputfile = &run($hash_t, $hash_n, $f, $r, $v, $wa, $wb, $wo);
    my $hash;
    # for removing duplication
    my ($hash1, $hash2);
        $DB::single = 1;
    open fh_, "<$outputfile" or die $!;
    while(<fh_>){
        my @a = split(/\t/, $_);
        my ($key, $value);
        # to avoid multiple to one and one to multiple
        if(defined $a[5]){
            my $key1 = join(":", @a[0 .. 2]);
            my $key2 = join(":", @a[4 .. 6]);
            next if(defined $hash1->{$key1});
            next if(defined $hash2->{$key2});
            # setting threshold of at least 500bp to known start and end, respectively
            next if(abs($a[1] - $a[5]) > 500 || abs($a[2] - $a[6]) > 500);
            $hash1->{$key1} = 1;
            $hash2->{$key2} = 1;
        }
        $key = &print_specific($a[3], $a[7]) if($a[3] =~ /:P/ || $a[7] =~ /:P/);
        $key = join("\t", @a[0 .. 1], "+", $a[0], $a[2], "-", $hash_t->{$a[0]}->{"$a[1].$a[2]"}, $a[2] - $a[1], "99") . "\n" if($a[3] !~ /:P/ && $a[7] !~ /:P/);
        $hash->{$key} = 1;
    }
    close fh_;
    foreach my $key (sort keys %$hash){
        print $key;
    }
}


# overlap ctx, hash in a form of hash->{chr1:pos1} = chr1:pos1:chr2:pos2, the two ends should be both overlapped, a flanking region is placed
sub run_ctx_somatic{
    my ($hash_t, $hash_n, $flank) = @_;
    # prepare hash from hash_a and hash_b
    my ($hash_1, $hash_2);
    $hash_1 = &prepare_ctx($hash_t, $flank);
    $hash_2 = &prepare_ctx($hash_n, $flank);
    # run overlap
    my ($f, $r, $v, $wa, $wb, $wo) = (0.5, 1, 1, 2, 2, 1);
    my $outputfile = &run($hash_1, $hash_2, $f, $r, $v, $wa, $wb, $wo);

    my $somatic = 1;
    &print_ctx($hash_1, $hash_2, $f, $r, $v, $wa, $wb, $wo, $somatic);
}

sub prepare_ctx{
    my ($h, $flank) = @_;
    my $hash;
    foreach my $key (keys %$h){
        my ($chr, $pos) = split(/\:/, $key);
        my ($s, $e) = ($pos - $flank, $pos + $flank);

        $s = 1 if($s <= 0);
        my $str = "$s.$e";
        $hash->{$chr}->{$str} = $h->{$key};
    }
    return $hash;
}

sub write_tmpbed{
    my ($hash, $file) = @_;
    open file_fh, ">$file" or die $!;
    foreach my $chr (keys %$hash){
        foreach my $pos (keys %{$hash->{$chr}}){
            my ($a, $b) = split(/\./, $pos);
            print file_fh join("\t", ($chr, $a, $b, $hash->{$chr}->{$pos}))."\n";
        }
    }
    close file_fh;
}


