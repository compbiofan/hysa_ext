use warnings;
use strict;
#require List::Util::sum;
#require POSIX::ceil;
package asm;

# Some functions for assembly evaluation.

sub calculate_N50{
    my ($seq_len) = @_;
    my @new_l;
    foreach (@$seq_len){
        foreach (1 .. $_){
            push @new_l, $_;
        }
    }
    $DB::single = 1;
    return median(\@new_l);
}

sub median {
    my ($len) = @_;
    my @len_sort = sort {$a <=> $b} @$len;
    return $len_sort[int(scalar(@len_sort)/2)] if(scalar(@len_sort) % 2 == 1);
    return ($len_sort[scalar(@len_sort)/2] + $len_sort[scalar(@len_sort)/2 - 1])/2 if(scalar(@len_sort) % 2 == 0);
#    sum( ( sort { $a <=> $b } @_ )[ int( $#_/2 ), ceil( $#_/2 ) ] )/2;
}
1;
