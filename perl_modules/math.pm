use warnings;
use strict;
package math;

sub rolling_sum{
    my ($a, $num) = @_;
    my @b;
    my $s = 0;
    my $n = 0;
    foreach (@$a){
        if($n < $num){
            $n ++;
            $s += $_;
        }
        else{
            push @b, $s;
            $s = $_;
            $n = 1;
        }
    }
    return \@b;
}

sub sum_ab {
    my ($a) = @_;
    my $s = 0;
    foreach (@$a){
        $s += abs($_);
    }
    return $s;
}


sub sum {
    my ($a) = @_;
    my $s = 0;
    foreach (@$a){
        $s += $_;
    }
    return $s;
}

sub avg_ab {
    my ($a) = @_;
    my $t = scalar(@$a);
    if($t == 0){
        return "NA";
    }
    return &sum_ab($a)/$t;
}


sub avg {
    my ($a) = @_;
    my $t = scalar(@$a);
    if($t == 0){
        return "NA";
    }
    return &sum($a)/$t;
}

sub max {
    my ($a) = @_;
    my $s = 0;
    foreach (@$a){
        if($_ > $s){
            $s = $_;
        }
    }
    return $s;
}

# return max and the index
sub max_w_index {
    my ($a) = @_;
    my $s = 0;
    my $index;
    foreach my $i (0 .. scalar(@$a) - 1){
        if($a->[$i] > $s){
            $index = $i;
            $s = $a->[$i];
        }
    }
    return ($s, $index);
}

sub std {
    my ($a) = @_;
    my $t = scalar(@$a);
    my $s = 0;
    my $avg = &avg($a);
    foreach (@$a){
        $s = $s + ($_ - $avg)**2;
    }
    return sqrt($s/$t);
}

1;
