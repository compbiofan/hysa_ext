use warnings;
use strict;

package string;

# This package deal with strings of ACGT

sub to_capital{
    my $str = shift;
    $str =~ tr/a-z/A-Z/;
    return $str;
}

sub to_small{
    my $str = shift;
    $str =~ tr/A-Z/a-z/;
    return $str;
}

sub reverse_comp{
    my ($str) = shift;
    my $str_rev = scalar reverse("$str");
    $str_rev =~ tr/ACGT/TGCA/;
    return $str_rev;
}

# check if the two strings are the same, after converting both to capital, considering reverse complement
# check if string1 is a subset of string2, off by off1 and off2 at left and right, after converting both to capital, considering reverse complement
sub check_same{
    my ($str1, $str2, $off1, $off2) = @_;
    my $str1_ = &to_capital($str1);
    my $str2_ = &to_capital($str2);
    return 1 if($str1_ eq $str2_ || &reverse_comp($str1_) eq $str2_);
    return 0 if(!defined $off1 || !defined $off2);

    my ($str1_1, $str1_2) = (substr($str1_, abs($off1), length($str1) - abs($off2) - abs($off1)), substr($str1_, abs($off2), length($str1) - abs($off1) - abs($off2)));

    return 2 if($str1_1 eq $str2_ || &reverse_comp($str1_2) eq $str2_);
 
    return 0;
}

1;
