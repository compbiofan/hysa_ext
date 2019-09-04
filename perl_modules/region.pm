use warnings;
use strict;
package region;


# a package about regions

# convert from chr:start-end to (chr, start, end)
sub convert_region{
    my ($str) = @_;
    my @a = split(/:/, $str);
    my @b = split(/\-/, $a[1]);
    return ($a[0], $b[0], $b[1]);
}

sub overlap_w_flank{
    my ($r1, $r2, $flank) = @_;
    my ($c1, $s1, $e1) = ($r1 =~ /^(.+):(\d+)-(\d+)$/);
    my ($c2, $s2, $e2) = ($r2 =~ /^(.+):(\d+)-(\d+)$/);
    $s1 -= $flank;
    $s2 -= $flank;
    $e1 += $flank;
    $e2 += $flank;
    $s1 = 0 if($s1 < 0);
    $s2 = 0 if($s2 < 0);
    return 0 if($c1 ne $c2);
    return 0 if($e1 < $s2 || $e2 < $s1);
    return 1;
}

# overlap two regions in the format of chr:start-end, return 1 if yes, 0 if no
sub overlap{
    my ($r1, $r2) = @_;
    my ($c1, $s1, $e1) = ($r1 =~ /^(.+):(\d+)-(\d+)$/);
    my ($c2, $s2, $e2) = ($r2 =~ /^(.+):(\d+)-(\d+)$/);
    return 0 if($c1 ne $c2);
    return 0 if($e1 < $s2 || $e2 < $s1);
    return 1;

}

# the second argument can be an array of regions, return 1 as long as at least has overlap
sub overlap1{
    my ($r1, $s) = @_;
    my @r2s = @$s;
    foreach my $r2 (@r2s){
        if(&overlap($r1, $r2) == 1){
            return 1;
        }
    }
    return 0;
}
1;
