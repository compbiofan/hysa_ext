use warnings;
use strict;
package base;

# deal with bases
# check the base quality of a soft-clipped read, to see if it should be included. $min_clipped is the minimum length of soft clip sequence to check base
my $version = "0.0_10212015";
sub get_version{
    return $version;
}
sub check_basequal{
    my ($bq, $cigar, $baseThreshold, $num, $min_clipped) = @_;
    #if(defined $cigar && $cigar ne "*" && ($cigar =~ /^\d+M/ && $cigar =~ /\d+M$/)){
    #    return &checkBQ($bq, 1, length($bq), $baseThreshold, $num);
    #}
    if(!defined $cigar || $cigar eq "*" || ($cigar =~ /^\d+M/ && $cigar =~ /\d+M$/)){
        return &checkBQ($bq, 1, length($bq), $baseThreshold, $num);
    }
    if($cigar =~ /^(\d+)S.+\D(\d+)S$/){
        # if both are clipped
        if($1 > $min_clipped && $2 > $min_clipped){
            return &checkBQ($bq, 1, $1, $baseThreshold, $num) * &checkBQ($bq, 2, $2, $baseThreshold, $num);
        }
        elsif($1 > $min_clipped && $2 <= $min_clipped){
            return &checkBQ($bq, 1, $1, $baseThreshold, $num) * &checkBQ($bq, 1, length($bq), $baseThreshold, $num);
        }
        elsif($1 <= $min_clipped && $2 > $min_clipped){
            return &checkBQ($bq, 1, length($bq), $baseThreshold, $num) * &checkBQ($bq, 2, $2, $baseThreshold, $num);
        }
        else{
            return &checkBQ($bq, 1, length($bq), $baseThreshold, $num);
        }
    }
    elsif($cigar =~ /^(\d+)S/){
        if($1 > $min_clipped){
            return &checkBQ($bq, 1, $1, $baseThreshold, $num);
        }
        else{
            return &checkBQ($bq, 1, length($bq), $baseThreshold, $num);
        }
    }
    elsif($cigar =~ /(\d+)S$/){
        if($1 > $min_clipped){
            return &checkBQ($bq, 2, $1, $baseThreshold, $num);
        }
        else{
            return &checkBQ($bq, 1, length($bq), $baseThreshold, $num);
        }
    }
}

sub checkBQ{
    my ($bq, $tag, $len, $baseThreshold, $num) = @_;
    # check two ends no matter what end is soft clipped
    my $this;
    $this = substr($bq, 0, $num);
    # changed on 09292015, from mean to median (NA19240: chr1:53099221)
    #return 0 if(&mean(split(//, $this)) <= $baseThreshold);
    return 0 if(&median(split(//, $this)) <= $baseThreshold);
    $this = substr($bq, length($bq) - $num);
    #return 0 if(&mean(split(//, $this)) <= $baseThreshold);
    return 0 if(&median(split(//, $this)) <= $baseThreshold);

    # now look at the soft clipped end
    if($tag == 1){
        # clip from head
        $this = substr($bq, 0, $len);
    }
    elsif($tag == 2){
        $this = substr($bq, length($bq) - $len);
    }

    my @quals = split(//, $this);
    if(&median(@quals) > $baseThreshold){
        return 1;
    }
    return 0;
}

sub mean {
    my @vals = @_;
    my $s = 0;
    foreach my $v (@vals){
        $s += ord($v);
    }
    return $s/scalar(@vals);
}

sub median {
    my @vals = @_;
    my @vals_;
    foreach (@vals){
        push @vals_, ord($_);
    }
    @vals = sort{$a <=> $b} @vals_;
    my $len = @vals;
    if ($len % 2) {
        return $vals[int($len/2)];
    } else {
        return ($vals[$len/2] + $vals[$len/2-1])/2;
    }
}

1;
