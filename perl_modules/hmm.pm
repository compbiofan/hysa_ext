use warnings;
use strict;

package hmm;

# this takes the alignment result from pairHMM (with series of 0 and 1), report the cigar (always ctg at first, ref at second)
# 0: the same (match)
# 1: ctg has it but not ref
# 2: ref has it but not ctg
sub analyze_csv{
    my ($csv_file) = @_;
    my $tag = 0;
    my $num = 0;
    my $h_tag = {0 => "M",
        1 => "I",
        2 => "D"
    };
    my $previous = -1;
    my $str = "";
    open fh_, "<$csv_file" or die $!;
    while(<fh_>){
        next if($_ =~ /^#/);
        if($tag == 0){
            # the 0, 1, 2 string
            chomp;
            my @a = split(//, $_);
            foreach my $a_ (@a){
                if($a_ == $previous || $previous == -1){
                    $previous = $a_;
                    $num ++;
                }
                elsif($previous != -1){
                    $str .= $num . $h_tag->{$previous};
                    $num = 1;
                    $previous = $a_;
                }
            }
            $str .= $num . $h_tag->{$previous};
            last;
        }
    }
    close fh_;
    return $str;
}

# apart from the cigar string, also read the fourth line corresponding to the match or mismatches for every M
# match: 0; mismatch: 1; irrelevant: -
# this takes the alignment result from pairHMM (with series of 0 and 1), report the cigar (always ctg at first, ref at second)
# 0: the same (match)
# 1: ctg has it but not ref
# 2: ref has it but not ctg
sub analyze_csv_w_matches{
    my ($csv_file) = @_;
    my $tag = 0;
    my $tag1 = 0;
    my $num = 0;
    my $h_tag = {0 => "M",
        1 => "I",
        2 => "D"
    };
    my $previous = -1;
    my $str = "";
    my $line_num = 0;
    my $mismatch = 0;
    my @mismatches;
    my @a;
    my @b;
    open fh_, "<$csv_file" or die $!;
    while(<fh_>){
        next if($_ =~ /^#/);
        $line_num ++;
        if($line_num == 1){
            # the 0, 1, 2 string
            chomp;
            @a = split(//, $_);
            foreach my $a_ (@a){
                if($a_ == $previous || $previous == -1){
                    $previous = $a_;
                    $num ++;
                }
                elsif($previous != -1){
                    $str .= $num . $h_tag->{$previous};
                    $num = 1;
                    $previous = $a_;
                }
            }
            $str .= $num . $h_tag->{$previous};
        }
        elsif($line_num == 4){
            chomp;
            @b = split(//, $_);
            for(my $i = 0; $i < scalar(@b); $i ++){
                if($a[$i] == 0){
                    $tag = 1;
                    $tag1 = 0;
                    if($b[$i] == 1){
                        $mismatch ++;
                    }
                }
                elsif($tag != 0){
                    if($tag1 == 0){
                        push @mismatches, $mismatch;
                        push @mismatches, "-";
                        $mismatch = 0;
                        $tag1 = 1;
                    }
                }
                elsif($tag == 0){
                    if($tag1 == 0){
                        push @mismatches, "-";
                        $tag1 = 1;
                    }
                }
            }
            if($a[scalar(@a) - 1] == 0){
                push @mismatches, $mismatch;
            }
        }
    }
    close fh_;
    return ($str, join(":", @mismatches));
}







1;
