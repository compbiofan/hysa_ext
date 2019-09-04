use warnings;
use strict;
package seq;

my %h = (
    "A" => "T",
    "C" => "G",
    "G" => "C",
    "T" => "A",
    "N" => "N",
    "a" => "t",
    "c" => "g",
    "g" => "c",
    "t" => "a",
    "n" => "n",
);

# deal with the basic sequence
sub revcom{
    # reverse complement
    my ($s) = @_;
    my $r = "";
    my @a = split(//, $s);
    foreach (@a){
        if(!defined $h{$_}){
            print "Error: in sequence $s, base $_ is not in our database. Will use N instead.\n";
            $_ = "N"; 
        }
        $r .= $h{$_};
    }
    $r = scalar reverse $r;
    return $r;
}

# deal with the fasta file, reverse com every read in the fa file and print to stdout. 
sub revcom_file{
    # reverse complement
    my ($fa) = @_;
    my $s = "";
    my $s_a = "";
    open fa_fh, "<$fa" or die $!;
    while(<fa_fh>){
        if($_ =~ /^>/){
            # changed 08282016
            $s_a .= &print_by_len(&revcom($s), 60) if($s ne "");
            $s_a .= $_;
            $s = "";
        }
        else{
            chomp;
            $s .= $_;
        }
    }
    $s_a .= &print_by_len(&revcom($s), 60) if($s ne "");
    close fa_fh;
    open fa_fh, ">$fa" or die $!;
    print fa_fh $s_a;
    close fa_fh;
}

# print the string with each line at most max_len
sub print_by_len{
    my ($s, $max_len) = @_;
    my $len = length($s);
    my $line = int($len/$max_len);
    my $start = 0;
    my $s_a = "";
    foreach my $i (0 .. $line - 1){
        # changed 08282016
        $s_a .= substr($s, $start, $max_len) . "\n";
        #print substr($s, $start, $max_len) . "\n";
        $start += $max_len;
    }
    # changed 08282016
    $s_a .= substr($s, $start, ) . "\n" if($len % $max_len != 0);
    return $s_a;
    #print substr($s, $start, ) . "\n" if($len % $max_len != 0);
}
