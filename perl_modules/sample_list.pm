use warnings;
use strict;

package sample_list;

sub get_samples {
    # given a hash of samples and a sample:location file, get the location
    my ($samples, $file) = @_;
    $DB::single = 1;
    open fh_, "<$file" or die $!;
    while(<fh_>){
        chomp;
        my ($sample, $location) = split(/:/, $_);
        if(defined $samples->{$sample}){
            $samples->{$sample} = $location;
        }
    }
    close fh_;
    return $samples;
}
1;


