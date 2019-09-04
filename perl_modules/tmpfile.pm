use warnings;
use strict;
use Time::localtime;
package tmpfile;

# generate time
sub time {
    my $tm=localtime(time);
    my $datestring = join("_", split(/\s/, localtime()));
    return $datestring;
}

1;
