use warnings;
use strict;

if(@ARGV == 0){
    die "This combines one-end mapped and both unmapped sam files into one sam file. Change the index of unmapped so that it will be appended to one-end mapped file. Generate a new index file; cp the split and drp index into this new directory. \nUsage: $0 <one_mapped_dir> <one_mapped_prefix> <both_unmapped_dir> <both_unmapped_prefix> <new_dir> <new_file_prefix>\n";
}

my ($one_mapped_dir, $one_mapped_prefix, $both_unmapped_dir, $both_unmapped_prefix, $new_dir, $prefix) = @ARGV;
my $new_sam = "$new_dir/$prefix.sam";
my $new_index = "$new_dir/$prefix.id";
`cp $one_mapped_dir/discordant.sam $new_sam`;
`cp $one_mapped_dir/index.txt $new_index`;
my @tmp = split(/\t/, `tail -n 1 $new_index`);
my $last_line_mapped = $tmp[1];
my $id_unmapped = $last_line_mapped;
my $tag = 1;
open index_fh, ">>$new_index" or die $!;
open sam_fh, ">>$new_sam" or die $!;
open fh_unmapped, "<$both_unmapped_dir/$both_unmapped_prefix.filtered.sam" or die $!;
while(<fh_unmapped>){
    chomp;
    my @a = split(/\t/, $_);
    if($tag == 0){
        $tag = 1;
    }
    elsif($tag == 1){
        $id_unmapped ++;
        $tag = 0;
        my $index_line = join("\t", $a[0], $id_unmapped) . "\n";
        print index_fh $index_line;
    }
    my $line = join("\t", $id_unmapped, @a[1 .. $#a]) . "\n";
    print sam_fh $line;
}
close fh_unmapped;
close sam_fh;
close index_fh;

# prepare fa for paired reads
`~/bin/sam2fa_pair.sh $new_dir/$prefix`;



