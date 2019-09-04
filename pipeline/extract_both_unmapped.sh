if [ $# -lt 2 ]; then
    echo "Function: extract both unmapped reads from a bam, generating a sam called unmapped.both.sam, then filter out low quality reads. "
    echo "Usage: $0 <bam> <output_dir> <output_prefix>"
    exit
fi
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
bam=$1
dir=$2
prefix=$3
samtools view -f 12 $bam > $dir/$prefix.sam
perl $DIR/filterPairedUnmapped.pl $dir/$prefix.sam $prefix.filtered > run.out
