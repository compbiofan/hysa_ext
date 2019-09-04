if [ $# -lt 2 ]; then
    echo "Function: given the alignment of the assembled contigs to the reference, generate commands to infer the SV confirmed by short reads. ";
    echo "Usage: $0 <output_dir> <sam> <reference> <cluster_prefix> <reordered_fa> <ctg_fasta> <hmm_trans_matrix> <blasr_binary> <whether_for_complex>"; 
    exit
fi
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "1" > tmp;

output_dir=$1
sam=$2
ref=$3
cluster_prefix=$4
reordered_prefix=$5
ctg_fasta=$6
trans_matrix=$7
blasr=$8


ctg=`grep -c ">" $ctg_fasta | awk '{print $1}'`
job_num=$(($ctg/100))

if [ ! -d $output_dir ]; then
    mkdir $output_dir;
fi

if [ -e $output_dir/run_inferSV.sh ]; then
    rm $output_dir/run_inferSV.sh;
fi

if [ $# = 8 ]; then
    for i in `seq 0 $job_num`; do awk -v x=$i -v DIR=$DIR -v sam=$sam -v output_dir=$output_dir -v ref=$ref -v cluster_prefix=$cluster_prefix -v reordered_prefix=$reordered_prefix -v ctg_fasta=$ctg_fasta -v trans_matrix=$trans_matrix -v blasr=$blasr '{b=x*100+1;e=(x+1)*100;print "perl "DIR"/infer_SV.pl -n 1 -K "sam" -d "output_dir"/out."b"_"e" -P out -R "ref" -a "cluster_prefix" -b "blasr" -J "trans_matrix" -g "b":"e" -C 1 -X "reordered_prefix" "ctg_fasta" > "output_dir"/out."b"_"e".txt";}' tmp >> $output_dir/run_inferSV.sh; echo "#" >> $output_dir/run_inferSV.sh; done 
elif [ $# = 9 ]; then 
    for i in `seq 0 $job_num`; do awk -v x=$i -v DIR=$DIR -v sam=$sam -v output_dir=$output_dir -v ref=$ref -v cluster_prefix=$cluster_prefix -v reordered_prefix=$reordered_prefix -v ctg_fasta=$ctg_fasta -v trans_matrix=$trans_matrix -v blasr=$blasr '{b=x*100+1;e=(x+1)*100;print "perl "DIR"/infer_SV.pl -n 1 -K "sam" -d "output_dir"/out."b"_"e" -P out -R "ref" -a "cluster_prefix" -b "blasr" -J "trans_matrix" -g "b":"e" -C 1 -X "reordered_prefix" -D bwa -M 20 "ctg_fasta" > "output_dir"/out."b"_"e".txt";}' tmp >> $output_dir/run_inferSV.sh; echo "#" >> $output_dir/run_inferSV.sh; done 
fi
rm tmp;
