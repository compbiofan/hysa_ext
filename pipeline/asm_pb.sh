if [ $# -lt 6 ]; then
    echo "Function: This script runs the assembly for the long reads in a cluster. \nUsage: $0 <pb_prefix> <tmp_dir> <cluster> <assembled_ctg_prefix> <cluster_index> <group> <time_limit> <thread>\n"; 
    exit
fi
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# this is exactly the same as that in ~/d/, except the directory is now CHM1.
#source ~/.profile
pb_prefix=$1
tmp_dir=$2
cluster=$3
overall_f_prefix=$4
group=$5
time_limit=$6
thread=$7
if [ ! -d $tmp_dir ]; then
    mkdir -p $tmp_dir
fi
pushd $tmp_dir
if [ -d asm_g$group ]; then
    rm -r asm_g$group
fi
# generate the group file
echo "perl $DIR/generate_group_file.pl $cluster $group ./out.g$group.txt"
perl $DIR/generate_group_file.pl $cluster $group ./out.g$group.txt
# generate the reads
perl $DIR/make_group_fa_dir.pl ./out.g$group.txt NA $pb_prefix ref.fa PB_only
# make fastq file
perl $DIR/fa2fq.pl out.g$group.txt.PB.fa > out.g$group.txt.PB.fq
# generate frg file
fastqToCA -libraryname NA -technology pacbio-raw -reads out.g$group.txt.PB.fq > out.g$group.txt.PB.frg
# assembly
timeout $time_limit $DIR/AssembleCA_wFrg.sh asm_g$group $thread out.g$group.txt.PB.frg 
# infer breakpoints
fasta_f=asm_g$group/asm_g$group.ctg.fasta
if [ -e $fasta_f ] && [ -s $fasta_f ]; then
    # change readname
    perl $DIR/change_readname.pl $fasta_f g$group head 

    # combine both SV inference and INDEL inference to output a breakdancer formatted file
    #perl ~/u/infer_breakpoints.pl $fasta_f > inferred.sv
    #perl ~/u/check_breakpoint.pl inferred.sv out.g$2.txt.IL out.g$2.txt.PB $fasta_f > checked.sv
    # clean up
    cat $fasta_f >> $overall_f_prefix.fa 
    #cat inferred.sv >> $overall_f_prefix.inferred.sv
    #cat checked.sv >> $overall_f_prefix.checked.sv
else
    echo "Cluster $2 does not assemble successfully within $time_limit seconds." >> $overall_f_prefix.err
fi
popd
rm -r $tmp_dir
echo "Cluster $2 has been finished. Result in $overall_f_prefix.[fa|err]";
