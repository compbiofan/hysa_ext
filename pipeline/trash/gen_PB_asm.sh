perl ~/combinePBIL/scripts/make_PB_assembly.pl /scratch/1KGENOME/GIAB/NA12878/analysis3/pipe/union_find/out ./NANA 40 2400 /scratch/1KGENOME/GIAB/NA12878/analysis3 /scratch/1KGENOME/GIAB/NA12878/analysis3/pipe/assembly NA NA > run_PB_asm.sh
pb_prefix=$1
tmp_dir=$2
cluster=$3
overall_f_prefix=$4
group=$5
time_limit=$6
