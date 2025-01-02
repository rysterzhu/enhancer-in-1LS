#!/bin/sh

#*****************************************************************************************
#    Filename:  thread_contact.sh
#     Creator:  Ryster Zhu
# Create Time:  2022/12/2
# Description:  build_matrix, from HiC-Pro
#     Version:  2
#*****************************************************************************************
set -e

help()
{
    cat << HELP
This is a shell script used for build_matrix and ice normalize, extract from HiC-Pro
	-w work dir, output dir in hicpro align; [require]
	-k keyword for filter input data in data dir; [default *]
	-r resolution for build matrix; [default 50000 100000]
	-c chromSize file; [default /home/qszhu/ann/hic-pro/chrom_mm10.sizes]
	-p Number of data running at the same time; [default 4]
	--hicpro hicpro script dir; [default /usr/local/software/HiC-Pro_2.11.0/scripts]

ice:
	--filter_low_counts_perc  <Percentage of reads to filter out; [default 0.02]
	--filter_high_counts_perc <Percentage of reads to filter out; [default 0]
	--remove_all_zeros_loci <If provided, all non-interacting loci will be removed prior to the filtering strategy chosen;0 or 1> [default 1]
	--max_iter <Maximum number of iterations> [default 100]
	--eps <Precision> [default 0.01]
	--output_bias <Output the bias vector;0 or 1> [default 1]
	-v 0 or 1 [default 1]
HELP
    exit 0
}

p=4
HiCPro=/usr/local/software/HiC-Pro-3.1.0/scripts
python=/usr/bin/python
chromSize=/home/qszhu/ann/hic-pro/chrom_mm10.sizes
k="*"
ress="50000 100000"
verbose=1
filter_low_counts_perc=0.02
filter_high_counts_perc=0
remove_all_zeros_loci=1
max_iter=100
eps=0.01
output_bias=1
while [ -n "$1" ]; do
case "$1" in
    -h) help; shift 1;exit 1;;
	--help) help; shift 1;exit 1;;
	-w) wdir=$2; shift 2;;
    -p) p=$2; shift 2;;
	-k) k=$2; shift 2;;
	-r) ress=$2; shift 2;;
	-c) chromSize=$2; shift 2;;
	--hicpro) HiCPro=$2; shift 2;;
	--python) python=$2; shift 2;;
	--filter_low_counts_perc) filter_low_counts_perc=$2; shift 2;;
	--filter_high_counts_perc) filter_high_counts_perc=$2; shift 2;;
#	--remove_all_zeros_loci) remove_all_zeros_loci=$2; shift 2;;
	--max_iter) max_iter=$2; shift 2;;
	--eps) eps=$2; shift 2;;
	--output_bias) output_bias=$2; shift 2;;
	-v) verbose=$2; shift 2;;
	*) echo "error: no such option $1. -h for help";exit 1;;
esac
done


mkfifo $wdir/temp.$$.fifo
exec 1000<>$wdir/temp.$$.fifo
for((i=0;i<$p;i++));do echo >&1000; done


for name in $wdir/*$k*pairs; do key=$(basename $name .pairs);for res in ${ress[@]}; do
read -u 1000
{
cat $name | $HiCPro/build_matrix --binsize $res --chrsizes $chromSize --ifile /dev/stdin --oprefix $wdir/${key}.${res} --matrix-format upper --progress 2> $wdir/logs/6.$key.$res.matrix.log && echo "$key build matrix $res done." || echo "$key build matrix $res error.";

$python $HiCPro/ice --results_filename $wdir/${key}.${res}.iced --filter_low_counts_perc $filter_low_counts_perc --filter_high_counts_perc $filter_high_counts_perc --max_iter $max_iter --eps $eps --remove-all-zeros-loci --output-bias $output_bias --verbose $verbose $wdir/${key}.${res}.matrix > $wdir/logs/6.$key.$res.ice.log && echo "$key ICE normalize $res done." || echo "$key ICE normalize $res error.";

	echo >&1000
	echo $key" done";
}&
done;done
wait
echo "All threads done."
/bin/rm -rf $wdir/temp.$$.fifo
