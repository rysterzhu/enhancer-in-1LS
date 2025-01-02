#!/bin/sh

#*****************************************************************************************
#    Filename:  thread_contact.sh
#     Creator:  Ryster Zhu
# Create Time:  2021/12/2
# Description:  thread for bam to valid pairs, from HiC-Pro
#     Version:  2
#*****************************************************************************************
set -e

help()
{
    cat << HELP
This is a shell script used for merge contact and build matrix, extract from HiC-Pro
	-w work dir, output dir in hicpro align; [require]
	-k keyword for filter input data in data dir; [default *]
	-r resolution for build matrix; [default 50000 100000]
	-c chromSize file; [default /home/qszhu/ann/hic-pro/chrom_mm10.sizes]
	-l ligation site file; [default /home/qszhu/ann/hic-pro/mm10_MboI.bed]
	-p Number of data running at the same time; [default 4]
	-t threads for sort; [default 8]
	--hicpro hicpro script dir; [default /usr/local/software/HiC-Pro_2.11.0/scripts]
	--cutoffShort short and long cis contact cutoff; [default 20000]

mapped_2hic_fragments:
	--shortestInsertSize] <Shortest insert size of mapped reads to consider> [default 100]
	--longestInsertSize] <Longest insert size of mapped reads to consider> [default 1000]
	--shortestFragmentLength] <Shortest restriction fragment length to consider> [default 50]
	--longestFragmentLength] <Longest restriction fragment length to consider> [default 100000]
	--minCisDist] <Minimum distance between intrachromosomal contact to consider> [default 0]
#	--gtag] <Genotype tag. If specified, this tag will be reported in the valid pairs output for allele specific classification> #default none
#	--all] <Write all additional output files, with information about the discarded reads (self-circle, dangling end, etc.)> #default
#	--sam] <Output an additional SAM file with flag 'CT' for pairs classification #default
#	-v 0 or 1 [default 1]

HELP
    exit 0
}

thread=8
p=4
HiCPro=/usr/local/software/HiC-Pro_2.11.1/scripts
LIGATION_FILE=/home/qszhu/ann/hic-pro/mm10_MboI.bed
python=/usr/bin/python
k="*"
shortestInsertSize=100;
longestInsertSize=1000;
shortestFragmentLength=50;
longestFragmentLength=100000;
minCisDist=50;
cutoffShort=20000;
memSize=5;
while [ -n "$1" ]; do
case "$1" in
    -h) help; shift 1;exit 1;;
	--help) help; shift 1;exit 1;;
	-w) wdir=$2; shift 2;;
    -p) p=$2; shift 2;;
    -t) thread=$2; shift 2;;
	-l) LIGATION_FILE=$2; shift 2;;
	-k) k=$2; shift 2;;
	-m) memSize=$2; shift 2;;
	--hicpro) HiCPro=$2; shift 2;;
	--python) python=$2; shift 2;;
	--shortestInsertSize) shortestInsertSize=$2; shift 2;;
	--longestInsertSize) longestInsertSize=$2; shift 2;;
	--shortestFragmentLength) shortestFragmentLength=$2; shift 2;;
	--longestFragmentLength) longestFragmentLength=$2; shift 2;;
	--minCisDist) minCisDist=$2; shift 2;;
	--cutoffShort) cutoffShort=$2; shift 2;;
	*) echo "error: no such option $1. -h for help";exit 1;;
esac
done


mkfifo $wdir/temp.$$.fifo
exec 1000<>$wdir/temp.$$.fifo
for((i=0;i<$p;i++));do echo >&1000; done


for name in $wdir/*$k*.bam; do
read -u 1000
{
	key=$(basename $name .bam)
	$python $HiCPro/mapped_2hic_fragments.py -v -t $shortestFragmentLength -m $longestFragmentLength -s $shortestInsertSize -l $longestInsertSize -d $minCisDist -a -f $LIGATION_FILE -r $wdir/$key.bam -o $wdir > $wdir/logs/5.$key.mapped_2hic_fragments.log && mv $key.RSstat $wdir/stats/ && mv $key.pairstat $wdir/stats/ && mv $key*Pairs $wdir/tmp/ && echo "$key mapped_2hic_fragments done." || echo "$key mapped_2hic_fragments error."

	LANG=en sort --parallel=$thread -S ${memSize}% -k2,2V -k3,3n -k5,5V -k6,6n $wdir/tmp/$key.validPairs > $wdir/$key.sorted.pairs && awk -F"\t" 'BEGIN{c1=0;c2=0;s1=0;s2=0} (c1!=$2||c2!=$5||s1!=$3||s2!=$6){print;c1=$2;c2=$5;s1=$3;s2=$6}' $wdir/$key.sorted.pairs > $wdir/${key}.rmdup.pairs && echo "$key sort & remove PCR duplicate done." || echo "$key sort & remove PCR duplicate error."

	(echo -e "valid_interaction\t$(cat $wdir/$key.sorted.pairs | wc -l)" > $wdir/stats/$key.validPairs.stat;
	echo -e "valid_interaction_rmdup\t$(cat $wdir/${key}.rmdup.pairs | wc -l)" >> $wdir/stats/$key.validPairs.stat;
	awk -v cutoff=$cutoffShort '$2==$5{cis=cis+1; d=$6>$3?$6-$3:$3-$6; if(d<=cutoff){sr=sr+1}else{lr=lr+1}} $2!=$5{trans=trans+1} END{print  "trans_interaction\t"trans"\ncis_interaction\t"cis"\ncis_shortRange\t"sr"\ncis_longRange\t"lr}' $wdir/${key}.rmdup.pairs >> $wdir/stats/$key.validPairs.stat;) && echo "$key validPairs stats done." || echo "$key validPairs stats error."

	echo >&1000
	echo $key" done";
}&
done
wait
echo "All threads done."
rm -rf $wdir/temp.$$.fifo
