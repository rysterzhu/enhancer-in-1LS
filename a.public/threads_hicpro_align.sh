#!/bin/sh

#*****************************************************************************************
#    Filename:  thread_align.sh
#     Creator:  Ryster Zhu
# Create Time:  2021/03/27
# Description:  thread for alignment Hi-C to bam, from HiC-Pro
#     Version:  1
#*****************************************************************************************
set -e

help()
{
    cat << HELP
This is a shell script used for alignment Hi-C and merge, extract from HiC-Pro
	-d data dir; [require]
	-o output dir; [require]
	-k keyword for filter input data in data dir; [default *]
	-p Number of data running at the same time; [default 4]
	-t threads for bowtie2; [default 8]
	-g genome; [default mm10]
	-q min mapping quality; [default 0]
	-l ligation site sequence; [default GATC]
	--hicpro hicpro script dir; [default /usr/local/software/HiC-Pro_2.11.1/scripts]
#	--bowtie2indext bowtie2 index dir; [default /home/share/bowtie2_index]
#default other bowtie2 arguement:
#	global: --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder --rg-id BMG --phred33-quals
#	local: --very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder --rg-id BML --phred33-quals
HELP
    exit 0
}

thread=8
p=4
HiCPro=/usr/local/software/HiC-Pro_2.11.1/scripts
python=/usr/bin/python
minQual=0
LIGATION_SITE="GATCGATC"
k=*

while [ -n "$1" ]; do
case "$1" in
    -h) help; shift 1;exit 1;;
	--help) help; shift 1;exit 1;;
	-o) odir=$2; shift 2;;
    -d) ddir=$2; shift 2;;
    -p) p=$2; shift 2;;
	-g) g=$2; shift 2;;
	-q) minQual=$2; shift 2;;
    -t) thread=$2; shift 2;;
	-l) LIGATION_SITE=$2; shift 2;;
	-k) k=$2; shift 2;;
	--hicpro) HiCPro=$2; shift 2;;
	--python) python=$2; shift 2;;
    -*) echo "error: no such option $1. -h for help";exit 1;;
     *) break;
esac
done

bwt2index=/home/share/bowtie2_index/$g


mkdir -p $odir/logs $odir/tmp $odir/stats
mkfifo $odir/temp.$$.fifo
exec 1000<>$odir/temp.$$.fifo
for((i=0;i<$p;i++));do echo >&1000; done


for name in $ddir/*$k*.R1.fastq.gz; do
read -u 1000
{
	key=$(basename $name .R1.fastq.gz)
	for i in 1 2; do
		bowtie2 --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder --rg-id BMG --phred33-quals -p $thread --rg SM:$key.R${i} --un $odir/tmp/$key.bwt2glob.unmap.R${i}.fastq -x $bwt2index -U $ddir/$key.R${i}.fastq.gz 2> $odir/logs/1.$key.R${i}.bwt2glob.log | samtools view -F 4 -bS - > $odir/tmp/$key.R${i}.bwt2glob.bam && echo "$key bwtglob done." || echo "$key bwtglob error."

		$HiCPro/cutsite_trimming --fastq $odir/tmp/$key.bwt2glob.unmap.R${i}.fastq --cutsite ${LIGATION_SITE} --out $odir/tmp/$key.bwt2glob.unmap_trimmed.R${i}.fastq > $odir/logs/2.$key.R${i}.trim.global.log 2>&1 && echo "$key cutsite done." || echo "$key cutsite error."

		bowtie2 --very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder --rg-id BML --phred33-quals -p $thread --rg SM:$key.R${i}.bwt2glob.unmap -x $bwt2index -U $odir/tmp/$key.bwt2glob.unmap_trimmed.R${i}.fastq 2>$odir/logs/3.$key.R${i}.bwt2loc.log | samtools view -bS - > $odir/tmp/$key.R${i}.bwt2loc.bam && echo "$key bwt2loc done." || echo "$key bwt2loc error."

		samtools merge -@ $thread -n -f $odir/tmp/$key.R${i}.bwt2merged.bam $odir/tmp/$key.R${i}.bwt2glob.bam $odir/tmp/$key.R${i}.bwt2loc.bam && echo "$key bwt2merge done." || echo "$key bwt2merge error."

		samtools sort -@ $thread -n $odir/tmp/$key.R${i}.bwt2merged.bam -T $odir/tmp/$key.R${i} -o $odir/tmp/$key.R${i}.bwt2merged.sorted.bam && mv $odir/tmp/$key.R${i}.bwt2merged.sorted.bam $odir/tmp/$key.R${i}.bwt2merged.bam && echo "$key bwt2sort done." || echo "$key bwt2sort error."

		(echo -e "total\t$(samtools view -@ $thread -c $odir/tmp/$key.R${i}.bwt2merged.bam)" > $odir/stats/$key.R${i}.bwt.stat;
		echo -e "mapped\t$(samtools view -@ $thread -c -F 4 $odir/tmp/$key.R${i}.bwt2merged.bam)" >> $odir/stats/$key.R${i}.bwt.stat;
		echo -e "global\t$(samtools view -@ $thread -c -F 4 $odir/tmp/$key.R${i}.bwt2glob.bam)" >> $odir/stats/$key.R${i}.bwt.stat;
		echo -e "local\t$(samtools view -@ $thread -c -F 4 $odir/tmp/$key.R${i}.bwt2loc.bam)" >> $odir/stats/$key.R${i}.bwt.stat;) && echo "$key stats done." || echo "$key stats error."
	done

	$python $HiCPro/mergeSAM.py -q $minQual -t -v -f $odir/tmp/$key.R1.bwt2merged.bam -r $odir/tmp/$key.R2.bwt2merged.bam -o $odir/$key.bam > $odir/logs/4.$key.mergeSAM.log && echo "$key mergeSAM done." || echo "$key mergeSAM error."

	echo >&1000
	echo $key" done"
}&
done

wait
echo "All threads done."

/bin/rm -rf $odir/temp.$$.fifo


