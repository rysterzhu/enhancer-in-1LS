ddir=~/workspace/e.Blastocyst-K27ac/1.align/3.BWA-human/c.merge_rep
wdir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/2.bamCoverage-merge/a.zscore-human-K27ac
mkcd $wdir;mkd logs
for i in $ddir/hICM-D6*redun.bam; do k=$(basename $i .redun.bam)
(nohup bamCoverage -b $i -of bedgraph -o $k.bg --normalizeUsing RPKM -p 8 --ignoreDuplicates -bs 10 --smoothLength 30 --samFlagInclude 2 --samFlagExclude 3852 --minMappingQuality 20 > logs/$k.log
sort -S 5% --parallel=8 -k1,1 -k2n,2 $k.bg > $k.sorted.bg
awk -v bs=10 '{for(i=$2;i<=$3-bs;i+=bs){print $1,i,i+bs,$4}}' $k.sorted.bg > $k.full
nohup Rscript ~/codes/bg2scale.R $k.full > logs/${k}.scale.log
awk 'NR==1{pre1=$1;pre2=$2;pre3=$3;pre4=$4} NR>1{if($4==pre4&&$1==pre1){pre3=$3}else{print pre1,pre2,pre3,pre4;pre1=$1;pre2=$2;pre3=$3;pre4=$4}} END{print pre1,pre2,pre3,pre4}' $k.full.scale > $k.full.scale2
bedGraphToBigWig ${k}.full.scale2 ~/ann/short_chromSizes/chrom_hg38.sizes $k.bw)&
done;wait


########################mouse ATAC
#merge rep
ddir=~/workspace/e.Blastocyst-K27ac/1.align/6.xiewei-ATAC/c.merge_rep
ddir=~/workspace/e.Blastocyst-K27ac/1.align/2.BWA-mouse/c.merge_rep
wdir=~/workspace/e.Blastocyst-K27ac/2.deeptools/2.bamCoverage-merge/4.rpkm-ATAC
mkcd $wdir;mkd logs
for i in $ddir/8cell*ATAC*redun.bam; do k=$(basename $i .redun.bam)
nohup bamCoverage -b $i -o $k.bw --normalizeUsing RPKM -p 8 --ignoreDuplicates -bs 50 --smoothLength 150 --samFlagInclude 2 --samFlagExclude 3852 --minMappingQuality 20 > logs/$k.log &
done

ddir=~/workspace/e.Blastocyst-K27ac/1.align/2.BWA-mouse/c.merge_rep
wdir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/2.bamCoverage-merge/9.zscore-ATAC
mkcd $wdir;mkd logs
for i in $ddir/[IMT]*ATAC*redun.bam; do k=$(basename $i .redun.bam)
(nohup bamCoverage -b $i -of bedgraph -o $k.bg --normalizeUsing RPKM -p 8 --ignoreDuplicates -bs 10 --smoothLength 30 --samFlagInclude 2 --samFlagExclude 3852 --minMappingQuality 20 > logs/$k.log
sort -S 5% --parallel=8 -k1,1 -k2n,2 $k.bg > $k.sorted.bg
awk -v bs=10 '{for(i=$2;i<=$3-bs;i+=bs){print $1,i,i+bs,$4}}' $k.sorted.bg > $k.full
nohup Rscript ~/codes/bg2scale.R $k.full > logs/${k}.scale.log
awk 'NR==1{pre1=$1;pre2=$2;pre3=$3;pre4=$4} NR>1{if($4==pre4&&$1==pre1){pre3=$3}else{print pre1,pre2,pre3,pre4;pre1=$1;pre2=$2;pre3=$3;pre4=$4}} END{print pre1,pre2,pre3,pre4}' $k.full.scale > $k.full.scale2
bedGraphToBigWig ${k}.full.scale2 ~/ann/short_chromSizes/chrom_mm10.sizes $k.bw)&
done;wait



##############################mouse K27ac
#merge rep
ddir=~/workspace/e.Blastocyst-K27ac/1.align/2.BWA-mouse/c.merge_rep
wdir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/2.bamCoverage-merge/1.mouse-K27ac
mkcd $wdir;mkd logs
for i in $ddir/TE_K27ac*redun.bam; do k=$(basename $i .redun.bam)
nohup bamCoverage -b $i -o $k.bw --normalizeUsing RPKM -p 8 --ignoreDuplicates -bs 50 --smoothLength 150 --samFlagInclude 2 --samFlagExclude 3852 --minMappingQuality 20 > logs/$k.log &
done

#zscore #fixed or variable step, bin size not effect the result
ddir=~/workspace/e.Blastocyst-K27ac/1.align/2.BWA-mouse/c.merge_rep
wdir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/2.bamCoverage-merge/8.zscore-K27ac/a.bs50
mkcd $wdir;mkd logs
conda deactivate
for i in $ddir/*[lMaE]_K27ac.redun.bam; do k=$(basename $i .redun.bam)
(nohup bamCoverage -b $i -of bedgraph -o $k.bg --normalizeUsing RPKM -p 8 --ignoreDuplicates -bs 10 --smoothLength 30 --samFlagInclude 2 --samFlagExclude 3852 --minMappingQuality 20 > logs/$k.log
sort -S 5% --parallel=8 -k1,1 -k2n,2 $k.bg > $k.sorted.bg
awk -v bs=10 '{for(i=$2;i<=$3-bs;i+=bs){print $1,i,i+bs,$4}}' $k.sorted.bg > $k.full
nohup Rscript ~/codes/bg2scale.R $k.full > logs/${k}.scale.log
awk 'NR==1{pre1=$1;pre2=$2;pre3=$3;pre4=$4} NR>1{if($4==pre4&&$1==pre1){pre3=$3}else{print pre1,pre2,pre3,pre4;pre1=$1;pre2=$2;pre3=$3;pre4=$4}} END{print pre1,pre2,pre3,pre4}' $k.full.scale > $k.full.scale2
bedGraphToBigWig ${k}.full.scale2 ~/ann/short_chromSizes/chrom_mm10.sizes $k.bw)&
done;wait
