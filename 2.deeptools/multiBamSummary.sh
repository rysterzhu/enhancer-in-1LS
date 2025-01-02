ddir=/home/qszhu/workspace/e.Blastocyst-K27ac/1.align/2.BWA-mouse/c.merge_rep
wdir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/4.multiBamSummary/1.K27ac-merge_rep/
mkd $wdir/logs;cd $wdir
conditions=(8cell Morula ICM TE)
for i in ${conditions[@]}; do for k in promoter_1k-1k promoter_10k-10k promoter_2k-2k; do key=$i.$k
nohup multiBamSummary BED-file --BED ~/ann/makewindows/mm10.$k.bed -b $ddir/${i}_K27ac.redun.bam -o $wdir/logs/$key.gz --smartLabels -p 8 --centerReads --outRawCounts $wdir/$key.tab --ignoreDuplicates --samFlagInclude 2 --samFlagExclude 1024 > $wdir/logs/$key.log &
done;done;wait-m;for i in $wdir/*chrom.tab; do sed -i "1s/['#]//g" $i; sed -i "1s/.sorted//g" $i; sed -i "1s/.redun//g" $i; done

for i in *tab;do c=${i##*/};c=${c%%.*}
awk 'NR==1{print;next} NR==FNR{for(i=4;i<=NF;i++)a[i]+=$i} NR>FNR&&FNR>1{for(i=4;i<=NF;i++)$i=($i*1e9)/(($3-$2)*a[i]);print $0}' $wdir/$c.chrom.tab $i > ${i/tab/rpkm} &
done;

