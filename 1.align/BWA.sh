wdir=~/workspace/e.Blastocyst-K27ac/1.align/1.BWA-mESC
ddir=~/workspace/e.Blastocyst-K27ac/0.prepare/2.cutadapt
bash ~/workspace/e.Blastocyst-K27ac/1.align/thread_bwa.sh -w $wdir -d $ddir -t 8 -p 16 &

#mouse K27ac
wdir=~/workspace/e.Blastocyst-K27ac/1.align/2.BWA-mouse
ddir=~/workspace/e.Blastocyst-K27ac/0.prepare/2.cutadapt/temp
bash ~/workspace/e.Blastocyst-K27ac/1.align/thread_bwa.sh -w $wdir -d $ddir -t 8 -p 16 -k K27ac &

wdir=~/workspace/e.Blastocyst-K27ac/1.align/2.BWA-mouse
ddir=~/workspace/e.Blastocyst-K27ac/0.prepare/2.cutadapt/temp
bash ~/workspace/e.Blastocyst-K27ac/1.align/thread_bwa.sh -w $wdir -d $ddir -t 8 -p 16 -k ATAC &

wait
#human
wdir=~/workspace/e.Blastocyst-K27ac/1.align/3.BWA-human
ddir=~/workspace/e.Blastocyst-K27ac/0.prepare/2.cutadapt/H3K27ac-data/human
bash ~/workspace/e.Blastocyst-K27ac/1.align/thread_bwa.sh -w $wdir -d $ddir -x hg38 -t 8 -p 16 &



#################merge rep
ddir=~/workspace/e.Blastocyst-K27ac/1.align/2.BWA-mouse/b.filtered
wdir=~/workspace/e.Blastocyst-K27ac/1.align/2.BWA-mouse/c.merge_rep

ddir=~/workspace/e.Blastocyst-K27ac/1.align/3.BWA-human/b.filtered
wdir=~/workspace/e.Blastocyst-K27ac/1.align/3.BWA-human/c.merge_rep
mkd $wdir/logs
for j in $ddir/Morula-A485_K27ac_rep1.sorted.bam; do i=$(basename $j _rep1.sorted.bam)
(echo $i $j
samtools merge -@ 8 -f $wdir/${i}.merge.bam $ddir/${i}*.sorted.bam;
samtools flagstat -@ 8 $wdir/${i}.merge.bam > $wdir/logs/${i}.merge.flag.txt;
samtools view -@ 8 -bh -F 1024 $wdir/${i}.merge.bam -o $wdir/${i}.redun.bam;
samtools flagstat -@ 8 $wdir/${i}.redun.bam > $wdir/logs/${i}.redun.flag.txt;
samtools index -@ 8 $wdir/${i}.redun.bam) &
done;
wait


#downsample
ddir=~/workspace/e.Blastocyst-K27ac/1.align/2.BWA-mouse/c.merge_rep
wdir=~/workspace/e.Blastocyst-K27ac/1.align/2.BWA-mouse/d.downsample
mkdir -p $wdir/logs
SEED=10001; downSize=80000000 #downsample的目标数据量和种子 #K27 80M
conda deactivate
for i in $ddir/*.redun.bam; do (k=$(basename $i .redun.bam);
P=`awk -v d=$downSize '/in total/{if($1>d){p=d/$1}else{p=1.0};printf("%g",p)}' $ddir/logs/$k.redun.flag.txt`
nohup java -jar /usr/local/bin/picard.jar DownsampleSam I=$i O=$wdir/$k.downsample.bam R=$SEED P=$P A=0.00001 S=ConstantMemory > $wdir/logs/$k.log
samtools sort -@ 8 $wdir/$k.downsample.bam -o $wdir/$k.sorted.bam
samtools index -@ 8 $wdir/$k.sorted.bam
samtools flagstat -@ 8 $wdir/$k.sorted.bam > $wdir/logs/$k.flag.txt
) &
done

ddir=~/workspace/e.Blastocyst-K27ac/1.align/2.BWA-mouse/c.merge_rep
wdir=~/workspace/e.Blastocyst-K27ac/1.align/2.BWA-mouse/d.downsample
mkdir -p $wdir/logs
SEED=10001; downSize=60000000 #downsample的目标数据量和种子 #ATAC 60M
conda deactivate
for i in $ddir/*_ATAC*.redun.bam; do (k=$(basename $i .redun.bam);
P=`awk -v d=$downSize '/in total/{if($1>d){p=d/$1}else{p=1.0};printf("%g",p)}' $ddir/logs/$k.redun.flag.txt`
nohup java -jar /usr/local/bin/picard.jar DownsampleSam I=$i O=$wdir/$k.downsample.bam R=$SEED P=$P A=0.00001 S=ConstantMemory > $wdir/logs/$k.log
samtools sort -@ 8 $wdir/$k.downsample.bam -o $wdir/$k.sorted.bam
samtools index -@ 8 $wdir/$k.sorted.bam
samtools flagstat -@ 8 $wdir/$k.sorted.bam > $wdir/logs/$k.flag.txt
) &
done
#human
ddir=~/workspace/e.Blastocyst-K27ac/1.align/3.BWA-human/c.merge_rep
wdir=~/workspace/e.Blastocyst-K27ac/1.align/3.BWA-human/d.downsample
mkdir -p $wdir/logs
SEED=10001; downSize=13000000
conda deactivate
for i in $ddir/*.redun.bam; do (k=$(basename $i .redun.bam);
P=`awk -v d=$downSize '/in total/{if($1>d){p=d/$1}else{p=1.0};printf("%g",p)}' $ddir/logs/$k.redun.flag.txt`
nohup java -jar /usr/local/bin/picard.jar DownsampleSam I=$i O=$wdir/$k.downsample.bam R=$SEED P=$P A=0.00001 S=ConstantMemory > $wdir/logs/$k.log
samtools sort -@ 8 $wdir/$k.downsample.bam -o $wdir/$k.sorted.bam
samtools index -@ 8 $wdir/$k.sorted.bam
samtools flagstat -@ 8 $wdir/$k.sorted.bam > $wdir/logs/$k.flag.txt
) &
done
