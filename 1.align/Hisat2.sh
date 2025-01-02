#mouse hisat2
ddir=/home/qszhu/workspace/e.Blastocyst-K27ac/0.prepare/2.cutadapt/H3K27ac-data/mouse-RNA
ddir=~/workspace/e.Blastocyst-K27ac/0.prepare/2.cutadapt/temp
wdir=~/workspace/e.Blastocyst-K27ac/1.align/4.hisat2-mouse
mkdir -p $wdir/logs
for i in $ddir/*RNA*R1.fastq.gz; do o=`basename $i .R1.fastq.gz`
nohup hisat2 --dta --no-discordant --no-mixed -p 2 --no-unal -x /home/share/hisat2_index/new_index/mm10_tran -1 $i -2 ${i/R1./R2.} -S $wdir/$o.sam > $wdir/logs/$o.hisat2.log 2>&1 &
done
wait-m hisat2


for i in $wdir/*sam; do
(samtools view -Shb $i -o ${i/sam/bam}
samtools sort ${i/sam/bam} -o ${i/sam/sorted.bam}
samtools index ${i/sam/sorted.bam}
rm $i ${i/sam/bam}) &
done

wait