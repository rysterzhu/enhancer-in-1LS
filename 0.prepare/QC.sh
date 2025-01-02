wdir=~/workspace/e.Blastocyst-K27ac/0.prepare/0.data/20230908-human-K27ac
mkdir -p $wdir ${wdir/0.data/1.qc} ${wdir/0.data/2.cutadapt} ${wdir/0.data/3.qc_cutadapt}
nohup fastqc -o ${wdir/0.data/1.qc} $wdir/*gz -t 64 &
for i in $wdir/*R1.fastq.gz; do o=${i/0.data/2.cutadapt}
nohup trim_galore -j 8 --length 20 --paired -q 20 --trim-n -o ${wdir/0.data/2.cutadapt} $i ${i/R1/R2} > ${o/R1.fastq.gz/trim.log} &
done
wait-m trim
for i in ${wdir/0.data/2.cutadapt}/*.gz; do
mv $i ${i/_val_[12].fq.gz/.fastq.gz}
done
wait
nohup fastqc -o ${wdir/0.data/3.qc_cutadapt} ${wdir/0.data/2.cutadapt}/*.gz -t 64 &
wait
cd ${wdir/0.data/1.qc}
multiqc --no-data-dir -f -dd 10 ./ -n $wdir/raw_H3K27ac.qc.html
cd ${wdir/0.data/2.cutadapt}
multiqc --no-data-dir -f -dd 10 ./ -n $wdir/cutadapt_H3K27ac.html
cd ${wdir/0.data/3.qc_cutadapt}
multiqc --no-data-dir -f -dd 10 ./ -n $wdir/clean_H3K27ac.qc.html
wait-m

