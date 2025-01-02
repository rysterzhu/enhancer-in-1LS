##############################################################################################################################
######################narrow
##############################################################################################################################
ddir=~/workspace/e.Blastocyst-K27ac/1.align/2.BWA-mouse/d.downsample
wdir=/home/qszhu/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5
mkd $wdir/logs;cd $wdir
for i in $ddir/Morula-A485_K27ac*sorted.bam; do k=$(basename $i .sorted.bam)
nohup macs2 callpeak --nomodel --nolambda -t $i --outdir $wdir -f BAMPE -g mm -n $k --keep-dup all -q 1e-5 > $wdir/logs/$k.log &
done;wait
key=K27ac
echo -e "stage\thm\tnumber\tratio\tlength" > $key.peaks.stat
for i in *[lMaE]_$key*narrowPeak; do cut -f 1,2,3 $i | uniq | awk -v k=${i%%_peaks*} '{a+=1;b+=$3-$2} END{split(k,temp,"_");print temp[1],temp[2],a,b/31e8,b/a}'  >> $key.peaks.stat; done
#distal and promoter peaks
echo -e "stage\tpromoter\ttotal\tratio" > K27ac.promoter.stat
for i in 8cell Morula ICM TE;do
intersectBed -a ${i}_K27ac_peaks.narrowPeak -b ~/ann/Uniq.mm10/all.gene/promoter_2k-2k.bed -wao |awk -v i=$i '{a+=$3-$2;b+=$NF} END{print i,b,a,b/a}' >> K27ac.promoter.stat
done

mkd annStats
for i in *.peaks.bed; do k=$(basename $i .bed)
nohup annotatePeaks.pl $i mm10 -annStats annStats/$k.stat > annStats/$k.txt 2> annStats/$k.log &
done;wait

#jaccard
bedtools jaccard -a ICM_K27ac_peaks.narrowPeak -b TE_K27ac_peaks.narrowPeak #0.35
bedtools jaccard -a ICM_K27ac_peaks.narrowPeak -b Morula_K27ac_peaks.narrowPeak #0.36
bedtools jaccard -a Morula_K27ac_peaks.narrowPeak -b TE_K27ac_peaks.narrowPeak #0.32

###################################################################
###6.K27ac-ICM-TE
awk 'FNR==1{split(FILENAME,t,"_")} {print $1,$2,$3,t[1]}' I*K27ac*narrowPeak T*K27ac*narrowPeak | sort -k1,1 -k2n,2 | mergeBed -i - -d 0 -c 4 -o distinct -delim "-" | awk '$4=="ICM-TE"{$4="both"} {print $1,$2,$3 > "6.K27ac-ICM-TE/"$4".bed"}'

#jaccard
sort -k1,1 -k2n,2 ~/ann/Uniq.mm10/all.gene/promoter_2k-2k.bed|bedtools jaccard -a ICM.bed -b - #0.008
sort -k1,1 -k2n,2 ~/ann/Uniq.mm10/all.gene/promoter_2k-2k.bed|bedtools jaccard -a TE.bed -b - #0.027
sort -k1,1 -k2n,2 ~/ann/Uniq.mm10/all.gene/promoter_2k-2k.bed|bedtools jaccard -a both.bed -b - #0.104

echo -e "flag\tchr\tratio" > promoter.ratio.tab
for i in *bed; do
intersectBed -a $i -b ~/ann/Uniq.mm10/all.gene/promoter_2k-2k.bed -u -r -f 0.5 | cut -f 1 | uniqc > temp
cut -f 1 $i | uniqc | awk -v i=${i%.*} 'NR==FNR{a[$1]=$2} NR>FNR{print i,$1,a[$1]/$2}' temp - >> promoter.ratio.tab
done

j=CGI
echo -e "flag\tchr\tratio" > $j.ratio.tab
for i in *bed; do
intersectBed -a $i -b ~/ann/mm10.$j.bed -u -r -f 0.5 | cut -f 1 | uniqc > temp
cut -f 1 $i | uniqc | awk -v i=${i%.*} 'NR==FNR{a[$1]=$2} NR>FNR{print i,$1,a[$1]/$2}' temp - >> $j.ratio.tab
done


for i in *bed; do
(awk '{print $0,1}' $i > ${i/bed/bg}
 bedGraphToBigWig ${i/bed/bg} ~/ann/short_chromSizes/chrom_mm10.sizes ${i/bed/bw})&
done
echo -e "sample\tflag\tvalue\tratio" > retain.ratio.tab && for i in *bed;do intersectBed -a $i -b ../Morula_K27ac_peaks.narrowPeak -wao | awk -v i=${i%.*} '{a+=1} ($NF/($3-$2))>0.5{b+=1} END{print i,"new",a-b,(a-b)/a;print i,"retain",b,b/a}' >> retain.ratio.tab;done

#进一步根据morula是否继承自，分类
cd ~/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5/6.K27ac-ICM-TE/morula-cluster
for i in *bed;do intersectBed -a $i -b ../Morula_K27ac_peaks.narrowPeak -wao | awk -v i=${i%.*} '{if(($NF/($3-$2))>0.5){print $1,$2,$3>"morula-cluster/"i".retain.bed"}else{print $1,$2,$3>"morula-cluster/"i".new.bed"}}' ;done
#进一步根据A485是否影响K27ac，分类
cd ~/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5/6.K27ac-ICM-TE/morula-A485-cluster
for i in *bed;do intersectBed -a $i -b ../../Morula-A485_K27ac_peaks.narrowPeak -wao | awk -v i=${i%.*} '{if(($NF/($3-$2))>0.5){print $1,$2,$3>"../morula-A485-cluster/"i".normal.bed"}else{print $1,$2,$3>"../morula-A485-cluster/"i".abnorm.bed"}}';done

cd ~/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5/6.K27ac-ICM-TE/morula-A485-cluster
mkd annStats
for i in *.bed; do k=$(basename $i .bed)
nohup annotatePeaks.pl $i mm10 -annStats annStats/$k.stat > annStats/$k.txt 2> annStats/$k.log &
done;wait
key=ann
#Chr Start End DetailedAnnotation DistancetoTSS Gene Name # $2,$3,$4,$9,$10,$16
awk -v FS="\t" 'FNR==1{gsub("annStats/","",FILENAME);gsub(".txt","",FILENAME);next} (NF>=16)&&$10<500000{print FILENAME,$16}' annStats/*txt | sort -k1,1 -k2,2 | uniq > $key.genes
cut -f 1 $key.genes | uniqc
Rscript ~/R/clusterProfiler/clusterProfiler-symbol-multiple.R $key.genes $key GO BP 12 && mv Rplots.pdf $key.pdf

mkd motifs
for i in *.bed; do k=$(basename $i .bed)
nohup findMotifsGenome.pl $i mm10 motifs/$k.motif -size given -mask -p 8 &
done
wait

#######Morula analysis A485
cd ~/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5/7.K27ac-A485
intersectBed -a ../Morula_K27ac_peaks.narrowPeak -b ../Morula-A485_K27ac_peaks.narrowPeak > both.bed
intersectBed -a ../Morula_K27ac_peaks.narrowPeak -b ../Morula-A485_K27ac_peaks.narrowPeak -v > Morula.bed
intersectBed -b ../Morula_K27ac_peaks.narrowPeak -a ../Morula-A485_K27ac_peaks.narrowPeak -v > A485.bed

mkd motifs
for i in *bed ; do k=$(basename $i .bed);
nohup findMotifsGenome.pl $i mm10 motifs/$k.motif -size given -mask -p 8 &
done

intersectBed -a ../Morula_K27ac_peaks.narrowPeak -b ../8cell_K27ac_peaks.narrowPeak -v > Morula.new.bed
intersectBed -a ../Morula_K27ac_peaks.narrowPeak -b ../8cell_K27ac_peaks.narrowPeak -u > Morula.retain.bed
intersectBed -b ../Morula_K27ac_peaks.narrowPeak -a ../8cell_K27ac_peaks.narrowPeak -v > Morula.lose.bed

###########################################################################################################
####human
ddir=~/workspace/e.Blastocyst-K27ac/1.align/3.BWA-human/d.downsample
wdir=/home/qszhu/workspace/e.Blastocyst-K27ac/4.callPeaks/b.MACS2-narrow-q5-human
mkd $wdir/logs;cd $wdir
for i in $ddir/*sorted.bam; do k=$(basename $i .sorted.bam)
nohup macs2 callpeak --nomodel --nolambda -t $i --outdir $wdir -f BAMPE -g hs -n $k --keep-dup all -q 1e-3 > $wdir/logs/$k.log &
done;wait
key=K27ac
echo -e "stage\thm\tnumber\tratio\tlength" > $key.peaks.stat
for i in *$key*narrowPeak; do cut -f 1,2,3 $i | uniq | awk -v k=${i%%_peaks*} '{a+=1;b+=$3-$2} END{split(k,temp,"-");print temp[1],temp[2],a,b/31e8,b/a}'  >> $key.peaks.stat; done


