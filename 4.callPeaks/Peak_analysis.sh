###########################K27ac peaks analysis giggle
pdir=/home/qszhu/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5
wdir=/home/qszhu/workspace/e.Blastocyst-K27ac/4.callPeaks/7.K27ac-peaks-analysis/1.K27ac-motifs
mkcd $wdir
for i in $pdir/*_K27ac_peaks.narrowPeak ; do k=$(basename $i _K27ac_peaks.narrowPeak);
nohup findMotifsGenome.pl $i mm10 $k.motif -size given -mask -p 8 &
done

pdir=/home/qszhu/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5
wdir=/home/qszhu/workspace/e.Blastocyst-K27ac/4.callPeaks/7.K27ac-peaks-analysis/2.K27ac-ATAC-motifs
mkcd $wdir
for i in 8cell Morula ICM TE; do
intersectBed -a $pdir/${i}_K27ac_peaks.narrowPeak -b $pdir/${i}_ATAC_peaks.narrowPeak |awk '$3-$2>200' > $wdir/$i.bed & done; wait
wc -l *
for i in *bed ; do k=$(basename $i.bed);
nohup findMotifsGenome.pl $i mm10 $k.motif -size 200 -mask -p 8 &
done

cd ~/workspace/e.Blastocyst-K27ac/4.callPeaks/7.K27ac-peaks-analysis/3.giggle

for i in /home/qszhu/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5/*_K27ac_peaks.narrowPeak; do awk '{print $1,$2,$3}' $i | bgzip -c > $(basename $i _K27ac_peaks.narrowPeak).bed.gz & done;wait

docker run -i -t -v /home:/home kubor/giggle-docker /bin/bash
cd /home/qszhu/workspace/e.Blastocyst-K27ac/4.callPeaks/7.K27ac-peaks-analysis/3.giggle
for i in *bed.gz; do
giggle search -s -i /home/share/giggle_bed_index/mm10/giggle_index -q $i -g 2725537669 > ${i/bed.gz/giggle.tab} &
done;wait;
exit

for i in *giggle.tab; do
awk 'NR==1{h1=$0} NR==FNR&&FNR>1{a[$1]=$0} NR>FNR&&FNR==1{$1=h1;print $0}
NR>FNR&&FNR>1{sub("^sorted_bed/","",$1);sub(".bed.gz$","",$1);$1=a[$1];print $0}' /home/share/giggle_bed_index/mm10/metadata.tab $i > ${i/giggle.tab/result.tab} &
done

grep -i "embryonic stem cell" 8cell.result.tab |awk -v FS="\t" '{print toupper($3),"8cell",$6,$12}' > temp
grep -i "embryonic stem cell" Morula.result.tab |awk -v FS="\t" '{print toupper($3),"Morula",$6,$12}' >> temp
grep -i "embryonic stem cell" ICM.result.tab |awk -v FS="\t" '{print toupper($3),"ICM",$6,$12}' >> temp
grep -i "embryonic stem cell" TE.result.tab |awk -v FS="\t" '{print toupper($3),"TE",$6,$12}' >> temp

#取file size最大的TF data
sort -k1,1 -k2,2 -k3n,3 temp | awk -v FS="\t" '{a[$1"\t"$2]=$4} END{for(i in a){print i,a[i]}}' |sort -k1,1 -k2,2 | datamash -s crosstab 1,2 mean 3 | awk -v FS="\t" 'NR==1{$0="TF"$0} $0!~"N/A"{print $1,$2,$4,$3,$5}' > all.combo_score.tab

sort -k1,1 -k2,2 -k4gr,4 temp | datamash -s crosstab 1,2 max 4 | awk 'NR==1{$0="TF"$0} $0!~"N/A"{print}' > all.combo_score.tab




