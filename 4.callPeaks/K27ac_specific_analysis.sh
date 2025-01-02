pdir=/home/qszhu/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5
wdir=/home/qszhu/workspace/e.Blastocyst-K27ac/4.callPeaks/9.K27ac-peaks-specific-analysis
mkcd $wdir;
subtractBed -a $pdir/ICM_K27ac_peaks.narrowPeak -b $pdir/8cell_K27ac_peaks.narrowPeak -A > ICM-build.bed
subtractBed -a $pdir/TE_K27ac_peaks.narrowPeak -b $pdir/8cell_K27ac_peaks.narrowPeak -A > TE-build.bed

#ICM、TE去除8cell的K27ac再去除TE、ICM的K27ac,再去除promoter，最后去除小于200bp的peak
subtractBed -a $pdir/ICM_K27ac_peaks.narrowPeak -b $pdir/8cell_K27ac_peaks.narrowPeak | subtractBed -a - -b $pdir/TE_K27ac_peaks.narrowPeak | subtractBed -a - -b ~/ann/Uniq.mm10/all.isoform/promoter_2k-2k.bed | awk '$3-$2>=200{print}' > ICM-specific.bed
subtractBed -a $pdir/TE_K27ac_peaks.narrowPeak -b $pdir/8cell_K27ac_peaks.narrowPeak | subtractBed -a - -b $pdir/ICM_K27ac_peaks.narrowPeak | subtractBed -a - -b ~/ann/Uniq.mm10/all.isoform/promoter_2k-2k.bed | awk '$3-$2>=200{print}'> TE-specific.bed
#ICM和TE相同的peak，去除8cell、promoter最后去除小于200bp的peak
intersectBed -a $pdir/ICM_K27ac_peaks.narrowPeak -b $pdir/TE_K27ac_peaks.narrowPeak | subtractBed -a - -b $pdir/8cell_K27ac_peaks.narrowPeak | subtractBed -a - -b ~/ann/Uniq.mm10/all.isoform/promoter_2k-2k.bed | awk '$3-$2>=200{print}' > ICM-TE.bed

subtractBed -a $pdir/ICM_K27ac_peaks.narrowPeak -b $pdir/TE_K27ac_peaks.narrowPeak -A | wc -l
subtractBed -a $pdir/TE_K27ac_peaks.narrowPeak -b $pdir/ICM_K27ac_peaks.narrowPeak -A | wc -l

subtractBed -a ICM-specific.bed -b ~/ann/Uniq.mm10/all.isoform/promoter_2k-2k.bed -A | wc -l

#有ATAC的peak
intersectBed -a ICM-specific.bed -b $pdir/ICM_ATAC_peaks.narrowPeak -wa -u > ICM-spe.ATAC.bed
intersectBed -a TE-specific.bed -b $pdir/TE_ATAC_peaks.narrowPeak -wa -u > TE-spe.ATAC.bed
intersectBed -a ICM-TE.bed -b $pdir/ICM_ATAC_peaks.narrowPeak -wa -u | intersectBed -a - -b $pdir/TE_ATAC_peaks.narrowPeak -wa -u > ICM-TE-spe.ATAC.bed

mkd motif-specific
for i in *specific*bed ; do k=$(basename $i .bed);
nohup findMotifsGenome.pl $i mm10 motif-specific/$k.motif -size 200 -mask -p 8 &
done
wait;

mkd motif-speATAC
for i in ICM-TE-spe.ATAC*bed ; do k=$(basename $i .bed);
nohup findMotifsGenome.pl $i mm10 motif-speATAC/$k.motif -size 200 -mask -p 8 &
done
wait;
mkd motif-build
for i in *build*bed ; do k=$(basename $i .bed);
nohup findMotifsGenome.pl $i mm10 motif-build/$k.motif -size 200 -mask -p 8 &
done
wait;

cd motif-build
awk 'BEGIN{print "Cluster\tMotif\tPvalue\tTarget"} FNR>1{split(FILENAME,temp,".");split($3,temp2,"e");gsub("%","",$7);print temp[1],$1,temp2[2],$7}' */knownResults.txt > merge.motif


mkcd ~/workspace/e.Blastocyst-K27ac/4.callPeaks/9.K27ac-peaks-specific-analysis/giggle
for i in ../ICM-TE-spe.ATAC*bed; do awk '{print $1,$2,$3}' $i | bgzip -c > $(basename $i .bed).bed.gz & done;wait
docker run -i -t -v /home:/home kubor/giggle-docker /bin/bash
cd /home/qszhu/workspace/e.Blastocyst-K27ac/4.callPeaks/9.K27ac-peaks-specific-analysis/giggle
for i in ICM-TE-spe.ATAC*bed.gz; do
giggle search -s -i /home/share/giggle_bed_index/mm10/giggle_index -q $i -g 2725537669 > ${i/bed.gz/giggle.tab} &
done;wait;
exit

for i in *giggle.tab; do
awk 'NR==1{h1=$0} NR==FNR&&FNR>1{a[$1]=$0} NR>FNR&&FNR==1{$1=h1;print $0}
NR>FNR&&FNR>1{sub("^sorted_bed/","",$1);sub(".bed.gz$","",$1);$1=a[$1];print $0}' /home/share/giggle_bed_index/mm10/metadata.tab $i > ${i/giggle.tab/result.tab} &
done

grep -i "embryonic stem cell" 8cell.result.tab |awk -v FS="\t" '{print toupper($3),"8cell",$6,$12}' > temp
grep -i "embryonic stem cell" Morula.result.tab |awk -v FS="\t" '{print toupper($3),"Morula",$6,$12}' >> temp
grep -i "embryonic stem cell" ICM-spe.ATAC.result.tab |awk -v FS="\t" '{print toupper($3),"ICM",$6,$12}' > temp
grep -i "embryonic stem cell" TE-spe.ATAC.result.tab |awk -v FS="\t" '{print toupper($3),"TE",$6,$12}' >> temp

#取file size最大的TF data
sort -k1,1 -k2,2 -k3n,3 temp | awk -v FS="\t" '{a[$1"\t"$2]=$4} END{for(i in a){print i,a[i]}}' |sort -k1,1 -k2,2 | datamash -s crosstab 1,2 mean 3 | awk -v FS="\t" 'NR==1{$0="TF"$0} $0!~"N/A"{print $0}' > all.combo_score.tab

sort -k1,1 -k2,2 -k4gr,4 temp | datamash -s crosstab 1,2 max 4 | awk 'NR==1{$0="TF"$0} $0!~"N/A"{print}' > all.combo_score.tab


pdir=/home/qszhu/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5
wdir=/home/qszhu/workspace/e.Blastocyst-K27ac/4.callPeaks/9.K27ac-peaks-specific-analysis/1.8cell-ICM-TE
mkcd $wdir;
ln -s $pdir/[8IT]*_K27*narrowPeak ./
awk 'FNR==1{split(FILENAME,t,"_")} {print $1,$2,$3,t[1]}' [8IT]*narrowPeak | sort -k1,1 -k2n,2 | mergeBed -i - -d 0 -c 4 -o distinct -delim "-" | awk '{print $1,$2,$3 > $4".bed"}'

mkd motifs
for i in *bed ; do k=$(basename $i .bed);
nohup findMotifsGenome.pl $i mm10 motifs/$k.motif -size 200 -mask -p 8 &
done
#
pdir=/home/qszhu/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5
wdir=/home/qszhu/workspace/e.Blastocyst-K27ac/4.callPeaks/9.K27ac-peaks-specific-analysis/2.8cell-Morula-ICM-TE
mkcd $wdir;
for i in 8cell Morula ICM TE; do
intersectBed -a $pdir/${i}_K27ac_peaks.narrowPeak -b $pdir/${i}_ATAC_peaks.narrowPeak | awk -v m=200 '{a=int($2/m+0.5);b=int($3/m+0.5)
for(i=a;i<b;i++){print $1,i*m,i*m+m}
}' | sort -k1,1 -k2n,2 | uniq > ${i}.fraction &
done;wait
awk 'FNR==1{split(FILENAME,t,".")} {print $1,$2,$3,t[1]}' 8cell.fraction Morula.fraction ICM.fraction TE.fraction | sort -k1,1 -k2n,2 | mergeBed -i - -d 0 -c 4 -o distinct -delim "-" | awk '{print $1,$2,$3 > $4".bed"}'
wc -l *bed

mkd motifs
for i in *bed ; do k=$(basename $i .bed);
nohup findMotifsGenome.pl $i mm10 motifs/$k.motif -size given -mask -p 8 &
done;wait;cd motifs
awk 'BEGIN{print "Cluster\tMotif\tPvalue\tTarget"} FNR>1{split(FILENAME,temp,".");split($3,temp2,"e");gsub("%","",$7);print temp[1],$1,temp2[2],$7}' */knownResults.txt > merge.motif
