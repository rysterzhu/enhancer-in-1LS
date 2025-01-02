ddir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/2.bamCoverage-merge/8.zscore-K27ac/
wdir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/6.map2peaks/3.K27ac-zscore/3.8cell-Morula-ICM-TE
pdir=/home/qszhu/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5
cd $wdir && mkdir -p logs
key=merge
computeMatrix reference-point -S $ddir/8cell_K27ac.bw $ddir/Morula_K27ac.bw $ddir/ICM_K27ac.bw $ddir/TE_K27ac.bw -R $pdir/K27ac.peaks.bed -p 32 -a 5000 -b 5000 --referencePoint center -bs 50 --skipZeros -o $key.peaks.gz &
wait
plotHeatmap -m $key.peaks.gz -out $key.peaks.png --dpi 100 --refPointLabel center \
--kmeans 5 --sortRegions descend --sortUsing mean --colorList "white,blue" --whatToShow "plot, heatmap and colorbar"  --samplesLabel 8cell Morula ICM TE -T "" --zMax 5
plotProfile -m $key.peaks.gz -out $key.peaks.profile.png --perGroup --regionsLabel "" --refPointLabel center -T " " --legendLocation upper-right --samplesLabel 8cell Morula ICM TE

key=merge.scale-regions
computeMatrix scale-regions -S $ddir/8cell_K27ac.bw $ddir/Morula_K27ac.bw $ddir/ICM_K27ac.bw $ddir/TE_K27ac.bw -R $pdir/K27ac.peaks.bed -p 32 -a 3000 -b 3000 -m 2000 -bs 50 --skipZeros -o $key.peaks.gz &
wait
plotHeatmap -m $key.peaks.gz -out $key.peaks.png --dpi 100 --kmeans 5 --sortRegions descend --sortUsing mean --colorList "white,blue" --whatToShow "plot, heatmap and colorbar" -T "" --zMax 3 --zMin 0 --outFileSortedRegions $key.sorted.bed
plotProfile -m $key.peaks.gz -out $key.peaks.profile.png --perGroup --regionsLabel "" -T " " --legendLocation upper-right


key=merge.kmeansInPeaks
computeMatrix scale-regions -S $ddir/8cell_K27ac.bw $ddir/Morula_K27ac.bw $ddir/ICM_K27ac.bw $ddir/TE_K27ac.bw -R $pdir/K27ac.peaks.bed -p 32 -a 0 -b 0 -m 2000 -bs 50 --skipZeros -o $key.peaks.gz && \
plotHeatmap -m $key.peaks.gz -out $key.peaks.png --dpi 100 --kmeans 6 --sortRegions descend --sortUsing mean --colorList "white,blue" --whatToShow "plot, heatmap and colorbar" -T "" --zMax 3 --zMin 0 --outFileSortedRegions $key.sorted.bed
#cluster1只有10个去掉
awk 'NR>1&&$NF!="cluster_1"{print $1,$2,$3 > $NF".bed"}' $key.sorted.bed

computeMatrix scale-regions -S $ddir/8cell_K27ac.bw $ddir/Morula_K27ac.bw $ddir/ICM_K27ac.bw $ddir/TE_K27ac.bw -R cluster*bed -p 32 -a 3000 -b 3000 -m 2000 -bs 50 --skipZeros -o $key.peaks.gz && \
plotHeatmap -m $key.peaks.gz -out $key.peaks.pdf --dpi 100 --sortRegions descend --sortUsing mean --colorList "white,blue" --whatToShow "plot, heatmap and colorbar" -T "" --zMax 3 --zMin 0 &

#Chr Start End DetailedAnnotation DistancetoTSS Gene Name # $2,$3,$4,$9,$10,$16
awk -v FS="\t" 'FNR==1{next} NR==FNR&&(NF>=16){a[$2"\t"$3"\t"$4]=$16} NR>FNR&&($1"\t"($2+1)"\t"$3 in a)&&(a[$1"\t"($2+1)"\t"$3]!=""){print $NF,a[$1"\t"($2+1)"\t"$3]}' $pdir/annStats/K27ac.peaks.txt $key.sorted.bed | grep -v "cluster_1" | sort -k1,1 -k2,2 |uniq > $key.genes
cut -f 1 $key.genes | uniqc
Rscript ~/R/clusterProfiler/clusterProfiler-symbol-multiple.R $key.genes $key GO BP 12 && mv Rplots.pdf $key.pdf

#expression
mkd motifs
for i in *.bed; do k=$(basename $i .bed)
nohup findMotifsGenome.pl $i mm10 motifs/$k.motif -size given -mask -p 8 &
done

ddir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/2.bamCoverage-merge/8.zscore-K27ac/
wdir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/6.map2peaks/3.K27ac-zscore
pdir=~/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5/6.K27ac-ICM-TE
mkdir -p logs;cd $wdir
key=ICM_TE.scale-regions
computeMatrix scale-regions -S $ddir/8cell_K27ac.bw $ddir/Morula_K27ac.bw $ddir/ICM_K27ac.bw $ddir/TE_K27ac.bw -R $pdir/ICM-TE.bed $pdir/ICM.bed $pdir/TE.bed -p 32 -a 3000 -b 3000 -m 1000 -bs 50 --skipZeros -o $key.peaks.gz &
wait
plotHeatmap -m $key.peaks.gz -out $key.peaks.png --dpi 100 --sortRegions keep --sortUsing mean --colorList "white,blue" --whatToShow "plot, heatmap and colorbar" -T "" --outFileSortedRegions $key.sorted.bed --zMax 3

key=ICM_TE.center
computeMatrix reference-point -S $ddir/8cell_K27ac.bw $ddir/Morula_K27ac.bw $ddir/ICM_K27ac.bw $ddir/TE_K27ac.bw -R $pdir/ICM-TE.bed $pdir/ICM.bed $pdir/TE.bed -p 32 -a 3000 -b 3000 -bs 50 --skipZeros -o $key.peaks.gz --referencePoint center --sortRegions no &
wait
plotHeatmap -m $key.peaks.gz -out $key.peaks.pdf --dpi 100 --sortRegions no --sortUsing mean --colorList "white,blue" --whatToShow "plot, heatmap and colorbar" -T "" --zMax 3 --outFileSortedRegions $key.sorted.bed
#distance to promoter
awk 'NR>1{sub(".bed","",$13);m=int(($2+$3)/2);print $1,m-1000,m+1000,$13,NR-1}' $key.sorted.bed | sort -k1,1 -k2n,2 > temp
awk '{m=int(($2+$3)/2);print $1,m,m+1}' ~/ann/Uniq.mm10/all.gene/promoter_2k-2k.bed | sort -k1,1 -k2n,2 | closestBed -a temp -b - -d -t first | cut -f 1-5,9 | sort -k5n,5 > $key.distance.txt
#distance to LTR
awk '{m=int(($2+$3)/2);print $1,m,m+1}' ~/ann/mm10.LTR.bed | sort -k1,1 -k2n,2 | closestBed -a temp -b - -d -t first | cut -f 1-5,9 | sort -k5n,5 > $key.LTR.distance.txt

awk '{m=int(($2+$3)/2);print $1,m,m+1,$4}' ~/ann/Uniq.mm10/all.gene/promoter_2k-2k.bed | sort -k1,1 -k2n,2 | closestBed -a temp -b - -d -t first | awk 'BEGIN{print "chr\tstart\tend\tcluster\tindex\tdistance\tICM.exp\tTE.exp\tgene"} NR==FNR{icm[$1]=$2;te[$1]=$10} NR>FNR{print $1,$2,$3,$4,$5,$10,icm[$9],te[$9],$9}' ~/workspace/e.Blastocyst-K27ac/5.RNA/1.Stringtie/all.samples.TPM.txt - | body sort -k5n,5 > $key.distance.txt # discard
#add expression foldchange
awk '{m=int(($2+$3)/2);print $1,m,m+1,$4}' ~/ann/Uniq.mm10/all.gene/promoter_2k-2k.bed | sort -k1,1 -k2n,2 | closestBed -a temp -b - -d -t first | awk 'BEGIN{print "chr\tstart\tend\tcluster\tindex\tdistance\tlog2FC\tgene"} NR==FNR{fc[$7]=$2} NR>FNR{print $1,$2,$3,$4,$5,$10,fc[$9],$9}' ~/workspace/e.Blastocyst-K27ac/5.RNA/4.DESeq2-WT/res.ICM2TE.padj005.FC2.tab - | body sort -k5n,5 > $key.distance.txt


ddir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/2.bamCoverage-merge/2.zscore-mouse-K27ac/
wdir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/6.map2peaks/3.K27ac-zscore
pdir=~/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5/5.K27ac-merge
mkdir -p logs;cd $wdir
key=all.center
computeMatrix reference-point -S $ddir/8cell_K27ac.bw $ddir/Morula_K27ac.bw $ddir/ICM_K27ac.bw $ddir/TE_K27ac.bw -R $pdir/8cell-Morula-ICM-TE.bed $pdir/8cell.bed $pdir/Morula.bed $pdir/ICM.bed $pdir/TE.bed $pdir/8cell-ICM.bed $pdir/Morula-ICM.bed $pdir/8cell-TE.bed $pdir/Morula-TE.bed $pdir/ICM-TE.bed $pdir/8cell-Morula.bed $pdir/8cell-Morula-ICM.bed $pdir/8cell-Morula-TE.bed $pdir/8cell-ICM-TE.bed $pdir/Morula-ICM-TE.bed -p 32 -a 3000 -b 3000 -bs 50 --skipZeros -o $key.peaks.gz --referencePoint center &
wait
plotHeatmap -m $key.peaks.gz -out $key.peaks.png --dpi 100 --sortRegions descend --sortUsing mean --colorList "white,blue" --whatToShow "plot, heatmap and colorbar" -T "" --zMax 3

#
key=promoter-distal
computeMatrix reference-point -S $ddir/8cell_K27ac.bw $ddir/Morula_K27ac.bw $ddir/ICM_K27ac.bw $ddir/TE_K27ac.bw -R $pdir/promoter.peaks.bed $pdir/distal.peaks.bed  -p 32 -a 3000 -b 3000 --referencePoint center -bs 50 --skipZeros -o $key.gz
plotHeatmap -m $key.gz -out $key.png --dpi 100 --refPointLabel center \
--zMax 12 --sortRegions descend --sortUsing mean --colorList "white,blue" --whatToShow "plot, heatmap and colorbar" --samplesLabel 8cell Morula ICM TE -T ""
computeMatrix scale-regions -S $ddir/8cell_K27ac.bw $ddir/Morula_K27ac.bw $ddir/ICM_K27ac.bw $ddir/TE_K27ac.bw -R ~/ann/mm10.RefSeq.bed -p 32 -a 5000 -b 5000 -m 10000 -bs 50 --skipZeros -o $key.gz
#plotHeatmap -m $key.gz -out $key.png --dpi 100 --kmeans 5 --zMax 12 --sortRegions descend --sortUsing mean --colorList "white,blue" --whatToShow "plot, heatmap and colorbar"  --samplesLabel 8cell Morula ICM TE -T "" &
plotHeatmap -m $key.gz -out $key.png --dpi 100  --sortRegions descend --sortUsing mean --colorList "white,blue" --whatToShow "plot, heatmap and colorbar" --zMax 20 --samplesLabel 8cell Morula ICM TE -T ""

#reps
ddir=~/workspace/e.Blastocyst-K27ac/2.deeptools/1.bamCoverage/mouse-K27ac
wdir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/6.map2peaks/3.K27ac-zscore
pdir=/home/qszhu/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5
mkdir -p logs
key=reps
computeMatrix scale-regions -S $ddir/8cell_K27ac*.bw $ddir/Morula_K27ac*.bw $ddir/ICM_K27ac_rep*.bw $ddir/TE_K27ac_rep[7-9abced].bw -R $pdir/blastocyst.K27ac.peaks.bed -p 32 -a 3000 -b 3000 -m 1000 -bs 10 --skipZeros -o $key.peaks.gz &
wait
plotHeatmap -m $key.peaks.gz -out $key.peaks.png --dpi 100 \
--kmeans 5 --sortRegions descend --sortUsing mean --colorList "white,blue" --whatToShow "plot, heatmap and colorbar" -T ""
plotProfile -m $key.peaks.gz -out $key.peaks.profile.png --perGroup --regionsLabel "" -T " " --legendLocation upper-right


#atac to H3K27ac peaks
ddir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/2.bamCoverage-merge/9.zscore-ATAC/
wdir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/6.map2peaks/3.ATAC-zscore
pdir=/home/qszhu/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5
mkdir -p logs;cd $wdir
key=K27acPeaks
computeMatrix reference-point -S $ddir/Morula_ATAC.bw $ddir/ICM_ATAC.bw $ddir/TE_ATAC.bw -R $pdir/K27ac.peaks.bed -p 32 -a 5000 -b 5000 --referencePoint center -bs 50 --skipZeros -o $key.peaks.gz &
wait
plotHeatmap -m $key.peaks.gz -out $key.peaks.png --dpi 100 --refPointLabel center \
--kmeans 5 --sortRegions descend --sortUsing mean --colorList "white,blue" --whatToShow "plot, heatmap and colorbar"  --samplesLabel Morula ICM TE -T ""
plotProfile -m $key.peaks.gz -out $key.peaks.profile.png --perGroup --regionsLabel "" --refPointLabel center -T " " --legendLocation upper-right --samplesLabel Morula ICM TE


#####################################################Morula A485
ddir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/2.bamCoverage-merge/*zscore*/
wdir=~/workspace/e.Blastocyst-K27ac/2.deeptools/6.map2peaks/3.K27ac-zscore/2.Morula-A485
pdir=~/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5/7.K27ac-A485
mkdir -p logs;cd $wdir
key=all.center
computeMatrix reference-point -S $ddir/Morula_K27ac.bw $ddir/Morula-A485_K27ac.bw $ddir/ICM_K27ac.bw $ddir/TE_K27ac.bw -R $pdir/both.bed $pdir/Morula.bed $pdir/A485.bed -p 32 -a 3000 -b 3000 -bs 50 --skipZeros -o $key.peaks.gz --referencePoint center &
wait
plotHeatmap -m $key.peaks.gz -out $key.peaks.png --dpi 100 --sortRegions descend --sortUsing mean --colorList "white,blue" --whatToShow "plot, heatmap and colorbar" -T "" --zMax 3


ddir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/2.bamCoverage-merge/*zscore*/
wdir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/6.map2peaks/3.K27ac-zscore/1.ICM-TE
pdir=~/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5/6.K27ac-ICM-TE
mkdir -p logs;cd $wdir
key=ICM_TE.center
computeMatrix reference-point -S $ddir/8cell_K27ac.bw $ddir/Morula_K27ac.bw $ddir/Morula-A485_K27ac.bw $ddir/ICM_K27ac.bw $ddir/TE_K27ac.bw -R $pdir/ICM-TE.bed $pdir/ICM.bed $pdir/TE.bed -p 32 -a 3000 -b 3000 -bs 50 --skipZeros -o $key.peaks.gz --referencePoint center --sortRegions no &
wait
plotHeatmap -m $key.peaks.gz -out $key.peaks.png --dpi 120 --sortRegions no --sortUsing mean --colorList "white,blue" --whatToShow "plot, heatmap and colorbar" -T "" --zMax 3 --outFileSortedRegions $key.sorted.bed

ddir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/2.bamCoverage-merge/*zscore*/
wdir=~/workspace/e.Blastocyst-K27ac/2.deeptools/6.map2peaks/3.K27ac-zscore/2.Morula-A485
pdir=~/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5/7.K27ac-A485
mkdir -p logs;cd $wdir
key=Morula.new2
computeMatrix reference-point -S $ddir/8cell_K27ac.bw $ddir/Morula_K27ac.bw $ddir/Morula-A485_K27ac.bw $ddir/ICM_K27ac.bw $ddir/TE_K27ac.bw -R $pdir/Morula.new.bed $pdir/Morula.lose.bed $pdir/Morula.retain.bed -p 32 -a 3000 -b 3000 -bs 50 --skipZeros -o $key.peaks.gz --referencePoint center &
wait
plotHeatmap -m $key.peaks.gz -out $key.peaks.png --dpi 120 --sortRegions descend --sortUsing mean --colorList "white,blue" --whatToShow "plot, heatmap and colorbar" -T ""

key=Morula.new2.cluster
computeMatrix reference-point -S $ddir/8cell_K27ac.bw $ddir/Morula_K27ac.bw $ddir/Morula-A485_K27ac.bw $ddir/ICM_K27ac.bw $ddir/TE_K27ac.bw -R $pdir/Morula.new.bed -p 32 -a 3000 -b 3000 -bs 50 --skipZeros -o $key.peaks.gz --referencePoint center &
wait
plotHeatmap -m $key.peaks.gz -out $key.peaks.png --dpi 120 --sortRegions descend --sortUsing mean --colorList "white,blue" --whatToShow "plot, heatmap and colorbar" -T "" --kmeans 3 --outFileSortedRegions $key.kmeans3.bed

key=Morula.new.cluster
computeMatrix reference-point -S $ddir/8cell_K27ac.bw $ddir/Morula_K27ac.bw $ddir/Morula-A485_K27ac.bw -R $pdir/Morula.new.bed -p 32 -a 3000 -b 3000 -bs 50 --skipZeros -o $key.peaks.gz --referencePoint center &
wait
plotHeatmap -m $key.peaks.gz -out $key.peaks.png --dpi 120 --sortRegions descend --sortUsing mean --colorList "white,blue" --whatToShow "plot, heatmap and colorbar" -T "" --kmeans 3 --outFileSortedRegions $key.kmeans3.bed
awk 'NR>1{print $1,$2,$3 > $NF".bed"}' Morula.new2.kmeans3.bed

mkd motifs-retain
for i in cluster*bed ; do k=$(basename $i .bed);
nohup findMotifsGenome.pl $i mm10 motifs-retain/$k.motif -size given -mask -p 8 &
done

key=Morula.retain.cluster
computeMatrix reference-point -S $ddir/8cell_K27ac.bw $ddir/Morula_K27ac.bw $ddir/Morula-A485_K27ac.bw -R $pdir/Morula.retain.bed -p 32 -a 3000 -b 3000 -bs 50 --skipZeros -o $key.peaks.gz --referencePoint center &
wait
plotHeatmap -m $key.peaks.gz -out $key.peaks.png --dpi 120 --sortRegions descend --sortUsing mean --colorList "white,blue" --whatToShow "plot, heatmap and colorbar" -T "" --kmeans 3 --outFileSortedRegions $key.kmeans3.bed
awk 'NR>1{print $1,$2,$3 > $NF".bed"}' Morula.retain.cluster.kmeans3.bed

#ICM TE specific peaks, and retain or new from morula
ddir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/2.bamCoverage-merge/*zscore*/
wdir=~/workspace/e.Blastocyst-K27ac/2.deeptools/6.map2peaks/3.K27ac-zscore/2.Morula-A485
pdir=~/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5/6.K27ac-ICM-TE/morula-cluster
mkdir -p logs;cd $wdir
for i in ICM TE both; do key=CWT-$i
(computeMatrix reference-point -S $ddir/Morula_K27ac.bw $ddir/Morula-A485_K27ac.bw -R $pdir/$i.retain.bed $pdir/$i.new.bed -p 32 -a 3000 -b 3000 -bs 50 --skipZeros -o $key.peaks.gz --referencePoint center
plotHeatmap -m $key.peaks.gz -out $key.peaks.png --dpi 100 --sortRegions descend --sortUsing mean --colorList "white,blue" --whatToShow "plot, heatmap and colorbar" -T "" --zMax 3)&
done
wait

(key=CWT-all;computeMatrix reference-point -S $ddir/8cell_K27ac.bw $ddir/Morula_K27ac.bw $ddir/Morula-A485_K27ac.bw $ddir/ICM_K27ac.bw $ddir/TE_K27ac.bw -R $pdir/both.retain.bed $pdir/both.new.bed $pdir/ICM.retain.bed $pdir/ICM.new.bed $pdir/TE.retain.bed $pdir/TE.new.bed -p 32 -a 3000 -b 3000 -bs 50 --skipZeros -o $key.peaks.gz --referencePoint center
plotHeatmap -m $key.peaks.gz -out $key.peaks.png --dpi 100 --sortRegions descend --sortUsing mean --sortUsingSamples 2 --colorList "white,blue" --whatToShow "plot, heatmap and colorbar" -T "" --zMax 3)&


