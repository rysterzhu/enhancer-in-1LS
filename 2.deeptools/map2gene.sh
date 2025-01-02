ddir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/2.bamCoverage-merge/2.zscore-mouse-K27ac/
wdir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/7.map2gene/3.K27ac-zscore
pdir=/home/qszhu/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5
mkdir -p $wdir/logs
key=map2gene
computeMatrix scale-regions -S $ddir/8cell_K27ac.bw $ddir/Morula_K27ac.bw $ddir/ICM_K27ac.bw $ddir/TE_K27ac.bw -R ~/ann/mm10.RefSeq.bed -p 32 -a 5000 -b 5000 -m 10000 -bs 50 --skipZeros -o $wdir/$key.gz
#plotHeatmap -m $wdir/$key.gz -out $wdir/$key.png --dpi 360 --kmeans 5 --zMax 12 --sortRegions descend --sortUsing mean --colorList "white,blue" --whatToShow "plot, heatmap and colorbar"  --samplesLabel 8cell Morula ICM TE -T "" &
plotHeatmap -m $wdir/$key.gz -out $wdir/$key.png --dpi 360  --sortRegions descend --sortUsing mean --colorList "blue,black,yellow" --whatToShow "plot, heatmap and colorbar" --zMax 3 --samplesLabel 8cell Morula ICM TE -T ""


############ATAC
wdir=~/workspace/e.Blastocyst-K27ac/2.deeptools/7.map2gene/1.ATAC
ddir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/2.bamCoverage-merge/mouse-K27ac/
key=map2gene
computeMatrix scale-regions -S $ddir/Morula_ATAC.bw $ddir/ICM_ATAC.bw $ddir/TE_ATAC.bw -R ~/ann/mm10.RefSeq.bed -p 32 -a 5000 -b 5000 -m 10000 -bs 50 --skipZeros -o $wdir/$key.gz
#plotHeatmap -m $wdir/$key.gz -out $wdir/$key.png --dpi 360 --kmeans 5 --zMax 12 --sortRegions descend --sortUsing mean --colorList "white,blue" --whatToShow "plot, heatmap and colorbar"  --samplesLabel 8cell Morula ICM TE -T "" &
plotHeatmap -m $wdir/$key.gz -out $wdir/$key.png --dpi 360  --sortRegions descend --sortUsing mean --colorList "white,blue" --whatToShow "plot, heatmap and colorbar" --zMax 20 --samplesLabel Morula ICM TE -T "ATAC"
