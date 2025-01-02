ddir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/2.bamCoverage-merge/8.zscore-K27ac/
wdir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/9.map2TAD/1.zscore-K27ac
pdir=~/workspace/e.Blastocyst-K27ac/3.HiC/8.IS/3.cat_tad
mkdir -p $wdir/logs && cd $wdir
for key in ICM TE;do
(computeMatrix scale-regions -S $ddir/8cell_K27ac.bw $ddir/Morula_K27ac.bw $ddir/ICM_K27ac.bw $ddir/TE_K27ac.bw -R $pdir/$key.25000.filter.bed -p 32 -a 25000 -b 25000 --unscaled5prime 25000 --unscaled3prime 25000 -m 100000 -bs 100 --skipZeros -o $key.25k.gz &&\
plotHeatmap -m $key.25k.gz -out $key.25k.png --dpi 100 --sortRegions descend --sortUsing mean --colorList "white,blue" --whatToShow "plot, heatmap and colorbar" -T "" --zMax 3 --zMin 0 &
plotProfile -m $key.25k.gz -out $key.25k.profile.png --perGroup --regionsLabel "" -T " " --legendLocation upper-right) &
done;wait

ddir=~/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5/6.K27ac-ICM-TE
wdir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/9.map2TAD/2.peak-bw
pdir=~/workspace/e.Blastocyst-K27ac/3.HiC/8.IS/3.cat_tad
mkdir -p $wdir/logs && cd $wdir
for key in ICM TE;do
(computeMatrix scale-regions -S $ddir/both.bw $ddir/ICM.bw $ddir/TE.bw -R $pdir/$key.25000.filter.bed -p 32 -a 25000 -b 25000 --unscaled5prime 25000 --unscaled3prime 25000 -m 100000 -bs 1000 --skipZeros -o $key.25k.gz &
wait
plotHeatmap -m $key.25k.gz -out $key.25k.png --dpi 100 --sortRegions descend --sortUsing mean --colorList "white,blue" --whatToShow "plot, heatmap and colorbar" -T "" #--zMax 3 --zMin 0
plotProfile -m $key.25k.gz -out $key.25k.profile.png --perGroup --regionsLabel "" -T " " --legendLocation upper-right) &
done

