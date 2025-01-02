###################ATAC
ddir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/1.bamCoverage/mouse-ATAC
wdir=~/workspace/e.Blastocyst-K27ac/2.deeptools/3.plotCorrelation/3.mouse-ATAC
mkcd $wdir
#key=all #bins
key=xiewei-all
multiBigwigSummary BED-file --BED ~/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5/ATAC.peaks.bed -b $ddir/*bw -o $wdir/${key}.npz --outRawCounts $wdir/${key}.tab -p 32 --smartLabels &
wait
plotPCA -in $wdir/${key}.npz -o $wdir/${key}-transpose.png --transpose --plotWidth 10 --outFileNameData $wdir/${key}-transpose.tab &
plotCorrelation -in $wdir/${key}.npz -o $wdir/${key}.pearson.heatmap.pdf -c pearson -p heatmap --skipZeros --removeOutliers --plotNumbers -T "Correlation of ATAC" &

#####################in peaks correlation
ddir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/1.bamCoverage/1.mouse-K27ac
wdir=~/workspace/e.Blastocyst-K27ac/2.deeptools/3.plotCorrelation/5.mouse-K27ac-peaks
mkcd $wdir
key=20240410
multiBigwigSummary BED-file --BED ~/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5/merge.peaks.bed -b $ddir/!(*sc*)K27ac*bw -o $wdir/${key}.npz --outRawCounts $wdir/${key}.tab -p 32 --smartLabels &
wait
plotPCA -in $wdir/${key}.npz -o $wdir/${key}-transpose.png --transpose --plotWidth 10 --outFileNameData $wdir/${key}-transpose.tab &
plotPCA -in $wdir/${key}.npz -o $wdir/${key}-rowCenter.png --rowCenter --plotWidth 10 --outFileNameData $wdir/${key}-rowCenter.tab &
plotCorrelation -in $wdir/${key}.npz -o $wdir/${key}.pearson.heatmap.pdf -c pearson -p heatmap --skipZeros --removeOutliers --plotNumbers -T "Correlation of K27ac" &

###reps peaks
ddir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/1.bamCoverage/mouse-K27ac
wdir=~/workspace/e.Blastocyst-K27ac/2.deeptools/3.plotCorrelation/6.mouse-K27ac-peaks-reps
mkcd $wdir
key=WT4
cat ~/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5/4.K27ac-reps/*[laM]_K27ac_rep*narrowPeak ~/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5/4.K27ac-reps/TE_K27ac_rep[78ced]*narrowPeak | sort -k1,1 -k2n,2 | mergeBed -i - -d 1000 > merge.$key.peaks.bed
multiBigwigSummary BED-file --BED merge.$key.peaks.bed -b $ddir/*[laM]_K27ac_rep*bw $ddir/TE_K27ac_rep[78ced]*bw -o ${key}.npz --outRawCounts ${key}.tab -p 32 --smartLabels &
wait
plotPCA -in $wdir/${key}.npz -o $wdir/${key}-transpose.png --transpose --plotWidth 10 --outFileNameData $wdir/${key}-transpose.tab &
plotCorrelation -in $wdir/${key}.npz -o $wdir/${key}.pearson.heatmap.pdf -c pearson -p heatmap --skipZeros --removeOutliers --plotNumbers -T "Correlation of K27ac" &

ddir=/home/qszhu/workspace/e.Blastocyst-K27ac/2.deeptools/1.bamCoverage/mouse-ATAC
pdir=~/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5/4.K27ac-reps
wdir=~/workspace/e.Blastocyst-K27ac/2.deeptools/3.plotCorrelation/3.mouse-ATAC/
mkcd $wdir
key=mergePeaks
cat $pdir/[ITM]*ATAC*Peak $pdir/8cell_ATAC_rep[2-5]*Peak | sort -k1,1 -k2n,2 | mergeBed -i - -d 0 > merge.ATAC.peaks.bed
multiBigwigSummary BED-file --BED merge.ATAC.peaks.bed -b $ddir/8cell_ATAC_rep[2-5].bw $ddir/[ITM]*ATAC_rep*bw -o ${key}.npz --outRawCounts ${key}.tab -p 32 --smartLabels &
wait
plotPCA -in $wdir/${key}.npz -o $wdir/${key}-transpose.png --transpose --plotWidth 10 --outFileNameData $wdir/${key}-transpose.tab &
plotCorrelation -in $wdir/${key}.npz -o $wdir/${key}.pearson.heatmap.pdf -c pearson -p heatmap --skipZeros --removeOutliers --plotNumbers -T "Correlation of ATAC" &
