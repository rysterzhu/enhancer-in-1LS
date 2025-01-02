wdir=~/workspace/e.Blastocyst-K27ac/2.deeptools/8.multiBigwigSummary/1.K27ac_zscore
ddir=~/workspace/e.Blastocyst-K27ac/2.deeptools/2.bamCoverage-merge/8.zscore-K27ac
pdir=~/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5
mkcd $wdir
key=fraction
multiBigwigSummary BED-file --BED $pdir/merge.K27ac.fraction -b $ddir/8cell_K27ac.bw $ddir/Morula_K27ac.bw $ddir/ICM_K27ac.bw $ddir/TE_K27ac.bw -o $wdir/$key.npz --outRawCounts $wdir/$key.tab --smartLabels --numberOfProcessors 16

sed -i "1s/[#']//g" $key.tab
