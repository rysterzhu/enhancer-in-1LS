for i in Cluster*bed; do k=$(basename $i .bed)
nohup annotatePeaks.pl $i mm10 -annStats annStats/$k.txt > annStats/$k.log &
done

for i in Cluster*bed ; do k=$(basename $i .bed);
nohup findMotifsGenome.pl $i mm10 motif/$k.motif -size given -mask -p 8 &
done


#######K27ac PCA
pdir=~/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5
wdir=~/workspace/e.Blastocyst-K27ac/2.deeptools/5.PCA/2.zscore-K27ac-reps-A485
ddir=~/workspace/e.Blastocyst-K27ac/2.deeptools/1.bamCoverage/*K27ac*
mkcd $wdir
key=rpkm-merge-peaks
multiBigwigSummary BED-file --BED $pdir/merge.peaks.bed -b $ddir/8cell_K27ac*bw $ddir/Morula_K27ac*bw $ddir/ICM_K27ac*bw $ddir/TE_K27ac_rep[9abcde]*bw $ddir/Morula-A485_K27ac*bw -out $key.npz -p 64 --outRawCounts $key.tab
plotPCA -in $key.npz -o $key.pdf --plotTitle "PCA of K27ac" --plotHeight 6 --plotWidth 8 -T "PCA of K27ac" --transpose

for i in Nfe2; do
Rscript ~/R/e.blastocyst-K27ac/1.RNA/plotExp-mouse.R $i &
done


