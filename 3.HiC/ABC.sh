res=25000
conda-init
odir=~/workspace/e.Blastocyst-K27ac/3.HiC/5.ABC-ATAC-MACS2
for sample in 8cell ICM TE;do
wdir=$odir/${sample}
mkd $wdir/hic/raw
ln -s ~/workspace/e.Blastocyst-K27ac/3.HiC/3.ABC-analysis/${sample}/hic/raw/* $wdir/hic/raw/
done

for sample in 8cell ICM TE;do
(wdir=$odir/${sample}
src=~/git/ABC-Enhancer-Gene-Prediction-0.2.2/src
bamdir=/home/qszhu/workspace/e.Blastocyst-K27ac/1.align/2.BWA-mouse/c.merge_rep
expdir=~/workspace/e.Blastocyst-K27ac/5.RNA/1.Stringtie/4.samples
juicedir=/home/qszhu/workspace/e.Blastocyst-K27ac/3.HiC/0.matrix/juice
peakdir=/home/qszhu/workspace/e.Blastocyst-K27ac/4.callPeaks/4.MACS2-narrow-q5
chrom_sizes=~/ann/chrom_sizes/mm10.chrom.sizes
clean_gene=~/ann/Uniq.mm10/mm10.clean-gene.bed
blacklist=~/ann/mm10.blacklist.bed
juicer_jar=~/software/juicer_tools_1.22.01.jar
mkcd $wdir
rm $wdir/peaks/*
python $src/makeCandidateRegions.py \
--narrowPeak $peakdir/${sample}_ATAC_peaks.narrowPeak \
--bam $bamdir/${sample}_ATAC.redun.bam \
--outDir $wdir/peaks \
--chrom_sizes $chrom_sizes \
--peakExtendFromSummit 0 \
--nStrongestPeaks 175000 \
--ignoreSummits --minPeakWidth 500 \
--regions_blacklist $blacklist
#
rm $wdir/neighborhoods/*
python $src/run.neighborhoods.py \
--candidate_enhancer_regions $wdir/peaks/${sample}_ATAC_peaks.narrowPeak.candidateRegions.bed \
--genes $clean_gene \
--H3K27ac $bamdir/${sample}_K27ac.redun.bam \
--ATAC $bamdir/${sample}_ATAC.redun.bam \
--expression_table $expdir/${sample}-control.TPM.txt \
--tss_slop_for_class_assignment 500 \
--chrom_sizes $chrom_sizes \
--cellType ${sample} \
--outdir $wdir/neighborhoods > $wdir/neighborhoods.log

rm $wdir/Predictions/*
python $src/predict.py \
--enhancers $wdir/neighborhoods/EnhancerList.txt \
--genes $wdir/neighborhoods/GeneList.txt \
--score_column ABC.Score \
--threshold 0.02 \
--window 5000000 \
--tss_slop 500 \
--expression_cutoff 1 \
--cellType ${sample} \
--outdir $wdir/Predictions/ \
--make_all_putative \
--scale_hic_using_powerlaw \
--HiCdir $wdir/hic/raw \
--hic_resolution $res > $wdir/predict.log)&


cd 2.map2enhancer/
ddir=~/workspace/9.NT-ChIP/2.public/5.deeptools/4.bamCoverage/*/
for i in 8cell ICM TE;do
(tail -n +2 ../${i}/Predictions/EnhancerPredictionsFull.txt | sort -k12,12g | cut -f 1-3 > ${i}.sortDistance.bed
key=${i}.sortDistance.K4me3
computeMatrix scale-regions -S $ddir/${i}_K4me3.bw -R ${i}.sortDistance.bed -p 32 -a 5000 -b 5000 -m 2000 -bs 50 --skipZeros -o $key.peaks.gz --sortRegions keep
plotHeatmap -m $key.peaks.gz -out $key.peaks.pdf --dpi 360 --sortRegions keep --colorList "white,blue" --whatToShow "plot, heatmap and colorbar" -T "" --outFileSortedRegions $key.sorted.bed )&
done

for k in 8cell ICM TE; do
awk -v OFS="\t" 'FNR!=1{$1=$1;print $1,$2,$3}' ../$k/Predictions/EnhancerPredictions.txt | sort -k1,1 -k2n,2 | uniq > $k.enhancer.bed
awk -v OFS="\t" 'FNR!=1{$1=$1;if($6<1000){t=0}else{t=$6-1000};print $1,t,$6+1000}' ../$k/Predictions/EnhancerPredictions.txt | sort -k1,1 -k2n,2 | uniq > $k.promoter.bed
done
for i in 8cell ICM TE; do
multiBigwigSummary BED-file --BED $i.enhancer.bed -b $ddir/${i}_K4me3.bw -p 32 -o $i.K4me3-enhancer.npz --outRawCounts $i.K4me3-enhancer.tab &
multiBigwigSummary BED-file --BED $i.promoter.bed -b $ddir/${i}_K4me3.bw -p 32 -o $i.K4me3-promoter.npz --outRawCounts $i.K4me3-promoter.tab &
done


for i in 8cell ICM TE; do
awk 'BEGIN{print "chr\tstart\tend\tgene\tdistance\tPromoterK4me3\tEnhancerK4me3\tactivityBase"}
ARGIND==1&&FNR>1{promoter[$1"\t"$2"\t"$3]=$4} ARGIND==2{enhancer[$1"\t"$2"\t"$3]=$4}
ARGIND==3&&FNR>1{if($8<1000){t=0}else{t=$8-1000};p2=$8+1000;
print $1,$2,$3,$7,$12,promoter[$1"\t"t"\t"p2],enhancer[$1"\t"$2"\t"$3],$6}' $i.K4me3-promoter.tab $i.K4me3-enhancer.tab ../${i}/Predictions/EnhancerPredictionsFull.txt > $i.K4me3-enhancer-gene.tab &
done
