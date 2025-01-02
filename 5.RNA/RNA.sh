#stringtie
wdir=~/workspace/e.Blastocyst-K27ac/5.RNA/1.Stringtie
ddir=~/workspace/e.Blastocyst-K27ac/1.align/4.hisat2-mouse
mkd $wdir/1.expressions $wdir/2.gtfs $wdir/3.denovo $wdir/0.logs
for i in $ddir/8cell*.sorted.bam; do k=$(basename $i .sorted.bam);
nohup stringtie $i -p 8 -G /home/share/gtf/mm10.gtf -l $k -A $wdir/1.expressions/$k.exp -C $wdir/2.gtfs/$k.gtf -o $wdir/3.denovo/$k.gtf > $wdir/0.logs/$k.log 2>&1 &
done
wait


#DESeq2
ddir=~/workspace/e.Blastocyst-K27ac/1.align/4.hisat2-mouse
wdir=~/workspace/e.Blastocyst-K27ac/5.RNA/4.DESeq2-WT
mkcd $wdir
ln -s $ddir/*control*bam* ./
nohup featureCounts -T 64 -p -B -C -M -t exon -g gene_id -a /home/share/gtf/mm10.gtf -o all.txt *bam > all.log &
awk 'NR==FNR{a[$4]=$0} NR>FNR&&($1 in a){print a[$1]>$2".gene.bed"}' ~/ann/Uniq.mm10/mm10.uniq-gene.bed cluster-8MIT.tab

awk 'NR==FNR{a[$4]=$0} NR>FNR&&($2 in a){print a[$2]>"cluster"$1".gene.bed"}' ~/ann/Uniq.mm10/mm10.uniq-gene.bed cluster-kmeans7.tab

############morula
ddir=~/workspace/e.Blastocyst-K27ac/1.align/4.hisat2-mouse
wdir=~/workspace/e.Blastocyst-K27ac/5.RNA/3.DESeq2-Morula
mkcd $wdir
nohup featureCounts -T 64 -p -B -C -M -t exon -g gene_id -a /home/share/gtf/mm10.gtf -o Morula.txt $ddir/Morula-[AST]*bam > Morula.log &

mkcd clusterPM
awk 'NR>1&&$9!="Unchange"{print $8"_"$9,$7}' ../res.all.padj005.FC2.tab > res.all.padj005.FC2.genes
clusterPM res.all.padj005.FC2.genes res.all.padj005.FC2 GO BP 12 &


##plot exp
for i in `grep Kdm ~/ann/mm10.gene.bed| cut -f 4 `;do
Rscript ~/R/e.blastocyst-K27ac/1.RNA/plotExp-mouse.R $i &
done
#Tfap2a Tfap2c Sox2 Brd3 Brd4 Ctcf Hoxb13 Gata3 Sox17 Thap11 Brd4 Brd2 Brd1 Brd7 Brd8 Brd9 Brdt
for i in Junb Elf1; do
Rscript ~/R/e.blastocyst-K27ac/1.RNA/plotExp-mouse.R $i &
done
