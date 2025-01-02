ddir=~/workspace/e.Blastocyst-K27ac/3.HiC/0.matrix/iced
wdir=~/workspace/e.Blastocyst-K27ac/3.HiC/8.IS
conda activate py2.7
bash ~/codes/hic_insulation_score.sh -w $wdir -d $ddir -r 25000 --is 500000 --ids 250000 --nt 0.1 --bmoe 3
for i in 8cell ICM TE; do
awk '{a+=$3-$2;b+=1} END{print a,b,a/b}' $i.25000.iced.bed
done

sort -k1,1 -k2n,2 $bdir/ICM.25000.iced.bed | awk '{print $1,($2+$3-1)/2,($2+$3-1)/2+1}' | bedtools reldist -a - -b $pdir/ICM.bed | awk 'NR==1{print $0,"boundary","peak"} NR>1{print $0,"ICM","ICM"}' > K27ac.reldist.tab
sort -k1,1 -k2n,2 $bdir/ICM.25000.iced.bed | awk '{print $1,($2+$3-1)/2,($2+$3-1)/2+1}' | bedtools reldist -a - -b $pdir/both.bed | awk 'NR>1{print $0,"ICM","both"}' >> K27ac.reldist.tab
sort -k1,1 -k2n,2 $bdir/TE.25000.iced.bed | awk '{print $1,($2+$3-1)/2,($2+$3-1)/2+1}' | bedtools reldist -a - -b $pdir/TE.bed | awk 'NR>1{print $0,"TE","TE"}' >> K27ac.reldist.tab
sort -k1,1 -k2n,2 $bdir/TE.25000.iced.bed | awk '{print $1,($2+$3-1)/2,($2+$3-1)/2+1}' | bedtools reldist -a - -b $pdir/both.bed | awk 'NR>1{print $0,"TE","both"}' >> K27ac.reldist.tab

