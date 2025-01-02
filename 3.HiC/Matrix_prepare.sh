res=5000
wdir=/home/qszhu/workspace/e.Blastocyst-K27ac/3.HiC/0.matrix/$res
ddir=/mnt/pangu/home2/qszhu/workspace/3.TE-HiC/3.processed/4.matrix-merged
mkd $wdir/logs
ln -s $ddir/*.$res.matrix $wdir/
cd $wdir
for i in *.$res.matrix; do
hicConvertFormat -m $i --inputFormat hicpro --bedFileHicpro $ConfigHP/${res}_mm10.bed --outputFormat h5 -o $wdir/${i/.matrix} --resolutions $res &
done
wait
hicNormalize -m NF-ICM.$res.h5 NF-TE.$res.h5 -n smallest -o NF-ICM.normed.$res.h5 NF-TE.normed.$res.h5

for i in ICM TE; do
hicCorrectMatrix diagnostic_plot --matrix NF-$i.normed.$res.h5 -o NF-$i.normed.$res.diagnostic.pdf &
done
wait
for i in ICM TE; do
hicCorrectMatrix correct --matrix NF-$i.normed.$res.h5 --correctionMethod ICE --iterNum 500  --outFileName NF-$i.iced.$res.h5 --filterThreshold -2 3 &
done
wait
for i in ICM TE; do
hicConvertFormat -m NF-$i.iced.$res.h5 --inputFormat h5 --outputFormat cool -o NF-$i.iced.$res.cool &
done
for i in ICM TE; do
hicPlotMatrix -m NF-$i.iced.$res.h5 -o NF-$i.hic_plot.png --region 7:45760000-47500000 --dpi 360 --log1p &
done
