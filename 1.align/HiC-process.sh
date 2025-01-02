ddir=~/workspace/e.Blastocyst-K27ac/0.prepare/2.cutadapt/temp
wdir=~/workspace/e.Blastocyst-K27ac/1.align/7.HiC/temp
mkcd $wdir;mkd logs
nohup bash ~/workspace/e.Blastocyst-K27ac/a.public/threads_hicpro_align.sh -o $wdir -d $ddir -k HiC -p 8 -t 8 -g mm10 -q 20 -l GATCGATC --hicpro /usr/local/software/HiC-Pro-3.1.0/scripts/ > $wdir/hicpro.align.log &
wait-m align

mv temp/logs/* logs/*
mv temp/stats/* stats/*
mv temp/tmp/* tmp/*
mv temp/* ./

nohup bash ~/workspace/e.Blastocyst-K27ac/a.public/threads_hicpro_pairs.sh -w $wdir --hicpro /usr/local/software/HiC-Pro-3.1.0/scripts/ --python /usr/bin/python -k HiC -m 5 -p 14 -t 8 -l $ConfigHP/mm10_MboI.bed > $wdir/contact.log &
wait-m hicpro


####merge rep
wdir=~/workspace/e.Blastocyst-K27ac/a.public/6.Liujiang/4.HiC-Pro
for i in TE; do
(cat $wdir/${i}_HiC_rep*rmdup.pairs > ${i}_HiC.temp.pairs
sort -S 50% --parallel=64 -k2,2V -k3,3n -k5,5V -k6,6n ${i}_HiC.temp.pairs > ${i}_HiC.merged.pairs
)&
done;
wait
rm *temp.pairs

nohup bash ~/workspace/e.Blastocyst-K27ac/a.public/threads_hicpro_buildMatrix.sh -w $wdir -k merged -r "10000 25000 50000 100000 1000000" --hicpro /usr/local/software/HiC-Pro-3.1.0/scripts/ --python /usr/bin/python --output_bias 0 -c $ConfigHP/chrom_mm10.sizes -p 10 > $wdir/logs/buildMatrix.log