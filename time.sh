genome=Lmexicana.fasta
illumina1=SRR19895146_1.fastq
illumina2=SRR19895146_2.fastq

lines=$(/usr/bin/time -v ls | wc -l)
cpu=$(nproc)
/usr/bin/time -v bwa index $genome 2> >(tail -n $lines > bwa-index.time.txt)
/usr/bin/time -v bwa mem -a -Y -t $cpu $genome $illumina1 > bwa_1.sam 2> >(tail -n $lines > bwa-mem.time.txt)
/usr/bin/time -v polypolish filter --in1 bwa_1.sam --in2 bwa_2.sam --out1 bwa_1_filtered.sam --out2 bwa_2_filtered.sam 2> polypolish-filter.time.txt
/usr/bin/time -v polypolish polish $genome bwa_1_filtered.sam bwa_2_filtered.sam > genome-polypolishfilter.fasta 2> >(tail -n $lines > polypolish-polish.time.txt)

/usr/bin/time -v python3 -m polyalign splitfiltered $genome $illumina1 $illumina2 polyalign-splitfiltered 2> polyalign-splitfiltered.time.txt
/usr/bin/time -v python3 -m polyalign splitfasta $genome polyalign-splitfiltered 2> polyalign-splitfasta.time.txt
if [ -f polyalign-polish.time.txt ]; then
  rm polyalign-polish.time.txt
fi
for fasta in polyalign-splitfiltered/*.fasta; do
  filename=$(basename "$fasta")
  filename="${filename%.*}"
  /usr/bin/time -v polypolish polish $fasta polyalign-splitfiltered_1/$filename.sam polyalign-splitfiltered_2/$filename.sam >> genome-polyalign-splitfiltered.fasta 2> >(tail -n $lines >> polyalign-polish.time.txt)
done