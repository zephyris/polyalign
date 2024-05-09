export PATH=$PATH:~/sratoolkit/bin/

# get illumina data
function fetchsra() {
  prefetch ${1}
  fasterq-dump ${1}
  rm -r ${1}
}
illumina1=SRR19895146_1.fastq
illumina2=SRR19895146_2.fastq
if [ ! -f $illumina1 ]; then
  fetchsra SRR19895146
fi

if [ ! -f genome.fasta ]; then
  curl -L https://tritrypdb.org/common/downloads/release-68/LmexicanaMHOMGT2001U1103/fasta/data/TriTrypDB-68_LmexicanaMHOMGT2001U1103_Genome.fasta > genome.fasta
fi

python3 -m polyalign genome.fasta $illumina1 $illumina2 polyalign