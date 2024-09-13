# install dependencies
BASEDIR=$(pwd)
cd ~

if [ ! -d ~/sratoolkit ]; then
  curl -L https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.1.0/sratoolkit.3.1.0-ubuntu64.tar.gz | tar jxf -
  mv sratoolkit.3.1.0-ubuntu64 ~/sratoolkit
fi

if ! command -v cargo &> /dev/null; then
  apt update && sudo apt install -y cargo
fi

if [ ! -d polypolish ]; then
  git clone https://github.com/rrwick/Polypolish polypolish
  cd polypolish && cargo build --release && cd ..
fi

if ! command -v bwa &> /dev/null; then
  apt update && sudo apt install -y bwa
fi

if ! command -v git &> /dev/null; then
  apt update && sudo apt install -y git
fi

cd $BASEDIR

python3 -m pip install git+https://github.com/zephyris/polyalign

# add tools to path
export PATH=$PATH:~/sratoolkit/bin/:~/polypolish/target/release/

# get illumina and genome data
# Leishmania mexicana genome and Illumina genome sequencing data
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

genome=Lmexicana.fasta
if [ ! -f $genome ]; then
  curl -L https://tritrypdb.org/common/downloads/release-68/LmexicanaMHOMGT2001U1103/fasta/data/TriTrypDB-68_LmexicanaMHOMGT2001U1103_Genome.fasta > $genome
fi

# do the tests

# with polypolish
# polyalign in splitfiltered mode
python3 -m polyalign splitfiltered $genome $illumina1 $illumina2 polyalign-splitfiltered
python3 -m polyalign splitfasta $genome polyalign-splitfiltered
if [ -f genome-polyalign-splitfiltered.fasta ]; then
  rm genome-polyalign-splitfiltered.fasta
fi
for fasta in polyalign-splitfiltered/*.fasta; do
  filename=$(basename "$fasta")
  filename="${filename%.*}"
  polypolish polish $fasta polyalign-splitfiltered_1/$filename.sam polyalign-splitfiltered_2/$filename.sam >> genome-polyalign-splitfiltered.fasta
done

# polyalign in filtered mode
python3 -m polyalign filtered $genome $illumina1 $illumina2 polyalign-filtered
polypolish polish $genome polyalign-filtered_1.sam polyalign-filtered_2.sam > genome-polyalign-filtered.fasta

# bwa and polypolish filter
cpu=$(nproc)
bwa index $genome
bwa mem -a -Y -t $cpu $genome $illumina1 > bwa_1.sam
bwa mem -a -Y -t $cpu $genome $illumina2 > bwa_2.sam
polypolish filter --in1 bwa_1.sam --in2 bwa_2.sam --out1 bwa_1_filtered.sam --out2 bwa_2_filtered.sam
polypolish polish $genome bwa_1_filtered.sam bwa_2_filtered.sam > genome-polypolishfilter.fasta

# polyalign in paired mode (no polish)
python3 -m polyalign paired $genome $illumina1 $illumina2 polyalign-paired
