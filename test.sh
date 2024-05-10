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

export PATH=$PATH:~/sratoolkit/bin/:~/polypolish/target/release/

python3 -m pip install git+https://github.com/zephyris/polyalign

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

python3 -m polyalign filtered genome.fasta $illumina1 $illumina2 polyalign
polypolish polish genome.fasta polyalign_1.sam polyalign_2.sam > genome-polyalign.fasta

bwa index genome.fasta
bwa mem -a -Y -t 40 genome.fasta $illumina1 > bwa_1.sam
bwa mem -a -Y -t 40 genome.fasta $illumina2 > bwa_2.sam
polypolish filter --in1 bwa_1.sam --in2 bwa_2.sam --out1 bwa_1_filtered.sam --out2 bwa_2_filtered.sam
polypolish polish genome.fasta bwa_1_filtered.sam bwa_2_filtered.sam > genome-polypolish_filter.fasta
