function fetchref() {
    if [ -f $2.fasta ]; then
        echo "Reference genome for ${2} already exists, skipping download."
        return
    fi
    version=$(wget -qO- https://tritrypdb.org/common/downloads/Current_Release/Build_number)
    echo "TriTrypDB version ${version}"
    wget "https://tritrypdb.org/common/downloads/release-${version}/${1}/fasta/data/TriTrypDB-${version}_${1}_Genome.fasta" -O ${2}.fasta
}

function fetchsra() {
    if [ -f ${1}_1.fastq ] && [ -f ${1}_2.fastq ]; then
        echo "FASTQ files for ${1} already exist, skipping download."
        return
    fi
    prefetch ${1}
    fasterq-dump ${1}
    rm -r ${1}
}

fetchref LmexicanaMHOMGT2001U1103 Lmexicana
fetchsra SRR19895146

python3 -m polyalign splitfiltered Lmexicana.fasta SRR19895146_1.fastq SRR19895146_2.fastq genome-polyalign-splitfiltered
for fasta in genome-polyalign-splitfiltered/*.fasta; do
  filename=$(basename "$fasta")
  filename="${filename%.*}"
  polypolish polish $fasta genome-polyalign-splitfiltered_1/${filename}_1.sam genome-polyalign-splitfiltered_2/${filename}_2.sam >> genome-polyalign-splitfiltered-polishedgenome.fasta
done

#python3 -m polyalign filtered Lmexicana.fasta SRR19895146_1.fastq SRR19895146_2.fastq genome-polyalign-filtered
#python3 -m polyalign paired Lmexicana.fasta SRR19895146_1.fastq SRR19895146_2.fastq genome-polyalign-paired
