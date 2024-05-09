from .polyalign import *

# python3 polyalign.py genome.fasta reads_1.fastq reads_2.fastq
if __name__ == "__main__":
    print("Polyalign")
    print("-PA-", "indicates Polyalign output")
    if len(sys.argv) < 5:
        print("Usage: python3 polyalign.py genome.fasta reads_1.fastq reads_2.fastq output_basename")
        sys.exit(1)
    polyalign = Polyalign(sys.argv[1], sys.argv[2], sys.argv[3], output_basename=sys.argv[4])
    polyalign.polyalign()
