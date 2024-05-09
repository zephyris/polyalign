from .polyalign import *

if __name__ == "__main__":
    print("Polyalign")
    print("-PA-", "indicates Polyalign output, other output is from the aligner.")
    if len(sys.argv) < 5:
        print("Usage: python3 -m polyalign [filter|paired] <reference.fasta> <reads_1.fastq> <reads_2.fastq> <output_basename>")
        print("", "Output will be <output_basename>_1.sam and <output_basename>_2.sam in filter mode or <output_basename>.sam in paired mode.")
        print("", "If in paired mode, <output_basename> can be - to output to stdout.")
        sys.exit(1)
    polyalign = Polyalign(sys.argv[2], sys.argv[3], sys.argv[4], output_basename=sys.argv[5], output_type=sys.argv[1])
    polyalign.polyalign()