from .polyalign import *

if __name__ == "__main__":
    def print_help():
        print("Usage: python3 -m polyalign [filtered|paired|splitfiltered] <reference.fasta> <reads_1.fastq> <reads_2.fastq> <output_basename>")
        print("", "In filtered mode, output will be <output_basename>_1.sam and <output_basename>_2.sam")
        print("", "In splitfiltered mode, output will two sam files per contig in <output_basename>_1/<contig>.sam and <output_basename>_2/<contig>.sam")
        print("", "In paired mode, output will be <output_basename>.sam.")
        print("", "", "In paired mode, <output_basename> can be - to output to stdout.")
        print("or")
        print("Usage: python3 -m polyalign [splitfasta] <reference.fasta> <output_basename>")
        print("", "In splitfasta mode, output will be <output_basename>/<contig>.fasta")

    print("Polyalign")
    print("[PA::*]", "indicates Polyalign output, other output is from the aligner.")
    if len(sys.argv) < 5:
        if len(sys.argv) == 4 and sys.argv[1] == "splitfasta":
            splitfasta = Splitfasta(sys.argv[2], output_basename=sys.argv[3])
            splitfasta.splitfasta()
            exit(0)
        else:
            print_help()
            sys.exit(1)
    elif len(sys.argv) == 6:
        polyalign = Polyalign(sys.argv[2], sys.argv[3], sys.argv[4], output_basename=sys.argv[5], output_type=sys.argv[1])
        polyalign.polyalign()
    else:
        print_help()
        sys.exit(1)
