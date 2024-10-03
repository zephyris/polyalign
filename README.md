# PolyAlign

A simple Python and BWA MEM-based aligner with improved read pairing.

`polyalign` aims to generate `sam`-format alignments which can be used for high accuracy polishing of repetitive sequences.

It is heavily inspired by `polypolish align` from [Polypolish](https://github.com/rrwick/Polypolish), but aims to reduce memory usages. This allows its use on eukaryote genomes.

## Installation

### Dependencies

PolyAlign requires `python`, and is easiest to install using `pip` and `git`. It also requires `bwa` to be installed and available to the system command line.

On Linux, install using your normal package manager, for example:
``` shell
sudo apt update
sudo apt install bwa python3 python3-pip git
```

Alternatively, `bwa` can be installed in a conda environment by running
```shell
conda install -c bioconda bwa
```
Please see the [Conda documentation](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) on how to install conda.

### Polyalign

Install using `pip` and `git`:
``` shell
pip install git+https://github.com/zephyris/polyalign
```

To reinstall and upgrade use `pip` and `git`:
``` shell
pip install --upgrade --force-reinstall git+https://github.com/zephyris/polyalign
```

To uninstall use `pip`
``` shell
pip uninstall polyalign
```

## Standalone usage

``` shell
python3 -m polyalign [filtered|filteredsplit|paired] <reference.fasta> <reads_1.fastq> <reads_2.fastq> <output_basename>
```

Reads have to be `fastq` format, `fasta` is not supported. In `filtered` mode, the output is two `sam` files, `<output_basename>_1.sam` and `<output_basename>_2.sam`.
In `filteredsplit` mode, the output is two `sam` files per sequence in the input `<reference.fasta>` file, output to the directories `<output_basename>_1` and `<output_basename>_2`.
In `paired` mode, the output is one `sam` file, `<output_basename>.sam`. In `paired` mode, `<output_basename>` of `-` will give `sam` output on `stdout`.

`polyalign` first samples a subset of reads from the start of both `fastq` files to identify orientation of read pairs and typical insert size from reads pairs aligned as a unique pair.
Next, it aligns the entire set of reads. Read pairs where one is not aligned and one is aligned to a single place are retained. Read pairs where both are aligned to a single place are retained.
Reads alignments to multiple places are retained, if a pair can be formed which gives typical insert size and correct orientation. For `filtered` and `filteredsplit` outputs, all such alignments are retained. For `paired`, one good pairing is randomly selected as the output.

In `filtered` mode, this behaviour broadly matches matches [Polypolish](https://github.com/rrwick/Polypolish) `polypolish filter`. The resulting output can be used for subsequent polishing using `polypolish polish`.

In `filteredsplit` mode, each individual `sam` file can be used for `polypolish polish` against the appropriate reference sequence.

In `paired` mode, this outputs a `sam` file similar to normal `bwa mem` paired alignments, and can be used for general subsequent analyses.

## Python module usage

You can use `polyalign` in your Python scripts - however it is subject to change. The `Polyalign` class carries out high-level operation, ouputting using the `Output` class.
`BwaMem` is to run `bwa mem` alignments. `Alignment` and `AlignmentClass` are used to parse alignments and alignment pairs.

## History and citing

This module was written to optimise polishing of small eukaryotic genomes assembled from noisy Nanopore data. The Polypolish strategy for polishing repetitive sequences is very promising, but designed for small (bacterial) genomes. Polyalign allows application of the same method to larger, eg. eukarotic genomes, without requiring enormous computational resources.

I haven't ultimately used it for a published genome assembly, but I've made it available in case it is useful. Please send me a message and cite this Github repository if you are publishing anything using this as a tool. Please also cite Polypolish, as this is _very_ closely modeled on that work.
