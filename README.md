# OrfFinder

A simple Python and BWA MEM-based aligner with improved read pairing.

## Installation

### Dependencies

OrfFinder requires `python`, and is easiest to install using `pip` and `git`. It also requires `bwa` to be installed and available to the system command line.

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
python3 -m polyalign <genome.fasta> <reads_1.fastq> <reads_2.fastq> <output_basename>
```

Reads have to be `fastq` format, `fasta` is not supported. The output is two `sam` files, `<output_basename>_1.sam` and `<output_basename>_2.sam`.

## Python module usage

You can use `polyalign` in your Python scripts - however it is subject to change. The `Polyalign` class carries out high-level operation. `BwaMem` is to run `bwa mem` alignments. `Alignment` and `AlignmentClass` are used to parse alignments and alignment pairs.
