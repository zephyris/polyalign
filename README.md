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
python3 -m polyalign [filter|paired] <reference.fasta> <reads_1.fastq> <reads_2.fastq> <output_basename>
```

Reads have to be `fastq` format, `fasta` is not supported. In `filter` mode, the output is two `sam` files, `<output_basename>_1.sam` and `<output_basename>_2.sam`.
In `paired` mode, the output is one `sam` file, `<output_basename>.sam`. In `paired` mode, `<output_basename>` of `-` will give `sam` output on `stdout`.

`polyalign` first samples a subset of reads from the start of both `fastq` files to identify orientation of read pairs and typical insert size from reads pairs aligned as a unique pair.
Next, it aligns the entire set of reads. Read pairs where one is not aligned and one is aligned to a single place are retained. Read pairs where both are aligned to a single place are retained.
Read pairs which aligned to multiple places are retained, if a pair can be formed which gives typical insert size and correct orientation.

In `filter` mode, this behaviour broadly matches matches (Polypolish)[https://github.com/rrwick/Polypolish] `polypolish filter`, and can be used for subsequent `polypolish polish`. There are a few exceptions:
Reads which are not retained are discarded from the SAM instead of marked with the `ZP:Z:fail` flag, and read pairs where both are not aligned are discarded.

In `paired` mode, this outputs a `sam` file similar to normal `bwa mem` paired alignments, and can be used for general subsequent analyses.

## Python module usage

You can use `polyalign` in your Python scripts - however it is subject to change. The `Polyalign` class carries out high-level operation.
`BwaMem` is to run `bwa mem` alignments. `Alignment` and `AlignmentClass` are used to parse alignments and alignment pairs.
