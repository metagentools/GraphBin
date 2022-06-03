<p align="center">
  <img src="https://raw.githubusercontent.com/metagentools/GraphBin/master/GraphBin_logo.png" width="400" title="Final Labelling" alt="Final Labelling">
</p>

# GraphBin: Refined Binning of Metagenomic Contigs using Assembly Graphs

[![CI](https://github.com/metagentools/GraphBin/actions/workflows/testing_python_app.yml/badge.svg)](https://github.com/metagentools/GraphBin/actions/workflows/testing_python_app.yml)
[![codecov](https://codecov.io/gh/metagentools/GraphBin/branch/develop/graph/badge.svg?token=0S310F6QXJ)](https://codecov.io/gh/metagentools/GraphBin)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/Vini2/GraphBin.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/Vini2/GraphBin/context:python)
[![Documentation Status](https://readthedocs.org/projects/graphbin/badge/?version=latest)](https://graphbin.readthedocs.io/en/latest/?badge=latest)

**GraphBin** is a NGS data-based metagenomic contig bin refinment tool that makes use of the contig connectivity information from the assembly graph to bin contigs. It utilizes the binning result of an existing binning tool and a label propagation algorithm to correct mis-binned contigs and predict the labels of contigs which are discarded due to short length.

**For detailed instructions on installation, usage and visualisation, please refer to the [documentation hosted at Read the Docs](https://graphbin.readthedocs.io/).**

## Dependencies

GraphBin installation requires python 3 (tested on Python 3.6 and 3.7). The following dependencies are required to run GraphBin and related support scripts.
* [python-igraph](https://igraph.org/python/) - version 0.7.1
* [biopython](https://biopython.org/) - version 1.74
* [cairocffi](https://pypi.org/project/cairocffi/)

## Installing GraphBin

### Using Conda

You can install GraphBin via [Conda](https://docs.conda.io/en/latest/). You can download [Anaconda](https://www.anaconda.com/distribution/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) which contains Conda.

```
# create conda environment and install
conda create -n graphbin -c bioconda graphbin

# activate conda environment
conda activate graphbin

# check graphbin installation
graphbin -h
```

### Using pip

You can install GraphBin using pip.

```
pip install graphbin
```

For ***development*** purposes, please clone the repository and install via [flit](https://pypi.org/project/flit/).

```
# clone repository to your local machine
git clone https://github.com/metagentools/GraphBin.git

# go to repo direcotry
cd GraphBin

# install flit
pip install flit

# install graphbin via flit
flit install -s --python `which python`
```

## Example Usage

```
# SPAdes version
graphbin --assembler spades --graph /path/to/graph_file.gfa --contigs /path/to/contigs.fasta --paths /path/to/paths_file.paths --binned /path/to/binning_result.csv --output /path/to/output_folder

# SGA version
graphbin --assembler sga --graph /path/to/graph_file.asqg --contigs /path/to/contigs.fa --binned /path/to/binning_result.csv --output /path/to/output_folder

# MEGAHIT version
graphbin --assembler megahit --graph /path/to/graph_file.gfa --contigs /path/to/contigs.fa --binned /path/to/binning_result.csv --output /path/to/output_folder
```

## Visualization of the Assembly Graph of ESC+metaSPAdes Test Dataset

### Initial Assembly Graph
<p align="center">
  <img src="https://raw.githubusercontent.com/metagentools/GraphBin/master/images/3G_SPAdes_graph_plot.png" width="400" title="Initial assembly graph" alt="Initial assembly graph">
</p>

### TAXAassign Labelling
<p align="center">
  <img src="https://raw.githubusercontent.com/metagentools/GraphBin/master/images/3G_SPAdes_taxaassign_graph_plot.png" width="400" title="TAXAassign Labelling" alt="TAXAassign Labelling">
</p>

### Original MaxBin Labelling with 2 Mis-binned Contigs
<p align="center">
  <img src="https://raw.githubusercontent.com/metagentools/GraphBin/master/images/3G_SPAdes_maxbin_graph_plot_edit.png" width="400" title="MaxBin Labelling" alt="MaxBin Labelling">
</p>

### Refined Labels
<p align="center">
  <img src="https://raw.githubusercontent.com/metagentools/GraphBin/master/images/3G_SPAdes_maxbin_graph_plot_correct.png" width="400" title="Refined Labels" alt="Refined Labels">
</p>

### Final Labelling of GraphBin
<p align="center">
  <img src="https://raw.githubusercontent.com/metagentools/GraphBin/master/images/3G_SPAdes_after_label_prop_graph_plot.png" width="400" title="Final Labelling" alt="Final Labelling">
</p>


## Citation
If you use GraphBin in your work, please cite GraphBin as,

Vijini Mallawaarachchi, Anuradha Wickramarachchi, Yu Lin. GraphBin: Refined binning of metagenomic contigs using assembly graphs. Bioinformatics, Volume 36, Issue 11, June 2020, Pages 3307–3313, DOI: [10.1093/bioinformatics/btaa180](http://dx.doi.org/10.1093/bioinformatics/btaa180)

```bibtex
@article{10.1093/bioinformatics/btaa180,
    author = {Mallawaarachchi, Vijini and Wickramarachchi, Anuradha and Lin, Yu},
    title = "{GraphBin: refined binning of metagenomic contigs using assembly graphs}",
    journal = {Bioinformatics},
    volume = {36},
    number = {11},
    pages = {3307-3313},
    year = {2020},
    month = {03},
    abstract = "{The field of metagenomics has provided valuable insights into the structure, diversity and ecology within microbial communities. One key step in metagenomics analysis is to assemble reads into longer contigs which are then binned into groups of contigs that belong to different species present in the metagenomic sample. Binning of contigs plays an important role in metagenomics and most available binning algorithms bin contigs using genomic features such as oligonucleotide/k-mer composition and contig coverage. As metagenomic contigs are derived from the assembly process, they are output from the underlying assembly graph which contains valuable connectivity information between contigs that can be used for binning.We propose GraphBin, a new binning method that makes use of the assembly graph and applies a label propagation algorithm to refine the binning result of existing tools. We show that GraphBin can make use of the assembly graphs constructed from both the de Bruijn graph and the overlap-layout-consensus approach. Moreover, we demonstrate improved experimental results from GraphBin in terms of identifying mis-binned contigs and binning of contigs discarded by existing binning tools. To the best of our knowledge, this is the first time that the information from the assembly graph has been used in a tool for the binning of metagenomic contigs.The source code of GraphBin is available at https://github.com/Vini2/GraphBin.vijini.mallawaarachchi@anu.edu.au or yu.lin@anu.edu.auSupplementary data are available at Bioinformatics online.}",
    issn = {1367-4803},
    doi = {10.1093/bioinformatics/btaa180},
    url = {https://doi.org/10.1093/bioinformatics/btaa180},
    eprint = {https://academic.oup.com/bioinformatics/article-pdf/36/11/3307/33329097/btaa180.pdf},
}
```
