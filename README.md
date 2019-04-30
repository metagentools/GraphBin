# GraphBin: Improved Binning of Metagenomic Contigs using Assembly Graphs

**GraphBin** is a metagenomic contig binning tool that makes use of the contig connectivity information from the assembly graph to bin contigs. It utilizes the binning result of an existing binning tool and a label propagation algorithm to correct mis-binned contigs and predict the labels of contigs which are discarded due to short length.

## Dependencies
GraphBin is coded in Python 3.6. To run GraphBin, you will need to install the following python modules.
* [python-igraph](https://igraph.org/python/)

The [python-labelpropagation](https://github.com/ZwEin27/python-labelpropagation) module supporting Python 3 is provided in the source code.

You can go to these links and follow the instructions to download these modules.

## Downloading GraphBin
To download GraphBin, you have to clone the GraphBin repository to your machine.

```
git clone https://github.com/Vini2/GraphBin.git
```
## Assembly
The assembly of contigs can be done using 2 software

### SPAdes
Use [**SPAdes**](http://cab.spbu.ru/software/spades/) software to assemble reads into contigs. Use the metagenomics mode for assembly.

### SGA
Use [**SGA**](https://github.com/jts/sga) (String Graph Assembler) software to assemble reads into contigs.

Once you have obtained the assembly output, you can run GraphBin.

## Using GraphBin
You can see the usage options of GraphBin by typing ```python graphbin_SPAdes.py -h``` or ```python graphbin_SGA.py -h``` on the command line.

```
usage: graphbin_SGA.py [-h] --graph GRAPH --paths PATHS --n_bins N_BINS --binned BINNED --output OUTPUT

optional arguments:
  -h, --help         show this help message and exit
  --graph GRAPH      path to the assembly graph file
  --paths PATHS      path to the contigs.paths file
  --binned BINNED    path to the .csv file with the initial binning output
                     from an existing tool
  --output OUTPUT    path to the output folder
```
```
usage: graphbin_SGA.py [-h] --graph GRAPH --n_contigs N_CONTIGS --binned BINNED --output OUTPUT

optional arguments:
  -h, --help             show this help message and exit
  --graph GRAPH          path to the assembly graph file
  --binned BINNED        path to the .csv file with the initial binning output
                         from an existing tool
  --output OUTPUT        path to the output folder
```
## Input Format

graphbin_SPAdes.py takes in 3 files as inputs.
* Assembly graph file (in .gfa format)
* Paths of contigs (in .paths format)
* Binning output from an existing tool (in .csv format)

graphbin_SGA.py takes in 1 file and 1 value as inputs.
* Assembly graph file (in .asqg format)
* Binning output from an existing tool (in .csv format)

**Note:** The binning output file should have comma separated values ```(node_number, bin_number)``` for each contig. The contents of the binning output file should look similar to the example given below. The numbering of contigs starts from 0 and numbering of bins starts from 1.

```
0,1
1,1
2,1
3,2
4,2
5,1
...
```

## Example Usage

```
python graphbin_SPAdes.py --graph /path/to/graph_file.gfa --paths /path/to/paths_file.paths --binned /path/to/binning_result.csv --output /path/to/output_folder
```
```
python graphbin_SGA.py --graph /path/to/graph_file.gfa --binned /path/to/binning_result.csv --output /path/to/output_folder
```

## Test Data

The data used to test GraphBin can be found in the `test data` folder. You can try running GraphBin using this test data.

## References
[1] Wu, Y.W., Tang, Y.H., Tringe, S.G., Simmons, B.A., Singer, S.W.: MaxBin: an automated binning method to recover individual genomes from metagenomes using an expectation-maximization algorithm. Microbiome 2(26), (2014)

[2] Wu, Y.W., Simmons, B.A., Singer, S.W.: MaxBin 2.0: an automated binning algorithm to recover genomes from multiple metagenomic datasets. Bioinformatics 32(4), (2016)

[3] Sczyrba, A., et. al : Critical Assessment of Metagenome Interpretation a Benchmark of Metagenomics Software. Nature Methods 14, 1063-1071 (2017)

[4] Barnum, T.P., et al.: Genome-resolved metagenomics identifies genetic mobility, metabolic interactions, and unexpected diversity in perchlorate-reducing communities. The ISME Journal 12, 1568-1581 (2018)

[5] Bankevich, A., et al.: SPAdes: A New Genome Assembly Algorithm and Its Applications to Single-Cell Sequencing. Journal of Computational Biology 19(5), 455-477 (2012)

[6] Nurk, S., Meleshko, D., Korobeynikov, A., Pevzner, P.A.: metaSPAdes: a new versatile metagenomic assembler. Genome Researcg 5, 824-834 (2017)

[7] Zhu, X., Ghahramani, Z.: Learning from Labeled and Unlabeled Data with Label Propagation. Technical Report CMU-CALD-02, Carnegie Mellon University, (2002)

[8] python-labelpropagation: Python implementation of label propagation, [https://github.com/ZwEin27/python-labelpropagation](https://github.com/ZwEin27/python-labelpropagation).
