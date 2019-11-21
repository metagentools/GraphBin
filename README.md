# GraphBin: Improved Binning of Metagenomic Contigs using Assembly Graphs

**GraphBin** is a metagenomic contig binning tool that makes use of the contig connectivity information from the assembly graph to bin contigs. It utilizes the binning result of an existing binning tool and a label propagation algorithm to correct mis-binned contigs and predict the labels of contigs which are discarded due to short length.

## Dependencies
GraphBin is coded in Python 3.6. To run GraphBin and support scripts, you will need to install the following python modules.
* [python-igraph](https://igraph.org/python/)
* [Biopython](https://biopython.org/)

The [python-labelpropagation](https://github.com/ZwEin27/python-labelpropagation) module supporting Python 3 is provided in the source code.

You can go to these links and follow the instructions to download these modules.

## Downloading GraphBin
To download GraphBin, you have to clone the GraphBin repository to your machine.

```
git clone https://github.com/Vini2/GraphBin.git
```
Now go in to the source folder using the command
```
cd GraphBin/src/
```

## Assembly
The assembly of contigs can be done using 2 assembly software.

### SPAdes
[**SPAdes**](http://cab.spbu.ru/software/spades/) is an assembler based on the de Bruijn graph approach. Use SPAdes software to assemble reads into contigs. Use the metagenomics mode for assembly.

### SGA
[**SGA**](https://github.com/jts/sga) (String Graph Assembler) is an assembler based on the overlap-layout-consensus (more recently string graph) approach. Use SGA software to assemble reads into contigs.

Once you have obtained the assembly output, you can run GraphBin.

## Using GraphBin
You can see the usage options of GraphBin by typing ```python3 graphbin_SPAdes.py -h``` or ```python3 graphbin_SGA.py -h``` on the command line.

```
usage: graphbin_SPAdes.py [-h] --graph GRAPH --paths PATHS --binned BINNED
                          --output OUTPUT [--max_iteration [MAX_ITERATION]]
                          [--diff_threshold [DIFF_THRESHOLD]]

optional arguments:
  -h, --help            show this help message and exit
  --graph GRAPH         path to the assembly graph file
  --paths PATHS         path to the contigs.paths file
  --binned BINNED       path to the .csv file with the initial binning output
                        from an existing tool
  --output OUTPUT       path to the output folder
  --prefix [PREFIX]     prefix for the output file
  --max_iteration [MAX_ITERATION]
                        maximum number of iterations for label propagation
                        algorithm. [default: 100]
  --diff_threshold [DIFF_THRESHOLD]
                        difference threshold for label propagation algorithm.
                        [default: 0.1]
```
```
usage: graphbin_SGA.py [-h] --graph GRAPH --binned BINNED --output OUTPUT
                       [--max_iteration [MAX_ITERATION]]
                       [--diff_threshold [DIFF_THRESHOLD]]

optional arguments:
  -h, --help            show this help message and exit
  --graph GRAPH         path to the assembly graph file
  --binned BINNED       path to the .csv file with the initial binning output
                        from an existing tool
  --output OUTPUT       path to the output folder
  --prefix [PREFIX]     prefix for the output file
  --max_iteration [MAX_ITERATION]
                        maximum number of iterations for label propagation
                        algorithm. [default: 100]
  --diff_threshold [DIFF_THRESHOLD]
                        difference threshold for label propagation algorithm.
                        [default: 0.1]
```
`max_iteration` and `diff_threshold` parameters are set by default to `100` and `0.1` respectively. However, the user can specify them when running GraphBin.

## Input Format

`graphbin_SPAdes.py` takes in 3 files as inputs (required).
* Assembly graph file (in `.gfa` format)
* Paths of contigs (in `.paths` format)
* Binning output from an existing tool (in `.csv` format)

`graphbin_SGA.py` takes in 2 files as inputs (required).
* Assembly graph file (in `.asqg` format)
* Binning output from an existing tool (in `.csv` format)

**Note:** The binning output file should have comma separated values ```(contig_identifier, bin_number)``` for each contig. The contents of the binning output file should look similar to the example given below. The numbering of bins starts from 1.

SPAdes binned input
```
NODE_1,1
NODE_2,1
NODE_3,1
NODE_4,2
NODE_5,2
...
```
SGA binned input
```
contig-0,1
contig-1,2
contig-2,1
contig-3,1
contig-4,2
...
```
GraphBin provides a support script to generate binning result files. You can refer to [SupportREADME.md](https://github.com/Vini2/GraphBin/blob/master/support%20scripts/SupportREADME.md) file for more details.

## Example Usage

```
python3 graphbin_SPAdes.py --graph /path/to/graph_file.gfa --paths /path/to/paths_file.paths --binned /path/to/binning_result.csv --output /path/to/output_folder
```
```
python3 graphbin_SGA.py --graph /path/to/graph_file.asqg --binned /path/to/binning_result.csv --output /path/to/output_folder
```

## Support Scripts

GraphBin provides support scripts to format an initial binning result and visualise binning results in the assembly graph. Details about support scripts and how to execute them are provided in [SupportREADME.md](https://github.com/Vini2/GraphBin/blob/master/support%20scripts/SupportREADME.md) file.

## Test Data

The data used to test GraphBin can be found in the `test data` folder. The test data for each of the datasets include the following files.
* Contigs file
* Assembly graph file
* Paths file for the assembly graph (for the datasets assembled using SPAdes)
* Initial binning result from MaxBin 2.0
* Initial binning result from MetaWatt
* Ground truth labelling of contigs from TAXAassign

You can try running GraphBin using these test data files.

## Visualization of the Assembly Graph of ESC+SPAdes Test Dataset

### Initial Assembly Graph
<p align="center">
  <img src="images/3G_SPAdes_graph_plot.png" width="400" title="Initial assembly graph" alt="Initial assembly graph">
</p>

### TAXAassign Labelling
<p align="center">
  <img src="images/3G_SPAdes_taxaassign_graph_plot.png" width="400" title="TAXAassign Labelling" alt="TAXAassign Labelling">
</p>

### Original MaxBin Labelling with 2 Mis-binned Contigs
<p align="center">
  <img src="images/3G_SPAdes_maxbin_graph_plot_edit.png" width="400" title="MaxBin Labelling" alt="MaxBin Labelling">
</p>

### Refined Labels
<p align="center">
  <img src="images/3G_SPAdes_maxbin_graph_plot_correct.png" width="400" title="Refined Labels" alt="Refined Labels">
</p>

### Final Labelling of GraphBin
<p align="center">
  <img src="images/3G_SPAdes_after_label_prop_graph_plot.png" width="400" title="Final Labelling" alt="Final Labelling">
</p>

## References
[1] Wu, Y.W., Tang, Y.H., Tringe, S.G., Simmons, B.A., Singer, S.W.: MaxBin: an automated binning method to recover individual genomes from metagenomes using an expectation-maximization algorithm. Microbiome 2(26), (2014)

[2] Wu, Y.W., Simmons, B.A., Singer, S.W.: MaxBin 2.0: an automated binning algorithm to recover genomes from multiple metagenomic datasets. Bioinformatics 32(4), (2016)

[3] Strous, M. _et al._: The Binning of Metagenomic Contigs forMicrobial Physiology of Mixed Cultures. Frontiers in Microbiology 3, 410, (2012).

[4] Sczyrba, A., _et. al_: Critical Assessment of Metagenome Interpretation a Benchmark of Metagenomics Software. Nature Methods 14, 1063-1071 (2017)

[5] Barnum, T.P., _et al._: Genome-resolved metagenomics identifies genetic mobility, metabolic interactions, and unexpected diversity in perchlorate-reducing communities. The ISME Journal 12, 1568-1581 (2018)

[6] Bankevich, A., _et al._: SPAdes: A New Genome Assembly Algorithm and Its Applications to Single-Cell Sequencing. Journal of Computational Biology 19(5), 455-477 (2012)

[7] Nurk, S., Meleshko, D., Korobeynikov, A., Pevzner, P.A.: metaSPAdes: a new versatile metagenomic assembler. Genome Researcg 5, 824-834 (2017)

[8] Zhu, X., Ghahramani, Z.: Learning from Labeled and Unlabeled Data with Label Propagation. Technical Report CMU-CALD-02, Carnegie Mellon University, (2002)

[9] python-labelpropagation: Python implementation of label propagation, [https://github.com/ZwEin27/python-labelpropagation](https://github.com/ZwEin27/python-labelpropagation).

## Publication
The manuscript for GraphBin is currently under review at the _Bioinformatics_ (Oxford Academic) journal.
