# Using GraphBin

## GraphBin Options
You can see the usage options of GraphBin by typing `graphbin -h` on the command line. For example,

```
usage: graphbin [-h] [--version] [--graph GRAPH] [--binned BINNED]
                [--output OUTPUT] [--prefix PREFIX]
                [--max_iteration MAX_ITERATION]
                [--diff_threshold DIFF_THRESHOLD] [--assembler ASSEMBLER]
                [--paths PATHS] [--contigs CONTIGS] [--delimiter DELIMITER]

GraphBin Help. GraphBin is a metagenomic contig binning tool that makes use of
the contig connectivity information from the assembly graph to bin contigs. It
utilizes the binning result of an existing binning tool and a label
propagation algorithm to correct mis-binned contigs and predict the labels of
contigs which are discarded due to short length.

optional arguments:
  -h, --help            show this help message and exit
  --version
  --graph GRAPH         path to the assembly graph file
  --binned BINNED       path to the .csv file with the initial binning output
                        from an existing tool
  --output OUTPUT       path to the output folder
  --prefix PREFIX       prefix for the output file
  --max_iteration MAX_ITERATION
                        maximum number of iterations for label propagation
                        algorithm. [default: 100]
  --diff_threshold DIFF_THRESHOLD
                        difference threshold for label propagation algorithm.
                        [default: 0.1]
  --assembler ASSEMBLER
                        name of the assembler used (SPAdes, SGA or MEGAHIT).
                        GraphBin supports Flye, Canu and Miniasm long-read
                        assemblies as well.
  --paths PATHS         path to the contigs.paths file, only needed for SPAdes
  --contigs CONTIGS     path to the contigs.fa file.
  --delimiter DELIMITER
                        delimiter for input/output results. Supports a comma
                        (,), a semicolon (;), a tab ($'\t'), a space (" ") and
                        a pipe (|) [default: , (comma)]
```

`max_iteration` and `diff_threshold` parameters are set by default to `100` and `0.1` respectively. However, the user can specify them when running GraphBin.

## Inputs

For the SPAdes version, `graphbin` takes in 3 files as inputs (required).

* Assembly graph file (in `.gfa` format)
* Contigs file (in FASTA format)
* Paths of contigs (in `.paths` format)
* Binning output from an existing tool (in `.csv` format)

For the SGA version, `graphbin` takes in 2 files as inputs (required).

* Assembly graph file (in `.asqg` format)
* Contigs file (in FASTA format)
* Binning output from an existing tool (in `.csv` format)

For the MEGAHIT version, `graphbin` takes in 3 files as inputs (required).

* Assembly graph file (in `.gfa` format. To convert fastg to gfa refer [here](https://github.com/Vini2/GraphBin/blob/master/support/README.md#fastg2gfa))
* Contigs file (in FASTA format)
* Binning output from an existing tool (in `.csv` format)

**Note:** Make sure that the initial binning result consists of contigs belonging to only one bin. GraphBin is designed to handle initial contigs which belong to only one bin. Multiple bins for the initial contigs are not supported.

**Note:** You can specify the delimiter for the initial binning result file and the final output file using the delimiter paramter. Enter the following values for different delimiters; `,` for a comma, `;` for a semicolon, `$'\t'` for a tab, `" "` for a space and `|` for a pipe.

**Note:** The binning output file should have comma separated values ```(contig_identifier, bin_identifier)``` for each contig. The contents of the binning output file should look similar to the example given below. Contigs are named according to their original identifier and bin identifier.

Example metaSPAdes binned input
```
NODE_1_length_458813_cov_136.660185,bin_1
NODE_2_length_409135_cov_127.776630,bin_1
NODE_3_length_346431_cov_35.887787,bin_2
NODE_4_length_333835_cov_129.134427,bin_1
NODE_5_length_282512_cov_36.193003,bin_2
...
```
Example SGA binned input
```
contig-0,bin_1
contig-1,bin_2
contig-2,bin_1
contig-3,bin_1
contig-4,bin_2
...
```
Example MEGAHIT binned input
```
k99_10059,bin_1
k99_9367,bin_1
k99_15595,bin_2
k99_18709,bin_1
k99_15596,bin_2
...
```
Make sure to follow the steps provided in the [preprocessing section](https://graphbin.readthedocs.io/en/latest/preprocess/) to prepare the initial binning result.

## Example Usage

```
graphbin --assembler spades --graph /path/to/graph_file.gfa --contigs /path/to/contigs.fasta --paths /path/to/paths_file.paths --binned /path/to/binning_result.csv --output /path/to/output_folder
```
```
graphbin --assembler sga --graph /path/to/graph_file.asqg --contigs /path/to/contigs.fa --binned /path/to/binning_result.csv --output /path/to/output_folder
```
```
graphbin --assembler megahit --graph /path/to/graph_file.gfa --contigs /path/to/contigs.fa --binned /path/to/binning_result.csv --output /path/to/output_folder
```

## Output

The output from GraphBin will be a `.csv` file with comma separated values ```(contig_identifier, bin_identifier)``` for the refined binning result and the `.fasta` files of the refined bins.