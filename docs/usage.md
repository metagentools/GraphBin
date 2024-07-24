# Using GraphBin

## GraphBin Options
You can see the usage options of GraphBin by typing `graphbin -h` on the command line. For example,

```
Usage: graphbin [OPTIONS]

  GraphBin: Refined Binning of Metagenomic Contigs using Assembly Graphs

Options:
  --assembler [spades|sga|megahit|flye|canu|miniasm]
                                  name of the assembler used (SPAdes, SGA or
                                  MEGAHIT). GraphBin supports Flye, Canu and
                                  Miniasm long-read assemblies as well.
                                  [required]
  --graph PATH                    path to the assembly graph file  [required]
  --contigs PATH                  path to the contigs file  [required]
  --paths PATH                    path to the contigs.paths (metaSPAdes) or
                                  assembly.info (metaFlye) file
  --binned PATH                   path to the .csv file with the initial
                                  binning output from an existing tool
                                  [required]
  --output PATH                   path to the output folder  [required]
  --prefix TEXT                   prefix for the output file
  --max_iteration INTEGER         maximum number of iterations for label
                                  propagation algorithm  [default: 100]
  --diff_threshold FLOAT RANGE    difference threshold for label propagation
                                  algorithm  [default: 0.1; 0<=x<=1]
  --delimiter [,|;|$'\t'|" "]     delimiter for input/output results. Supports
                                  a comma (,), a semicolon (;), a tab ($'\t'),
                                  a space (" ") and a pipe (|)  [default: ,]
  -v, --version                   Show the version and exit.
  --help                          Show this message and exit.
```

`max_iteration` and `diff_threshold` parameters are set by default to `100` and `0.1` respectively. However, the user can specify them when running GraphBin.

## Inputs

For the SPAdes version, `graphbin` takes in 3 files as inputs (required).

* Assembly graph file (in `.gfa` format)
* Contigs file (`contigs.fasta` file in FASTA format)
* Paths of contigs (`contigs.paths` file)
* Binning output from an existing tool (in `.csv` format)

For the SGA version, `graphbin` takes in 2 files as inputs (required).

* Assembly graph file (in `.asqg` format)
* Contigs file (in FASTA format)
* Binning output from an existing tool (in `.csv` format)

For the MEGAHIT version, `graphbin` takes in 3 files as inputs (required).

* Assembly graph file (in `.gfa` format. To convert fastg to gfa refer [here](https://github.com/Vini2/GraphBin/blob/master/support/README.md#fastg2gfa))
* Contigs file (in FASTA format)
* Binning output from an existing tool (in `.csv` format)

For the Flye version, `graphbin` takes in 3 files as inputs (required).

* Assembly graph file (in `.gfa` format)
* Contigs file (`assembly.fasta` file in FASTA format)
* Paths of contigs (`assembly_info.txt` file)
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
Example Flye binned input
```
contig_1,bin_1
contig_2,bin_1
contig_3,bin_2
contig_4,bin_1
contig_5,bin_2
...
```

Make sure to follow the steps provided in the [preprocessing section](https://graphbin.readthedocs.io/en/latest/preprocess/) to prepare the initial binning result.

## Example Usage

```
# metaSPAdes assembly
graphbin --assembler spades --graph /path/to/graph_file.gfa --contigs /path/to/contigs.fasta --paths /path/to/paths_file.paths --binned /path/to/binning_result.csv --output /path/to/output_folder
```
```
# SGA assembly
graphbin --assembler sga --graph /path/to/graph_file.asqg --contigs /path/to/contigs.fa --binned /path/to/binning_result.csv --output /path/to/output_folder
```
```
# MEGAHIT assembly
graphbin --assembler megahit --graph /path/to/graph_file.gfa --contigs /path/to/contigs.fa --binned /path/to/binning_result.csv --output /path/to/output_folder
```
```
# metaFlye assembly
graphbin --assembler flye --graph /path/to/assembly_graph.gfa --contigs /path/to/assembly.fasta --paths /path/to/assembly_info.txt --binned /path/to/binning_result.csv --output /path/to/output_folder
```

## Output

The output from GraphBin will be a `.csv` file with comma separated values ```(contig_identifier, bin_identifier)``` for the refined binning result and the `.fasta` files of the refined bins.