# GraphBin Tutorial

This tutorial walks you through the steps and commands used to setup GraphBin, prepare results, run GraphBin and visualise results.

## Installing GraphBin

Let's create a new conda environment and install GraphBin from bioconda using the following command.
```
conda create -n graphbin -c bioconda graphbin
```

## Preprocessing

Let's set the path to our data as follows
```
my_path=/path/to/data/folder
```

### Assembly

We can assemble our set of reads into contigs. For this purpose, we will use **metaSPAdes** as follows. [metaSPAdes](https://genome.cshlp.org/content/27/5/824) is the dedicated metagenomic assembler of [SPAdes](http://cab.spbu.ru/software/spades/)
```
spades --meta -1 $my_path/Reads_1.fastq -2 $my_path/Reads_2.fastq -o $my_path/ -t 8
```

### Initial binning

Any contig binning tool can be used to get an initial binning result. We will be using [MaxBin 2](https://sourceforge.net/projects/maxbin2/) in this example.


### Prepare binning results

`prepResult.py` is a support script that allows you to format an initial binning result in to the .csv format with contig identifiers and bin ID. Contigs are named according to their original identifier and bins are numbered according to the fasta file name. You can run `prepResult.py` as follows.

```
python support/prepResult.py --binned $my_path/maxbin_bins --output $my_path/
```

## Using GraphBin

You can run the metaSPAdes version of GraphBin as follows.
```
graphbin --assembler spades --graph $my_path/assembly_graph_with_scaffolds.gfa --contigs $my_path/contigs.fasta --paths $my_path/contigs.paths --binned $my_path/initial_contig_bins.csv --output $my_path/
```

The final binning result from GraphBin can be found in the file `graphbin_output.csv`.

## Visualising the binning results

`visualiseResult_SPAdes.py` allows you to visualize the metaSPAdes binning result by denoting coloured contigs in the assembly graph according to their corresponding bins. You can visualise the initial binning result obtained from an existing binning tool and the final binning result obtained from GraphBin and compare.

```
python support/visualiseResult_SPAdes.py --graph $my_path/assembly_graph_with_scaffolds.gfa --paths $my_path/contigs.paths --initial $my_path/initial_contig_bins.csv --final $my_path/graphbin_output.csv --output $my_path/
```
