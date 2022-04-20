# GraphBin Tutorial - FAME Fortnightly Bioinf Meeting

This tutorial walks through the steps and commands used to set up GraphBin, prepare results for input, run GraphBin and visualise the final results.

## Step 1 - Installing GraphBin

Let's create a new conda environment and install GraphBin from bioconda using the following command.
```
conda create -n graphbin -c bioconda graphbin
```

We can check if GraphBin is working properly using the following command.
```
graphbin -h
```

Now we can clone the GraphBin repository to our local machine.

```
git clone https://github.com/Vini2/GraphBin.git
```

Make sure you go into the GraphBin folder using the `cd` command.

```
cd GraphBin/
```

## Step 2 - Preprocessing

Let's set the path to our data as follows.
```
my_path=/path/to/data/folder
```

### Step 2a - Assembly

We can assemble our reads into contigs using any metagenomic assembler. For this purpose, we will use [**metaSPAdes**](https://genome.cshlp.org/content/27/5/824) (available from [SPAdes](http://cab.spbu.ru/software/spades/)) as follows.
```
spades --meta -1 $my_path/Reads_1.fastq -2 $my_path/Reads_2.fastq -o $my_path/ -t 8
```

### Step 2b - Obtain the initial binning result

Any contig binning tool can be used to get an initial binning result. We will be using [MaxBin 2](https://sourceforge.net/projects/maxbin2/) in this example.


### Step 2c - Prepare the initial binning result

`prepResult.py` is a support script that allows you to format an initial binning result into the .csv format with contig identifiers and bin ID. Contigs are named according to their original identifier and bins are numbered according to the fasta file name. We can run `prepResult.py` as follows.

```
python support/prepResult.py --binned $my_path/maxbin_bins --output $my_path/
```

## Step 3 - Using GraphBin

We can run the **metaSPAdes** version of GraphBin as follows.
```
graphbin --assembler spades --graph $my_path/assembly_graph_with_scaffolds.gfa --contigs $my_path/contigs.fasta --paths $my_path/contigs.paths --binned $my_path/initial_contig_bins.csv --output $my_path/
```

The final binning result from GraphBin can be found in the file `graphbin_output.csv`.

## Step 4 - Visualising the binning results

`visualiseResult_SPAdes.py` allows you to visualise the metaSPAdes binning result by denoting coloured contigs in the assembly graph according to their corresponding bins. We can visualise the initial binning result obtained from an existing binning tool and the final binning result obtained from GraphBin, and compare the differences.

```
python support/visualiseResult_SPAdes.py --graph $my_path/assembly_graph_with_scaffolds.gfa --paths $my_path/contigs.paths --initial $my_path/initial_contig_bins.csv --final $my_path/graphbin_output.csv --output $my_path/
```
