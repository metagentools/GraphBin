# GraphBin Tutorial

This tutorial walks through the steps and commands used to set up GraphBin, prepare results for input, run GraphBin and visualise the final results. 

## Prerequisites

Make sure you have installed the following.

* [Conda](https://docs.conda.io/en/latest/miniconda.html)
* [Git](https://github.com/git-guides/install-git)

## Step 1 - Installing GraphBin

Let's create a new conda environment and install GraphBin from bioconda using the following command.
```
conda create -n graphbin -c bioconda graphbin
```

Activate the conda environment using,
```
conda activate graphbin
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

Let's set the path to our data as follows. You can use the path to our test data in [`tests/data`](https://github.com/Vini2/GraphBin/tree/master/tests/data).
```
mypath=/path/to/data/folder
```

### Step 2a - Assembly

We can assemble our reads into contigs using any metagenomic assembler. For this purpose, we will use [**metaSPAdes**](https://genome.cshlp.org/content/27/5/824) (available from [SPAdes](http://cab.spbu.ru/software/spades/)) as follows.
```
spades --meta -1 $mypath/Reads_1.fastq -2 $mypath/Reads_2.fastq -o $mypath/ -t 8
```

### Step 2b - Obtain the initial binning result

Any contig binning tool can be used to get an initial binning result. We will be using [MaxBin 2](https://sourceforge.net/projects/maxbin2/) in this example.


### Step 2c - Prepare the initial binning result

`prepResult.py` is a support script that allows you to format an initial binning result into the .csv format with contig identifiers and bin ID. Contigs are named according to their original identifier and bins are numbered according to the fasta file name. We can run `prepResult.py` as follows.

```
python support/prepResult.py --binned $mypath/maxbin_bins --output $mypath/
```

## Step 3 - Using GraphBin

We can run the **metaSPAdes** version of GraphBin as follows.
```
graphbin --assembler spades --graph $mypath/assembly_graph_with_scaffolds.gfa --contigs $mypath/contigs.fasta --paths $mypath/contigs.paths --binned $mypath/initial_contig_bins.csv --output $mypath/
```