# Preprocessing

Before running GraphBin, we have to assemble our read data into contigs and bin the contigs.

## Assembly
Reads can be assembled into contigs using 3 assembly software.

**metaSPAdes**

[**SPAdes**](http://cab.spbu.ru/software/spades/) is an assembler based on the de Bruijn graph approach. [**metaSPAdes**](https://genome.cshlp.org/content/27/5/824) is the dedicated metagenomic assembler of SPAdes. Use metaSPAdes (SPAdes in metagenomics mode) software to assemble reads into contigs.

**SGA**

[**SGA**](https://github.com/jts/sga) (String Graph Assembler) is an assembler based on the overlap-layout-consensus (more recently string graph) approach. Use SGA software to assemble reads into contigs.

**MEGAHIT**

[**MEGAHIT**](https://github.com/voutcn/megahit) is an assembler based on the de Bruijn graph approach. Use MEGAHIT software to assemble reads into contigs.

## Initial Binning

Once you have obtained the assembly output, you can run a metagenomic binning tool such as [MaxBin2](https://sourceforge.net/projects/maxbin2/), [CONCOCT](https://concoct.readthedocs.io/en/latest/), [MetaBAT2](https://bitbucket.org/berkeleylab/metabat) or [VAMB](https://github.com/RasmussenLab/vamb) to get an initial binning result.

You can use the [`prepResult.py` support script](https://github.com/Vini2/GraphBin/blob/master/support/prepResult.py) to format an initial binning result in to the .csv format with contig identifiers and bin ID. Contigs are named according to their original identifier and bins are numbered according to the fasta file name. You can run `prepResult.py` as follows.

```
python prepResult.py --binned /path/to/folder_with_binning_result --output /path/to/output_folder
```
You can see the usage options of `prepResult.py` by typing `python prepResult.py -h` on the command line.

Formatted binning result will be stored in a file named `initial_contig_bins.csv` in the output folder provided. Bin IDs and corresponding fasta files for each bin will be recorded in a file named `bin_ids.csv` in the output folder provided.

Now we are all set to run GraphBin.