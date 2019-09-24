`prepResult.py` is a support script that allows you to format an initial binning result in to the .csv format with contig ID and bin ID. Contigs are numbered starting from 0 and bins are numbered starting from 1. You can run `prepResult.py` as follows.

```
python3 prepResult.py --binned /path/to/folder_with_binning_result --assembler assembler_type_(SPAdes/SGA) --output /path/to/output_folder
```
You can see the usage options of `prepResult.py` by typing `python3 prepResult.py -h` on the command line.

Formatted binning result will be stored in a file named `initial_contig_bins.csv` in the output folder provided. Bin IDs and corresponding fasta files for each bin will be recorded in a file named `bin_ids.csv` in the output folder provided.