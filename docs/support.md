# Support scripts for GraphBin

GraphBin provides a set of support scripts on [GitHub](https://github.com/Vini2/GraphBin/tree/master/support). Details of these support scripts are as follows.

## prep_result.py

`prep_result.py` is a support script that allows you to format an initial binning result in to the .csv format with contig identifiers and bin ID. Contigs are named according to their original identifier and bins are numbered according to the fasta file name. You can run `prep_result.py` as follows.

```
python prep_result.py --binned /path/to/folder_with_binning_result --output /path/to/output_folder
```
You can see the usage options of `prep_result.py` by typing `python prep_result.py -h` on the command line.

Formatted binning result will be stored in a file named `initial_contig_bins.csv` in the output folder provided. Bin IDs and corresponding fasta files for each bin will be recorded in a file named `bin_ids.csv` in the output folder provided.

## Before using MEGAHIT assemblies

After getting the MEGAHIT assemblies, please make sure to run `contig2fastg` from the MEGAHIT toolkit to build the assembly graph as follows.

```
megahit_toolkit contig2fastg 119 final.contigs.fa > final.graph.fastg
```

The MEGAHIT toolkit will result in a FASTG file which you can convert to GFA using [fastg2gfa](https://github.com/lh3/gfa1/blob/master/misc/fastg2gfa.c).

```
fastg2gfa final.fastg > final.gfa
```

If you want to run GraphBin on an assembly from a different `k` value output found in the MEGAHIT output folder `intermediate_contigs`, please make sure to build the `.fastg` file from the `.fa` file with the correct `k` value. For example, if you want to bin contigs from `k79.final.contigs.fa`, you should first build the corresponding `k79.final.graph.fastg` file and then run `fastg2gfa` as follows.

```
fastg2gfa 79.final.graph.fastg > 79.final.graph.gfa
```

## Before using Miniasm assemblies

Please note that, if you are using Miniasm assemblies, you should provide the edge sequences for the initial binning tool (not the contigs output from Miniasm). To get the edge sequences from the GFA file, you can use [`gfa2fasta.py` script](https://github.com/Vini2/GraphBin/blob/master/support/gfa2fasta.py) as the assembly graph consists of these edge sequences and not contigs.

If you come across a Python3 error when plotting graphs, please refer to [this thread](https://github.com/igraph/python-igraph/issues/88) to fix it.

## Visualising Assemblies

`visualise_result_SPAdes.py`, `visualise_result_SGA.py`, `visualise_result_MEGAHIT.py` and `visualise_result_Flye_Canu_Miniasm.py` allows you to visualize the binning result by denoting coloured contigs in the assembly graph according to their corresponding bins. You can visualise the initial binning result obtained from an existing binning tool and the final binning result obtained from GraphBin and compare.

You can see the usage options by typing `python visualise_result_SPAdes.py -h` or `python visualise_result_SGA.py -h` or `python visualise_result_MEGAHIT.py -h` or `python visualise_result_Flye_Canu_Miniasm.py -h` on the command line.

```
usage: visualise_result_SPAdes.py [-h] --initial INITIAL --final FINAL --graph
                                 GRAPH --paths PATHS --output OUTPUT
                                 [--prefix PREFIX] [--type TYPE]
                                 [--width WIDTH] [--height HEIGHT]
                                 [--vsize VSIZE] [--lsize LSIZE]
                                 [--margin MARGIN] [--dpi DPI]

optional arguments:
  -h, --help         show this help message and exit
  --initial INITIAL  path to the file containing the initial binning result
                     from an existing tool
  --final FINAL      path to the file containing the final GraphBin binning
                     result
  --graph GRAPH      path to the assembly graph file
  --paths PATHS      path to the contigs.paths file
  --output OUTPUT    path to the output folder
  --prefix PREFIX    prefix for the output image files
  --type TYPE        type of the image (jpg, png, eps, svg)
  --width WIDTH      width of the image in pixels
  --height HEIGHT    height of the image in pixels
  --vsize VSIZE      size of the vertices
  --lsize LSIZE      size of the vertex labels
  --margin MARGIN    margin of the figure
  --dpi DPI          dpi value

```
```
usage: visualise_result_SGA.py [-h] --initial INITIAL --final FINAL --graph
                              GRAPH --output OUTPUT [--prefix PREFIX]
                              [--type TYPE] [--width WIDTH] [--height HEIGHT]
                              [--vsize VSIZE] [--lsize LSIZE]
                              [--margin MARGIN] [--dpi DPI]

optional arguments:
  -h, --help         show this help message and exit
  --initial INITIAL  path to the file containing the initial binning result
                     from an existing tool
  --final FINAL      path to the file containing the final GraphBin binning
                     result
  --graph GRAPH      path to the assembly graph file
  --output OUTPUT    path to the output folder
  --prefix PREFIX    prefix for the output image files
  --type TYPE        type of the image (jpg, png, eps, svg)
  --width WIDTH      width of the image in pixels
  --height HEIGHT    height of the image in pixels
  --vsize VSIZE      size of the vertices
  --lsize LSIZE      size of the vertex labels
  --margin MARGIN    margin of the figure
  --dpi DPI          dpi value
```
```
usage: visualise_result_MEGAHIT.py [-h] --initial INITIAL --final FINAL --graph
                                  GRAPH --output OUTPUT [--prefix PREFIX]
                                  [--type TYPE] [--width WIDTH]
                                  [--height HEIGHT] [--vsize VSIZE]
                                  [--lsize LSIZE] [--margin MARGIN]
                                  [--dpi DPI]

optional arguments:
  -h, --help         show this help message and exit
  --initial INITIAL  path to the file containing the initial binning result
                     from an existing tool
  --final FINAL      path to the file containing the final GraphBin binning
                     result
  --graph GRAPH      path to the assembly graph file
  --output OUTPUT    path to the output folder
  --prefix PREFIX    prefix for the output image files
  --type TYPE        type of the image (jpg, png, eps, svg)
  --width WIDTH      width of the image in pixels
  --height HEIGHT    height of the image in pixels
  --vsize VSIZE      size of the vertices
  --lsize LSIZE      size of the vertex labels
  --margin MARGIN    margin of the figure
  --dpi DPI          dpi value
```

### Example visualisation

**Original MaxBin Labelling of the ESC+SPAdes dataset**

![](images/3G_initial_binning_result.png)

**Final Labelling of the of the ESC+SPAdes dataset produced from GraphBin**

![](images/3G_final_GraphBin_binning_result.png)
