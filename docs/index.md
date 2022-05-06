![](images/GraphBin_logo.png)

**GraphBin** is a NGS data-based metagenomic contig bin refinment tool that makes use of the contig connectivity information from the assembly graph to bin contigs. It utilizes the binning result of an existing binning tool and a label propagation algorithm to correct mis-binned contigs and predict the labels of contigs which are discarded due to short length.

**NEW:** The conda package of GraphBin is now available on bioconda at [https://anaconda.org/bioconda/graphbin](https://anaconda.org/bioconda/graphbin).

**Note:** Due to recent requests from the community, we have added support for long-read assemblies produced from Flye, Canu and Miniasm. Please note that GraphBin has not been tested extensively on long-read assemblies. We originally developed GraphBin for short-read assemblies. Long-read assemblies might have sparsely connected graphs which can make the label propagation process less effective and may not result in improvements.

