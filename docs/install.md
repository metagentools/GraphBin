# Setting up GraphBin

## Dependencies

GraphBin installation requires python 3 (tested on Python 3.6 and 3.7). The following dependencies are required to run GraphBin and related support scripts.

* [python-igraph](https://igraph.org/python/) - version 0.7.1
* [biopython](https://biopython.org/) - version 1.74
* [cairocffi](https://pypi.org/project/cairocffi/)

## Setting up GraphBin

### Method 1: Conda Install

You can install GraphBin via [Conda](https://docs.conda.io/en/latest/). You can download [Anaconda](https://www.anaconda.com/distribution/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) which contains Conda.

Once you have installed Conda, you can install conda directly from the bioconda distribution using the command

```
conda install -c bioconda graphbin
```

You can also create a new conda environment and install GraphBin from bioconda using the following command and activate it.

```
conda create -n graphbin -c bioconda graphbin
conda activate graphbin
```


### Method 2: Setting up from GitHub repo
You can clone the GraphBin repository to your machine for development purposes.

```
git clone https://github.com/Vini2/GraphBin.git
```

Now go in to the GraphBin folder using the command

```
cd GraphBin/
```

We can install GraphBin via [flit](https://pypi.org/project/flit/).

```
# install flit
pip install flit

# install graphbin via flit
flit install -s --python `which python`
```

After setup, check if GraphBin is properly installed by typing `graphbin -h` on the command line. You should see the usage options as shown in section [Using GraphBin](https://github.com/Vini2/GraphBin#using-graphbin)

Now let's prepare our results to run GraphBin.