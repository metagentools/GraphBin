# Setting up GraphBin

## Dependencies

GraphBin installation requires python 3 (tested on Python 3.6 and 3.7). The following dependencies are required to run GraphBin and related support scripts.

* [python-igraph](https://igraph.org/python/) - version 0.7.1
* [biopython](https://biopython.org/) - version 1.74
* [cairocffi](https://pypi.org/project/cairocffi/)

## Setting up GraphBin

### Method 1: conda install

You can install GraphBin via [Conda](https://docs.conda.io/en/latest/). You can download [Anaconda](https://www.anaconda.com/distribution/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) which contains Conda.

Now let's add conda channels so we know the locations where packages are stored.

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Once you have added the channels, you can install GraphBin directly from the bioconda distribution using the command

```
conda install -c bioconda graphbin
```

You can also create a new conda environment and install GraphBin from bioconda using the following commands.

```
# create conda environment
conda create -n graphbin

# activate conda environment
conda activate graphbin

# install graphbin
conda install -c bioconda graphbin
```


### Method 2: pip install

You can install GraphBin using pip.

```
pip install graphbin
```

After setup, check if GraphBin is properly installed by typing `graphbin -h` on the command line. You should see the usage options as shown in section [Using GraphBin](https://github.com/Vini2/GraphBin#using-graphbin)

Now let's prepare our results to run GraphBin.