# Setting up GraphBin

## Dependencies

GraphBin installation requires python 3 (tested on Python 3.6 and 3.7). The following dependencies are required to run GraphBin and related support scripts.

* [python-igraph](https://igraph.org/python/) - version 0.7.1
* [biopython](https://biopython.org/) - version 1.74
* [cairocffi](https://pypi.org/project/cairocffi/)

## Downloading GraphBin

### Method 1: Conda Install

You can install GraphBin via [Conda](https://docs.conda.io/en/latest/). You can download [Anaconda](https://www.anaconda.com/distribution/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) which contains Conda.

Once you have installed Conda, you can install conda directly from the bioconda distribution using the command

```
conda install -c bioconda graphbin
```

Check if GraphBin is properly installed by typing `graphbin -h` on the command line. You should see the usage options as shown in section [Using GraphBin](https://github.com/Vini2/GraphBin#using-graphbin)

### Method 2: Setting up from GitHub repo
You can download the latest release of GraphBin from [Releases](https://github.com/Vini2/GraphBin/releases) or clone the GraphBin repository to your machine.

```
git clone https://github.com/Vini2/GraphBin.git
```

If you have downloaded a release, you will have to extract the files using the following command.

```
unzip [file_name].zip
```

Now go in to the GraphBin folder using the command

```
cd GraphBin/
```

Now we have to setup and install GraphBin.

## Setting up GraphBin

You can use Conda to setup an environemnt to run GraphBin **OR** you can use pip3 to install GraphBin.

### Option 1: Setup and environment using Conda

For this option, you need to have [Conda](https://docs.conda.io/en/latest/) installed on your machine. Then make sure you are in the GraphBin folder. Now run the following commands to create a Conda environment and activate it to run GraphBin.

```
conda env create -f environment.yml
conda activate graphbin
```

Add GraphBin to your PATH variable.
```
echo export PATH=$PATH:$(pwd) >> ~/.bashrc   # adding path to bashrc file to be available on login
export PATH=$PATH:$(pwd)                     # enabling the PATH for current terminal session
```

Now you are ready to run GraphBin.

If you want to switch back to your normal environment, run the following command.

```
conda deactivate
```

### Option 2: Install using pip3

You can install GraphBin globally or per user depending on your privileges to the system.

### Installing as admin
```
pip3 install .
```

### Installing for the active user
```
pip3 install . --user
```

**Note for Ubuntu users**

If you come across an error as `Failed building wheel for python-igraph` when installing GraphBin, you can install python-igraph as shown in [this thread](https://stackoverflow.com/questions/34962410/igraph-failed-to-install-through-pip).

Now let's prepare our results to run GraphBin.