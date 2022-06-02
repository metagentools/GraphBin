# Contributing to GraphBin project

We love to have your contributions to the GraphBin, whether it's:
* Reporting a bug
* Discussing the current state of the code
* Submitting a fix
* Proposing new features

## Clone and install GraphBin onto your machine

On GitHub, [fork](https://docs.github.com/en/get-started/quickstart/fork-a-repo) the GraphBin repository and clone it to your machine.

```
# clone repository to your local machine
git clone https://github.com/metagentools/GraphBin.git
```

Move to the GraphBin directory, install via [flit](https://pypi.org/project/flit/).

```
# go to repo direcotry
cd GraphBin

# install flit
pip install flit

# install graphbin via flit
flit install -s --python `which python`
```

## Test GraphBin installation

Run the following command and the all the tests should pass.

```
pytest
```

## Coding Style

We adhere to the [PEP 8](https://peps.python.org/pep-0008/) style guide. 

Before committing, run [`black`](https://pypi.org/project/black/) and [`isort`](https://pypi.org/project/isort/) before committing.

## Report bugs using Github's issues

We use GitHub issues to track public bugs. Report a bug by opening a new issue in GitHub [issues](https://github.com/metagentools/GraphBin/issues). You will get to select between templates for bug report and feature request.

## Committing code

Once you have finished coding and all the tests pass, commit your code and make a pull request. Make sure to follow the commit style of [c3dev](https://github.com/cogent3/c3dev/wiki#style-for-commit-messages).

```
git commit -m "<commit message>"
git push
```

Your contribution will be reviewed before accepting it. 

## License

By contributing, you agree that your contributions will be licensed under the BSD-3 License.

## References

This document was adapted from the open-source contribution guidelines for [Transcriptase](https://github.com/briandk/transcriptase-atom/blob/master/CONTRIBUTING.md) and [c3dev](https://github.com/cogent3/c3dev/wiki/How-to-Contribute-Code).