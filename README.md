# Pandora Simulations

## Setup initial gene alignment data

The first thing that needs to be done before running the pipeline is to get and organise
the gene alignment data. This is done with a provided script.  

Run the following from the project root directory

```sh
bash scripts/setup_panx_data.sh
```

## Install requirements

### Snakemake

For a basic installation, run:

```sh
pip3 install --user snakemake
```

Other installation options can be found in the [snakemake docs](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Conda

A basic installation (Linux) that also prevents `conda` from being your default python installation
can be done with:

```sh
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda3.sh
conda_prefix=/opt/miniconda
bash miniconda3.sh -b -p "$conda_prefix"
# add conda to the end of your path
echo "export PATH=$PATH:$conda_prefix" >> "$HOME"/.bashrc
source "$HOME"/.bashrc
conda activate base
```

For more installation options, see the [conda docs](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html#).

## Run pipeline

From the project root directory, run

```sh
snakemake --use-conda
```
