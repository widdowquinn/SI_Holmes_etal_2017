<img src="notebooks/images/JHI_STRAP_Web.png" style="width: 150px; float: right;">

# README.md - `SI_Holmes_etal_2017`

This repository contains files detailing the process of fitting the model of the enrichment array experiment described in [Holmes *et al.* (2017)](). These files are intended to enable independent reproduction, exploration and extension of the analysis reported in the paper.

To encourage exploration and reproduction, we have tried to make these notebooks compatible, so far as is possible, with [MyBinder](http://mybinder.org/), to enable you to run them in the cloud without having to install software on your own machine. To use these notebooks, click on [this link](http://mybinder.org:/repo/widdowquinn/si_holmes_etal_2017), or the button below.

[![Binder](http://mybinder.org/badge.svg)](http://mybinder.org:/repo/widdowquinn/si_holmes_etal_2017)

## Table of Contents

1. [Files and Directories](#files)
2. [Quick Start](#quickstart)
3. [Replicating the manuscript model](#replicate)


<a id="files"></a>
## Files and Directories

* `data/`: directory containing the raw microarray data, and the genomic data used in the analysis
* `models/`: directory containing Stan models in plan text format
* `multiplexing/`: directory containing scripts used to generate multiplexed data for *k*-fold cross-validation, and to fit the cross-validation models
* `notebooks/`: directory containing Jupyter notebooks describing and enabling reproduction of the data QA, model fitting and model validation
* `requirements.txt`: file describing the Python dependencies of the notebooks and scripts, which can be used to create a virtual environment for replication of the analysis from the paper.
* `LICENCE`: a copy of the MIT licence that governs the code contained in this repository
* `README.md`: this file

### How to get help for the code/analyses in this repository

Please raise any issues at the GitHub issues page for this repository:

* [GitHub issues page](https://github.com/widdowquinn/SI_Holmes_etal_2017/issues)

<a id="quickstart"></a>
## Quick Start

### Set up the environment

We would like our analysis to be reproducible, and for this we recommend using a Python virtual environment to ensure compatibility of dependencies and to replicate the environment used for the analysis. The virtual environment separates installation of Python packages from your system Python installation, enabling the running of these analyses without interfering with the system Python.

Using `pip` to install the required dependencies listed in `requirements.txt` should ensure that the code in this repository runs as expected.

#### Create and start the virtual environment

**NOTE:** You will need to have installed `virtualenv`[[*](http://docs.python-guide.org/en/latest/dev/virtualenvs/)] for your system.

```bash
virtualenv venv-SI_Holmes_2016 -p python3.6
source venv-SI_Holmes_2016/bin/activate
pip install -r requirements.txt
```

**Read more**

* `virtualenv`: [The Hitchhiker's Guide to Python](http://docs.python-guide.org/en/latest/dev/virtualenvs/)
* `pip`: [Installing Python Modules](https://docs.python.org/3/installing/)

### Use the notebooks

We have used the [Jupyter notebook](http://jupyter.org/) environment to facilitate [literate programming](https://en.wikipedia.org/wiki/Literate_programming), and to encourage exploration of and experimentation with the code and model. These notebooks have sections of explanatory text that are punctuated by code snippets. In the Jupyter notebook environment, all of these code snippets are editable and runnable.

#### Start the notebook environment

From the top level directory of this repository, start the Jupyter notebook server by issuing the command:

```bash
jupyter notebook
```

A new browser window or tab should open, containing the Jupyter homepage, which will show a listing of files and directories in the top level of the repository.

**Read more**

* `jupyter notebook`: [Quick-start guide](https://jupyter-notebook-beginner-guide.readthedocs.io/en/latest/)
* `jupyter notebook`: [Tutorial](https://www.datacamp.com/community/tutorials/tutorial-jupyter-notebook)

#### Opening the notebooks

From the Jupyter homepage in your browser window, click on the link to the `notebooks/` subdirectory. Then click on the link for the notebook you wish to use. The selected notebook should then open in a new browser tab.

When they were committed to the repository, the notebooks contained output from the original runs, so they can be read and understood without needing to rerun the models. If you would like to rerun/reproduce/modify these outputs, we recommend restarting the kernel and clearing all output before beginning. This can be done by clicking on `Kernel -> Restart & Clear Output`, in the notebook window.

<a id="replicate"></a>
## Replicating the manuscript model

To replicate the manuscript model from scratch: start the virtual environment, then run the notebooks and scripts as described below (remembering to use `Kernel -> Restart & Clear Output` in each notebook to remove the original/existing outputs):

### Data processing, QA and normalisation

* `01-data_qa.ipynb`: this will take the input data from the `notebooks/data/` directory, and process it into the output files:
  *  `notebooks/datasets/normalised_array_data.tab`: used for the full model fit, and to produce the multiplexed output datasets.
  *  `reduced_probe_data.tab`: a subset of `notebooks/datasets/normalised_array_data.tab`, used for testing code
  *  `reduced_locus_data.tab`: a subset of `notebooks/datasets/normalised_array_data.tab`, used for testing code


### Fitting the model on the full dataset

* `02-full_model_fit.ipynb`: this will fit the Stan model described in the notebook to the `notebooks/datasets/normalised_array_data.tab` processed data file, and conduct analyses to produce the final set of genes for which the estimated effect on enrichment due to passage (treatment) was positive, and render the figures used in the paper.

**NOTE:** the complete fit takes between 5 and 9 hours on my laptop (2013 MacBook Pro, 2.8GHz i7 16GB RAM).

### 10-fold crossvalidation

The crossvalidation dataset construction and model fit were conducted in the `multiplexing` directory, using Python scripts rather than Jupyter notebooks. To reproduce the dataset construction and fits, first change directory to `multiplexing`:

```bash
cd multiplexing
```

then build the input datasets with the `multiplex_data.py` script:

```bash
./multiplex_data.py -v -d ../notebooks/datasets/normalised_array_data.tab \
                    -k 10 -o 10-fold_CV --seed 123456789 \
                    -l 10-fold_CV_multiplex.log
```

This will create a new directory called `10-fold_CV`, containing one new subdirectory for each training/test split of the input dataset.

Next, use the `run_multiplex_models.py` script to fit the Stan model to each of the multiplexed training/test sets.

```bash
./run_multiplex_models.py -v -i 10-fold_CV --script ./run_model.py \
                          --seed 123456789 \
                          -l 10-fold_CV_run_models.log
```

**NOTE:** the `run_multiplex_models.py` has a dependency on the [`pysge` module](https://github.com/widdowquinn/pysge) for submission of jobs to our local cluster. This is not included in the `requirements.txt` file, so the script will fail at this point. The command-lines that this script produces in the log file can, however, be copied for execution on any system available to you. If you happen to be running on a cluster with SGE scheduling, then installation of `pysge` in the virtual environment will enable use of the cluster to fit the multiplexed models.

Finally, use the `join_multiplexed_data.py` script to combine prediction output from each of the 10 test sets into a single `.tab` file. This will contain predictions for each of the probes from the input dataset, using the model fit to the remaining training data.

```bash
./join_multiplexed_data.py -v -i 10-fold_CV -o 10-fold_CV.tab \
                           -l 10-fold_CV_join_data.log
```

The combined data produced in this way can then be used as input for the notebook `03-model_validation.ipynb`.

* `03-model_validation.ipynb`: this will conduct analyses on the combined output from 10-fold crossvalidation on the input dataset in `normalised_array_data.tab`. These analyses estimate the ability of the model to predict unseen 'output' array intensities by training it on 90% of the data at any one time, and testing it on the remaining 10% of the dataset.

### NOTE: PRNG seeds

All random processes in the model building and fitting can take a seed value for the pseudorandom number generator. For replication of the values in the paper, this seed should be set to `123456789` for all processes:

* the seed used for the main Stan fit (in notebook `02-full_model_fit`)
* the seed for splitting the input data into multiplexed sets (`multiplex_data.py`)
* the seed used for the multiplexing Stan fits (`run_multiplex_models.py`)
