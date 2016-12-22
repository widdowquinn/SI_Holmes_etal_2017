# README.md

This repository contains files detailing the process of fitting the model of the enrichment array experiment described in [Holmes *et al.* (2017)](). These files are intended to enable independent reproduction and extension of the analysis reported in the paper.

## Files and Directories

* `data/`: directory containing the raw microarray data, and the genomic data used in the analysis
* `notebooks/`: directory containing Jupyter notebooks describing the data QA, model fitting and model validation
* `multiplexing/`: directory containing scripts used to generate multiplexed data for 10-fold cross-validation, and to fit the cross-validation models
* `LICENCE`: a copy of the licence that governs the code and data contained in this repository
* `README.md`: this file
* `requirements.txt`: file describing the Python dependencies of the notebooks and scripts, used to create a virtual environment for replication of the analysis, and experimentation.

## Quick Start

### Set up the environment

We would like our analysis to be reproducible, so recommend using a Python virtual environment to ensure compatibility of dependencies. The virtual environment separates the running of these analyses from your system Python installation and enables installation of Python packages, without interfering with the system Python.

Using `pip` to install the required dependencies listed in `requirements.txt` should ensure that the code in this repository runs as expected.

#### Create and start the virtual environment

```bash
virtualenv venv-SI_Holmes_2016
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


## Replicating the manuscript model

To replicate the manuscript model from scratch: start the virtual environment, then run the notebooks and scripts as described below (remembering to use `Kernel -> Restart & Clear Output` in each notebook to remove the original/existing outputs):

### Data processing, QA and normalisation

* `01-data_qa.ipynb`: this will take the input data from the `notebooks/data/` directory, and process it into the `notebooks/datasets/normalised_array_data.tab` file that was used for the full model fit, and to produce the multiplexed output datasets.

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
./multiplex_data.py -v -d ../notebooks/datasets/normalised_array_data.tab -k 10 -o 10-fold_CV --seed 123456789 -l 10-fold_CV_multiplex.log
```

This will create a new directory called `10-fold_CV`, containing one new subdirectory for each training/test split of the input dataset.

Next, use the `run_multiplex_models.py` script to fit the Stan model to each of the multiplexed training/test sets.

```bash
./run_multiplex_models.py -v -i 10-fold_CV --script ./run_model.py --seed 123456789 -l 10-fold_CV_run_models.log
```

**NOTE:** the `run_multiplex_models.py` has a dependency on a `pysge` module for submission of jobs to our local cluster. This is not included in the `requirements.txt` file, so the script will fail at this point. The command-lines that this script produces can, however, be executed on any system available to you.


### NOTE: PRNG seeds

All random processes in the model building/fitting can take a seed value for the pseudorandom number generator. For replication of the values in the paper, this seed should be set to `123456789` for all processes:

* the seed used for the main Stan fit (in notebook `02-full_model_fit`)
* the seed for splitting the input data into multiplexed sets (`multiplex_data.py`)
* the seed used for the multiplexing Stan fits (`run_multiplex_models.py`)
