#!/usr/bin/env python3.5
#
# run_model.py
#
# Script to apply a single PyStan model to a test/training sets.
# The training set is loaded in and the model parameters estimated. The
# resulting fitted model is then used to make predictions for the test set.
# The test set predicitons are then written to file.
#
#
# (C) The James Hutton Institute 2016
# Author: Leighton Pritchard

"""
run_model.py

A script to run a single PyStan model script on training/test data and
write predictions to a new output file
"""

import logging
import logging.handlers
import os
import pickle
import pystan
import random
import sys
import time

import pandas as pd

from argparse import ArgumentParser


# The fitting/prediction PyStan model
MODEL = """
data {
  int<lower=0> M;  # number of training datapoints
  int<lower=0> N;  # number of test datapoints
  int<lower=0> J;  # number of unique probes
  int<lower=1, upper=J> tagidx_train[M];  # probe indices (training)
  int<lower=1, upper=J> tagidx_test[N];   # probe indices (test)
  vector[M] t_train;
  vector[N] t_test;
  vector[M] x_train;  # training input datapoints
  vector[N] x_test;   # test input datapoints
  vector[M] y;
}
parameters {
  vector[J] a;
  vector[J] b;
  vector[J] g;
  vector[J] d;
  real mu_a;
  real mu_b;  
  real mu_g;
  real mu_d;  
  real<lower=0> sigma_y;
  real<lower=0,upper=100> sigma_a;
  real<lower=0,upper=100> sigma_b;
  real<lower=0,upper=100> sigma_g;  
  real<lower=0,upper=100> sigma_d;  
}
transformed parameters{
  vector[M] y_hat;
  vector[N] mu_pred;  

  for (i in 1:M)
    y_hat[i] = a[tagidx_train[i]] + b[tagidx_train[i]] * x_train[i] +
               g[tagidx_train[i]] * t_train[i] +
               d[tagidx_train[i]] * t_train[i] * x_train[i];
    
  for (j in 1:N)
    mu_pred[j] = a[tagidx_test[j]] + b[tagidx_test[j]] * x_test[j] +
                 g[tagidx_test[j]] * t_test[j] +
                 d[tagidx_test[j]] * t_test[j] * x_test[j];
}
model {
  sigma_a ~ uniform(0, 100);
  a ~ cauchy(mu_a, sigma_a);

  sigma_b ~ uniform(0, 100);
  b ~ cauchy(mu_b, sigma_b);

  sigma_g ~ uniform(0, 100);
  g ~ cauchy(mu_g, sigma_g);

  sigma_d ~ uniform(0, 100);
  d ~ cauchy(mu_d, sigma_d);

  y ~ normal(y_hat, sigma_y);
}
generated quantities {
  vector[N] y_pred;
  
  for (i in 1:N)
    y_pred[i] = normal_rng(mu_pred[i], sigma_y);
}
"""


class RunModelException(Exception):
    """Exception raised when script fails"""
    def __init__(self, message):
        super().__init__(message)


# Parse command-line
def parse_cmdline():
    """ Parse command-line arguments
    """
    parser = ArgumentParser(prog=__file__)
    parser.add_argument("-i", "--indir", dest="indirname",
                        action="store", default='multiplexed_data', type=str,
                        help="Parent directory for multiplexed subfolders")
    parser.add_argument("-l", "--logfile", dest="logfile",
                        action="store", default=None, type=str,
                        help="Logfile location")
    parser.add_argument("-v", "--verbose", dest="verbose",
                        action="store_true",
                        help="Give verbose output")
    parser.add_argument("--seed", dest="seedval",
                        action="store", default=None, type=int,
                        help="Seed random values for testing")
    return parser.parse_args()


# Set up logger
def logger_setup(logfilename=None, verbose=False):
    """Return logger for script"""
    logger = logging.getLogger(__file__)
    logger.setLevel(logging.DEBUG)

    # Add handler for STDERR
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)

    # Add logfile handler
    if logfilename:
        try:
            logstream = open(logfilename, 'w')
            err_handler_file = logging.StreamHandler(logstream)
            err_handler_file.setFormatter(err_formatter)
            err_handler_file.setLevel(logging.INFO)
            logger.addHandler(err_handler_file)
        except:
            msg = "Could not open {0} for logging".format(logfilename)
            logger.error(msg)
            raise RunModelException(msg)

    # Set loglevel
    if verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)

    return logger


# Main method for script
def main():
    """Main function for script."""
    args = parse_cmdline()

    # Set up logger
    logger = logger_setup(args.logfile, args.verbose)
    logger.info('# %s logfile', __file__)
    logger.info('# %s', time.asctime())
    logger.info(args)

    # Build filenames and load data
    traindatafile = os.path.join(args.indirname, "train.tab")
    testdatafile = os.path.join(args.indirname, "test.tab")
    preddatafile = os.path.join(args.indirname, "prediction.pkl")

    logger.info("Loading training data from %s", traindatafile)
    traindata = pd.read_csv(traindatafile, sep="\t")
    logger.info("Training data shape: %d x %d", *traindata.shape)
    logger.info("Training data headers: %s", traindata.columns)

    logger.info("Loading test data from %s", testdatafile)
    testdata = pd.read_csv(testdatafile, sep="\t")
    logger.info("Test data shape: %d x %d", *testdata.shape)
    logger.info("Test data headers: %s", testdata.columns)
    
    tag_ids = traindata['locus_tag'].unique()
    ntags = len(tag_ids)

    # Create data dictionary for model
    datadict = {'M': len(traindata),
                'N': len(testdata),
                'J': ntags,
                'tagidx_train': traindata['locus_tag_index'] + 1,
                'tagidx_test': testdata['locus_tag_index'] + 1,
                't_train': traindata['treatment'],
                't_test': testdata['treatment'],                        
                'x_train': traindata['log_input'],
                'x_test': testdata['log_input'],
                'y': traindata['log_output']}

    # Run the fit and prediction
    logger.info("Running Stan fit: %s", time.asctime())
    try:
        fit = pystan.stan(model_code=MODEL,
                          data=datadict,
                          iter=1000, chains=2,
                          seed=args.seedval)    
    except:
        msg = "Could not fit model (exiting)"
        logger.error(msg)
        raise RunModelException(msg)
    logger.info("Fit concluded: %s", time.asctime())

    # Extract the fit/prediction
    logger.info("Writing prediction/fit output to %s", preddatafile)
    unpermutedChains = fit.extract()
    unpermutedChains_df = pd.DataFrame([dict(unpermutedChains)])
    pickle.dump(unpermutedChains_df, open(preddatafile, 'wb'))
    logger.info("To open this file, use pickle.load(data, open(%s, 'rb))",
                preddatafile)
    logger.info("Script exiting: %s", time.asctime())


# SCRIPT
if __name__ == '__main__':
    main()
