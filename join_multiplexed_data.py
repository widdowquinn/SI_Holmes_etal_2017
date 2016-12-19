#!/usr/bin/env python3.5
#
# join_multiplexed_data.py
#
# Script to read pickled output from multiplexed runs generated by
# run_multiplex_models.py, and produce a single output .tab file
# containing predictions, prediction credibility intervals, and the
# observed output for the probes.
#
# (C) The James Hutton Institute 2016
# Author: Leighton Pritchard

"""
join_multiplexed_data.py

A script to join the pickled PyStan output from multiplexed data into a
single .tab file
"""


import logging
import logging.handlers
import os
import pickle
import sys
import time

from argparse import ArgumentParser

import numpy as np
import pandas as pd


class JoinMultiplexException(Exception):
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
    parser.add_argument("-o", "--outfile", dest="outfilename",
                        action="store", default="multiplexed_predictions.tab",
                        type=str,
                        help="Path to model script")
    parser.add_argument("-l", "--logfile", dest="logfile",
                        action="store", default=None, type=str,
                        help="Logfile location")
    parser.add_argument("-v", "--verbose", dest="verbose",
                        action="store_true",
                        help="Give verbose output")
    return parser.parse_args()


# Set up logger
def logger_setup(logfilename=None, verbose=False):
    """Return logger for script"""
    logger = logging.getLogger(__file__)
    logger.setLevel(logging.DEBUG)

    # Add handler for STDERR
    err_handler = logging.StreamHandler(sys.stderr)
    logger.addHandler(err_handler)
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
            raise MultiplexException(msg)

    # Set loglevel
    if verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)

    return logger


def get_pickle_filenames(dirname):
    """Returns a list of paths to the pickled output files"""
    pickles = []
    for subdir in [dname for dname in os.listdir(dirname) if
                   os.path.isdir(dname)]:
        for fname in [ifn for ifn in
                      os.listdir(os.path.join(dirname, subdir)) if
                      os.path.splitext(ifn)[-1] == '.pkl']:
            pickles.append(os.path.join(dirname, subdir, fname))
    return pickles


def extract_variable_summaries(df, varname):
    """Extracts summary variables for a variable in the passed dataframe

    Calculates mean, std, median, and 5%, 25%, 75% and 95% percentiles
    for the passed variable, returning them as a dataframe.
    """
    # Using Pandas methods
    mean = pd.Series(data[varname][0].mean(0))
    se = pd.Series(data[varname][0].std(0))

    # Need to use numpy functions
    median = pd.Series(np.median(prediction_fit[varname][0], 0))
    perc_5 = pd.Series(np.percentile(prediction_fit[varname][0], 5, 0))
    perc_25 = pd.Series(np.percentile(prediction_fit[varname][0], 25, 0))
    perc_75 = pd.Series(np.percentile(prediction_fit[varname][0], 75, 0))
    perc_95 = pd.Series(np.percentile(prediction_fit[varname][0], 95, 0))



def load_pickle_data(datafilename, picklename, df=None):
    """Returns a Pandas dataframe containing summary of pickled data

    datafilename takes the path to input data corresponding to the pickled
    results
    picklename takes the path to the pickled results file
    If a dataframe is passed in df, data from the pickled file is
    appended to the dataframe.
    """
    # load pickled results and input data
    with open(picklename, 'rb') as ifh:
        data = pickle.load(ifh)
    indata = pd.read_csv(datafilename, sep="\t")

    # Obtain summary of pickled data predictions and join to the input data
    y_pred = extract_variable_summaries(df, 'y_pred')  # Predictions
    outdata = indata.reset_index(drop=True)
    outdata = outdata.join(y_pred, how="outer")

    # Append this dataframe to the passed df
    if df is not None:
        outdata = df.append(outdata)

    return outdata


# Main method for script
def main():
    """Main function for script."""
    args = parse_cmdline()

    # Set up logger
    logger = logger_setup(args.logfile, args.verbose)
    logger.info('# %s logfile', __file__)
    logger.info('# %s', time.asctime())
    logger.info(args)

    # Identify pickle files
    logger.info("Looking for test data/pickle files in %s", args.indirname)
    infiles = get_result_filenames(args.indirname)
    logger.info("Identified input files:\n%s", "\n\t".join(infiles))

    # Load each pickle file, and add results to a growing Pandas dataframe
    logger.info("Loading pickle files into dataframe")
    df = None
    for datafile, results in pickles:
        logger.info("Loading test data from %s, results from %s",
                    datafile, results)
        df = load_pickle_data(datafile, results, df)
    loger.info("Loaded pickled data. Dataframe size: %d x %d", *dfshape)



# SCRIPT
if __name__ == '__main__':
    main()
