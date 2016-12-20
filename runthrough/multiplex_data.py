#!/usr/bin/env python3.5
#
# multiplex_data.py
#
# Script to take input datasets for the Holmes et al. paper, and split them
# into several test and training sets for cross-validation, writing the
# datasets to several subdirectories.
#
# (C) The James Hutton Institute 2016
# Author: Leighton Pritchard

"""
multiplex_data.py

A script to multiplex input data for Stan models into several training and
test sets.
"""

import logging
import logging.handlers
import os
import random
import sys
import time

from argparse import ArgumentParser

import pandas as pd


class MultiplexException(Exception):
    """Exception raised when script fails"""
    def __init__(self, message):
        super().__init__(message)


# Parse command-line
def parse_cmdline():
    """ Parse command-line arguments
    """
    parser = ArgumentParser(prog=__file__)
    parser.add_argument("-o", "--outdir", dest="outdirname",
                        action="store", default='multiplexed_data', type=str,
                        help="Parent directory for output subfolders")
    parser.add_argument("-k", dest="kfold",
                        action="store", default=10, type=int,
                        help="Number of test/training datasets to create")
    parser.add_argument("-d", "--data", dest="datafile",
                        action="store", default=None, type=str,
                        help="Path to input data in tab-separated format")
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


def load_data(filename):
    """Return input data and unique locus tags.

    Assumes that data is tab-separated.
    """
    data = pd.read_csv(filename, sep="\t")
    return (data, data['locus_tag'].unique())


def chunks(iterable, chunksize):
    """Return the passed iterable in chunks of given size"""
    for idx in range(0, len(iterable), chunksize):
        yield iterable[idx:idx + chunksize]


def split_and_write_data(dataset, kfold, outdir, seedval):
    """Write k-fold CV datasets to subdirectories

    We insist that the parent directory is created anew, to avoid problems
    with overwriting some, but not all, previous k-fold CV datasets.
    """
    os.makedirs(outdir, exist_ok=False)

    # seed the PRNG?
    if seedval:
        random.seed(seedval)

    # Shuffle dataset indices and calculate k-fold step size
    indices = list(dataset.index)
    random.shuffle(indices)
    step = int(len(dataset)/kfold)

    # Shuffled index chunks of size step are the holdout sets.
    holdouts = chunks(indices, step)

    #
    for idx, test_indices in enumerate(holdouts):
        mplexdir = os.path.join(outdir, "multiplex{:05d}".format(idx))
        os.makedirs(mplexdir, exist_ok=False)
        dataset.iloc[test_indices].to_csv(os.path.join(mplexdir, "test.tab"),
                                          sep="\t",
                                          index=False)
        dataset.drop(test_indices).to_csv(os.path.join(mplexdir, "train.tab"),
                                          sep="\t",
                                          index=False)


# Main method for script
def main():
    """Main function for script."""
    args = parse_cmdline()

    # Set up logger
    logger = logger_setup(args.logfile, args.verbose)
    logger.info('# %s logfile', __file__)
    logger.info('# %s', time.asctime())
    logger.info(args)

    # Load input data
    try:
        logger.info("Loading data file %s", args.datafile)
        data, tag_ids = load_data(args.datafile)
        ntags = len(tag_ids)
        logger.info("Loaded data from %s", args.datafile)
        logger.info("Data shape: %d x %d", *data.shape)
        logger.info("Data contains %d unique locus tags", ntags)
        logger.info("Data headers: %s", list(data.columns))
    except:
        msg = "Could not load data file {0}".format(args.datafile)
        logger.error(msg)
        raise MultiplexException(msg)

    # Split data into k subsets, writing each to its own subdirectory
    if args.kfold < 2 or args.kfold > data.shape[0]:
        msg = "Value of k not valid (got %d, should be in range [2,%d])"
        logger.error(msg, args.kfold, data.shape[0])
        raise MultiplexException(msg % (args.kfold, data.shape[0]))
    try:
        logger.info("Writing %d k-fold CV datasets to subdirectories",
                    args.kfold)
        logger.info("Main subdirectory: %s", args.outdirname)
        split_and_write_data(data, args.kfold, args.outdirname, args.seedval)
        logger.info("New CV data subdirectories:")
        logger.info("%s", args.outdirname)
        for dirname in os.listdir(args.outdirname):
            logger.info("\t|-%s", dirname)
            for fname in os.listdir(os.path.join(args.outdirname, dirname)):
                logger.info("\t|\t|-%s", fname)
            logger.info("\t|")
    except:
        msg = "Writing CV datasets failed (exiting)"
        logger.error(msg)
        raise MultiplexException(msg)

    logger.info("Exiting cleanly %s", time.asctime())


# SCRIPT
if __name__ == '__main__':
    main()
