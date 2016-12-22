#!/usr/bin/env python3.5
#
# run_multiplex_models.py
#
# Script to apply a PyStan model across several multiplexed test/training sets.
# The training set is loaded in and the model parameters estimated. The
# resulting fitted model is then used to make predictions for the test set.
# The test set predicitons are then written to file.
#
# The script uses pysge to create command-lines that are submitted to the
# SGE scheduler. Each command-line runs the model on a single training/test
# set pair.
#
# (C) The James Hutton Institute 2016
# Author: Leighton Pritchard

"""
run_multiplex_models.py

A script to run a PyStan model script on multiplexed data
"""

import logging
import logging.handlers
import os
import sys
import time

from argparse import ArgumentParser

import pysge


class RunMultiplexException(Exception):
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
    parser.add_argument("--script", dest="scriptfile",
                        action="store", default="run_model.py", type=str,
                        help="Path to model script")
    parser.add_argument("-l", "--logfile", dest="logfile",
                        action="store", default=None, type=str,
                        help="Logfile location")
    parser.add_argument("--seed", dest="seedval",
                        action="store", default=None, type=int,
                        help="Seed for PRNG when running model")
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


# Main method for script
def main():
    """Main function for script."""
    args = parse_cmdline()

    # Set up logger
    logger = logger_setup(args.logfile, args.verbose)
    logger.info('# %s logfile', __file__)
    logger.info('# %s', time.asctime())
    logger.info(args)

    # Compose command-lines/joblist
    logger.info("Composing jobs for submission")
    if not os.path.isfile(args.scriptfile):
        msg = "Could not find script %s (exiting)"
        logger.error(msg, args.scriptfile)
        raise RunMultiplexException(msg % args.scriptfile)

    joblist = []
    try:
        for dirname in [dname for dname in os.listdir(args.indirname) if
                        dname.startswith('multiplex')]:
            cmd = "{0} -i {1} -l {2}".format(args.scriptfile,
                                             os.path.join(args.indirname,
                                                          dirname),
                                             os.path.join(args.indirname,
                                                          dirname, 'run.log'))
            if args.seedval is not None:
                cmd += " --seed {0}".format(args.seedval)
            joblist.append(pysge.Job(dirname,
                                     cmd))
    except:
        msg = "Could not create command-line list (exiting)"
        logger.error(msg)
        raise(RunMultiplexException(msg))
    logger.info("Jobs created:\n\t" + 
                "\n\t".join(["{0}: {1}".format(job.name, job.command) for
                             job in joblist]))

    # Run commands via pysge
    logger.info("Submitting jobs to SGE")
    pysge.build_and_submit_jobs(os.path.curdir, joblist)
    logger.info("Script exiting: %s", time.asctime())


# SCRIPT
if __name__ == '__main__':
    main()
