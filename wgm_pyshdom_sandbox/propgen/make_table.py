#!/usr/bin/python
# 
# Author: William GK Martin       
# Date: February 11 2015
# 


"""
# Module with an API to command line calls for property file generation

This Module contains class definitions for interacting with the following
parts of the SHDOM code:
(1) make_mie_table      [Running, not tested 6/15]
(2) make_tmatrix_table  [Not implemented 6/15]
(3) propgen             [Not implemented 6/15]  

"""

# Imports
# ------------------------------------------------------------------------
from __future__ import print_function

import warnings 
import os
from os.path import join, abspath, isfile
from subprocess import Popen, PIPE, STDOUT
# from sys_tools import PATH_ABS_SHDOM, get_abs_path, script_switch



# Define a folder to use for storing files
DIR = os.path.dirname(__file__)
DATA_FOLDER = os.path.join(DIR, "data")
BIN_FOLDER = os.path.join(DIR, "bin")

# Make mie table interface
# --
mie_flags = {'partype': ['W', 'I', 'A'], 
             'aveflag': ['A', 'C'], 
             'distflag': ['G', 'M', 'L']}

def _check_make_mie_flags(partype, aveflag, distflag, flags=mie_flags):
    "Check the input flags for a call to make mie table."
    # Checking input flags for make_mie table
    if not (partype in flags['partype']):
        msg = "Got partype = {0}.  Expected, \n\t{1}."
        raise ValueError(msg.format(partype, flags['partype']))
    if not (aveflag in flags['aveflag']):
        msg = "Got aveflag = {0}.  Expected, \n\t{1}."
        raise ValueError(msg.format(aveflag, flags['aveflag']))
    if not (distflag in flags['distflag']):
        msg = "Got distflag = {0}.  Expected, \n\t{1}."
        raise ValueError(msg.format(distflag, flags['distflag']))
    
def _check_filename(filename_out):
    "Force the output file to end in '.part'."
    if not (filename_out[-5:] == ".part"):
        msg = "Filename, '{0}' should end in '.part'."
        raise ValueError(msg.format(filename_out))

def _check_options(polarized_table, logre):
    "These options need to be true or false."
    # Polarized table is build when True is given
    if not (polarized_table in {True, False}):
        msg = "Argument polarized_table expected True or False, got {}."
        raise ValueError(msg.format(polarized_table))
    # Option to use log spacing for effective radii
    if not (logre in {True, False}):
        msg = "Choose logre True for log space radius False for linear, got {}."
        raise(ValueError(msg.format(logre)))
            
def _check_value(wavelen1, wavelen2, deltawave,
                 alpha, gamma, nretab, sretab, eretab, maxradius):
    "Check sanity of incoming values for Mie Scattering."
    if not (wavelen1 <= wavelen2):
        msg = "wavelen1 must be less than wavelen2.  Got wavelen1={1} and wavelen2={1} "
        raise ValueError(msg.format(wavelen1, wavelen2))
    if not deltawave >= 0:
        msg = "Wavelength increment for average must be positive. Got, deltawave={0}."
        raise ValueError(msg.format(deltawave))
    if not alpha > 0:
        msg = "Parameter alpha should be positive. Got, alpha={0}."
        raise ValueError(msg.format(alpha))
    if not alpha > 0:
        msg = "Parameter gamma should be positive. Got, gamma={0}."
        raise ValueError(msg.format(gamma))
    if not nretab >= 1:
        msg = "Number of effective radii must be positive. Got, nretab={0}."
        raise ValueError(msg.format(nretab))
    if not (sretab > 0 and sretab <= eretab):
        msg = "Starting and ending effective radii must be positive/ordered. Got, {0}."
        raise ValueError(msg.format((sretab, eretab)))
    if not (maxradius > eretab):
        msg = ("Size average must go larger than effective radius." + 
               " ertab={0} and maxradius={1}.")
        raise ValueError(msg.format(ertab, maxradius))


def make_mie_table(poltab, wavelen1, wavelen2, 
                   partype, aveflag, deltawave, rindex, pardens, 
                   distflag, alpha, gamma, nretab, sretab, eretab, logre, maxradius, 
                   mietabfile, overwrite=False):
    "Make a mie table and write to the specified file."
    # Check arguments
    _check_filename(mietabfile)
    _check_value(wavelen1, wavelen2, deltawave, alpha, gamma, nretab, sretab, eretab, maxradius)
    _check_options(poltab, logre)
    _check_make_mie_flags(partype, aveflag, distflag)
    
    # Make a new name for the file
    file_already_exists = os.path.isfile(mietabfile)
    if (file_already_exists) and (overwrite is False):
        msg = "File, {0}, exists and overwrite = {1}."
        raise ValueError(msg.format(mietabfile, overwrite))

    # For non-aerosol particles
    if not (partype is "A"):
        inputs = ["T" if poltab else "F", 
                  "{0} {1}".format(str(wavelen1), 
                                   str(wavelen1 if (wavelen2 is None) else wavelen2)),
                  str(partype), str(aveflag), str(distflag), str(alpha), 
                  "{0} {1} {2}".format(nretab, sretab, eretab), 
                  "T" if logre else "F", 
                  str(maxradius), mietabfile]
        make_mie_table_input = os.linesep.join(inputs)
        # protram_path = os.path.join()


        # Calculation
        # ============
        
        # get absolute path
        program_path = os.path.join(BIN_FOLDER, "make_mie_table")
        
        # Call the program with an open pipe and quited output
        p = Popen(program_path, stdin=PIPE, stdout=open(os.devnull, 'wb'))
        
        # Pass use input via the open pipe
        p.communicate(make_mie_table_input)
        p.wait()                    # wait for the process to finish
     
    
    else:                   # if aerosol calculations are needed . . . 
        msg = "Add aerosol suport to make_mie_table."
        raise(NotImplementedError(msg))
        

# Some convenience routines
def makeMieWaterGamma(fname, wavlen, reffmin, reffmax, reffnum, 
                      alpha=7.0,  # shape parameter estimate [Evans propgen.txt]
                      reff_cutoff=None, logspacing=True, polarized=True, 
                      overwrite=False,
                      make_mie_table=make_mie_table):
    """
    Call Frank Evan's `make_mie_table` program to generate the single-
    scattering properties of water droplets.  Mie scattering is for spherical 
    particles.  The output table is used with the function `propgen` to make a 
    (optical) property file for SHDOM.  

    Some arguments are hidden for simplicity.  Look carefully to make sure 
    the variables give the calculations you want.

    For example, we assume a gamma distribution of droplet number concentration 
    over droplet size (specified by particle radius).
      n(r) = a * r**alpha * exp(-b*r).
    
    Note that reff and alpha are the only parameters specified.  These are 
    related to effective variance by the following formula,
      reff = (alpha + 3) / b,
      veff = 1.0 / (3.0 + alpha).

    ## input
    fname   - filename prefix or full filename ending in .part
    wavelen - single wavelength used for the table
    reffmin - minimum effective radius available in the table
    reffmax - maximum effective radius available in the table
    reffnum - number of effective radii in between
    alpha=7 - shape parameter for Gamma distribution 

    ## output
    fname   - returns the filename of the table (may have a extra .part)

   """

    # Define fname to end in .part it doesn't already
    if not (fname[-5:] == ".part"):
        fname += ".part"

    # Integrate Reff up to 3 effective variances (if not given a cutoff)
    if reff_cutoff is None:
        effective_variance = 1.0 / (alpha + 3.0)
        reff_cutoff = reffmax + 3 * effective_variance

    # Make the call to make_mie_table
    make_mie_table(polarized, wavlen, wavlen, 'W', 'C', 0.0, 0.0, 1.0, 
                   'G', alpha, 0.1, reffnum, reffmin, reffmax, logspacing,
                   reff_cutoff, fname, overwrite=overwrite)

    # Return the filename
    return fname



