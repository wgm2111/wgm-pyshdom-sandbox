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
DATA_FOLDER = "data"
BIN_FOLDER = os.path.abspath('bin')


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
        raise ValueError(msg.format{filename_out})

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
    if not (waveleng1 <= wavelen2):
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
                   distflag, alpha, gamma, nretanb, sretab, eretab, logre, maxradius, 
                   mietabfile, overwrite=False):
    "Make a mie table and write to the specified file."
    # Check arguments
    _check_filename(mietablefile)
    _check_value(wavelen1, wavelen2, deltawave, alpha, gamma, nretab, sretab, eretab, maxradius)
    _check_options(polarized_table, logre)
    _check_make_mie_flags(partype, aveflag, distflag)
    
    # Make a new name for the file
    file_already_exists = os.path.isfile(mietabfile)
    if (file_already_exists) and (overwrite is False):
        msg = "File, {0}, exists and overwrite = {1}."
        raise ValueError(msg.format(mietabfile, overwrite))

    # For non-aerosol particles
    if not (partype is "A"):
        inputs = ["T" if polarized_table else "F", 
                  "{0} {1}".format(str(wavelen1), 
                                   str(wavelen1 if (wavelen2 is None) else wavelen2)),
                  str(partype), str(avgflag), str(distflag), str(alpha), 
                  "{0} {1} {2}".format(nretanb, sretab, eretab), 
                  "T" if logre else "F", 
                  str(maxradius), self.mietabfile]
        
        protram_path = os.path.join()


            # Calculation
            # ============
            
            # get absolute path
            self.program_path = os.path.join('bin', "make_mie_table") # get absolute path

            # Call the program with an open pipe and quited output
            p = Popen(self.program_path, stdin=PIPE, stdout=open(os.devnull, 'wb'))
            
            # Pass use input via the open pipe
            p.communicate(self.make_mie_table_input)
            p.wait()                    # wait for the process to finish

            
            
        else:                   # if aerosol calculations are needed . . . 
            msg = "Add aerosol suport to {}."
            raise(NotImplementedError(msg.format(self.__class__)))



# Make a class for calling a mie table simulation
class MakeMieTable(object):
    """
    Class providing and OO interface to SHDOM's 'make_mie_table' routine.

    * INPUT PARAMETERS *
    POLTAB      logical flag for polarized scattering table output 
    .               (T=polarized, F=unpolarized)
    WAVELEN1    wavelength range (microns) for this band
    WAVELEN2    for monochromatic choose WAVELEN1=WAVELEN2
    PARTYPE     particle type: W=water, I=ice, A=aerosol
    .               if PARTTYPE='A' then the index of refraction is input,
    .               otherwise tables for water and ice index are used.
    AVGFLAG     'A' for spectral average over the wavelength range (for
    .               PARTYPE='W' or 'I'), 'C' to use the central wavelength.
    DELTAWAVE   wavelength interval for averaging (micron)
    RINDEX      aerosol complex index of refraction (negative imaginary part)
    PARDENS     aerosol particle bulk density (g/cm^3)
    DISTFLAG    'G' for gamma distribution, 'M' for modified gamma distribution,
    .               or 'L' for lognormal distribution
    ALPHA       distribution shape parameter (either alpha in gamma distribution
    .               or sigma in lognormal distribution).  Effective variance
    .               = 1/(alpha+3) for gamma, exp(alpha^2)-1 for lognormal.
    GAMMA       modified gamma distribution parameter
    NRETANB     number of effective radii entries in Mie table
    SRETAB      starting effective radius (micron) in Mie table
    ERETAB      ending effective radius (micron) in Mie table
    LOGRE       logical flag for log-spaced effective radius ouput 
    .              (T=log-spaced r_e, F=evenly spaced r_e)
    MAXRADIUS   maxium particle radius in size distribution (micron)
    MIETABFILE  output Mie scattering table file name

    * Ref * 
    Taken from propgen.txt (shdom)

    """

    # Data folder 
    data_folder = DATA_FOLDER

    # Construct a MakeMieTable object
    def __init__(self, poltab, wavelen1, wavelen2, 
                 partype, aveflag, deltawave, rindex, pardens, 
                 distflag, alpha, gamma, nretanb, sretab, eretab, logre, maxradius, 
                 mietabfile):
        "Construct a mie table with an OOP interface to the table data."
        
        # Check arguments
        _check_filename(mietabfile)
        _check_dist_flags(partype, aveflag, distflag)
        

        # Check incoming arguments and make a list
        self._check_args(mietabfile, polarized_table,
                         wavelen1, wavelen2, 
                         partype, avgflag, deltawa, rindex, 
                         pardens, distflag, alpha, gamma, 
                         nretanb,  sretab, eretab, 
                         logre, maxradius)


        # Make a new name for the file
        mietabfile = os.path.join(self.data_folder, mietabfile)
        self.mietabfile = mietabfile


        # prepare arguments
        if not(partype is "A"):      # Aerosols have refractive index
        
            # List of inputs
            make_mie_table_inputs = [
                "T" if polarized_table else "F", 
                "{0} {1}".format(str(wavelen1), 
                                 str(wavelen1 if (wavelen2 is None) else wavelen2)),
                str(partype), str(avgflag), str(distflag), str(alpha), 
                "{0} {1} {2}".format(nretanb, sretab, eretab), 
                "T" if logre else "F", 
                str(maxradius), self.mietabfile]

            # print('The input filename is, {}.'.format(make_mie_table_inputs))
            

            # Store the list of inputs
            self.make_mie_table_inputs = make_mie_table_inputs


            # String of inputs separated by lines
            self.make_mie_table_input = os.linesep.join(make_mie_table_inputs)
            

            # Calculation
            # ============
            
            # get absolute path
            self.program_path = os.path.join('bin', "make_mie_table") # get absolute path

            # Call the program with an open pipe and quited output
            p = Popen(self.program_path, stdin=PIPE, stdout=open(os.devnull, 'wb'))
            
            # Pass use input via the open pipe
            p.communicate(self.make_mie_table_input)
            p.wait()                    # wait for the process to finish

            
            
        else:                   # if aerosol calculations are needed . . . 
            msg = "Add aerosol suport to {}."
            raise(NotImplementedError(msg.format(self.__class__)))
                





    # Check the incoming arguments and make a list in the right order
    def _check_args(self, mietabfile, polarized_table,
                    wavelen1, wavelen2, 
                    partype, avgflag, deltawa, rindex, 
                    pardens, distflag, alpha, gamma, 
                    nretanb,  sretab, eretab, 
                    logre, maxradius):
        """
        Check incoming arguments and prepare for the call to make_mie_table.
        """
        
        # Filename
        if not (mietabfile.split('.')[-1] == "part"):
            msg = "The output file must end in .part, got {0}."
            raise ValueError(msg.format(
                mietabfile.split('.')[-1]))

        # polarized table
        if not (polarized_table in {True, False}):
            msg = "Argument polarized_table expected True or False, got {}."
            raise ValueError(msg.format(polarized_table))
            
        # Wavelengths
        wavelengths_are_ok = (wavelen1 > 0) 
        if not (wavelen2 is None):
            wavelengths_are_ok *= (wavelen2 > wavelen1)
        if not wavelengths_are_ok:
            msg = "Wavelengths unacceptable, min={0}, max={1}."
            raise ValueError(msg.format(wavelen1, wavelen2))
    
        # Particle type
        if not (partype in {"W", "I", "A"}):
            msg = "Supported partype 'W'ater, 'I'ce, and 'A'erosol, got {0}"
            raise ValueError(msg.format(partype))
          
        else:
            if (partype=="A") and (rindex is None):
                msg = "Specify 'rindex' for aerosol when using partype={0}"
                raise(ValueError(msg.format(partype)))
                
        # Averaging parameters
        if not(avgflag in {'A', 'C'}):
            msg = "Expected avgflag 'A'verage or 'C'enter, got {0}"
            raise(ValueError(msg.format(avgflag)))
            
        # Size distribution parameters
        if not (distflag in {'G', 'M', 'L'}):
            msg = "distflag can be (G)amma, (M)odified or (L)ognormal, got {}."
            raise(ValueError(msg.format(distflag)))
            
        # logre Use log spacing for the size distribution 
        if not (logre in {True, False}):
            msg = "Choose logre True for log space radius False for linear, got {}."
            raise(ValueError(msg.format(logre)))
            







if __name__ == "__main__":
    

    # Make a test_file.part
    test_file = 'test_file.part'
    mie_table  = MakeMieTable(test_file,             
                              polarized_table=True,
                              wavelen1 = 1.6, wavelen2 = None, 
                              partype = 'W', 
                              avgflag = 'C', 
                              distflag = "G",
                              alpha = 7, 
                              nretanb = 15, 
                              sretab = 0.1,
                              eretab = 1.0, 
                              logre = False, 
                              maxradius = 15.)
