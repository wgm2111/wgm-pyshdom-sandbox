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
import sys, os
from subprocess import Popen, PIPE, STDOUT
from sys_tools import PATH_ABS_SHDOM, get_abs_path, script_switch


# Define a folder to use for storing files 
DATA_FOLDER = "data_files"



# Input processing
# --


class MakeMieTable(object):
    """
    Class providing and OO interface to SHDOM's 'make_mie_table' routine.

    * INPUT PARAMETERS *
    POLTAB      logical flag for polarized scattering table output 
    -               (T=polarized, F=unpolarized)
    WAVELEN1    wavelength range (microns) for this band
    WAVELEN2    for monochromatic choose WAVELEN1=WAVELEN2
    PARTYPE     particle type: W=water, I=ice, A=aerosol
    -               if PARTTYPE='A' then the index of refraction is input,
    -               otherwise tables for water and ice index are used.
    AVGFLAG     'A' for spectral average over the wavelength range (for
    -               PARTYPE='W' or 'I'), 'C' to use the central wavelength.
    DELTAWAVE   wavelength interval for averaging (micron)
    RINDEX      aerosol complex index of refraction (negative imaginary part)
    PARDENS     aerosol particle bulk density (g/cm^3)
    DISTFLAG    'G' for gamma distribution, 'M' for modified gamma distribution,
    -               or 'L' for lognormal distribution
    ALPHA       distribution shape parameter (either alpha in gamma distribution
    -               or sigma in lognormal distribution).  Effective variance
    -               = 1/(alpha+3) for gamma, exp(alpha^2)-1 for lognormal.
    GAMMA       modified gamma distribution parameter
    NRETANB     number of effective radii entries in Mie table
    SRETAB      starting effective radius (micron) in Mie table
    ERETAB      ending effective radius (micron) in Mie table
    LOGRE       logical flag for log-spaced effective radius ouput 
    -              (T=log-spaced r_e, F=evenly spaced r_e)
    MAXRADIUS   maxium particle radius in size distribution (micron)
    MIETABFILE  output Mie scattering table file name

    """

    # Data folder 
    self.data_folder = DATA_FOLDER


    # Construct a MakeMieTable object
    def __init__(self, mietablefile,    
            polarized_table=None,
            wavelen1 = None, wavelen2 = None, 
            partype = None, 
            avgflag = None, 
            deltawa = None, 
            rindex = None, 
            pardens = None, 
            distflag = None, 
            alpha = None, 
            gamma = None,
            nretanb = None, 
            sretab = None,
            eretab = None, 
            logre = None, 
            maxradius = None):
        "Construct a mie table with an OOP interface to the table data."
        
        # Check incoming arguments and make a list
        self._check_args(mietablefile, polarized_table,
                         wavelen1, wavelen2, 
                         partype, avgflag, deltawa, rindex, 
                         pardens, distflag, alpha, gamma, 
                         nretanb,  sretab, eretab, 
                         logre, maxradius)


        # Make a new name for the file
        mietablefile = os.join(self.data_folder, mietablefile)
        self.mietablefile = mietablefile 


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
                str(maxradius), self.mietablefile]
            
            # Store the list of inputs
            self.make_mie_table_inputs = make_mie_table_inputs


            # String of inputs separated by lines
            make_mie_table_input = os.linesep.join(make_mie_table_inputs)


            # Calculation
            # ============
            
            # get absolute path
            program_path = get_abs_path("make_mie_table") # get absolute path

            # Call the program with an open pipe and quited output
            p = Popen(program_path, stdin=PIPE, stdout=open(os.devnull, 'wb'))
            
            # Pass use input via the open pipe
            p.communicate(make_mie_table_input)
            p.wait()                    # wait for the process to finish

            
            
        else:                   # if aerosol calculations are needed . . . 
            msg = "Add aerosol suport to {}."
            raise(NotImplementedError(msg.format(self.__class__)))
                





    # Check the incoming arguments and make a list in the right order
    def _check_args(self, mietablefile, polarized_table,
                    wavelen1, wavelen2, 
                    partype, avgflag, deltawa, rindex, 
                    pardens, distflag, alpha, gamma, 
                    nretanb,  sretab, eretab, 
                    logre, maxradius):
        """
        Check incoming arguments and prepare for the call to make_mie_table.
        """
        
        # Filename
        if not (mietablefile.split('.')[-1] == "part"):
            msg = "The output file must end in .part, got {0}."
            raise ValueError(msg.format(
                mietablefile.split('.')[-1]))

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
                              wavelen1 = .67, wavelen2 = None, 
                              partype = 'W', 
                              avgflag = 'C', 
                              distflag = "G",
                              alpha = 7, 
                              nretanb = 50, 
                              sretab = 0.5,
                              eretab = 15.0, 
                              logre = False, 
                              maxradius = 75.0)
