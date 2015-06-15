#!/usr/bin/python
# 
# Author: William GK Martin       
# Date: February 11 2015
# 


"""
# Module with an API to command line calls for property file generation

This Module contains class definitions for interacting with the following
parts of the SHDOM code:
(1) make_mie_table
(2) make_tmatrix_table
(3) propgen

"""

# Imports
# ------------------------------------------------------------------------
from __future__ import print_function

import warnings 
import sys, os
from subprocess import Popen, PIPE, STDOUT
from sys_tools import PATH_ABS_SHDOM, get_abs_path, script_switch


# Define a folder to use for storing files 
DATA_FOLDER = "mono_les_files"



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

    # Define a dictionary of acceptable arguments
    ok_values = {'poltab': [True, False], 
                 'partype': ['W', 'I', 'A'], 
                 'avgflag': ['A', 'C'], 
                 'distflag': ['G', 'M', 'L'], 
                 'logre': [True, False]}
                 

    # Define bounds
    


        
    def __init__(
            self, 
            mietablefile,    
            polarized_table=None,
            wavelen1 = None, 
            wavelen2 = None, 
            partype = None, 
            avgflag = None, 
            deltawave = None, 
            rindex = None, 
            pardens = None, 
            alpha = None, 
            gamma = None
            nretanb = None, 
            sretab = None
            eretab = None, 
            logre = None, 
            maxradius = None)
        "Construct a mie table with an OOP interface to the table data."
        
        # Check incoming arguments and make a list
        argument_list = self._get_arg_list(mietablefile, polarized_table,
                                           wavelen1, wavelen2, 
                                           partype, avgflag, deltawa, rindex, 
                                           pardens, alpha, gamma, 
                                           nretanb,  sretab, eretab, 
                                           logre, maxradius)

        

    # Check the incoming arguments and make a list in the right order
    def _get_arg_list(self, mietablefile, polarized_table,
                      wavelen1, wavelen2, 
                      partype, avgflag, deltawa, rindex, 
                      pardens, alpha, gamma, 
                      nretanb,  sretab, eretab, 
                      logre, maxradius):
        "Check incoming arguments and prepare for the call to make_mie_table."
        
        # polarized table
        if not (polarized_table in {True, False}):
            msg = "Argument polarized_table expected True or False, got {}."
            raise ValueError(msg.format(polarized_table))
            
        # Wavelengths
        if (wavelen2 < wavelen1) or (min(wavelen2, wavelen1)<0): 
            msg = "Wavelengths unacceptable, min={0}, max={1}."
            raise ValueError(msg.format(wavelen1, wavelen2))
        if not (partype in {""})
    
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
            
        # 

        arglist = [{True:'T', False:'F'}[polarized_table], 
                   str(wavelen1), str(wavelen2), ]
        # Make inputs to pass to the propgen command line tool
        self.make_mie_table_inputs = [
            polarized_table, "{0} {0}".format(wavelength), 
            partype, avgflag, distflag, alpha, 
            "{0} {1} {2}".format(Nretab, Sretab, Eretab), logre, maxradius, scattable]
        
        # String of inputs separated by lines
        self.make_mie_table_input = os.linesep.join(self.make_mie_table_inputs)



        # Define the incoming parameters
        return arg_list

        
        # 
        

        

        # Define the incoming parameters
        (scattable, polarized_table, 
         distflag, alpha, 
         Nretab, Sretab, Eretab, logre, maxradius,
         partype, avgflag) = args[:]
        
        
