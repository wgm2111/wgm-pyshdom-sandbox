

"""
Script to make a pancake cloud test file.

This script is a prototype.  It defines module variable and constructs
the output cloud properties by calling at the command line.  

** Next steps **
(1) Testing with SHDOM simulation
(2) Automate the construction. Make .part files with a function call. 


"""


# import 
import numpy as np
from numpy import ma

# my imports
from atmos import get_atm








def make_propfile(nscattab, scattabfiles, scatnums, poltabs, parfile,
                  maxnewphase, asymtol, fracphasetol, raywavelen, raysfcpres, 
                  nzo, zother, tempother, polout, propfile):
    """
    Make a property file from the given `.scat` and `.part` files. 
    """
    # Check arguments
    _check_filenames(scattabfiles, scatnums, poltabs, parfile, propfile)
    _check_extra_layers(nzo, zother, tmpother)
    _check_tolerances(maxnewphase, asymtol)
    

    msg = "Need a basic wrapper for calling 'propgen' with Python args."
    raise NotImplementedError(msg)
