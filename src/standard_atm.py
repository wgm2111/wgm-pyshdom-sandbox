

# 
# Python written by William G. K. Martin 
# Copywrite 2015
# GPL 2
# 
                                                                               
# ------------------------------------------------------------------        
# -------------     MCCLATCHY (1972) ATMOSPHERE DATA     -----------        
# ------------------------------------------------------------------        

"""
Standard atmosphere calculator based on Fortran source code and 
data from [MCCLATCHY 1972].  Contains routines for relating pressure, 
height, density and temperature for differnent lattitudes. 

"""


# Array operations
import numpy as np

# Source code import
from atm_mod import atmos




# Dictionary of acceptable lattitude keywords
# --
lats = {'tropical':1, 
       'midlat_summer':2, 'midlat_winter':3, 
       'subarc_summer':4, 'subarc_winter':5, 
       'standard':6, 
       'standard_rh0':7, 'standard_rh17':8, 'standard_rh87':9}

vcoords = {'pressure': 1, 
           'height': 2, 
           'density':3}    

def _argcheck(lat, vcoord):
    """
    Check the incoming arguments for compatibility with interface.  
    Print the acceptable options.
    """
    if not (lat in lats): 
        msg = "Got kwrd lat='{0}'. Expected one of the following: \n{1}."
        raise ValueError(msg.format(lat, sorted(lats.keys())))
    if not (vcoord in vcoords): 
        msg = "Got kwrd vcoord='{0}'. Expected one of the following: \n{1}."
        raise ValueError(msg.format(lat, sorted(vcoords.keys())))



# Define routines for computing armophere
# --
def get_atm(phd, lat='standard', vcoord='height'):
    """
    Routine for defining atmospheric data

    ## Input 
    phd    - float or vector (N,) of vertical coordinate values
    lat    - key for different latitude regions
    vcoord - specify the type of verticle coordinate 
             ('pressure', 'height' or 'density')
    
    ## Output
    atm_array - Array of shape (N, 9) or (9,).

    """
    
    # Check the arguments
    _argcheck(lat, vcoord)
    if np.isscalar(phd):
        phd = np.array([phd,])
    else: 
        if len(phd.shape) != 1:
            msg = "Expected phd to be scalar or 1-dimensional, phd.shape = {}."
            raise ValueError(msg.format(phd.shape))            

    # Define the output array 
    atm_dtype = np.dtype([('pressure', 'f8'), 
                          ('height', 'f8'), 
                          ('density', 'f8'), 
                          ('temperature', 'f8'), 
                          ('ozone_mixing_ratio', 'f8'), 
                          ('specific_humidity', 'f8'), 
                          ('saturation_ratio', 'f8'), 
                          ('ozone_above', 'f8'), 
                          ('water_vapor_above', 'f8')])
    atm_array = np.zeros((phd.size,), dtype=atm_dtype)
    
    # Define arguments for the call to atmos
    nphd = vcoords[vcoord]
    natm = lats[lat]

    # Call atmos
    atm_array[:] = [atmos(_phd, _phd, _phd, nphd, natm) for _phd in np.nditer(phd)]
    
    # for _phd, _atm_array in zip(np.nditer(phd),
    #                             np.nditer(atm_array)):

    #     # Call for each value of phd
    #     _atm_array[:] = atmos(phd, phd, phd, nphd, natm)
    
        
        
    return atm_array.squeeze()
    # if np.isphd
    # p, h, d, t, o, q, s, ocm, wcm, = atmos(p, h, d, nphd, natm)
    
    # return p, h, d, t, o, q, s, ocm, wcm
    


# Make a Atmosphere class
# --







if __name__ == "__main__":
    

    # Test call to atmos
    h = 0.1                     # hight in kilometers
    std_atm_h = get_atm(h)
    
    # Test call with array 
    hvec = np.linspace(0, 3, 50)
    
    std_atm_hvec = get_atm(hvec)
