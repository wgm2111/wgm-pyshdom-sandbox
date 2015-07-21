

"""
Script to make a pancake cloud test file.
"""


# import 
import numpy as np

# my imports
from atmos import get_atm



# Globals
FILEOUT = "pancake.part"

# Domain parameters
NX = NY = 64
XMAX = YMAX = 10.0
DX = XMAX / (NX - 1)
DY = YMAX / (NY - 1)

NZ = 32
ZMAX = 5.0
P = np.linspace(get_atm(0.0, vcoord='height')['pressure'], 
                get_atm(ZMAX, vcoord='height')['pressure'], NZ)
atm_array = get_atm(P, vcoord='pressure')
Z = atm_array['height']
T = atm_array['temperature']

# Define a cloud
CLOUD_CENTER = np.array([4.0, 4.0, 3.0])
CLOUD_DIMS = np.array([3.0, 3.0, .5])

def massdense(x, y, z, mass=1.0, cen=CLOUD_CENTER, dim=CLOUD_DIMS):
    """
    Define the mass density of liquid water in cloud droplets
    as a function of possition (x,y,z).
    """

    # Define the location of the cloud
    bool = ((abs(x - cen[0]) < dim[0]) *
            (abs(y - cen[1]) < dim[1]) * 
            (abs(z - cen[2]) < dim[2]))

    # Define the extinction where the cloud is
    if np.isscalar(x):
        out = bool * mass
    else:
        out = bool * mass * np.ones_like(x)

    # Return the float or float-array with mass density 
    return out

def cloudreff(x, y, z, reff=0.1, massdense=massdense):
    """
    Define the effective radius of liquid water droplets in cloud 
    as a function of possition (x,y,z).
    """
    return reff * massdense(x, y, z, mass=1.0)




# Script
if __name__ == "__main__":
    
    with open(FILEOUT, 'w') as fout:
        fout.write('{}\n'.format(3))
        fout.write('{} {} {}\n'.format(NX, NY, NZ))
        fout.write('{} {}\n'.format(DX, DY))
        fout.write('{}\n'.format(" ".join([str(z) for z in Z])))
        fout.write('{}\n'.format(" ".join([str(t) for t in T])))
        fout.flush()
        fout.close()

    print(open(FILEOUT,'r').read())
