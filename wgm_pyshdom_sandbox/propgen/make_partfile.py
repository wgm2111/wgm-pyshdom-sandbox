

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




# Globals
FILEOUT = "pancake.part"
DECIMALS = 5
NCOMP = 1
MASSTOL = 1e-5

# Domain parameters
NX = NY = 64                    # xy - horizontal parameters
XMAX = YMAX = 7.0
X = np.linspace(0, XMAX, NX)
Y = np.linspace(0, YMAX, NY)
DX = X[1] - X[0]
DY = Y[1] - Y[0]                
NZ = 32                         # z - verticle parameters
ZMAX = 5.0
P = np.linspace(get_atm(0.0, vcoord='height')['pressure'], 
                get_atm(ZMAX, vcoord='height')['pressure'], NZ)
atm_array = get_atm(P, vcoord='pressure')
Z = np.round(atm_array['height'], DECIMALS)
T = np.round(atm_array['temperature'], DECIMALS)

# Define a cloud
CLOUD_CENTER = np.array([4.0, 4.0, 3.0])
CLOUD_DIMS = np.array([2.0, 2.0, .5])

def cloudmass(x, y, z, mass=1.0, cen=CLOUD_CENTER, dim=CLOUD_DIMS):
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

def cloudreff(x, y, z, reff=0.1, cloudmass=cloudmass):
    """
    Define the effective radius of liquid water droplets in cloud 
    as a function of possition (x,y,z).
    """
    return reff * cloudmass(x, y, z, mass=1.0)




# Script
if __name__ == "__main__":
    



    with open(FILEOUT, 'w') as fout:

        # Writing the header
        fout.write('{}   propgen particle file\n'.format(3))
        fout.write('{} {} {}\n'.format(NX, NY, NZ))
        fout.write('{} {}\n'.format(DX, DY))
        fout.write('{}\n'.format(
            " ".join([("{:."+str(DECIMALS)+"f}").format(z) for z in Z])))
        fout.write('{}\n'.format(
            " ".join([("{:."+str(DECIMALS)+"f}").format(t) for t in T])))
        fout.flush()

        # Writing out particle positions in sparse format
        ix = np.arange(X.size) + 1
        iy = np.arange(Y.size) + 1
        iz = np.arange(Z.size) + 1
        ixx, iyy, izz = np.array(np.meshgrid(ix, iy, iz))
        xx, yy, zz = np.array(np.meshgrid(X, Y, Z))
        
        mass = cloudmass(xx, yy, zz)
        reff = cloudreff(xx, yy, zz)
        pointdata_list = [(_ix, _iy, _iz, _x, _y, _z, _mass, _reff) 
                          for _ix, _iy, _iz, _x, _y, _z, _mass, _reff
                          in zip(
                              ixx.flatten(), iyy.flatten(), izz.flatten(),
                              xx.flatten(), yy.flatten(), zz.flatten(),
                              mass.flatten(), reff.flatten()) if _mass >= (MASSTOL)]
        
        # Loop over the points with liquid mass
        for _ix, _iy, _iz, _x, _y, _z, _mass, _reff in pointdata_list:
            temp = ("{ix:3d} {iy:3d} {iz:3d} {ncomp:2d} " + 
                    "2 {mass:.8f} {reff:.8f}\n")
            line = temp.format(
                ix=_ix, iy=_iy, iz=_iz, ncomp=NCOMP, mass=_mass, reff=_reff)

            fout.write(line)
            
            if _iz == 1: 
                fout.flush()

           # end = "2 {mass:.8f} {reff:.8f}".format(
           #     mass=cloudmass(_x, _y, _z), reff=cloudreff(_x, _y, _z))


           # start = ('{ix:3s} {iy:3s} {iz:3s} {ncomp:2s} ').format(


        # Finishing up
        fout.flush()
        fout.close()

    print(open(FILEOUT,'r').read())
