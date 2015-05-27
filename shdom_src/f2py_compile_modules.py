
# 
# Author: William GK Martin
# 




import numpy.f2py as f2py


# Open and compile miewig
with open('miewig.f', 'r') as file_miewig:
    f2py.compile(file_miewig.read(), modulename="py_miewig")

# # Open and compile shdom.f90
# with open('shdom.f90', 'r') as file_shdom:
#     f2py.compile(file_shdom.read(), modulename="py_shdom")
    
# # ocean brdf
# with open('ocean_brdf.f', 'r') as file_ocean_brdf:
#     f2py.compile(file_ocean_brdf.read(), modulename="py_ocean_brdf")

# with open('shdomsub2.f', 'r') as file_shdomsub2:
#     f2py.compile(file_shdomsub2.read(), modulename='py_shdomsub2')





if __name__ == "__main__":
    

    # Test mie wig
    import py_miewig
    # import py_shdom
    # import py_ocean_brdf
    # import py_shdomsub2
    
