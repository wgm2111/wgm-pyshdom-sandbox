#!/usr/bin/env python

# Author: William G.K. Martin (wgm2111@cu where cu=columbia.edu)
# copyright (c) 2011
# liscence: BSD style


# future
from __future__ import division, absolute_import, print_function


# Import the numpy version of distutils (which has support for f2py)
import os
from os.path import join
import sys

from setuptools.command.install import install
from setuptools.command.develop import develop
# from distutils.command.develop import develop
# from numpy.distutils.command.develop import develop
from numpy.distutils.core import Extension

# from numpy.distutils.command.install import install as DistutilsInstall
# from numpy.distutils.command.develop import develop as DistutilsDevelop

import platform
import subprocess as sbp
import tempfile
import stat
import shlex
import shutil


# Global names
# ================================================================

NAME = "pyshdom"
VERSION = '0.0.1.dev1'
DESCRIPTION = (
    "Pre-development work on Python wrappers for SHDOM.")
AUTHOR = "William G.K. Martin"
AUTHOR_EMAIL = "wgm2111@columbia.edu"
LICENSE = "GPL2"



# Targe location for developement
# --
PYTHON_SOURCE_DIR = 'wgm_pyshdom_sandbox' # python bindings
# TARGET_SO_LIB = join(PYTHON_SOURCE_DIR, 'lib')  # .so files
# TARGET_BIN = join(PYTHON_SOURCE_DIR, 'bin')     # compiled binaries

# Souces locations for Fortran codes
SHDOM_ORIG_SOURCE_DIR = 'shdom'    # Unpacking Frank's SHDOM
SHDOM_MOD_SOURCE_DIR = 'shdom_src' # Modified SHDOM source files
FORTRAN_SOURCE_DIR = 'src'         # Other Fortran sources 

# List of propgen scripts to make
TARGET_PROPGEN_DIR = join(PYTHON_SOURCE_DIR, 'propgen')
TARGET_PROPGEN_BIN = join(PYTHON_SOURCE_DIR, 'propgen','bin')
SHDOM_PROPGEN_EXECS = [
    'propgen', 'make_mie_table', 'make_tmatrix_table', 'plotscattab']


# Define the structure of the compiled package
PACKAGES = [PYTHON_SOURCE_DIR, 
            PYTHON_SOURCE_DIR+'.propgen', 
            PYTHON_SOURCE_DIR+'.shdom']

PACKAGE_DIR = None


# Begin the f2py stuf for building the package
# ================================================================

# Define the fortran extensions to make
ext_modules = [
    Extension(TARGET_PROPGEN_DIR+'.'+'atm_mod', 
              sources = [join(FORTRAN_SOURCE_DIR, 'atm_mod.f'),]),]




# Routines for making SHDOM programs for property file generation
# =====================================================================

def make_shdom(src_path=SHDOM_MOD_SOURCE_DIR, 
               propgen_bin_dir=TARGET_PROPGEN_BIN, 
               shdom_execs=SHDOM_PROPGEN_EXECS):
    """
    Run make the FORTRAN programs that will be called directly. This is 
    mainly for 'propgen' and related routines. 
    """
    
    # Get the absolute path to 
    shdom_dir = os.path.abspath(src_path)
    propgen_bin_dir = os.path.abspath(propgen_bin_dir)

    # Change directory to the current working directory
    current_path = os.getcwd()
    os.chdir(shdom_dir)

    # Windows needs executables
    ext = '.exe' if (platform.system() == 'Windows') else ''
        
    # Make routines for property generation 
    for name in shdom_execs:#('propgen', 'make_mie_table'):
        os.system('make {exe}'.format(exe=name))
        src_path = os.path.join(shdom_dir, name+ext)
        dst_path = os.path.join(propgen_bin_dir, name+ext)
        shutil.copyfile(src_path, dst_path)
       
        # Change permissions to make the programs executable
        if not platform.system() == 'Windows':
            os.chmod(dst_path, 
                     stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | stat.S_IWUSR)

    # Return to original path
    os.chdir(current_path)


class SrcInstall(install):#(DistutilsInstall):
    """
    Class which is passed to setup to define what actions
    are perfomed on call to 

      $ python setup.py install
    
    """
    def run(self):
        install.run(self)
        make_shdom()#F2PY_SRC_PATH, self.script_dir)


class SrcDevelop(develop):#DistutilsDevelop):
    """
    Class which is passed to setup to define what actions
    are perfomed on call to 

      $ python setup.py develop
    
    """
    def run(self):
        'running make'

        develop.run(self)
        make_shdom()#F2PY_SRC_PATH, self.script_dir)


# Run setup when called as a program
# --
if __name__ == "__main__":
    
    # import setup
    # from setuptools import setup
    from numpy.distutils.core import setup
    # from numpy.distutils.core import setup

    # Call to setup the package
    setup(name = NAME,
          version = VERSION,
          description = DESCRIPTION,
          author = AUTHOR,
          author_email = AUTHOR_EMAIL,
          packages=PACKAGES, 
          package_dir=PACKAGE_DIR,
          license = LICENSE, 
          ext_modules=ext_modules, 
          cmdclass={
              'install':SrcInstall,
              'develop':SrcDevelop})


