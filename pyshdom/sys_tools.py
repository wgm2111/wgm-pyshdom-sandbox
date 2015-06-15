#!/usr/bin/python
# 
# Author: William GK Martin       
# Date: February 11 2015
# 


"""
System tools for often used routines needed to interact with 
SHDOM from python.  

*DATA*
PATH_REL_SHDOM
PATH_ABS_SHDOM



This module defines the path to the source code of SHDOM.  Later it will keep other useful information that is required for running SHDOM scripts.
"""

# Basic imports
from __future__ import print_function
import os
from contextlib import contextmanager

# Get path to the SHDOM source code
PATH_REL_SHDOM = os.path.join('..', '..', 'shdom')     # relative path from examples
PATH_ABS_SHDOM = os.path.abspath(PATH_REL_SHDOM)

# Define a function to add the source directory to the pythonpath
def get_abs_path(relative_path):
    "Get the absolute path to shdom source" 
    return os.path.join(PATH_ABS_SHDOM, relative_path)
    
# Define a context for running a script
@contextmanager
def script_switch(name, runtrue):
    if runtrue:
        print("Running {0} ...".format(name))
        yield runtrue
        print("Done")
    if not runtrue:
        print("Skipping {0}.".format(name))
        exit
        yield runtrue
