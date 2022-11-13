# magnetar.py

This contains a python package and model files to calculate the total polarized emission from the surface of a neutron star.

Depending on the version of numba installed, you may encounter errors.  You can replace 

  import numba
  
with

  import jnumba
  
to force the code to fall back regular python

Some key routines are written in C, so go to the Magnetar directory and type make to compile these.  This should work fine on
Mac and Linux systems, but Windows users will have to make some alterations to the python loading routines and the Makefile.
