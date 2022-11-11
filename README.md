# magnetar.py

This contains a python package and model files to calculate the total polarized emission from the surface of a neutron star.

Depending on the version of numba installed, you may encounter errors.  You can replace 

  import numba
  
with

  import jnumba
  
to force the code to fall back regular python
