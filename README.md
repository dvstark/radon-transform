# Radon Transform Codes. See Stark et al. 2018 for more detailed
  description of methodology. Send questions/comments to David
  V. Stark at david.stark@ipmu.jp

ds_radon: calculates variety of radon transforms for an input map

calc_radon_covar.pro: code to calculate the covariance matrix of the
radon transform given the covariance matrix of the input map.  This is
an IDL wrapper to python code. It requires IDL 8.5 or above and
Python. The Python bridge for IDL will need to be configured (see
http://www.harrisgeospatial.com/docs/Python.html)

example.pro: A simple tutorial showing how to run calc_radon_covar.pro