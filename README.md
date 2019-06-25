# Radon Transform Codes. See Stark et al. 2018 for more detailed
  description of methodology. Send questions/comments to David
  V. Stark at david.stark@ipmu.jp

Many of the codes rely on outside packages.  These include:
1) The IDL astronomy library
2) The Coyote Library
3) The MPFIT library

IDL 8 or greater is needed for most codes (or whichever version first
introduced indexing using -1 to select the last entry in an array).
IDL 8.5 needed for the internal calls to python (just
calc_radon_covar.pro when calculating the covariance matrix of the
radon transform).

###Quick descriptions of codes###

add_corr_noise.pro -- adds random correlated noise to a data set using
its covariance matrix

calc_radon_covar.pro: code to calculate the covariance matrix of the
radon transform given the covariance matrix of the input map.  This is
an IDL wrapper to python code. It requires IDL 8.5 or above and
Python. The Python bridge for IDL will need to be configured (see
http://www.harrisgeospatial.com/docs/Python.html)

convert_radon.pro -- a simple rescaling routine

ds_radon: calculates variety of radon transforms for an input map

example.pro A simple tutorial showing how to run ds_radon

find_zeros.pro -- used internally to find zeros in a smooth function

logger.pro -- just a little helper routine to log failures

peakfinder.pro -- identifies peaks in a 1D array

radon_wrapper_working.pro -- a wrapper specific for MaNGA maps
files. Loads and processeds the maps files, calculates the RT, then
traces it. Outputs information into a fits file, and plots results

trace_radon.pro -- traces the ridge in the absolute radon transform
defined by where the transform is minimized as a function of radius

translate_radon.pro -- translates a radon transform map along the
theta axis, flipping as necessary at 0 and 180.

upwrap_pa.pro -- unwraps an array of postion angles so there's no
dramatic jump at +/- 180

von_mises_line.pro -- a von-mises + linear function

wrapper_example.pro -- example showing how to run radon_wrapper.pro