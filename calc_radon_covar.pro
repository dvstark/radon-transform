function calc_radon_covar,weights_arr,im_covar

  ;; NAME : CALC_RADON_COVAR
  ;;  
  ;; PURPOSE: Calculates tqhe covariance matrix of the Radon
  ;; transform. Meant to be used in conjunction with DS_RADON.pro
  ;;
  ;; CALLING SEQUENCE:
  ;; rt_covar = calc_radon_covar(weights_arr,im_covar)
  ;;
  ;; INPUTS
  ;;
  ;; WEIGHTS_ARR: A weights array of size MxN, where M is the number
  ;; of pixels in the input velocity map, and N is the number of
  ;; pixels in calculated radon transform map. The weights array is
  ;; 1/0, where values of 1 indicate a certain spaxel was used to
  ;; calculate the radon transform for a given index in the Radon
  ;; transform map.  This weights array is calculated in DS_RADON.pro.
  ;;
  ;; IM_COVAR: covariance matrix of velocity field used to calculate
  ;; radon transform
  ;;
  ;; OPTIONAL INPUTS: None.
  ;;
  ;; KEYWORD PARAMETERS: None.
  ;;
  ;; OUTPUTS: Returns the covariance matrix of the Radon
  ;; transform. Note: this can get BIG
  ;;  
  ;; PROCEDURE:
  ;;
  ;; Calculation described in Stark et al. 2018 (submitted), but
  ;; relatively straightforward matrix multiplication. Radon transform
  ;; covariance matrix is just:
  ;;
  ;; C_RT = W * C_V * W^T
  ;;
  ;; where W is the weights map and C_V is the velocity field
  ;; covariance matrix.
  ;;  
  ;; The actual calcuation is done in Python since scipy can handle
  ;; non-square sparse matrices, and sparse matrix speed things up
  ;; dramatically.  IDL seems to calculate radon transform faster,
  ;; while Python calculates covariance matrices faster.  I could
  ;; optimize python radon transform code, but for now, unholy
  ;; marriage of IDL/Python
  ;;
  ;; Requires IDL 8.5 or above and Python. The Python bridge will need
  ;; to be configured (see
  ;; http://www.harrisgeospatial.com/docs/Python.html)
  ;;
  ;; UPDATE HISTORY:
  ;;
  ;;     May 15 2018 (D. V. Stark): finally documented. Written a
  ;;     while back.

                            
  Python.w_arr = weights_arr
  Python.im_covar = im_covar
  >>>import numpy as np
  >>>import scipy.sparse as s
  >>>from scipy.sparse.linalg import inv
  >>>from scipy.sparse import find, triu
  >>>from astropy.io import fits
  >>>import sys
  >>>import time
  
  >>>sp_w_arr = s.csr_matrix(w_arr)
  >>>sp_im_covar = s.csr_matrix(im_covar)
  >>>sp_radon_covar=sp_w_arr.dot(sp_im_covar.dot(sp_w_arr.transpose()))
  >>>radon_covar = sp_radon_covar.toarray()

  radon_covar = Python.radon_covar
  return,radon_covar
  
end
