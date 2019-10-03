function ds_radon, im, normal = normal, theta = theta, rho = rho, $
                   mask = mask, weight = weight, error = error, $
                   xmin = xmin, ymin = ymin, dx = dx, dy = dy, $
                   median = median, robust=robust, $
                   aperture=aperture, do_covar=do_covar, $
                   im_covar=im_covar,radon_covar = radon_covar
;
; NAME: DS_RADON
;
; PURPOSE: Carry out a Radon-inspired transform on an image  
;
; CATEGORY: Image analysis routines.
;
; CALLING SEQUENCE: 
;       res = DS_RADON(IM, [KEYWORDS=KEYWORDS])
;
; INPUTS: 
;
;    IM: An MxN image - the routine will work in the precision of this
;        image. The routine will probably work if given long or int
;        input but float or double is recommended for proper
;        operation. 
; 
;
; OPTIONAL INPUTS: None.
;
; KEYWORD PARAMETERS: 
;
;              DX: The resolution to apply in the x-direction. The
;                  default is 1.
;              DY: Ditto for the y-direction
;           ERROR: An error image corresponding the IM image input
;                  with estimates of the error on each pixel. Gaussian
;                  statistics is assumed throughout and the error is
;                  expected to be given as the standard deviation.
;            MASK: A mask image with the same dimension as IM. If
;                  given the pixels that are set to zero on this image
;                  are excluded from the calculation. 
;          NORMAL: If this keyword is set, the operation of the
;                  routine is like the built-in RADON transform. This
;                  is mostly intended for debugging purposes and
;                  should not be used for production code - use RADON
;                  instead. 
;             RHO: [Output] - the values for the rho axis of the
;                  output image (this is the y-axis). This is also
;                  returned in the result structure.
;           THETA: [Output] - the theta axis (x-axis) of the output
;                  image. This is also returned in the result
;                  structure. 
;          WEIGHT: A weight image with the same dimensions as IM. The
;                  integral along each ray is weighted by this image.
;            XMIN: The value of the minimum x value. The default is
;                  the same as for RADON: -(m-1)/2 if IM is MxN
;            YMIN: The value of the minimum y value. The default is
;                  the same as for RADON: -(n-1)/2 if IM is MxN
;          MEDIAN: If keyword set, use the median instead of the
;                  average to normalize the absolute radon transform
;         APERTURE: Use only pixels within radius r of given rho/theta    
;         DO_COVAR: Set to calculate the covariance matrix of the
;                   Radon transform
;         IM_COVAR: Covariance matrix of input velocity field
;  
;           ROBUST: Keyword to remove outliers when calcalating RT at
;                   given location. This is experimental.
;  
;         IM_COVAR: Covariance matrix of input velocity field   
  
; OUTPUTS: The routine returns a structure with various transforms
;          included. Each transform is discussed in detail under
;          PROCEDURE below.  The keys in the result structure are as
;          follows (when the NORMAL keyword is set only MAP, RHO and
;          THETA contains anything useful):
;
;      MAP: The standard transformed image. If the NORMAL keyword is
;           set this is the standard RADON transform.
;    THETA: The theta axis.
;      RHO: The rho axis.
;     DIFF: The velocity diffence transform.
;     VMAX: The maximum velocity transform.
;     VMIN: The minimum velocity transform.
;  N_ZEROC: The zero-crossings transform.
;   LENGTH: The length of integration along each line.
;    DDIFF: If present, this contains the uncertainty estimate for
;           DIFF. 
;    ERROR: If present, this contains the uncertainty estimate for
;           MAP.
; MASKFRAC: The fraction of pixels in a line integral that were
;           removed due to masking; can be used to identify regions
;           where measurements may be compromised by missing/bad data.
;     EDGE: Indicates how many pixels included in an aperture extend
;           beyond the edge of the map.  Values with 999 indicate no data.
;           This improves upon the maskfrac flag which does not account for
;           This improves upon hte maskfrac flag which does not account for
;           edges explicitly
;RADON_COVAR: calculated radon transform covariance matrix
;  
; SIDE EFFECTS: For large images a considerable amount of memory might
;               be consumed as at least six transforms are calculated
;               (eight if an error image is given). 
;
; RESTRICTIONS: Very few error checks are made. 
;
; PROCEDURE: 
;   The standard Radon transform of a 2D function f(x,y) involves
;   carrying out line integrals along each possible line crossing the
;   array, parametrised by the normal distance to the origin, rho, and
;   the angle of this normal with the x-axis, theta. Here we provide
;   some alternatives to this. We still use the same parametrisation
;   and algorithm, but we modify the operator on f. 
;
;   The basic approach, of chosing the nearest neighbour etc. are
;   exactly as discussed for RADON so the reader is directed there for
;   more information on those aspect. This is also the case if you
;   want to use the NORMAL keyword (but you shouldn't - use RADON
;   instead!) 
;
;   When a given line is chosen, we then figure out whether there are
;   any points of this within the image - if not a zero is placed in
;   the output MAP. We then extract the function along this line. If a
;   mask is given we use this at this point and mask out any points
;   that have MASK = 0. We then calculate the following transforms on
;   the masked g(s). (here labelled by the name of the key in the
;   output structure) 
;  
;      MAP: We subtract the average of the function and then integrate
;           the absolute value of the resulting function. In formuli
;           this is: 
;   
;                MAP(t,r) = int |g(s)-<g(s)>| ds
;           
;           The intention of this is to highlight directions that have
;           large changes in the function value without having to
;           calculate derivatives.
;
;     VMAX: The value of the transform is the maximum function values
;           along the line. The name comes from the fact that the
;           function is expected to be a velocity map.
;     VMIN: The minimum function value (velocity) along the line.
;     DIFF: The difference between VMAX and VMIN. The intention of
;           this transform is to provide an alternative measure of the
;           change of the function along a given line. This is likely
;           to be more sensitive to outliers and noise than MAP.
;
; The effect of the weight and error arrays makes relatively small
; changes to the above. The VMAX, VMIN and DIFF arrays are not
; affected, but for the main transform we change the integration to be
;
;     int g(s) w(s) ds/int w(s) ds
;
; in the case of weights.  For the errors only MAP and DIFF are given
; as ERROR and DDIFF respectively.
;
; We then calculate the number of points we integrate over and
; store this as LENGTH - this is the effective length and takes into
; account masking, but not points with zero weight. 
;
; Finally we estimate a monotonicity index. This is given by the
; number of zero-crossings of the derivative of the function along
; each line. This is constructed as follows:
;
;    - Calculate Delta = f(i+1)-f(i) 
;    - Delta -> Delta/abs(Delta) [now only -1 or 1 or NaN]
;    - Where Delta is undefined set it to zero.
;    - Use UNIQ to remove all identical consecutive numbers and find
;    the final number of objects - this gives the number of
;    zero-crossing.
;
; The number of zero crossings of the derivative is stored as
; N_ZEROC. 
; 
; EXAMPLE:
;
; MODIFICATION HISTORY:
;
;       Oct 4, 2005, J. Brinchmann (jarle@astro.up.pt)
;	   Documented routine.  
;
;       Sep 28, 2015, D. Stark (david.stark@ipmu.jp)
;         Added option to use median to normalize absolute radon
;         transform rather than avg which could be strongly affected by
;         outliers
;         Added outlier identification (probably needs tweaking)
;
;       June 7, 2017, D. Stark
;         Added LMASK to output
;
;       Aug 1, 2017 D. Stark
;         Added EDGE flag to output
;         Added option to calculate covariance matrix using python  
;  
;-
  

   ;;-----------------------------
   ;; Get dimensions of image.
   ;;-----------------------------

  do_covar=keyword_set(do_covar)
  radon_covar = -1
  
   dims = size(im, /dimen)
   m = dims[0]
   n = dims[1]

   ;create array of x and y indices (only needed if including covariance)
   xind = intarr(m,n)
   yind = intarr(m,n)
   for ll = 0,n-1 do xind[*,ll]=indgen(m)
   for ll = 0,m-1 do yind[ll,*]=indgen(n)
   
   ;; Notice that we use now the same default as RADON because this
   ;; makes it easier to deal with rotation.
   if (n_elements(xmin) eq 0) then xmin = -(m-1)/2.
   if (n_elements(ymin) eq 0) then ymin = -(n-1)/2.
   if (n_elements(dx) eq 0) then dx = 1.0
   if (n_elements(dy) eq 0) then dy = 1.0

   do_mask = (n_elements(mask) gt 0)
   do_weight = (n_elements(weight) gt 0)
   do_error = (n_elements(error) gt 0)

   xmax = (xmin + m)
   ymax = (ymin + n)

   ;;---------------------------------------------------------------
   ;; Set the number of rho, theta etc variables using the default
   ;; settings of the IDL routine.
   ;;---------------------------------------------------------------
   theta_min = 0.0
   theta_max = !pi
   drho = 0.5*sqrt(dx^2 + dy^2)

   n_theta = ceil(!pi*sqrt(0.5*(m^2 + n^2)))
   ;; This gives 
   n_rho = 2L*ceil(sqrt(max([xmax^2 + ymax^2, xmin^2+ymin^2]))/drho)+1L

   dtheta = (theta_max-theta_min)/n_theta

   if (n_elements(rmin) eq 0) then rmin = -0.5*(n_rho-1)*drho

   
   ;; 0-> !pi-dtheta
   theta = findgen(n_theta)*dtheta ;+ 0.01
   cost = cos(theta)
   sint = sin(theta)
   sint[0] = sint[1]

   rho = findgen(n_rho)*drho + rmin
   
   ;;-------------------------
   ;; Create transform image
   ;;-------------------------
   tim = make_array(n_theta, n_rho, type=size(im, /type))
   diff = tim
   vmax = tim
   vmin = tim
   crosszero = tim
   n_zeroc = long(tim)-1
   length = tim
   maskfrac = length*0.+1
   edge = maskfrac*0.+999
   if (do_error) then begin
      dtim = tim
      ddiff = tim
   endif

   inds_radon = make_array(n_theta, n_rho,/index) ;1d indices of radon transform
   inds_map = make_array(m,n,/index)  ;1d indices of input map
   weights_arr = intarr(n_elements(inds_map),n_elements(inds_radon)) ;1 if it's included in calculation at inds_map[i]
   
   ;; Calculate some variables for the paths
   a = -(dx/dy)*cost/sint
   ;; When calcuating (rho - b_part[i])/norm[i] this gives b... 
   b_part = xmin*cost + ymin*sint
   norm = dy*sint

   ii = where(abs(sint) le sqrt(2)/2., n_low)
   if (n_low gt 0) then begin
      a[ii] = 1.0/a[ii]
      norm[ii] = dx*cost[ii]
   endif
   
   ;; Carry out calculation
   xinds = findgen(m)
   yinds = findgen(n)

                                ;define center index
   xind0=max(xinds)/2.
   yind0=max(yinds)/2.

   for i = 0L, n_theta-1 do begin
      if (abs(sint[i]) gt sqrt(2.)/2.) then begin
         path = 1 
         pre_factor = dx/abs(sint[i])
      endif else begin 
         path = 0
         pre_factor = dy/abs(cost[i])
      endelse
      b = (rho-b_part[i])/norm[i]
      
      for j = 0L, n_rho-1 do begin
         ;; Construct the ray for this rho, theta - depending on the
         ;; path 
         if (path eq 1) then begin 
            iy = round(a[i]*xinds + b[j])
            ix = xinds
         endif else begin
            ix = round(a[i]*yinds + b[j])
            iy = yinds
         endelse

         
         
        

         xcent = rho[j]*cos(theta[i])
         ycent = rho[j]*sin(theta[i])
         rad = sqrt(((ix-xind0)-xcent)^2 + ((iy-yind0)-ycent)^2)
         
                   ;  if abs(rho[j]) lt 10 and rho[j] gt 0 then stop

         ;; Ensure that the indices are within the image and do not
         ;; wrap. We want at least 3 pixels
         use = where(ix ge 0 and ix lt m and iy ge 0 and iy lt n, $
                     n_use)
         edgeflag = 0
         
         if keyword_set(aperture) then begin
            use = where(ix ge 0 and ix lt m and iy ge 0 and iy lt n and rad lt aperture, n_use)
            if (n_use le 2) then continue ;if not enough data, just skip
            
            edgeflag = total((ix[use]+aperture gt m) or (ix[use]-aperture lt 0) $
                             or (iy[use]+aperture gt n) or (iy[use]-aperture lt 0))
         endif 
         
                                ;if filterr is set, do the remove any
                                ;that have radii </> the current rho         
         
         ix = ix[use]
         iy = iy[use]

         ;; Carry out the evaluation of the integrand (to do: make
         ;; "normal" version accept masks/set maskfrac variable properly
         if (keyword_set(normal)) then begin
            func = im[ix, iy]
            if (do_mask) then begin
               t_mask = mask[ix, iy]
               use2 = where(t_mask eq 1, n_use2)
               frac_mask = float(n_elements(t_mask) - n_use2)/n_elements(t_mask)
               if (n_use2 le 1) then continue
               func = func[use2]
            endif
            func_norm = func
         endif else begin
            ;; Here we have make our required adjustments.
            func = im[ix, iy]
            if (do_mask) then begin
               t_mask = mask[ix, iy]
               use2 = where(t_mask eq 1, n_use2)
               frac_mask = float(n_elements(t_mask) - n_use2)/n_elements(t_mask)
               if (n_use2 le 1) then continue
               func = func[use2]
            endif
            if 1-keyword_set(median) then func_norm = (func-avg(func)) else $
               func_norm = func - median(func)
            func_raw = func
            func = abs(func_norm)
         endelse
         func_raw = func
         
                                ;find outliers
         outlier=func_norm*0.
         if keyword_set(robust) then begin
            resistant_mean,func,4,mean,sigma,num_rej,goodvec=goodvec
            outlier[*]=1
            outlier[goodvec]=0
         endif
         
         
         
         
;         print, 'THE VMAX & VMIN AND DIFF operators should operate on  the unadulterated function!'
;         print, 'Should also check carefully how the WEIGHT array is  used for VMAX, VMIN etc.'
;         


         mx = max(func_norm, i_mx)
         mn = min(func_norm, i_mn)

         
         if (do_weight) then begin
            w = weight[ix, iy]
            wnorm = total(w)
            if (wnorm eq 0.0) then continue
            if (do_mask) then w = w[use2]
            if keyword_set(robust) then tim[i, j] = pre_factor*total(func*w*(1-outlier))/wnorm else $
               tim[i, j] = pre_factor*total(func*w)/wnorm                       
            if (do_error) then begin
               err = error[ix, iy]
               if (do_mask) then err = err[use2]
               dtim[i, j] = pre_factor*sqrt(total((err*w)^2)/total(w^2))
            endif
         endif else begin
            ;no weights
            if keyword_set(robust) then tim[i, j] = pre_factor*total(func*(1-outlier)) else $
               tim[i, j] = pre_factor*total(func) 
            if (do_error) then begin
               err = error[ix, iy]
               if (do_mask) then err = err[use2]

               if 1-do_covar then begin
                  dtim[i, j] = pre_factor*sqrt(total(err^2))
                  ddiff[i, j] = sqrt(err[i_mx]^2+err[i_mn]^2)
               endif else begin
                  ;need to properly account for covariance!
                                ;for each pixel in i_x,i_y, need to
                                ;find all covariance pixels
                  iind_1d = inds_map[ix,iy]
                  npix = n_elements(iind_1d)
                  ;create covariacne matrix for just this subset
                  covar = fltarr(npix,npix)
                  for xx = 0,npix-1 do begin
                     covar[*,xx]=im_covar[iind_1d,iind_1d[xx]]
                  endfor
                  weights=intarr(npix)+1 ;generalize this for if the "weights" keyword is actually used
                  err_sq = weights##covar
                  err_sq = err_sq##transpose(weights)
                  dtim[i,j] = pre_factor*sqrt(err_sq)
                                 
               endelse
               
               
                  
            endif
               
         endelse
         diff[i, j] = mx-mn
         vmax[i,j ] = mx
         vmin[i, j] = mn
         length[i, j] = n_elements(func)
         edge[i,j] = edgeflag
         radon_pix = inds_radon[i,j]
         map_spaxes = inds_map[ix,iy]
         weights_arr[map_spaxes,radon_pix] = 1;./n_elements(map_spaxes)
         if do_mask then maskfrac[i,j] = frac_mask

                                ;find whether it crosses the zero velocity contour
         lowv=total(func_raw lt 0)
         highv=total(func_raw gt 0)
         crosszero[i,j] = lowv ne 0 and highv ne 0

         ;; Now calculate a monotonicity index. This is simply the
         ;; number of zero-crossings of the derivative.
;;==         delta_func = func_norm-shift(func_norm,1)
         ;; Discard the first since it comes from wrapping
;;==         delta_func = delta_func[1:*]
;;==         delta_func = delta_func/abs(delta_func)
;;==         izero = where(finite(delta_func) eq 0, n_zero)
;;==         if (n_zero gt 0) then delta_func[izero] = 0.0
         
;;==         ui = uniq(delta_func)
;;==         n_zeroc[i, j] = n_elements(ui)
;;==         if (tim[i, j] ne 0.0 ) then stop

      endfor

   endfor
 
   if do_covar then begin
   
      weights_file = 'weights_arr.fits'
      im_covar_file = 'im_covar.fits'
      radon_covar_file = 'radon_covar.fits'
      radon_covar_diag_file = 'radon_covar_diag.fits'
      

      if keyword_set(covar_file_root) then begin
         weights_file = covar_file_root + '_' + weights_file
         im_covar_file = covar_file_root + '_' + im_covar_file
         radon_covar_file = covar_file_root + '_' + radon_covar_file
         radon_covar_diag_file = covar_file_root + '_' + radon_covar_diag_file
      endif
      
      
      print,'creating covariance matrix for radon transform'
      radon_covar = calc_radon_covar(weights_arr,im_covar)
      covar_dim = size(radon_covar,/dim)
      err_1d = radon_covar[indgen(covar_dim[0]),indgen(covar_dim[0])]

           
      dtim_alt = reform(sqrt(err_1d),n_theta,n_rho)
      dtim = dtim_alt
   endif
   

      
   res = {map: tim, theta: theta, rho: rho, diff: diff, vmax: vmax, $ 
          vmin: vmin, n_zeroc: float(n_zeroc), length: length, xmin: xmin, $ 
          ymin: ymin, crosszero:crosszero, maskfrac: maskfrac, edge:edge}
   
   if (do_error) then res = create_struct(res, $ 
                                          'error', dtim, $ 
                                          'ddiff', diff)
 
  return, res

end
