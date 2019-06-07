
function trace_radon_vm,im,rho,theta,mask=mask,smo=smo,error=error,ploton=ploton,inspect=inspect,silent=silent,covar=covar,mc_iter=mc_iter,invert=invert,center_only=center_only,tolerance=tolerance
  
  if keyword_set(inspect) then silent=0

;Program to trace the "ridge" in rho vs theta in an inverted absolute
;radon transform.  This is accomplished by making an initial guess at
;the position of the ridge (using a peak finder) and then fitting this
;feature with a von Mises function (i.e., a guassian for polar
;coordinates).
;
;INPUTS
;
;im -- the absolute radon transformed velocity field. The program
;      traces where this is minimized. Should be arranged so rho is
;      the y axis, theta is the x axis
;  
;rho -- radius array (units unimportant)
; 
;theta -- array of angular coordinates. Should be in radons and extend
;         from 0 to pi (I haven't yet built in routines to check for
;         this!
; 
;mask -- array of regions to ignore in im
;
;smo  -- (optional) smoothing factor given in normalized units e.g.,
;        0.1 means use a smoothing kernel equal to 0.1x the length of
;        the theta axis.  This is highly recommended to enable
;        reliable peak/trough detection, especially if there is noise
;        in the data (0.15 seems to work well)
;  
;error -- (optional)  error on im
;  
;ploton -- (optional) keyword to plot intermediate and final results
;          of fitting
;
;inspect -- (optional) keyword to halt the program at each step so
;           that the fitting can be inspected
;  
;silent -- (optional) keyword that turns off all text output in the
;          terminal
;  
;
;covar -- (optional) covariance matrix of input map
;
;mc_iter -- (optional) number of monte-carlo iterations when adding
;           random noise with covariance matrix in order to estimate
;           true uncertainty on fit parameters
;  
;invert -- keyword to invert the radon transform before calculation
;
;center_only -- keyword to only calculate centroid at rho=0
;
;tolerance -- set to ensure estimate centroid is not more than
;             "tolerance" degrees from previous measurement. Defaults
;              to 30 degrees
  
  
;OUTPUTS
;
;trace -- a structure containing the following information:
;
;         rho -- array of radii
;         theta -- best fitting there where im is minimized
;         etheta -- uncertainty on theta
;         status -- output from mpfitfun indicating fitting status
;                   (e.g., fit succeeded or failed). See mpfitfun
;                   documentation for details. Status=-999 means that
;                   the fit ended at one of the limiting values set in
;                   the code and is no good
;        slice_length -- fraction of data at given rho that was usable
;        bias_flag -- indicates that the estimated theta value may be
;                     compromised because of missing data.
;
;DEPENDENCIES
;
;        IDL astro library,
;        coyote library 
;        peakfinder
;        mpfit package
;        add_corr_noise
;        von_mises_line
;        unwrap_pa  
;  
;UPDATE HISTORY
;
;Oct 25, 2017 -- Code finally documented
;
;Jan 16, 2019 -- added linear term to fitting function to improve
;                performance. Fixed bug when taking median of best-fit
;                centers when bad fits exist.

  
  dim=size(im,/dim)
  if 1-keyword_set(theta) then theta=findgen(dim[0])
  if 1-keyword_set(rho) then rho=findgen(dim[1])-(dim[1]-1)/2.
  if 1-keyword_set(mask) then mask=1-finite(im)
  if 1-keyword_set(error) then error = im*0.+0.01
  if 1-keyword_set(silent) then silent = 0
  if 1-keyword_set(tolerance) then tolerance = 30.
  if 1-keyword_set(mc_iter) then mc_iter=1

  csz = size(covar,/dim)
  if keyword_set(covar) and n_elements(csz) eq 1 then begin
     print,'invalid covariance matrix, ignoring'
     covar = 0
  endif

  
  if keyword_set(ploton) then set_plot,'x'
  
                                ;define the 1d indices within the RT image
  sz=size(im,/dim)
  im_inds = make_array(sz[0],sz[1],/index)
  
  dtheta = abs(theta[1]-theta[0])

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;put flipped version of transform on either side, helps avoid;;;;;;;;
  ;;;;;wrapping issues. kind of brute force way of doing this, see if;;;;;;
  ;;;;;something;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;better;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  translate_radon,im,left,-1*dim[0]
  translate_radon,im,right,dim[0]
  im_ext = [left,im,right]
  theta_ext = [theta - (max(theta)+dtheta),theta,theta+max(theta)+dtheta]

  translate_radon,error,left,-1*dim[0]
  translate_radon,error,right,dim[0]
  error_ext = [left,error,right]
  
  translate_radon,mask,left,-1*dim[0]
  translate_radon,mask,right,dim[0]
  mask_ext = [left,mask,right]

  translate_radon,im_inds,left,-1*dim[0]
  translate_radon,im_inds,right,dim[0]
  im_inds_ext = [left,im_inds,right]
  
                                ;replace old arrays with these ones
  im_orig=im
  theta_orig=theta
  convert_radon,im_orig,im_inv,2,mask=mask
  
  im=im_ext
  theta=theta_ext
  error = error_ext
  mask=mask_ext
  im_inds = im_inds_ext

  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;initialize output;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;
  rho_arr = []
  theta_arr = []
  etheta_arr = []
  slice_length = []
  flag_arr = []
  status_arr = []                      
  lastparms = [-9999,-9999,-9999] ;store dummy parameters from previous row (these will eventually be used as guesses for the next row)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;iterate over each side, one row at a time start at rho=0, then;;;
  ;;;top rho>0 and rho<0;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  for side=0,1 do begin

     rho0_ind=where(rho eq 0)
     if rho0_ind eq -1 then rho0_ind = (dim[1]-1)/2.
     
     if side eq 0 then begin
        rho_start = rho0_ind[0]
        rho_end = dim[1]-1
        if keyword_set(center_only) then rho_end = rho0_ind[0]
        step=1
     endif else begin
        rho_start = rho0_ind[0]
        rho_end = 0
        if keyword_set(center_only) then rho_end = rho0_ind[0]
        step=-1
     endelse
     
     for j=rho_start,rho_end,step do begin

        if not silent then print,'rho=',rho[j]
        
        ;may want to bring this back
        ;catch,error_status
        ;; if error_status ne 0 then begin
        ;;    print,'error, skipping iteration'
        ;;    catch,/cancel
        ;;    goto,loop_end
        ;; endif
        

        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;;;pull out the slice, renormalize;;;
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        
        slice = im[*,j]
        if not keyword_set(invert) then slice=max(slice)-slice  
        slice_err = error[*,j]
        slice_mask = mask[*,j]
        slice_inds = im_inds[*,j] ;1d indices of the 2d radon transform
        sel=where(1-slice_mask) ;identify good regions
        norm = total(slice[where(1-slice_mask)]*dtheta) ;renormalizing
        slice=slice/norm
        slice_err = slice_err/norm
        slice_orig=slice        ;save for later

;        print,n_elements(slice),n_elements(slice_mask)
        
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;;;;;smooth the slice if necessary. Identify the "good" chunks;;
        ;;;;;and smooth them separately.  Also checking if;;;;;;;;;;;;;;
        ;;;;;there's enough data to proceed.  If not, skip to;;;;;;;;;;;
        ;;;;;next row;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        good_inds = get_chunks(slice,slice_mask,nchunks=nchunks) ;good regions
              
        ;if data doesn't cover 10 degrees, just ignore this row
        if total(1-slice_mask) lt 10*!dtor/dtheta then nchunks=0                      
        if nchunks eq 0 then continue
        
        if keyword_set(smo) then begin
                                ;smoothing each region separately

           smo_pix = dim[0]*smo
           for mm=0,nchunks-1 do begin
              if good_inds[1,mm]-good_inds[0,mm] ge smo_pix then $
                 slice[good_inds[0,mm]:good_inds[1,mm]] = $
                 smooth(slice[good_inds[0,mm]:good_inds[1,mm]],smo_pix,/edge_truncate)
           endfor
        endif


        
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;;;next step, find peaks and troughs in;;;;;
        ;;;data, use to define localized fitting;;;;
        ;;;region. This helps the algorithm from;;;;
        ;;;getting screwed up by other local;;;;;;;;
        ;;;minima/maxima in the radon transform;;;;;
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        
        if keyword_set(ploton) then begin
           set_plot,'x'
           window,0,xsize=1500,ysize=700
           !p.multi=[0,2,1]
           image_plot,im_inv,theta_orig,rho
           cgoplot,[-1000,1000],[rho[j],rho[j]],thick=3,color=cgcolor('yellow')
           p1=!p
           x1=!x
           y1=!y
           cgplot,theta,slice_orig,err_yhigh=slice_err,err_ylow=slice_err,xrange=[-!pi/2,1.5*!pi],/err_clip
           p2=!p
           x2=!x
           y2=!y
           
           cgoplot,theta,slice,thick=5,color=cgcolor('red')
        endif
           

        ;run tool to find peaks in the data
        p=peakfinder(slice[sel],theta[sel],/optimize,climits=climits,/widget_off,silent=silent)
        if p[-1] eq 0 then begin
           if not silent then print,'no peaks found at rho = '+strtrim(string(rho[j]),2)+'. Skipping'
           continue
        endif
        
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;;;filter out junk peaks. This is tricky. The program has;;;
        ;;;a built in routine to define a reliability limit based;;;
        ;;;on where the number of peaks stabilizes as a function;;;;
        ;;;of different reliability cuts.  However, I've found it;;;
        ;;;can be too restrictive, so I force it to not go above;;;;
        ;;;0.2.  0.2 is just found via trial and error;;;;;;;;;;;;;;
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              
        keep = where(p[-2,*] gt (climits[2] < 0.2))
                                ;if first iteration, keep focus on primary region between 0 and pi
        if j eq rho0_ind then keep=where(p[-2,*] gt (climits[2] < 0.2) $
                                         and p[1,*] ge -!pi/4 and p[1,*] le !pi+!pi/4)
        p=p[*,keep]
        
                                ;same process on the flipped data to
                                ;find troughs. If no troughs
                                ;found, it's handled below
        
        t = peakfinder(max(slice[sel])-slice[sel],theta[sel],/optimize,climits=climits,/widget_off,silent=silent)
        if t[0] ne 0 then begin 
           keep = where(t[-2,*] gt (climits[2] < 0.1))
           t=t[*,keep]
        endif           

        if keyword_set(ploton) then begin
           cgoplot,theta[sel[p[0,*]]],slice[sel[p[0,*]]],psym=4,color=cgcolor('red'),symsize=3,thick=3
           cgoplot,theta[sel[t[0,*]]],slice[sel[t[0,*]]],psym=4,color=cgcolor('blue'),symsize=3,thick=3
        endif
        

        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;;;If first iteration, find best peak (i.e., one with;;;;;;;;;;;;
        ;;;highest value). Do this also if there is no previous;;;;;;;;;;
        ;;;theta guess (but this shuld be extremely rare).;;;;;;;;;;;;;;;
        ;;;Otherwise, use output from previous iteration as guess;;;;;;;;
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;;;then we fit von mises function to the data for the data out;;;
        ;;;to the nearest troughs, or betwen the peak +/-!pi/4,;;;;;;;;;;
        ;;;whichever is a smaller range;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        if (j eq rho0_ind) or n_elements(theta_arr) eq 0 or (n_elements(lastparms) gt 0 and lastparms[0] eq -9999) then begin
           maxval = max(p[2,*],max_p_ind)
           theta_guess = p[1,max_p_ind]
           rt_guess=  p[2,max_p_ind]
           width=p[-1,max_p_ind]/2*dtheta/sqrt(2*alog(2))
           guess = [theta_guess,1./width^2,rt_guess,0]
        endif else begin
           guess = lastparms    ;values from previous successful iteration
        endelse
        guess = double(guess)
        
        parinfo={value:0.d,limited:[0,0],limits:[-9999.d,9999.d]} ;partinfo structure which is input to mpfit
        parinfo = replicate(parinfo,5)

        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;;;;setting up different initial guesses and ranges;;;;
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                                ;1) mean
        parinfo[0].limited=[0,0]
        parinfo[0].limits = [0d,2d*!pi]
                                ;if we're beyond the first iteration, ensure mean stays
                                ;within 'tolerance' degrees of previous value
        if j ne rho0_ind and lastparms[0] ne -9999 then begin
           parinfo[0].limited=[1,1]
           parinfo[0].limits=[guess[0]-tolerance*!dtor,guess[0]+tolerance*!dtor]
        endif            
        parinfo[0].value=guess[0]
        
                                ;2) width.  These limits are very non restrictive and just
                                ;ensure the fits don't become unrealistic
        parinfo[1].limited=[1,1]
        parinfo[1].limits=[(1./(90.*!dtor)^2),(1./(5.*!dtor)^2)] 
        parinfo[1].value=(guess[1] > 1.05*parinfo[1].limits[0]) < 0.95*parinfo[1].limits[1]
        
                                ;3) scale factor. Again, just ensuring they aren't
                                ;unrealistic but not restricting heavily
        parinfo[2].limited=[1,1]
        parinfo[2].limits=[0,10]
        parinfo[2].value=(guess[2] > 1.05*parinfo[2].limits[0]) < 0.95*parinfo[2].limits[1]
        parinfo[3].value=guess[3]
        
                                ;re-adjust guesses if needed using the
                                ;previously determined troughs
        guess = parinfo[*].value

                                ;default to +/- !pi/4
        min_fit_theta = guess[0]-!pi/4.
        max_fit_theta = guess[0]+!pi/4.

              ;find closest troughs
        if t[0] ne 0 then begin
                 
           theta_trough = t[1,*]
           
           diff = (guess[0] - theta_trough)
           diff = abs(diff*(diff gt 0)) + 9999*(diff lt 0)
           if total(diff lt 9999) gt 0 then begin
              mindiff = min(diff,ind)
              min_fit_theta = theta_trough[ind]
           endif
           
           diff = (guess[0] - theta_trough)
           diff = abs(diff*(diff lt 0)) + 9999*(diff gt 0)
           if total(diff lt 9999) gt 0 then begin
              mindiff = min(diff,ind)
              max_fit_theta = theta_trough[ind]
           endif 
        endif
        
              ;quick sanity check on ranges
        if abs(guess[0] - min_fit_theta) gt !pi/4. then min_fit_theta = guess[0] - !pi/4.
        if abs(guess[0] - max_fit_theta) gt !pi/4. then max_fit_theta = guess[0] + !pi/4.
        
        if keyword_set(ploton) then begin
           cgoplot,[min_fit_theta,min_fit_theta],[-100,100],color='blue',thick=5,linestyle=2
           cgoplot,[max_fit_theta,max_fit_theta],[-100,100],color='blue',thick=5,linestyle=2
        endif
        
              ;now select data for fitting and fit!
              
        selfit = where(theta ge min_fit_theta and $
                       theta le max_fit_theta and $
                       slice_err gt 0 and $
                       1-slice_mask,count)

        if count lt 5 then begin
           if not keyword_set(silent) then print,'too few points to do fitting. Skipping.'
           continue
        endif

                                ;going back to the unsmoothed data
                                ;with proper uncertainties (smoothing
                                ;helps with finding hte peaks)
        theta_fit = theta[selfit]
        slice_fit = slice_orig[selfit] 
        slice_err_fit = slice_err[selfit]
        slice_ind_fit = slice_inds[selfit]

        if keyword_set(ploton) then begin
           cgoplot,theta_fit,slice_fit,color=cgcolor('orange'),psym=4,thick=5
           cgoplot,theta[where(slice_mask)],slice[where(slice_mask)],psym=7,color='magenta',symsize=2,thick=4
        endif
        
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;;;;;set up the covariance matrix for this slice.  If the;;;;;;
        ;;;;;covar keyword is set, we extract the relevant indices;;;;;
        ;;;;;from the full 2d covariance matrix.  If not, we just;;;;;;
        ;;;;;create a diagonal covariance matrix with the assumed;;;;;;
        ;;;;;errors (if they are given);;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        slice_covar = fltarr(count,count)
        slice_cor = fltarr(count,count)
        if keyword_set(covar) then begin
        
           for xx = 0,count-1 do slice_covar[*,xx]=covar[slice_ind_fit,slice_ind_fit[xx]]


           for aa = 0,count-1 do begin
              for bb=0,count-1 do begin
                 slice_cor[aa,bb] = slice_covar[aa,bb]/sqrt(slice_covar[aa,aa]*slice_covar[bb,bb])
              endfor
           endfor

                                ;delete completely correlated pixels
                                ;out=delcorpix(slice_cor,xdata=theta_fit,ydata=slice_fit,covar=slice_covar,eydata=slice_err_fit)
           ;renormalize the same way as the slice
           slice_covar = slice_covar/norm^2

        endif else slice_covar[indgen(count),indgen(count)]=slice_err_fit^2
 
        
                                ;correct covariance matrix if needed
        blah = add_corr_noise(slice_fit,slice_covar,return_covar=return_covar)
        slice_covar_orig = slice_covar
        slice_covar = return_covar

        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;;;do monte-carlo method to add;;;;;;;;;;;
        ;;;correlated noise and estimate;;;;;;;;;;
        ;;;uncertainties.  Note to self;;;;;;;;;;;
        ;;;consider multithreading this process;;;
        ;;;to speed up;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        func_name='von_mises_line'

                                ;will keep track of fit parameters,
                                ;errors, and fit status for mc mode
        parms_stack = []
        eparms_stack = []
        status_stack = []

        
        for ll = 0,mc_iter-1 do begin

           if mc_iter gt 1 then newslice_fit = add_corr_noise(slice_fit,slice_covar) $
           else newslice_fit = slice_fit

           parms = mpfitfun(func_name,theta_fit,newslice_fit,slice_err_fit,yfit=out,$
                            bestnorm=chisq,perror=eparms,status=status,parinfo=parinfo,$
                            quiet=1)

                                ;the true maximum is not always at
                                ;parms[0]. Use derivative to get
                                ;better estimate
           
           df = deriv(theta_fit,out)
 ;          if rho[j] gt 11 then print,rho[j]
           z = find_zeros(theta_fit,df)
;           if rho[j] gt 11 then print,z
           if z ne !NULL then begin
              sep = abs(parms[0] - z)
              minsep = min(sep,index)
              parms[0] = z[index]
           endif else status = 0

           
           ;set all parameters to -1 if there was a failure
           if status eq 0 then begin
              eparms = [-1,-1,-1,-1,-1]
              parms = [-1,-1,-1,-1,-1]
           endif
                    
           ;keep track of parameters, error, and status of fit
           parms_stack = [[parms_stack],[parms]]
           eparms_stack = [[eparms_stack],[eparms]]
           status_stack = [status_stack,status]
        endfor

        ;estimate fit parameters and uncertainties if doing mc
        if mc_iter gt 1 then begin
           newparms = fltarr(n_elements(parms))
           neweparms = fltarr(n_elementS(parms))
           for ll = 0,n_elements(parms)-1 do begin
              sel=where(status_stack gt 0,count)
              newstatus = 1
              if count gt 0 then begin
                 newparms[ll] = median(parms_stack[ll,sel])
                 neweparms[ll] = robust_sigma(parms_stack[ll,sel])
              endif else begin
                 newparms[ll] = -999
                 neweparms[ll]=-999
                 newstatus = -999
              endelse
           endfor
           if 1-keyword_set(silent) then print,'rho, parameters, errors,fit status: ',rho[j],newparms[0],neweparms[0],status
           parms = newparms
           eparms = neweparms
           status = newstatus
        endif

        if (n_elements(eparms) eq 0) or (total(1-finite(parms)) gt 0) then status=-999

        if keyword_set(ploton) then begin
           cgoplot,[parms[0],parms[0]],[-100,100],color='turquoise',thick=1
           cgoplot,[parms[0],parms[0]]-eparms[0],[-100,100],color='turquoise',thick=1,linestyle=2
           cgoplot,[parms[0],parms[0]]+eparms[0],[-100,100],color='turquoise',thick=1,linestyle=2
        endif
        



;        if keyword_set(ploton) and status ne -999 then cgoplot,theta_fit,von_mises_line(theta_fit,parms),color=cgcolor('green'),linestyle=2,thick=5
        
        if total(1-finite(parms)) eq 0 then begin
                                ;save the fit output
           rho_arr = [rho_arr,rho[j]]
           theta_arr = [theta_arr,parms[0]]
           etheta_arr = [etheta_arr,eparms[0]]
           status_arr = [status_arr,status]
           slice_length = [slice_length,total(1-slice_mask)/n_elements(slice_mask)]
                 
                 
                                ;see if earlier theta value +/- 95% confidence interval
                                ;was in masked region this time. This indicates our new
                                ;estimate may be biased by missing data.
                 
           flag=0
           if j ne rho0_ind and n_elements(theta_arr) gt 1 then begin
              lastgood = where(flag_arr ne 1 and status_arr gt 0)
              lastgood = lastgood[-1]
              
              nmask = total(slice_mask[where(theta gt (theta_arr[lastgood]-3*etheta_arr[lastgood]) and theta lt (theta_arr[lastgood]+3*etheta_arr[lastgood]))])
              if nmask gt 0 then flag=1
              if flag and not keyword_set(silent) then begin
                 print,''
                 print,'***FLAG***'
                 print,''
              endif

              if status le 0 and not keyword_set(silent) then print,'***BAD FIT***'
           endif
                 
           flag_arr = [flag_arr,flag]


           if 1-flag and status gt 0 and eparms[0] ne 0 then lastparms=parms $
           else if (flag or status le 0 or eparms[0] eq 0) and lastparms[0] ne -9999 then  lastparms = lastparms $
           else lastparms=parms*0-9999
                 
        endif

                 
        if keyword_set(inspect) then begin
           print,'guess: ',guess
           print,parms[0]
           print,eparms[0]
           print,status
           print,flag
           input=''
           read,input
        endif
        loop_end:
     endfor
  endfor
  
    
  trace = {rho:-1,theta:-1,etheta:-1,status:-1,slice_length:-1,bias_flag:-1}
  if n_elements(rho_arr) gt 0 then begin

                                ;fix jumps around 180/360 if present
     unwrap_pa,rho_arr,theta_arr,srtarr=srt
     status_arr = status_arr[srt]
     flag_arr = flag_arr[srt]
     etheta_arr = etheta_arr[srt]
     slice_length = slice_length[srt]
     
     if keyword_set(ploton) then begin
     
        cgplot,rho_arr,theta_arr/!dtor,err_yhigh=etheta_arr/!dtor,err_ylow=etheta_arr/!dtor,/err_clip,psym=1,/ynozero
        sel=where(flag_arr,count)
        if count gt 0 then oplot,rho_arr[sel],theta_arr[sel]/!dtor,psym=4,color=cgcolor('red'),thick=5
        sel=where(status_arr lt 1,count)
        if count gt 0 then cgoplot,rho_arr[sel],theta_arr[sel]/!dtor,psym=1,color=cgcolor('green'),thick=5
     endif
   
     srt = sort(rho_arr)
     rho_arr = rho_arr[srt]
     theta_arr = theta_arr[srt]
     etheta_arr = etheta_arr[srt]
     status_arr = status_arr[srt]
     slice_length = slice_length[srt]
     flag_arr = flag_arr[srt]

  ;we'll have two rho = 0 entries, remove 1
     u = uniq(rho_arr)
     
     trace = {rho:rho_arr[u],theta:theta_arr[u],etheta:etheta_arr[u],status:status_arr[u],slice_length:slice_length[u],bias_flag:flag_arr[u]}


     
  endif
  
  return,trace
;find rho=0
end
