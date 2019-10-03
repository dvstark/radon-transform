pro radon_wrapper_minor,name,output,outfile,xshift=xshift,yshift=yshift,plotfile=plotfile, radon_ap=radon_ap,err_log=err_log,do_covar=do_covar,mc_iter=mc_iter,quiet=quiet,eh_on=eh_on,mapstype=mapstype,drpallpath=drpallpath,mapspath=mapspath,stars=stars,sncut=sncut, maps_suffix = maps_suffix

  print,'analyzing '+name
  
                                ;This is a wrapper
                                ;to run the radon transform and
                                ;tracing algorithm on MaNGA maps
                                ;files. Specialized for my purposes so
                                ;may need significant editing.
  
                                ;to do: keyword inheritance

;  INPUTS:
;
; name -- galaxy plate-ifu
;  
; OUTPUT
;  
; output -- just returns 1 if worked or -1 if failed  
;
; outfile  -- output fits file which holds various data products.  It
;             contains 13 extensions. See radon_wrapper.readme for
;             more info
;
;
;
;
; KEYWORDS:
;
; xshift/yshift -- optional pixel shifts to apply in x/y direction
;                 (otherwise center assumed to be at center of array)
; 
; plotfile -- name of plot file (will be a png, but don't put
;             .png at the end)
;
; radon_ap -- the 'radon aperture' in pixels to use. This is
;             technically half of the full length of the line segment
;             we integrate over during the transform
;
;err_log -- log file where we'll store error messages
;
;do_covar -- keyword to turn on calculation of the full covariance
;            matrix of the radon transform (warning, this step can
;            increase runtime significantly)
;
;mc_iter -- number of monte-carlo iterations to run when estimating
;           errors during the trace step (warning: also increases run
;           time a lot)
;
;quiet -- turn off output (may not fully work)
;
;eh_on -- turn on error handling. Good if running over many galaxies
;         because code won't crash if there's an error
;
;maps_type -- type of maps file (SPX, HYB10, etc)
;
;drpalpath -- where to find the drpall file
;
;mapspath -- hwere to find the maps files
;
;stars -- set to 1 to run on stellar velocity fields rather than
;         ionized gas
;
; sncut -- required signal-to-noise per spaxel (default = 3)  
;
;update log
;
;July 5, 2019 -- added SN keyword and ability to mask low S/N regions
;                for stellar velocity fields

  
;;;;;Some initial setup;;;;
  
  if 1-keyword_set(eh_on) then eh_on = 1
  if 1-keyword_set(xshift) then xshift = 0
  if 1-keyword_set(yshift) then yshift = 0
  if 1-keyword_set(radon_ap) then radon_ap = 9e9
  if 1-keyword_set(do_covar) then do_covar = 0
  if 1-keyword_set(sncut) then sncut=3
  if not keyword_set(mapstype) then mapstype='HYB10'
  if not keyword_set(mapspath) then mapspath='~/manga/fits/mpl7/dap/'+mapstype+'/';HYB10/'
  if not keyword_set(drpallpath) then drpallpath = '/home/dstark/manga/fits/drpall-v2_5_3.fits'
  if 1-keyword_set(maps_suffix) then maps_suffix = 'MILESHC-MILESHC'


                                ;if there's an error, print it
                                ;to the log, but then just
                                ;exist script rather than crashing
  catch,error_status
  if (error_status ne 0) and eh_on then begin
                                ;write down error
     if 1-keyword_set(err_log) then err_log = './logs/radon_wrapper_log_'+timestamp()
     if 1-file_test(err_log) then spawn,'touch '+err_log ;just create file
     openu,1,err_log
     print,error_status,!error_state.msg
     printf,1,''
     printf,1,name+' failed'
     printf,1,'Error index: ',error_status
     printf,1,'Error message: ',!error_state.msg
     help,/last_message,output=message
     printf,1,message
     close,1
     output=-1
     goto,the_end
  endif

                                ;default output in case of failure
  output = -1

  ;read drpall file
  db=mrdfits(drpallpath,1,/silent)
  ii = where(db.plateifu eq name,count) ;index of this galaxy in the database
  dbind = ii
  if count eq 0 then begin
     logger,err_log,name+' not found in drpall file'
;     openw,1,err_log,/append
;     printf,1,''
;     printf,1,name+' not found in drpall file'
;     close,1
     goto,the_end
  endif
  

                                ;read in velocity field
  plate = db[dbind].plate     
  ifudsgn = db[dbind].ifudsgn  
 
  mapsfile = mapspath+'manga-'+strtrim(string(plate),2)+'-'$
             +strtrim(ifudsgn,2)+'-MAPS-'+mapstype+'-'+maps_suffix+'.fits.gz'
  
;             +strtrim(ifudsgn,2)+'-MAPS-'+mapstype+'-GAU-MILESHC.fits.gz'


  if 1-file_test(mapsfile) then begin
     logger,err_log,'maps file not found for '+name
     goto,the_end
  endif 
  
     ;;;;first, read in the coordinates. We will initially remap the
     ;;;;velocity field so that the optical center is exactly at the
     ;;;;center of the 2d map. Calculating new grid below
     
  dummy = mrdfits(mapsfile,0,shdr,/silent)
  coords = mrdfits(mapsfile,'SPX_SKYCOO',/silent)
  xgrid = coords[*,*,0]
  ygrid = coords[*,*,1]
  
  dx = 0.5
  maxx = round(abs(min(xgrid)) < abs(max(xgrid)))
  maxy = round(abs(min(ygrid)) < abs(max(ygrid)))
  newxarr = findgen((maxx < maxy)*2/dx+1)*0.5 - (maxx < maxy)
  newyarr = newxarr
  newxarr = -1*newxarr       ;just to keep orientation the same. Note that this is in arcseconds, not pixels!
  newxgrid = rebin(newxarr,n_elements(newxarr),n_elements(newyarr))
  newygrid = transpose(rebin(newyarr,n_elements(newyarr),n_elements(newxarr)))
  
     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;;;;;read in velocies now;;;;;;;
     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if 1-keyword_set(stars) then begin
  
     vels = mrdfits(mapsfile,'EMLINE_GVEL',hdr,/silent)
     vels_ivar = mrdfits(mapsfile,'EMLINE_GVEL_IVAR',/silent)
     vels_masks=mrdfits(mapsfile,'EMLINE_GVEL_MASK',/silent) ;bad pixel mask
     gflux = mrdfits(mapsfile,'EMLINE_GFLUX',/silent)
     gflux_ivar=mrdfits(mapsfile,'EMLINE_GFLUX_IVAR',/silent)
     
     pxscl=abs(sxpar(hdr,'PC2_2')*3600)
     
     halpha_ind = 18
     vel = vels[*,*,halpha_ind]
     vel_ivar = vels_ivar[*,*,halpha_ind]
     evel = 1/sqrt(vel_ivar)
     evel[where(1-finite(evel))]=0
     gf = gflux[*,*,halpha_ind]
     gf_ivar = gflux_ivar[*,*,halpha_ind]
     vel_mask = vels_masks[*,*,halpha_ind]
     sn = gf*sqrt(gf_ivar)
     
     vel_mask = (vel_mask gt 0) or (sn lt sncut)

  endif else begin

     vel = mrdfits(mapsfile,'STELLAR_VEL',hdr,/silent)
     vel_ivar = mrdfits(mapsfile,'STELLAR_VEL_IVAR',/silent)
     vel_mask = mrdfits(mapsfile,'STELLAR_VEL_MASK',/silent)
     bin_sn = mrdfits(mapsfile,'BIN_SNR',/silent)
     
     pxscl=abs(sxpar(hdr,'PC2_2')*3600)
     evel=1./sqrt(vel_ivar)
     evel[where(1-finite(evel))]=0
     vel_mask = (vel_mask gt 0) or (bin_sn lt sncut)

  endelse
  
     
  if total(1-vel_mask) lt 10 then begin
     ;to do, come up with more reasonable flag for too much masked data
     logger,err_log,name+' has too much masked data'
     goto,the_end
  endif

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ;; When regridding the algorithm tends to behave poorly near the
 ;; edges of regions with good data.  Let's first run a regridding of
 ;; the data where I simply replace missing regions (i.e. masked
 ;; values) with their nearest neighbor. This keeps future regridding
 ;; when we actually shift the velocity around from going haywire near
 ;; the edges
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;    

  sel=where(1-vel_mask)
  triangulate,reform(xgrid[sel]),reform(ygrid[sel]),triangles
  vel_replace = griddata(reform(xgrid[sel]),reform(ygrid[sel]),reform(vel[sel]),/grid,xout=xgrid[*,0],yout=ygrid[0,*],/nearest_neighbor,triangles=triangles)
  vel_ivar_replace = griddata(reform(xgrid[sel]),reform(ygrid[sel]),reform(vel_ivar[sel]),/grid,xout=xgrid[*,0],yout=ygrid[0,*],/nearest_neighbor,triangles=triangles)
  
  vel = vel_replace
  vel_ivar = vel_ivar_replace
  

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;put onto new grid where optical center is at exactly;;;;;
  ;;;;zero using the coordinates calculated above. Use only;;;;
  ;;;;reliable data for the velocity and error;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     
  sel=where(1-vel_mask,count)
  if count eq 0 then begin
     logger,err_log,'data all masked'
     goto,the_end
  endif

  triangulate,reform(xgrid),reform(ygrid),triangles
  vel_shift = griddata(reform(xgrid),reform(ygrid),reform(vel),/grid,xout=newxarr,yout=newyarr,/linear,triangles=triangles)
  vel_ivar_shift = griddata(reform(xgrid),reform(ygrid),reform(vel_ivar),/grid,xout=newxarr,yout=newyarr,/linear,triangles=triangles)
  triangulate,reform(xgrid),reform(ygrid),triangles
  vel_mask_shift = griddata(reform(xgrid),reform(ygrid),reform(vel_mask),/grid,xout=newxarr,yout=newyarr,triangles=triangles,/nearest_neighbor)

       ;;;;replace old values
     
  vel = vel_shift
  vel_mask = vel_mask_shift
  evel = 1/sqrt(vel_ivar_shift)
  evel[where(1-finite(evel))]=0
  xgrid = rebin(newxarr,n_elements(newxarr),n_elements(newyarr))
  ygrid = transpose(rebin(newyarr,n_elements(newyarr),n_elements(newxarr)))

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;If an additional shift is needed, apply it now.  To Do:
  ;;;;;recalculate shift so that there's only one resampling
  ;;;;;applied.  This will be reduce interpolation errors
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  sz = size(vel,/dim)
  xind = indgen(sz[0])
  yind = indgen(sz[1])
  xout = xind
  yout = yind
  
  if (xshift ne 0) or (yshift ne 0) then begin
  

     xgrid = rebin(xind,sz[0],sz[1])
     ygrid = transpose(rebin(yind,sz[1],sz[0]))
     xout = xind + xshift
     yout = yind + yshift
     
     sel=where(1-vel_mask)
     triangulate,reform(xgrid),reform(ygrid),triangles
     vel_shift = griddata(reform(xgrid),reform(ygrid),reform(vel),/grid,xout=xout,yout=yout,/linear,triangles=triangles)
     evel_shift = griddata(reform(xgrid),reform(ygrid),reform(evel),/grid,xout=xout,yout=yout,/nearest_neighbor,triangles=triangles)
     
     triangulate,reform(xgrid),reform(ygrid),triangles
     vel_mask_shift = griddata(reform(xgrid),reform(ygrid),reform(vel_mask),/grid,xout=xout,yout=yout,triangles=triangles,/nearest_neighbor)

     ;replace old values again
     vel = vel_shift
     vel_mask = vel_mask_shift
     evel = evel_shift
     
  endif 
  
  if total(1-vel_mask) lt 10 then begin
     logger,err_log,name+': too much data masked'
     goto,the_end
  endif

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;Now we calculate the vfield covariance matrix;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  
  sz = size(vel)
  xind = intarr(sz[1],sz[2])
  yind = intarr(sz[1],sz[2])
  for ll = 0,sz[2]-1 do xind[*,ll]=indgen(sz[1])
  for ll = 0,sz[1]-1 do yind[ll,*]=indgen(sz[2])
  
  ind_arr = make_array(sz[1],sz[2],/index)
  ind_1d = reform(ind_arr,sz[1]*sz[2])
  sep = ind_1d*0.
  vel_cor = dblarr(sz[1]*sz[2],sz[1]*sz[2]) ;correlation matrix
  vel_covar = dblarr(sz[1]*sz[2],sz[1]*sz[2])

  for xx=0,max(ind_1d) do begin
     sep=sqrt((xind[xx]-xind[ind_1d])^2+(yind[xx]-yind[ind_1d])^2)
     err_term = evel[xx]*evel[ind_1d] ;*(1-vel_mask[ind_1d])
     vel_cor[ind_1d,xx]=exp(-0.5*(sep/1.9)^2)*(sep lt 6.4)
     vel_covar[ind_1d,xx]=vel_cor[ind_1d,xx]*err_term
  endfor

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;and apply absolute radon transform;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  radon = ds_radon(vel,mask=1-vel_mask,/median,aperture=radon_ap,error=evel,do_covar=do_covar,im_covar=vel_covar,radon_covar = radon_covar)

  radon_norm = ds_radon(vel,mask=1-vel_mask,/median,aperture=radon_ap,error=evel,do_covar=do_covar,im_covar=vel_covar,radon_covar = nradon_covar,/normal)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;run trace, including account for covariance matrix (if keyword;;;;
  ;;;;;;set
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  radon_mask = (1-finite(radon.map)) or (radon.map eq 0) or (radon.maskfrac gt 0) or (radon.length lt radon_ap)
  map = abs(radon.map)
  rho = radon.rho
  theta = radon.theta
  error = radon.error
  mask = radon_mask
  covar = radon_covar
  smo= 0.15

  trace=trace_radon_vm(map,rho,theta,smo=smo,error=error,mask=mask,covar=covar,mc_iter=mc_iter,/silent);,/ploton,/inspect,tolerance=45,fit_range=45.) ;,/silent,/ploton,/inspect)

  nradon_mask = (1-finite(radon_norm.map)) or (radon_norm.map eq 0) or (radon_norm.maskfrac gt 0) or (radon_norm.length lt radon_ap)
  nmap = abs(radon_norm.map)
  nrho = radon_norm.rho
  ntheta = radon_norm.theta
  nerror = radon_norm.error
  nmask = nradon_mask
  ncovar = nradon_covar
  nsmo= 0.15 

  trace_minor=trace_radon_vm(nmap,nrho,ntheta,smo=nsmo,error=enrror,mask=nmask,covar=ncovar,mc_iter=mc_iter,guess_mode='peak',/ploton,/inspect,tolerance=90) ;,/silent,/ploton,/inspect)

  
   
  convert_radon,radon.map,radon_rescale,1,mask=radon_mask

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;create a multi-extension fits file to hold all the information
  ;; ext0 =vel
  ;; ext1 = evel
  ;; ext2 = vel_mask
  ;; ext3 = radon
  ;; ext4 = radon_rescale
  ;; ext5 = radon_mask
  ;; ext6 = error
  ;; ext7 = theta
  ;; ext8 = rho
  ;; ext9 = length
  ;; ext10 = maskfrac
  ;; ext11 = edge
  ;; ext12 = trace

  if file_test(outfile) then spawn,'rm '+outfile
  writefits,outfile,0
  mwrfits,vel,outfile
  mwrfits,evel,outfile
  mwrfits,vel_mask,outfile
  mwrfits,radon.map,outfile
  mwrfits,radon_rescale,outfile
  mwrfits,radon_mask,outfile
  mwrfits,radon.error,outfile
  mwrfits,radon.theta,outfile
  mwrfits,radon.rho,outfile
  mwrfits,radon.length,outfile
  mwrfits,radon.maskfrac,outfile
  mwrfits,radon.edge,outfile
  mwrfits,trace,outfile

  ;add extension names
  extnames = ['velocity','velocity_error','velocity_mask','radon','radon_rescale','radon_mask','radon_error','radon_theta','radon_rho','radon_length','radon_maskfrac','radon_edge','trace']
  for i=1,n_elements(extnames) do begin
     t=mrdfits(outfile,i,hd)
     sxaddpar,hd,'extname',extnames[i-1]
     modfits,outfile,0,hd,exten_no=i
  endfor
  t=mrdfits(outfile,0,hd)
  sxaddpar,hd,'plateifu',db[ii].plateifu
  sxaddpar,hd,',mangaid',db[ii].mangaid
  sxaddpar,hd,'r_ap_arcsec',radon_ap*pxscl
  sxaddpar,hd,'r_ap',radon_ap
  if keyword_set(stars) then sxaddpar,hd,'component','stars' else sxaddpar,hd,'component','Halpha'
  sxaddpar,hd,'xshift',xshift
  sxaddpar,hd,',yshift',yshift
  modfits,outfile,0,hd,exten_no=0

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;now make diagnostic plots;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if keyword_set(plotfile) then begin
     cleanplot,/silent
     set_plot,'ps'
     !p.font=0
     !p.charsize=1
     device,filename=plotfile+'.eps',encapsulated=1,/color,bits=8,xsize=11.69,ysize=8.27,/inch
     pos = cglayout([3,2])      ;,oxmargin=[15,2])

                                ;panel 1, velocity field
     pc = percentiles(vel[where(1-vel_mask)],conflimit=0.95)
     nlevels=15
     dlevel = 2*max(abs(pc))/nlevels
     levels=findgen(nlevels+1)*dlevel - max(abs(pc))
     cgloadct,70,/silent,/reverse

     sz = size(vel,/dim)
     xgrid = rebin(xout,sz[0],sz[1])
     ygrid = transpose(rebin(yout,sz[1],sz[0]))
     cgcontour,vel*(1-vel_mask),xgrid,ygrid,missingvalue=0,/cell_fill,levels=levels,aspect=1,xtitle='x index',ytitle='y index',title='velocity',position=pos[*,2]
     x0=(max(xout)-min(xout))/2+min(xout)
     y0=(max(yout)-min(yout))/2+min(yout)
     xopt = x0-xshift
     yopt = y0-yshift
     cgoplot,[x0,x0],[-1000,1000],linestyle=2
     cgoplot,[-1000,1000],[y0,y0],linestyle=2
     cgoplot,[xopt,xopt],[-1000,1000],color='lime green',linestyle=2
     cgoplot,[-1000,1000],[yopt,yopt],color='lime green',linestyle=2

     xmin = min(xgrid)
     xmax = max(xgrid)
     ymin = min(ygrid)
     ymax = max(ygrid)
     image_size = (xmax-xmin)
     offset_perc = 0.05
     x0 = xmin+offset_perc*image_size
     x1 = xmin+offset_perc*image_size + radon_ap*2
     y0 = ymax - offset_perc*image_size
     y1 = ymax - offset_perc*image_size
     
     plots,[x0,x1],[y0,y1],thick=10,color=cgcolor('orange')
     
     drawifu,strtrim(string(db[ii].IFUDESIGNSIZE),2),xoff=xopt,yoff=yopt,pxscl=pxscl,coords=coords

     sel=where((trace.status gt 0 and 1-trace.bias_flag),count)
     trace_x = trace.rho*cos(trace.theta) + max(xgrid)/2
     trace_y = trace.rho*sin(trace.theta) + max(ygrid)/2
     cgoplot,trace_x[sel],trace_y[sel],psym=4,symsize=0.5,thick=3,color='magenta'
     
     sel=where((trace_minor.status gt 0 and 1-trace_minor.bias_flag),count)
     trace_x = trace_minor.rho*cos(trace_minor.theta) + max(xgrid)/2
     trace_y = trace_minor.rho*sin(trace_minor.theta) + max(ygrid)/2
     cgoplot,trace_x[sel],trace_y[sel],psym=4,symsize=0.5,thick=3,color='turquoise'

     
     sdss_pxscl = 0.396127
     height = 0
     buffer = abs(xgrid[0,0] - xgrid[-1,0])*pxscl*0.4
     while height lt 64 do begin
        height = round((max(ygrid*pxscl+buffer) - min(ygrid*pxscl-buffer))/sdss_pxscl)
        width = round((max(xgrid*pxscl+buffer) - min(xgrid*pxscl-buffer))/sdss_pxscl)
        if height lt 64 then sdss_pxscl = sdss_pxscl/2.
        endwhile
     
     sdss_ra = db[dbind].objra
     sdss_dec = db[dbind].objdec
     url = '"http://skyserver.sdss.org/dr13/SkyServerWS/ImgCutout/getjpeg?ra='+strtrim(string(sdss_ra),2)+'&dec='+strtrim(string(sdss_dec),2)+'&scale='+strtrim(string(sdss_pxscl),2)+'&width='+strtrim(string(width),2)+'&height='+strtrim(string(height),2)+'"'
;        print,url
     jpegname = (strsplit(url,'?',/extract))[-1]
     jpegname = (strsplit(jpegname,'"',/extract))[0]
     jpegname = './cutouts/'+name+'_'+jpegname+'.jpg'
     file_info = file_info(jpegname)

                                ;if cutouts directory doesnt exit, create
     if 1-file_test('./cutouts/',/directory) then spawn,'mkdir cutouts'
     ;if already exists don't redownload
     if file_info.size eq 0 then spawn,'wget -O "'+jpegname+'" '+url
     file_info = file_info(jpegname)
     if file_info.size gt 0 then begin
        READ_JPEG, jpegname, a, TRUE=1
        cgimage,a,/noerase,position=pos[*,1],/keep_aspect_ratio,/axes,xtitle='X',ytitle='Y',color=cgcolor('black'),title='SDSS gri'
     endif
     imsz = size(a,/dim)
     x0 = (imsz[1]-1)/2
     y0 = (imsz[2]-1)/2
     drawifu,strtrim(string(db[ii].IFUDESIGNSIZE),2),xoff=x0,yoff=y0,pxscl=sdss_pxscl,coords=coords

                                ;normal radon map;;;
     sz = size(radon_rescale,/dim)
     thetagrid = rebin(radon.theta,sz[0],sz[1])
     rhogrid = transpose(rebin(radon.rho,sz[1],sz[0]))

     cgcontour,abs(radon_norm.map)*(1-nradon_mask),thetagrid/!dtor,rhogrid,nlevels=15,/cell_fill,position=pos[*,3],/noerase,missingvalue=0,xtitle='theta [deg]',ytitle='rho',title='RT'
     sel=where((trace_minor.status gt 0 and 1-trace_minor.bias_flag),count)
     cgoplot,(trace_minor.theta[sel]*(trace_minor.theta[sel] gt 0 and trace_minor.theta[sel] lt !pi) + (trace_minor.theta[sel] + !pi)*(trace_minor.theta[sel] lt 0) + (trace_minor.theta[sel]-!pi)*(trace_minor.theta[sel] gt !pi))/!dtor,trace_minor.rho[sel],psym=4,color='turquoise',symsize=0.5,thick=2
     
     ;;;;radon map;;;;
     sz = size(radon_rescale,/dim)
     thetagrid = rebin(radon.theta,sz[0],sz[1])
     rhogrid = transpose(rebin(radon.rho,sz[1],sz[0]))
     
     cgcontour,radon_rescale*(1-radon_mask),thetagrid/!dtor,rhogrid,nlevels=15,/cell_fill,position=pos[*,4],/noerase,missingvalue=0,xtitle='theta [deg]',ytitle='rho',title='RT'
     cgcolorbar,range=[0,1],/vertical,/right,position=[pos[2,4]+0.01,pos[1,4],pos[2,4]+0.03,pos[3,4]]
     sel=where((trace.status gt 0 and 1-trace.bias_flag),count)
     cgoplot,(trace.theta[sel]*(trace.theta[sel] gt 0 and trace.theta[sel] lt !pi) + (trace.theta[sel] + !pi)*(trace.theta[sel] lt 0) + (trace.theta[sel]-!pi)*(trace.theta[sel] gt !pi))/!dtor,trace.rho[sel],psym=4,color='magenta',symsize=0.5,thick=2

     
     ;;trace

     rho_arr = trace.rho
     theta_arr = trace.theta
     etheta_arr = trace.etheta
     bias_flag_arr = trace.bias_flag
     status_arr = trace.status
     
     ;; unwrap_pa,rho_arr,theta_arr,srtarr=srtarr
     ;; bias_flag_arr = bias_flag_arr[srtarr]
     ;; status_arr = status_arr[srtarr]
     ;; etheta_arr = etheta_arr[srtarr]
     
     sel=where(trace.status ne 0 and 1-trace.bias_flag)
     sel=where(status_arr gt 0 and 1-bias_flag_arr)
     cgplot,rho_arr,theta_arr/!dtor,psym=3,err_ylow = etheta_arr/!dtor,err_yhigh=etheta_arr/!dtor,xtitle='rho',ytitle='theta',position=pos[*,5],/noerase,/err_clip,/ynozero,title='trace',yrange=minmax(theta_arr[sel]/!dtor)+[-10,10]
     sel=where(1-(status_arr gt 0 and 1-bias_flag_arr),count)
     if count gt 0 then cgoplot,rho_arr[sel],theta_arr[sel]/!dtor,psym=7,color=cgcolor('red'),thick=2,/ynozero
     

                                ;write a few useful numbers
        sep = 0.05
        xtoff = -0.0
        xyouts,pos[0,0]+xtoff,pos[3,0]-0.1,'PLATE-IFU: '+strtrim(string(name),2),charsize=1.5,/norm,color=cgcolor('black')
        xyouts,pos[0,0]+xtoff,pos[3,0]-0.1-sep*1,'MaNGA ID: '+strtrim(string(db[ii].mangaid),2),charsize=1.5,/norm,color=cgcolor('black')
        xyouts,pos[0,0]+xtoff,pos[3,0]-0.1-sep*2,'log M!d*!n='+strtrim(string(alog10(db[ii].nsa_sersic_mass),F='(F6.2)'),2),/norm,color=cgcolor('black'),charsize=1.5
        xyouts,pos[0,0]+xtoff,pos[3,0]-0.1-sep*3,'Radon R!dap!n (")='+strtrim(string(radon_ap*pxscl,F='(F6.2)'),2),/norm,color=cgcolor('black'),charsize=1.5
        xyouts,pos[0,0]+xtoff,pos[3,0]-0.1-sep*4,'RA='+strtrim(string(db[ii].objra,F='(F8.4)'),2),/norm,color=cgcolor('black'),charsize=1.5
        xyouts,pos[0,0]+xtoff,pos[3,0]-0.1-sep*5,'DEC='+strtrim(string(db[ii].objdec,F='(F8.4)'),2),/norm,color=cgcolor('black'),charsize=1.5

        if keyword_set(stars) then ststr = 'Component: stars' else ststr = "Component: gas" ;else ststr="Component:stars"
        xyouts,pos[0,0]+xtoff,pos[3,0]-0.1-sep*6,ststr,/norm,color=cgcolor('black'),charsize=1.5

        device,/close

                                ;comment these lines if you
                                ;don't want to convert the eps
                                ;file to a png file
        spawn,'convert -density 125 '+plotfile+'.eps'+' '+plotfile+'.png'
;        spawn,'rm '+plotfile+'.eps'
        
        set_plot,'x'
        
      
  endif
  output=1
  
        
the_end:print,' '
  
close,1

  
end

;; testing

;; db=mrdfits('/home/dstark/manga/fits/drpall-v2_5_3.fits',1)

;; name='8082-6102'
;; dbind=where(db.plateifu eq name)
;; step=1;db[dbind].nsa_elpetro_th50_r/4/0.5 ;1
;; xshift=0;-0.5;1.75;.25
;; yshift=0;1.25
;; radon_ap=db[dbind].nsa_elpetro_th50_r/0.5*db[dbind].nsa_elpetro_ba
;; plotfile='blah'
;; sncut=3
;; filter_radon=0
;; usediff=0
;; stars=1
;; smooth=0
;; do_covar=1
;; radon_wrapper_working,name,output,'test_radon_working.fits',xshift=xshift,yshift=yshift,radon_ap=radon_ap,err_log='test_radon_log.log',do_covar=do_covar,eh_on=0,plotfile='test_plotfile',/stars,mapstype='SPX',mapspath='~/manga/fits/mpl8/dap/SPX/'

;; end

