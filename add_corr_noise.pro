function add_corr_noise,vel,covar,return_covar = return_covar
  sz = size(vel,/dim)
  if n_elements(sz) eq 2 then vel_inds = make_array(sz[0],sz[1],/index) $
  else vel_inds = make_array(sz[0],/index)
  
  
  csz = size(covar,/dim)
  csz=csz[0]
  
                                ;test to see if the covariance matrix
                                ;is positive definite
  covar_sub_adj = covar
  covar_inds = make_array(csz[0],csz[0],/index)

  ;removing empty rows
  goodrow = []
  for i=0,csz-1 do if covar[i,i] ne 0 then goodrow = [goodrow,i]
  covar_sub = covar[goodrow,*]
  covar_sub = covar_sub[*,goodrow]
  covar_inds = covar_inds[goodrow,*]
  covar_inds = covar_inds[*,goodrow]
  
  ;test to see if positive definite
  test = covar_sub
  la_choldc,test,status=status
  ;print,'status check: ',status
  if status gt 0 then begin
     print,'WARNING: covariance matrix not positive definite! Removing diagonal values set to zero and adding small factor into remaining diagonal elements'

                                ;now add some value to diagonal to make positive def
     eps = 0.00
     step=0.01
     status=999
     subsz = size(covar_sub,/dim)
     inds=indgen(subsz[0])
     
     while status gt 0 do begin
        print,'eps = ',eps
        covar_sub_adj = covar_sub
        covar_sub_adj[inds,inds]=(1+eps)*covar_sub_adj[inds,inds]
        test = covar_sub_adj
        
        la_choldc,test,status=status
        ;print,'status: ',status
        if status gt 0 then eps = eps+step
     endwhile
  endif else covar_sub_adj = covar_sub
     
                                ;now get new velocity values
  offset = mrandomn(seed,covar_sub_adj)
  newvel = vel
  newvel[goodrow] = newvel[goodrow]+offset

  return_covar = covar;covar_sub_adj
  return_covar[reform(covar_inds)] = reform(covar_sub_adj)
  
  return,newvel
  
end
