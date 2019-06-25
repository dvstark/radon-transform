function get_chunks,arr,mask,nchunks=nchunks

                                ;identifies indices of chuncks of good
                                ;data. To be used with trace_radon_vm
  sel=where(mask,count)
  if count eq 0 then begin

     nchunks=1
     return,[0,n_elements(arr)-1]
     
  endif else if count eq n_elements(mask) then begin

     nchunks=0
     return,[-1,-1]
     
  endif else begin

     arr = [arr[0],arr,arr[-1]]
     mask = [mask[0],mask,mask[-1]]

     cond1 =  abs(mask - shift(mask,1)) eq 0
     cond2 = abs(mask - shift(mask,-1)) eq 0

     arr = arr[1:-2]
     mask = mask[1:-2]
     cond1=cond1[1:-2]
     cond2=cond2[1:-2]

     ind1 = where(cond1 and not cond2 and not mask)
     ind0 = where(cond2 and not cond1 and not mask)   
     
     if mask[0] eq 0 then ind0 = [0,ind0]
     if mask[-1] eq 0 then ind1 = [ind1,n_elements(mask)-1]
     
     ind0=ind0[where(ind0 ne -1)]
     ind1=ind1[where(ind1 ne -1)]

     nchunks=n_elements(ind0)
     return,transpose([[ind0],[ind1]])
          
  endelse
  

end
