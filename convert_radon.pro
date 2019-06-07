pro convert_radon,input,output,mode,mask=mask,error=error,output_error=output_error

                                ;converts radon transform into
                                ;alternate form (e.g., renormalized,
                                ;inverted)

                                ;mode: (1) normalized radon transform
                                ;(2) inverted, renormalized radon transform
  

  
  if keyword_set(mask) then input[where(mask)]=-1
  sz = size(input)

  output = fltarr(sz[1],sz[2])
  output_error = fltarr(sz[1],sz[2])
  
  for mm=0,sz[2]-1 do begin
     row = input[*,mm]
     if keyword_set(error) then row_error = error[*,mm]
                                ;set to range from 0 --> 1
     goodpix = where(row ge 0,count)
     if count gt 1 then begin
        minval = min(row[goodpix])
        outrow = row - minval
        maxval = max(outrow[goodpix])
        outrow = outrow/maxval

        if keyword_set(error) then outrow_error = row_error/maxval

        if mode eq 2 then outrow = 1-outrow
        
        output[*,mm]=outrow
        if keyword_set(error) then output_error[*,mm] = outrow_error

        
     endif
  endfor
  
                                ;make sure all masked values are still masked  
  if keyword_set(mask) then output[where(mask)]=-1
end

  
