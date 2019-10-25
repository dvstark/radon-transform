function find_zeros,x,y

                                ;simple function to find where a data
                                ;set crosses zero. Useful in cases
                                ;where there may be multiple zero
                                ;crossings. Meant for smooth data.
  xl = x[0:-2]
  xh = y[1:-1]
  yl = y[0:-2]
  yh = y[1:-1]

  dsign = yl*yh lt 0
  sel=where((dsign eq 1) or (yl eq 0),count)
;  print,count,' zero crossings found'

  if count eq 1 then begin
     zero = interpol(x,y,0)
  endif else begin
     sec = xl*0
     for i=0,count-1 do begin
        if i eq 0 then sec[0:sel[i]] = i
        if i gt 0 and i le (count-1) then sec[sel[i-1]+1:sel[i]]=i
        if (i eq (count-1)) and (sel[i] ne n_elements(sec)-1) then sec[sel[i]+1:-1]=i+1
     endfor

     zero = []
     for i=0,count-1 do begin
        sel=where(sec eq i or sec eq i+1)
        z = interpol(xl[sel],yl[sel],0)
        zero = [zero,z]
     endfor
     if n_elements(zero) eq 0 then zero=[-999]

;     print,''
  endelse

;  cgoplot,zero,intarr(n_elements(zero)),psym=4,color='red',symsize=2,thick=3
;  print,'zeros at: ',zero
  return,zero

end

;; x=findgen(101)-50
;; y=(x)^3+8*x^2-x*5-20

;; plot,x,y,xrange=[-10,5]
;; oplot,[-1000,1000],[0,0]

;; z = find_zeros(x,y)

;; end
  
