pro drawifu,ifudsgn,noplot=noplot,coords=coords,xoff=xoff,yoff=yoff,pxscl=pxscl

  if 1-keyword_set(xoff) or 1-keyword_set(yoff) then begin
     xoff = 0
     yoff = 0
  endif

  if 1-keyword_set(pxscl) then pxscl = 0.5
  
  
                                ;simple code to draw the MaNGA IFU
                                ;outline.  Thanks go to Edmond Cheung
                                ;for sending this
  
  case ifudsgn of 
     '19': r_arc=12.5/2.
     '37': r_arc=17.5/2.
     '61': r_arc=22.5/2.
     '91': r_arc=27.5/2.
     '127':r_arc=32.5/2
     else: begin
        print,'Invalid IFU design provided'
        r_arc=0
        return
     end
     
  endcase
  
;this statement plots the polygon over the image. This assumes that the image is centered at 0, 0, and the axis are in arcseconds. So make sure that each spaxel has the right dimensions. 

  outx=replicate(0, 7) + [r_arc/2,r_arc, r_arc/2, -r_arc/2, -r_arc,-r_arc/2, r_arc/2] ;+ xoff

  outy=replicate(0, 7) + [sqrt(3)*r_arc/2., 0, -sqrt(3)*r_arc/2., -sqrt(3)*r_arc/2.,0, sqrt(3)*r_arc/2.,sqrt(3)*r_arc/2.] ;+ yoff

  if keyword_set(pxscl) then begin
;     print,'adjusting by pxscl = ',pxscl
     outx = outx/pxscl
     outy = outy/pxscl
  endif

  outx = outx + xoff
  outy = outy + yoff
  
    coords = [[outx],[outy]]
  
  
;  if 1-keyword_set(noplot) then cgPolygon, replicate(0, 7) + [r_arc/2, r_arc, r_arc/2, -r_arc/2, -r_arc, -r_arc/2, r_arc/2], replicate(0, 7) + [sqrt(3)*r_arc/2., 0, -sqrt(3)*r_arc/2., -sqrt(3)*r_arc/2., 0, sqrt(3)*r_arc/2., sqrt(3)*r_arc/2.], COLOR='magenta',thick=2,/data
  if 1-keyword_set(noplot) then cgPolygon, outx, outy, COLOR='magenta',thick=2,/data
end
