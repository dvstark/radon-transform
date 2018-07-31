;;;The following is a simple script showing how to calculate the
;;;absolute Radon transform as described in Stark et al. 2018. This
;;;example uses the IDL Astro Library and the Coyote Library

;Let's generate a simple model velocity field for this example.

;some parameters:
dim = [71,71]       ;dimensions
x0=(dim[0]-1)/2.    ;center in x direction
y0=(dim[1]-1)/2.    ;center in y direction
inclin=30.          ;inclination
PA=45               ;position angle
v0=200.             ;rotation velocity (in flat part of rotation curve)
re = (dim[0]-1)/4   ;effective radius
h_rot=re/2.         ;scale length of rotation curve (vrot = v0*tanh(r/h_rot))

;create a grid of the true galacto-centric radius
b_a = cos(inclin*!dtor) ;axis ratio assuming an infinitely thin disk
dist_ellipse,R_true,dim,x0,y0,1/b_a,PA

;calculate cosine of phi (angle relative to PA in plane of sky)
xarr = findgen(dim[0])-x0
yarr = findgen(dim[1])-y0
xgrid = rebin(xarr,dim[0],dim[1]) ;x coordinate in plane of sky
ygrid = transpose(rebin(yarr,dim[1],dim[0])) ;y coordinate in plane of sky
rgrid = sqrt(xgrid^2+ygrid^2) ;radius in plane of sky

phi=atan(ygrid/xgrid)
phi[where(xgrid lt 0)]=phi[where(xgrid lt 0)]+!pi
phi[where(xgrid gt 0 and ygrid lt 0)] = phi[where(xgrid gt 0 and ygrid lt 0)]+2*!pi
phi[where(xgrid eq 0 and ygrid lt 0)]=3*!pi/2.
cosphigrid = cos(phi-PA*!dtor)

;calculate costheta (angle with respect to PA in plane of galaxy)
costheta = rgrid*cosphigrid/r_true

;true rotation velocity assuming a typical rotation curve
v_rot = v0*tanh(r_true/h_rot)

;Now let's plot the velocity field, calculate the RT (with and
;without bounds), and plot the results

window,0,xsize=900,ysize=300
!p.multi=[0,3,1]
;observed velocity
v_obs = sin(inclin*!dtor)*(v0*tanh(r_true/h_rot)*costheta)
v_obs[x0,y0]=0 ;fixing that point in the middle
cgloadct,70,/reverse
cgcontour,v_obs,xgrid,ygrid,xrange=reverse(minmax(xgrid)),yrange=minmax(ygrid),nlevels=15,/fill,xtitle='X',ytitle='Y'
;(I've reversed the x axis to follow astro convention)

;calculate RT with no bounds on the integral and plot

rt = ds_radon(v_obs)
cgloadct,68,/silent
rt_sz = size(rt.map,/dim)
thetagrid = rebin(rt.theta,rt_sz[0],rt_sz[1])
rhogrid = transpose(rebin(rt.rho,rt_sz[1],rt_sz[0]))

;In this case, the RT will depend not only on the velocities
;themselves, but how many data points are integrated for a given
;[theta,rho].  We can find and flag any points outside our input
;velocity field by finding where LENGTH=0.

rt.map[where(rt.length eq 0)]=-9999
cgcontour,rt.map,thetagrid/!dtor,rhogrid,nlevels=15,/fill,xrange=minmax(thetagrid/!dtor),yrange=minmax(rhogrid),xtitle='theta (degrees)',ytitle='rho',missingvalue=-9999

;do same thing, but place bounds on the integral

rt_ap = ds_radon(v_obs,aperture=10)

;The use of an aperture means the RT will not depend on the number of
;pixels being integrated over (which, aside from slight variations due
;to finite pixel sizes, is constant) except near the edges. We can
;[theta,rho] values where the integration extended beyond the edge of
;our map using the EDGE parameter.

rt_ap.map[where(rt_ap.edge gt 0)]=-9999
cgcontour,rt_ap.map,thetagrid/!dtor,rhogrid,nlevels=15,/fill,xrange=minmax(thetagrid/!dtor),yrange=minmax(rhogrid),xtitle='theta (degrees)',ytitle='rho',missingvalue=-9999

;The extra noise is due to slight variations in # of pixels being
;integrated over at a given [theta,rho].

;reset the color table
loadct,0,/silent

end
