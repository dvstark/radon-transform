;example run of radon_wrapper.pro

drpallpath = '/home/david-stark/manga/fits/drpall-v2_5_3.fits'
db=mrdfits(drpallpath,1)
name = '8082-3702'
dbind=where(db.plateifu eq name)
pxscl=0.5                       ;pixel scale

output_file = 'test_radon_SPX.fits'
xshift=0
yshift=0

;I typically define the radon aperture as r_e*b_a. Needs to be in spaxels
radon_ap=db[dbind].nsa_elpetro_th50_r/pxscl*db[dbind].nsa_elpetro_ba

err_log='test_radon_HBY10.log'
do_covar = 0        ;let's not calculate the full covariance matrix just yet...

eh_on = 0                       ;set to 1 to turn on error handling
plotfile='test_radon_SPX'           ;final output file will be test_radon.png
stars=1                         ;will analyze stellar velocity field
mapstype = 'HYB10'              ;'SPX'
mapspath = './'   ;'~/manga/fits/mpl8/dap/SPX/'
maps_suffix ='GAU-MILESHC'  ;'MILESHC-MILESHC'

radon_wrapper,name,output,output_file,xshift=xshift,yshift=yshift,radon_ap=radon_ap,err_log=err_log,do_covar=do_covar,eh_on=0,plotfile=plotfile,stars=stars,mapstype=mapstype,mapspath=mapspath,mc_iter=0,drpallpath='/home/david-stark/manga/fits/drpall-v2_5_3.fits',maps_suffix=maps_suffix

end
