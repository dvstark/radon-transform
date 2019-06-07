Function sdep, DS=ds, PS=ps, VF=vf, VS=vs, W=w

;+
; NAME:
;	SDEP
;
; PURPOSE:
;	Returns several System DEPendent parameters, like OS family (default), 
;	directory separator, path separator, etc. 
;
; CATEGORY:
;	General utilities
;
; CALLING SEQUENCE:
; 
;	Result = SDEP()
;
; INPUTS:
;	None
;
; KEYWORD PARAMETERS:
;	DS:	(Directory Separator) When set, SDEP returns the directory 
;		separator (i.e. 
;		/ (slash) under Unix or \ (backslash) under Windows).
;	PS:	(Path Separator) When set, SDEP returns the path separator 
;		(i.e. :  under Unix or ; (backslash) under Windows).
;	VS:	(Version Short) returns the idl short-version (i.e. 
;		'5' for idl 5.0.2)
;	VF:	(Version Full) returns the idl full versiom (i.e. 
;		'5.0.2' for idl 5.0.2)
;	W:	(Widget allowed?) Returns 1 when widget are allowed (i.e. 
;		!d.nane is 'WIN','X' or 'MAC') otherwise returns 0 (i.e. 'PS').
;
; OUTPUTS:
;	This function returns (with no keywords) the !version.os_family
;	value in uppercase. If any keyword is set, the returned value
;	is changed to the one described below.
;
; RESTRICTIONS:
;	Never used under Mac
;
; PROCEDURE:
;	Straightforward
;
; EXAMPLE:
;	print,sdep()
;
; MODIFICATION HISTORY:
; 	Written by:	srio@esrf.fr and dejus@aps.anl.gov
;	Sept, 1997	
;	97/10/16 srio@esrf.fr adds /w keyword.
;	98/12/23 srio@esrf.fr tentative update for Mac
;-
; returns system dependent values
; returns osversion if no keywords are set, otherwise returns the
; requested separator (ds = directory separator, ps = path separator)
; srio@esrf.fr 97-09-15 and dejus@aps.anl.gov 09/15/97.

osversion = StrUpCase(!version.os_family)

if keyword_set(ds) then begin		; directory separator
  CASE osversion OF
    'UNIX': 	return,'/'
    'WINDOWS': 	return,'\'
    'MACOS': 	return,':'
    else:		return,''
  ENDCASE
endif

if keyword_set(ps) then begin		; path separator
  CASE osversion OF
    'UNIX': 	return,':'
    'WINDOWS':	return,';'
    'MACOS': 	return,','
     else:		return,''
  ENDCASE
endif

if keyword_set(vf) then return,!version.release

if keyword_set(vs) then return,strmid(!version.release,0,1)

if keyword_set(w) then begin		; path separator
  CASE !d.name OF
    'WIN': 	return,1
    'X':	return,1
    'MAC': 	return,1
     else:	return,0
  ENDCASE
endif

return, osversion
end ; sdep

FUNCTION Make_Set,x,y1,y2,y3,y4,y5,y6,y7,y8, ToFile=tofile, group=group, $
  _Extra=extra

;+
; NAME:
;	MAKE_SET
;
; PURPOSE:
;       This function returns a multicolumn array (array(ncols,npoints))
;	build from separate arrays with the column data.
;
; CATEGORY:
;	General utilities
;
; CALLING SEQUENCE:
;       Result = Make_Set(x,y1 [,y2,...,y8])
;
; INPUTS:
;       
;	x: The array with abscissas
;	y1: The array with ordinates (must be of the same type as x and
;		with the same number of points)
;
; OPTIONAL INPUTS:
;	y2,...,y8 additional ordinates arrays 
;	
; KEYWORD PARAMETERS
;	ToFile: a string cointaining a file name where to write
;		(optionally) the data. If this strinf is '?' then 
;		the Dialog_PickFile() window is started to receive
;		the file name.
;	Group: The widget if of the caller.
;	_Extra: any other keyword to ba passed to Dialog_PickFile()
;
; OUTPUTS:
;       This function return the global array and optionally writes 
;	a disk file.
;
; PROCEDURE:
;       
;       Easy.
;
; EXAMPLE:
;	help,Make_Set(FindGen(100),FindGen(100))
;	<Expression>    FLOAT     = Array[2, 100]
;
; MODIFICATION HISTORY:
; 	Written by:	Manuel Sanchez del Rio (srio@esrf.fr), 98-12-21
;
;-

Catch, error_status
IF error_status NE 0 THEN BEGIN
   Message,/Info,'error caught: '+!err_string
   itmp = Dialog_Message(/Error,Dialog_Parent=group, $
     'MAKE_SET: error caught: '+!err_string)
   Catch, /Cancel
   On_Error,2
   IF N_Elements(out) NE 0 THEN RETURN,out ELSE RETURN,0
ENDIF

nn = N_Params()
nx = N_Elements(x)
IF nn LT 2 THEN Message,'Usage: result = Make_Set(x,y [,y2,y3,y4,y5,y6,y7,y8])'
CASE nn OF
 2: out = reform([x,y1],nx,nn)
 3: out = reform([x,y1,y2],nx,nn)
 4: out = reform([x,y1,y2,y3],nx,nn)
 5: out = reform([x,y1,y2,y3,y4],nx,nn)
 6: out = reform([x,y1,y2,y3,y4,y5],nx,nn)
 7: out = reform([x,y1,y2,y3,y4,y5,y6],nx,nn)
 8: out = reform([x,y1,y2,y3,y4,y5,y6,y7],nx,nn)
 9: out = reform([x,y1,y2,y3,y4,y5,y6,y7,y8],nx,nn)
 else: Message,'Bad inputs'
ENDCASE
out = Transpose(out)
IF Keyword_Set(tofile) THEN BEGIN
  IF StrCompress(tofile,/Rem) EQ '?' THEN BEGIN
    tofile = Dialog_PickFile(/Write,group=group,file='tmp.dat', $
	_Extra=extra)
    IF tofile EQ '' THEN RETURN,out
  ENDIF
  Openw,unit,tofile,/Get_Lun
  FOR i=0L,nx-1 DO PrintF,unit,out[*,i]
  Free_Lun,unit
  Message,/Info,'File '+tofile+' written to disk.'
  IF SDep(/W) THEN $
     itmp = Dialog_Message(Dialog_Parent=group, $
	/Info,'File '+tofile+' written to disk.')
ENDIF
RETURN,out
END

FUNCTION PeakFinder,y,x, Group=group, Sort=sort, $
  PCutoff=pcutoff, CLimits=cLimits, NPeaks=nPeaks, Error=error, $
  Optimize=optimize,widget_off=widget_off,silent=silent

;+
; NAME:
;	PEAKFINDER
;
; PURPOSE:
;	This function finds the peaks of a 1-D data. 
;
;       It returns a multicolumn array with the following data:
;	Column 0: The indices of the peaks.
;	Column 1: The abscissa values of the found peaks
;	Column 2: The ordinate values of the found peaks.
;	Column 3: A weight indicating how "important" is the peak.
;	Column 4: The "normalized peak weight" (i.e., col 3 normalized 
;		between 0 and 1.
;	Column 5: The "peak width" in number of points.
;
;	If no peaks are found, it returns 0.
;
; CATEGORY:
;	Maths.
;
; CALLING SEQUENCE:
;	result = PeakFinder(y[,x])
;
; INPUTS:
;	y:	  An array with input (ordinates) data.
;
; OPTIONAL INPUTS:
;	x:	  An array with input abscissas data (must be of
;		the same dimension of y)
;
; KEYWORD PARAMETERS:
;
;	INPUTS: 
;	Group:	The widget id of the caller. This is used to center
;	 	the Dialog_Message window(s).
;	Sort:	If set, the returned array is sorted by weight. By
;		default is it sorted by index.
;	Optimize: If set, for each found peak looks if the immediate
;		leaft and right neighmourhood is higher. If so,
;		substitute by this value.
;
;
;	OUTPUTS:
;	PCutoff: A named variable where the data (FltArr(2,npts))
;		 of the number of peaks versus cutoff value is returned.
;	CLimits:	A named variable where an array is returned. It
;		contains [cutoff_a,cutoff_w,cutoff_1,cutoff_2]
;		being: 
;		cutoff_1: The minimum cutoff value that produces 
;			a "stabilized" number of peaks. 
;		cutoff_2: The maximum cutoff value that produces 
;			a "stabilized" number of peaks. 
;		cutoff_a: The average value 0.5*(cutoff_1+cutoff_2)
;		cutoff_w: The width value cutoff_2-cutoff_1
;	NPeaks:  A named variable where the "stabilized" or "reasonable"
;		number of peaks is stored. 
;	Error:  A named variable where to put error flag (0=No Error,
;		1=Error). (To find zero peaks is not considered an error.)
;	
; OUTPUTS:
;	A multicolumn array. 
;
; SIDE EFFECTS:
;
; PROCEDURE:
;	This IDL procedure implements a peak search procedure that 
;	I invented. 
;
;       A peak is a local maximum  in your data. A local maximum 
;	has derivative zero. I calculate the derivative of the y
;	data (using IDL's Deriv) and I consider "a peak" any 
;	value which derivative is positive (or zero) and the 
;	point after has a negative derivative. All peaks are returned.
;	In the case thet "Optimize" keyword is set, each peak is checked
;	against its immidiate left and right neighbour. If the neighbour is 
;	larger, we take the neighbour instead.
;
;	However, real data contains noise and the returned peaks
;	are much more than the "intrinsic", "reasonable" or supposely
;	"real" peaks. The routine performs some kind of evaluation
;	of the peaks to see is the peaks is "good" or not. For that
;	purpose we assign a weight (col 3) value to each peak. The larger
;	is the weight, the peak is more "important". 
;
;	How to get the weight: looking at the derivetive of some
;	data one realizes that the in the neighbourhood of the 
;	"important" peaks the absolute values of the derivative are
;	larger. In addition, an "important" peak has meny points 
;	with positive derivtive to its left and many points with 
;	negative derivative to its right. We evaluate, for each peak,
;	the number of points that have positive derivative at its 
;	left plus the number of points with negative derivative
;	at its right. This is the "peak width" (col 5).
;	The sum of the absolute values of the derivatives for all
;	the points inside the "peak width" constitutes the "peak
;	weight". The "normalized peak weight" is the peak weight
;	divided by the maximum peak weight. 
;	The "peak width", "peak weight" and "normalized peak weight"
;	are returned for eack peak, in columns 5, 3 and 4, respectively.
;
;	One can reduce the number of peaks setting a "cutoff" value,
;	i.e., considering only the peaks which normalized weight
;	is larger than this cutoff value. 
;	What we do next is calculate the number of peaks as a function
;	of the cutoff value. This information is optionally returned 
;	in the PCutoff array. 
;	From this plot one can estimate the number of peaks that 
;	are "reasonable". I obtain this when the plot "stabilizes",
;	i.e, I get the same number of points for a wide internal of
;	cutoff values. The minimum and maximum cutoff values that 
;	provide with the "stabilized" number of points are optionally 
;	returned in the CLimits keyword. The difference between 
;	these two numbers (CLimits[1]) is an estimator of the 
;	confidence of the "resonable" of "stabilized" number of peaks
;	(nPeaks). The "reasonable" or "stabilized" number of points is 
;	returned in the NPeaks keyword.
;
;	The "reasonable" NPeaks are easily obtainmed from the NPeak
;	first lines of the returned array when the Sort keyword is used.
;	(See example).
;	 
; EXAMPLE: 
;	;This small program illustrates the use of PeakFinder
;
;	;
;	; create some data
;	;
;	nn=200
;	x=FIndGen(nn)/Float(nn-1)*4
;	x=x-2
;	
;	y0 = 0.11*x + 1.0
;	y1 = Voigt1(x,[10,-0.5,0.1,0.5])
;	y2 = Voigt1(x,[10,0,0.2,0.3])
;	y3 = Voigt1(x,[10,0.5,0.4,0.0])
;	y=y0+y1+y2+y3
;	yran=max(y) - min(y)
;	yran = (yran * 0.05)*randomu(seed,nn)
;	y=y+yran
;	
;	;
;	; calls PeakFinder
;	;
;	pcutoff=0 & CLimits=0 & npeaks=0
;	a=PeakFinder(y,x,PCutoff=pcutoff,CLim=CLimits,NPeaks=npeaks,/Sort,/Opt)
;	
;	;
;	; Display/plot results
;	;
;	Print,'Peaks found: ',N_Elements(a[0,*])
;	Print,'Good Peaks: ',nPeaks
;	Print,'cutoff value: ',climits[0]
;	Print,'confidence value: ',climits[1]
;	Plot,x,y,linestyle=1
;	OPlot,a[1,*],a[2,*],PSym=1
;	OPlot,a[1,0:npeaks-1],a[2,0:npeaks-1],PSym=6
;	
;	;
;	; This part overplots Lorentzian functions estimated from 
;	; the PeakFinder results
;	;
;	pause
;	FOR i=0,nPeaks-1 DO BEGIN
;	  ypeak = a[2,i]
;	  xpeak = a[1,i]
;	  imin = a[0,i]+(a[5,i]/2)
;	  imax = a[0,i]-(a[5,i]/2)
;	  xfwhm = x[imax]-x[imin]
;	  OPlot,x,Voigt1(x,[ypeak,xpeak,xfwhm,1])
;	ENDFOR
;	
;	;
;	; plot "number of peaks" versus "cutoff value"
;	;
;	pause
;	Plot,pcutoff[0,*],pcutoff[1,*],PSym=10,XTitle='cutoff value',$
;	  YTitle='Number of peaks'
;	END
;	
; MODIFICATION HISTORY:
; 	Written by:	M. Sanchez del Rio, srio@esrf.fr
;	22 January , 1999
;-

Catch, error_status
IF error_status NE 0 THEN BEGIN
   Message,/Info,'error caught: '+!err_string
   IF SDep(/w) THEN itmp = Dialog_Message(/Error,Dialog_Parent=event.top, $
     'PEAKFINDER: error caught: '+!err_string)
   error=1
   Catch, /Cancel
   On_Error,2
   RETURN,0
ENDIF

if 1-keyword_set(silent) then silent=0

error=0
IF N_Elements(y) EQ 0 THEN BEGIN
  nn=200
  x=FIndGen(nn)/Float(nn-1)*4
  x=x-2
  y0 = 0.11*x + 1.0
  y1 = Voigt1(x,[10,-0.5,0.1,0.5])
  y2 = Voigt1(x,[10,0,0.2,0.3])
  y3 = Voigt1(x,[10,0.5,0.4,0.0])
  y=y0+y1+y2+y3
  y=y+(0.5*randomu(seed,nn))
  set = Make_Set(x,y)
ENDIF

npts = N_Elements(y)
IF N_Elements(x) EQ 0 THEN x=FindGen(npts)

;
; calculates derivative
;
yp = Deriv(x,y)   ; calculates y prime

;
; flags (100) the points with positive derivative
;
i_plus = Where(yp GE 0)

yp_plus_flag = yp*0
yp_plus_flag[i_plus] = 100

;
; flags the peaks (their derivative switches from >0 to <0)
; and compute the indices of all peaks.
;
ychange = yp*0
FOR i=0L,N_Elements(ychange)-2 DO BEGIN
 IF yp_plus_flag[i] GT 0 AND yp_plus_flag[i+1] EQ 0 THEN ychange[i]=1
ENDFOR
i_peaks = Where(ychange EQ 1) ; indices of peak points
ni_peaks = N_Elements(i_peaks)

IF Keyword_Set(optimize) THEN BEGIN
  FOR i=0L,ni_peaks-1 DO BEGIN
    ii = i_peaks[i]
    mmax = max( [y[(ii-1)>0],y[ii],y[(ii+1)<(npts-1)]], idx)
    IF idx NE 1 THEN i_peaks[i]=((i_peaks[i]-1+idx)>0)<(npts-1)
  ENDFOR
ENDIF

IF i_peaks[0] EQ -1 and 1-keyword_set(widget_off) THEN BEGIN
  Message,/Info,'No peaks found'
  IF SDep(/w) THEN itmp = Dialog_Message('PEAKFINDER: No peaks found.',$
	Dialog_Parent=group)
  pcutoff=0
  climits=0  
  npeaks=0
  error=0
  RETURN,0
ENDIF


;
; for each peak found compute:
; i) The number of points with derivarive positive sitting at its
; left plus the number of points with negative derivarive sitting at 
; its right. This is stored in ynpts
; ii) A weight for each peak, consisting of the addition of absolute value 
; of the derivative for the points calculated in i)
;
yweight = FltArr(ni_peaks)
ynpts = FltArr(ni_peaks)

FOR i=0L,ni_peaks-1 DO BEGIN
  ynpts[i]=1.0
  yweight[i]=0.0
  ; left search
  good =1
  iii=i_peaks[i]
  WHILE good EQ 1 DO BEGIN
    iii=(iii-1)>0
    IF yp[iii] GE 0 THEN BEGIN
       ynpts[i]=ynpts[i]+1 
       yweight[i]=yweight[i]+abs(yp[iii]) 
    ENDIF ELSE good=0
    IF iii EQ 0 THEN good=0
  ENDWHILE
  ; right search
  good =1
  iii=i_peaks[i]
  WHILE good EQ 1 DO BEGIN
    iii=(iii+1)<(ni_peaks-1)
    IF yp[iii] LT 0 THEN BEGIN
      ynpts[i]=ynpts[i]+1 
      yweight[i]=yweight[i]+abs(yp[iii]) 
    ENDIF ELSE good=0
    IF iii EQ ni_peaks-1 THEN good=0
  ENDWHILE
ENDFOR

; 
; Now calculate the number of peaks for nnn diffent cutoff values. 
; (the good peaks are those which have weight larger than cutoff value)
; 
nnn = 100
xtmp = FindGen(nnn)/Float(nnn-1)
ytmp = FltArr(nnn)
FOR i=0,nnn-1 DO BEGIN
  cutoff = xtmp[i]*max(yweight)
  itmp = Where(yweight GE cutoff)
  IF itmp[0] EQ -1 THEN nn=0 ELSE nn=N_Elements(itmp)
  ;i_change = i_change[itmp]
  ;yweight = yweight[itmp]
  ytmp[i]=nn
ENDFOR
pCutoff = Make_Set(xtmp,ytmp) ; optional keyword (returned)
;
; Now calculates the histogram of the number pf peaks as a function 
; of the cutoff value
; 
htmp = Histogram(ytmp)
index=0
hmax=max(htmp,index)
nPeaks = index+1
if not silent then print,'PEAKFINDER: Most probable number of maxima and weight: ',nPeaks,hmax
cutoff = Where(ytmp EQ (index+1)) 
cutoff_1 = xtmp(cutoff[0])
cutoff_2 = xtmp(cutoff[N_Elements(cutoff)-1])
cutoff_a = 0.5*(cutoff_1+cutoff_2)
cutoff_w = cutoff_2-cutoff_1
if not silent then print,'PEAKFINDER: Stabilized cutoffs (min, max, average): ',$
	cutoff_1,cutoff_2,cutoff_a
if not silent then print,'PEAKFINDER: Cutoff width: ',cutoff_w
cLimits=[cutoff_a,cutoff_w,cutoff_1,cutoff_2]

xx=x[i_peaks]
yy=y[i_peaks]
cc=yweight/max(yweight)
IF Keyword_Set(sort) THEN BEGIN
 tmp = Reverse(Sort(yweight))
 i_peaks = i_peaks[tmp]
 yweight = yweight[tmp]
 ynpts = ynpts[tmp]
 xx = xx[tmp]
 yy = yy[tmp]
 cc = cc[tmp]
ENDIF

RETURN,Make_Set(i_peaks,xx,yy,yweight,cc,ynpts)

END
