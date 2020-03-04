;+
;; COS_LSF.pro - return normalized COS LSF at nearest tabulated wavelength value
;;               based on data at 
;; http://www.stsci.edu/hst/cos/performance/spectral_resolution/ 
;; 
;; inputs:
;;   lambda - central wavelength for LSF 
;; 
;; outputs: 
;;   xarray - output vector of wavelengths (in angstroms or pixels)
;;
;; options:
;;   /pixel - dispersion units in pixels, else in angstroms (default) 
;;   chan= - Specifies data channel (g130m,g160m,g140l,g225m). Defaults to selecting based on
;;            wavelength, where chan='g130m' if (lambda le 1450), 'g160m' if 
;;            (lambda gt 1450 and lambda lt 1800), else chan='g225m'.
;;   vercos= - Scalar or 1D integer array specifying which LSF(s) to use.  
;;              Acceptable values are 1 for COS ISR 2009-01, 2 for COS ISR 2011-01,
;;              3 for the middle LTP2 central wavelength for the selected channel,
;;              or the central wavelength setting of a lifetime position 2 LSF.
;;              Valid values are: 1, 2, 3, 1055, 1096, 1222, 1291, 1300, 1309, 1318, 1327, 1105, 
;;              1230, 1280, 1577, 1589, 1600, 1611, 1623. If VERCOS undefined, defaults
;;              to VERCOS=2. If invalid value specified, closest valid value will be chosen.
;;              If (n_elements(VERCOS) gt 1), the LSFs at those values will be combined
;;              as a weighted mean according to keyword WEIGHTS. Note that if you mix channels,
;;		the pixel scales, and probably the whole LSF, will be wrong, but no checking
;;              is currently done to prevent that.
;;   weights= - Scalar or 1D array of weights used when combining LSFs.
;;              Must have (n_elements(WEIGHTS) eq n_elements(VERCOS)). Equal weighting is default.
;;   /nopad   - By default, routine pads LSF with zeros to the maximum number of pixels in any LSF.
;;              If this keyword is set, all elements of the LSF that are zero at the end of the routine
;;              are removed.
;;
;; VERSION NOTE: The 5/30/14 version and later is incompatible with calls to the older versions
;;  if those calls specify the OLD or NEW keywords. If the calls implicitly use
;;  the default behavior for those keywords, then this version is probably compatible.
;;  The addition of forced normalization may break some uses.
;;
;; POSSIBLE FUTURE AREAS FOR IMPROVEMENT:
;;  -Prevent mixing channels.
;;  -Are the present pixel scales precise enough?
;;  -Should we interpolate LSFs to get the wavelength of interest?
;;  -Improve centering of the LSF on the central pixel?
;;
;; C. Danforth 10/12/09 --> 12/2/09
;; B. Keeney  03/26/11: added the NEW keyword [Seemingly this was changed to OLD keyword? -ET, 5/30/14]
;; E. Tilton  05/30/14: Removed OLD keyword. Expanded to include LTP2 LSFs and composite LSFs.
;;			Routine now forces normalization of the LSF to allow the composites.
;; E. Tilton  06/02/14: Added NOPAD keyword. Additional minor improvements.
;; E. Tilton  09/22/14: Fixed vercos indexing bug.
;; E. Tilton  12/29/2016: add ltp3 but its badly programmed due to rush
;-

function cos_lsf,lambda,xarray,chan=chan,pixel=pixel,vercos=vercos,weights=weights,nopad=nopad,ltp3=ltp3,fixltp31055=fixltp31055
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;+
	;USER CHANGEABLE SETTINGS
        ;
	;save path for LSFs
	folder='~/idl/ltp2/'        ;path to folder containing LSFs
	;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;-
	if n_elements(fixltp31055) eq 0 then fitxltp31055=0
	if not keyword_set(nopad) then nopad=0
        if n_elements(ltp3) eq 0 then ltp3=0
	;if no LSF specified, use the COS ISR 2011-01 LSF
	if n_elements(vercos) eq 0 then vercos=[2]
	
	;if weights unspecified, use equal weighting
	if n_elements(weights) eq 0 then weights=make_array(n_elements(vercos),value=1.0)
	
	;cheack for right number of weights
	if n_elements(weights) ne n_elements(vercos) then begin
		message,'WEIGHTS must have same number of elements as VERCOS.'
	endif

	;sort vercos
	sind=sort(vercos)
	vercos=vercos[sind]
	weights=weights[sind]

	;the 1055 lsf at longer wavelengths doesn't seem to match observations based on EMT's qualitative observations.
	;this is a hack to deal with it. it just uses the generic LTP2 LSF, which qualitatively seems to match better.
	;beware. this may be incorrect for most purposes.
	if  fixltp31055 then begin
		ind1055=where(vercos eq 1055,found1055)
		if found1055 gt 0 then begin
			vercos[ind1055]=2
		endif
	endif
	

	;; assume medium resolution grating and pick according to lambda
	if not keyword_set(chan) then begin
		if lambda lt 1800. then begin ; assume FUV
			if lambda gt 1450. then chan='g160m' else chan='g130m'
		endif else begin
			chan='g225m' 
			print,'WARNING: Wavelength is NUV, but no channel specified. Assuming G225M.'
			print,'All NUV LSFs from year 2009 are the same.'
			print,'Shape will be correct, but pixel scale may be wrong.'
		endelse
	endif 
	chan = strtrim(strlowcase(chan),2)
	
	;choose closest valid vercos values
	valid=[1,2,3,4,1055,1096,1222,1291,1300,1309,1318,1327,1105,1230,1280,1577,1589,1600,1611,1623]
	ltp3valid=[1222,1291,1300,1309,1318,1327,1105,1230,1280,1577,1589,1600,1611,1623]
	;valid=[valid,ltp3valid]
	for i=0,n_elements(vercos)-1 do begin
		abdiff=abs(valid-vercos[i])
		mindiff=min(abdiff,mindex)
		vercos[i]=valid[mindex] 
	endfor

	;if LTP2 desired but cenwave not known, select middle cenwave in channel
	threeind=where(vercos eq 3,threecount)
	if threecount gt 0 then begin
		if chan eq 'g130m' then begin
			vercos[threeind]=1291
		endif else if chan eq 'g140l' then begin
			vercos[threeind]=1230
		endif else if chan eq 'g160m' then begin
			vercos[threeind]=1600
		endif else begin
			print,'WARNING: Channel invalid. Setting vercos=1 in these locations.'
			vercos[threeind]=1
		endelse
	endif


	;create array that will hold the LSFs to be used
	maxpixels=331 ;331 is the largest number of pixels that the present LSFs contain
	lsfarr=dblarr(n_elements(vercos),maxpixels) 
	lsfpixfull=indgen(maxpixels)-165

	;; load LSFs and pad to maxpixels pixels if necessary
	for i=0,n_elements(vercos)-1 do begin
		if (vercos[i] eq 1 or vercos[i] eq 2) or $
		  (vercos[i] ge 3 and ((chan ne 'g130m') and (chan ne 'g160m') and (chan ne 'g140l'))) $
		  then begin
			if vercos[i] eq 2 and ((chan eq 'g130m') or (chan eq 'g160m')) then begin
				savefile=folder+'cos_lsf_new.idl'
			endif else if vercos[i] eq 1 then begin
				savefile=folder+'cos_lsf.idl'
			endif else begin
				print,'WARNING:'
				print,'Selected channel not available with selected LSF version.'
				print,'Using 2009 LSFs instead.'
				savefile=folder+'cos_lsf.idl'
			endelse
			restore, savefile
			chanind=where(chan eq lsfchan)
			chanind=chanind[0]
			;; pick nearest wavelength point (LSF varies slowly enough with lambda
			;; that this should be good enough)
			junk=min(abs(lsfwave[chanind,*]-lambda),lamind)
			tmplsf=lsf[chanind,lamind,*]
		endif else begin
		    if ltp3 and vercos[i] gt 1097 then begin
			;lsfwave=read_table(folder+'ltp3/'+string(vercos[i],format='(I4)'),nrows=1)
			;lsf=read_table(folder+'ltp3/'+string(vercos[i],format='(I4)'),head=1)
		        filename=folder+'ltp3/'+string(vercos[i],format='(I4)')+'.idl'
			restore,filename

		    endif else begin
		        filename=folder+string(vercos[i],format='(I4)')+'.idl'
			restore,filename

		    endelse
		    ;; pick nearest wavelength point (LSF varies slowly enough with lambda
		    ;; that this should be good enough)
		    junk=min(abs(lsfwave-lambda),lamind)
		    tmplsf=lsf[lamind,*]
		endelse
		if n_elements(tmplsf) lt maxpixels then begin
			halfnum=(maxpixels-n_elements(tmplsf))/2
			tmplsf=[intarr(halfnum),transpose(tmplsf),intarr(halfnum)]
		endif

		lsfarr[i,*]=tmplsf/total(tmplsf)
	endfor

	;grating pixel scales from 
        ;http://www.stsci.edu/hst/cos/documents/handbooks/current/ch05.COS_Spectroscopy02.html#427421
	;retreived 2014-06-02
	lsfchan=    ['g130m','g160m','g140l','g185m','g225m','g285m','g230l']
	lsfpixscale=[9.97000,12.2300,80.3000,37.0000,33.0000,40.0000,390.000] 
	
	chanind=where(chan eq lsfchan)
	chanind=chanind[0]	
	if keyword_set(pixel) then xarray=lsfpixfull else xarray=lambda+lsfpixfull*0.001*lsfpixscale[chanind]
  
	;create mean lsf according to weights
	lsf=dblarr(n_elements(lsfpixfull))
	for i=0,n_elements(lsfpixfull)-1 do begin
		lsf[i]=total(weights*lsfarr[*,i])/total(weights)
	endfor
	
	;remove zeros if desired
	if nopad then begin
		zeroind=where(lsf ne 0)
		lsf=lsf[zeroind]
		xarray=xarray[zeroind]
	endif

	;check normalization
	lsf=lsf/total(lsf)

	return,lsf
end
