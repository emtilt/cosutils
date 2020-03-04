;+
;; Version 2.1
;;
;; Version 2.1 of COADD_X1D adds new functionality and changes one of
;; the default behaviors from v2.0.
;;
;; (1) The DQ array of the X1D files is now explicitly parsed to throw
;;     out some regions (e.g., invalid wavelengths) and de-weight
;;     others (bad pixels, grid wires).  Previously, anything
;;     flagged with DQ_WGT=0 was thrown out.  As part of this change,
;;     regions near the edges of the detectors are automatically
;;     de-weighted, so "flanging" the edges is no longer necessary.
;;     Users can still "flange" the data using the /FLANGE keyword,
;;     but it is no longer the default behavior.
;;
;; (2) Adam Jensen has added the /aonly and /bonly keywords to coadd an
;;     individual segment of the detector while ignoring the other.
;;
;;
;; Version 2.0 of COADD_X1D differs substantially from the versions
;; previously released, primarily in its default behaviors and
;; corresponding keywords.  
;;
;; (1) We have found that care must be taken in choosing an
;;     interpolation method to minimize the non-Poissonian noise
;;     introduced by coaddition (see Keeney et al. 2012, PASP for
;;     details). Therefore, we now use the /nointerp option to LINTERP
;;     to interpolate individual exposures onto a common wavelength
;;     scale. This has a small effect on the coadded flux and error
;;     values, but a profound effect on the noise properties of the
;;     final arrays and, by extension, the limiting equivalent widths
;;     achievable in the coadded data. 
;;
;; (2) Flat fielding is no longer performed by default, so the /no_flat 
;;     option has been removed.  If the user wishes to have COADD_X1D 
;;     flat field their data, they should ensure that no flat fielding 
;;     has been applied as part of standard CalCOS processing and then 
;;     choose either the /gridflat (the old default) or /iterflat 
;;     (previously, /newflat) options.
;;
;; (3) COADD_X1D now tries to scale flux levels in individual
;;     exposures relative to one another automatically.  This behavior
;;     can be turned off using the new /no_scale option, in which case
;;     the scale factors all default to 1.0, or explicitly overwritten
;;     using the familiar scale= keyword.
;;
;; (4) The dispersion of the output spectrum can now be specified
;;     explicitly using the dl= keyword, but defaults to the nominal
;;     values for each grating from the COS Instrument Handbook.
;-

pro COADD_X1D, wave, flux, err, $
               files=files, path=path, chan=chan, method=method, $
               no_align=no_align, indivshift=indivshift, ubershift=ubershift, $
               no_scale=no_scale, scale=scale, $
               gridflat=gridflat, iterflat=iterflat, $
               flange=flange, update_flux=update_flux, $
               aonly=aonly, bonly=bonly, dl=dl, plot=plot, $
               xrange=xrange, yrange=yrange, savefile=savefile

;+
;; COADD_X1D - generalized coaddition routine for aligning and
;;             coadding COS FUV observations 
;;
;; INPUT:  files = list of x1d.fits files to process (e.g. 'lb4r10z4q')   
;;                 If file list is not provided, code will search for
;;                 *x1d.fits in the path directory.
;;
;;         path = path to x1d.fits files.  If not provided, current
;;                working directory is assumed.
;;
;; OUTPUT:  wave,flux,err = coadded vectors (optional)
;; 
;; PROCESSING OPTIONS:
;;          chan = 1 G130M (default)
;;               = 2 G160M
;;               = 3 both M gratings
;;               = 4 G140L
;;          
;;          method = -1 simple mean (previously /lowsn)
;;                 =  0 
;;                 =  1 modified exposure weighting (default)
;;                 =  2 err^-2 weighting
;;                 =  3 (S/N)^2 weighting
;;
;;          /no_align - do not cross-correlate exposures
;;
;;          indivshift = exposure/side specific shift (in angstroms).
;;                       Vector must have dimensions of [nfiles,2]
;;
;;          ubershift = constant shift (in angstroms) to shift the
;;                      wavelength vector by
;;
;;          /no_scale - do not scale flux levels of individual
;;                      exposures relative to one another
;;
;;          scale = set of scaling factors to multiply input flux,
;;                  errors by.  Should be one scale per input file.
;;
;;          /gridflat - use the original flat field files for
;;                      gridwire removal (previous default behavior)
;;
;;          /iterflat - use the iterative flatfield files (Fall, 2010)
;;                      (previously called /newflat)
;;
;;          /flange - de-weight pixels at the edges of detectors in 
;;                    each exposure. As of v2.1 this should no longer
;;                    be necessary since de-weighting of these regions
;;                    is performed when parsing the DQ array
;;
;;          /update_flux - correct early flux calibration; generally
;;                         unnecessary with data downloaded >2010
;;
;;          /aonly - process only segment A
;;
;;          /bonly - process only segment B
;;
;;
;; OUTPUT OPTIONS:
;;          dl = dispersion in output wavelength array (mA/pix). 
;;               defaults to [9.97, 12.23, 80.3] mA/pix for G130M,
;;               G160M, and G140L, respectively.
;;
;;          plot = 1 basic flux/wavelength plot
;;               = 2 two panels: flux/wavelength and exptime/wavelength
;;
;;          [x,y]range = wavelength/flux axis ranges
;;      
;;          savefile = name for IDL save file, if keyword is set and
;;                     name not supplied, 'targ_grating_coadd.idl' will
;;                     be used. 
;;
;;
;; REQUIRES: 
;;          cos_sens_update.pro - Steve Penton's routine to correct
;;                                early CalCOS output. Only needed
;;                                when /update_flux keyword is set.
;;          wireflat_a.idl, wireflat_b.idl - pseudo-flatfields used
;;                                           in wire removal with
;;                                           /gridflat command
;;          FUV_FF_ITER.dat - iterative flatfield triggered with
;;                            /iterflat command
;;
;; NOT YET IMPLEMENTED: 
;;          rigorous wavelength solutions (particularly for G140L data)
;;          NUV support (probably won't ever happen)
;;
;; KNOWN BUGS
;;   Requires the version of LINTERP.pro distributed in the ASTROLIB
;;   package.  The common version in GHRSLIB differs slightly, but
;;   importantly.  An error of "Keyword "missing" not allowed in
;;   linterp" means you're default is probably the GHRS
;;   version. (12/14/10) 
;; 
;; HISTORY:  written by Charles Danforth (danforth@casa.colorado.edu)
;; substantial input and testing by Brian Keeney, Kevin France,
;; Yangsen Yao, Hao Yang, Steve Penton (U. Colorado) and Anand
;; Narayanan (Wisconsin) 
;;
;; 11/16/09 - beta release 1
;; 11/20/09 - grid wire removal, flanging, release 2
;; 11/30/09 - removed many external code dependencies
;;            added weighting options and processing options
;; 12/03/09 - release 3
;; 12/04/09 - modified exposure time weighting added, adopted as
;;            default, release 4
;; 12/07/09 - added scaling option
;; 01/13/10 - linear (ubershift) option
;; 01/15/10 - fixed coalignment bug 
;; 01/19/10 - fixed reference exposure bug
;; 02/01/10 - corrected method=1 error coaddition method,
;;            added multi-ref-exposure warning
;; 11/22/10 - re-corrected method=1 error coaddition
;;          - /newflat option to use newest iterative flat field
;;            files (later changed to /iterflat)
;; 12/14/10 - documented linterp bug
;; 03/09/12 - added DL keyword (BAK)
;;          - switched to using the DQ_WGT information from the x1d
;;            files to define good regions, and made the wavein,
;;            fluxin, errin arrays double precision (BAK)
;; 04/05/12 - significant changes to default behavior and keyword
;;            structures, summarized above (BAK)
;; 06/18/12 - added prop_id and pi_name into save data structure to be
;;            compatible with proc_cos and plotting functions (CD)
;; 08/13/12 - parsed the DQ array ourselves instead of relying on the
;;            DQ_WGT array to determine which data to throw out (BAK)
;-

; Keywords for if only one segment (or "side") is to be processed
sideend   = (keyword_set(aonly)) ? 0 : 1
sidestart = (keyword_set(bonly)) ? 1 : 0

;; path, variable, and parameter definitions
wirefile='~/bin/my_idlpro/wireflat_'+['a','b']+'.idl'   ; change path to reflect local structure
angstrom='!sA!r!u!9 %!x!n'                         ; makes an Angstrom symbol the correct way
channame=['G130M','G160M','G130M + G160M','G140L'] ; name used on plots
savefilechanname='_'+['g130m','g160m','both','g140l']+'_' ; name used in .idl savefile generation
segname=['FUVA','FUVB']                            ; detector names
disp=[0.00997d,0.01223d,-1,0.0803d]                ; grating dispersion (A/pix) 
xcorfeaturename=[['CII 1334.5','AlII 1670.8','n/a','CII 1334.5'], $
                 ['SiII 1260.4','SiII 1526.7','n/a','n/a']]
;minxcorwave    =[[1330,1664,-1,1300],[1255,1520,-1,300]]      
;maxxcorwave    =[[1340,1676,-1,1370],[1266,1533,-1,900]]
minxcorwave    =[[1330,1664,-1,1300],[1255,1520,-1,1000]]      
maxxcorwave    =[[1340,1676,-1,1370],[1266,1533,-1,1100]]

flange_size=100.                                   ; error-flanging scale height, pixels
snthresh=5.                                        ; local S/N threshhold for dynamic weighting
speedoflight=2.9979e5

;; parse weighting method
if not keyword_set(method) then method=1 ; modified exposure time weighting default
if method gt 3 or method lt -1 then method=1

;; if no path supplied, assume current working directory
if not keyword_set(path) then path=''

;; determine channel index to deal with various grating-specific
;; settings (dispersion, etc)
if not keyword_set(chan) then begin
  chan=1 ; G130M default
  message,/info,' no data channel specified.  Assuming '+channame[chan-1]
endif 
if chan eq 3 then begin
  chanind=[0,1] ; G130M+G160M
endif else begin
  chanind=chan-1 ; single-grating mode
endelse

;; determine input files
if keyword_set(files) then nfiles=n_elements(files) else begin
;  files=findfile2(path+'*_x1d.fits',count=nfiles)
  files=findfile(path+'*_x1d.fits',count=nfiles)
  if files[0] eq '' then begin
    message,/info,'  !!! no files found !!!'
    stop
  endif 
endelse

;; parse file names to extract roots (strip off path and '_x1d.fits',
;; if any)
for i=0,nfiles-1 do begin
  if strpos(files[i],'/') ne -1 then begin
    filetext=strsplit(files[i],'/',/extract) ; strip path out of filename
    files[i]=filetext[n_elements(filetext)-1]
  endif 
  posn=strsplit(files[i],'_') ; locate the positions of '_' in the filename
  files[i]=strmid(files[i], 0, posn[n_elements(posn)-1]-1) ; keep everything 
                                                           ; before the last '_'
endfor

;; error checking on scale keyword
if ((n_elements(scale) gt 0) and (n_elements(scale) ne nfiles)) then begin
  message,/info,'SCALE must have NFILES elements; using /NO_SCALE instead...'
  no_scale = 1
endif

;; treat all exposures equally if /no_scale is set
if keyword_set(no_scale) then scale=replicate(1d,nfiles) 

;; define exposure quantities
grating=strarr(nfiles) ; G130M, G160M, etc
cenwave=intarr(nfiles) ; centra wavelength (blue edge of A-side grating)
fppos=intarr(nfiles)   ; FP position (1-4)
filenote=strarr(nfiles)
pi_name=strarr(nfiles) ; PI name
prop_id=intarr(nfiles) ; HST proposal ID

;; loop through files, reading header information
for i=0,nfiles-1 do begin
;; grating and central wavelength are in fits extension 0
  junk=mrdfits(path+files[i]+'_x1d.fits',0,hdr0,/silent)
  grating[i]=strtrim(sxpar(hdr0,'OPT_ELEM'),2)
  cenwave[i]=fix(sxpar(hdr0,'CENWAVE'))
  fppos[i]=fix(sxpar(hdr0,'FPPOS'))
  pi_name[i]=sxpar(hdr0,'PR_INV_L')
  prop_id[i]=sxpar(hdr0,'PROPOSID')
endfor

;; pick only files matching the selected grating(s)
if n_elements(chanind) eq 1 $
  then g=where(grating eq channame[chanind]) $
  else g=where(grating eq channame[chanind[0]] or grating eq channame[chanind[1]])
if g ne [-1] then begin
  grating=grating[g]
  cenwave=cenwave[g]
  fppos=fppos[g]
  files=files[g]
  nfiles=n_elements(g)
endif else begin
  message,/info,'!!! WARNING !!!'
  message,/info,' no files match requested grating: '+channame[chan-1]
  message,/info,' exiting!!!'
  return
endelse

;; Determine some meta info
spec=mrdfits(path+files[0]+'_x1d.fits',1,fitshdr,/silent)
targname=strtrim(sxpar(hdr0,'TARGNAME'),2)
nwave = n_elements(spec[0].wavelength)

;; loop through files and load spectral data into structures
wavein=dblarr(nwave,nfiles,2) ; wavelength
fluxin=dblarr(nwave,nfiles,2) ; flux
errin=dblarr(nwave,nfiles,2) ; error
pix2waveshift=dblarr(nfiles,2) ; pixel to wavelength shift (pixels, used in flatfielding)
exptin=dblarr(nwave,nfiles,2) ; exposure time (seconds)
dateobs=strarr(nfiles) ; Observation date (UT)
timeobs=strarr(nfiles) ; observation time (UT)

;if keyword_set(newflat) then restore,'FUV_FF_ITER.dat' ; fall 2010 flatfield files
if keyword_set(iterflat) then restore,'~/bin/my_idlpro/wireflats_iter.idl' ; processed 2010 flats

for i=0,nfiles-1 do begin
;; spectral data as well as exposure time is stored in fits extension 1
  spec=mrdfits(path+files[i]+'_x1d.fits',1,hdr1,/silent)
  dateobs[i]=strtrim(sxpar(hdr1,'DATE-OBS'),2)
  timeobs[i]=strtrim(sxpar(hdr1,'TIME-OBS'),2)
  pix2waveshift[i,*]=-1.*[sxpar(hdr1,'SHIFT1A'),sxpar(hdr1,'SHIFT1B')]

;; parse segments A, B into wave, flux, err vectors
  for j=sidestart,sideend do begin ; FUVA, FUVB
    waveintmp = spec[j].wavelength


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; UPDATE FLUX CALIBRATION (cos_sens_update can't (yet) handle G140L data)
;; Data calibrated under older versions of CalCOS contained
;; dramatically-incorrect flux calibration as opposed to the only
;; slightly-incorrect flux calibration now.  This step shouldn't be
;; necessary for more current reductions. 
    if keyword_set(update_flux) and chanind[0] ne 3 then begin
      fluxintmp=cos_sens_update(spec[j].wavelength,spec[j].flux,seg=segname[j],grat=strtrim(grating[i],2),/lv)
      errintmp=cos_sens_update(spec[j].wavelength,spec[j].error,seg=segname[j],grat=strtrim(grating[i],2),/lv)
    endif else begin
      fluxintmp=spec[j].flux
      errintmp=spec[j].error
    endelse 
    exptintmp=replicate(spec[j].exptime,n_elements(fluxintmp))


;; limit errors - in error-weighted coaddition, undue weight is given
;;                to points with ridiculously small errors.  For data
;;                where there are significant regions of F~0, use the
;;                methods=+-1 keyword.
;minerr = median(errin)-2.*stddev(errin)>1d-16
    minerr=1d-17
    errin  = errin>minerr       ; abs(errin) > 1d-17


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; CORRECT FOR GRID WIRES (wire flat)
;;
;; Kevin France has determined the best-available flat fields for the
;; G130M detectors.  This includes both sensitivity changes and grid
;; wire locations in pixel space.  For the moment, I don't bother with
;; the <10% flat-fielding issues and concentrate on the major 
;; instrumental features which are narrow and >10% opaque (in
;; particular, shadows from grid wires).  These profiles are
;; detector-specific and thus should work the same for
;; G130M/G160M/G140L observations made on the same detector (FUVA/B).
;; The corrections are made in detector coordinates, not 
;; wavelength, so the 1-D flat is shifted by an amount listed in the
;; header keywords 'SHIFT1A' and 'SHIFT1B'.  Flux vectors are divided
;; by the 'wire flat' while errors and local exposure time vectors are
;; divided by a fattened version of the wire flat, squared, where each
;; wire region is broadened by ~2 resel.  A 15% opacity in a wire
;; should yield a 38% increase in error which in turn will give it
;; only ~62% of the weight in exposure time weighted coaddition.  
    if (keyword_set(gridflat) or keyword_set(iterflat)) then begin
      if keyword_set(iterflat) then begin
;      	message,/info,'USING ITERATIVE FLATS!'
      	case grating[i] of 
          'G130M' : begin
;            pseudoflat=(j eq 0  ? G130MA_FF : G130MB_FF) < 1.2
            pseudoflat=(j eq 0 ? G130MA_F0 : G130MB_F0)
          end
          'G160M' : begin
;            pseudoflat=(j eq 0  ? G160MA_FF : G160MB_FF) < 1.2
            pseudoflat=(j eq 0 ? G160MA_F0 : G160MB_F0)
          end
          else: message,/info,'G140L iterative flats not supported'
        endcase
      endif else begin
        restore,wirefile[j]
      endelse 

;; flatten wires out of flux vector
      flat=shift(pseudoflat,round(pix2waveshift[i,j]))
;      flat=shift_interp(pseudoflat,pix2waveshift[i,j]) ; non-integral pixel shift
      fluxintmp=fluxintmp/flat
;; broaden wire marks to properly flag error vector with wire locations
      fatflat=flat*0
      errbox=7 ; pixels, half-width
      for k=errbox,n_elements(flat)-errbox-1 do $
            fatflat[k]=min(flat[k-errbox:k+errbox])
      errintmp=errintmp/fatflat^2
      exptintmp=exptintmp*fatflat^2
    endif
    

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; PARSE THE DQ ARRAY TO DETERMINE WHICH (OSTENSIBLY BAD) PIXELS TO
;; THROW OUT ALTOGETHER AND WHICH TO DE-WEIGHT
;;
;;     1 = 2^0  = DQ Softerr        --> ignore
;;     2 = 2^1  = DQ Brush Mark     --> ignore
;;     4 = 2^2  = DQ Grid Shadow    --> de-weight
;;     8 = 2^3  = DQ Near Edge      --> de-weight
;;    16 = 2^4  = DQ Dead           --> de-weight
;;    32 = 2^5  = DQ Hot            --> de-weight
;;    64 = 2^6  = DQ Burst          --> throw out???
;;   128 = 2^7  = DQ Out of Bounds  --> throw out
;;   256 = 2^8  = DQ Data Fill      --> throw out???
;;   512 = 2^9  = DQ PH Low         --> throw out
;;  1024 = 2^10 = DQ PH High        --> throw out
;;  2048 = 2^11 = DQ Bad Time       --> throw out
;;  4096 = 2^12 = DQ Bad Wavelength --> ignore
;;  8192 = 2^13 = DQ Divots         --> de-weight
;; 16384 = 2^14 = DQ Sdistortion    --> de-weight
;;
;; In principle, we would throw out DQ Bad Wavelength, which is meant
;; to flag regions where the wavelength is below 900 A (COS Data
;; Handbook), but this flag is set for some G130M/G160M data points so
;; it seems untrustworthy.  Perhaps it's set properly for G140L data?

    ; Throw out appropriate values
    aa = where((spec[j].dq and 2^6)  or (spec[j].dq and 2^7)  or $
               (spec[j].dq and 2^8)  or (spec[j].dq and 2^9)  or $
               (spec[j].dq and 2^10) or (spec[j].dq and 2^11),   $
               naa, comp=bb, ncomp=nbb)
    if (naa gt 0) then begin
      fluxintmp[aa] = 0d
      errintmp[aa]  = 0d
      exptintmp[aa] = 0d
    endif


    ; De-weight other values as appropriate
    if (nbb gt 0) then begin
      cc = where((spec[j].dq[bb] and 2^2)  or (spec[j].dq[bb] and 2^3)  or $
                 (spec[j].dq[bb] and 2^4)  or (spec[j].dq[bb] and 2^5)  or $
                 (spec[j].dq[bb] and 2^13) or (spec[j].dq[bb] and 2^14), ncc)
      if (ncc gt 0) then begin
        if (method eq 1) then exptintmp[bb[cc]] /= 2d else $
                              errintmp[bb[cc]]  *= 2d
      endif
    endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 'FLANGE' the outer few angstroms of error and exposure time
;; vectors.   This is a kludge to give contributions to the sum from
;; the edges of detectors where flux calibration is less certain less
;; weight than those from the middle.  Flange scale height in pixels
;; is set by flange_size= above. 
    if keyword_set(flange) then begin
      edge=where((spec[j].dq and 2^7),comp=active,ncomp=npix)
      edgedist=findgen(npix)
      flange=1.+((exp(-1.*(edgedist+1)/flange_size)) > $
                 (exp(-1.*(npix-2-edgedist)/flange_size)))
      errintmp[active]*=flange
      exptintmp[active]/=flange^2
    endif 

;; copy temporary vectors into wavein/fluxin/errin arrays
    wavein[*,i,j]=waveintmp
    fluxin[*,i,j]=fluxintmp
    errin[*,i,j]=errintmp
    exptin[*,i,j]=exptintmp

;; filter out single-pixel down-spikes in error by comparing values to
;; seven-pixel medians.  This takes care of some of it, but there are
;; still problems. 2/2/10
    smootherrin=median(errin[*,i,j],7)
    for k=3,nwave-4 do if errin[k,i,j] le .1*smootherrin[k] then errin[k,i,j]=median(errin[k-3:k+3,i,j])

  endfor                          ; side A/B loop
endfor                          ; exposure loop


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; shift G140L side B data (added 3/19/10)
;;
;; The wavelength solution for the G140L side B data is grossly
;; wrong.  Since no good wavelength solution exists, I assume that the
;; gap between A and B detectors (1589+-25 pixels) is representative
;; and shift the B-side data such that it is this many pixels to the
;; blue of the blue edge of the side A data for each exposure. 

;g140lgap=123. ; angstroms

;for i=0,nfiles-1 do begin
;  if grating[i] eq 'G140L' then begin
;    rawgap=min(wavein[*,i,0])-max(wavein[*,i,1])
;    wavein[*,i,1]=wavein[*,i,1]-g140lgap+rawgap
;    print,'*** Applying G140L-B gap correction of '+string(rawgap-g140lgap,'(F6.1)')+' Angstroms'
;  endif 
;endfor


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; shift diagnostics
;; this was a section of code used when I was figuring out the sense
;; of the shifts required in the X-correlation.  It's left in the code
;; for completeness.  Feel free to ignore it.
;;
;; replace actual data with fake data
;fakeshift=20.*randomu(systime(/seconds),nfiles)
;vel1=([1260.4,1335.5,1526,1670.8]-1215.67)/1215.67*speedoflight
;for i=0,nfiles-1 do begin
;  for j=0,1 do begin
;    fluxin[*,i,j]=1.
;    line=where((wavein[*,i,j] ge 1335.4 and wavein[*,i,j] le 1335.6) or (wavein[*,i,j] ge 1260.3 and wavein[*,i,j] le 1260.5) or (wavein[*,i,j] ge 1525.9 and wavein[*,i,j] le 1526.1) or (wavein[*,i,j] ge 1670.7 and wavein[*,i,j] le 1670.9))
;    fluxin[line,i,j]=.1
;    wavein[*,i,j]=wavein[*,i,j]+i*disp[0]
;  endfor
;endfor
;
;print,fakeshift
;
;stop


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; COALIGN DATA
;; For the moment, the code assumes the CalCOS wavelength scale is
;; relatively correct.  As of Nov '09, the best wavelength 
;; solutions were determined from FPPOS=3 at the middle grating
;; position in each band (1309 and 1600).  Coalign on these positions
;; if available, and on FPPOS=3 if not.  If neither is available, the
;; first exposure is used as the reference.  The final coadded data
;; should be relatively aligned, but may need a global shift into
;; whatever rest frame you like.  This can be done using the
;; ubershift= keyword or in post-processing.
;;
;; Data is cross-correlated based on the spectral region around a
;; single strong feature in each segment.  The given data is based on
;; strong ISM features which appear in most ISM and IGM continuum
;; targets.  If this doesn't work for you, substitute your own regions
;; of interest (emission lines, etc.).  The vector xcorfeaturename is
;; purely for your own edification and not used by the code in any way.
;; 
;; If you wish to add a specific shift for each exposure/side, use the
;; keyword indivshift=[] (added 5/27/10)

xshift=dblarr(nfiles,2) ; wavelength shift for each exposure

if keyword_set(indivshift) then no_align=1.

if nfiles gt 1 and not keyword_set(no_align) then begin
  bestgrate=[1309,1600,-1,1230]
;; pick reference exposure(s)
  for i=0,n_elements(chanind)-1 do begin ; do this twice if G130M+G160M, once otherwise
    ref=where(grating eq channame[chanind[i]] and cenwave eq bestgrate[chanind[i]] and fppos eq 3) ; center grating setting, FPPos=3 (best option)
    if ref eq [-1] then ref=where(grating eq channame[chanind[i]] and fppos eq 3) ; anything at FPpos=3 (2nd choice)
    if ref eq [-1] then ref=where(grating eq channame[chanind[i]]) ; use first file 
    if n_elements(ref) gt 1 then begin
      print,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      print,'!!! NOTE: multiple reference exposures available for grating '+grating[ref[0]]+':'
      for j=0,n_elements(ref)-1 do begin
        print,strtrim(j,2)+'- '+files[ref[j]]+'  '+string(grating[ref[j]],cenwave[ref[j]],fppos[ref[j]],max(exptin[*,ref[j],0]),median(fluxin[*,ref[j],0]),max(exptin[*,ref[j],1]),median(fluxin[*,ref[j],1]),'(A5,I5,2x,I1,x,F8.1,E9.2,F8.1,E9.2)')+'  '+dateobs[ref[j]]+' '+timeobs[ref[j]]
      endfor 
      print,'  Using exposure 0 above as alignment and scaling reference.'
      print,'  Re-run code with explicit file list if this is non-optimal.'
      print,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    endif 
    ref=ref[0]
    filenote[ref]=filenote[ref]+'*';alignment reference exposure '

;; identify other exposures from the same grating
    thisgrat=where(grating eq grating[ref])
    
;; x-correlate exposures to that grating, by side
;    if grating[ref] eq 'G140L' then xcor_width=300 else 
    xcor_width=30               ; 30 pixels = 3A (G140L) or 0.3A (M gratings)
    for j=sidestart,sideend do begin ; side A, B
      for k=0,n_elements(thisgrat)-1 do begin
        
; define equivalent wavelength ranges to compare
        refrange=where(wavein[*,ref,j] ge minxcorwave[chanind[i],j] $
                       and wavein[*,ref,j] le maxxcorwave[chanind[i],j])
        comprange=where(wavein[*,thisgrat[k],j] ge minxcorwave[chanind[i],j] $
                        and wavein[*,thisgrat[k],j] le maxxcorwave[chanind[i],j])
; cross-correlate comparison with reference data
        refx=wavein[refrange,ref,j] 
        refy=fluxin[refrange,ref,j]
        compx=wavein[comprange,thisgrat[k],j] 
        compy=fluxin[comprange,thisgrat[k],j]
 
; invoke Herczeg Clause (non-standard central wavelength, FP position
; check)
      herczeg:
        if abs(n_elements(compx)-n_elements(refx)) gt 2 then begin
          if n_elements(compx) lt n_elements(refx) then begin
            print,'!!! comparison data '+files[thisgrat[k]]+' does not cover entire reference cross-correlation range'
            print,'!!!    reference range:'+string(minxcorwave[chanind[i],j],maxxcorwave[chanind[i],j],'(2F10.2)')
            print,'!!!    comparison data:'+string(minmax(compx),'(2F10.2)')
            filenote[thisgrat[k]]=filenote[thisgrat[k]]+'cross-correlation range truncated! '
          endif 
          if n_elements(refx) lt n_elements(compx) then begin
            print,'!!! reference data '+files[ref]+' does not cover entire reference cross-correlation range.'
            print,'!!!    reference range:'+string(minxcorwave[chanind[i],j],maxxcorwave[chanind[i],j],'(2F10.2)')
            print,'!!!    reference data:'+string(minmax(refx),'(2F10.2)')
            filenote[ref]=filenote[ref]+'reference range truncated! '
          endif 
          print,'!!! Truncating reference range to attempt a fix.'
          refxrange=minmax(refx) & compxrange=minmax(compx)
          userange=[refxrange[0]>compxrange[0],refxrange[1]<compxrange[1]]
          compgood=where(compx ge userange[0] and compx le userange[1])
          refgood =where(refx ge  userange[0] and refx  le userange[1])
          compx=compx[compgood] & compy=compy[compgood]
          refx=refx[refgood] & refy=refy[refgood]
          goto,herczeg
        endif 

; diagnostics
;        print,i,j,k
;        print,'reference= '+grating[ref]+string(cenwave[ref],fppos[ref],minxcorwave[chanind[i],j],maxxcorwave[chanind[i],j],'(I5,I2,2F10.2)')
;        print,'comparison='+grating[thisgrat[k]]+string(cenwave[thisgrat[k]],fppos[thisgrat[k]],'(I5,I2)')

        CROSS_CORRELATE,refy,compy,shift,corr,width=xcor_width
        
        xshift[thisgrat[k],j]=shift*disp[chanind[i]]
        wavein[*,thisgrat[k],j]=wavein[*,thisgrat[k],j]+xshift[thisgrat[k],j] 
        
      endfor ; exposure loop
    endfor ; side A/B loop
  endfor ; grating loop
ENDIF else begin
  if not keyword_set(indivshift) then begin
    xshift=replicate(0.,nfiles,2) ; /no_align keyword set
  endif else begin 
    xshift=indivshift
    for i=0,nfiles-1 do begin ; exposure
      for j=sidestart,sideend do begin ; side
        print,'Individual wavelength shift for ',i,j,xshift[i,j],indivshift[i,j]
        wavein[*,i,j]=wavein[*,i,j]+xshift[i,j]
      endfor 
    endfor 
  endelse
endelse 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; SCALE FLUX, ERROR - Scale factors determined automatically by
;;                     default but can be overridden by the scale= and
;;                     /no_scale keywords

if (n_elements(scale) eq 0) then begin
  scale = dblarr(nfiles)

  ;; pick reference exposure(s)
  bestgrate=[1309,1600,-1,1230]

  ; do this twice if G130M+G160M, once otherwise
  if (chan eq 3) then begin
    ref1 = 0
    ref2 = 0
  endif

  for i=0,n_elements(chanind)-1 do begin 
    ; center grating setting, FPPos=3 (best option)
    ref=where((grating eq channame[chanind[i]]) and $
              (cenwave eq bestgrate[chanind[i]]) and (fppos eq 3))

    ; anything at FP-POS = 3 (2nd choice)
    if ref eq [-1] then ref=where((grating eq channame[chanind[i]]) and $
                                  (fppos eq 3))

    ; use first file (3rd choice)
    if ref eq [-1] then ref=where(grating eq channame[chanind[i]])

    ref=ref[0]


    ; Save reference exposure info for use later
    if (chan eq 3) then $
      if (i eq 0) then ref1 = ref else ref2 = ref


    ; identify other exposures from the same grating
    thisgrat=where(grating eq grating[ref])
    
    for k=0,n_elements(thisgrat)-1 do begin
      ; Determine wavelength regions for comparison
      if (grating[thisgrat[k]] eq 'G130M') then begin
        aa = WHERE((wavein[*,thisgrat[k],0] ge 1378.0) and $
                   (wavein[*,thisgrat[k],0] le 1408.0)) 
        bb = WHERE((wavein[*,thisgrat[k],1] ge 1240.0) and $
                   (wavein[*,thisgrat[k],1] le 1270.0))

        dla = wavein[1,thisgrat[k],0] - wavein[0,thisgrat[k],0]
        dlb = wavein[1,thisgrat[k],1] - wavein[0,thisgrat[k],1]

        scale[thisgrat[k]] = total(fluxin[aa,thisgrat[k],0])*dla + $
                             total(fluxin[bb,thisgrat[k],1])*dlb
      endif else if (grating[thisgrat[k]] eq 'G160M') then begin
        aa = WHERE((wavein[*,thisgrat[k],1] ge 1510.0) and $
                   (wavein[*,thisgrat[k],1] le 1538.0), naa) 
        bb = WHERE((wavein[*,thisgrat[k],1] ge 1452.0) and $
                   (wavein[*,thisgrat[k],1] le 1484.0), nbb)

        dla = wavein[1,thisgrat[k],1] - wavein[0,thisgrat[k],1]
        dlb = wavein[1,thisgrat[k],1] - wavein[0,thisgrat[k],1]

        scale[thisgrat[k]] = total(fluxin[aa,thisgrat[k],1])*dla + $
                             total(fluxin[bb,thisgrat[k],1])*dlb
      endif else begin
        message, /info, 'Autoscaling not yet supported for G140L.'
      endelse
    endfor ; exposure loop

    ; Normalize so that reference exposure has scale=1
    scale[thisgrat] = scale[ref] / scale[thisgrat]
  endfor ; grating loop


  ; Snippet of code to normalize G160M relative to G130M, if necessary.
  ; Assumes that the reference exposures have some overlapping wavelengths.
  if (chan eq 3) then begin
    aa = where(wavein[*,ref1,0] ge min(wavein[*,ref2,1]), naa)
    bb = where(wavein[*,ref2,1] le max(wavein[*,ref1,0]), nbb)

    dla  = wavein[1,ref1,0] - wavein[0,ref1,0]
    dlb  = wavein[1,ref2,1] - wavein[0,ref2,1]
    scl  = total(fluxin[aa,ref1,0])*dla / (total(fluxin[bb,ref2,1])*dlb)
    scl *= (wavein[nbb-1,ref2,1]-wavein[0,ref2,1]) / $
           (wavein[naa-1,ref1,0]-wavein[0,ref1,0])

    lastgrat = where(grating eq grating[ref2])
    scale[lastgrat] *= scl
  endif
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; OBSERVATION INFO: list input files and parameters
print,'COADD_X1D: Exposure Summary'
print,'                       cen  FP   ------ side A -------   ------ side B -------' 
print,'Exposure name grating wave pos   t_exp   <F>     shift   t_exp   <F>     shift  Scale    Date      Time   '
print,'------------------------------------------------------------------------------------------------------------'
for i=0,nfiles-1 do print,string(i,'(I2)')+' - '+files[i]+'  '+ string(grating[i],cenwave[i],fppos[i],max(exptin[*,i,0]),median(fluxin[*,i,0]),xshift[i,0], max(exptin[*,i,1]),median(fluxin[*,i,1]),xshift[i,1],scale[i],'(A5,I5,2x,I1,x,F8.1,E9.2,F7.3,F8.1,E9.2,2F7.3)')+'  '+dateobs[i]+' '+timeobs[i]+' '+filenote[i]
print,' (* alignment+scaling reference exposure)'
print,' '


;; apply scale factors to input data...
for i=0,nfiles-1 do begin
  fluxin[*,i,*] *= scale[i]
  errin[*,i,*]  *= scale[i]
endfor


;; collapse sides A and B into the same vector (remove third dimension
;; from the arrays)
IF (keyword_set(aonly) OR keyword_set(bonly)) THEN BEGIN
   wavein=[wavein[*,*,sidestart]]
   fluxin=[fluxin[*,*,sidestart]]
   errin=[errin[*,*,sidestart]]
   exptin=[exptin[*,*,sidestart]]
   sidenames=['A','B']
   sidename=sidenames[sidestart]
   side=replicate(sidename,nfiles)
ENDIF ELSE BEGIN
   wavein=[[wavein[*,*,0]],[wavein[*,*,1]]]
   fluxin=[[fluxin[*,*,0]],[fluxin[*,*,1]]]
   errin=[[errin[*,*,0]],[errin[*,*,1]]]
   exptin=[[exptin[*,*,0]],[exptin[*,*,1]]]
   grating=[grating,grating]
   cenwave=[cenwave,cenwave]
   fppos=[fppos,fppos]
   scale=[scale,scale]
   side=[replicate('A',nfiles),replicate('B',nfiles)]
   nfiles=nfiles*2
ENDELSE


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Define array for coadded wavelength
;; At the moment, I assume a linear dispersion for the data.  This is
;; approximately correct for the medium-resolution gratings, but
;; definitely not correct for G140L.  The G130M+G160M data get a
;; dual-dispersion solution broken at 1450A. 
goodwave=where(fluxin ne 0. and wavein gt 0.) ; changed 4/7/10
if (n_elements(dl) gt 0) then begin
    dl *= 1d-3   ; input dispersion is in mA/pix
    minwave = round(min(wavein[goodwave]));/disp)) * disp[chanind]
    maxwave = round(max(wavein[goodwave]));/disp)) * disp[chanind]
    wave = minwave + dindgen((maxwave-minwave)/dl + 1.)*dl
endif else begin
  if chan eq 3 then begin ; create dual-dispersion G130M+G160M vector
    wave=disp[0]*dindgen(33099)+min(wavein[goodwave]) ; G130M section
    good=where(wave le 1460.) & wave=wave[good] ; truncate at 1460A
    wave=[wave,max(wave)+disp[1]+disp[1]*dindgen(27800)] ; G160M section
    good=where(wave ge min(wavein[goodwave]) and wave le max(wavein[goodwave]))
    wave=wave[good]
  endif else begin
    minwave = round(min(wavein[goodwave]))
    maxwave = round(max(wavein[goodwave]))
    wave = minwave + dindgen((maxwave-minwave)/disp[chanind] + 1.)*disp[chanind]
  endelse
endelse


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Interpolate the individual exposures onto the coadded wavelength
;; vector.  (Note, this code requires the ASTROLIB version of
;; LINTERP.pro.  The common GHRSLIB version produces subtlely wrong
;; results).
fluxint=dblarr(n_elements(wave),nfiles) ; interpolated flux
errint=dblarr(n_elements(wave),nfiles) ; interpolated error
exptint=dblarr(n_elements(wave),nfiles) ; interpolated exposure time array

for i=0,nfiles-1 do begin
  good=where(wavein[*,i] gt min(wave)) ; use only real wavelength values
  LINTERP, wavein[good,i], fluxin[good,i], wave, fluxtmp, missing=0d, /nointerp
  fluxint[*,i]=fluxtmp[*]
  LINTERP, wavein[good,i], errin[good,i],  wave, errtmp,  missing=0d, /nointerp
  errint[*,i]=errtmp[*]
  LINTERP, wavein[good,i], exptin[good,i], wave, exptmp,  missing=0d, /nointerp
  exptint[*,i]=exptmp[*]
endfor


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Perform the coaddition
;;
;; COADDITION: flux elements are coadded pixel by pixel according to
;; several methods set by the keyword method.  In our extensive
;; testing, modified exposure time weighting with flanging seems to
;; produce the best coadditions with the least spurious features.
;; Thus we have made method=1 the default.  
;;   method=-1 - simple mean of pixel values, gives too much weight to
;;               short exposures
;;   method= 1 - modified exposure weighting:  Exposure time is
;;               modified at each pixel location by flanging and wire
;;               shadows (if selected).
;;   method= 2 - err^-2 weighting, allows error array tuning, but
;;               tends to weight toward lower-flux pixels.  
;;   method= 3 - (S/N)^2 weighting, allows for error array tuning, but
;;               tends to weight toward higher-flux pixels.

;; define output vectors
flux = dblarr(n_elements(wave))
err  = dblarr(n_elements(wave))
exptime = dblarr(n_elements(wave))

for i=0L,n_elements(wave)-1 do begin ; loop through by pixel
  good=where(fluxint[i,*] ne 0.) ; which exposures have data at this location?
  if good ne [-1] then begin
    exptime[i]=total(exptint[i,good]) ; track total exposure time per pixel

    if n_elements(good) eq 1 then begin ; only one exposure at this wavelength
      flux[i]=fluxint[i,good]
      err[i]=errint[i,good]

    endif else begin
      localsn=median(fluxint[i,good]/errint[i,good]) ; local pixel S/N
      
;; simple mean
      if method eq -1 then begin ; simple mean
        flux[i]=total(fluxint[i,good])/double(n_elements(good))
        err[i]=stddev(fluxint[i,good])/sqrt(double(n_elements(good)))
      endif 
        
;; modified exposure-weighted mean
      if method eq 1 then begin
        flux[i]=total(fluxint[i,good]*exptint[i,good])/total(exptint[i,good])
;        err[i]=stddev(fluxint[i,good])/sqrt(double(n_elements(good)))
;        err[i]=sqrt(1d/total(1d/errint[i,good]^2))
	err[i]=sqrt(total((exptint[i,good]*errint[i,good])^2))/total(exptint[i,good]) ; SVP method, fall 2010

      endif  

;; error-weighted mean
      if method eq 2 then begin
        wts=(1d/errint[i,good]^2)
        flux[i]=total(fluxint[i,good]*wts)/total(wts)
        err[i]=sqrt(1d/total(wts))
      endif

;; (flux/err)^2-weighted mean (flux(lam)/err(lam) should be equivalent to sqrt(t))
      if method eq 3 then begin
        wts=(fluxint[i,good]/errint[i,good])^2
        flux[i]=total(fluxint[i,good]*wts)/total(wts)
        err[i]=sqrt(1d/total(wts))
      endif 

    endelse
  endif 

  ; Enforce minimum flux, error in final coaddition
  flux[i] = flux[i] > 0d
  err[i]  = err[i] > minerr
endfor ; coaddition loop

if keyword_set(ubershift) then wave+=ubershift 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; PLOTting
if keyword_set(plot) then begin
  if plot eq 2 then begin ; fancy plot with exposure time sub-graph underneath
    !p.multi=[0,1,2]
    !y.margin=[-10,2]
  endif 

if not keyword_set(xrange) then xrange=[min(wave)-10.,max(wave)+10]
if xrange[1]-xrange[0] gt 100 then nsum=7 else nsum=1
  
;; flux, error vs wavelength
  nonairglow=where(abs(wave-1215.67) gt 20. and abs(wave-1305) gt 20. and wave gt min(wave)+2. and wave lt max(wave)-2.)
  if not keyword_set(yrange) then yrange=minmax(smooth(flux[nonairglow],49))
  plot,wave,flux,xrange=xrange,/xsty,yrange=[(-0.1*yrange[1])>yrange[0],yrange[1]],/ysty,nsum=nsum,xtit='Wavelength ('+angstrom+')',ytit='Flux (erg cm!u-2!n s!u-1!n '+angstrom+'!u-1!n)',tit=targname+' (COS/'+channame[chan-1]+')'
  oplot,wave,err,psym=3,color=128,nsum=nsum
  oplot,!x.crange,[0.,0.],lines=1
  
;; exposure time vs wavelength
  if plot eq 2 then begin
    !y.margin=[4,14]
    plot,wave,exptime/1000.,xtit='Wavelength ('+angstrom+')',ytit='Exposure (ksec)',yr=[0,1.1*max(exptime)/1000.],/ysty,xr=!x.crange,/xsty,nsum=nsum
    !y.margin=[4,2]
    !p.multi=0
  endif 
endif 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; SAVE everything as IDL binary file
if keyword_set(savefile) then begin
  if strtrim(savefile,2) eq '1' then savefile=strlowcase(targname)+savefilechanname[chan-1]+'coadd.idl'
  readme='Generated by COADD_X1D.pro: '+systime()
  save,readme,targname,fitshdr,files,grating,side,cenwave,fppos,dateobs,prop_id,pi_name,exptin,wavein,fluxin,errin,wave,flux,err,exptime,scale,file=savefile
  message,/info,' data saved as '+savefile
endif 

end