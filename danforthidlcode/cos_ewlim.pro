;; Determine an equivalent width limit given Nsig, wavelength, 
;; b-value, and S/N using Equations 4-5,7,9-10 of Keeney et al. (2012).
;;
;; INPUTS
;;   Nsig   - significance level of limit (number of sigma)
;;   lambda - observed wavelength of line (A)
;;   b      - b-value (km/s)
;;   snpix  - S/N per pixel
;;   snx    - S/N measured at optimal resolution
;;
;; OPTIONAL PARAMETERS
;;   disp   - dispersion in mA/pix (defaults to G130M value for 
;;            lambda <= 1425 A and G160M value otherwise)
;;   bin    - the number of native COS pixels that the data have been
;;            pre-binned by, used to convert SNPIX to SN1
;;   xopt   - optimal integration width, in (binned) pixels
;;
;; B. Keeney 02/16/12
;; B. Keeney 06/12/14

function cos_ewlim, Nsig, lambda, b, snpix=snpix, snx=snx, $
                    disp=disp, bin=bin, xopt=xopt

if ((n_elements(snpix) eq 0) and (n_elements(snx) eq 0)) then begin
  print, 'EWLIM: either SNPIX or SNX must be specified!!!'
  return, -1
endif

if (n_elements(disp) eq 0) then disp = (lambda le 1425) ? 9.97d : 12.23d
if (n_elements(bin)  eq 0) then bin  = 1

dx = 1d3*b*lambda / (2.998d5*disp)

xopt = 1.605*dx + 5.1d*dx^(-0.25d)
eta  = 0.15d + xopt^(0.37d)
fcx  = 0.743d - 0.185d*exp(-dx/11.6d)

sn1 = (bin eq 1) ? snpix : snpix / (0.15 + bin^(0.37d))

ewlim = (n_elements(snx) gt 0) ? (Nsig*disp/snx) * xopt/fcx         : $
                                 (Nsig*disp/sn1) * xopt/(eta*fcx)

xopt /= bin

return, ewlim

end



