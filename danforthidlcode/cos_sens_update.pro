function cos_sens_update,wave,flux,segment=segment,lv=lv,grating=grating,p=p,sensvector=sensvector
if n_elements(grating) ne 1 then grating='G130M'
if n_elements(segment) ne 1 then segment='FUVA'
if n_elements(lv) ne 1 then lv=0

; pass in wave, segment, and grating
; and get out the smov sensitivity update
;
; USE the lv option for data taken after Aug 13, 2009
;
;
; divide the returned sensvector value times your SMOV flux to get updated values
;
nw=n_elements(wave)
nf=n_elements(flux)

if nw ne nf then message,'ERROR: the wavelength and flux vectors must have the same size'

useg=strupcase(segment)
ugrat=strupcase(strtrim(grating,2))

case ugrat of 
	'G130M' : begin
			pa=[2736.3667,-6.5024309,0.0040286034,1.2969438e-06,-2.0995578e-09,5.5387341e-13]
			pb=[7385.5333,-13.271672,-0.0048998948,2.4903635e-05,-1.8383126e-08,4.2931604e-12]
			pa_lv=[913.33789,5.3526649,-0.020921567,2.5117444e-05,-1.2848511e-08,2.4239719e-12]
			pb_lv=[6566.5360,-12.140919,-0.0032876628,2.0892240e-05,-1.5698942e-08,3.6928369e-12]
			end
	'G160M': begin
			pa=[-3279.1462,5.9503510,-0.0025139167,-1.2436785e-06,1.1917261e-09,-2.4044604e-13]
			pb=[-5080.8995,9.1067530,-0.0015959567,-6.0123991e-06,4.3706627e-09,-8.9840028e-13]
			pa_lv=[-1378.5339,2.6580205,-0.0014300817,-1.8917433e-07,3.6880377e-10,-8.1514001e-14]
			pb_lv=[-7073.9153,12.539628,-0.0019500167,-8.5522165e-06,6.1290053e-09,-1.2525614e-12]
			end
	'G140L' : begin
				message,'Sorry G140L data are not yet supported.'
			end
			
	else: message,'Bummer you input a bad grating -> '+ugrat
endcase

p=double((useg eq 'FUVA' ? (lv ? pa_lv :pa) : (lv ? pb_lv : pb)))

sensvector=dblarr(nw)

for i=0,5 do sensvector+=p[i]*wave^i

return,double(flux)/sensvector

end

