; stokes vectors the Tinbergen way
function calcstokes_tin, flux
	Rq = sqrt((flux[0]/flux[4])/(flux[2]/flux[6]))
	Ru = sqrt((flux[1]/flux[5])/(flux[3]/flux[7]))
	q = (Rq - 1)/(Rq + 1)
	u = (Ru - 1)/(Ru + 1)
	stokes = [q,u]
	return, stokes
end
; stokes vectors the Terry Jones way
function calcpol_tjj, flux
	q = (flux[0]-flux[4])-(flux[2]-flux[6])/(flux[0]+flux[4]+flux[2]+flux[6])
	u = (flux[1]-flux[5])-(flux[3]-flux[7])/(flux[1]+flux[5]+flux[3]+flux[7])
	return, sqrt(q^2+u^2)
end
; calculate polarization the ESO/Claudia/Matt way
function calcpol_eso, flux
	q = 0.5*(flux[0]-flux[4])/(flux[0]+flux[4])	- $
		0.5*(flux[2]-flux[6])/(flux[2]+flux[6])
	u = 0.5*(flux[1]-flux[5])/(flux[1]+flux[5]) - $
		0.5*(flux[3]-flux[7])/(flux[3]+flux[7])
	return, sqrt(q^2+u^2)
end

; I want this to return EVERYTHING -- make it a struct!!
; INPUTS: 	Flux values of the 8 ord/ext beams, total intensity, err per PIXEL 
;			of the total intensity image, bin size: [x,y], elements to loop 
;			over: [x,y]
; NOTE: 	Flux of the ord/ext beams must have same binning as Tot Intensity img
; OUTPUTS: 	Struct containing Q, U, P, P err, P STN, Theta, Theta err
function polarization, flux, int, sig, bin,TWO_D=two_d,MULTIBIN=multibin
	; check to see if the 2D keyword is set 
	; create arrays of appropriate sizes
	sf = size(flux)
	if KEYWORD_SET(two_d) then p = fltarr(sf[1],sf[2]) else p = fltarr(sf[1])

	; check to see if the MULITBIN keyword is set
	; create an err array with the same # of rows as the bin sizes
	sb = size(bin)
	if KEYWORD_SET(multibin) then err = fltarr(sb[2])

	; if we have a 2D array ... 
	if KEYWORD_SET(two_d) then begin
		for i=0,sf[1]-1 do begin
			for j=0,sf[2]-1 do begin
				p[i,j] = calcpol_eso(flux[i,j,*])		; calculate stokes 
				;p = sqrt(s[*,*,0]^2+s[*,*,1]^2) 			; polarization
				;t = .5*atan(s[*,*,1]/s[*,*,0])*180/3.14		; theta in degrees	
			endfor
			; Determine BIN ERROR -- ASSUMING THE BINS WERE SUMMED
			if KEYWORD_SET(multibin) then begin
;				; error per bin - many bin sizes
				err[i] = sqrt(bin[0,i]*bin[1,i])*sig		
			endif else begin
				; error per bin - 1 bin size for all
				if  i eq 0 then err = sqrt(bin[0]*bin[1])*sig
			endelse
		endfor
	; if we have a 1D array ... 
	endif else begin
		for i=0,sf[1]-1 do begin
			p[i] = calcpol_eso(flux[i,*])		; calculate stokes 
			;t = .5*atan(s[*,1]/s[*,0])*180/3.14		; theta in degrees
;			; Determine BIN ERROR
			if KEYWORD_SET(multibin) then begin
				; error per bin - many bin sizes
				err[i] = sqrt(bin[0,i]*bin[1,i])*sig		
			endif else begin
				; error per bin - 1 bin size for all
				if i eq 0 then err = sqrt(bin[0]*bin[1])*sig
			endelse
		endfor	
	endelse
	stn = int/err			; tot int S/N
	pstn = sqrt(2)*stn*p	; polarization S/N
	pe = 1/(sqrt(2)*stn)	; polarization error
;	te = pe/(2*p)			; theta error
;stop
	if KEYWORD_SET(two_d) then pol=create_struct('p',p,'pe',pe,'pstn',pstn,'stn',stn) $
	else pol=create_struct('p',p,'pe',pe,'pstn', pstn,'stn',stn)
	return, pol
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro pol2d, date
	; read in ord/ext beams
	spawn,'ls stacked/'+date+'/[A,B]*_stack_*'+date+'_2d.fits > stacked_sci.list'
	readcol,'stacked_sci.list',format='a',sci_files

	; read in total intensity images
	spawn,'ls stacked/'+date+'/I*_stack_*'+date+'_2d.fits > stacked_int.list'
	readcol,'stacked_int.list',format='a',int_files	

	; read in a test file to determine wavelength range(s)
	temp=mrdfits(sci_files[0],0,h)
	ss = size(temp)
	disp=sxpar(h,'CD1_1')
	pix=sxpar(h,'CRVAL1')
	wlen = dindgen(ss[1])*disp+pix
	y = dindgen(ss[2])
	; x range of 2d spectrum to look at
	xrange=where(wlen gt 4930. and wlen lt 5057.5)
	lx = xrange[0]
	; y range of 2d spectrum to look at
	yrange=where(y gt 3 and y lt 74)
	ly = yrange[0]
	; portion of the spectrum to use as the background
	srange=where(wlen gt 5025. and wlen lt 5500.)	
	; bin sizes in the x and y directions for total image binning
	bx = 8.
	by = 10.	
	; number of total bins in x and y directions
	num_x = (n_elements(xrange))/bx		
	num_y = (n_elements(yrange))/by

	; declare variables for total image binning
	spec=fltarr(ss[1],ss[2],n_elements(sci_files))		; 8 2d full spectra
	spec1d=fltarr(ss[1],n_elements(sci_files))
	spec1d_f = spec1d
	spec_c=fltarr(num_x*bx,num_y*by,n_elements(sci_files))
	spec_b=fltarr(num_x,num_y,n_elements(sci_files))	; 8 binned cropped 2d spectra
	spec1d_b=fltarr(num_x,n_elements(sci_files))
	stokes = fltarr(num_x,num_y,2)						; q,u for each bin/spectra
	pol = fltarr(num_x,num_y)							; p for each bin/spectra
	pol1d = fltarr(num_x)

	tot=fltarr(ss[1],ss[2],n_elements(int_files))		; 5 2d full tot int spectra
	tot1d=fltarr(ss[1],n_elements(int_files))
	tot_c=fltarr(num_x*bx,num_y*by,n_elements(int_files))
	stn_pp = tot_c
	tot_b=fltarr(num_x,num_y,n_elements(int_files))		; 5 binned cropped 2d spectra
	tot1d_b=fltarr(num_x,n_elements(sci_files))
	stn_pb = tot_b										; stn for each bin in tot int spec
	sig = fltarr(n_elements(int_files))					; rms for each tot int spectrum	

	; make 2d binned versions of the total intensity images
	for i=0,n_elements(int_files)-1 do begin
		; read in the total intensity spectra
		tot[*,*,i] = mrdfits(int_files[i],0,h)
		st=size(tot)
;		name=strsplit(int_files[i],'/.',/extract)
		; make 1D extracted total intensity spectrum
;		for j=0,ss[1]-1 do tot1d[*,i]=total(tot[j,*,i])
		; select the portion of the entire spectrum to be rebinned
;		tot_c[*,*,i] = tot[xrange,yrange,i]
;		tot1d_c[*,i] = tot1d[xrange,i]

		; REBIN FOR 2D
;		tot_b[*,*,i] = rebin(tot[xrange,yrange,i],num_x,num_y)
		; REBIN FOR 1D
		for j=0,st[1]-1 do tot1d[j,i] = total(tot[j,yrange,i])
		for j=0,num_x-1 do tot1d_b[j,i] = total(tot1d[lx+j*bx:lx+(j+1)*bx-1,i])		

;		writefits,name[0]+'/crop/'+name[1]+'_2dcrop.fits',tot_c[*,*,i]
;		writefits,name[0]+'/rebin/'+name[1]+'_2drebin.fits',tot_b[*,*,i]

		; Pixel Error: RMS PER PIXEL of each science frame (RMS of the BKGD)
		sig[i] = sigma(tot[srange,yrange,i]) 
		; S/N PER PIXEL of each science frame
;		stn_pp[*,*,i] = tot_c[*,*,i]/sig[i]			;;; write to file??
stop
	endfor

	for i=0,n_elements(sci_files)-1 do begin
		; read in the science spectra
		spec[*,*,i] = mrdfits(sci_files[i],0,h)
;		name = strsplit(sci_files[i],'.',/extract)
		; CROP SPECTRA
;		spec_c[*,*,i] = spec[xrange,yrange,i]

		; FLATTEN the 1d spectrum
		for j=0,ss[1]-1 do spec1d[j,i] = total(spec[j,yrange,i])
		mask = where((wlen gt 4700 and wlen lt 4900) or $
					(wlen gt 5050 and wlen lt 5400))		
		x=dindgen(ss[1])
		coeff=poly_fit(x[mask],spec1d[mask,i],1)
		y=coeff[0]+coeff[1]*x	
		spec1d_f[*,i] = spec1d - y
		plot, wlen,spec1d[*,i],xr=[4800,5500]
		oplot,wlen,spec1d_f[*,i],color=240
;stop
		; REBIN FOR 1D
		;spec1d_b[*,i] = rebin(spec[xrange,*,i],num_x)
		for j=0,num_x-1 do spec1d_b[j,i] = total(spec1d_f[lx+j*bx:lx+(j+1)*bx-1,i])		

;		writefits,name[0]+'/'+date+'/crop/'+name[1]+'_2dcrop.fits',spec_c[*,*,i],h
		; REBIN FOR 2D
;		spec_b[*,*,i] = rebin(spec[xrange,yrange,i],num_x,num_y)
;		writefits,name[0]+'/'+date+'/rebin/'+name[1]+'_2drebin.fits',spec_b[*,*,i],h
	endfor

;	for i=0,num_x-1 do pol1d[i] = calcpol_eso(spec1d_b[i,*])
;	err = sqrt(bx*n_elements(yrange))*sig[0]
;	stn = tot1d_b[*,0]/err		; tot int S/N
;	pstn = sqrt(2)*stn*pol1d	; polarization S/N
;	pe = 1/(sqrt(2)*stn)		; polarization error
;stop

	; Calculate POLARIZATION and ERROR for 2D BINNED SPECTRA
;	pol = polarization(spec_b,tot_b[*,*,0],sig[0],[bx,by],/two_d)
	; Calculate POLARIZATION and ERROR for 1D BINNED SPECTRA
	pol1d = polarization(spec1d_b,tot1d_b[*,0],sig[0],[bx,n_elements(yrange)])
stop
;	openw,out,'poldata2d_'+date,/get_lun
;	for i=0,num_x-1 do begin
;		for j=0,num_y-1 do begin
;			printf,out,pol.p[i,j],pol.pe[i,j]
;		endfor
;	endfor
;	close,out
;	free_lun,out

	openw,out,'pol1d_stack7814_2_fitsub',/get_lun
	for i=0,num_x-1 do begin
		printf,out,pol1d.p[i],pol1d.pe[i]
	endfor
	close,out
	free_lun,out
stop
	writefits,'pol_'+strtrim(fix(bx),2)+'_'+date+'.fits',pol1d.p
	writefits,'perr_'+strtrim(fix(bx),2)+'_'+date+'.fits',pol1d.pe
	writefits,'pstn_'+strtrim(fix(bx),2)+'_'+date+'.fits',pol1d.pstn

	writefits,'pol_'+strtrim(fix(bx),2)+'_'+strtrim(fix(by),2)+'_'+date+'.fits',pol.p
	writefits,'perr_'+strtrim(fix(bx),2)+'_'+strtrim(fix(by),2)+'_'+date+'.fits',pol.pe
	writefits,'pstn_'+strtrim(fix(bx),2)+'_'+strtrim(fix(by),2)+'_'+date+'.fits',pol.pstn
stop

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; polarization calculations
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; region information 
c1 = [[93,50],[93,40],[103,45],[110,46],[83,23],[91,20],[99,20],[107,20]];,[105,53]		; centers
s1 = [[12,8],[12,12],[8,10],[6,8],[8,10],[8,16],[8,16],[8,16]]	;,[8,8]	; spatial extent (pixels)

ss1 = size(c1)
bbins = fltarr(ss1[2],n_elements(sci_files))
tbins = fltarr(ss1[2])
for i=0,ss1[2]-1 do begin
	for j=0,n_elements(sci_files)-1 do begin
		bbins[i,j] = total(spec_c[c1[0,i]-s1[0,i]/2:c1[0,i]+s1[0,i]/2-1, $
									c1[1,i]-s1[1,i]/2:c1[1,i]+s1[1,i]/2-1, j])
	endfor
	tbins[i] = total(tot_c[c1[0,i]-s1[0,i]/2:c1[0,i]+s1[0,i]/2-1, $
								c1[1,i]-s1[1,i]/2:c1[1,i]+s1[1,i]/2-1, 0])
endfor
bpol = polarization(bbins,tbins,sig[0],s1)
;stop

; create a polarization "image" from scratch
bpol_img = spec_c[*,*,0]*0.
for i=0,ss1[2]-1 do begin
	if bpol.pstn[i] gt 1. then begin
		bpol_img[c1[0,i]-s1[0,i]/2:c1[0,i]+s1[0,i]/2-1, $
					c1[1,i]-s1[1,i]/2:c1[1,i]+s1[1,i]/2-1] = bpol.p[i]
	endif
endfor
writefits,'bpol.fits',bpol_img
;stop

; print region files
openw,out,'pol_spec.reg',/get_lun
printf,out, 'global color=red dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
printf,out,'image'
for i=0,ss1[2]-1 do begin
	printf,out,'box('+strtrim(c1[0,i],2)+','+strtrim(c1[1,i],2)+',' $
					 +strtrim(s1[0,i],2)+','+strtrim(s1[0,i],2)+',0)'
	printf,out,'# text('+strtrim(c1[i],2)+','+strtrim(c1[i],2)+') text={'+strtrim(i+1,2)+'}'
endfor
close,out
free_lun,out

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro fig_boxes

int=mrdfits('stacked/I_all_2dstack2_6514_crop.fits',0,h)
ssi = size(int)

int_s = gauss_smooth(int,1.5,kernel=3,/edge_truncate)

; region information for upper lya shits 
c1 = [[93,50],[93,40],[103,45],[110,46],[83,23],[91,20],[99,20],[107,20]];,[105,53]		; centers
s1 = [[12,8],[12,12],[8,10],[6,8],[8,10],[8,16],[8,16],[8,16]]	;,[8,8]	; spatial extent (pixels)

XC=(c1[0,*])*sxpar(h,'CD1_1') - sxpar(h,'CRPIX1')*sxpar(h,'CD1_1')+sxpar(h,'CRVAL1')
XS1=(c1[0,*]-(s1[0,*]+1)/2)*sxpar(h,'CD1_1') - sxpar(h,'CRPIX1')*sxpar(h,'CD1_1')+sxpar(h,'CRVAL1')
XS2=(c1[0,*]+(s1[0,*]+1)/2)*sxpar(h,'CD1_1') - sxpar(h,'CRPIX1')*sxpar(h,'CD1_1')+sxpar(h,'CRVAL1')
YC=c1[1,*]+4
YS=(s1[1,*])/2

;pos1=[0.1,0.1,0.9,0.9]
pos1=[0.1,0.35,0.9,0.65]
pos2=[0.1,0.65,0.9,0.95]

x=dindgen(ssi[1])*sxpar(h,'CD1_1') - sxpar(h,'CRPIX1')*sxpar(h,'CD1_1')+sxpar(h,'CRVAL1')
y=dindgen(ssi[2])+4 

;sc=dindgen(10)/10.
sc=[0.0,.05,.1,.15,.2]
sci=dblarr(5,3)

sci[*,0]=sc
sci[*,1]=sc
sci[*,2]=sc

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

set_plot,'ps'
device,/encap,filename='special_boxes2.ps',/color

loadct,0
tvimage,bytscl(int_s,min=0.,max=.003),POSITION=pos2,/keep_aspect_ratio
loadct,0
cgloadct,3,/reverse;, ncolors=100,clip=100,
tvimage,bytscl(bpol_img,min=0.,max=.45),POSITION=pos1,/keep_aspect_ratio
loadct,0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

plot,x,y,POSITION=pos2,/noerase,/nodata,thick=4,charthick=4,xthick=4,ythick=4,charsize=1.5,xstyle=1,xtickformat='(a2)',ystyle=13

for i=0,ss1[2]-1 do begin
oplot,[XS1[i],XS2[i],XS2[i],XS1[i],XS1[i]], $
      [YC[i]-YS[i],YC[i]-YS[i],YC[i]+YS[i],YC[i]+YS[i],YC[i]-YS[i]], $
      color=cgcolor('red'),thick=2
;xyouts,XC[i]-1,YC[i]-1,strtrim(i+1,2),charsize=.8,charthick=3,color=cgcolor('red'),/data
endfor

axis,4921.6,yaxis=0,yrange=[0.001,17.5],ystyle=1,charthick=4,charsize=1.5,ythick=4,/save
axis,5048.2,yaxis=1,ystyle=1,ythick=4,/nodata,ytickformat='(a2)',/save

oplot,dindgen(100)*0. + 4984.117,dindgen(100),thick=5,linestyle=2,color=1

xyouts,4904,-7,'Arcseconds',charsize=1.5,charthick=4,/data,orientation=90

;stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

plot,x,y,POSITION=pos1,/noerase,/nodata,thick=4,charthick=4,xthick=4,ythick=4,xtitle='Wavelength [Angstroms]',charsize=1.5,xstyle=1,ystyle=1,yr=[0.001,19.99]
;axis,4922,yaxis=0,yrange=[0.001,17.5],ystyle=1,charthick=4,charsize=1.5,ythick=4,/save
;axis,5048,yaxis=0,ystyle=1,charthick=4,charsize=1.5,ythick=4,/save

oplot,dindgen(100)*0. + 4984.117,dindgen(100),thick=5,linestyle=2,color=1
xyouts,5010,2,'S/N > 1.0',charthick=4,charsize=1.2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

cgloadct,3,/reverse;clip=100,
tvimage,bytscl(sci,min=0.,max=.45),POSITION=[0.15,0.12,0.9,0.2];,/noerase
loadct,0

plot,sc,sc/sc,POSITION=[0.15,0.12,0.9,0.2],/noerase,/nodata,thick=4,xtitle='P fraction',charsize=1.2,xstyle=1,ystyle=4,charthick=4,ytickformat='(a2)',xr=[0,.25]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

device,/close
set_plot,'X'
stop

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function polcalc2, a1,a2,a3,a4,b1,b2,b3,b4
	f1 = (a1-b1)/(a1+b1)
	f2 = (a2-b2)/(a2+b2)
	f3 = (a3-b3)/(a3+b3)
	f4 = (a4-b4)/(a4+b4)
	q = 0.5*f1 - 0.5*f3
	u = 0.5*f2 - 0.5*f4
	return, sqrt(q^2 + u^2)
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mc
	; read in ord/ext beams
	spawn,'ls stacked/crop/[A,B]*stack_7614_2dcrop.fits > stacked_sci.list'
	readcol,'stacked_sci.list',format='a',sci_files
	
	; read in the sig maps for each ord/ext beam
	spawn,'ls stacked/crop/[A,B]*2dstack2_6514_sig_crop.fits > stacked_sig.list'
	readcol,'stacked_sig.list',format='a',sig_files

	; read in test image for array sizes
	test=mrdfits(sci_files[0],0,h)
	ss=size(test)
	; let's start with the fully binned image : 
	bx = 8.
	by = 10.	
	; number of total bins in x and y directions
	num_x = ss[1]/bx		
	num_y = ss[2]/by

	spec = fltarr(ss[1],ss[2],n_elements(sci_files))
	spec_b = fltarr(num_x,num_y,n_elements(sci_files))
	sigs = fltarr(ss[1],ss[2],n_elements(sci_files))
	sigs_b = fltarr(num_x,num_y,n_elements(sci_files))

	; read in all the ord/ext beams and their corresponding sigma maps
	a1 = mrdfits(sci_files[0],0,h)
	a2 = mrdfits(sci_files[1],0,h)
	a3 = mrdfits(sci_files[2],0,h)
	a4 = mrdfits(sci_files[3],0,h)
	b1 = mrdfits(sci_files[4],0,h)
	b2 = mrdfits(sci_files[5],0,h)
	b3 = mrdfits(sci_files[6],0,h)
	b4 = mrdfits(sci_files[7],0,h)

	sa1 = mrdfits(sig_files[0],0,h)
	sa2 = mrdfits(sig_files[1],0,h)
	sa3 = mrdfits(sig_files[2],0,h)
	sa4 = mrdfits(sig_files[3],0,h)
	sb1 = mrdfits(sig_files[4],0,h)
	sb2 = mrdfits(sig_files[5],0,h)
	sb3 = mrdfits(sig_files[6],0,h)
	sb4 = mrdfits(sig_files[7],0,h)

	; rebin the data
	ra1 = rebin(a1,num_x,num_y)
	ra2 = rebin(a2,num_x,num_y)
	ra3 = rebin(a3,num_x,num_y)
	ra4 = rebin(a4,num_x,num_y)
	rb1 = rebin(b1,num_x,num_y)
	rb2 = rebin(b2,num_x,num_y)
	rb3 = rebin(b3,num_x,num_y)
	rb4 = rebin(b4,num_x,num_y)

	rsa1 = sqrt(rebin(sa1^2,num_x,num_y))/sqrt(bx*by)
	rsa2 = sqrt(rebin(sa2^2,num_x,num_y))/sqrt(bx*by)
	rsa3 = sqrt(rebin(sa3^2,num_x,num_y))/sqrt(bx*by)
	rsa4 = sqrt(rebin(sa4^2,num_x,num_y))/sqrt(bx*by)
	rsb1 = sqrt(rebin(sb1^2,num_x,num_y))/sqrt(bx*by)
	rsb2 = sqrt(rebin(sb2^2,num_x,num_y))/sqrt(bx*by)
	rsb3 = sqrt(rebin(sb3^2,num_x,num_y))/sqrt(bx*by)
	rsb4 = sqrt(rebin(sb4^2,num_x,num_y))/sqrt(bx*by)

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; Run Simulation
	n = 10000
	; calc actual polarization 
	orig_p = polcalc2(ra1,ra2,ra3,ra4,rb1,rb2,rb3,rb4)
	
	sim_p = dblarr(num_x,num_y,n)
	sigma1 = fltarr(num_x,num_y)
	sigma2=sigma1

	ff=dindgen(10)*1000.
	for j=1,n do begin
		; simulate a NEW value of each of the ord/ext beams by adding to that value a small deviation 
		; based on the a value of 
		ra1_s = ra1 + (randomu(seed,num_x,num_y)-0.5d)*rsa1*2
		ra2_s = ra2 + (randomu(seed,num_x,num_y)-0.5d)*rsa2*2
		ra3_s = ra3 + (randomu(seed,num_x,num_y)-0.5d)*rsa3*2
		ra4_s = ra4 + (randomu(seed,num_x,num_y)-0.5d)*rsa4*2
		rb1_s = rb1 + (randomu(seed,num_x,num_y)-0.5d)*rsb1*2
		rb2_s = rb2 + (randomu(seed,num_x,num_y)-0.5d)*rsb2*2
		rb3_s = rb3 + (randomu(seed,num_x,num_y)-0.5d)*rsb3*2
		rb4_s = rb4 + (randomu(seed,num_x,num_y)-0.5d)*rsb4*2

		sim_p[*,*,j-1] = polcalc2(ra1_s,ra2_s,ra3_s,ra4_s,rb1_s,rb2_s,rb3_s,rb4_s)
 
		fr=where(ff eq j)
        if fr[0] gt 0. then print,'done ',ff[fr]/10000.,'%'
	endfor

;stop


	for i=0,num_x-1 do begin
		for j=0,num_y-1 do begin 
;			set_plot,'PS'
;			device,filename='plots/hist_'+strtrim(i,2)+'_'+strtrim(j,2)+'.ps',/color,/encap

;			plothist,sim_p[i,j,*],xhist,yhist,bin=0.01,/fill,fcolor=230,xr=[0,1.],title= $
;						'P Hist for bin [' +strtrim(i,2)+', '+strtrim(j,2)+']', xtit = $
;						'Polarization fraction'
;			oplot,yhist*0.+ orig_p[i,j],yhist,thick=2,color=cgcolor('black') 
;
;			avg = mean(sim_p[i,j,*])
;			med = median(sim_p[i,j,*])
;			std = stddev(sim_p[i,j,*])
 ;    
	;		oplot,yhist*0.+ avg,yhist,thick=2, color=cgcolor('red')
	;		oplot,yhist*0.+ med,yhist,thick=2,color=cgcolor('green')
;
;			oplot,yhist*0.+(avg+std),yhist,thick=2,linestyle=1,color=cgcolor('blue')
;			oplot,yhist*0.+(avg-std),yhist,thick=2,linestyle=1,color=cgcolor('blue')
;			oplot,yhist*0.+(avg+2*std),yhist,thick=2,linestyle=1,color=cgcolor('magenta')
;			oplot,yhist*0.+(avg-2*std),yhist,thick=2,linestyle=1,color=cgcolor('magenta')
;	
			percent = PERCENTILES(sim_p[i,j,*],VALUE=[0.025,0.16,0.5,0.84,0.975])
;
 ;     ;oplot,yhist*0.+percent[2],yhist,thick=2,color=cgcolor('orange'),linestyle=3
  ;    oplot,yhist*0.+percent[1],yhist,thick=2,color=cgcolor('blue'),linestyle=2
   ;   oplot,yhist*0.+percent[3],yhist,thick=2,color=cgcolor('blue'),linestyle=2
    ;  oplot,yhist*0.+percent[0],yhist,thick=2,color=cgcolor('magenta'),linestyle=2
     ; oplot,yhist*0.+percent[4],yhist,thick=2,color=cgcolor('magenta'),linestyle=2
;
	  		sigma1[i,j] = (orig_p[i,j] - percent[1]) >  $
									(percent[3] - orig_p[i,j])
	  		sigma2[i,j] = (orig_p[i,j] - percent[0]) >  $
									(percent[4] - orig_p[i,j])

;stop
		;if i eq 12 and j eq 3 then stop

;		device,/close
;stop
		endfor
	endfor

set_plot,'X'

writefits,'MC_polstn1.fits',orig_p/sigma1
writefits,'MC_polstn2.fits',orig_p/sigma2
writefits,'MC_1sigma.fits',sigma1
writefits,'MC_2sigma.fits',sigma2
writefits,'ORIG_pol.fits',orig_p

stop

end

