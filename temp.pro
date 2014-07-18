

; NAME:		master_flat
; PURPOSE:	create master flat from bias-subtracted flats
;			i don't think i'm going to make my own
; INPUTS:	
pro master_flat, top, bot, list
	readcol,list,format='a',dirs
	for i=0,n_elements(dirs)-1 do begin
		cd,dirs[i]+'/calibration/flats/'
		; read in the five flats for each directory
		f1=readfits('flat1_bcor.fits',h1,exten_no=0)
		f2=readfits('flat2_bcor.fits',h2,exten_no=0)
		f3=readfits('flat3_bcor.fits',h3,exten_no=0)
		f4=readfits('flat4_bcor.fits',h4,exten_no=0)
		f5=readfits('flat5_bcor.fits',h5,exten_no=0)
		big_arr=[[[f1]],[[f2]],[[f3]],[[f4]],[[f5]]]
		; median combine them
		medarr, big_arr, medflat
		ss=size(medflat)
		temp = make_array(ss[1],ss[2],value=0.0)
		; look at each slit in the median flat
		;for j=0,n_elements(top)-1 do begin
		for j=0,ss[2]-1 do begin
			; fit the spectral response
			;row = top[i] - (top[i]-bot[i])/2
			x = dindgen(ss[1])
			;prof = medflat[*,row]
			prof = medflat[*,j]
			fit = poly_fit(x,prof,4)
			y = fit[0]+fit[1]*x+fit[2]*x^2+fit[3]*x^3
			plot,x,y
			;stop
			; divide out the spectral response
			;slit_h=top[j]-bot[j]+1
			;yslit = fltarr(ss[1],slit_h)
			;for k=0,slit_h-1 do begin
			;	yslit[*,k]=y
			;endfor
			;stop
			;temp[*,bot[j]:top[j]] = medflat[*,bot[j]:top[j]]/yslit
			;stop
			temp[*,j] = medflat[*,j]/y
		;stop
		endfor
		dir = strsplit(dirs[i],'=/',/extract)
		writefits,'masterflat3.fits',temp,h1
		i+=4
		cd,'../../../'
		stop
	endfor
end







; NAME:		find_xshift
; PURPOSE:	
; INPUTS:	list of reduced sci data -- for now accepts 1D spectra
; OUTPUTS:	The avg x value that 
function find_xshift, filelist
	readcol,filelist,format='a',files
	avgs = fltarr(n_elements(files))
	xvals = fltarr(n_elements(files))	
	for i=0,n_elements(files)-1 do begin		
		image=mrdfits(files[i],0,hdr)
		dir = strsplit(files[i],'/',/extract)
		name = strsplit(files[i],'_',/extract)
		ss=size(image)	
		x=dindgen(ss[1])

		; build the 4 profiles along the x-direction 
		;skyA = fltarr(ss[1])
		;skyB = fltarr(ss[1])
		;labA = fltarr(ss[1])
		;labB = fltarr(ss[1])

		; generate profiles of each slit along x direction
		;for p=0,ss[1]-1 do skyA[p]=median(image[p,570-50:570+50])
		;for p=0,ss[1]-1 do skyB[p]=median(image[p,478-50:478+50])
		;for p=0,ss[1]-1 do labA[p]=median(image[p,390-50:390+50])
		;for p=0,ss[1]-1 do labB[p]=median(image[p,300-50:300+50])

		; cut out the most likely location of the central-most sky line (in x-dir)
		; if slit1 or slit2 then it's a lab slit -- appropriate range for the skyline is 
		; between 950 and 1030
		if name[3] eq 'slit1' or name[3] eq 'slit2' then begin
			range = [950,1050]
			linecut = where(x gt 950 and x lt 1040)
		endif 
		; if slit1 or slit2 then it's a sky slit -- appropriate range for the skyline varies
		if name[3] eq 'slit3' or name[3] eq 'slit4' then begin
			; if the directory is N1_1, N1_2, or N2_1 then the range is 1050 to 1125
			if dir[0] eq 'N1_1' or dir[0] eq 'N1_2' or dir[0] eq 'N2_1' then begin
				range = [1050, 1150]
				linecut = where(x gt 1030 and x lt 1120)
			; otherwise the range is somewhere betwen 850 and 950
			endif else begin
				range = [850, 950]
				linecut = where(x gt 860 and x lt 960)
			endelse
		endif
		; fit gaussian to the sky lines
		result = gaussfit(x[linecut],image[linecut],params,nterms=6)
		;rlabB = gaussfit(x[linecut_lab],labB[linecut_lab],lBparams,nterms=6)
		;rskyA = gaussfit(x[linecut_sky],skyA[linecut_sky],sAparams,nterms=6)
		;rskyB = gaussfit(x[linecut_sky],skyB[linecut_sky],sBparams,nterms=6)

		plot,x,image,xr=range,psym=2
		oplot,x[linecut],result,color=250
		center=params[1]
		if i eq 0 then xcomp = center
		xvals[i] = xcomp-center
		print,files[i]
stop
	endfor
		openw,out,'xshifts.out',/get_lun
		printf,out,transpose(xvals)
		close,out	
		free_lun,out
	return,xvals
end



; NAME:		xalign
; PURPOSE:	Align the sci data in the x direction 
; INPUTS:	
pro xalign, filelist, avgs
	readcol, filelist, format='a', files
	name=strsplit(files[0],'/',/extract)
	for i=0,n_elements(files)-1 do begin
		image=mrdfits(files[i],0,h)
		ss = size(image)
		if ss[0] gt 1 then test=subpix_shift(image,avgs[i],0) else test=subpix_shift(image,avgs[i])
		sxaddpar,h,'XALIGN','T',' x-aligned to '+name[2],after='ALIGN',format='(a1)'
		dir = strsplit(files[i],'.',/extract)	
;stop
		writefits,dir[0]+'_x.fits',test,h
	endfor
end



pro sky_subtract, skylist, datlist
	readcol, skylist, format='a', sky
	readcol, datlist, format='a', lab
	for i=0,n_elements(lab)-1 do begin
		; read in the interpolated sky data
		dir = strsplit(lab[i],'/',/extract)
		names = strsplit(lab[i],'_',/extract)
		skyname = strsplit(sky[i],'.',/extract)
		readcol,skyname[0]+'_interp.out',wlength,skyval
		; read in the science data
		dat = mrdfits(lab[i],0,h)

		pix=sxpar(h,'CRVAL1')
		disp=sxpar(h,'CD1_1')
		x=dindgen(n_elements(dat))*disp+pix

		;plot,x,dat,xr=[5500,5700]
		;oplot,wlength,skyval,color=240
		result = dat - skyval
;		plot,x,result,xr=[4700,5100]	
		writefits,dir[0]+'/'+dir[1]+'/'+dir[0]+'_'+names[3]+'_'+names[4]+'_skysub.fits',result, h
	;stop
	endfor 
end

pro	create_tot_intens
	outfile='tot_i_calint.fits'
	spawn,'ls stacked/*calint.fits > stacked.list'
	readcol,'stacked.list',format='a',files
	test = mrdfits(files[0])
	ss = size(test)
	if ss[0] eq 1 then begin
		spec=fltarr(ss[1],n_elements(files))
		tot_is = fltarr(ss[1],4)
		tot = fltarr(ss[1])
		for i=0,n_elements(files)-1 do spec[*,i] = mrdfits(files[i],0,h)
		for i=0,3 do tot_is[*,i] = spec[*,i] + spec[*,i+4]
		for i=0,2047 do tot[i] = mean(tot_is[i,*])
		writefits,'tot_intens_1d.fits',tot
	endif else begin
		spec=fltarr(ss[1],ss[2],n_elements(files))
		tot_is = fltarr(ss[1],ss[2],4)
		tot = fltarr(ss[1],ss[2])
		for i=0,n_elements(files)-1 do spec[*,*,i]=mrdfits(files[i],0,h)
		for i=0,3 do tot_is[*,*,i] = spec[*,*,i] + spec[*,*,i+4]
		for i=0,ss[1]-1 do begin
			for j=0,ss[2]-1 do begin
				tot[i,j] = mean(tot_is[i,j,*])
			endfor
		endfor
stop
		writefits,outfile,tot,h
	endelse

	;x = dindgen(2048)
	;plot,x[500:800],tot[500:800],yr=[-0.004,0.008]
stop
end



;;;;;;;;;;;;;NOTE: USED IRAF INSTEAD;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NAME:		stack
; PURPOSE:	averages together a set of 1D or 2D spectra that have been coarsely aligned
;			(aligned by shifting them all wrt a bright sky line)
;			This means all their columns are already roughly aligned
; INPUTS:	none
; OUTPUTS:	8 stacked spectra: 1 for each position of the HWP and each side of the Wolly
pro stack
	prefix=['A00','B00','A22','B22','A45','B45','A67','B67']
	readcol,'multiplot.list',format='a',list
	filename='temp.list'
	for i=0,n_elements(list)-1 do begin
		spawn,'ls '+list[i]+' > '+filename
		readcol, filename, format='a', lis
		test = mrdfits(lis[0],0,h)
		ss = size(test)
		if ss[1] eq 1 then begin
			test=fltarr(ss[1],10)
			combo=fltarr(ss[1])
			x=dindgen(ss[1])
			for j=0,n_elements(lis)-1 do begin
				test[*,j]=mrdfits(lis[j],0,h)
			endfor
			for j=0,ss[1]-1 do combo[j]=mean(test[j,*])
			writefits,prefix[i]+'_1dstacked.fits',combo,h
			plot,x,combo,xr=[550,800]
		endif else begin
			test=fltarr(ss[1],ss[2],10)
			combo=fltarr(ss[1],ss[2])
			for j=0,n_elements(lis)-1 do begin
				test[*,*,j]=mrdfits(lis[j],0,h)
			endfor
			for j=0,ss[1]-1 do begin
				for k=0,ss[2]-1 do begin
					combo[j,k]=mean(test[j,k,*])
				endfor
			endfor
			writefits,prefix[i]+'_2dstacked_calib.fits',combo,h
		endelse
	;stop 
	endfor
end

pro stn_maps
	spawn,'ls stacked/*2dstack2_6514.fits > stacked.list'
	spawn,'ls stacked/*2dstack2_6514_sig.fits > stacksig.list'

	readcol,'stacked.list',format='a',sci
	readcol,'stacksig.list',format='a',sig

	for i=0,n_elements(sci)-1 do begin
		name=strsplit(sci[i],'.',/extract)
		tempsci=mrdfits(sci[i],0,h)
		tempsig=mrdfits(sig[i])
		stn=tempsci/tempsig
		writefits,name[0]+'_ppstn.fits',stn,h
	endfor
end


pro pol_2
	A00=mrdfits('A00_stacked.fits')
	A22=mrdfits('A22_stacked.fits')
	A45=mrdfits('A45_stacked.fits')
	A67=mrdfits('A67_stacked.fits')
	B00=mrdfits('B00_stacked.fits')
	B22=mrdfits('B22_stacked.fits')
	B45=mrdfits('B45_stacked.fits')
	B67=mrdfits('B67_stacked.fits')
	q = 0.5*((A00-B00)/(A00+B00))-0.5*((A45-B45)/(A45+B45))
	u = 0.5*((A22-B22)/(A22+B22))-0.5*((A67-B67)/(A67+B67))
	x=dindgen(2048)
	;plot,x,q,xr=[550,800]
	;oplot,x,u,color=240
	p= sqrt(q^2 + u^2)
	;oplot,x,p,color=cgcolor('green')

	set_plot,'PS'
	spawn,'ls *_stacked.fits > stacked.list'
	readcol,'stacked.list',format='a',list
	;device,filename='stacked_AB.eps',/encapsulated
	cgerase & multiplot,[2,4],mTitle="All data stacked by A/B & HWP angle",/doyaxis
	FOR i=0,n_elements(list)-1 DO BEGIN
		temp=mrdfits(list[i],0,h,/silent)
		x=dindgen(n_elements(temp))
		cgplot,x,temp,xr=[500,1299];,yr=[-10.,10.]
		multiplot
		;stop
	ENDFOR
	;endfor
	multiplot,/reset
	device,/close
stop
	;device, filename='stacked_QU2.eps',/encapsulated
	cgerase & multiplot,[1,3],mTitle='Stacked Q, U and P',/doyaxis
	cgplot,x,q,xr=[550,700];,yr=[-1.,1.]
	multiplot
	cgplot,x,u,xr=[550,700];,yr=[-1.,1.]
	multiplot
	cgplot,x,p,xr=[550,700];,yr=[0.,1.]
	multiplot,/reset
	;device,/close
	;set_plot,'X'
stop
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Polarization Calculation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; make Q and U for each night, for each set of images
; Let 'slit1' be A and 'slit2' be B

pro polarization
	dirlist = 'directories.list'
	readcol,dirlist,format='a',dir
	spawn,'pwd', wdir
	bin = 300
	
	for i=0,n_elements(dir)-1 do begin
		cd,dir[i]+'/1dspectra'
		; Read in files for Q
		A00 = mrdfits(dir[i]+'_slit1_00_1d_x_sub.fits',0,h1q)
		B00 = mrdfits(dir[i]+'_slit2_00_1d_x_sub.fits',0,h2q)
		A45 = mrdfits(dir[i]+'_slit1_45_1d_x_sub.fits',0,h3q)
		B45 = mrdfits(dir[i]+'_slit2_45_1d_x_sub.fits',0,h4q)
		; Rebin
		A00r = rebin(A00[0:1499],bin)
		B00r = rebin(B00[0:1499],bin)
		A45r = rebin(A45[0:1499],bin)
		B45r = rebin(B45[0:1499],bin)
		; Read in files for U	
		A22 = mrdfits(dir[i]+'_slit1_22_1d_x_sub.fits',0,h1u)
		B22 = mrdfits(dir[i]+'_slit2_22_1d_x_sub.fits',0,h2u)
		A67 = mrdfits(dir[i]+'_slit1_67_1d_x_sub.fits',0,h3u)
		B67 = mrdfits(dir[i]+'_slit2_67_1d_x_sub.fits',0,h4u)
		; Rebin
		A22r = rebin(A22[0:1499],bin)
		B22r = rebin(B22[0:1499],bin)
		A67r = rebin(A67[0:1499],bin)
		B67r = rebin(B67[0:1499],bin)
		; calculate Q and U
		q = 0.5*((A00-B00)/(A00+B00))-0.5*((A45-B45)/(A45+B45))
		u = 0.5*((A22-B22)/(A22+B22))-0.5*((A67-B67)/(A67+B67))
		qr = 0.5*((A00r-B00r)/(A00r+B00r))-0.5*((A45r-B45r)/(A45r+B45r))
		ur = 0.5*((A22r-B22r)/(A22r+B22r))-0.5*((A67r-B67r)/(A67r+B67r))

		x1=dindgen(n_elements(A00))
		x2=dindgen(bin)
		plot,x1,q,xr=[500,800]
		oplot,x1,u,color=240
	;stop
		plot,x2,qr,xr=[100,160]
		oplot,x2,ur,color=240	
	;stop

		writefits,dir[i]+'_coarse_q.fits',q
		writefits,dir[i]+'_coarse_u.fits',u
		writefits,dir[i]+'_coarse_qr.fits',qr
		writefits,dir[i]+'_coarse_ur.fits',ur
		cd,wdir		
	endfor
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Course Wavelength Alignment -- Using Gaussians on the Sky Lines 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro analyze_1d

; Terry says to test the 00 HWP for initial polarization signal
; 1. roughly align the sky and lab slits by the bright sky-line closest to ly-alpha
; 2. subtract the sky slit from the lab slit for both A/B sides

;spawn,'ls N*/1dspectra/N*slit*_1d.fits > 1d_slits.list'
;x = find_xshift('1d_slits.list')
;stop
;readcol,'xshifts.out',xshifts1
;spawn, 'ls N*/1dspectra/N*_1d.fits > 1d_slits.list'
;xalign,'1d_slits.list',xshifts1

spawn,'ls N*/1dspectra/N*slit[1,2]*_x.fits > labspec.list'
spawn,'ls N*/1dspectra/N*slit[3,4]*_x.fits > skyspec.list'

readcol,'labspec.list',format='a',labs
readcol,'skyspec.list',format='a',skys
for i=0,n_elements(labs)-1 do begin
	; coarse sky subtraction for 1D spectra
	name=strsplit(labs[i],'.',/extract)
	t1 = mrdfits(labs[i],0,h)
	t2 = mrdfits(skys[i],0,h)
	test = t1-t2

	x=dindgen(2048)
	plot,x,t1,xr=[0,1500],yr=[-0.02,0.1]
	oplot,test
;stop
	writefits,name[0]+'_sub.fits',test,h
endfor

stop

end

pro analyze_2d
; this is the same type of functionality as analyze_1d
; except for 2d spectra --
; align in the x directin and then sky subtract
; sky subtraction : 
;	take the median of the columns of the sky slits
;	multiply the median by the total number of rows in the slit (scaling factor)


;	readcol,'xshifts.out',xshifts1
	; use 'alignment values' found in the 1D case and apply them to the 2D spectra
;	spawn, 'ls N*/2dspectra/* > 2d_slits.list'
;	xalign,'2d_slits.list',xshifts1

	; need to sky subtract
	spawn, 'ls N*/2dspectra/N*slit[3,4]*_x.fits > skyslits.list'
	spawn, 'ls N*/2dspectra/N*slit[1,2]*_x.fits > labslits.list'
	readcol,'skyslits.list',format='a',skys
	readcol,'labslits.list',format='a',labs
	for i=0,n_elements(labs)-1 do begin
		s=mrdfits(skys[i])
		l=mrdfits(labs[i])
		lab_name=strsplit(labs[i],'.',/extract)
		sky_name=strsplit(skys[i],'.',/extract)
		ss = size(l)
		; extract the sky slit, s, using the median
		s_ext = fltarr(ss[1],ss[2])
		for j=0,ss[1]-1 do s_ext[j,*] = median(s[j,*])
		writefits,sky_name[0]+'_sub.fits',s_ext
		
		test1 = l-s_ext
		test2 = l-s

		writefits,lab_name[0]+'_subsky.fits',test1	
		writefits,lab_name[0]+'_submed.fits',test2	
	;stop
	endfor
stop
end

pro plotstuff
	set_plot,'PS'
	prefix=['A00','B00','A22','B22','A45','B45','A67','B67']
	; read in a list of all the different subsets of spectra
	readcol,'multiplot.list',format='a',list
	filename='temp.list'
	;for i=0,n_elements(list)-1 do begin
	; generate the subset of spectra to plot
	;spawn,'ls '+list[i]+' > '+filename
	spawn,'ls N*/1dspectra/N*_coarse_q.fits > '+filename
	; read in those file names for plotting
	readcol, filename, format='a', lis
	;device,filename=prefix[i]+'_x_sub.eps',/encapsulated
	device,filename='q.eps',/encapsulated
	;cgerase & multiplot,[2,5],mTitle=prefix[i]+" HWP pos angle",/doyaxis
	cgerase & multiplot,[2,5],mTitle="q (unbinned) for each set of data",/doyaxis
		FOR j=0,n_elements(lis)-1 DO BEGIN
			spec=mrdfits(lis[j],0,h,/silent)
			x=dindgen(n_elements(spec))
			cgplot,x,spec,xr=[550,800],yr=[-10.,10.]
			;cgplot,x,spec,xr=[110,160],yr=[-3.,3.]
			multiplot
			;stop
		ENDFOR
	;stop
	;endfor
	multiplot,/reset
	device,/close
	set_plot,'X'	
end


pro pol1d
	; define a bin size in the x/lambda direction
	bx = 8.
	; initial estimates to use in the gaussfit of the ly-a line
;	gest = [0.2,4981.,15,0.,0.,0.]	
	gest =[0.2,650.,10.,0.,0.,0.]

	; read in file lists
	spawn,'ls stacked/[A,B]*2dstack2_6514.fits > stacked_sci.list'
	readcol,'stacked_sci.list',format='a',sci_files	

	spawn,'ls stacked/I*2dstack2_6514.fits > stacked_int.list'
	readcol,'stacked_int.list',format='a',int_files	

	; read in a test science file to declare variables
	temp=mrdfits(sci_files[0],0,h)
	ss = size(temp)
	disp=sxpar(h,'CD1_1')
	pix=sxpar(h,'CRVAL1')
	wlen = dindgen(ss[1])*disp+pix
	y = dindgen(ss[2])
	; x range of 2d spectrum to look at
	xrange=where(wlen gt 4930. and wlen lt 5057.5)
	la = xrange[0]
	; y range of 2d spectrum to look at
	yrange=where(y gt 3 and y lt 74)
	; portion of the spectrum to use as the background
	srange=where(wlen gt 5025. and wlen lt 5500.)	
	; determine the number of bins based on wavelength range and bin size
	num_x = (n_elements(xrange))/bx	
	; declare a bunch of variables
	tot=fltarr(ss[1]+2,ss[2],n_elements(int_files))
	tot1d=fltarr(ss[1]+1,n_elements(int_files))		; 5 1d full Total Intensity spectra
	tot1d_b=fltarr(num_x,n_elements(int_files))		; 5 1d rebinned 1d spectra (cut)
	stn_b=tot1d_b									; stn for each tot_int bin
	sig = fltarr(n_elements(int_files))				; per pixel rms for each tot_int spectrum
	err = sig										; per bin error for each tot_int spectrum
	spec=fltarr(ss[1],ss[2],n_elements(sci_files))	; 8 2D full spectra 
	spec1d=fltarr(ss[1],n_elements(sci_files))		; 8 1D full spectra 
	spec1d_rb=fltarr(num_x,n_elements(sci_files))	; 8 rebinned 1D spectra (cut)
	spec1d_sb=spec1d_rb								; 8 sum-binned 1d spectra (cut)
	lya = fltarr(n_elements(sci_files))			; total integrated Lya line flux
	pix = lya
	center = pix
	stokes = fltarr(num_x,2)						; q,u pairs for each bin in the spectra
	pol = fltarr(num_x)								; polarization for each bin
	theta = pol										; theta for each bin
	qerr = pol										; polarization error for each bin

	avg=sig

	; make a binned versions of the total intensity images
	for i=0,n_elements(int_files)-1 do begin
		tot[*,*,i] = mrdfits(int_files[i],0,h)
		name=strsplit(int_files[i],'/.',/extract)
		for j=0,ss[1]-1 do tot1d[j,i]=total(tot[j,yrange])
		writefits,name[0]+'/'+name[1]+'_1d.fits',tot[*,i]
;		plot,wlen,tot1d[*,i],xr=[4930,5057.5]

		x = dindgen(ss[1])
		; mask out the ly-a line 
		mask = where((wlen lt 4960 and wlen gt 4700) or (wlen gt 5020 and wlen lt 5600))
		; fit the background with a linear polynomial
		coeffs = poly_fit(x[mask],tot1d[mask,i],1)
		y = coeffs[0] + coeffs[1]*x
		; subtract the linear fit from the spectrum
		flat = tot1d[*,i] - y
;		oplot,wlen,flat,color=240
		; rebin the 1d data
		tot1d_b[*,i] = rebin(flat[xrange],num_x,1)
		writefits,name[0]+'/rebin/'+name[1]+'1d_rebin.fits',tot1d_b[*,i]
		; calculate the RMS PER PIXEL of each science frame
		sig[i] = sigma(tot[srange,yrange,i])
		avg[i] = mean(tot[srange,yrange,i])
		; Bin Error: err = sqrt(number of pixels)*RMS PER PIXEL
		err[i] = sqrt(bx*n_elements(yrange))*sig[i]/bx
		; S/N per bin per angle for the total intensity images
		stn_b[*,i] = tot1d_b[*,i]/err[i] 
;stop
	endfor

;stop
	; read in all the stacked science spectra
	for i=0,n_elements(sci_files)-1 do begin
		spec[*,*,i] = mrdfits(sci_files[i],0,h)
		name = strsplit(sci_files[i],'/.',/extract)
		; extract the spectra -- SUM each column
		for j=0,ss[1]-1 do spec1d[j,i]=total(spec[j,yrange,i])
;stop
		x = dindgen(ss[1])
		; mask out the ly-a line 
		mask = where((wlen lt 4960 and wlen gt 4700) or (wlen gt 5020 and wlen lt 5600))
		; fit the background with a linear polynomial
		coeffs = poly_fit(x[mask],spec1d[mask,i],1)
		y = coeffs[0] + coeffs[1]*x
		; subtract the linear fit from the spectrum
		flat = spec1d[*,i] - y
;		plot,wlen,flat_t,xr=[4900,5400],yr=[-0.05,0.15]
;		oplot,wlen,flat_m*70,color=240
;stop

		; compute a gaussian for the lyman alpha line (centered at ~650 pixels
		result = gaussfit(x[xrange],flat[xrange],p,est=gest,nterms=6)
		plot,x[xrange],flat[xrange];,xr=[4900,5100];,yr=[-0.005,0.005]
	;	oplot,wlen,spec1d[*,i];,xr=[4900,5100];,yr=[-0.005,0.005]
		oplot,x[xrange],result,color=cgcolor('red')
		; integrate over entire ly-a line (2 sigma from center of line profile)
		lya[i] = total(flat[640-12:640+12])
		center[i] = p[1]
		pix[i] = (p[1]+p[2]) - (p[1]-p[2])
;stop
		
		; rebin the 1d data
		spec1d_rb[*,i] = rebin(flat[xrange],num_x,1)
		writefits,name[0]+'/rebin/'+name[1]+'1d_rebin.fits',spec1d_rb[*,i]
;stop
		; *sum* bin instead of rebinning (averaging)
		for j=0,num_x-1 do begin
			spec1d_sb[j,i] = total(flat[la+j*bx:la+(j+1)*bx])
		endfor
		writefits,name[0]+'/rebin/'+name[1]+'1d_sumd.fits',spec1d_sb[*,i]
	endfor

Rq= sqrt(lya[0]/lya[4])/sqrt(lya[2]/lya[6])
Ru= sqrt(lya[1]/lya[3])/sqrt(lya[5]/lya[7])
q=(Rq-1)/(Rq+1)
u=(Ru-1)/(Ru+1)
pola=sqrt(q^2+u^2)
print,pola
err=sqrt(24*n_elements(yrange))*sig[0]
print,err

q2=0.5*((lya[0]-lya[4])/(lya[0]+lya[4])) - 0.5*((lya[2]-lya[6])/(lya[2]+lya[6]))
u2=0.5*((lya[1]-lya[5])/(lya[1]+lya[5])) - 0.5*((lya[3]-lya[7])/(lya[3]+lya[7]))
pol2=sqrt(q2^2+u2^2)
print,pol2
stop	
stokes2=stokes
pol2=pol
qerr2=qerr
theta=pol
theta2=theta
;stop
	; calculate polarization for each bin
	for i=0,num_x - 1 do begin
		; REBINNED -- Per Pixel Errors -- Tinbergen Pol Calc
		stokes[i,*] = calc_stokes1(spec1d_rb[i,*])
		pol[i] = sqrt(stokes[i,0]^2+stokes[i,1]^2)
		theta[i] = 0.5*atan(stokes[i,1]/stokes[i,0])
;		qerr[i] = calc_error(stokes[i,0],err_t,spec1d_rb[i,*])/bx
;		stn[i] = calc_stn(stn_b[i,0],pol[i])

		; REBINNED -- Per Pixel Errors -- "Method 2"  
		stokes2[i,*] = calc_stokes3(spec1d_rb[i,*])		
		pol2[i] = sqrt(stokes2[i,0]^2+stokes2[i,1]^2)
		theta2[i] = 0.5*atan(stokes2[i,1]/stokes2[i,0])	
;		qerr2[i] = calc_error(stn_b[i,0])
	endfor

stn=sqrt(2)*pol*stn_b[*,3]
stn2=sqrt(2)*pol2*stn_b[*,3]

stncut=3.0

pol_p=pol*0. - 100
pol_p[where(stn gt stncut)]=pol[where(stn gt stncut)]

perr=pol_p/stn

pol_p2=pol2*0.-100
pol_p2[where(stn2 gt stncut)]=pol2[where(stn2 gt stncut)]

perr2=pol_p2/stn2

perr[where(perr lt 0.)]=0.
perr[where(perr gt 1.)]=0.

perr2[where(perr2 lt 0.)]=0.
perr2[where(perr2 gt 1.)]=0.

; plot total everything image 
; (the last one in the intensity list is the combination of all 80 frames)
plot,wlen,flat,xr=[4930,5057.5],yr=[0.,.2]
for k = 0,num_x-1 do begin
	oplot,dindgen(100)*0. + wlen(xrange(k*bx)),dindgen(100),linestyle=5
endfor
oplot,dindgen(200)*0. + 4981.5 ,dindgen(200),thick=1,linestyle=2;,color=cgcolor('blue')
axis,4930,yaxis=0,yrange=[0,1],ystyle=1,charthick=1,charsize=1.5,ythick=1,ytitle='Intensity [arbitrary units]',/save
axis,5057.5,yaxis=1,yrange=[0,0.4],ystyle=1,charthick=1,charsize=1.5,ythick=1,ytitle='Polarization Fraction',color=cgcolor('dark green'),/save
; overplot the polarization fraction 
circsym
;stop
xp=dindgen(n_elements(pol))*0.636*bx + wlen(xrange[0]-bx/2);
oploterror,xp,pol_p,perr,linestyle=6,psym=8,symsize=4,errthick=1,color=cgcolor('blue')
oploterror,xp,pol_p2,perr2,linestyle=6,psym=8,symsize=4,errthick=1,color=cgcolor('dark green')
;oploterror,xp,pol_p,perr,linestyle=6,psym=8,symsize=4,color=cgcolor('blue')
;oploterr0r,xp,pol_p2,perr2,linestyle=6,psym=8,symsize=4,color=cgcolor('dark green')


stop
	
;	writefits,'q2_'+strtrim(fix(bx),2)+'_1d.fits',stokes2[*,0]
;	writefits,'u2_'+strtrim(fix(bx),2)+'_1d.fits',stokes2[*,1]
	writefits,'pol_products/pol_'+strtrim(fix(bx),2)+'_R_cstn_Tin.fits',pol
;	writefits,'pol_products/perr_'+strtrim(fix(bx),2)+'_RppM2.fits',qerr2
	writefits,'pol_products/stn_'+strtrim(fix(bx),2)+'_R_cstn_Tin.fits',stn
	writefits,'pol_products/pol_'+strtrim(fix(bx),2)+'_R_cstn_ESO.fits',pol2
	writefits,'pol_products/stn_'+strtrim(fix(bx),2)+'_R_cstn_ESO.fits',stn2
stop

end
