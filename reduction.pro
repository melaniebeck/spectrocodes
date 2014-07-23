; NAME: 	prep_scidata
; PURPOSE: 	create a ROTANG key word in the header, trim overscan,
;			rename files appropriately
; INPUTS: 	List = list of file names to adjust
pro prep_scidata, datalist
; read in list of file names
	readcol, datalist, format='a', dat
	for i = 0,n_elements(dat)-1 do begin
		; read in the image/header
		img = readfits(dat[i],h,exten_no=0,/silent)
		; find the rotation angle in the messed up header file
		rotang = h(152)
		rotang = strsplit(rotang,'=/',/extract)
		; rewrite it as an actual KEYWORD at top of header
		sxaddpar,h,'ROTANG',double(rotang[1]),rotang[2], $
				 after='CD2_2',format='(d20.1)'
		; trim overscan region
		img = img[*,5:1028]
		; rename file appropriately
		dir = strsplit(dat[i],'/',/extract)
		writefits,dir[0]+'/'+dir[0]+'_'+strtrim(rotang[1],2)+'.fits',img, h
	endfor
end

; NAME: 	master_bias
; PURPOSE: 	creates median combined master bias for each data set (overscan trimmed)
; INPUTS: 	List = list of bias frames
pro master_bias, biaslist
	readcol,biaslist,format='a', bias
	for i=0,n_elements(bias)-1 do begin
		b1=readfits(bias[i],h,exten_no=0)
		b2=readfits(bias[i+1],h,exten_no=0)
		b3=readfits(bias[i+2],h,exten_no=0)
		b4=readfits(bias[i+3],h,exten_no=0)
		b5=readfits(bias[i+4],h,exten_no=0)
		big_arr=[[[b1]],[[b2]],[[b3]],[[b4]],[[b5]]]
		;master = median([b1,b2,b3,b4,b5],dimension=3)
		medarr, big_arr, master
		; trim overscane region
		master = master[*,5:1028]
		dir = strsplit(bias[i],'/',/extract)
		writefits,dir[0]+'/'+dir[1]+'/'+dir[2]+'/masterbias.fits',master
		i+=4
;	stop
	endfor
end

; NAME:		sci_bias_subtract
; PURPOSE:	subtract master bias from science data
; INPUTS:	list of science data to bias subtract
; NOTES:	if SCIENCE keyword is set, routine runs on science data;
;			if FLATS keyword is set, routine runs on flat data. duh. 
pro bias_subtract, dirlist, SCIENCE=science, FLATS=flats
	; read in directory list
	readcol,dirlist,format='a',dir
	; go into each directory and bias subtract
	for i=0,n_elements(dir)-1 do begin
		; descend into appropriate directory
		cd,dir[i]
		; determine whether to bias subtract science data or flats
		if keyword_set(science) then begin
			; read in the appropriate masterbias
			bias=readfits('calibration/sci_bias/masterbias.fits',hb,exten_no=0,/silent)
			; list the .fits files to apply the bias to
			spawn,'ls N*_{00,22,45,67}.fits > sci_data.lis'
			readcol,'sci_data.lis', format='a',imgs
		endif else begin 
			; read in the appropriate masterbias
			bias=readfits('calibration/flat_bias/masterbias.fits',0,hb)
			; list the .fits files to apply the bias to
			spawn,'ls calibration/flats/*FLAT*.fits > flats.lis'
			readcol,'flats.lis', format='a',imgs
		endelse
		for j=0,n_elements(imgs)-1 do begin
			; read in the image
			img=readfits(imgs[j],h,exten_no=0,/silent)
			; subtract the master bias
			sub=img-bias
			;adjust the header
			sxaddpar,h,'BIAS','T',' Bias corrected', $
					 after='ROTANG',format='(a1)'
			if keyword_set(flats) then begin
				name=strsplit(imgs[j],'/',/extract)
				writefits,'calibration/flats/flat'+strtrim(j+1,2) + $
						  '_bs.fits',sub,h
			endif else begin
				name=strsplit(imgs[j],'.',/extract)
				writefits,name[0]+'_bs.fits',sub,h
			endelse
		endfor
		cd,'../'
	endfor
end

; NAME:		flatfield
; PURPOSE:	divides the master flat into the science imgs
; INPUTS:	directory list
pro flatfield, dirlist
	readcol,dirlist,format='a',dir
	for i=0,n_elements(dir)-1 do begin
		; descend into proper directory
		cd,dir[i]		
		; read in the master flat
		flat = readfits('calibration/master_norm_flat_pmos.fits',hf,exten_no=0,/silent)
		; list the .fits files to apply the bias to
		spawn,'ls N*_{00,22,45,67}_bs.fits > sci_data.lis'
		readcol,'sci_data.lis', format='a',imgs						

		; loop through all four hwp angles in the directory
		for j=0,n_elements(imgs)-1 do begin 	
			; read in each science img
			img = readfits(imgs[j],h,exten_no=0,/silent)
			; divide the flat into each science img
			temp = float(img)/flat
			; adjust sci header to reflect flatfielding
			sxaddpar,h,'FLAT','T',' Flatfielded',after='BIAS',format='(a10)'
			; save flatfielded sci img
			name = strsplit(imgs[j],'.',/extract)
			writefits,name[0]+'_ff.fits',temp,h
		endfor
		; back up into main working directory
		cd,'../'
;stop
	endfor
end

; NAME:		add_headers
; PURPOSE:	Adds headers to the cosmic ray removed science data because I can't
;			figure out how to do this in the already existing python code
; INPUTS:	list of cosmic ray removed science data, list of flatfielded data 
;			(data that has appropriate headers to transfer over)
pro fix_crc_headers, tofixlist, gethdrlist
	readcol, tofixlist, format='a', files_to_fix
	readcol, gethdrlist, format='a', matching_hdrs
	for i=0,n_elements(files_to_fix)-1 do begin
		; read in GOOD HEADER
		get_hdr = mrdfits(matching_hdrs[i],0,hgood,/silent)
		; adjust GOOD HEADER to reflect cosmic ray removal
		sxaddpar,hgood,'CRREM','T',' Cosmic rays removied with L.A.Cosmic', $
				 after='FLAT', format='(a10)'
		; read in cosmic ray removed science data with BAD HEADER
		give_hdr = mrdfits(files_to_fix[i],0,hbad,/silent)
		; write CRREM  science data with GOOD HEADER
		writefits,files_to_fix[i],give_hdr,hgood
	endfor
end

; NAME:		normalize
; PURPOSE:	divide each science image by its exposure time
; INPUTS:	list of files to normalize
pro normalize, filelist
	readcol, filelist, format='a', files
	for i=0,n_elements(files)-1 do begin
		img=mrdfits(files[i],0,h)
		exptime = sxpar(h,'EXPTIME')
		temp = img/exptime
		; adjust sci header to reflect flatfielding
		sxaddpar,h,'NORM','T',' Normalized to EXPTIME', $
				 after='CRREM',format='(a10)'
		; save flatfielded sci img
		name = strsplit(files[i],'.',/extract)
		writefits,name[0]+'_n.fits',temp,h
	endfor
end

; NAME:		find_yshift
; PURPOSE:	Fits gaussians to the 4 stellar continua in each img (not night 1 
;			imgs --	they don't have stellar continua). Finds the center of 
;			each gaussian. Subtracts these central pixel values from those 
;			obtained in the first img on the list. Takes a mean of these four 
;			values for each img as an averge y shift necessary to align the 
;			sci data in the y direction.
; INPUTS:	list of reduced sci data with stellar continua
; OUTPUTS:	The avg y value that each img must be shifted wrt the first img on the list.
function find_yshift, filelist
	readcol,filelist,format='a',files
	avgs = fltarr(n_elements(files))
	for i=0,n_elements(files)-1 do begin		
		image=mrdfits(files[i],0,hdr,/silent)
		name = strsplit(files[i],'/',/extract)
		; if this is the first night's data, skip it! (NO STELLAR CONTINUA)
		if name[0] eq 'N1_1' or name[0] eq 'N1_2' then begin
			avgs[i] = 0		; should be 0 for N1_1, N1_2
		endif else begin
			ss=size(image)	
			; placeholder for the central row pixels of the 4 stellar
			; continuua 
			;per sci img
			yvals = fltarr(4)
			; profile array with elements = # of rows
			prof=dblarr(ss[2])
			; fill the profile array: 
			; for each row in the image take the median of 200 columns
			; around 
			; x = 650 pixels
			for p=0,ss[2]-1 do prof[p]=mean(image[650-150:650+150,p])	
			; make a second array with values of the number of rows		
			xx1=dindgen(ss[2])
			; for each of the stellar continua in the image ...
			for j=0,3 do begin
				case j of
				; define the pixel rows in which the stellar continuum is
				; most likely to be
				0:	cont=where(xx1 gt 80 and xx1 lt 150)	; continuum 
															; in slit 10
				1:	cont=where(xx1 gt 175 and xx1 lt 245)	; in slit 9
				2: 	cont=where(xx1 gt 625 and xx1 lt 700)	; in slit 4
				3:	cont=where(xx1 gt 725 and xx1 lt 795)	; in slit 3
					ENDCASE
				; plot the row values against the median of the columns
				; defined by "cont"
				plot,xx1[cont],prof[cont],psym=5
				; fit a gaussian to these data 
				result = gaussfit(xx1[cont],prof[cont],yparams, nterm=6, $
						 sigma=sig, yerror=yerr)
				oplot,xx1[cont],result
				; save central pixel of the gaussian to an array
				yvals(j)=yparams[1]
				; compare everything to the first img in the list 
				if i eq 8 then ycomp = yvals
			endfor
			; compare every other img to the first img
			shift = ycomp - yvals
			; take an avg of the 4 shifts
			avgs[i] = mean(shift)		; should be 0 for i = 0
;stop
		endelse
	endfor
stop
	return, avgs
end

; NAME:		yalign
; PURPOSE:	Align the sci data in the y direction 
; INPUTS:	list of reduced sci data to shift, y value each img should be
;			shifted by. These two must have equal elements - 1 y value for
;			EACH data frame in the list	in the same order as the data list
pro yalign, filelist, shifts
	readcol,filelist,format='a',files
	for i=0,n_elements(files)-1 do begin
		img=mrdfits(files[i],0,h,/silent)
		; sub pixel shift each img in y direction
		temp=shift_sub(img,0,shifts[i])		
		; adjust header to reflect the alignment process
		sxaddpar,h,'ALIGN','T',' y-aligned to '+files[0], after='NORM', $
				 format='(a1)'
		name = strsplit(files[i],'.',/extract)
		writefits,name[0]+'_y.fits',temp,h
	endfor
end

; NAME:		cut_slits
; PURPOSE:	Cut out individual slits from the multi-spectra image 
; INPUTS:	List of sci data to be manipulated, list of y-values (paired:
;			bottom, top) for each slit. This list of y-vals is applied to ALL
;			the sci images in the list -- different	cut values for different
;			sci data is not allowed in this procedure. The y-val list can
;			contain up to 11 rows of y values (there are 11 slits in this
;			science data) but the same list will be applied to all images. 
;			Saves each cut slit: NX_X_slitX_HWP.fits
pro cut_slits, filelist, cutvals
	readcol, filelist, format='a', files
	for i=0,n_elements(files)-1 do begin
		img=mrdfits(files[i],0,h,/silent)
		for j=0,n_elements(cutvals[0,*])-1 do begin
			slit = img[*,cutvals[0,j]:cutvals[1,j]]
			; adjust header to reflect the cut
			dir = strsplit(files[i],'./',/extract)
			ang = strsplit(files[i],'_',/extract)
			; if it's a comp lamp then who cares
			if n_elements(ang) lt 3 then begin
				writefits,dir[0]+'/lampslit'+strtrim(j+1,2)+'.fits',slit,h
			; if it's science data then adjust the filename/header info
			endif else begin
				sxaddpar,h,'CUT','T',' Rows '+strtrim(cutvals[0,j],2)+ $
						 ' to '+ strtrim(cutvals[1,j],2)+' of '+dir[1]+ $
						 '.fits'
				writefits,dir[0]+'/2dspectra/'+dir[0]+'_slit'+ $
						  strtrim(j+1,2)+'_'+ang[3]+'.fits',slit,h
			endelse
		endfor
	endfor
end

; NAME:		adjust_headers
; PURPOSE:	If many .fits files need to have a the same new keyword added to
;			their headers, this will do it!
; INPUTS:	list of .fits files that need headers updated, the keyword toadd 
;			to the header, the value for that keyword
pro adjust_headers, listoffiles, keyword, value
	readcol,listoffiles,format='a',files
	for i=0,n_elements(files)-1 do begin
		dat=mrdfits(files[i],0,h,/silent)
		sxaddpar,h,keyword,value
		writefits,files[i],dat,h
	endfor
end

; NAME:		extract_spectra
; PURPOSE:	 
; INPUTS:	
pro extract_spectra, filelist
	readcol,filelist,format='a',files
	for i=0,n_elements(files)-1 do begin
		img=mrdfits(files[i],0,h,/silent)
		ss=size(img)
		oned=fltarr(ss[1])
		for j=0,ss[1]-1 do begin
			oned[j]=median(img[j,*])
		endfor
		dir = strsplit(files[i],'.',/extract)
		writefits,dir[0]+'_1d.fits',oned,h
	endfor
end


; NAME:		build_dither_lists
; PURPOSE:	For 1D spec, build "dither" lists of spectra that have the same
;			wavelength pattern to apply the proper IRAF functions to them
; INPUTS:
pro build_dither_lists, dirlist
	readcol, dirlist, format='a', dirs
	for i=0,n_elements(dirs)-1 do begin
		spawn,'ls '+dirs[i]+'/1dspectra/* > '+dirs[i]+'/'+dirs[i]+'.list'
	endfor
	spawn,'cat N1_1/N1_1.list N1_2/N1_2.list N2_1/N2_1.list >> dither1.list'
	spawn,'cat N3_1/N3_1.list N4_1/N4_1.list N5_2/N5_2.list >> dither2.list'
	spawn,'cat N2_2/N2_2.list N3_2/N3_2.list N4_2/N4_2.list N5_1/N5_1.list >> dither3.list'
end

; NAME:		build_dispcor_cl
; PURPOSE:	Builds .cl file to run in IRAF that will execute a series of
;			dispcor commands: Puts the file name as the file to excute
;			dispcor on, tacks '_calib' onto the file name for the output of
;			dispcor
; INPUTS:	List of .fits files to perform dispcor on; the database directory 
; 			appropriate for the file set - all files must have the same
;			database directory!
pro build_dispcor_cl,listoffiles,data_dir
	readcol,listoffiles,format='a',files
	dith=strsplit(listoffiles,'.',/extract)
	openw,out,'dispcor_'+dith[0]+'.cl',/get_lun
	for i=0,n_elements(files)-1 do begin
		name=strsplit(files[i],'.',/extract)
		printf,out,'dispcor '+files[i]+' '+name[0]+'_calib.fits databas='+ $
			data_dir+'/database/'
	endfor
	close,out
	free_lun, out
;stop
end

; NAME:		interp_1d
; PURPOSE:	Based on a python script that Mike wrote: 
;			Interpolates appropriate sky values onto the wavelength range of
;			the	corresponding lab slit. For example, Ordinary beam of the sky
;			and lab on a given CCD image are not perfectly wavelength
;			correlated. Use this to	determine appropriate sky values at the
;			exact wavelengths of the lab slit
; INPUTS:	List of 1d sky spectra; list of 1d sci spectra. These lists must
;			be correlated: the first sky spectrum will be interpolated onto
;			the wavelength range of the first science spectrum, etc. 
; OUTPUTS:	An output file with naming format "input_sky_file_name_interp.out"
;			for each spectra pair: the first column contains wavelengths
;			the second column contains the interpolated sky values for those
;			wavelengths
pro interp_1d, skylist, datlist
	readcol, skylist, format='a', sky
	readcol, datlist, format='a', lab
	for i=0,n_elements(sky)-1 do begin
		dir=strsplit(sky[i],'.',/extract)
		spawn,'python interp.py '+sky[i]+' '+lab[i]+' '+dir[0]+'_interp.out 0'
	;stop
	endfor
end

; NAME:		oned_spectra
; PURPOSE:	Run the the final calibration steps on the 1d spectra:
;			Adjust headers -- add REFSPEC1
;			Create dither.cl lists for quick IRAF use
;			Interpolate the 1d sky spectra onto the wavelength range of the 
;			science data
; INPUTS:	
pro oned_spectra
; THIS SHIT NEEDS TO BE REWRITTEN BEFORE USING IT AGAIN. SUPER FUCKED. ><
	spawn,'ls N*/1dspectra/N*slit[1,2]*.fits > 1dlab_spec.list'
	spawn,'ls N*/1dspectra/N*slit[3,4]*.fits > 1dsky_spec.list'
	adjust_headers,'2dspec_master.list','DISPAXIS',1
	adjust_headers,'1dsky_spec.list','REFSPEC1','lampslit2'
	;stop 
	
	; Build IRAF script to run DISPCOR, applying calibration to sci data
	; Dither1: N1_1, N1_2, N2_1
	; Dither2: N2_2, N3_2, N4_2, N5_1
	; Dither3: N3_1, N4_1, N5_1
	if FILE_TEST('dither1.lis') eq 0 then build_dither_lists, dirlist
	
	build_dispcor_cl,'dither1.list','N1_1'
	build_dispcor_cl,'dither2.list','N2_2'
	build_dispcor_cl,'dither3.list','N3_1'
	
	print,'Go run DISPCOR in IRAF now. Then come back.'
	stop
	
	; Now Interpolate
	;spawn,'ls N*/1dspectra/N*slit[1,2]*1d_calib.fits > 1dlab_spec.list'
	;spawn,'ls N*/1dspectra/N*slit[3,4]*1d_calib.fits > 1dsky_spec.list'
	interp_1d,'1dsky_spec.list','1dlab_spec.list'
end


; NAME:		build_transform_cl
; PURPOSE:	Builds a .cl script to run TRANSFORM in an IRAF xgterm
; INPUTS:
function transform_list
	; if the list of imgs to be transformed hasn't yet been made -- make it!
 	if FILE_TEST('transform.lis') eq 0 then begin
		spawn, 'ls '+dir+'/2dspectra/N*00.fits > 00.list'
		spawn, 'ls '+dir+'/2dspectra/N*22.fits > 22.list'
		spawn, 'ls '+dir+'/2dspectra/N*45.fits > 45.list'
		spawn, 'ls '+dir+'/2dspectra/N*67.fits > 67.list'
		spawn, 'cat 00.list 22.list 45.list 67.list >> transform.lis'
	endif
	return,'transform.lis'
end

; NAME:		build_transform_cl
; PURPOSE:	Builds a .cl script to run TRANSFORM in an IRAF xgterm
; INPUTS:
pro build_transform_cl, trans_list
	; read in the list of imgs to transform
	readcol,trans_list,format='a',files
	openw,out,'transform.cl',/get_lun
	for i=0,n_elements(files)-1 do begin
		dir = strsplit(files[i],'/',/extract)
		name = strsplit(files[i],'.',/extract)
		slit = strsplit(files[i],'_',/extract)
		; need to associate each img with it's proper LAMPSLIT[1,2,3,4]
		; AND DATABASE
		case slit[3] of
			'slit1' : lampname='lampslit1'
			'slit2' : lampname='lampslit2'
			'slit3' : lampname='lampslit3'
			'slit4' : lampname='lampslit4'
		endcase
		case dir[0] of
			'N1_1'	: database='N1_1/database/'
			'N1_2'  : database='N1_1/database/'
			'N2_1'	: database='N2_1/database/'
			'N2_2'  : database='N2_2/database/'
			'N3_2'  : database='N2_2/database/'
			'N4_2'  : database='N2_2/database/'
			'N5_1' 	: database='N2_2/database/'
			'N3_1' 	: database='N3_1/database/'
			'N4_1'	: database='N3_1/database/'
			'N5_2'  : database='N3_1/database/'
		endcase
		; print the IRAF command to file
		printf,out,'transform '+files[i]+' '+name[0]+'_wc.fits databas='+ $
			  	   database+' fitname='+lampname
	endfor
	close,out
	free_lun,out
end

; NAME:		interp_2d
; PURPOSE:	performs 2d interpolation of spectra. Interpolates flux values
;			of all spectra onto the wavelength & spatial range of the first
;			spectrum in the list. Use this after you've sky subtracted to
;			ensure that all the spectrum have flux values for the exact same
;			wavelength values before stacking them.
; INPUTS:	none
; OUTPUTS:	the same 2d spectrum interpolated onto different wavelength values
;			specified by the first spectrum in the list
pro interp_2dcd
	spawn,'ls N*/2dspectra/N*_skysub.fits > 2d_calib_skysub.list'
	readcol,'2d_calib_skysub.list',format='a',list
	; read in the first spectra on the list -- use it's wavelength to run 
	; all other interpolations against
	comp=mrdfits(list[0],0,hc,/silent)
	ss=size(comp)
	disp=sxpar(hc,'CD1_1')
	pix=sxpar(hc,'CRVAL1')
	; define the x and y values to interpolate onto (wavelength and index)
	x1 = dindgen(ss[1])*sxpar(hc,'CD1_1')+sxpar(hc,'CRVAL1')
	y1=dindgen(ss[2])
	; run through the rest of the list
	for i=1,n_elements(list)-1 do begin
		temp=mrdfits(list[i],0,ht,/silent)
		; definte the x and y values to interpolate from
		x0=dindgen(ss[1])*sxpar(ht,'CD1_1')+sxpar(hc,'CRVAL1')
		y0=y1
		; call interp2d to perform binlinear interpolation on temp
		result=interp2d(temp,x0,y0,x1,y1,grid=1)
		if i eq 1 then begin
;			comp=comp[100:1949,*]
			name=strsplit(list[0],'.',/extract)
;			sxaddpar,hc,'CRVAL1',100*disp+pix
			writefits,name[0]+'_interp2.fits',comp,hc
		endif
;		result=result[100:1949,*]
		; save result as a new .fits
		name=strsplit(list[i],'.',/extract)
		; update each header with the new wavelength information -- 
		; that which each spectra has been interpolated onto
;		pix2=sxpar(ht,'CRVAL1')
		; set the new CRVAL to account for the chopping off of 100 pixels
;		sxaddpar,ht,'CRVAL1',100*sxpar(ht,'CD1_1')+pix2
		writefits,name[0]+'_interp2.fits',result,ht
	;stop
	endfor
end

function extract_airm, header
	m = intarr(n_elements(header))
	for i=0,n_elements(header)-1 do begin
		m[i] = strpos(header[i],"AIRM ")
	endfor
	matches = where(m ne -1)
	temp = strsplit(header[matches[0]],'=/',/extract)
	airm_start = temp[1]
	temp = strsplit(header[matches[1]],'=/',/extract)
	airm_end = temp[1]
	return, [airm_start, airm_end]
end

pro correct_atmos, filelist, TWOD=twod
	; read in airmass coefficient data
	readcol,'airmass_coeff.dat',kwlen,kmag
	; read in the list of spectra to correct
	readcol, filelist, format='a', files
	airm = fltarr(n_elements(files))
	for i=0,n_elements(files)-1 do begin
		img = mrdfits(files[i],0,h,/silent)
		ss = size(img)
		if KEYWORD_SET(twod) then temp = fltarr(ss[1], ss[2]) $
			else temp = fltarr(ss[1])
		; extract airmass values from header 
		airmass = extract_airm(h)
		; take the median airmass value
		airm[i] = mean(float(airmass))
		; build the wavelength values of the spectrum (interpolating ONTO)
		wlen = dindgen(ss[1])*sxpar(h,'CD1_1') +sxpar(h,'CRVAL1')
		; interpolate the airm coefficients onto these wlen vals
		k = interpol(kmag,kwlen,wlen,/spline)
;		img_m = -2.5*ALOG10(img)
		; m_orig = m_obs - k*airm
;		if KEYWORD_SET(twod) then for j=0,ss[2]-1 do temp1[*,j] = img_m[*,j] - k1*airm $
;			else temp1 = img_m - k1*airm
		; convert spectrum back into flux 
;		temp1 = 2.512^(-temp1)
		if KEYWORD_SET(twod) then for j=0,ss[2]-1 do temp[*,j] = img[*,j]*10^(.4*k*airm[i]) $
			else temp = img*10^(.4*k*airm[i])  
		sxaddpar,h,'ATMOS','T',' Corrected for atmospheric extinction', $
				 after='ALIGN', format='(a10)'
		name = strsplit(files[i],'.',/extract)
;stop
		writefits, name[0]+'_acorr3.fits',temp,h	
	endfor
	openw,out,'airmass.dat',/get_lun
	printf,out,transpose(airm)
	close,out
	free_lun,out
end

pro sky_subtract, scilist, skylist
	readcol, scilist, format='a', labs
	readcol, skylist, format='a', skys
	
	; Check to see if the 1d spectra have already been extracted
	if not FILE_TEST('1dsky_data.lis') then begin 
		; extract the wavelength calibrated sky spectra
		extract_spectra, skylist
		spawn, 'ls N*/2dspectra/N*wc_1d.fits > 1dsky_data.lis'
	endif
	; read in the 1d wlen calib sky spectra
	readcol, '1dsky_data.lis', format='a', skys1d
	; interpolate the 1d sky spectra onto the wavelength range of the
	; corresponding lab spectrum
	for i=0,n_elements(skys1d)-1 do begin
		; read in the 1d median extracted sky and the lab spectra
		temps=mrdfits(skys1d[i],0,hs,/silent)
		templ=mrdfits(labs[i],0,hl,/silent)
		; extract their filenames for later
		sky_name=strsplit(skys[i],'.',/extract)
		lab_name=strsplit(labs[i],'.',/extract)
		; determine the sizes of the arrays
		ss_sky=size(temps)
		ss_lab=size(templ)

		; create the abscissa that we're interpolating FROM
		x = dindgen(ss_sky[1])*sxpar(hs,'CD1_1')+sxpar(hs,'CRVAL1')
		; create the abscissa that we're interpolating ONTO
		u = dindgen(ss_lab[1])*sxpar(hl,'CD1_1')+sxpar(hl,'CRVAL1')
		; interpolate the sky wlens onto the lab wlens
		result = interpol(temps,x,u,/spline)
		; write the interpolated sky spectrum to file
		writefits,sky_name[0]+'_interp.fits',result,hs

		; test to see how well it worked!
;		plot,x,temps,xr=[4900,5100],yr=[-.01,.05]	; this is the sky
;		;oplot,u,result						; results of interpolation
;		spawn,'ls '+labs[i]+' > shit.list'
;		extract_spectra,'shit.list'
;		labname = strsplit(labs[i],'.',/extract)
;		test = mrdfits(labname[0]+'_1d.fits',0,h)		
;		sub = test - result
;		oplot, u, test, color = cgcolor('blue')
;		oplot, u, sub, color=cgcolor('red')			

		; create a 2D median, interpolated sky spectra	
		med_sky = fltarr(ss_lab[1],ss_lab[2])	
		for j=0,ss_lab[1]-1 do med_sky[j,*] = result[j]
		; save the 2D median sky
		writefits,sky_name[0]+'_med.fits',med_sky,hs
		; now subtract this 2d interpolated median sky spectrum from the
		; appropriate lab spectrum
		sky_sub = templ - med_sky
		writefits,lab_name[0]+'_ss.fits',sky_sub,hl

		final = fltarr(ss_lab[1],ss_lab[2])
		; mask out the lyman alpha line
		mask = where((u gt 4700 and u lt 4900) or (u gt 5100 and u lt 5500))
		xx = dindgen(ss_lab[1])
		y = fltarr(ss_lab[1])
		; subtract a linear from each row of the sky_subtracted spectrum 
		for j=0,ss_lab[2]-1 do begin
			skyfit=sky_sub[mask,j]
			rfit = poly_fit(xx[mask],skyfit,1,STATUS=stat)
			y = rfit[0]+rfit[1]*xx
			final[*,j] = sky_sub[*,j] - y
		endfor
		writefits,lab_name[0]+'_ssfit2.fits',final,hl
		; Create 1D spectra
		lab1d=fltarr(ss_lab[1])
		for j=0,ss_lab[1]-1 do lab1d[j]=total(final[j,4:73])
		plot, u, lab1d, xr=[4800,5500]
		name=strsplit(labs[i],'/.',/extract)
		writefits,name[0]+'/1dspectra/'+name[2]+'_1d_ssfit2.fits',lab1d,hl
;stop
;		writefits,lab_name[0]+'_skyfitsub.fits',flat,h
		; test to see how well sky-subtracted the lab spectra are by
		; extracting and plotting 1d spec
;		temp = fltarr(ss_lab[1])
;		for j=0,ss_lab[1]-1 do temp[j]=mean(flat[j,*])
;		plot,u,temp,xr=[4700,5400]
	endfor
end

pro build_stack_cl, suffix, date, TWOD=twod
	if KEYWORD_SET(twod) then dir = '2dspectra' else dir = '1dspectra'
	blah = strsplit(dir,'s',/extract)
	angle = ['00','22','45','67']
	openw,out,'stack_'+blah[0]+'.cl',/get_lun
	printf,out,'# this is 8 groups of 10 images each to create the ord and ext at 4 angles'
	for i=0,3 do begin
		printf,out,'!ls N*/'+dir[0]+'/N*slit1*'+angle[i]+'*'+suffix+' > '+ $
				   'combine_temp.list'
		printf,out,'imcombine @combine_temp.list stacked/A'+angle[i]+'_stack_'+ $
				   date+'_'+blah[0]+'.fits combine=average offsets=wcs '+ $
				   'reject=sigclip nlow=1 nhigh=1 sigmas=stacked/A'+angle[i]+ $
				   '_stack_'+date+'_'+blah[0]+'_sig.fits'
		printf,out,''
	endfor
	for i=0,3 do begin
		printf,out,'!ls N*/'+dir[0]+'/N*slit2*'+angle[i]+'*'+suffix+' > '+ $
				   'combine_temp.list'
		printf,out,'imcombine @combine_temp.list stacked/B'+angle[i]+'_stack_'+ $
				   date+'_'+blah[0]+'.fits combine=average offsets=wcs '+ $
				   'reject=sigclip nlow=1 nhigh=1 sigmas=stacked/B'+angle[i]+ $
				   '_stack_'+date+'_'+blah[0]+'_sig.fits'
		printf,out,''

	endfor
	printf,out,''
	printf,out,'# this is 4 groups of 20 images to create individual intensity images at each angle'
	for i=0,3 do begin
		printf,out,'!ls N*/'+dir[0]+'/N*slit[1,2]*'+angle[i]+'*'+suffix+' > '+ $
				   'combine_temp.list'
		printf,out,'imcombine @combine_temp.list stacked/I'+angle[i]+'_stack_'+date+ $
				   '_'+blah[0]+'.fits combine=average offsets=wcs reject=sigclip '+ $
				   'nlow=1 nhigh=1 sigmas=stacked/I'+angle[i]+'_stack_'+date+'_'+ $
				   blah[0]+'_sig.fits'
		printf,out,''
	endfor
	printf,out,''
	printf,out,'# this is 1 group of 80 images '
	printf,out,'!ls N*/'+dir[0]+'/N*slit[1,2]*'+suffix+' > combine_temp.list'
	printf,out,'imcombine @combine_temp.list stacked/Iall_stack_'+date+'_'+ $
			   blah[0]+'.fits combine=average offsets=wcs reject=sigclip '+ $
			   'nlow=1 nhigh=1 sigmas=stacked/Iall_stack_'+date+'_'+blah[0]+'_sig.fits'
	close,out
	free_lun,out
end

pro fix_cols_A, filelist, TWOD=twod
	readcol, filelist, format='a',files
	if KEYWORD_SET(twod) then begin
		openw,out,'adjust_pixels2d.cl',/get_lun
		for i=0,n_elements(files)-1 do begin
			name=strsplit(files[i],'.',/extract)
			printf,out,'imcopy '+files[i]+'[2:2069,*] '+files[i]+'[0,overwrite]'
			printf,out,'imcopy '+name[0]+'_sig.fits[2:2069,*] '+ $
				   name[0]+'_sig.fits[0,overwrite]'
		endfor
	endif else begin
		openw,out,'adjust_pixels1d.cl',/get_lun
		for i=0,n_elements(files)-1 do begin
			name=strsplit(files[i],'.',/extract)
			printf,out,'imcopy '+files[i]+'[2:2069] '+files[i]+'[0,overwrite]'
			printf,out,'imcopy '+name[0]+'_sig.fits[2:2069] '+ $
				   name[0]+'_sig.fits[0,overwrite]'
		endfor		
	endelse
	close,out
	free_lun,out
end
pro fix_cols_I, filelist, TWOD=twod
	readcol, filelist, format='a',files
	if KEYWORD_SET(twod) then begin
		openw,out,'adjust_pixels2d.cl',/get_lun	
		for i=0,n_elements(files)-1 do begin
			name=strsplit(files[i],'.',/extract)
			printf,out,'imcopy '+files[i]+'[3:2070,*] '+files[i]+'[0,overwrite]'
			printf,out,'imcopy '+name[0]+'_sig.fits[3:2070,*] '+ $
				   name[0]+'_sig.fits[0,overwrite]'
		endfor
	endif else begin
		openw,out,'adjust_pixels1d.cl',/get_lun	
		for i=0,n_elements(files)-1 do begin
			name=strsplit(files[i],'.',/extract)
			printf,out,'imcopy '+files[i]+'[3:2070] '+files[i]+'[0,overwrite]'
			printf,out,'imcopy '+name[0]+'_sig.fits[3:2070] '+ $
				   name[0]+'_sig.fits[0,overwrite]'
		endfor
	endelse
	close,out
	free_lun,out
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;  BEGINNING OF OFFICIAL REDUCTION   ;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro reduction
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	spawn,'pwd', wdir
	
	date = '7814'
	; specify which directories/data to reduce
	dirs = 'N*'	; this specifies the science data; 
				; starN* specifies standard star dirs
	to_reduce = dirs+'/N*.fits'
	sci_bias =  dirs+'/calibration/sci_bias/*BIAS*.fits'
	flat_bias = dirs+'/calibration/flat_bias/*BIAS*.fits'
	dirlist = 'directories.list'
	
	; y pixel positions of the top and bottom of slits
	lab_slit=[[263,340],[357,434]] 
	sky_slit=[[445,522],[533,610]]	
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; Initial Reduction: Bias Subtraction, Flatfielding, Normalizing, 
	; Spatial Alignment
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	; creates ROTANG keyword in the header (rotation angle of HWP)
	; renames files to this format: NX_X/NX_X_HWP.fits
;	spawn,'ls '+to_reduce+' > sci_data.lis'
;	prep_scidata,'sci_data.lis'
	
	; create the master bias frames (both for science imgs and flats)
;	spawn,'ls '+sci_bias+' > sci_bias.lis
;	master_bias, 'sci_bias.lis'
;	spawn,'ls '+flat_bias+' > flat_bias.lis
;	master_bias, 'flat_bias.lis'
	
	; subtract the master bias frames from the appropriate images
;	bias_subtract,dirlist,/SCIENCE
;	bias_subtract,dirlist,/FLATS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	print, 'BIAS SUBTRACTION COMPLETE'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;stop
	; Use the master flats created by the FORS pipeline
;	flatfield, dirlist
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	print, 'FLAT FIELDING COMPLETE'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;stop
	; remove cosmic rays with rmCosmic.py
;	spawn,'rmCosmic.py'
;	spawn,'ls '+dirs+'/*ff.fits > sci_data.lis'
;	spawn,'rmCosmic	-l sci_data.lis -gain 1.89 -rn 3.2 -maxiter 4 -sigclip 8'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	print, 'COSMIC RAYS REMOVED'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; the cosmic ray remover doesn't save the headers properly and I can't 
	; figure out how to change the python code appropriately so ... a quick
	; procedure to fix them:
;	spawn,'ls '+dirs+'/*crc.fits > fixheaders.list'
;	spawn,'ls '+dirs+'/*ff.fits > getheaders.list'
;	fix_crc_headers, 'fixheaders.list', 'getheaders.list'
	
	; normalize all sci data by exposure time
;	spawn,'ls '+dirs+'/*crc.fits > sci_data.lis'
;	normalize, 'sci_data.lis'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	print, 'DATA NORMALIZED'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; align the data in the y direction by fitting the stellar continua in slits
	; 2 and 5
;	spawn,'ls '+dirs+'/*_crc_n.fits > sci_data.lis'
	; find the shift in the y-direction for each science image (except for N1_1
	; and N1_2
;	yshifts = find_yshift('sci_data.lis')
;	yalign,'sci_data.lis',yshifts
	
	; cut out main slits 3 and 4 (actual slits 5-8)
	; Slit 3 are sky slits -- Slit 4 are LAB1 slits
;	spawn,'ls '+dirs+'/*n_y.fits > sci_data.lis'
;	cut_slits,'sci_data.lis',[[lab_slit],[sky_slit]]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	print, 'ALIGNED IN SPATIAL DIRECTION'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;stop
	
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Wavelength Calibration
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	; load in comparison lamps
;	spawn,'ls '+dirs+'/calibration/arclamp* > lamps.lis'
	; cut out slits in the comp lamps the same size as the science data
;	cut_slits,'lamps.lis',[[lab_slit],[sky_slit]] 
	
	; Wavelength calibration done in IRAF:
	; 1. epar IDENTIFY (run on the 2D comp lamp slits for each unique dither)
	;		images = lampslit1(2).fits
	;		section = middle line
	;		coordli = ../linelist
	;  		all other parameters are defaults
	;	identification fit with legendre/chebyshev functions
	; 	wavelength solution with RMS ~ 0.06 achieved for all slits keeping at
	; 	least 10 lines each
	
	; 2. epar REIDENTIFY (run on the 2D comp lamp slits for each unique dither
	;			(using itself as a reference)
	;			find a 2D solution for each 2D comp lamp slit 
	;		reference = lampslit[1,2,3,4].fits
	;		images = 	lampslit[1,2,3,4].fits
	;		interac = 	yes	(interactive fitting)
	;		section = 	middle line
	;		trace = yes (use the previous fit for each subsequent fitting)
	;		step = 3	(number of lines to skip between each subsequent
	;				 	fitting)
	
	; IF WORKING WITH 2D SPECTRA:
	; 3. epar FITCOORDS (apply to comp lamp slits that were run in REIDENTIFY --
	;			this creates a new files called database/fclampslit[1,2,3,4]
	;		images = 	lampslit[1,2,3,4]
	;		interac = 	yes (fit interactively)
	
	; 4. epar TRANSFORM (apply the coordinate fit to the science data)
	; 		If steps 1-3 have already been completed, they don't need to be
	;		run a second time -- Only the TRANSFORM step is required
	;		(Don't forget to add DISPAXIS = 1 to header file!!)
	
	; build the list of files to TRANSFORM
;	trans_file = transform_list()
	; Adjust their headers to include the DISPAXIS keyword
;	adjust_headers, trans_file, 'DISPAXIS', 1
	; Build the .cl IRAF script to perform TRANSFORM
;	build_transform_cl, trans_file
	; run this script in an IRAF xgterm
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	print,'Run transform.cl in an IRAF xgterm. Then come back.'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;stop
	; perform SKY SUBTRACTION
;	spawn,'ls '+dirs+'/2dspectra/N*slit[1,2]*wc.fits > 2dsci_data.lis'
;	spawn,'ls '+dirs+'/2dspectra/N*slit[3,4]*wc.fits > 2dsky_data.lis'
;	sky_subtract, '2dsci_data.lis', '2dsky_data.lis'

	; Correct for ATMOSPHERIC EXTINCTION
;	spawn,'ls '+dirs+'/2dspectra/N*slit[1,2]*wc_ssfit2.fits > sci_data.lis'
;	correct_atmos, 'sci_data.lis',/twod
;	spawn,'ls '+dirs+'/1dspectra/N*slit[1,2]*wc_1d_ssfit2.fits > sci_data.lis'
;	correct_atmos, 'sci_data.lis'
;stop
	suffix = '715'
	; STACK the spectra in IRAF
	build_stack_cl, 'wc_ssfit2_acorr3.fits', suffix, /twod
	build_stack_cl, 'wc_1d_ssfit2_acorr3.fits', suffix
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	print,"Run stack.cl in IRAF then return." 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
stop
	; For some reason when IRAF stacks these, the As and Bs have a ~1 pix diff
	; in column length and the Is and Bs have ~2 pixels difference! 
	; Fix them in IRAF so that it will automatically adjust the headers to 
	; preserve the wavelength calibration

	; fix the As - 2D
	spawn,'ls stacked/'+suffix+'/A*'+suffix+'*2d.fits > fix_cols2d.lis
	fix_cols_A, 'fix_cols2d.lis',/twod
	; fix the As - 1D
	spawn,'ls stacked/'+suffix+'/A*'+suffix+'*1d.fits > fix_cols1d.lis
	fix_cols_A, 'fix_cols1d.lis'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	print,"Run adjust_pixels.cl in IRAF (round 1)" 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
stop
	; fix the Is - 2D
	spawn,'ls stacked/'+suffix+'/I*'+suffix+'*2d.fits > fix_cols2d.lis
	fix_cols_I, 'fix_cols2d.lis',/twod
	; fix the Is - 1D
	spawn,'ls stacked/'+suffix+'/I*'+suffix+'*1d.fits > fix_cols1d.lis
	fix_cols_I, 'fix_cols1d.lis'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	print,"Run adjust_pixels.cl in IRAF (round 2)" 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; READY TO GO SCIENCE PRODUCT HURRAY!
	; MOVE ON TO POLARIMETRY!
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	print,"REDUCTION COMPLETE." 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	stop
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;   END OF OFFICIAL REDUCTION   ;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;spawn,'ls N*/1dspectra/N*slit[1,2]*1d_calib.fits > 1dlab_spec.list'
;spawn,'ls N*/1dspectra/N*slit[3,4]*1d_calib.fits > 1dsky_spec.list'
;sky_subtract, '1dsky_spec.list', '1dlab_spec.list'

; create smaller samples of the spectra averaged over a given pixel size
;spawn,'ls N*/1dspectra/*skysub.fits > 1d_reduced_skysub.list'
;rebin, '1d_reduced_skysub.list'

	
;		; fit residual background with a linear polynomial
;		; don't fit to the lyman alpha line itself -- create a mask
;		x_fit = fltarr(ss_lab[1],ss_lab[2])
;		y_fit = fltarr(ss_lab[1],ss_lab[2])
;		for k=0,n_elements(ss_lab[2]) do x_fit[*,k]=x
;		for k=0,n_elements(ss_lab[1]) do y_fit[k,*]=dindgen(ss_lab[1])


