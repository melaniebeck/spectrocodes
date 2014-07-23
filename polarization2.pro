; calculate polarization the ESO/Claudia/Matt way
function polcalc1, flux
	q = 0.5*(flux[0]-flux[4])/(flux[0]+flux[4])	- $
		0.5*(flux[2]-flux[6])/(flux[2]+flux[6])
	u = 0.5*(flux[1]-flux[5])/(flux[1]+flux[5]) - $
		0.5*(flux[3]-flux[7])/(flux[3]+flux[7])
	return, sqrt(q^2+u^2)
end

function polcalc2, a1,a2,a3,a4,b1,b2,b3,b4
	f1 = (a1-b1)/(a1+b1)
	f2 = (a2-b2)/(a2+b2)
	f3 = (a3-b3)/(a3+b3)
	f4 = (a4-b4)/(a4+b4)
	q = 0.5*f1 - 0.5*f3
	u = 0.5*f2 - 0.5*f4
	return, sqrt(q^2 + u^2)
end

; NAME: 	mc_sim
; PURPOSE:	runs a monte carlo simulation on the measured polarization and 
;			uses the statistics of the resulting histogram to compute the 
;			polarization error
; INPUTS:	FLUX VALUES for each of the HWP positions and both halves of the
;			Wolly (8 flux values in total),
;			SIGMA VALUES for each of the 8 flux values,
;			Original POLARIZATION measurement
; OUTPUTS: 	Various error calculations (1sig, 2sig, percentiles?)
function mc_sim, hwps, sigs, pol, coord, date, PLOT=PLOT
	n = 10000
	sims = fltarr(n_elements(hwps))
	sim_p = fltarr(n)
	ff=dindgen(10)*1000.
	for j=1,n do begin
		; simulate a NEW value of each of the ord/ext beams by adding to that
		; value a small deviation based on a random fluctuation of the sigma map
		sims[0] = hwps[0] + (randomu(seed)-0.5d)*sigs[0]*2
		sims[1] = hwps[1] + (randomu(seed)-0.5d)*sigs[1]*2
		sims[2] = hwps[2] + (randomu(seed)-0.5d)*sigs[2]*2
		sims[3] = hwps[3] + (randomu(seed)-0.5d)*sigs[3]*2
		sims[4] = hwps[4] + (randomu(seed)-0.5d)*sigs[4]*2
		sims[5] = hwps[5] + (randomu(seed)-0.5d)*sigs[5]*2
		sims[6] = hwps[6] + (randomu(seed)-0.5d)*sigs[6]*2
		sims[7] = hwps[7] + (randomu(seed)-0.5d)*sigs[7]*2
		; calculate the simulated polarization
		sim_p(j-1) = polcalc1(sims)
		fr=where(ff eq j)
       ; if fr[0] gt 0. then print,'done ',ff[fr]/10000.,'%'
	endfor
	; calculate statistics typical of a symmetric distribution
	avg = mean(sim_p)
	std = stddev(sim_p)
	; calculate statistics for a skewed distribution
	med = median(sim_p)
	percent = PERCENTILES(sim_p,VALUE=[0.025,0.16,0.5,0.84,0.975])
	; define the error as the larger of the two absolute deviations from
	; the measurement 
	sigma1 = (pol - percent[1]) > (percent[3] - pol)
	sigma2 = (pol - percent[0]) > (percent[4] - pol)
	; If requested, plot the histogram and various statistics
	if KEYWORD_SET(PLOT) then begin
		set_plot,'PS'
		if n_elements(coord) eq 2 then begin
			device,filename='plots/'+date+'/hist_'+date+'_'+strtrim(coord[0],2)+'_'+ $
						strtrim(coord[1],2)+'.ps',/color,/encap
			title = 'Pol Histrogram for bin '+strtrim(coord[0],2)+','+ $
				 		strtrim(coord[1],2)
		endif else begin
			device,filename='plots/'+date+'/hist_'+date+'_'+strtrim(coord[0],2)+'.ps'
			title = 'Pol Histrogram for bin '+strtrim(coord[0],2)		
		endelse		
;		device,filename=name,/color,/encap
		; plot the histogram of simulated polarizations
		plothist,sim_p,xhist,yhist,bin=0.01,/fill,fcolor=230,xr=[0,1.], $
				 xtit='Polarization fraction'
		; overplot the mean and median of the simulated polarization dist
		oplot,yhist*0.+avg,yhist,thick=2,color=cgcolor('red')
		oplot,yhist*0.+med,yhist,thick=2,color=cgcolor('green')
		; overplot the original measured polarization
		oplot,yhist*0.+pol,yhist,thick=2,color=cgcolor('black') 
		; overplot the avg sim pol +- std of the histogram
		oplot,yhist*0.+(avg+std),yhist,thick=2,linestyle=1,color=cgcolor('blue')
		oplot,yhist*0.+(avg-std),yhist,thick=2,linestyle=1,color=cgcolor('blue')
		; overplot the avg sim pol +- 2*std of the histogram
		oplot,yhist*0.+(avg+2*std),yhist,thick=2,linestyle=1, $
			  color=cgcolor('magenta')
		oplot,yhist*0.+(avg-2*std),yhist,thick=2,linestyle=1, $
			  color=cgcolor('magenta')
		; overplot the percentiles of the histogram
;		oplot,yhist*0.+percent[2],yhist,thick=2,linestyle=3, $
;			  color=cgcolor('orange')
    	oplot,yhist*0.+percent[1],yhist,thick=2,linestyle=2, $
			  color=cgcolor('blue')
   		oplot,yhist*0.+percent[3],yhist,thick=2,linestyle=2, $
			  color=cgcolor('blue')
 		oplot,yhist*0.+percent[0],yhist,thick=2,linestyle=2, $
			  color=cgcolor('magenta')
 		oplot,yhist*0.+percent[4],yhist,thick=2,linestyle=2, $
			  color=cgcolor('magenta')

		items = ['Raw P', 'Mean', 'Median', 'Geom (68%)', 'Geom (95%)', $
				 'PDF (68%)', 'PDF (95%)']
		lines = [0,0,0,1,1,2,2]
		colors = [cgcolor('black'), cgcolor('red'), cgcolor('green'), cgcolor('blue'), $
				 cgcolor('magenta'), cgcolor('blue'), cgcolor('magenta')]
		al_legend,items, linestyle=lines,color=colors,linsize=.5,/right
		device,/close
		set_plot,'X
;stop
	endif
	; return SIGMA1 and SIGMA2 
	return, [sigma1, sigma2]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro pol2d, date, TWOD=twod
	;if KEYWORD_SET(twod) then date = date+'_2d' else date=date+'_1d'
	; read in ord/ext beams
	spawn,'ls stacked/'+date+'/[A,B]*_stack_*'+date+'_2d.fits > stacked_sci.list'
;	spawn,'ls stacked/'+date+'/noacorr/[A,B]*_stack_noacorr.fits > stacked_sci.list'
	readcol,'stacked_sci.list',format='a',sci_files
	; read in the matching sigma maps
	spawn,'ls stacked/'+date+'/[A,B]*_stack_*'+date+'_2d_sig.fits > stacked_sig.list'
;	spawn,'ls stacked/'+date+'/noacorr/[A,B]*_stack_noacorr_sig.fits > stacked_sig.list'
	readcol,'stacked_sig.list',format='a',sig_files

	; choose bin sizes
	bx = 8.
	by = 10.
	start_wlen = 4930.
	end_wlen = 5032.
	; read in a test file to determine wavelength range(s)
	temp=mrdfits(sci_files[0],0,h)
	ss = size(temp)
	wlen = dindgen(ss[1])*sxpar(h,'CD1_1')+sxpar(h,'CRVAL1')- $
		   sxpar(h,'CRPIX1')*sxpar(h,'CD1_1')
	y = dindgen(ss[2])
	; x range of 2d spectrum to look at
	xrange=where(wlen gt start_wlen and wlen lt end_wlen)
;stop
	lx = xrange[0]
	; y range of 2d spectrum to look at
	yrange=where(y gt 3 and y lt 74)
	ly = yrange[0]
	; number of total bins in x and y directions
	num_x = (n_elements(xrange))/bx		
	num_y = (n_elements(yrange))/by
stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; declare variables for science spectra binning
	spec=fltarr(ss[1],ss[2],n_elements(sci_files))
	spec1d_t=fltarr(num_x*bx,n_elements(sci_files))	
	spec1d_b=fltarr(num_x*bx,n_elements(sci_files))

	spec_c=fltarr(num_x*bx,num_y*by,n_elements(sci_files))
	spec_b=fltarr(num_x,num_y,n_elements(sci_files))
	; declare variables for sigma binning	
	sigs=fltarr(ss[1],ss[2],n_elements(sci_files))	
	sigs1d_t=fltarr(num_x*bx,n_elements(sci_files))	
	sigs1d_b=fltarr(num_x*bx,n_elements(sci_files))	

	sigs_c=fltarr(num_x*bx,num_y*by,n_elements(sci_files))
	sigs_b=fltarr(num_x,num_y,n_elements(sci_files))
	sigs_b2=sigs_b
	; declare variables for polarization
	pol = fltarr(num_x,num_y)
	perr = fltarr(num_x,num_y,2)	
	pol1d = fltarr(num_x)
	perr1d = fltarr(num_x,2)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; CROP & REBIN 2D TOTAL INTENSITY SPECTRUM 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	tot_a = mrdfits('stacked/'+date+'/Iall_stack_'+date+'_2d.fits',0,hta)
	ss = size(tot_a)
	tot1d_a = fltarr(ss[1])
	tot_c = fltarr(num_x*bx,num_y*by)
	tot_d = fltarr(num_x,num_y)

	; wavelength range of 2d spectrum
	wlen1 = dindgen(ss[1])*sxpar(hta,'CD1_1')+sxpar(hta,'CRVAL1')- $
		   	   sxpar(hta,'CRPIX1')*sxpar(hta,'CD1_1')
	; x range of 2d spectrum to look at
	xrange1=where(wlen1 gt start_wlen and wlen1 lt end_wlen)

	; extract 1d spectrum for use later
	for i=0,ss[1]-1 do tot1d_a[i]=total(tot_a[i,4:73])
	writefits,'stacked/'+date+'/Iall_stack_'+date+'_1d_2.fits',tot1d_a,hta
;stop
	; test wlen range with some plots
	plot,wlen1,tot1d_a,xr=[4800,5200]
	oplot,dindgen(100)*0. + 4981.5,dindgen(100),thick=2;,linestyle=2,color=1

	; crop the 2d total intensity spectrum
	tot_c = tot_a[xrange1,4:73]
	; rebin the 2d tot int spectrum
	tot_b = rebin(tot_a[xrange1,4:73],num_x,num_y)

	writefits,'stacked/'+date+'/crop/Iall_stack_'+date+'_2d_crop_2.fits',tot_c,hta
	writefits,'stacked/'+date+'/rebin/Iall_stack_'+date+'_2d_rebin_2.fits',tot_b,hta
;stop
;	writefits,'stacked/'+date+'/noacorr/crop/Iall_stack_noacorr_2d_crop.fits',tot_c,ht
;	writefits,'stacked/'+date+'/noacorr/rebin/Iall_stack_noacorr_2d_rebin.fits',tot_b,ht

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; CREATE 1D TOTAL INTENSITY SPECTRA -- TOP AND BOTTOM
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	ssc=size(tot_c)
	tot_bot = fltarr(ssc[1])
	tot_top = tot_bot

	; extract 2 1d spectra
	for i=0,ssc[1]-1 do begin
		; 1d tot intensity for the bottom bright line
		tot_bot[i] = total(tot_c[i,11:28])
		; 1d tot intensity for the top bright line
		tot_top[i] = total(tot_c[i,35:54])
	endfor

	; test plot some shits
	crop = where(wlen1 gt start_wlen and wlen1 lt end_wlen)
	plot,wlen1[crop],tot1d_a[crop];,xr=[4930,5057.5];,yr=[-0.1,.5]
	oplot,wlen1[crop],tot_bot,color=cgcolor('blue')
	oplot,wlen1[crop],tot_top,color=cgcolor('red')
	oplot,dindgen(100)*0. + 4981.5,dindgen(100),thick=2;,linestyle=2,color=1

	; save data to files
	writefits,'stacked/'+date+'/crop/Iall_stack_'+date+'_1d_top_crop_2.fits', tot_top, hta
	writefits,'stacked/'+date+'/crop/Iall_stack_'+date+'_1d_bot_crop_2.fits', tot_bot, hta
;stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; CROP AND REBIN ALL SCIENCE SPECTRA
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	for i=0,n_elements(sci_files)-1 do begin
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		; READ IN SCIENCE SPECTRA
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		spec[*,*,i] = mrdfits(sci_files[i],0,h)	
		; the wlen could be slightly different because each header could have
		; a different value of CRPIX1
		wlen2 = dindgen(ss[1])*sxpar(h,'CD1_1')+sxpar(h,'CRVAL1')- $
		   sxpar(h,'CRPIX1')*sxpar(h,'CD1_1')	
		xrange2=where(wlen2 gt start_wlen and wlen2 lt end_wlen) ;5057.3
		name = strsplit(sci_files[i],'/.',/extract)
;stop
		; CROP 2D SCIENCE SPECTRA
		spec_c[*,*,i] = spec[xrange2,yrange,i]
		writefits,name[0]+'/'+name[1]+'/crop/'+name[2]+'_crop_2.fits',spec_c[*,*,i],h
;		writefits,name[0]+'/'+date+'/noacorr/crop/'+name[3]+'_2d_crop.fits',spec_c[*,*,i],h
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		; CREATE 1D SCIENCE PRODUCTS -- TOP AND BOTTOM SPECTRA
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		; extract the science data 
		name1d=strsplit(sci_files[i],'_',/extract)
		for j=0,ssc[1]-1 do spec1d_t[j,i] = total(spec_c[j,11:28,i])
		for j=0,ssc[1]-1 do spec1d_b[j,i] = total(spec_c[j,34:54,i])
		writefits,name[0]+'/'+date+'/crop/'+name1d[2]+'_stack_'+ $
				  date+'_1d_top_crop_2.fits',spec1d_t,h
		writefits,name[0]+'/'+date+'/crop/'+name1d[2]+'_stack_'+ $
				  date+'_1d_bot_crop_2.fits',spec1d_b,h
		; rebin the 1D science spectra
;		spec1d_b[*,i] = rebin(spec[xrange,i],num_x)
;		writefits,name[0]+'/'+name[1]+'/rebin/'+name[2]+'_1d_rebin.fits',spec1d_b,h

		; REBIN 2D SCIENCE SPECTRA
		spec_b[*,*,i] = rebin(spec[xrange2,yrange,i],num_x,num_y)
		writefits,name[0]+'/'+date+'/rebin/'+name[2]+'_rebin_2.fits',spec_b[*,*,i],h
;		writefits,name[0]+'/'+date+'/noacorr/rebin/'+name[3]+'_2d_rebin.fits',spec_b[*,*,i],h
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		; READ IN SIGMA MAPS
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		sigs[*,*,i] = mrdfits(sig_files[i],0,h)
		name = strsplit(sig_files[i],'/.',/extract)
		; CROP 2D SIGMA MAPS TO SAME RANGE AS SCIENCE FRAMES
		sigs_c[*,*,i] = sigs[xrange2,yrange,i]
		writefits,name[0]+'/'+date+'/crop/'+name[2]+'_crop_2.fits',sigs_c[*,*,i],h
;		writefits,name[0]+'/'+date+'/noacorr/crop/'+name[3]+'_2d_crop.fits',sigs_c[*,*,i],h
		; REBIN 2D SIGMA MAPS
		sigs_b[*,*,i] = sqrt(rebin(sigs_c[*,*,i]^2,num_x,num_y))/sqrt(bx*by)
		writefits,name[0]+'/'+date+'/rebin/'+name[2]+'_rebin_2.fits',sigs_b[*,*,i],h
;		writefits,name[0]+'/'+date+'/noacorr/rebin/'+name[3]+'_2d_rebin.fits',sigs_b[*,*,i],h
		for j=0,num_x-1 do begin
			for k=0,num_y-1 do begin
				sigs_b2[j,k,i] = sqrt(total(sigs_c[0+j*bx:0+(j+1)*bx-1, $
										0+k*by:0+(k+1)*by-1,i]^2))
			endfor
		endfor
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		; CREATE 1D SIGMA PRODUCT 
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;sigs1d_b2=sigs1d_b
;		sigs1d_b[*,i] = sqrt(rebin(sigs_c[*,*,i]^2,num_x))/ $
;						sqrt(bx*n_elements(yrange))
;		writefits,name[0]+'/'+date+'/rebin/'+name[1]+'_1drebin.fits',sigs1d_b,h
;		for j=0,num_x-1 do begin
;			sigs1d_b2[j,i] = sqrt(total(sigs_c[0+j*bx:0+(j+1)*bx-1,*,i]^2))
;		endfor

	endfor	
;stop
;	openw,out,'MC_error68_'+date,/get_lun
	; Calculate the measured POLARIZATION
	for i=0,num_x-1 do begin
		; calculate the polarization for the 1D data
;		pol1d[i] = polcalc1(spec1d_b[i,*])
		;if (i eq 10) or (i eq 11) or (i eq 12) then 
;		perr1d[i,*]= $
;		   mc_sim(spec1d_b[i,*],sigs1d_b[i,*],pol1d[i],i,date,/plot) 
		;else perr1d[i,*]=mc_sim(spec1d_b[i,*],sigs1d_b[i,*],pol1d[i])
		for j=0,num_y-1 do begin
			; calculate the polarization for the 2D data
			pol[i,j] = polcalc1(spec_b[i,j,*])
			if (i eq 10 and j eq 5) or (i eq 11 and (j eq 4 or j eq 2)) or $
			   (i eq 12 and (j eq 1 or j eq 3)) or (i eq 5 and j eq 4) then perr[i,j,*] =  $
			   mc_sim(spec_b[i,j,*],sigs_b[i,j,*],pol[i,j],[i,j],date,/plot) $
			else perr[i,j,*] = mc_sim(spec_b[i,j,*], sigs_b[i,j,*], pol[i,j])
;			printf,out,perr[i,j,0]
		endfor
	endfor
;	close,out
;	free_lun,out
;stop
	; WRITE 1D DATA TO FILE
;	writefits,'pol_products/MC_polstn68_'+date+'_1d.fits',pol1d/perr1d[*,0]
;	writefits,'pol_products/MC_polstn95_'+date+'_1d.fits',pol1d/perr1d[*,1]
;	writefits,'pol_products/MC_sigma68_'+date+'_1d.fits',perr1d[*,0]
;	writefits,'pol_products/MC_sigma95_'+date+'_1d.fits',perr1d[*,1]
;	writefits,'pol_products/pol_'+date+'_1d.fits',pol1d
	; WRITE 2D DATA TO FILE
	writefits,'pol_products/MC_polstn68_'+date+'_2.fits',pol/perr[*,*,0]
	writefits,'pol_products/MC_polstn95_'+date+'_2.fits',pol/perr[*,*,1]
	writefits,'pol_products/MC_sigma68_'+date+'_2.fits',perr[*,*,0]
	writefits,'pol_products/MC_sigma95_'+date+'_2.fits',perr[*,*,1]
	writefits,'pol_products/pol_'+date+'_2.fits',pol
stop

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Must already have run pol2d FIRST so that you have certain files required for 1D
pro pol1d, date
	; read in ORD/EXT Beams (ORIGINALS)
	spawn,'ls stacked/[A,B]*'+date+'.fits > stacked_sci.list'
	readcol,'stacked_sci.list',format='a',sci_files
	; read in the matching sigma maps
	spawn,'ls stacked/[A,B]*'+date+'_sig.fits > stacked_sig.list'
	readcol,'stacked_sig.list',format='a',sig_files

	; choose bin size
	bx = 8.

	; read in a test file to determine wavelength range(s)
	temp=mrdfits(sci_files[0],0,h)
	ss = size(temp)
	wlen = dindgen(ss[1])*sxpar(h,'CD1_1')+sxpar(h,'CRVAL1')
	; x range of 1d spectrum to look at
	xrange=where(wlen gt 4930. and wlen lt 5057.5)
	lx = xrange[0]
	; number of total bins in x direction
	num_x = (n_elements(xrange))/bx		

	; declare variables for science spectra binning
	spec=fltarr(ss[1],n_elements(sci_files))		
	spec_c=fltarr(num_x*bx,n_elements(sci_files))
	spec_b=fltarr(num_x,n_elements(sci_files))
	; declare variables for sigma binning	
	sigs=fltarr(ss[1],n_elements(sci_files))		
	sigs_c=fltarr(num_x*bx,n_elements(sci_files))
	sigs_b=fltarr(num_x,n_elements(sci_files))
	pol = fltarr(num_x)
	perr = fltarr(num_x,2)

	; read in and rebin all the science spectra and sigma maps
	for i=0,n_elements(sci_files)-1 do begin
		; read in the science spectra
		spec[*,i] = mrdfits(sci_files[i],0,h)
		name = strsplit(sci_files[i],'/.',/extract)
		; extract the 2D spectra
		
		; rebin the 1d science spectra
		spec_b[*,i] = rebin(spec[xrange,i],num_x)
		writefits,name[0]+'/rebin/'+name[1]+'_1drebin.fits',spec_b[*,i],h
		; read in the sigma maps
		sigs[*,i] = mrdfits(sig_files[i],0,h)
		name = strsplit(sig_files[i],'.',/extract)
		; rebin the sigma maps
		sigs_b[*,i] = sqrt(rebin(sigs_c[*,i]^2,num_x))/sqrt(bx)
		writefits,name[0]+'/rebin/'+name[1]+'_1drebin.fits',sigs_b[*,i],h
	endfor
	openw,out,'MC_error68_'+date+'_1d',/get_lun
	; Calculate the measured POLARIZATION
	for i=0,num_x-1 do begin
		pol[i] = polcalc1(spec_b[i,*])
		if (i eq 10) or (i eq 11) or (i eq 12) then perr[i,*] =  $
		   mc_sim(spec_b[i,*],sigs_b[i,*],pol[i],i,date,/plot) $
		else perr[i,*] = mc_sim(spec_b[i,*], sigs_b[i,*], pol[i])
		printf,out,perr[i,0]
	endfor
	close,out
	free_lun,out

	print,pol
stop
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro pol2d_special, date
	; read in cropped 2d ord/ext beams
	spawn,'ls stacked/[A,B]*'+date+'_2dcrop.fits > stacked_sci.list'
	readcol,'stacked_sci.list',format='a',sci_files
	; read in the corresponding sigma maps
	spawn,'ls stacked/[A,B]*'+date+'_sig_2dcrop.fits > stacked_sig.list'
	readcol,'stacked_sig.list',format='a',sig_files
	; choose specific regions of interest
	c1 = [[93,50],[93,40],[103,45],[110,46],[83,23],[91,20],[99,20],[107,20]]
	s1 = [[12,8],[12,12],[8,10],[6,8],[8,10],[8,16],[8,16],[8,16]]	

	ss = size(c1)
	dbins = fltarr(ss[2],n_elements(sci_files))
	sbins = fltarr(ss[2],n_elements(sci_files))
	pol = fltarr(ss[2])
	perr = fltarr(ss[2],2)
	; rebin the data according to these bin sizes
	for i=0,n_elements(sci_files)-1 do begin
		spec = mrdfits(sci_files[i],0,h)
		sig = mrdfits(sig_files[i],0)
		for j=0,ss[2]-1 do begin
			; rebin the data according to these bin sizes
			dbins[j,i] = mean(spec[c1[0,i]-s1[0,i]/2:c1[0,i]+s1[0,i]/2-1, $
									c1[1,i]-s1[1,i]/2:c1[1,i]+s1[1,i]/2-1])
			; rebin the sigma maps according to these bin sizes
			sbins[i,j] = mean(sig[c1[0,i]-s1[0,i]/2:c1[0,i]+s1[0,i]/2-1, $
									c1[1,i]-s1[1,i]/2:c1[1,i]+s1[1,i]/2-1])
		endfor
	endfor
	; calculate the polarization and error in each bin
	for i=0,ss[2]-1 do begin
		pol[i] = polcalc1(dbins[i,*])
		perr[i,*] = mc_sim(dbins[i,*], sbins[i,*], pol[i])
	endfor
	; create stn maps
	stn68 = pol/perr[*,0]
	stn95 = pol/perr[*,1]
	; create a polarization "image" from scratch
	pol_img = spec[*,*]*0.
	for i=0,ss[2]-1 do begin
		if stn68[i] gt 1. then begin
			pol_img[c1[0,i]-s1[0,i]/2:c1[0,i]+s1[0,i]/2-1, $
						c1[1,i]-s1[1,i]/2:c1[1,i]+s1[1,i]/2-1] = pol[i]
		endif
	endfor
	; write data products to file
	sxaddpar, h, 'NOTE','T',' Pol w S/N(sig68) gt 1'
	writefits,'p_special_stn68_'+date+'.fits',pol_img,h

	;stop
	
	; print region files
	openw,out,'pol_special.reg',/get_lun
	printf,out, 'global color=red dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 	highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
	printf,out,'image'
	for i=0,ss[2]-1 do begin
		printf,out,'box('+strtrim(c1[0,i],2)+','+strtrim(c1[1,i],2)+',' $
						 +strtrim(s1[0,i],2)+','+strtrim(s1[0,i],2)+',0)'
		printf,out,'# text('+strtrim(c1[i],2)+','+strtrim(c1[i],2)+') text={'+strtrim(i+1,2)+'}'
	endfor
	close,out
	free_lun,out
	stop	
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

pro comp_errors
	readcol, 'MC_error68_6514', mc_old
	readcol, 'MC_error68_7814', mc_new
	readcol, 'Prop_error_6514', prop_old
	readcol, 'Prop_error_7814', prop_new

	x = dindgen(10)/10.
	set_plot,'PS'
	device,filename='comperr_mcprop.ps',/encap,/color
	plot, x, x, xr=[0,.3], yr=[0,.3],xtit='MC Error (sig68)',ytit='Propagated Error', $
		 tit='MC vs. Propagated error for old and new (acorr) data sets'

	oplot, mc_old, prop_old, psym=2, color=cgcolor('blue')
	oplot, mc_new, prop_new, psym=2, color=cgcolor('red')

	result1 = poly_fit(mc_old[where(mc_old lt 0.3)], prop_old[where(mc_old lt 0.3)],1)
	y1 = result1[0]+result1[1]*x
	oplot, x, y1,linestyle=2
	result2 = poly_fit(mc_new[where(mc_new lt 0.3)], prop_new[where(mc_new lt 0.3)],1)
	y2 = result2[0]+result2[1]*x
	oplot, x, y2,linestyle=3

	xyouts,.03,.27,'6514 Data', color=cgcolor('blue')
	xyouts,.03,.25,'7814 Data (acorr)', color=cgcolor('red')
	device,/close

	device,filename='comperr_oldnew.ps',/encap,/color
	plot, x, x, xr=[0,.15], yr=[0,.5],xtit='7814 Error',ytit='6514 Error',$
		  tit='Old vs. New (acorr) error for MC and propagated error calc'
	oplot, prop_new, prop_old, psym=2, color=cgcolor('blue')
	oplot, mc_new, mc_old, psym=2, color=cgcolor('red')
	xyouts,.1,.4,'Propagated Error', color=cgcolor('blue')
	xyouts,.1,.38,'MC sim Error (sig68)', color=cgcolor('red')
	device,/close

	set_plot,'X'
stop
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;		sigs_b[*,*,i] = sqrt(rebin(sigs_c[*,*,i]^2,num_x,num_y))/sqrt(bx*by)

;		for j=0,num_x-1 do begin
;			for k=0,num_y-1 do begin
;				spec_b[j,k,i] = total(spec[lx+j*bx:lx+(j+1)*bx-1, $
;										ly+k*by:ly+(k+1)*by-1,i])
;			endfor
;		endfor


