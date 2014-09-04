pro crop, imgs, lambda1, lambda2
	for i=0,n_elements(imgs)-1 do begin
		; read in a science frame
		temp=mrdfits(imgs[i],0,h)
		name=strsplit(imgs[i],'/.',/extract)
		ss = size(temp)
		; generate wavelength range for this science frame
		wlen = dindgen(ss[1])*sxpar(h,'CD1_1')+sxpar(h,'CRVAL1')- $
		   		sxpar(h,'CRPIX1')*sxpar(h,'CD1_1')
		; define x range to crop to
		xrange=where(wlen gt lambda1 and wlen lt lambda2)
		; crop the science frame
		img_crop = temp[xrange,4:73]
		; save this file for later 
		writefits,'stacked/715/crop/'+name[2]+'_crop.fits',img_crop,h
		;writefits,'stacked/'+name[1]+'_crop.fits',img_crop,h
	endfor
end

; Rebins the data into identical bins for the entire image
pro rebin, imgs, x, y, SIG=sig
	; define the proper wavelength range
	for i=0,n_elements(imgs)-1 do begin
		; read in a science frame
		temp=mrdfits(imgs[i],0,h)
		ss=size(temp)
		name=strsplit(imgs[i],'/.',/extract)
		; rebin the science frame
		if keyword_set(sig) then img_bin = sqrt(rebin(temp^2,ss[1]/x,ss[2]/y))/sqrt(x*y) $
			else img_bin = rebin(temp,ss[1]/x,ss[2]/y)
		; save this file for later 
		writefits,'stacked/715/rebin/'+name[3]+'_rebin'+strtrim(fix(x),2)+'x'+ $
				  strtrim(fix(y),2)+'.fits',img_bin,h
	endfor
end

; Rebins the data in bins where the wlen direction is identical but 
; the spatial direction is allowed to vary - 
;		y must be a vector of pixel values to bin over
;		x must be a scalar defining the bin size in the wlen direction
pro rebins, imgs, y, x, SIG=sig
	for i=0,n_elements(imgs)-1 do begin
		; read in a science frame
		temp=mrdfits(imgs[i],0,h)
		ss=size(temp)
		img_bin = fltarr(ss[1]/x,n_elements(y)-1)
		name=strsplit(imgs[i],'/.',/extract)

		; if these are sigma maps then use error propagation...
		;
		if KEYWORD_SET(SIG) then begin
			; loop through the wlen direction binning in fixed size
			for j=0,ss[1]/x-1 do begin
				; loop through the spatial direction binning in particular sizes
				for k=0,n_elements(y)-2 do begin
					if k eq 0 then begin
						pixels=x*(y[k+1]-y[k]+1)				
						img_bin[j,k] = $
						sqrt(total(temp[j*x:(j+1)*x-1,y[k]:y[k+1]]^2))/pixels				
					endif else begin	
						pixels=x*(y[k+1]-y[k])				
						img_bin[j,k] = $
						sqrt(total(temp[j*x:(j+1)*x-1,y[k]+1:y[k+1]]^2))/pixels
					endelse
				endfor
			endfor
		; ...otherwise they're regular spectra so just take the mean
		;
		endif else begin
			; loop through the wlen direction binning in fixed size
			for j=0,ss[1]/x-1 do begin
				; loop through the spatial direction binning in particular sizes
				for k=0,n_elements(y)-2 do begin
					if k eq 0 then $
						img_bin[j,k] = mean(temp[j*x:(j+1)*x-1,y[k]:y[k+1]]) $
					else img_bin[j,k] = mean(temp[j*x:(j+1)*x-1,y[k]+1:y[k+1]])
				endfor
			endfor
		endelse
		writefits,'stacked/715/rebin/'+name[3]+'_rebin'+strtrim(fix(x),2)+'x'+ $
				  'V_bb8.fits',img_bin,h
	endfor
end

pro makepolfigs, images, x, y;, suffix;, err, tot
	pimg=fltarr(160, 70)
	timg=pimg
	pe68=pimg
	pe95=pimg
	
	for i=0,n_elements(images)-1 do begin
		img=mrdfits(images[i])
		name=strsplit(images[i],'.',/extract)
		; loop through the wlen direction binning in fixed size
		for j=0,160./x-1 do begin
			; loop through the spatial direction binning in particular sizes
			for k=0,n_elements(y)-2 do begin
				if k eq 0 then begin
					pimg[j*x:(j+1)*x-1,y[k]:y[k+1]] = img[j,k]
					;timg[j*x:(j+1)*x-1,y[k]:y[k+1]] = tot[j,k]
					;pe68[j*x:(j+1)*x-1,y[k]:y[k+1]] = err[j,k,*,0]	
					;pe95[j*x:(j+1)*x-1,y[k]:y[k+1]] = err[j,k,*,1]
				endif else begin	
					pimg[j*x:(j+1)*x-1,y[k]+1:y[k+1]] = img[j,k]
					;timg[j*x:(j+1)*x-1,y[k]+1:y[k+1]] = tot[j,k]
					;pe68[j*x:(j+1)*x-1,y[k]+1:y[k+1]] = err[j,k,*,0]	
					;pe95[j*x:(j+1)*x-1,y[k]+1:y[k+1]] = err[j,k,*,1]
				endelse
			endfor
		endfor
		writefits,name[0]+'_resized.fits',pimg
	endfor	
	; WRITE 2D DATA TO FILE
;	writefits,'stacked/715/rebin/Iall_avg_2d_crop_rebin8xV_'+suffix+'.fits',timg
;	writefits,'pol_products/715/MC_pstn68_715_mysig_8xV_'+suffix+'.fits',pimg/pe68
;	writefits,'pol_products/715/MC_pstn95_715_mysig_8xV_'+suffix+'.fits',pimg/pe95
;	writefits,'pol_products/715/MC_sig68_715_mysig_8xV_'+suffix+'.fits',pe68
;	writefits,'pol_products/715/MC_sig95_715_mysig_8xV_'+suffix+'.fits',pe95
;	writefits,'pol_products/715/pol_715_mysig_8xV_'+suffix+'.fits',pimg
end

pro extract, imgs
	for i=0,n_elements(imgs)-1 do begin
		; read in a science frame
		temp=mrdfits(imgs[i],0,h)
		ss=size(temp)
		name=strsplit(imgs[i],'/.',/extract)
		
		bot = fltarr(ss[1])
		top = bot
		; extract 2 1d spectra: top "blob" and bottom "blob"
		for j=0,ss[1]-1 do begin
			; 1d tot intensity for the bottom bright line
			bot[i] = total(temp[i,11:28])
			; 1d tot intensity for the top bright line
			top[i] = total(temp[i,35:54])
		endfor
		; save this file for later 
		writefits,'stacked/715/oned/'+name[3]+'_1dtop.fits',top,h
		writefits,'stacked/715/oned/'+name[3]+'_1dbot.fits',bot,h
	endfor
end

; calculate polarization the ESO/Claudia/Matt way
function polcalc1, flux
	q = 0.5*(flux[0]-flux[4])/(flux[0]+flux[4])	- $
		0.5*(flux[2]-flux[6])/(flux[2]+flux[6])
	u = 0.5*(flux[1]-flux[5])/(flux[1]+flux[5]) - $
		0.5*(flux[3]-flux[7])/(flux[3]+flux[7])
	theta = 0.5*ATAN(u/q)
	return, [sqrt(q^2+u^2),theta]
end

function thetacalc, flux
	q = 0.5*(flux[0]-flux[4])/(flux[0]+flux[4])	- $
		0.5*(flux[2]-flux[6])/(flux[2]+flux[6])
	u = 0.5*(flux[1]-flux[5])/(flux[1]+flux[5]) - $
		0.5*(flux[3]-flux[7])/(flux[3]+flux[7])
;	arctan,q,u,a,a_deg;*180/3.14159 	; pol angle in degrees
	theta1 = 0.5*a
	theta2 = 0.5*ATAN(u/q);*180./!pi
	return, [theta1,theta2]
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
function mc_sim, hwps, sigs, pol, coord, name, PLOT=PLOT
	n = 10000
	sims = fltarr(n_elements(hwps))
	sim_p = fltarr(n,2)
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
		; calculate the simulated polarization and theta
		sim_p(j-1,*) = polcalc1(sims)
		fr=where(ff eq j)
        ;if fr[0] gt 0. then print,'done ',ff[fr]/10000.,'%'
	endfor
	; calculate statistics typical of a symmetric distribution
	avg = [mean(sim_p[*,0]), mean(sim_p[*,1])]
	std = [stddev(sim_p[*,0]), stddev(sim_p[*,1])]
	; calculate statistics for a skewed distribution
	med = [median(sim_p[*,0]), median(sim_p[*,1])]
	percent = [[PERCENTILES(sim_p[*,0],VALUE=[.025,.16,.5,.84,.975,.0015,.9985])], $ ;pol
				[PERCENTILES(sim_p[*,1],VALUE=[.025,.16,.5,.84,.975,.0015,.9985])]]	; theta
	; define the error as the larger of the two absolute deviations from
	; the measurement 
	sigma1 = [(pol[0] - percent[1,0]) > (percent[3,0] - pol[0]), $	; pol 1sig error
				(pol[1] - percent[1,1]) > (percent[3,1] - pol[1])]	; theta 1sig err
	sigma2 = [(pol[0] - percent[0,0]) > (percent[4,0] - pol[0]), $
				(pol[1] - percent[0,1]) > (percent[4,1] - pol[1])]

	; If requested, plot the polarization histogram and various statistics
	if KEYWORD_SET(PLOT) then begin
		set_plot,'PS'
		if n_elements(coord) eq 2 then begin
			device,filename='plots/715/hist_715_'+name+'_'+strtrim(coord[0],2)+'_'+ $
						strtrim(coord[1],2)+'.ps',/color,/encap
			title = 'Pol Histrogram for bin '+strtrim(coord[0],2)+','+ $
				 		strtrim(coord[1],2)
		endif else begin
			device,filename='plots/715/hist_715_'+name+'_'+strtrim(coord[0],2)+'.ps'
			title = 'Pol Histrogram for bin '+strtrim(coord[0],2)		
		endelse		
;		device,filename=name,/color,/encap
		; plot the histogram of simulated polarizations
		plothist, sim_p[*,0], xhist, yhist, bin=0.01, /fill, xr=[0,1.], $
				  xtit='P Frac', $
				  charsize=1.5, charthick=4, thick=1, xthick=4, ythick=4, $
				  fcolor=cgcolor('yellow'), peak=500.
		; overplot the mean and median of the simulated polarization dist
		oplot,yhist*0.+avg[0],yhist,thick=3,color=cgcolor('red')
		;oplot,yhist*0.+med[0],yhist,thick=3,color=cgcolor('green')

		; overplot the original measured polarization
		oplot,yhist*0.+pol[0],yhist,thick=4,color=cgcolor('black') 

		; overplot the avg sim pol +- std of the histogram
		oplot,yhist*0.+(avg[0]+std[0]),yhist,thick=3,linestyle=1,color=cgcolor('blue')
		oplot,yhist*0.+(avg[0]-std[0]),yhist,thick=3,linestyle=1,color=cgcolor('blue')

		; overplot the avg sim pol +- 2*std of the histogram
		oplot,yhist*0.+(avg[0]+2*std[0]),yhist,thick=3,linestyle=1, $
			  color=cgcolor('magenta')
		oplot,yhist*0.+(avg[0]-2*std[0]),yhist,thick=3,linestyle=1, $
			  color=cgcolor('magenta')

		; overplot the percentiles of the histogram
		oplot,yhist*0.+percent[2,0],yhist,thick=4,linestyle=0, $	; 50%
			  color=cgcolor('green')
    	oplot,yhist*0.+percent[1,0],yhist,thick=4,linestyle=2, $	; 16%
			  color=cgcolor('blue')
   		oplot,yhist*0.+percent[3,0],yhist,thick=4,linestyle=2, $	; 84%
			  color=cgcolor('blue')
 		oplot,yhist*0.+percent[0,0],yhist,thick=4,linestyle=2, $	; 2.5%
			  color=cgcolor('magenta')	
 		oplot,yhist*0.+percent[4,0],yhist,thick=4,linestyle=2, $	; 97.5%
			  color=cgcolor('magenta')
 		;oplot,yhist*0.+percent[5,0],yhist,thick=4,linestyle=2, $	; 2.5%
		;	  color=cgcolor('green')	
 		;oplot,yhist*0.+percent[6,0],yhist,thick=4,linestyle=2, $	; 97.5%
		;	  color=cgcolor('green')

		items = ['Raw P', 'Mean', 'Median', '1 Std Dev (68%)', '2 Std Dev (95%)', $
				 'Percentiles (68%)', 'Percentiles (95%)']

		;items = ['Raw P', 'Median', '1Sig (68%)', '2Sig (95%)']
;		items = ['Raw P', 'Mean', 'Median', '1 Sigma (68%)', '2 Sigma (95%)']
		lines = [0,0,0,1,1,2,2]
		colors = [cgcolor('black'), cgcolor('red'), cgcolor('green'), cgcolor('blue'), $
				  cgcolor('magenta'), cgcolor('blue'), cgcolor('magenta')] 

		;colors = [cgcolor('black'), cgcolor('red'), cgcolor('blue'), $
		;		  cgcolor('magenta'), cgcolor('green')]
		al_legend, items, linestyle=lines, color=colors, linsize=.5, /right, $
				   background_color=cgcolor('white'), bthick=4, charsize=1.2, $
				   charthick=3, thick=4
		device,/close
		set_plot,'X
;stop
	endif
	; return SIGMA1 and SIGMA2 
	return, [[sigma1], [sigma2]]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro pol2d, date, CROP=crop, REBIN=rebin, REBINS=rebins, EXTRACT=extract, MYSIG=mysig
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; INITIALIZE SOME SHITS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; choose bin sizes
	bx = 8.
	by = 7.
	
	; this is from the previous procedure: 
	; bin the spectra in the spatial direction the same way as in that figure
	;ts=[0,5,11,18,26,35,44,53,61,69]	
	ts=[0,11,18,26,35,42,50,57,69]	; 4th attempt -- labelled '_bb8.fits'

	; define the proper wavelength range
	start_wlen = 4930.
	end_wlen = 5032.

	; read in a test file to determine wavelength range(s)
	temp = mrdfits('stacked/'+date+'/Iall_stack_'+date+'_2d.fits',0,h)
	ss = size(temp)

	; DETERMINE X RANGE - NoooUMBER OF PIXELS IN THIS RANGE
	wlen = dindgen(ss[1])*sxpar(h,'CD1_1')+sxpar(h,'CRVAL1')- $
		   sxpar(h,'CRPIX1')*sxpar(h,'CD1_1')
	xrange=where(wlen gt start_wlen and wlen lt end_wlen)
	lx = xrange[0]
	num_x = (n_elements(xrange))/bx	
	
	; DETERMINE Y RANGE - NUMBER OF PIXELS IN THIS RANGE
	y = dindgen(ss[2])
	yrange=where(y gt 3 and y lt 74)
	ly = yrange[0]
	num_y = (n_elements(yrange))/by

	num_y = n_elements(ts)-1

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; DECLARE SOME MOTHERFUCKIN VARIABLES
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	science = fltarr(num_x, num_y, 8)
	sigmas = fltarr(num_x, num_y, 8)
	pol = fltarr(num_x, num_y, 2)
	perr = fltarr(num_x, num_y, 2, 2)

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; CROP 2D SPECTRA 
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	if KEYWORD_SET(CROP) then begin
;		spawn,'ls stacked/715/[I,A,B]*stack_715_2d_*sum2.fits > stacked.list'
;		spawn,'ls stacked/715/[I,A,B]*stack_715_2d_sc.fits > stacked.list'
		spawn,'ls stacked/Iall*.fits > stacked.list'
		readcol,'stacked.list',imgs,format='a'
		crop, imgs, start_wlen, end_wlen
;		spawn,'ls stacked/715/[I,A,B]*stack_715_2d_sc_sig.fits > stacked.list'
		spawn,'ls stacked/Iall*_sig.fits > stacked.list'
		readcol,'stacked.list',imgs,format='a'
		crop, imgs, start_wlen, end_wlen
	endif
;stop
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; REBIN 2D SPECTRA 
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	if KEYWORD_SET(REBIN) then begin
		; rebin science frames
;		spawn,'ls stacked/715/crop/[I,A,B]*2d_*sum2_crop.fits > stacked.list'
		spawn,'ls stacked/715/crop/[I,A,B]*2d_crop.fits > stacked.list'
		readcol,'stacked.list',imgs,format='a'
		rebin, imgs, bx, by

		; rebin IRAF sigma maps
;		spawn,'ls stacked/715/crop/[I,A,B]*2d_sig_sum2_crop.fits > stacked.list'
		spawn,'ls stacked/715/crop/[I,A,B]*2d_sig_crop.fits > stacked.list'
		readcol,'stacked.list',imgs,format='a'
		rebin, imgs, bx, by, /sig

		; read in cropped sigma frames - MY SIGMA MAPS
		spawn,'ls stacked/715/crop/[A,B]*_2d_MYsig_avg.fits > stacked.list'
		readcol,'stacked.list',imgs,format='a'
		rebin, imgs, bx, by, /sig
	endif

	if KEYWORD_SET(REBINS) then begin
		spawn,'ls stacked/715/crop/[I,A,B]*2d_crop.fits > stacked.list'
		readcol,'stacked.list',imgs,format='a'
		rebins, imgs, ts, bx

		spawn,'ls stacked/715/crop/[I,A,B]*_2d_MYsig_avg.fits > stacked.list'
		readcol,'stacked.list',imgs,format='a'
		rebins, imgs, ts, bx, /sig
	endif
;stop
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; CREATE 1D SPECTRA -- TOP AND BOTTOM
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	if KEYWORD_SET(EXTRACT) then begin
;		spawn,'ls stacked/715/crop/Iall_stack*sum2_crop.fits > stacked.list'
		spawn,'ls stacked/715/crop/Iall_stack_2d_sc_crop.fits > stacked.list'
		readcol,'stacked.list',totes,format='a'
		extract,totes

;		spawn,'ls stacked/715/crop/[A,B]*stack*crop.fits > stacked.list'
		spawn,'ls stacked/715/crop/[A,B]*stack_2d_sc_crop.fits > stacked.list'
		readcol,'stacked.list',imgs,format='a'
		extract,imgs
	endif
;stop
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; CALCULATE MEASURED POLARIZATION IN BINNED FRAMES
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	spawn,'ls stacked/715/rebin/[A,B]*2d_sum2*rebin.fits > stacked.list'
	spawn,'ls stacked/715/rebin/[A,B]*2d_crop_rebin8xV_bb8.fits > stacked.list'
	readcol,'stacked.list',imgs,format='a'

;	spawn,'ls stacked/715/rebin/[A,B]*2d_sig_sum2*rebin.fits > stacked.list'
	if KEYWORD_SET(MYSIG) then $
		spawn,'ls stacked/715/rebin/[A,B]*2d_MYsig_avg_rebin8xV_bb8.fits > stacked.list' $
	else spawn,'ls stacked/715/rebin/[A,B]*2d_sig_crop_rebin8xV_bb8.fits > stacked.list'
	readcol,'stacked.list',sigs,format='a'

	; read in science and sigma map frames
	for i=0,n_elements(imgs)-1 do begin
		science[*,*,i] = mrdfits(imgs[i])
		sigmas[*,*,i] = mrdfits(sigs[i])
	endfor

	; calculate polarization for each bin
	for i=0,num_x-1 do begin
		for j=0,num_y-1 do begin
			; calculate the polarization for the 2D data
			pol[i,j,*] = polcalc1(science[i,j,*])
			if i ge 8 and i le 13 then perr[i,j,*,*] = mc_sim(science[i,j,*], $
				sigmas[i,j,*], pol[i,j,*], [i,j], 'bb8', /plot) $
			else perr[i,j,*,*] = mc_sim(science[i,j,*], sigmas[i,j,*], $
				pol[i,j,*], [i,j])
		endfor
	endfor
	
	if KEYWORD_SET(MYSIG) then begin
	; WRITE 2D DATA TO FILE
	writefits,'pol_products/715/MC_pstn68_715_mysig_8xV_bb8.fits',pol[*,*,0]/perr[*,*,0,0]
	writefits,'pol_products/715/MC_pstn95_715_mysig_8xV_bb8.fits',pol[*,*,0]/perr[*,*,0,1]
	writefits,'pol_products/715/MC_sig68_715_mysig_8xV_bb8.fits',perr[*,*,0,0]
	writefits,'pol_products/715/MC_sig95_715_mysig_8xV_bb8.fits',perr[*,*,0,1]
	writefits,'pol_products/715/pol_715_mysig_8xV_bb8.fits',pol[*,*,0]
	endif else begin
	writefits,'pol_products/715/MC_pstn68_715_IRAF.fits',pol/perr[*,*,0]
	writefits,'pol_products/715/MC_pstn95_715_IRAF.fits',pol/perr[*,*,1]
	writefits,'pol_products/715/MC_sig68_715_IRAF.fits',perr[*,*,0]
	writefits,'pol_products/715/MC_sig95_715_IRAF.fits',perr[*,*,1]
	writefits,'pol_products/715/pol_715_IRAF.fits',pol
	endelse

	;tot=mrdfits('stacked/715/rebin/I*_avg_2d_crop_rebin8xV_bb8.fits')
	spawn,'ls stacked/715/rebin/I*_rebin8xV_bb8.fits > resize.list'
	readcol,'resize.list',imgs,format='a'
	; make a pol figure that looks like it has varying sized boxes in y dir
	makepolfigs,imgs, bx, ts

	spawn,'ls pol_products/715/*8xV_bb8.fits > resize.list'
	readcol,'resize.list',imgs,format='a'
	makepolfigs,imgs, bx, ts

	; WRITE 1D DATA TO FILE
;	writefits,'pol_products/MC_polstn68_'+date+'_1d.fits',pol1d/perr1d[*,0]
;	writefits,'pol_products/MC_polstn95_'+date+'_1d.fits',pol1d/perr1d[*,1]
;	writefits,'pol_products/MC_sigma68_'+date+'_1d.fits',perr1d[*,0]
;	writefits,'pol_products/MC_sigma95_'+date+'_1d.fits',perr1d[*,1]
;	writefits,'pol_products/pol_'+date+'_1d.fits',pol1d

stop

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro pol_bigboxes, MYSIG=mysig, SAVEDAT=savedat, ONED=oned
	; access the cropped science frames
	spawn,'ls stacked/715/crop/[A,B]*2d_crop.fits > stacked.list'
	readcol,'stacked.list',imgs,format='a'
	;if KEYWORD_SET(MYSIG) then $
	spawn,'ls stacked/715/crop/[A,B]*2d_MYsig_avg.fits > stacked.list' 
	;else spawn,'ls stacked/715/crop/[A,B]*2d_sig_crop.fits > stacked.list'
	readcol,'stacked.list',sigs,format='a'

	science=fltarr(160,70,8)
	sigmas=science

	; read in science and sigma map frames
	for i=0,n_elements(imgs)-1 do begin
		science[*,*,i] = mrdfits(imgs[i])
		sigmas[*,*,i] = mrdfits(sigs[i])
	endfor
	
	; read in I00 total intensity frames
	i00=mrdfits('stacked/715/crop/I00_stack_715_2d_crop.fits')
	i00s=mrdfits('stacked/715/crop/I00_stack_715_2d_MYsig_avg.fits')

	; SPECIFY THE LEFT AND RIGHT HAND WLEN RANGES (IN PIXELS) FOR EACH BOX 
	;xl=[55,55,55,55,65,65,65,65,65,65]
	;xr=[110,110,110,110,105,105,105,105,105,105]
	;xl=[58,58,58,58,58,58,58,58]
	xl=[55,55,55,55,55,55,55,55]
	xr=[110,110,110,110,110,110,110,110]
	;xr=[106,106,106,106,106,106,106,106]
	; AND THEN SPECIFY THE TOPS (IN PIXELS) OF THE BOXES FOR THE SPATIAL DIRECTION
	;ts=[7,13,19,26,34,44,54,61,69]
	;ts=[5,11,18,26,35,44,53,61,69]
	;ts=[5,11,18,25,33,40,48,55,62,69]	; adjusted to match matt's output
	;ts=[11,18,26,34,44,54,69]	; 3rd attempt after reading W10 paper
	ts=[11,18,26,35,42,50,57,69]	; 4th attempt 
	bs=[0,ts[0:n_elements(ts)-1]]

	; declare some variables
	tot=fltarr(n_elements(ts))
	tots=tot
	flux=fltarr(n_elements(ts),8)
	sig=flux
	; contains pol and theta measurements
	pol=fltarr(n_elements(ts),2)	
	theta=pol
	; containts 1 & 2 sigma errors for pol measurements
	perr=fltarr(n_elements(ts),2,2)

	; integrate over the entire bin
	for i=0,n_elements(imgs)-1 do begin
		for j=0,n_elements(ts)-1 do begin
			if j eq 0 then begin
				tot[j] = total(i00[xl[j]:xr[j],bs[j]:ts[j]])
				tots[j] = total(i00s[xl[j]:xr[j],bs[j]:ts[j]]^2)
				flux[j,i] = total(science[xl[j]:xr[j],bs[j]:ts[j],i])
				sig[j,i] = sqrt(total(sigmas[xl[j]:xr[j],bs[j]:ts[j],i]^2))
			; after doing the first box, you don't want to use a pixel value
			; from the previous box so add 1 to the lower bound in the y direction
			endif else begin
				tot[j] = total(i00[xl[j]:xr[j],bs[j]+1:ts[j]])
				tots[j] = total(i00s[xl[j]:xr[j],bs[j]+1:ts[j]]^2)
				flux[j,i] = total(science[xl[j]:xr[j],bs[j]+1:ts[j],i])
				sig[j,i] = sqrt(total(sigmas[xl[j]:xr[j],bs[j]+1:ts[j],i]^2))
			endelse
		endfor
	endfor
	stn00 = tot/tots
;stop

	; calculate polarization for each bin
	for i=0,n_elements(ts)-1 do begin
		pol[i,*] = polcalc1(flux[i,*])
		;theta[i,*] = thetacalc(flux[i,*])
		perr[i,*,*] = mc_sim(flux[i,*], sig[i,*], pol[i,*])
	endfor
	pol[*,1]=pol[*,1]*180./!pi

	if KEYWORD_SET(SAVEDAT) then begin
		openw,out,'spectropol_bb8_wlensame.dat',/get_lun
		printf,out,'#	P%	   Perr%      THETA	    Terr'
		for i=0,n_elements(pol[*,1])-1 do begin
		printf,out,[pol[i,0], perr[i,0,1],pol[i,1], $
			   perr[i,0,0]/(2*pol[i,0])*180./!pi], $
			   format='(4(d12.5))'
		endfor
		close,out
		free_lun,out
	endif

	LETTERS=['A','B','C','D','E','F','G','H','I','J']

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; MAKE A PRETTY FIGURE
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	tot=mrdfits('stacked/Iall_avg_crop.fits')
	ss=size(tot)

	tot1d=mrdfits('stacked/715/Iall_stack_715_2d.fits',0,h)
	ss=size(tot1d)
	wlen=dindgen(ss[1])*sxpar(h,'CD1_1') + sxpar(h,'CRVAL1')- $
		 sxpar(h,'CRPIX1')*sxpar(h,'CD1_1')
	wlen2=wlen[where(wlen gt 4930. and wlen lt 5032.)]

	; READ IN MATT'S POLARIZATION DATA
	readcol, 'poldat_list.txt', pix, p, dp68, dp95

	; READ IN MATT'S LYA FLUX DATA
	readcol, 'lab1_slitflux.dat', pixel, flux


	pos_i=[0.2,0.46,0.9,0.95]
	pos_p=[0.2,0.15,0.9,0.45]

	; display the cropped/transpose/smoothed total intensity image
	loadct,0
	tot2=tot[xl[0]:xr[0],*]
	ss=size(tot2)

	; FLUX TEST
	; for each pixel that matt sampled, sum up the the wlen range for the 
	; corresponding pixel in our slit
	lyaf=fltarr(n_elements(pixel))
	for i=0,n_elements(pixel)-1 do begin
		lyaf[i] = total(tot2[*,pixel[i]])
	endfor

	set_plot,'PS'

	;device,filename='lyaflux_test.ps',/color,/encap
	;plot, pixel, flux/max(flux), yr=[0,1.1], xr=[0,70]
	;oplot, pixel, lyaf/max(lyaf), color=cgcolor('red')
	;xyouts,50, 1.1, 'Imaging', charthick=2, charsize=1.5
	;xyouts,50, 1.0, 'Spectrum', charthick=2, charsize=1.5, color=cgcolor('red')
	;device,/close

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	device,filename='bigboxes8_wlensame_2sig.ps',/color,/encap

	totsm = transpose(gauss_smooth(tot2,.85,/edge_truncate))
;	totsm = transpose(tot2)
	cgimage,bytscl(-totsm,min=-0.006,max=0.0009),POSITION=pos_i,/noerase;,/keep_aspect_ratio
	loadct,0

	; create x and y axes for the lya image
	x1=dindgen(ss[2]+1)*.25
	y1=dindgen(ss[1])
	; plot the axes for the cropped lya image on it's side
	plot,x1,y1,POSITION=pos_i,/noerase,/nodata,thick=4,xthick=4,ythick=4, $
		 charthick=2,charsize=1.,xstyle=1,ystyle=1,yr=[0.,51], xr=[0,17.5], $
		 yticklen=0.01, yticks=4, xticklen=0.025,  $
		 ytickformat='(a2)', xtickformat='(a2)'

	; define the starting point for plotting the red lines that are the boxes i
	; summed over
	j=ts[0]
	; create an x axis for plotting the polarization values
	x = [0]
	; create an offsets array for plotting the pol in the middle of bins
	offsets = [0]

	; plot the red lines to delineate the boxes i summed over
	for i=0,n_elements(ts)-1 do begin
		x=[x,x1[j+1]]
		offsets=[offsets,float(x[i+1]-x[i])/2.]
		; if it's not the last red line, plot it
		if i lt n_elements(ts)-1 then begin 
			; here we do j+1 because there's a shift of one pixel between IDL 
			; and everything else (0-start indexing vs 1-start indexing)
			plots, [x1[j+1], x1[j+1]], [0, 51], /DATA,thick=4,color=cgcolor('red')
			xyouts,x1[j+1]-offsets[i+1], 45, LETTERS[i], color=cgcolor('red'), $
					charsize=1.5,charthick=4
		endif

		if i eq n_elements(ts)-1 then $
			xyouts,x1[j+1]-offsets[i+1], 45, LETTERS[i], color=cgcolor('red'), $
					charsize=1.5,charthick=4
		if i le n_elements(ts)-2 then j+=bs[i+2]-bs[i+1] else j=ts[i]
	endfor
	; the first offset isn't really 0 so get rid of that
	offsets = offsets[1:n_elements(offsets)-1]

	; overlay the wavelength axis on the cropped lya image	
	axis, yaxis=0, yr=[wlen2[xl[0]],wlen2[xr[0]]], xstyle=1, charthick=4, $
		  charsize=1.25, ythick=4, yticklayout=1, ytit='Wavelength';, yticks=4

	; plot the polarization axes below the lya image
	y2=dindgen(15)/100.
	plot,x,y2,POSITION=pos_p,/noerase,/nodata,thick=4,xthick=4,ythick=4, $
		 charthick=4,charsize=1.25,xstyle=1,ystyle=1,yr=[-0.1,.3], $
		 xr=[0,17.5], yticklen=0.01, yticks=5, xticklen=0.025, $
		 xtit='Arcseconds',ytit='Frac. P', ytickinterval=.1
	oplot,x,dindgen(100)*0, linestyle=2,thick=3

	;oplot, pix*25, p, psym=5, thick=3, color=cgcolor('blue')
	oploterror,pix*.25, p, dp68, linestyle=6, psym=5, errthick=3, thick=3

	; plot the fractional polarization points using the x axis created during 
	; the plotting of the red lines
	;oplot,x+offsets,pol,psym=6,thick=3, symsize=1, color=cgcolor('blue')
	oploterror,x+offsets,pol[*,0],perr[*,0,1],linestyle=6,psym=6,symsize=1, $
				errthick=4, thick=4, color=cgcolor('orange red')

	device,/close
	set_plot,'X'



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;; FUCKING AROUND WITH 1D STUFF FOR FUNSIES ;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	if KEYWORD_SET(ONED) then begin
		filename='stacked/715/crop/Iall_avg_2d_crop.fits'
		; extract 1D line profiles for EACH bin
		tot=mrdfits(filename)
		name=strsplit(filename,'.',/extract)
		ss=size(tot)
		tot1ds=fltarr(ss[1],n_elements(ts))
		for i=0,ss[1]-1 do begin
			for j=0,n_elements(ts)-1 do begin
				if j eq 0 then begin
					tot1ds[i,j]=total(tot[i,bs[j]:ts[j]])
				endif else begin
					tot1ds[i,j]=total(tot[i,bs[j]+1:ts[j]])
				endelse
			endfor
		endfor

		; read in 2D binned pol data
		pol_file='pol_products/715/pol_715_mysig_8xV_bb8.fits'
	    stn_file='pol_products/715/MC_pstn95_715_mysig_8xV_bb8.fits'
		err_file='pol_products/715/MC_sig68_715_mysig_8xV_bb8.fits'

		i00=mrdfits('stacked/715/rebin/I00_stack_715_2d_crop_rebin8xV_bb8.fits')
		i00s=mrdfits('stacked/715/rebin/I00_stack_715_2d_MYsig_avg_rebin8xV_bb8.fits')
		pola=mrdfits(pol_file)
		stna=mrdfits(stn_file)
		erra=mrdfits(err_file)	
	
		stn00 = i00/i00s
		pstncut = 1.0
		stn00cut = 0.0
	
		; remove any nans
		nan=finite(pol,/nan)
		pola[where(nan eq 1)]=0.
		; ceate array where only pol above stn cut remains
		pol_sig=pola*0.+10.
		err_sig=pola*0.+10.
		pol_sig[where(stna gt pstncut and stn00 gt stn00cut)]= $
			pola[where(stna gt pstncut and stn00 gt stn00cut)]
		err_sig[where(stna gt pstncut and stn00 gt stn00cut)]= $
			erra[where(stna gt pstncut and stn00 gt stn00cut)]
		;pol_sig[where(stn00 gt stn00cut)] = pol[where(stn00 gt stn00cut)]
		pol_sig[where(pol_sig gt 0.5)]=10.
		err_sig[where(pol_sig gt 0.5)]=0.
	
		x=dindgen(ss[1])	; x axis in pixel coordinates for plotting the 1d spec
		x2=dindgen(ss[1]/8.)*8.-5 ; x axis in bins of 8 pixels each
		x2=x2[1:n_elements(x2)-1]
		x2=[x2,x2[n_elements(x2)-1]+8.]
		wlen = x*0.636 + 4930. ; x axis in wavelength coordinates

		; SET POSITION ARGUMENTS
		pos = [[0.15,0.71,0.85,0.91], $
			   [0.15,0.51,0.85,0.71], $
			   [0.15,0.31,0.85,0.51], $
			   [0.15,0.11,0.85,0.31]]

		circsym

		set_plot,'PS'
		device,filename='multi1d_bb8_wlenfixed_kms.ps',/encap,/color
	
		for i=0,n_elements(ts)-1 do begin
			; plot the 1d spec IF there is appreciable polarization for that aperture!
			if total(pol_sig[*,i]) lt 200. then begin
				if i lt 4 then begin
				;print,total(pol_sig[*,i])
				plot, x[20:140], tot1ds[20:140,i]/max(tot1ds[*,i]), POSITION=$
						pos[*,i],/noerase, xtickformat='(a2)', thick=3, $ 
						charsize=1.2, charthick=3, xthick=3, ythick=3,  $
						yticks=3, yr=[-0.1, 1.1], ystyle=9, xstyle=9, /data, $
						ytickinterval=.4, ytickformat='(F3.1)'
				xyouts,30, .70, LETTERS[i], charthick=4, charsize=1.2
				; oplot line center
				oplot, dindgen(10)*0.+81., [-0.5,1.0], linestyle=2, thick=3
				oplot, [20.,140.],dindgen(10)*0., linestyle=2, thick=2
				; oplot integration lines
				;oplot, dindgen(10)*0.+xl[i], [-0.5,1.0], linestyle=3, thick=3, $
				;		color=cgcolor('blue')
				;oplot, dindgen(10)*0.+xr[i], [-0.5,1.0], linestyle=3, thick=3, $
				;		color=cgcolor('blue')

				; PLOT CORRESPONDING POLARIZATION 2D ;;;;;;;;;;;;;;;;;;;
				axis, 140, yaxis=1, yrange=[-0.06,0.6], ystyle=1, charthick=3,  $
						charsize=1.2, ythick=3, color=cgcolor('orange red'), $
						ytickinterval=.2, /save
				oploterror, x2, pol_sig[*,i], err_sig[*,i], linestyle=6, psym=8, $
						symsize=3, errthick=3, color=cgcolor('orange red');, $
				;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
				endif else begin
		 		plot, x[20:140], tot1ds[20:140,i]/max(tot1ds[20:140,i]), POSITION=$
						pos[*,i-4], /noerase, xtickformat='(a2)', thick=3, $ 
						charsize=1.2, charthick=3, xthick=3, ythick=3, $
						yticks=3, yr=[-0.1, 1.1], ystyle=9, xstyle=9, /data, $
						ytickinterval=.4, ytickformat='(F3.1)'
				xyouts,30, .70, LETTERS[i], charthick=4, charsize=1.2
				; oplot line center
				oplot, dindgen(10)*0.+81., [-0.5,1.0], linestyle=2, thick=3
				oplot, [20.,140.],dindgen(10)*0., linestyle=2, thick=2
				; oplot integration lines
				;oplot, dindgen(10)*0.+xl[i], [-0.5,1.0], linestyle=3, thick=3, $
				;		color=cgcolor('blue')
				;oplot, dindgen(10)*0.+xr[i], [-0.5,1.0], linestyle=3, thick=3, $
				;		color=cgcolor('blue')

				; PLOT CORRESPONDING POLARIZATION 2D ;;;;;;;;;;;;;;;;
				axis, 140, yaxis=1, yrange=[-0.06,.6], ystyle=1, charthick=3,  $
						charsize=1.2, ythick=4, color=cgcolor('orange red'), $
						ytickinterval=.2, /save
				oploterror, x2, pol_sig[*,i], err_sig[*,i], linestyle=6, psym=8, $
						symsize=3, errthick=3, color=cgcolor('orange red');, $
				;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
				endelse
			if i eq 0 or i eq 4 then axis, xaxis=1, xr=[-2499.,2319.], charthick=3, $
					  charsize=1.2, xtit='delta v [km/s]', xstyle=1, xthick=3
			if (i+1) mod 4 eq 0 then begin
				axis, xaxis=0, xr=[4940., 5020.],charthick=3, xthick=3, $
					  charsize=1.2,  xtit='Wavelength [Angstroms]', xstyle=1
				xyouts,7,.35,'Intensity [Arbitrary Units]', orientation=90, $
					  charsize=1.5, charthick=3
				xyouts,153,2.1,'Polarization Fraction', orientation=270, $
					  charsize=1.5, charthick=3, color=cgcolor('orange red')
				erase
			endif
			endif
		endfor
		device,/close
		set_plot,'X'
	endif
stop


	; define appropriate boxes in which to calculate polarization
;	cb1 = [65,12,103,19]			; lower bottom
;	cb2 = [58,19,103,27]			; upper bottom

;	cb2 = [58,27,103,33]			; low intensity stuff between U and L
			
;	cb_num = (cb[2]-cb[1]+1)*(cb[3]-cb[1]+1)
;	ct1 = [65,35,98,45]				; lower top
;	ct2 = [65,47,98,57]				; upper top

;	plots,[ct2[0],ct2[2]],[ct2[1],ct2[1]]*.25,/DATA,thick=4,color=cgcolor('red')
;	plots,[ct2[0],ct2[2]],[ct2[3],ct2[3]]*.25,/DATA,thick=4,color=cgcolor('red')
;	plots,[ct2[0],ct2[0]],[ct2[1],ct2[3]]*.25,/DATA,thick=4,color=cgcolor('red')
;	plots,[ct2[2],ct2[2]],[ct2[1],ct2[3]]*.25,/DATA,thick=4,color=cgcolor('red')
;
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


